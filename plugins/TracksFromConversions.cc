// -*- C++ -*-
//
// Package:    CMGTools/TracksFromConversions
// Class:      TracksFromConversions
// 
/**\class TracksFromConversions TracksFromConversions.cc CMGTools/TracksFromConversions/plugins/TracksFromConversions.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Mauro Verzetti
//         Created:  Mon, 11 Jun 2018 15:24:14 GMT
//
//


// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/stream/EDProducer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/StreamID.h"

#include "DataFormats/Common/interface/TriggerResults.h"
#include "DataFormats/Common/interface/ValueMap.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectron.h"
#include "DataFormats/HLTReco/interface/TriggerEvent.h"
#include "DataFormats/Math/interface/deltaR.h"

#include "MagneticField/Engine/interface/MagneticField.h"
#include "MagneticField/Records/interface/IdealMagneticFieldRecord.h"
#include "MagneticField/ParametrizedEngine/src/OAEParametrizedMagneticField.h"

#include "RecoVertex/KalmanVertexFit/interface/KalmanVertexFitter.h"

#include "TrackingTools/TransientTrack/interface/TransientTrack.h"
#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"

#include "DataFormats/EgammaCandidates/interface/Conversion.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "CommonTools/Statistics/interface/ChiSquaredProbability.h"

//
// class declaration
//

class TracksFromConversions : public edm::stream::EDProducer<> {
   public:
      explicit TracksFromConversions(const edm::ParameterSet&);
      ~TracksFromConversions();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

   private:
      virtual void beginStream(edm::StreamID) override;
      virtual void produce(edm::Event&, const edm::EventSetup&) override;
      virtual void endStream() override;

      //virtual void beginRun(edm::Run const&, edm::EventSetup const&) override;
      //virtual void endRun(edm::Run const&, edm::EventSetup const&) override;
      //virtual void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;
      //virtual void endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;

      // ----------member data ---------------------------
	const edm::EDGetTokenT< std::vector<reco::Conversion> > src_;
	const edm::EDGetTokenT< reco::BeamSpot > beamspotToken_;
	const edm::EDGetTokenT< std::vector<reco::Track> > tracks_;
	
  bool generalTracksOnly_;
  bool arbitratedMerged_;
  bool arbitratedEcalSeeded_;
  bool ecalalgotracks_;
  bool highPurity_;
  double minProb_;
  uint maxHitsBeforeVtx_;
  double minLxy_;
  double minVtxR_;
  double maxVtxR_;
	double minLeadPt_;
};

//
// constants, enums and typedefs
//


//
// static data member definitions
//

//
// constructors and destructor
//
TracksFromConversions::TracksFromConversions(const edm::ParameterSet& cfg):
	src_{consumes< std::vector<reco::Conversion> >(cfg.getParameter<edm::InputTag>("src"))},
	beamspotToken_{consumes< reco::BeamSpot >(cfg.getParameter<edm::InputTag>("beamspot"))},
	tracks_{consumes< std::vector<reco::Track> >(cfg.getParameter<edm::InputTag>("tracks"))},
  generalTracksOnly_{cfg.getParameter<bool>("generalTracksOnly")},
  arbitratedMerged_{cfg.getParameter<bool>("arbitratedMerged")},
  arbitratedEcalSeeded_{cfg.getParameter<bool>("arbitratedEcalSeeded")},
  ecalalgotracks_{cfg.getParameter<bool>("ecalalgotracks")},
  highPurity_{cfg.getParameter<bool>("highPurity")},
  minProb_{cfg.getParameter<double>("minProb")},
  maxHitsBeforeVtx_{cfg.getParameter<uint>("maxHitsBeforeVtx")},
  minLxy_{cfg.getParameter<double>("minLxy")},
	minVtxR_{cfg.getParameter<double>("minVtxR")},
  maxVtxR_{cfg.getParameter<double>("maxVtxR")},
	minLeadPt_{cfg.getParameter<double>("minLeadPt")}
{
   //register your products
	produces<reco::TrackRefVector>("electrons");
	produces<reco::TrackRefVector>("NOTelectrons");  
}


TracksFromConversions::~TracksFromConversions()
{
 
   // do anything here that needs to be done at destruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called to produce the data  ------------
void
TracksFromConversions::produce(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
   using namespace edm;
	 std::unique_ptr<reco::TrackRefVector> elecs(    new reco::TrackRefVector()); 
	 std::unique_ptr<reco::TrackRefVector> not_elecs(new reco::TrackRefVector()); 

	 edm::Handle< std::vector<reco::Conversion> > conversions;
	 iEvent.getByToken(src_, conversions);

	 edm::Handle<reco::BeamSpot> beamspot;
	 iEvent.getByToken(beamspotToken_, beamspot);

	 edm::Handle< std::vector<reco::Track> > tracks;
	 iEvent.getByToken(tracks_, tracks);

	 std::set<reco::TrackRef> electrons_set;
	 
	 //std::cout << "Starting from " << conversions->size() << " conversions" << std::endl;
	 size_t passed[4] = {0,0,0,0};
	 for(auto& conv : *conversions) {
		 if( arbitratedMerged_ && !conv.quality(reco::Conversion::arbitratedMerged)  ) continue;
		 if( generalTracksOnly_ && !conv.quality(reco::Conversion::generalTracksOnly) ) continue;
		 if( arbitratedEcalSeeded_ && !conv.quality(reco::Conversion::arbitratedEcalSeeded)  ) continue;
		 if( highPurity_ && !conv.quality(reco::Conversion::highPurity) ) continue;
		 
		 passed[0]++;
		 const reco::Vertex& vtx = conv.conversionVertex();
		 if (conv.tracks().size() !=2 || !(vtx.isValid())) continue;

		 if (ChiSquaredProbability( conv.conversionVertex().chi2(),  conv.conversionVertex().ndof() ) <= minProb_) continue;
		 passed[1]++;

		 if (conv.nHitsBeforeVtx().size()>1 && std::max(conv.nHitsBeforeVtx().at(0), conv.nHitsBeforeVtx().at(1)) > maxHitsBeforeVtx_ ) continue;
		 
		 //compute transverse decay length with respect to beamspot
		 math::XYZVectorF  themom = conv.refittedPairMomentum();
		 double dbsx = conv.conversionVertex().x() - beamspot->x0();
		 double dbsy = conv.conversionVertex().y() - beamspot->y0();
		 double lxy = (themom.x()*dbsx + themom.y()*dbsy)/themom.rho();
		 if (lxy < minLxy_) continue;
		 passed[2]++;

		 double vtx_r = std::sqrt(conv.conversionVertex().position().perp2());
		 if(vtx_r < minVtxR_ || vtx_r > maxVtxR_) continue;

		 //fill set for later
		 electrons_set.insert(conv.tracks().front().castTo<reco::TrackRef>());
		 electrons_set.insert(conv.tracks().back( ).castTo<reco::TrackRef>());

		 float lead_pt = (conv.tracks().front()->pt() > conv.tracks().back()->pt()) ? conv.tracks().front()->pt() : conv.tracks().back()->pt();
		 if(lead_pt < minLeadPt_) continue;
		 passed[3]++;

		 try {
			 elecs->push_back(conv.tracks().front().castTo<reco::TrackRef>());
		 } catch(edm::Exception& e) {
			 //std::cout << "Invalid reference found" << std::endl;
		 }
		 try {
			 elecs->push_back(conv.tracks().back( ).castTo<reco::TrackRef>());
		 } catch(edm::Exception& e) {
			 //std::cout << "Invalid reference found" << std::endl;
		 }			 
	 }
	 //std::cout << "Passed: " << passed[0] << ", " << passed[1] << ", " << passed[2] << ", " << passed[3] << std::endl;

	 for(size_t itrk = 0; itrk < tracks->size(); itrk++) {
		 reco::TrackRef tkref(tracks, itrk);
		 if(electrons_set.find(tkref) == electrons_set.end()) {
			 not_elecs->push_back(tkref);
		 }
	 }

	 iEvent.put(std::move(elecs), "electrons");
	 iEvent.put(std::move(not_elecs), "NOTelectrons"); 
}

// ------------ method called once each stream before processing any runs, lumis or events  ------------
void
TracksFromConversions::beginStream(edm::StreamID)
{
}

// ------------ method called once each stream after processing all runs, lumis and events  ------------
void
TracksFromConversions::endStream() {
}

// ------------ method called when starting to processes a run  ------------
/*
void
TracksFromConversions::beginRun(edm::Run const&, edm::EventSetup const&)
{
}
*/
 
// ------------ method called when ending the processing of a run  ------------
/*
void
TracksFromConversions::endRun(edm::Run const&, edm::EventSetup const&)
{
}
*/
 
// ------------ method called when starting to processes a luminosity block  ------------
/*
void
TracksFromConversions::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
*/
 
// ------------ method called when ending the processing of a luminosity block  ------------
/*
void
TracksFromConversions::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
*/
 
// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
TracksFromConversions::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(TracksFromConversions);
