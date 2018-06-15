// -*- C++ -*-
//
// Package:    CMGTools/TracksFromGenParticles
// Class:      TracksFromGenParticles
// 
/**\class TracksFromGenParticles TracksFromGenParticles.cc CMGTools/TracksFromGenParticles/plugins/TracksFromGenParticles.cc

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
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/RecoCandidate/interface/TrackAssociation.h"
//
// class declaration
//

class TracksFromGenParticles : public edm::stream::EDProducer<> {
   public:
      explicit TracksFromGenParticles(const edm::ParameterSet&);
      ~TracksFromGenParticles();

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
	const edm::EDGetTokenT< edm::View<reco::Track> > tracks_;
	const edm::EDGetTokenT< reco::RecoToSimCollection > association_;
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
TracksFromGenParticles::TracksFromGenParticles(const edm::ParameterSet& cfg):
	tracks_{consumes< edm::View<reco::Track> >(cfg.getParameter<edm::InputTag>("tracks"))},
	association_{consumes< reco::RecoToSimCollection >(cfg.getParameter<edm::InputTag>("association"))}
{
   //register your products
	produces<reco::TrackRefVector>();
}


TracksFromGenParticles::~TracksFromGenParticles()
{
 
   // do anything here that needs to be done at destruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called to produce the data  ------------
void
TracksFromGenParticles::produce(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
   using namespace edm;
	 std::unique_ptr<reco::TrackRefVector> elecs(    new reco::TrackRefVector()); 

	 edm::Handle< edm::View<reco::Track> > tracks;
	 iEvent.getByToken(tracks_, tracks);

	 edm::Handle<reco::RecoToSimCollection> matching;
	 iEvent.getByToken(association_, matching);

	 for(size_t idx=0; idx<tracks->size(); ++idx) {
		 RefToBase<reco::Track> key(tracks, idx);
		 auto match = matching->find(key);
		 if(match != matching->end()) { 
			 auto tracking_particle = match->val.front().first;
			 if(std::abs(tracking_particle->pdgId()) == 11) {
				 elecs->push_back(key.castTo<reco::TrackRef>());
			 }
		 }
	 }

	 iEvent.put(std::move(elecs));
}

// ------------ method called once each stream before processing any runs, lumis or events  ------------
void
TracksFromGenParticles::beginStream(edm::StreamID)
{
}

// ------------ method called once each stream after processing all runs, lumis and events  ------------
void
TracksFromGenParticles::endStream() {
}

// ------------ method called when starting to processes a run  ------------
/*
void
TracksFromGenParticles::beginRun(edm::Run const&, edm::EventSetup const&)
{
}
*/
 
// ------------ method called when ending the processing of a run  ------------
/*
void
TracksFromGenParticles::endRun(edm::Run const&, edm::EventSetup const&)
{
}
*/
 
// ------------ method called when starting to processes a luminosity block  ------------
/*
void
TracksFromGenParticles::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
*/
 
// ------------ method called when ending the processing of a luminosity block  ------------
/*
void
TracksFromGenParticles::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
*/
 
// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
TracksFromGenParticles::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(TracksFromGenParticles);
