#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/makeRefToBaseProdFrom.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Common/interface/TriggerNames.h"
#include "FWCore/ServiceRegistry/interface/Service.h"

#include "DataFormats/Common/interface/TriggerResults.h"
#include "DataFormats/Common/interface/ValueMap.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectron.h"
#include "DataFormats/HLTReco/interface/TriggerEvent.h"

#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "MagneticField/Engine/interface/MagneticField.h"
#include "MagneticField/Records/interface/IdealMagneticFieldRecord.h"
#include "MagneticField/ParametrizedEngine/src/OAEParametrizedMagneticField.h"

#include "RecoVertex/KalmanVertexFit/interface/KalmanVertexFitter.h"

#include "TrackingTools/TransientTrack/interface/TransientTrack.h"
#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"

#include <vector>
#include <array>
#include <memory>

#include "TTree.h"
#include "TH1F.h"

using namespace std;

class ConversionsNtuples: public edm::EDAnalyzer {
public:
	explicit ConversionsNtuples(const edm::ParameterSet& cfg);

private:
  virtual void beginJob() {};
  virtual void endJob() {};
	virtual void analyze(const edm::Event&, const edm::EventSetup&);

  virtual void beginRun(edm::Run const&, edm::EventSetup const&){};
  virtual void endRun(edm::Run const&, edm::EventSetup const&){};
  virtual void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&){};
  virtual void endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&){};

protected:
	const edm::EDGetTokenT< std::vector<reco::Track> > tracks_tok_;
	const edm::EDGetTokenT< std::vector<reco::GsfElectron> > electrons_tok_;
	const edm::EDGetTokenT< std::vector<reco::GsfTrack> > electron_tracks_tok_;
	const edm::EDGetTokenT< edm::ValueMap<bool> > eid_tok_;
	const edm::EDGetTokenT< edm::TriggerResults > trig_tok_;
	const edm::EDGetTokenT< trigger::TriggerEvent > trig_evt_tok_;
	//const edm::EDGetTokenT< > tracks_tok_;

	//vertex fitter
	KalmanVertexFitter fitter_;
	//unique_ptr<OAEParametrizedMagneticField> field_;

	//selected trigger bits
	int ele27_idx_ = -1;
	int ele32_idx_ = -1;

	//cuts
	float trk_min_pt_ = 0;
	float trk_max_pt_ = 999;
	float trk_max_dtheta_ = 10;
	float vtx_max_chi2_ = 10;
	float vtx_max_mass_ = 10;

	//output
	TTree *tree_ = 0;
	
	//branches
	unsigned int b_lumi_ = 0;
  unsigned int b_run_ = 0;
  unsigned long long b_evt_ = 0;

	unsigned int b_nele_ = 0.;
	float b_ele_pt_ = -1.;
	float b_ele_eta_ = -10.;

	float b_vtx_chi2_ = -1.;
	float b_vtx_mass_ = -10.;
	float b_vtx_dtheta_ = -10.;

	float b_trk_pt_ = -1.;
	float b_trk_eta_ = -10.;
}; 

ConversionsNtuples::ConversionsNtuples(const edm::ParameterSet& cfg): 
	tracks_tok_{consumes< std::vector<reco::Track> >(cfg.getParameter<edm::InputTag>("tracks"))},
	electrons_tok_{consumes< vector<reco::GsfElectron> >(cfg.getParameter<edm::InputTag>("electrons"))},
	electron_tracks_tok_(consumes< vector<reco::GsfTrack> >(cfg.getParameter<edm::InputTag>("electronTracks"))),
	eid_tok_(consumes< edm::ValueMap<bool> >(cfg.getParameter<edm::InputTag>("eid"))),
	trig_tok_(consumes< edm::TriggerResults >(cfg.getParameter<edm::InputTag>("trigger"))),
	trig_evt_tok_(consumes< trigger::TriggerEvent >(cfg.getParameter<edm::InputTag>("triggerEvent"))),
	fitter_{},
	//field_(new OAEParametrizedMagneticField("3_8T")),
	trk_min_pt_(cfg.getParameter<double>("trkMinPt")),
	trk_max_pt_(cfg.getParameter<double>("trkMaxPt")),
	trk_max_dtheta_(cfg.getParameter<double>("trkMaxDTheta")),
	vtx_max_chi2_(cfg.getParameter<double>("vtxMaxChi2")),
	vtx_max_mass_(cfg.getParameter<double>("vtxMaxMass")) {

	edm::Service<TFileService> fs;
	cout << "fs " << fs << " available " << fs.isAvailable() << endl;
	
	tree_ = fs->make<TTree>("electrons", "test");

	//book branches
	tree_->Branch("run",  &b_run_, "run/i");
  tree_->Branch("lumi", &b_lumi_, "lumi/i");
  tree_->Branch("evt",  &b_evt_, "evt/i");

	tree_->Branch("nele"   , &b_nele_   , "nele/f");
	tree_->Branch("ele_pt" , &b_ele_pt_ , "ele_pt/f");
	tree_->Branch("ele_eta", &b_ele_eta_, "ele_eta/f");

	tree_->Branch("vtx_chi2",	  &b_vtx_chi2_	, "vtx_chi2/f"	);
	tree_->Branch("vtx_mass",	  &b_vtx_mass_	, "vtx_mass/f"	);
	tree_->Branch("vtx_dtheta", &b_vtx_dtheta_, "vtx_dtheta/f");

	tree_->Branch("trk_pt",		  &b_trk_pt_		, "trk_pt/f"		);
	tree_->Branch("trk_eta",    &b_trk_eta_   , "trk_eta/f"   );
}

void ConversionsNtuples::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup){
	edm::ESHandle<TransientTrackBuilder> ttrack_builder;
	iSetup.get<TransientTrackRecord>().get(
		"TransientTrackBuilder", ttrack_builder);

	typedef edm::Ref< std::vector<reco::GsfElectron> > ElectronRef;
	typedef edm::Ref< std::vector<reco::Track> > TrackRef;

	edm::Handle< std::vector<reco::Track> > tracks;
	iEvent.getByToken(tracks_tok_, tracks);

	edm::Handle< std::vector<reco::GsfElectron> > electrons;
	iEvent.getByToken(electrons_tok_, electrons);

	edm::Handle< std::vector<reco::GsfTrack> > electron_tracks;
	iEvent.getByToken(electron_tracks_tok_, electron_tracks);

	edm::Handle< edm::ValueMap<bool> > eid;
	iEvent.getByToken(eid_tok_, eid);

	edm::Handle< edm::TriggerResults > trigger;
	iEvent.getByToken(trig_tok_, trigger);

	edm::Handle< trigger::TriggerEvent > trig_event;
	iEvent.getByToken(trig_evt_tok_, trig_event);
	
	//find trigger index
	if(iEvent.id().run() != b_run_) {
		b_run_ = iEvent.id().run();
		ele32_idx_ = -1;
		ele27_idx_ = -1;
		const edm::TriggerNames& names = iEvent.triggerNames(*trigger);
		for(size_t idx = 0 ; idx < trigger->size() ; ++idx) {
			if(names.triggerName(idx).find("HLT_Ele32_WPTight_Gsf_v") != string::npos) {
					ele32_idx_ = idx;
			} else if(names.triggerName(idx).find("HLT_Ele27_WPTight_Gsf_v") != string::npos) {
					ele27_idx_ = idx;
			}			
		}
	}

	if(ele32_idx_ == -1 && ele27_idx_ == -1) {
		std::cerr << "NO TRIGGER FOUND!" << std::endl;
		throw 42;
	}

	bool ele32 = (ele32_idx_ > -1) && (trigger->accept(ele32_idx_));
	bool ele27 = (ele27_idx_ > -1) && (trigger->accept(ele27_idx_));
	
	if(!(ele32 || ele27)) return;
	float ele_thr = (ele27) ? 30. : 35.;

  b_lumi_ = iEvent.id().luminosityBlock();
  b_evt_  = iEvent.id().event();

	std::vector<ElectronRef> selected_electrons;
	for(size_t eidx=0; eidx < electrons->size(); ++eidx) {
		auto& ele = electrons->at(eidx);
		edm::Ref< std::vector<reco::GsfElectron> > ele_ref(electrons, eidx);
		if(ele.pt() > ele_thr && (*eid)[ele_ref]) selected_electrons.push_back(ele_ref);
	}
	b_nele_ = selected_electrons.size();
	if(b_nele_ == 0) return;

	//filter tracks
	std::vector<TrackRef> selected_tracks;
	for(size_t idx=0; idx < tracks->size(); ++idx) {
		if(tracks->at(idx).pt() > trk_min_pt_ &&
			 tracks->at(idx).pt() < trk_max_pt_
			) selected_tracks.emplace_back(tracks, idx);
	}
	if(!selected_tracks.size()) return;

	for(auto& ele : selected_electrons) {
		b_ele_pt_ = ele->pt();
		b_ele_eta_= ele->pt();

		b_vtx_chi2_	  = 999;
		b_vtx_mass_	  = 999;
		b_vtx_dtheta_ = 999;
		for(auto& trk : selected_tracks) {
			if((trk->charge()+ele->charge()) != 0) continue;
			reco::Candidate::LorentzVector trk_p4{trk->px(), trk->py(), trk->pz(), trk->p()}; //assume massless (good enough)
			float mass = (ele->p4() + trk_p4).mass();
			if(mass > vtx_max_mass_) continue;
			float dtheta = fabs(trk->theta() - ele->theta());
			if(dtheta < trk_max_dtheta_) continue;
			std::vector<reco::TransientTrack> trks{
				ttrack_builder->build(trk), ttrack_builder->build(ele->gsfTrack())};

			TransientVertex vtx = fitter_.vertex(trks);
			if(vtx.isValid() && vtx.normalisedChiSquared() < vtx_max_chi2_ &&
				 vtx.normalisedChiSquared() < b_vtx_chi2_) {
				b_vtx_chi2_ = vtx.normalisedChiSquared();
				b_vtx_mass_ = mass;
				b_vtx_dtheta_ = dtheta;
			}
		}

		if(b_vtx_chi2_ >= vtx_max_chi2_) continue;
		
		tree_->Fill();
	}
}


DEFINE_FWK_MODULE(ConversionsNtuples);

