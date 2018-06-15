/*
	Ugly copy-paste of RecoParticleFlow/PFTracking/plugins/GoodSeedProducer.cc
	to dump ALL the values used by the seeding to create a new training
	the code should be refactored to extract all the configurable parts
	(matching, extraction of features, etc.)
 */

// system include files
#include <memory>

// user include files

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "DataFormats/ParticleFlowReco/interface/PFCluster.h"
#include "DataFormats/ParticleFlowReco/interface/PFClusterFwd.h"
#include "DataFormats/Math/interface/LorentzVector.h"
#include "DataFormats/Math/interface/deltaPhi.h"
#include "TrackingTools/PatternTools/interface/Trajectory.h"
#include "RecoParticleFlow/PFTracking/interface/PFGeometry.h"

#include "RecoTracker/TransientTrackingRecHit/interface/TkTransientTrackingRecHitBuilder.h"

#include "RecoParticleFlow/PFTracking/interface/PFTrackTransformer.h"

#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "TrackingTools/TrackFitters/interface/TrajectoryFitter.h"
#include "TrackingTools/PatternTools/interface/TrajectorySmoother.h"
#include "TrackingTools/PatternTools/interface/TrajTrackAssociation.h"
#include "TrackingTools/Records/interface/TransientRecHitRecord.h"  
//#include "TrackingTools/TrajectoryState/interface/TrajectoryStateTransform.h"
#include "FastSimulation/BaseParticlePropagator/interface/BaseParticlePropagator.h"
#include "MagneticField/Engine/interface/MagneticField.h"
#include "MagneticField/Records/interface/IdealMagneticFieldRecord.h"

#include "DataFormats/Math/interface/deltaR.h"
#include "TMVA/MethodBDT.h"
#include "TMVA/Reader.h"
#include "CondFormats/EgammaObjects/interface/GBRForest.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"

#include <fstream>
#include <string>
#include "TMath.h"
#include "TFile.h"
#include "TTree.h"
#include "TH1D.h"
#include "TRandom3.h"

using namespace edm;
using namespace std;
using namespace reco;
constexpr size_t kMaxWeights=9;

class TrackerDrivenSeedFeatures: public edm::EDAnalyzer {
  typedef TrajectoryStateOnSurface TSOS;
public:
  explicit TrackerDrivenSeedFeatures(const edm::ParameterSet&);
	~TrackerDrivenSeedFeatures() {
		file_->cd();
		tree_->Write();
		pt_hist_->Write();
		file_->Write();
		file_->Close();
	}

	struct Features { 
		unsigned int lumi = 0;
		unsigned int run = 0;
		unsigned long long evt = 0;

		float nhits = -1;
		int ibin = -1;
		int high_purity = 0;
		bool trk_ecal_match = false;
		float dxy = -1;
		float dxy_err = -1;
		float bdtout = -1.;
		float trk_pt = -1.;
		float trk_inp = -1.;
		float trk_outp = -1.;
		float trk_eta = -1.;
		float trk_ecal_Deta = -1.;
		float trk_ecal_Dphi = -1.;
		float trk_ecal_AbsDphi = -1.;
		float e_over_p = -1.;
		float trk_chi2red = -1.;

		//ctk residuals
		vector<float> trk_residuals = {};
		vector<float> trk_momenta = {};
		
		//stage 2, with GSF
		bool gsf_success = false;
		float gsf_dpt = -1.;
		float trk_gsf_chiratio = -1.;
		float gsf_chi2red = -1.;
		
		void reset() {
			lumi = 0;
			run = 0;
			evt = 0;

			nhits = -1;
			ibin = -1;
			high_purity = 0;
			trk_ecal_match = false;
			dxy = -1;
			dxy_err = -1;
			bdtout = -1.;
			trk_pt = -1.;
			trk_inp = -1.;
			trk_outp = -1.;
			trk_eta = -1.;
			trk_ecal_Deta = -1.;
			trk_ecal_AbsDphi = -1.;
			trk_ecal_Dphi = -1.;
			e_over_p = -1.;
			trk_chi2red = -1.;
			
			//ctk residuals
			trk_residuals.clear();
			trk_momenta.clear();

			//stage 2, with GSF
			gsf_success = false;
			gsf_dpt = -1.;
			trk_gsf_chiratio = -1.;
			gsf_chi2red = -1.;			
		}
	};

private:
	virtual void beginRun(const edm::Run & run,const edm::EventSetup&);
	virtual void analyze(const edm::Event&, const edm::EventSetup&);
		
	int getBin(float eta, float pt) {
		  int ie=0;
			int ip=0;
			if (fabs(eta)<0.8) ie=0;
			else{ if (fabs(eta)<1.479) ie=1;
				else ie=2;
			}
			if (pt<6) ip=0;
			else {  if (pt<12) ip=1;     
				else ip=2;
			}
			int iep= ie*3+ip;
			LogDebug("GoodSeedProducer")<<"Track pt ="<<pt<<" eta="<<eta<<" bin="<<iep;
			return iep;
	}
	// ----------member data ---------------------------

	Features trk_features_;
		
	///GSF Fitter & Smoother and ancillary stuff
	std::unique_ptr<TrajectoryFitter> fitter_;
	std::unique_ptr<TrajectorySmoother> smoother_;
	TkClonerImpl hitCloner_;
	
	///PFTrackTransformer (needed?)
	std::unique_ptr<PFTrackTransformer> pfTransformer_;
	
	///Minimum transverse momentum and maximum pseudorapidity
	double minPt_;
	double maxPt_;
      
	///Cut on the energy of the clusters
	double clusThreshold_;

	///Min and MAx allowed values forEoverP
	double minEp_;
	double maxEp_;
	double prescale_;

	// ----------access to event data
	const edm::EDGetTokenT<reco::PFClusterCollection> pfCLusTagECLabel_;
	const edm::EDGetTokenT<reco::TrackRefVector > tracks_;
	//const edm::EDGetTokenT< std::vector<Trajectory> > trajectories_;
	const edm::EDGetTokenT< reco::BeamSpot > beamspot_;

	std::string fitterName_;
	std::string smootherName_;
	std::string trackerRecHitBuilderName_;
      
	double Min_dr_;

	///B field
	math::XYZVector B_;

	//output
	TFile *file_ = 0;
	TTree *tree_ = 0;
	TH1D *pt_hist_;
	std::array<std::unique_ptr<const GBRForest>,kMaxWeights> gbrs_;    
};

TrackerDrivenSeedFeatures::TrackerDrivenSeedFeatures(const ParameterSet& iConfig):
	fitter_(nullptr),
	smoother_(nullptr),
	hitCloner_(),
  pfTransformer_(nullptr),
	minPt_{iConfig.getParameter<double>("MinPt")},
	maxPt_{iConfig.getParameter<double>("MaxPt")},
	clusThreshold_{iConfig.getParameter<double>("ClusterThreshold")},
  minEp_{iConfig.getParameter<double>("MinEOverP")},
  maxEp_{iConfig.getParameter<double>("MaxEOverP")},
	prescale_{1.-iConfig.getParameter<double>("prescale")},
  pfCLusTagECLabel_{consumes<reco::PFClusterCollection>(iConfig.getParameter<InputTag>("PFEcalClusterLabel"))},
	tracks_{consumes<reco::TrackRefVector>(iConfig.getParameter<edm::InputTag>("tracks"))},
	//trajectories_{consumes< std::vector<Trajectory> >(iConfig.getParameter<edm::InputTag>("tracks"))},	
	beamspot_{consumes<reco::BeamSpot>(iConfig.getParameter<edm::InputTag>("beamspot"))},
  fitterName_  {iConfig.getParameter<string>("Fitter")  },
  smootherName_{iConfig.getParameter<string>("Smoother")},
  trackerRecHitBuilderName_{iConfig.getParameter<std::string>("TTRHBuilder")},
	Min_dr_{iConfig.getParameter<double>("Min_dr")},
	B_()
{
  LogInfo("TrackerDrivenSeedFeatures")<<"Electron PreIdentification started  ";
    
	file_ = new TFile(iConfig.getParameter<string>("filename").c_str(), "RECREATE");
	file_->cd();
	pt_hist_ = new TH1D("pt", "pt", 100, 0, 100);
	tree_ = new TTree("tree", "test");

	//book branches
	//tree_->Branch("run",  &b_run_, "run/i");
	tree_->Branch("run",  &trk_features_.run, "run/i");
  tree_->Branch("lumi", &trk_features_.lumi, "lumi/i");
  tree_->Branch("evt",  &trk_features_.evt, "evt/i");

	tree_->Branch("nhits",					&trk_features_.nhits         , "nhits/f");
	tree_->Branch("ibin", 					&trk_features_.ibin         , "ibin/I");
	tree_->Branch("high_purity",		&trk_features_.high_purity	 , "high_purity/i");
	tree_->Branch("trk_ecal_match", &trk_features_.trk_ecal_match, "trk_ecal_match/O");
	tree_->Branch("bdtout",				 	&trk_features_.bdtout				 , "bdtout/f");
	tree_->Branch("dxy",			  	 	&trk_features_.dxy		  		 , "dxy/f");
	tree_->Branch("dxy_err",			 	&trk_features_.dxy_err			 , "dxy_err/f");
	tree_->Branch("trk_pt",				 	&trk_features_.trk_pt				 , "trk_pt/f");
	tree_->Branch("trk_inp",  			&trk_features_.trk_inp  		 , "trk_inp/f");
	tree_->Branch("trk_outp",	  		&trk_features_.trk_outp	  	 , "trk_outp/f");
	tree_->Branch("trk_eta",		 	  &trk_features_.trk_eta		 , "trk_eta/f");
	tree_->Branch("trk_ecal_Deta",	&trk_features_.trk_ecal_Deta , "trk_ecal_Deta/f");
	tree_->Branch("trk_ecal_Dphi",	&trk_features_.trk_ecal_Dphi , "trk_ecal_Dphi/f");
	tree_->Branch("e_over_p",			 	&trk_features_.e_over_p			 , "e_over_p/f");
	tree_->Branch("trk_chi2red",    &trk_features_.trk_chi2red   , "trk_chi2red/f"); 

	tree_->Branch("trk_momenta", &trk_features_.trk_momenta);
	/*tree_->Branch("trk_res0",	&trk_features_.trk_res0, "trk_res0/f");
	tree_->Branch("trk_res1",	&trk_features_.trk_res1, "trk_res1/f");
	tree_->Branch("trk_res2",	&trk_features_.trk_res2, "trk_res2/f");
	tree_->Branch("trk_res3",	&trk_features_.trk_res3, "trk_res3/f");
	tree_->Branch("trk_res4",	&trk_features_.trk_res4, "trk_res4/f");
	tree_->Branch("trk_res5",	&trk_features_.trk_res5, "trk_res5/f");*/

	tree_->Branch("gsf_success"			, &trk_features_.gsf_success     , "gsf_success/O");
	tree_->Branch("gsf_dpt"					, &trk_features_.gsf_dpt				 , "gsf_dpt/f");					
	tree_->Branch("trk_gsf_chiratio", &trk_features_.trk_gsf_chiratio, "trk_gsf_chiratio/f");
	tree_->Branch("gsf_chi2red"     , &trk_features_.gsf_chi2red     , "gsf_chi2red/f");     

	const std::string method = iConfig.getParameter<string>("TMVAMethod");
	std::array<edm::FileInPath,kMaxWeights> weights = {{ 
			edm::FileInPath(iConfig.getParameter<string>("Weights1")),
			edm::FileInPath(iConfig.getParameter<string>("Weights2")),
			edm::FileInPath(iConfig.getParameter<string>("Weights3")),
			edm::FileInPath(iConfig.getParameter<string>("Weights4")),
			edm::FileInPath(iConfig.getParameter<string>("Weights5")),
			edm::FileInPath(iConfig.getParameter<string>("Weights6")),
			edm::FileInPath(iConfig.getParameter<string>("Weights7")),
			edm::FileInPath(iConfig.getParameter<string>("Weights8")),
			edm::FileInPath(iConfig.getParameter<string>("Weights9")) }};
            
	for(UInt_t j = 0; j < gbrs_.size(); ++j){
		TMVA::Reader reader("!Color:Silent");
		
		reader.AddVariable("NHits", &trk_features_.nhits);
		reader.AddVariable("NormChi", &trk_features_.trk_chi2red);
		reader.AddVariable("dPtGSF", &trk_features_.gsf_dpt);
		reader.AddVariable("EoP", &trk_features_.e_over_p);
		reader.AddVariable("ChiRatio", &trk_features_.trk_gsf_chiratio);
		reader.AddVariable("RedChi", &trk_features_.gsf_chi2red);
		reader.AddVariable("EcalDEta", &trk_features_.trk_ecal_Deta);
		reader.AddVariable("EcalDPhi", &trk_features_.trk_ecal_AbsDphi);
		reader.AddVariable("pt", &trk_features_.trk_pt);
		reader.AddVariable("eta", &trk_features_.trk_eta);
    
		reader.BookMVA(method, weights[j].fullPath().c_str());
    
		gbrs_[j].reset( new GBRForest( dynamic_cast<TMVA::MethodBDT*>( reader.FindMVA(method) ) ) );
	}
}


//
// member functions
//

// ------------ method called to produce the data  ------------
void
TrackerDrivenSeedFeatures::analyze(const Event& iEvent, const EventSetup& iSetup)
{
  LogDebug("TrackerDrivenSeedFeatures")<<"START event: "<<iEvent.id().event()
															<<" in run "<<iEvent.id().run();

	edm::ESHandle<TrajectoryFitter> aFitter;
	edm::ESHandle<TrajectorySmoother> aSmoother;
	iSetup.get<TrajectoryFitter::Record>().get(fitterName_, aFitter);
	iSetup.get<TrajectoryFitter::Record>().get(smootherName_, aSmoother);
	smoother_.reset(aSmoother->clone());
	fitter_ = aFitter->clone();
	edm::ESHandle<TransientTrackingRecHitBuilder> theTrackerRecHitBuilder;
	iSetup.get<TransientRecHitRecord>().get(trackerRecHitBuilderName_,theTrackerRecHitBuilder);
	hitCloner_ = static_cast<TkTransientTrackingRecHitBuilder const *>(theTrackerRecHitBuilder.product())->cloner();
	fitter_->setHitCloner(&hitCloner_);
	smoother_->setHitCloner(&hitCloner_);

  //Magnetic Field
  ESHandle<MagneticField> magneticField;
  iSetup.get<IdealMagneticFieldRecord>().get(magneticField);

  //Handle input collections
  //ECAL clusters	      
  Handle<PFClusterCollection> theECPfClustCollection;
  iEvent.getByToken(pfCLusTagECLabel_,theECPfClustCollection);  

	edm::Handle<reco::BeamSpot> beam_spot;
	iEvent.getByToken(beamspot_, beam_spot);

	vector<PFCluster const *> basClus;
  for ( auto const & klus : *theECPfClustCollection.product() ) {
    if(klus.correctedEnergy()>clusThreshold_) basClus.push_back(&klus);
  }

  //Vector of track collections
	size_t ntrks = 0;
	//Track collection
	Handle<reco::TrackRefVector> tkRefCollection;
	iEvent.getByToken(tracks_, tkRefCollection);
	//const TrackRefVector&  Tk=*(tkRefCollection.product());

	// edm::Handle< std::vector<Trajectory> > trajectories;  
	// iEvent.getByToken(trajectories_, trajectories);	
	// assert(tkRefCollection->size() == trajectories->size());

	//loop over the track collection
	//cout << "Starting from " << tkRefCollection->size() << " tracks! " << endl;
	for(const auto& trackRef : *tkRefCollection) {
		trk_features_.reset();
      
		trk_features_.run  = iEvent.id().run();
		trk_features_.lumi = iEvent.id().luminosityBlock();
		trk_features_.evt  = iEvent.id().event();
		//const TrackRef& trackRef = tkRefCollection.at(i);
		//Trajectory const * trajectory = &(*trajectories)[i];

		math::XYZVectorF tkmom(trackRef->momentum());
		auto tketa= trackRef->eta();
		auto tkpt = trackRef->pt();
		auto const & Seed=(*trackRef->seedRef());

		pt_hist_->Fill((tkpt > 100) ? 99.9 : tkpt);
		if (tkpt>maxPt_ || tkpt<minPt_) continue;			
		if(gRandom->Rndm() < prescale_) continue;
		ntrks++;
		//cout << "pt: " << tkpt <<endl;
		//cout << "P: " << tkmom.mag2() << ", p(in): " << trackRef->innerMomentum().mag2() << ", p(out): " << trackRef->outerMomentum().mag2() << endl;

		// const auto& measurements = trajectory->measurements();
		// for(auto& measurement : measurements) {
		// 	trk_features_.trk_momenta.push_back(
		// 		trackRef->
		// 		);
		// }

		float oPTOB=1.f/std::sqrt(trackRef->innerMomentum().mag2()); // FIXME the original code was buggy should be outerMomentum...
		float nchi=trackRef->normalizedChi2();

		int nhitpi=trackRef->found();
		float EP=0;
      
		//CLUSTERS - TRACK matching
      
		auto pfmass=  0.0005;
		auto pfoutenergy=sqrt((pfmass*pfmass)+trackRef->outerMomentum().Mag2());

		XYZTLorentzVector mom =XYZTLorentzVector(
			trackRef->outerMomentum().x(),
			trackRef->outerMomentum().y(),
			trackRef->outerMomentum().z(),
			pfoutenergy);
		XYZTLorentzVector pos =   XYZTLorentzVector(
			trackRef->outerPosition().x(),
			trackRef->outerPosition().y(),
			trackRef->outerPosition().z(),
			0.);

		BaseParticlePropagator theOutParticle( RawParticle(mom,pos),
																					 0,0,B_.z());
		theOutParticle.setCharge(trackRef->charge());
      
		theOutParticle.propagateToEcalEntrance(false);
      
//const std::vector<TrajectoryMeasurement> &tmColl = trajectory->measurements();
      
		float toteta=10.f;
		float totphi=6.f;
		float dr=1000.f;
		GlobalPoint ElecTrkEcalPos(0,0,0);

		PFClusterRef clusterRef;
		math::XYZPoint meanShowerSaved;
		if(theOutParticle.getSuccess()!=0){
			ElecTrkEcalPos=GlobalPoint(theOutParticle.vertex().x(),
																 theOutParticle.vertex().y(),
																 theOutParticle.vertex().z()
				);

			constexpr float psLim = 2.50746495928f; // std::sinh(1.65f);
			bool isBelowPS= (ElecTrkEcalPos.z()*ElecTrkEcalPos.z()) > (psLim*psLim)*ElecTrkEcalPos.perp2();
			// bool isBelowPS=(std::abs(ElecTrkEcalPos.eta())>1.65f);	
	
			unsigned clusCounter=0;
			float max_ee = 0;
			for(auto aClus : basClus) {

				float tmp_ep=float(aClus->correctedEnergy())*oPTOB;
				if ((tmp_ep<minEp_)|(tmp_ep>maxEp_)) { ++clusCounter; continue;}
	    
				double ecalShowerDepth
					= PFCluster::getDepthCorrection(aClus->correctedEnergy(),
																					isBelowPS,
																					false);
				auto mom = theOutParticle.momentum().Vect();
				auto meanShower = ElecTrkEcalPos +
					GlobalVector(mom.x(),mom.y(),mom.z()).unit()*ecalShowerDepth;	
	  
				float etarec=meanShower.eta();
				float phirec=meanShower.phi();
	     

				float tmp_phi = deltaPhi(aClus->positionREP().phi(), phirec); 
				tmp_phi *= trackRef->charge();
	      
				float tmp_dr=std::sqrt(std::pow(tmp_phi,2.f)+
															 std::pow(aClus->positionREP().eta()-etarec,2.f));
	  
				if (tmp_dr<dr){
					dr=tmp_dr;
					if(dr < Min_dr_){ // find the most closest and energetic ECAL cluster
						if(aClus->correctedEnergy() > max_ee){

							toteta=aClus->positionREP().eta()-etarec;
							totphi=tmp_phi;
							EP=tmp_ep;
							clusterRef = PFClusterRef(theECPfClustCollection,clusCounter);
							meanShowerSaved = meanShower;
		    
						}
					}
				}
				++clusCounter;
			}
		}
		float trk_ecalDeta = fabs(toteta);
		float trk_ecalDphi = totphi;
			
		trk_features_.nhits = nhitpi;
		trk_features_.ibin = getBin(trackRef->eta(),trackRef->pt());
		trk_features_.high_purity = trackRef->quality(
			TrackBase::qualityByName("highPurity"));
		trk_features_.dxy		  = trackRef->dxy(*beam_spot);
		trk_features_.dxy_err	= trackRef->dxyError();
		trk_features_.trk_ecal_match = (toteta < 10.f);
		trk_features_.trk_pt = tkpt;
		trk_features_.trk_outp = sqrt(trackRef->outerMomentum().mag2());
		trk_features_.trk_inp = sqrt(trackRef->innerMomentum().mag2());
		trk_features_.trk_eta = tketa;
		trk_features_.trk_ecal_Deta = trk_ecalDeta;
		trk_features_.trk_ecal_Dphi = trk_ecalDphi;
		trk_features_.trk_ecal_AbsDphi = fabs(trk_ecalDphi);
		trk_features_.e_over_p = EP;
		trk_features_.trk_chi2red = nchi;

		float chired=1000;
		float chiRatio=1000;
		float dpt=0;

		Trajectory::ConstRecHitContainer tmp;
		auto hb = trackRef->recHitsBegin();
		for(unsigned int h=0;h<trackRef->recHitsSize();h++){
			auto recHit = *(hb+h); tmp.push_back(recHit->cloneSH());
		}

		GlobalVector gv(trackRef->innerMomentum().x(),trackRef->innerMomentum().y(),trackRef->innerMomentum().z());
		GlobalPoint  gp(trackRef->innerPosition().x(),trackRef->innerPosition().y(),trackRef->innerPosition().z());
		GlobalTrajectoryParameters gtps(gp,gv,trackRef->charge(),&*magneticField);
		TrajectoryStateOnSurface tsos(gtps,trackRef->innerStateCovariance(),*tmp[0]->surface());
		Trajectory  && FitTjs= fitter_->fitOne(Seed,tmp,tsos);
	
		if(FitTjs.isValid()){
			Trajectory && SmooTjs= smoother_->trajectory(FitTjs);
			if(SmooTjs.isValid()){
					
				//Track refitted with electron hypothesis
					
				float pt_out=SmooTjs.firstMeasurement().
					updatedState().globalMomentum().perp();
				float pt_in=SmooTjs.lastMeasurement().
					updatedState().globalMomentum().perp();
				dpt=(pt_in>0) ? fabs(pt_out-pt_in)/pt_in : 0.;
				// the following is simply the number of degrees of freedom
				chiRatio=SmooTjs.chiSquared()/trackRef->chi2();
				chired=chiRatio*nchi;

				trk_features_.gsf_success = true;
				trk_features_.gsf_dpt = dpt;
				trk_features_.trk_gsf_chiratio = chiRatio;
				trk_features_.gsf_chi2red = chired;
				
				float vars[10] = { 
					trk_features_.nhits,
					trk_features_.trk_chi2red,
					trk_features_.gsf_dpt,
					trk_features_.e_over_p,
					trk_features_.trk_gsf_chiratio,
					trk_features_.gsf_chi2red,
					trk_features_.trk_ecal_Deta,
					trk_features_.trk_ecal_AbsDphi,
					trk_features_.trk_pt,
					trk_features_.trk_eta
				};
              
				trk_features_.bdtout = gbrs_[trk_features_.ibin]->GetClassifier( 
					vars
					);
			}
		}
			  
		tree_->Fill();
	} //end loop on track collection
	//cout << "processed " << ntrks << " tracks! " << endl;
  
}


// ------------ method called once each job just before starting event loop  ------------
void 
TrackerDrivenSeedFeatures::beginRun(const edm::Run & run,
													 const EventSetup& es)
{
  //Magnetic Field
  ESHandle<MagneticField> magneticField;
  es.get<IdealMagneticFieldRecord>().get(magneticField);
  B_=magneticField->inTesla(GlobalPoint(0,0,0));
  
  pfTransformer_.reset( new PFTrackTransformer(B_) );
  pfTransformer_->OnlyProp();
}

DEFINE_FWK_MODULE(TrackerDrivenSeedFeatures);
