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

#include <fstream>
#include <string>
#include "TMath.h"
#include "TFile.h"
#include "TTree.h"
#include "TH1D.h"

using namespace edm;
using namespace std;
using namespace reco;

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
		int nhits = -1;
		int high_purity = 0;
		bool trk_ecal_match = false;
		float trk_pt = -1.;
		float trk_inp = -1.;
		float trk_outp = -1.;
		float trk_abseta = -1.;
		float trk_ecal_Deta = -1.;
		float trk_ecal_Dphi = -1.;
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
			nhits = -1;
			high_purity = 0;
			trk_ecal_match = false;
			trk_pt = -1.;
			trk_inp = -1.;
			trk_outp = -1.;
			trk_abseta = -1.;
			trk_ecal_Deta = -1.;
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

	// ----------access to event data
	edm::EDGetTokenT<reco::PFClusterCollection> pfCLusTagECLabel_;
	edm::EDGetTokenT<reco::TrackCollection > tracks_;
	edm::EDGetTokenT< std::vector<Trajectory> > trajectories_;

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
  pfCLusTagECLabel_{consumes<reco::PFClusterCollection>(iConfig.getParameter<InputTag>("PFEcalClusterLabel"))},
	tracks_{consumes<reco::TrackCollection>(iConfig.getParameter<edm::InputTag>("tracks"))},
	trajectories_{consumes< std::vector<Trajectory> >(iConfig.getParameter<edm::InputTag>("tracks"))},	
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

	tree_->Branch("nhits",					&trk_features_.nhits         , "nhits/I");
	tree_->Branch("high_purity",		&trk_features_.high_purity	 , "high_purity/i");
	tree_->Branch("trk_ecal_match", &trk_features_.trk_ecal_match, "trk_ecal_match/O");
	tree_->Branch("trk_pt",				 	&trk_features_.trk_pt				 , "trk_pt/f");
	tree_->Branch("trk_inp",  			&trk_features_.trk_inp  		 , "trk_inp/f");
	tree_->Branch("trk_outp",	  		&trk_features_.trk_outp	  	 , "trk_outp/f");
	tree_->Branch("trk_abseta",		 	&trk_features_.trk_abseta		 , "trk_abseta/f");
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
}


//
// member functions
//

// ------------ method called to produce the data  ------------
void
TrackerDrivenSeedFeatures::analyze(const Event& iEvent, const EventSetup& iSetup)
{
	trk_features_.reset();
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

	vector<PFCluster const *> basClus;
  for ( auto const & klus : *theECPfClustCollection.product() ) {
    if(klus.correctedEnergy()>clusThreshold_) basClus.push_back(&klus);
  }

  //Vector of track collections
	size_t ntrks = 0;
	//Track collection
	Handle<TrackCollection> tkRefCollection;
	iEvent.getByToken(tracks_, tkRefCollection);
	const TrackCollection&  Tk=*(tkRefCollection.product());

	// edm::Handle< std::vector<Trajectory> > trajectories;  
	// iEvent.getByToken(trajectories_, trajectories);	
	// assert(tkRefCollection->size() == trajectories->size());

	//loop over the track collection
	cout << "Starting from " << Tk.size() << " tracks! " << endl;
	for(unsigned int i=0;i<Tk.size();++i){		
      
		TrackRef trackRef(tkRefCollection, i);
		//Trajectory const * trajectory = &(*trajectories)[i];

		math::XYZVectorF tkmom(Tk[i].momentum());
		auto tketa= Tk[i].eta();
		auto tkpt = Tk[i].pt();
		auto const & Seed=(*trackRef->seedRef());

		pt_hist_->Fill((tkpt > 100) ? 99.9 : tkpt);
		if (tkpt>maxPt_ || tkpt<minPt_) continue;			
		ntrks++;
		//cout << "pt: " << tkpt <<endl;
		//cout << "P: " << tkmom.mag2() << ", p(in): " << Tk[i].innerMomentum().mag2() << ", p(out): " << Tk[i].outerMomentum().mag2() << endl;
v
		// const auto& measurements = trajectory->measurements();
		// for(auto& measurement : measurements) {
		// 	trk_features_.trk_momenta.push_back(
		// 		trackRef->
		// 		);
		// }

		float oPTOB=1.f/std::sqrt(Tk[i].innerMomentum().mag2()); // FIXME the original code was buggy should be outerMomentum...
		float nchi=Tk[i].normalizedChi2();

		int nhitpi=Tk[i].found();
		float EP=0;
      
		//CLUSTERS - TRACK matching
      
		auto pfmass=  0.0005;
		auto pfoutenergy=sqrt((pfmass*pfmass)+Tk[i].outerMomentum().Mag2());

		XYZTLorentzVector mom =XYZTLorentzVector(
			Tk[i].outerMomentum().x(),
			Tk[i].outerMomentum().y(),
			Tk[i].outerMomentum().z(),
			pfoutenergy);
		XYZTLorentzVector pos =   XYZTLorentzVector(
			Tk[i].outerPosition().x(),
			Tk[i].outerPosition().y(),
			Tk[i].outerPosition().z(),
			0.);

		BaseParticlePropagator theOutParticle( RawParticle(mom,pos),
																					 0,0,B_.z());
		theOutParticle.setCharge(Tk[i].charge());
      
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
				tmp_phi *= Tk[i].charge();
	      
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
		trk_features_.high_purity = Tk[i].quality(
			TrackBase::qualityByName("highPurity"));
		trk_features_.trk_ecal_match = (toteta < 10.f);
		trk_features_.trk_pt = tkpt;
		trk_features_.trk_outp = sqrt(Tk[i].outerMomentum().mag2());
		trk_features_.trk_inp = sqrt(Tk[i].innerMomentum().mag2());
		trk_features_.trk_abseta = fabs(tketa);
		trk_features_.trk_ecal_Deta = trk_ecalDeta;
		trk_features_.trk_ecal_Dphi = trk_ecalDphi;
		trk_features_.e_over_p = EP;
		trk_features_.trk_chi2red = nchi;

		float chired=1000;
		float chiRatio=1000;
		float dpt=0;

		Trajectory::ConstRecHitContainer tmp;
		auto hb = Tk[i].recHitsBegin();
		for(unsigned int h=0;h<Tk[i].recHitsSize();h++){
			auto recHit = *(hb+h); tmp.push_back(recHit->cloneSH());
		}
		auto const & theTrack = Tk[i]; 
		GlobalVector gv(theTrack.innerMomentum().x(),theTrack.innerMomentum().y(),theTrack.innerMomentum().z());
		GlobalPoint  gp(theTrack.innerPosition().x(),theTrack.innerPosition().y(),theTrack.innerPosition().z());
		GlobalTrajectoryParameters gtps(gp,gv,theTrack.charge(),&*magneticField);
		TrajectoryStateOnSurface tsos(gtps,theTrack.innerStateCovariance(),*tmp[0]->surface());
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
				chiRatio=SmooTjs.chiSquared()/Tk[i].chi2();
				chired=chiRatio*nchi;

				trk_features_.gsf_success = true;
				trk_features_.gsf_dpt = dpt;
				trk_features_.trk_gsf_chiratio = chiRatio;
				trk_features_.gsf_chi2red = chired;
			}
		}
			  
		tree_->Fill();
	} //end loop on track collection
	cout << "processed " << ntrks << " tracks! " << endl;
  
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