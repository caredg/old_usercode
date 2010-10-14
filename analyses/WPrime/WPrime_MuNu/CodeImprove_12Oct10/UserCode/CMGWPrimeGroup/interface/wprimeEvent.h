#ifndef _wprimeEvent_h_
#define _wprimeEvent_h_

// root stuff
#include <TROOT.h>
#include <TObject.h>
#include <TLorentzVector.h>
#include <TBranch.h>
#include <TClonesArray.h>

// std library
#include <string>

namespace wprime{

/// class containing per-job (or file) Wprime info
class JobInfo : public TObject {
 public: 
   /// RECO version
   std::string RECOversion;
   /// sample type (e.g. "QCD50_80" or "Spring10 collisions")
   std::string sample;
	  
   JobInfo();
   ~JobInfo();
	  
    ClassDef(JobInfo, 1)

};

/// class containing per-run Wprime info
class RunInfo : public TObject {
 public: 
   /// HLT version
   std::string HLTversion;
   /// HLT menu (ConfDB) name
   std::string HLTmenu;
   /// run number
   Int_t run_no;
   /// # of events processed
   Int_t Nproc_evt;

   RunInfo();
   ~RunInfo();
	  
    ClassDef(RunInfo, 1)
      
};

/// class containing info for MC-truth particles
class MCParticle : public TObject {
 public: 
  /// kinematic info
   TLorentzVector p;
   /// charge
   Int_t q;
   /// pdgId of mother particle
   Int_t momId;
   /// status of particle
   Int_t status;
   /// 
   MCParticle();
   MCParticle(TLorentzVector& p_,Int_t q_,Int_t momId_, Int_t status_);
   ~MCParticle();

   ClassDef(MCParticle, 1)
};


/// class with track information
class Track : public TObject {
 public: 
  /// kinematic info
   TLorentzVector p;
   /// charge
   Int_t q;
   /// chi^2
   Float_t chi2;
   /// d0 w.r.t. (0,0)
   Float_t d0_default;
   /// uncertainty on d0 
   Float_t dd0_default;
   /// d0 w.r.t. the beamspot
   Float_t d0;
   /// uncertainty on pt
   Float_t dpt;
   /// uncertainty on q/p
   Float_t dq_over_p;
   /// Ndof
   Int_t ndof;
   /// # of siStrip layers used by fit
   Int_t Nstrip_layer;
   /// # of pixel layers used by fit
   Int_t Npixel_layer;
   /// # of siStrip layers without measurement
   Int_t Nstrip_layerNoMeas;
   /// # of pixel layers without measurement
   Int_t Npixel_layerNoMeas;
   /// # of siStrip hits used by fit
   Int_t NsiStrip_hits;
   /// # of pixel hits used by fit
   Int_t Npixel_hits;
   /// # of muon hits used by fit
   Int_t Nmuon_hits;
   /// # of tracker-only hits used by fit
   Int_t Ntrk_hits;
   /// # of total (ie. tracker + muon) hits used by fit
   Int_t Ntot_hits;

   Track();
   ~Track();

   ClassDef(Track, 1)
};

const int N_CONESIZE = 9;
const float MIN_CONE = 0.20;
const float MAX_CONE = 0.60;
const float MUON_MASS = 0.105658369; // GeV

/// class with muon info
class Muon : public TObject {
 public: 
  /// # of standalone-muon hits
   Int_t Nmu_hits;
   /// # of chambers with matched segments 
   /// (ie. # of muon stations contributing segments to the track)
   Int_t Nmatches;
   /// isolation for different cone sizes;
   Float_t SumPtIso[N_CONESIZE];
   Float_t NtrkIso[N_CONESIZE];
   Float_t ECALIso[N_CONESIZE];
   Float_t HCALIso[N_CONESIZE];
   // tracking information obtained with 
   // tracker-only, global-fit, TPFMS (tracker plus first muon station)
   // picky, cocktail and TMR algorithms 
   Track tracker;
   Track global;
   Track tpfms;
   Track picky;
   Track dyt;
   Track cocktail;
   Track tmr;

   // quality flags
   Bool_t GlobalMuonPromptTight;
   Bool_t TMLastStationLoose;
   Bool_t TMLastStationTight;
   Bool_t TMLastStationAngTight;
   Bool_t AllGlobalMuons;
   Bool_t AllStandAloneMuons;
   Bool_t AllTrackerMuons;

   Muon();
   ~Muon();
    
   ClassDef(Muon, 1)
};



/// class containing per-event Wprime info
class Event : public TObject {
 public: 
  /// event #
   Int_t evt_no;
   /// run #
   Int_t run_no;
   /// luminosity section (aka block)
   Int_t LS_no;

   /// # of primary vertices
   Int_t Npv;
   /// # of primary vertices with beam-spot
   Int_t NpvBS;
     
   /// trigger decisions: -1: did not run, 0: fail, 1: accept
   signed char HLT_L1MuOpen;
   signed char HLT_L1Mu;
   signed char HLT_Mu3;
   signed char HLT_Mu5;
   signed char HLT_Mu7;
   signed char HLT_Mu9;
   signed char HLT_Mu11;
   signed char HLT_L2Mu5;
   signed char HLT_L2Mu9;
   signed char HLT_L2Mu11;
   signed char HLT_L2Mu15;
   signed char HLT_L2Mu25;
   
   
   void reset_triggers()
   {HLT_L1MuOpen = HLT_L1Mu = HLT_Mu3 = HLT_Mu5 = HLT_Mu7 = HLT_Mu9 = 
       HLT_Mu11 = HLT_L2Mu5 = HLT_L2Mu9 = HLT_L2Mu11 = HLT_L2Mu15 = 
       HLT_L2Mu25 = -1;}
   /// Particle-Flow MET
   TVector2 pfmet;
   //same as pfmet but with the pT (from global) of the hardest muon added back
   TVector2 pfmetaddmu; 
                        

   // Jets (class: TLorentzVector)
   //Calo jets
   TClonesArray * jet;
   //PF jets
   TClonesArray * pfjet;
   /// Muons (class: muon)
   TClonesArray * mu;
   // MC-truth muons (class: MCParticle)
   TClonesArray * mu_mc;
   // MC-truth neutrinos (class: MCParticle)
   TClonesArray * neu_mc;
   // MC-truth W (class: MCParticle)
   TClonesArray * w_mc;
   // MC-truth W' (class: MCParticle)
   TClonesArray * wp_mc;

   // constructor
   Event();  
   // destructor
   ~Event();
   
   ClassDef(Event, 1)
};


} // namespace wprime

#endif // #define _wprimeEvent_h_
