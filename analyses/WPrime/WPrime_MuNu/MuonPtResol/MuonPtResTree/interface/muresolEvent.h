#ifndef _muresolEvent_h_
#define _muresolEvent_h_

// root stuff
#include <TROOT.h>
#include <TObject.h>
#include <TLorentzVector.h>
#include <TBranch.h>
#include <TClonesArray.h>

// std library
#include <string>

namespace muresol{

/// class containing per-job (or file) Muresol info
class JobInfo : public TObject {
 public: 
  /// HLT version
   std::string HLTversion;
   /// RECO version
   std::string RECOversion;
   /// sample type (e.g. QCD50_80)
   std::string sample;
   /// # of events produced BEFORE filtering
   Int_t Nprod_evt;
	  
   JobInfo();
   ~JobInfo();
	  
    ClassDef(JobInfo, 1);

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
   /// uncertainty on pt
   Float_t dpt;
   /// uncertainty on q/p
   Float_t dq_over_p;
   /// Ndof
   Int_t ndof;
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
   /// isolation for different cone sizes;
   Float_t SumPtIso[N_CONESIZE];
   Float_t NtrkIso[N_CONESIZE];
   Float_t ECALIso[N_CONESIZE];
   Float_t HCALIso[N_CONESIZE];
   // tracking information obtained with 
   // tracker-only, global-fit or TeV-1st muon station algorithms 
   Track tracker;
   Track global;
   Track tev_1st;
   Track tev_oddhits;
   
    
   Muon();
   ~Muon();
    
   ClassDef(Muon, 1)
};


/// class containing per-event Muresol info
class Event : public TObject {
 public: 
  /// event #
   Int_t evt_no;
   /// run #
   Int_t run_no;
     
   /// trigger decisions
   Bool_t HLT_L1MuOpen;
   Bool_t HLT_L1Mu;
   Bool_t HLT_Mu3;
   Bool_t HLT_Mu5;
   Bool_t HLT_Mu9;
   
   /// Particle-Flow MET
   TVector2 pfmet;

   /// Jets (class: TLorentzVector)
   TClonesArray * jet;
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


} // namespace muresol

#endif // #define _muresolEvent_h_
