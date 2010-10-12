#include "UserCode/CMGWPrimeGroup/src/wprimeEvent_LinkDef.h"

using namespace wprime;

// constructor
Event::Event()
{
  calojet = new TClonesArray("TLorentzVector");
  pfjet = new TClonesArray("TLorentzVector");
  mu = new TClonesArray("wprime::Muon");
  mu_mc = new TClonesArray("wprime::MCParticle");
  neu_mc = new TClonesArray("wprime::MCParticle");
  w_mc = new TClonesArray("wprime::MCParticle");
  wp_mc = new TClonesArray("wprime::MCParticle");
}

// destructor
Event::~Event()
{
    delete calojet; delete pfjet; delete mu; 
    delete mu_mc; delete neu_mc; delete w_mc; delete wp_mc; 
}


RunInfo::RunInfo()
{
  HLTmenu = HLTversion = "invalid"; run_no = Nproc_evt = -999;
}
RunInfo::~RunInfo(){}

JobInfo::JobInfo()
{ 
  RECOversion = sample = "invalid";
}
JobInfo::~JobInfo(){}

MCParticle::MCParticle()
{ 
  q = momId = status = -999;
}

MCParticle::MCParticle(TLorentzVector & p_, Int_t q_, Int_t momId_, 
				   Int_t status_)
{
  p = p_; q = q_; momId = momId_; status = status_;
}


MCParticle::~MCParticle(){}

Track::Track()
{
  q = 0; ndof = Ntrk_hits = Ntot_hits = -999;
  
  dpt = dq_over_p = chi2 = -9999;
}
Track::~Track(){}

Muon::Muon(){}
Muon::~Muon(){}





