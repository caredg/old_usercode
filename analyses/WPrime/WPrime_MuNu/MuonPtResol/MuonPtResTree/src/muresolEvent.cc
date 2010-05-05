#include "MuonPtResol/MuonPtResTree/src/muresolEvent_LinkDef.h"

using namespace muresol;

// constructor
Event::Event()
{
  jet = new TClonesArray("TLorentzVector");
  mu = new TClonesArray("muresol::Muon");
  mu_mc = new TClonesArray("muresol::MCParticle");
  neu_mc = new TClonesArray("muresol::MCParticle");
  w_mc = new TClonesArray("muresol::MCParticle");
  wp_mc = new TClonesArray("muresol::MCParticle");
}

// destructor
Event::~Event()
{
  delete jet; delete mu; 
  delete mu_mc; delete neu_mc; delete w_mc; delete wp_mc; 
}


JobInfo::JobInfo()
{ 
  HLTversion = RECOversion = sample = "invalid"; Nprod_evt = 0;
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




