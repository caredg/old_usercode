// wprime stuff

#include "UserCode/CMGWPrimeGroup/root_macros/loadCutsAndThresholds.h"


// Calculate efficiencies
//------------------------------------------------------------------------
void getEff(float & eff, float & deff,float Num,float Denom)
{
  //------------------------------------------------------------------------
  eff = Num/Denom;
  deff = TMath::Sqrt(eff * (1-eff)/Denom);
}//---------------getEff

// Get the hardest muon (based on tracker-pt) [iMu is the index in ev->mu array]
//------------------------------------------------------------------------
void GetTheHardestMuon(const wprime::Event * ev, int & iMu)
{
//------------------------------------------------------------------------

    int nmuons = ev->mu->GetLast() + 1;
    float temp_muPT = -999;
    for(int j = 0; j != nmuons; ++j)
      {
	wprime::Muon * mu = (wprime::Muon *) ev->mu->At(j);
	float current_muPT = mu->tracker.p.Pt();
	if (current_muPT > temp_muPT) 
	  {
	    temp_muPT = current_muPT;
	    iMu = j;
	  }
      }
}//-----GetTheHardestMuon


// returns # of (global) muons with tracker-pt above <tracker_muon_pt>
//------------------------------------------------------------------------
unsigned NmuAboveThresh(float tracker_muon_pt, const wprime::Event * ev)
{
//------------------------------------------------------------------------
  unsigned N = 0;

  int nmuons = ev->mu->GetLast() + 1;
  for(int j = 0; j != nmuons; ++j)
    { // loop over muons
      wprime::Muon * mu = (wprime::Muon *) ev->mu->At(j);
      if(mu->tracker.p.Pt() > tracker_muon_pt)
	++N;
            
    } // loop over muons

  return N;
}//------------NumAboveThresh

// returns # of jets with Et above threshold
//------------------------------------------------------------------------
unsigned NjetAboveThresh(float threshold, const wprime::Event * ev)
{
  unsigned N = 0;
  
  int njets = ev->jet->GetLast() + 1;
  for(int j = 0; j != njets; ++j)
    { // loop over jets
      TLorentzVector * jet = (TLorentzVector *) ev->jet->At(j);
      if(jet->Et() > threshold)
	++N;
    } // loop over jets
  
  return N;
}//--------------NjetAboveThresh



// returns # of jets with Et above threshold 
// with angle greater than delta_phi from muon
//------------------------------------------------------------------------
unsigned NjetAboveThresh(float threshold, float delta_phi, 
			 const wprime::Muon * mu, const wprime::Event * ev)
{
  //------------------------------------------------------------------------
  unsigned N = 0;
  
  int njets = ev->jet->GetLast() + 1;
  for(int j = 0; j != njets; ++j)
    { // loop over jets
      TLorentzVector * jet = (TLorentzVector *) ev->jet->At(j);
      if(jet->Et() >threshold && jet->DeltaPhi(mu->tracker.p) > delta_phi)
	++N;
    } // loop over jets
  
  return N;

}//-----------------NjetAboveThresh


// true if HLT conditions are met
//-------------------------------------------------------------------
bool PassedHLT(const wprime::Event* ev)
{
  //-------------------------------------------------------------------
#if debugmemore
  cout << " Executing PassedHLT()" << endl;
#endif
  
  // here the triggers that are to be used
  bool HLT = (ev->HLT_Mu9 == 1);
  return HLT;
}//-------PassedHLT()


// check if muon is in pt-range for the different algorithms, fill isThere
//-------------------------------------------------------------------
void CheckMuonPtInRange(const wprime::Muon* mu, bool isThere[],
			float min_MuPt, float max_MuPt)
{
//-------------------------------------------------------------------
#if debugmemore
  cout << " Processing IsMuonPtInRange() " << endl;
#endif
  
  // Make sure there is a muon
  if(mu == 0) 
    {
      for(int i = 0; i != Num_trkAlgos; ++i)
	isThere[i] = false;

      return;
    }

  isThere[0] = mu->global.p.Pt() >=min_MuPt 
    && mu->global.p.Pt() <=max_MuPt;

  
  isThere[1] = mu->tracker.p.Pt() >=min_MuPt 
    && mu->tracker.p.Pt() <=max_MuPt;

  isThere[2] = mu->tev_1st.p.Pt() >=min_MuPt 
    && mu->tev_1st.p.Pt() <=max_MuPt;

}//-------IsMuonPtInRange



// true if only one muon with track pT > the threshold
//-------------------------------------------------------------------
bool OnlyOneHighTrackPtMuon(const wprime::Event* ev, float one_mu_pt_trkcut)
{
//-------------------------------------------------------------------
#if debugmemore
 cout << " Processing OnlyOneHighTrackPtMuon() " << endl;
#endif

 return (NmuAboveThresh(one_mu_pt_trkcut, ev) == 1);
}//----------ThereIsOneMuonOnly




// true if isolation requirements satisfied for muon
//-------------------------------------------------------------------
bool SumPtIsolation(const wprime::Muon* the_mu, 
		    unsigned detR_iso_index,
                    float sum_pt_cut)
{
//-------------------------------------------------------------------
#if debugmemore
  cout << " Processing SumPtIsolation() " << endl;
#endif

  return the_mu->SumPtIso[detR_iso_index] <= sum_pt_cut;
}// SumPtIsolation






// same as above for relative isolation
//-------------------------------------------------------------------
bool RelSumPtIsolation(const wprime::Muon* the_mu, 
		    unsigned detR_iso_index,
                    float rel_sum_pt_cut)
{
//-------------------------------------------------------------------
#if debugmemore
  cout << " Processing RelSumPtIsolation() " << endl;
#endif

  return (the_mu->SumPtIso[detR_iso_index]/the_mu->tracker.p.Pt())
    <= rel_sum_pt_cut;
}// SumPtIsolation





// combined Trk+ECAL+HCAL relative isolation
//-------------------------------------------------------------------
bool CombRelIsolation(const wprime::Muon* the_mu, 
                       unsigned detR_iso_index,
                       float rel_combiso_cut)
{
//-------------------------------------------------------------------
#if debugmemore
  cout << " Processing CombRelIsolation() " << endl;
#endif

  return (
      (the_mu->SumPtIso[detR_iso_index]+
       the_mu->ECALIso[detR_iso_index]+
       the_mu->HCALIso[detR_iso_index])/the_mu->tracker.p.Pt())
      <= rel_combiso_cut;
}// CombRelIsolation



// true if energetic Jet(s) found back to back with muon 
//-------------------------------------------------------------------
bool ExceedMaxNumJetsOpposedToMu(unsigned max_jets_aboveThresh,
				 float et_jet_cut,  
				 float delta_phi,
				 const wprime::Muon* the_mu,
				 const wprime::Event* ev)
{
//-------------------------------------------------------------------
#if debugmemore
  cout << " Processing ExeedMaxNumJetsOpposedToMu() " << endl;
#endif
  return (NjetAboveThresh(et_jet_cut, delta_phi, 
			  the_mu, ev) > max_jets_aboveThresh);

}//ExceedMaxNumEnergeticJets



// check if muon satisfies quality requirements for all tracking algorithms, fill goodQual
//-------------------------------------------------------------------
void CheckQuality(const wprime::Muon* mu, 
		  bool goodQual[], float pttrk_cut,
		  float chi2_cut, float muon_etacut)
{
//-------------------------------------------------------------------
#if debugmemore
  cout << " Processing HasQuality() " << endl;
#endif
  
  bool isHard = mu->tracker.p.Pt() >= pttrk_cut;
  bool checkqual = false;
  // recommendation on muon-ID by muon POG: 
  // - GlobalMuonPromptTight + Global + Tracker muon
  // - # of hits in tracker track >= 11
  // - |IP of tracker track| < 0.2 cm
  bool muonID = mu->GlobalMuonPromptTight && (mu->tracker.Ntrk_hits >= 11)
    && (TMath::Abs(mu->tracker.d0) < 0.2) && mu->AllTrackerMuons
    && mu->AllGlobalMuons;

  checkqual = ((mu->global.chi2 / mu->global.ndof)<chi2_cut)
    && TMath::Abs(mu->global.p.Eta()) < muon_etacut;

  // old value: muon's pt within range
  // new value: old value .AND. quality cuts

  goodQual[0] = goodQual[0] && checkqual && isHard && muonID;



  checkqual = ((mu->tracker.chi2 / mu->tracker.ndof) 
	       < chi2_cut)
    && TMath::Abs(mu->tracker.p.Eta()) < muon_etacut;
 
  goodQual[1] = goodQual[1] && checkqual && isHard && muonID;



  checkqual = ((mu->tev_1st.chi2 / mu->tev_1st.ndof) 
	       < chi2_cut) 
    && TMath::Abs(mu->tev_1st.p.Eta()) < muon_etacut;

  goodQual[2] = goodQual[2] && checkqual && isHard && muonID;

}//------HasQuality




