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




// Get the pT of the hardest global Muon regarless
// of its trackPT (ptMaxMu) and if it passes the 
// PtTrackCut (ptMaxMuTrackVeto)
//------------------------------------------------------------------------
void GetGlobalMuonPtMax(const wprime::Event * ev,
                        float& ptMaxMu, float& ptMaxMuTrackVeto)
{
//------------------------------------------------------------------------
    int nmuons = ev->mu->GetLast() + 1;
    for(int j = 0; j != nmuons; ++j)
	  { // loop over muons
	    wprime::Muon * mu = (wprime::Muon *) ev->mu->At(j);
	    float ptmu =  mu->global.p.Pt();
	    if(ptmu > ptMaxMu) ptMaxMu = ptmu;
        if( mu->tracker.p.Pt() > PtTrackCut ){
            if(ptmu > ptMaxMuTrackVeto) {ptMaxMuTrackVeto = ptmu;}
        }
	  } // loop over muons

}//------------GetGlobalMuonPtMax





// Get the hardest muon
//------------------------------------------------------------------------
void GetTheHardestMuon(const wprime::Event * ev,
                       wprime::Muon*& theMu)
{
//------------------------------------------------------------------------

    int nmuons = ev->mu->GetLast() + 1;
    float temp_muPT = 0;
    for(int j = 0; j != nmuons; ++j){
        wprime::Muon * mu = (wprime::Muon *) ev->mu->At(j);
        float current_muPT = mu->tracker.p.Pt();
        if (current_muPT > temp_muPT) {
            temp_muPT = current_muPT;
            theMu = mu;
        }
    }
}//-----GetTheHardestMuon





// returns # of (global) muons with tracker-pt above <tracker_muon_pt>
//------------------------------------------------------------------------
unsigned NmuAboveThresh(float tracker_muon_pt, 
                        const wprime::Event * ev,
                        wprime::Muon*& theMu)
{
//------------------------------------------------------------------------
  unsigned N = 0;

  int nmuons = ev->mu->GetLast() + 1;
  for(int j = 0; j != nmuons; ++j)
    { // loop over muons
      wprime::Muon * mu = (wprime::Muon *) ev->mu->At(j);
      if(mu->tracker.p.Pt() > tracker_muon_pt){
          theMu = mu;
          ++N;
      }
      
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






// returns # of jets with Et above 
// threshold with angle greater than
// delta_phi from muon
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
      if(jet->Et() > threshold && jet->DeltaPhi(mu->tracker.p) > delta_phi)
	++N;
    } // loop over jets

  return N;

}//-----------------NjetAboveThresh




//Check if HLT conditions are met.
//-------------------------------------------------------------------
bool PassedHLT(const wprime::Event* ev,wprime::Muon*& the_mu)
{
//-------------------------------------------------------------------
    if(debugme) cout<<"Executing PassedHLT()"<<endl;
    
    //here the triggers that are to be used
    bool HLT = ev->HLT_Mu9;
    //Pass the hardest muon in event
    GetTheHardestMuon(ev,the_mu);
    
    return HLT;
}//-------PassedHLT()





//Check to see if the hardest muon is in range for
//the different algorithms
//-------------------------------------------------------------------
bool IsMuonPtInRange(const wprime::Muon* mu,const string& algo_name,
                     const float& min_MuPt, const float& max_MuPt)
{
//-------------------------------------------------------------------
    if(debugme) cout<<"Processing IsMuonPtInRange()"<<endl;
    bool checkptrange = false;
    
    //Make sure there is a muon
    if(mu == 0) return checkptrange;
  
    if (!strcmp(algo_name.c_str(),"gbl")){
        checkptrange = mu->global.p.Pt() >=min_MuPt 
            && mu->global.p.Pt() <=max_MuPt;
    }
    else if (!strcmp(algo_name.c_str(),"trk")){
        checkptrange = mu->tracker.p.Pt() >=min_MuPt 
            && mu->tracker.p.Pt() <=max_MuPt;
    }
    else if (!strcmp(algo_name.c_str(),"tev")){
        checkptrange = mu->tev_1st.p.Pt() >=min_MuPt 
            && mu->tev_1st.p.Pt() <=max_MuPt;
    }
    else cout<<"Algo name not recognized"<<endl;

    return checkptrange;


}//-------IsMuonPtInRange



// Require only one muon with track pT > the threshold
//-------------------------------------------------------------------
bool OnlyOneHighTrackPtMuon(const wprime::Event* ev, 
                            wprime::Muon*& the_mu,
                            const float& one_mu_pt_trkcut)
{
//-------------------------------------------------------------------

    if(debugme) cout<<"Processing OnlyOneHighTrackPtMuon()"<<endl;
    bool OneMuon = (NmuAboveThresh(one_mu_pt_trkcut, ev, the_mu) == 1);
    return OneMuon;
    
}//----------ThereIsOneMuonOnly




//Require isolation for the Muon found
//-------------------------------------------------------------------
bool SumPtIsolation(const wprime::Muon* the_mu, 
                    const unsigned& detR_iso_index,
                    const float& sum_pt_cut)
{
//-------------------------------------------------------------------
    if(debugme) cout<<"Processing SumPtIsolation()"<<endl;
    bool iso = false;
    iso = the_mu->SumPtIso[detR_iso_index] <= sum_pt_cut;
    return iso;
}// SumPtIsolation



//Require no energetic Jet back to back with the muon 
//-------------------------------------------------------------------
bool ExeedMaxNumJetsOpposedToMu(const unsigned& max_jets_aboveThresh,
                              const float& et_jet_cut,  
                              const float& delta_phi,
                              const wprime::Muon* the_mu,
                              const wprime::Event* ev)
{
//-------------------------------------------------------------------
    if(debugme) cout<<"Processing ExeedMaxNumJetsOpposedToMu()"<<endl;
    bool JetActivity = (NjetAboveThresh(et_jet_cut, delta_phi, 
						the_mu, ev) > max_jets_aboveThresh);
    return JetActivity;

}//ExeedMaxNumEnergeticJets



//Require track quality for the muon.
//-------------------------------------------------------------------
bool HasQuality(const wprime::Muon* mu, 
                const string& algo_name, const float& pttrk_cut,
                const float& chi2_cut,const float& muon_etacut)
{
//-------------------------------------------------------------------
    if(debugme) cout<<"Processing HasQuality()"<<endl;
    bool isHard = mu->tracker.p.Pt() >= pttrk_cut;
    bool checkqual = false;
    if (!strcmp(algo_name.c_str(),"gbl")){
        checkqual = ((mu->global.chi2 / mu->global.ndof)<chi2_cut)
	      && TMath::Abs(mu->global.p.Eta()) < muon_etacut;
    }
    else if (!strcmp(algo_name.c_str(),"trk")){
        checkqual = ((mu->tracker.chi2 / mu->tracker.ndof) 
			      < chi2_cut)
	      && TMath::Abs(mu->tracker.p.Eta()) < muon_etacut;
    }
    else if (!strcmp(algo_name.c_str(),"tev")){
        checkqual = ((mu->tev_1st.chi2 / mu->tev_1st.ndof) 
				     < chi2_cut) 
            && TMath::Abs(mu->tev_1st.p.Eta()) < muon_etacut;
    }
    else cout<<"Algo name not recognized"<<endl;

    bool overallqual = isHard && checkqual;

    return overallqual;

}//------HasQuality




