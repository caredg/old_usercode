#ifndef _load_cuts_h__
#define _load_cuts_h__

#include <TROOT.h>
#include <TClonesArray.h>
#include <TSystem.h>
#include <TFile.h>
#include <TTree.h>
#include <TMath.h>
#include <TH1F.h>

#include <iostream>
#include <fstream>
#include <string>
#include <sstream>

#include "UserCode/CMGWPrimeGroup/interface/wprimeEvent.h"


// $$$$$$$$$$$$$$$$$$$$$$$ Main Study
// +++++++++++++++++++++++++++++++Analysis Thresholds:
const float EtJetCut = 100; 
const unsigned MaxNjetsAboveThresh = 0;
const float SumPtCut = 10; // Cone DeltaR =0.3; 
const float CombRelCut = 0.15; // Cone DeltaR =0.3; 
const unsigned deltaRIsoIndex = 2; //for the isolation container
const float PtTrackCut = 60 ;
const float OneMuPtTrackCut = 20; 
const float Chi2Cut = 10;
const float Delta_Phi = TMath::Pi() - 0.3;//min angle muon/jet for jet veto
const float Muon_Eta_Cut = 1.8;

// ++++++++++++++++++++++++++++++++Useful constants
const int Num_histo_sets = 6; // one new set of histograms after each cut
const int Num_trkAlgos = 3; // global, tracker, tev_1st

// use this for histogram names
const string algo_desc_short[Num_trkAlgos] = {"gbl","trk","tev"};
const string cuts_desc_short[Num_histo_sets] = {"hlt","1mu","ptrange","iso", "jet", "qual"};
// use this for histogram descriptions
const string algo_desc_long[Num_trkAlgos] = {"global", "tracker", "TPFMS"};
const string cuts_desc_long[Num_histo_sets]= {"HLT_Mu9,", "Pt within range,", 
					     "1 muon only,","Isolation,", 
					     "Jet Veto,", "Quality,"};

#define debugme  0
#define debugmemore 0
#define dumpHighPtMuons 0

// +++++++++++++++++++++++++++++++muon-pt histogram parameters
const unsigned  nBinPtMu = 45; // 400; // 45; // 18; 200; 380; 
const float minPtMu = 100; // 100;
const float  maxPtMu = 1500; // 800; 2000;
// +++++++++++++++++++++++++Declare histograms 
TH1F * hPT[Num_histo_sets][Num_trkAlgos] = {{0}};






// $$$$$$$$$$$$$$$$$$$$$$$ JetIso study (option = 2)

// ++++++++++++++++ JetIso study: jet distributions (# of jets and Et)
static const unsigned nBinNJets = 50;
static const float minNJets = -0.5;
static const float maxNJets = 49.5;
static const unsigned nBinEtJets = 150;
static const float minEtJets = 0;
static const float maxEtJets = 300;
// ++++++++++++++++ JetIso study: jet-activity veto (# of jets and Et)
static const unsigned nBinEtJets_veto = 10;
static const float minEtJets_veto = 0;
static const float maxEtJets_veto = 200;
static const unsigned nBinNJets_veto = 10;
static const float minNJets_veto = -0.5;
static const float maxNJets_veto = 9.5;
static const unsigned nBinSumPt = 120;
static const float minSumPt = 0;
static const float maxSumPt = 600;
// +++++++++++++++++++++++++Declare histograms
TH1F* hNMu, *hPtMaxMu, *hPtMaxMuTrackVeto, 
    *hPtMaxMuJetVeto, *hPtMaxMuTrackVetoJetVeto;






// $$$$$$$$$$$$$$$$$$$$$ Charge Asymmetry (option = 3)
// +++++++++++++++++++++++++Declare histograms
TH1F * hPTplus[Num_trkAlgos] = {0};
TH1F * hPTminus[Num_trkAlgos] = {0};
// ++++++++++++++++++++++++++++++++Useful constants
// if 1, create charge-asymmetry histograms after all cuts
// if 2, also create histograms before last cut, and so on
const int Num_histo_sets_chargePt = 1; 






// +++++++++++++++++++++++++Declare auxiliary methods

// Calculate efficiencies
void getEff(float & eff, float & deff, float Num, float Denom);

// Get the hardest muon (based on tracker-pt) [theMu is the index in ev->mu array]
void GetTheHardestMuon(const wprime::Event * ev, int & theMu);

// returns # of (global) muons with tracker-pt above <tracker_muon_pt>
unsigned NmuAboveThresh(float tracker_muon_pt, const wprime::Event * ev);

// returns # of jets with Et above threshold
unsigned NjetAboveThresh(float threshold, const wprime::Event * ev);

// returns # of jets with Et above 
// threshold with angle greater than
// delta_phi from muon
unsigned NjetAboveThresh(float threshold, float delta_phi, 
                         const wprime::Muon * mu, const wprime::Event * ev);


// +++++++++++++++++++++++++ Declare cut methods
// true if HLT conditions are met
bool PassedHLT(const wprime::Event* ev);

// check if muon is in pt-range for the different algorithms, fill isThere
void CheckMuonPtInRange(const wprime::Muon* mu, bool isThere[],
			float min_MuPt, float max_MuPt);

// true if only one muon with track pT > the threshold
bool OnlyOneHighTrackPtMuon(const wprime::Event* ev, float one_mu_pt_trkcut);

// true if isolation requirements satisfied for muon
bool SumPtIsolation(const wprime::Muon* the_mu, 
                    unsigned detR_iso_index,
                    float sum_pt_cut);

// same as above for relative isolation
bool RelSumPtIsolation(const wprime::Muon* the_mu, unsigned detR_iso_index,
		       float rel_sum_pt_cut);

// combined relative isolation
bool CombRelIsolation(const wprime::Muon* the_mu, 
                       unsigned detR_iso_index,
                      float rel_combiso_cut);


// true if energetic Jet(s) found back to back with muon 
bool ExceedMaxNumJetsOpposedToMu(unsigned max_jets_aboveThresh, float et_jet_cut,
				 float delta_phi, const wprime::Muon* the_mu,
				 const wprime::Event* ev);

// check if muon satisfies quality requirements for all tracking algorithms, 
// fill goodQual
void CheckQuality(const wprime::Muon* mu, bool goodQual[],
		  float pttrk_cut, float chi2_cut, float muon_etacut);


#endif // #define _load_cuts_h__
