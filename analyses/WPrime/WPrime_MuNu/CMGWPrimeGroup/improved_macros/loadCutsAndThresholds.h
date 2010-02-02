#ifndef _load_cuts_h__
#define _load_cuts_h__

#include <TROOT.h>
#include <TClonesArray.h>
#include <TSystem.h>
#include <TFile.h>
#include <TTree.h>
#include <TMath.h>
#include <TH1F.h>
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
const unsigned deltaRIsoIndex = 2; //for the isolation container
const float PtTrackCut = 60 ;
const float OneMuPtTrackCut = 20; 
const float Chi2Cut = 10;
const float Delta_Phi = TMath::Pi() - 0.3;//min angle muon/jet for jet veto
const float Muon_Eta_Cut = 1.8;

// ++++++++++++++++++++++++++++++++Useful constants
const int Num_histo_sets = 6; // one new set of histograms after each cut
const int Num_trkAlgos = 3; // global, tracker, tev_1st
const string algo_desc[Num_trkAlgos] = {"gbl","trk","tev"};
const string cuts_desc[Num_histo_sets] = {"hlt","ptrange","1mu","iso", "jet", "qual"};
const bool debugme = false;

// +++++++++++++++++++++++++++++++muon-pt histogram parameters
const unsigned  nBinPtMu = 45; // 400; // 45; // 18; 200; 380; 
const float minPtMu = 100; // 100;
const float  maxPtMu = 1500; // 800; 2000;
// +++++++++++++++++++++++++Declare histograms 
TH1F * hPT[Num_histo_sets][Num_trkAlgos] = {0};






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
// if 1, get just the final counts, if 2
const int Num_histo_sets_chargePt = 1; 





// $$$$$$$$$$$$$$$$$$$$$ Optimization (option = 4)
// +++++++++++++++++++++parameters
const string optimdir = "forOptimization";
const int o_maxNumCuts = 1035;
const int o_NmuTrkPt = 10;
const float o_muonTrackPt[o_NmuTrkPt] = {5,10,15,20,25,30,35,40,45,50};
const int o_indexDeltaR[9] = {0,1,2,3,4,5,6,7,8};
const float o_pTIso[5] = {6,8,10,12,14};
const float o_ETJetVeto[5] = {80,90,100,110,120};
const float o_DphiJetVeto[5] = { TMath::Pi() - 0.1, 
                                 TMath::Pi() - 0.2, 
                                 TMath::Pi() - 0.3, 
                                 TMath::Pi() - 0.4,
                                 TMath::Pi() - 0.5};




// +++++++++++++++++++++++++Declare auxiliary methods
void getEff(float & eff, float & deff, float Num, float Denom);
void GetTheHardestMuon(const wprime::Event * ev,
                       wprime::Muon*& theMu);
void GetGlobalMuonPtMax(const wprime::Event * ev,
                        float& ptMaxMu, float& ptMaxMuTrackVeto);
unsigned NmuAboveThresh(float tracker_muon_pt, 
                        const wprime::Event * ev,
                        wprime::Muon*& theMu);
unsigned NjetAboveThresh(float threshold, const wprime::Event * ev);
unsigned NjetAboveThresh(float threshold, float delta_phi, 
                         const wprime::Muon * mu, const wprime::Event * ev);


// +++++++++++++++++++++++++ Declare cut methods
bool PassedHLT(const wprime::Event* ev, wprime::Muon*& the_mu);
bool IsMuonPtInRange(const wprime::Muon* mu,const string& algo_name,
                     const float& min_MuPt, const float& max_MuPt);
bool OnlyOneHighTrackPtMuon(const wprime::Event* ev, 
                            wprime::Muon*& the_mu,
                            const float& one_mu_pt_trkcut);
bool SumPtIsolation(const wprime::Muon* the_mu, 
                    const unsigned& detR_iso_index,
                    const float& sum_pt_cut);
bool ExeedMaxNumJetsOpposedToMu(const unsigned& max_jets_aboveThresh,
                              const float& et_jet_cut,  
                              const float& delta_phi,
                              const wprime::Muon* the_mu,
                                const wprime::Event* ev);
bool HasQuality(const wprime::Muon* mu, 
                const string& algo_name, const float& pttrk_cut,
                const float& chi2_cut,const float& muon_etacut);


#endif // #define _load_cuts_h__
