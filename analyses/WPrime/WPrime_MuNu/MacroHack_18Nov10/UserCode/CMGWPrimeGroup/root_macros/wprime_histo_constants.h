#ifndef _wprime_histo_constants_h__
#define _wprime_histo_constants_h__

#include <string>

// ++++++++++++++++++++++++++++++++Useful constants
const int Num_selection_cuts = 11; // one new set of histograms after each cut
const int Num_trkAlgos = 7; // global, tracker, tpfms, cocktail, picky, tmr, dyt
// making distributions for two flavors: muon-pt and Mt
const int Num_flavors = 2;
const int PT_INDEX = 0;
const int MT_INDEX = 1;

const string FLAVOR_NAME[Num_flavors] = {" Muon-pt analysis", " Mt analysis"};

// lower, upper limits for algorithms to consider when processing data
// default: processing only cocktail tracking (3)
const int MuAlgo_MIN = 3;
const int MuAlgo_MAX = 3;

const int NumPtMtThresholds = 6;

//float PtThreshold[NumPtMtThresholds] = {100, 120, 200, 300, 350, 400};
float PtThreshold[NumPtMtThresholds] = {0,25, 40, 100, 120, 400};
//float MtThreshold[NumPtMtThresholds] = {200, 300, 350, 400, 500, 600};
float MtThreshold[NumPtMtThresholds] = {0,50, 100, 200, 300, 600};

// cut-indices for the <NumPtMtThresholds> thresholds
int ThreshIndices[NumPtMtThresholds] = {5, 6, 7, 8, 9, 10};

// use this for histogram names
const std::string algo_desc_short[Num_trkAlgos] = {"gbl","trk","tpfms","ckt","pic","tmr","dyt"};
// use this for histogram descriptions
const std::string algo_desc_long[Num_trkAlgos] = {"global", "tracker", "TPFMS","cocktail","picky","TMR", "DYT"};

// use this for histogram names
const std::string cuts_desc_short[Num_selection_cuts] = {"hlt","qual","1mu","iso", "jetmet", "thr1", "thr2", "thr3", "thr4", "thr5", "thr6"};
// use this for histogram descriptions
const std::string cuts_desc_long[Num_selection_cuts]= {"Single-muon HLT", "Quality", "1 muon only","Isolation", "Jet Veto OR kinematic cuts", "Threshold 1", "Threshold 2", "Threshold 3", "Threshold 4", "Threshold 5", "Threshold 6"}; 

// last histogram after all cuts to be used for fit
const std::string final_histo_desc = cuts_desc_short[ThreshIndices[0]];

#endif // #ifndef _wprime_histo_constants_h__
