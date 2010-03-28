#ifndef _ExecuteAnalysis_h_
#define _ExecuteAnalysis_h_

#include <stdio.h>
#include <math.h>
#include <string>
#include <iostream>
#include <fstream>
#include "TFile.h"
#include "TTree.h"
#include "TH1F.h"
#include "TRandom.h"
#include "TMath.h"
#include <vector>
#include <algorithm>
#include <limits>

//gROOT->Reset();

using namespace std;

// ++++++++++++++data structure to accumulate info about the file
struct InputFile{
  float x_sect; // cross-section in pb
  float Nprod_evt; // # of events produced
  float weight; // cross-section * integrated luminosity / (# of events produced)
  string pathname;
  string description; // sample description
  TTree * tree; // pointer to ROOT file
  //
  InputFile(float cross_sect)
  {
    x_sect = cross_sect; Nprod_evt = 0; weight = 0; tree = 0;
  }
};

// +++++++++++++++++++location of data files
string top_level_dir = "/localdata/ecarrera/analyses/wprime_to_wz/root_uples/";


// +++++++++++++++++++Variables to store Branch Addresses:
Int_t           W_flavor;
Int_t           Z_flavor;
Int_t           numberOfZs;
Float_t         WZ_invMassMinPz;

// +++++++++++++++++++useful constants
const double PI    = 2.0 * acos(0.);
const double TWOPI = 2.0 * PI;
const bool debugme = false; //print stuff if active


// +++++++++++++++++++ Histogram Definitions
const int Num_histo_sets = 3;
TH1F * hWZInvMass[Num_histo_sets];



// +++++++++++++++++++ Muon ID Cuts

const float cutZMuon_pt = 10.;
const float cutWMuon_pt = 20.;
const float cutWmuD0 = 8.;
const float cutWmuCombRelIso = 0.1;


// +++++++++++++++++++ Electron ID Cuts

const float cutBarrel = 1.479;
const float cutEndcap = 1.55;

const float cutZElectron_pt = 15;
const float cutWElectron_pt = 20;

const float cutWElCombRelIso = 0.1;

//const float NOCUT = numeric_limits<float>::max();
//just make a big number for now
const float NOCUT = 99999999;


// Cut arrays are for {barrel, endcap}
const float cutDeltaEtaIn[]  = {0.005, 0.007};
const float cutDeltaPhiIn[]  = {0.040, 0.040};
const float cutSigmaEtaEta[] = {0.011, NOCUT};
const float cutEOverP[]      = {0.760, 0.680};
const float cutHOverE[]      = {0.016, 0.025};
const float cutCalRelIso[]   = {0.100, 0.160};
const float cutTrackRelIso[] = {0.100, 0.100};

// ++++++++++++++++++++ Other Cuts

const float minSeparation = 0.1;




#endif//#define _ExecuteAnalysis_h_
