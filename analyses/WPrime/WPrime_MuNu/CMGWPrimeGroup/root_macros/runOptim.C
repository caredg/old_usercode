#include <TROOT.h>
#include <TFile.h>

#include <string>
#include <vector>
#include <iostream>
#include <fstream>

#include "UserCode/CMGWPrimeGroup/root_macros/input_file.h"

using std::cout; using std::endl;
using std::vector; using std::string;

extern void GetCutsOptimization(const wprime::InputFile& files,
                                TFile* fout,
                                const int option, const bool loopMuons,
                                const float o_min1, const float o_max1, 
                                const int o_nstep1, const float o_min2,
                                const float o_max2, const int o_nstep2);


extern void loadInputFiles(vector<wprime::InputFile> & files, float lumiPb);

void runOptim()
{
  //  float lumiPb = 100; // in pb^-1
  //float lumiPb = 0.255; // in pb^-1
    float lumiPb = 100; // in pb^-1

   vector<wprime::InputFile> all_files; 

 //Optimization option. Choose from the ones above
  //####### change both:XS
  const int option = 300;
  string namesopt[7] = {"1mu","1muiso","1mujet","1muisojet","1mujetiso","1muisojetqual","1mujetisoqual"};
  
  int optstrdim = int(option/100.)-1;
  const string optionstr = namesopt[optstrdim];
  const string optimdir = "forOptimization";
  string  histfile = optimdir +"/"+"h_" + optionstr+".root";
  TFile* fout = new TFile(histfile.c_str(),"recreate");
    //$$$$$$optimization
// Options for the optimization study.
// Current cuts suscentible of optimization are:
// A. One muon pT cut (vary muon pt)
// B. Isolation (vary deltaR and sum_pt)
// C. Jet Veto (vary deltaPhi and jet et)
// D. Track Quality cut

//In probability notation, the optimization study options are:
// option = 100 : A 
// option = 200 : <B|A>
// option = 300 : <C|A>
// option = 400 : <C|AB>
// option = 500 : <B|AC>
// option = 600 : <D|ABC>
// option = 700 : <D|ACB>

//Angle values are fixed in the header of the GetCutsOptimization macro.
//For the rest of the values, a min, max, and step value can be given,
//the program will then calculate an array of possible values.

//Muon track pt values
const float o_minMuonTrackPt = 0;
const float o_maxMuonTrackPt = 600;
const int o_nstepMuonTrackPt = 30;

//sum_pt for isolation
//const float o_minPtIso = 12;
//const float o_maxPtIso = 22;
//const int o_nstepPtIso = 5;
//combined relative isolation (use either sum_pt above or this one)
const float o_minPtIso = 0.05;
const float o_maxPtIso = 0.20;
const int o_nstepPtIso = 3;
//jet et for jet veto
const float o_minEtJet = 30;
const float o_maxEtJet = 90;
const int o_nstepEtJet = 2;
//delta phi for jet veto
//const float o_minDphiJet = TMath::Pi()-0.5;
//const float o_maxDphiJet = TMath::Pi()-0.1;
//const int o_nstepDphiJet = 5;
//Muon track pt values for qual
const float o_minQualPt = 20;
const float o_maxQualPt = 100;
const int o_nstepQualPt = 8;
//Chi2 cut values for qual
const float o_minChi2cut = 5;
const float o_maxChi2cut = 15;
const int o_nstepChi2cut = 10;


//Declare optimization parameter holders.  There are two sets because
//the program can perform on two variables. You shouldn't need to
//change anything below here.
float o_min1 = 0.0;
float o_max1 = 0.0;
int o_nstep1 = 0;
float o_min2 = 0.0;
float o_max2 = 0.0;
int o_nstep2 = 0;

if (option == 100){
  o_min1 =  o_minMuonTrackPt;
  o_max1 =  o_maxMuonTrackPt;
  o_nstep1 =  o_nstepMuonTrackPt;
}
else if (option == 200 || option == 500){
  o_min1 =  o_minPtIso;
  o_max1 =  o_maxPtIso;
  o_nstep1 =  o_nstepPtIso;
}

else if (option == 300 || option == 400){
  o_min1 =  o_minEtJet;
  o_max1 =  o_maxEtJet;
  o_nstep1 =  o_nstepEtJet;
}
else if (option == 600 || option == 700){
  o_min2 =  o_minQualPt;
  o_max2 =  o_maxQualPt;
  o_nstep2 =  o_nstepQualPt;
  o_min1 =  o_minChi2cut;
  o_max1 =  o_maxChi2cut;
  o_nstep1 =  o_nstepChi2cut;
}
  
  // switch on or off looping over muons;
  // if true, will only examine hardest muon in event (based on tracker-pt)
  bool highestPtMuonOnly = false; 

  loadInputFiles(all_files, lumiPb);
  
  vector<wprime::InputFile>::const_iterator it;
  for(it = all_files.begin(); it != all_files.end(); ++it)
      GetCutsOptimization(*it,fout,option, highestPtMuonOnly,
                          o_min1,o_max1,o_nstep1, o_min2, o_max2, o_nstep2);

  fout->Close();
}
