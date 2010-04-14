{

  gROOT->Reset();

  // Update the include path so that we can find wprimeEvent.cc
  TString incpath = gSystem->GetIncludePath();
  incpath.Append(" -I$CMSSW_BASE/src");
  gSystem->SetIncludePath(incpath.Data());

// flag for detector conditions 
  // options:
  // 1 -> 21x RECO, ideal conditions
  // 2 -> 2212 RECO, 50 pb-1 alignment for Wprime and W (all other bgd samples same as option 1)
  unsigned int detector_conditions = 2;

  // compile code
  gROOT->ProcessLine(".L UserCode/CMGWPrimeGroup/src/wprimeEvent.cc+");
  gROOT->ProcessLine(".L UserCode/CMGWPrimeGroup/root_macros/loadInputFiles.C+");
  gROOT->ProcessLine(".L UserCode/CMGWPrimeGroup/root_macros/loadCrossSections.C+");
  gROOT->ProcessLine(".L UserCode/CMGWPrimeGroup/root_macros/loadCutsAndThresholds.C+");
  gROOT->ProcessLine(".L UserCode/CMGWPrimeGroup/root_macros/GetCutsOptimization.C+");

  float lumiPb = 100;

  TFile *fout = new TFile("Wprime_analysis.root","recreate");

  vector<wprime::InputFile> qcd_files; vector<wprime::InputFile> z_files;
  vector<wprime::InputFile> w_files; vector<wprime::InputFile> top_files;
  vector<wprime::InputFile> wprime10_files;
  vector<wprime::InputFile> wprime15_files;
  vector<wprime::InputFile> wprime20_files;

  gROOT->ProcessLine("loadCrossSections(qcd_files, z_files, w_files,top_files, wprime10_files, wprime15_files, wprime20_files, detector_conditions)");

  string outfile("event_counts.txt");
  ofstream out(outfile.c_str());
  if(!out) { 
    cout << "Cannot open file " << outfile << endl; 
    abort();
  } 


//switch on or off looping over muons; 
//DO NOT USE, IT IS NOT FUNCTIONAL HERE YET
bool loopMuons = false;
//Optimization option. Choose from the ones above
const int option = 100;

//make optimization directory
const string optimdir = "forOptimization";
if(!gSystem->OpenDirectory(optimdir.c_str())){
    gSystem->MakeDirectory(optimdir.c_str());
  }
else cout<<"Directory \""<<optimdir.c_str()<<"\" exists"<<endl;

// Options for the optimization study.
// Current cuts suscentible of optimization are:
// A. One muon pT cut (vary muon pt)
// B. Isolation (vary deltaR and sum_pt)
// C. Jet Veto (vary deltaPhi and jet et)
// D. Track Quality cut

//In probability notation, the optimization study options are:
// option = 100 : A 
// option = 200 : <BC|A> 
// option = 300 : <CB|A>
// option = 400 : <B|A>
// option = 500 : <C|A>
// option = 600 : <C|AB>
// option = 700 : <B|AC>
// option = 800 : <D|ABC>
// option = 900 : <D|ACB>

//Angle values are fixed in the header of the GetCutsOptimization macro.
//For the rest of the values, a min, max, and step value can be given,
//the program will then calculate an array of possible values.

//Muon track pt values
const float o_minMuonTrackPt = 0;
const float o_maxMuonTrackPt = 800;
const int o_nstepMuonTrackPt = 40;
//sum_pt for isolation
const float o_minPtIso = 0;
const float o_maxPtIso = 15;
const int o_nstepPtIso = 60;
//jet et for jet veto
const float o_minEtJet = 0;
const float o_maxEtJet = 200;
const int o_nstepEtJet = 40;
//delta phi for jet veto
//const float o_minDphiJet = TMath::Pi()-0.5;
//const float o_maxDphiJet = TMath::Pi()-0.1;
//const int o_nstepDphiJet = 5;
//Muon track pt values for qual
const float o_minQualPt = 0;
const float o_maxQualPt = 600;
const int o_nstepQualPt = 120;
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
else if (option == 200 || option == 300){
  o_min1 =  o_minPtIso;
  o_max1 =  o_maxPtIso;
  o_nstep1 =  o_nstepPtIso;
  o_min2 =  o_minEtJet;
  o_max2 =  o_maxEtJet;
  o_nstep2 =  o_nstepEtJet;
}
else if (option == 400 || option == 700){
  o_min1 =  o_minPtIso;
  o_max1 =  o_maxPtIso;
  o_nstep1 =  o_nstepPtIso;
}

else if (option == 500 || option == 600){
  o_min1 =  o_minEtJet;
  o_max1 =  o_maxEtJet;
  o_nstep1 =  o_nstepEtJet;
}
else if (option == 800 || option == 900){
  o_min1 =  o_minQualPt;
  o_max1 =  o_maxQualPt;
  o_nstep1 =  o_nstepQualPt;
  o_min2 =  o_minChi2cut;
  o_max2 =  o_maxChi2cut;
  o_nstep2 =  o_nstepChi2cut;
}

  
   string dir = "QCD";
   gROOT->ProcessLine("loadInputFiles(dir, qcd_files, lumiPb, detector_conditions)");
   gROOT->ProcessLine("GetCutsOptimization(qcd_files, dir,option,loopMuons, o_min1, o_max1, o_nstep1, o_min2, o_max2, o_nstep2)");
  
   dir = "Z";
   gROOT->ProcessLine("loadInputFiles(dir, z_files, lumiPb, detector_conditions)");
   gROOT->ProcessLine("GetCutsOptimization(z_files, dir,option,loopMuons, o_min1, o_max1, o_nstep1, o_min2, o_max2, o_nstep2)");

   dir = "W";
   gROOT->ProcessLine("loadInputFiles(dir, w_files, lumiPb, detector_conditions)");
   gROOT->ProcessLine("GetCutsOptimization(w_files, dir,option,loopMuons, o_min1, o_max1, o_nstep1, o_min2, o_max2, o_nstep2)");

   dir = "Wsignal";
   gROOT->ProcessLine("loadInputFiles(dir, w_files, lumiPb, detector_conditions)");
   gROOT->ProcessLine("GetCutsOptimization(w_files, dir,option,loopMuons, o_min1, o_max1, o_nstep1, o_min2, o_max2, o_nstep2)");

   dir = "Wback";
   gROOT->ProcessLine("loadInputFiles(dir, w_files, lumiPb, detector_conditions)");
   gROOT->ProcessLine("GetCutsOptimization(w_files, dir,option,loopMuons, o_min1, o_max1, o_nstep1, o_min2, o_max2, o_nstep2)");

  
   dir = "Top";
   gROOT->ProcessLine("loadInputFiles(dir, top_files, lumiPb, detector_conditions)");
   gROOT->ProcessLine("GetCutsOptimization(top_files, dir,option,loopMuons, o_min1, o_max1, o_nstep1, o_min2, o_max2, o_nstep2)");  

   dir = "wprime10";
   gROOT->ProcessLine("loadInputFiles(dir, wprime10_files, lumiPb, detector_conditions)");
   gROOT->ProcessLine("GetCutsOptimization(wprime10_files, dir,option,loopMuons, o_min1, o_max1, o_nstep1, o_min2, o_max2, o_nstep2)"); 

   dir = "wprime15";
   gROOT->ProcessLine("loadInputFiles(dir, wprime15_files, lumiPb, detector_conditions)");
   gROOT->ProcessLine("GetCutsOptimization(wprime15_files, dir,option,loopMuons, o_min1, o_max1, o_nstep1, o_min2, o_max2, o_nstep2)");
 
   dir = "wprime20";
   gROOT->ProcessLine("loadInputFiles(dir, wprime20_files, lumiPb, detector_conditions)");
   gROOT->ProcessLine("GetCutsOptimization(wprime20_files, dir,option,loopMuons, o_min1, o_max1, o_nstep1, o_min2, o_max2, o_nstep2)"); 

  out.close(); 
  fout->Close();
}
