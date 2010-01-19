{

  gROOT->Reset();

  // Update the include path so that we can find wprimeEvent.cc
  TString incpath = gSystem->GetIncludePath();
  incpath.Append(" -I$CMSSW_BASE/src");
  gSystem->SetIncludePath(incpath.Data());

  // compile code
  gROOT->ProcessLine(".L UserCode/CMGWPrimeGroup/src/wprimeEvent.cc+");
  gROOT->ProcessLine(".L UserCode/CMGWPrimeGroup/root_macros/loadInputFiles.C+");
  gROOT->ProcessLine(".L UserCode/CMGWPrimeGroup/root_macros/loadCrossSections.C+");
  gROOT->ProcessLine(".L UserCode/CMGWPrimeGroup/root_macros/loadCuts.C+");
  gROOT->ProcessLine(".L UserCode/CMGWPrimeGroup/root_macros/util.C+");

  //gROOT->ProcessLine(".L UserCode/CMGWPrimeGroup/root_macros/GetDistributionGeneric.C+");
  gROOT->ProcessLine(".L UserCode/CMGWPrimeGroup/root_macros/GetCutsOptimization.C+");
  gROOT->ProcessLine(".L UserCode/CMGWPrimeGroup/root_macros/plotOptim.C+");

  float lumiPb = 100;

  TFile *fout = new TFile("Wprime_analysis_V54.root","recreate");

  vector<wprime::InputFile> qcd_files; vector<wprime::InputFile> z_files;
  vector<wprime::InputFile> w_files; vector<wprime::InputFile> top_files;
  vector<wprime::InputFile> wprime10_files;
  vector<wprime::InputFile> wprime15_files;
  vector<wprime::InputFile> wprime20_files;

  gROOT->ProcessLine("loadCrossSections(qcd_files, z_files, w_files,top_files, wprime10_files, wprime15_files, wprime20_files)");

  //put output in a new directory
  string outdir("forOptimization");
  gSystem->MakeDirectory(outdir.c_str());
  

  //GetCutsOptimization(type of files, outdirectory for the *.dat files,
  //sample dir, variable to be optimized,lower bound, upper bound, number
  //of points) 
  string dir = "QCD";
  gROOT->ProcessLine("loadInputFiles(dir, qcd_files, lumiPb)");
  gROOT->ProcessLine("GetCutsOptimization(qcd_files, outdir, dir,\"iso\",0.,20.,100)");

  string dir = "Z";
  gROOT->ProcessLine("loadInputFiles(dir, z_files, lumiPb)");
  gROOT->ProcessLine("GetCutsOptimization(z_files, outdir, dir,\"iso\",0.,20.,100)");

   dir = "W";
   gROOT->ProcessLine("loadInputFiles(dir, w_files, lumiPb)");
   gROOT->ProcessLine("GetCutsOptimization(w_files, outdir, dir,\"iso\",0.,20.,100)");
  
   dir = "Top";
   gROOT->ProcessLine("loadInputFiles(dir, top_files, lumiPb)");
  gROOT->ProcessLine("GetCutsOptimization(top_files, outdir, dir,\"iso\",0.,20.,100)");

   dir = "wprime10";
   gROOT->ProcessLine("loadInputFiles(dir, wprime10_files, lumiPb)");
   gROOT->ProcessLine("GetCutsOptimization(wprime10_files, outdir, dir, \"iso\", 0., 20., 100)"); 

//   dir = "wprime15";
//   gROOT->ProcessLine("loadInputFiles(dir, wprime15_files, lumiPb)");
//   gROOT->ProcessLine("GetDistributionGeneric(wprime15_files, fout, dir, out)");
 
//   dir = "wprime20";
//   gROOT->ProcessLine("loadInputFiles(dir, wprime20_files, lumiPb)");
//   gROOT->ProcessLine("GetDistributionGeneric(wprime20_files, fout, dir, out)"); 

   //if "all" is passed instead of a single background type,
   //all backgrounds are summed.
   //gROOT->ProcessLine("plotOptim(\"iso\", \"wprime10\", \"Top\")");
   gROOT->ProcessLine("plotOptim(\"iso\", \"wprime10\", \"all\")");
   

  //out.close(); 
  //fout->Close();
}
