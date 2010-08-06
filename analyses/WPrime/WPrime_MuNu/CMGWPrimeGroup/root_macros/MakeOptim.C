{

  gROOT->Reset();
  
  /*
  // Update the include path so that we can find wprimeEvent.cc
  TString incpath = gSystem->GetIncludePath();
  incpath.Append(" -I$CMSSW_BASE/src");
  gSystem->SetIncludePath(incpath.Data());
  */

  
  // compile code
  gROOT->ProcessLine(".L UserCode/CMGWPrimeGroup/src/wprimeEvent.cc+");
  gROOT->ProcessLine(".L UserCode/CMGWPrimeGroup/root_macros/loadInputFiles.C+");
  gROOT->ProcessLine(".L UserCode/CMGWPrimeGroup/root_macros/loadCutsAndThresholds.C+");
  gROOT->ProcessLine(".L UserCode/CMGWPrimeGroup/root_macros/GetCutsOptimization.C+");
  gROOT->ProcessLine(".x UserCode/CMGWPrimeGroup/root_macros/runOptim.C+");
}

