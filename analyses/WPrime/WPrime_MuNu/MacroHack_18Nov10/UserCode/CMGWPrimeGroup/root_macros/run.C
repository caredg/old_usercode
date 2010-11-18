#include <TROOT.h>
#include <TFile.h>
#include <TH1F.h>

#include <string>
#include <vector>
#include <iostream>
#include <fstream>

#include "UserCode/CMGWPrimeGroup/root_macros/input_file.h"
#include "UserCode/CMGWPrimeGroup/root_macros/loadCutsAndThresholds.h"

using std::cout; using std::endl;
using std::vector; using std::string;

extern void GetDistributionGeneric(wprime::InputFile & file, 
				   TFile *fout, ofstream & out, 
				   const int option = 1, 
				   const bool highestPtMuonOnly = true);

extern void loadInputFiles(vector<wprime::InputFile> & files, float & lumiPb);

void run()
{
  float lumiPb = -1; // in pb^-1, to be retrieved from samples_cross_sections.txt

  TFile *fout = new TFile("Wprime_analysis.root","recreate");

  vector<wprime::InputFile> all_files; 

  string outfile("event_counts.txt");
  ofstream out(outfile.c_str());
  if(!out) { 
    cout << "Cannot open file " << outfile << endl; 
    abort();
  } 

  // option = 1 : GetMuonPtDistribution
  // option = 2 : GetChargePtDistribution
  const int option = 1;
  
  // switch on or off looping over muons;
  // if true, will only examine hardest muon in event (based on tracker-pt)
  bool highestPtMuonOnly = false; 

  loadInputFiles(all_files, lumiPb);
  
  vector<wprime::InputFile>::iterator it;
  for(it = all_files.begin(); it != all_files.end(); ++it)
    GetDistributionGeneric(*it, fout, out, option, highestPtMuonOnly);

  TH1F * h = new TH1F("lumi_ipb", "Integrated luminosity in pb^{-1}", 1, 0, 1);
  h->SetBinContent(1, lumiPb);

  fout->cd();
  h->Write();

  fout->Close();

  cout << " Please find event-count summary after selection cuts in file " 
       << outfile << endl;

  out << "\n\n\n Event count summary " << endl;
  for(int f = 0; f != Num_flavors; ++f){ // loop over pt & mt

    out << "\n" << FLAVOR_NAME[f] << endl;
    
    for (int mual = MuAlgo_MIN; mual <= MuAlgo_MAX; ++mual){
      // loop over tracking algorithms
      out << " --------------------------------------------" << endl;
      out << " " << algo_desc_long[mual] << " muon reconstruction" << endl;
      out << " # of events after selection cuts for " << lumiPb << " pb^-1:"
	   << endl;

      int thresh_counter = -1;
      for(int i = 0; i != NumPtMtThresholds; ++i){
	int index = ThreshIndices[i];
	++thresh_counter;

	if (PT_INDEX == f)
	  out << "\n Muon-pt > " << PtThreshold[thresh_counter];
	else if(MT_INDEX == f)
	  out << "\n Mt > " << MtThreshold[thresh_counter];
	else
	  abort();
	out << " GeV " << endl;

	float N_SM = 0; 
	vector<wprime::InputFile>::const_iterator it;
	for(it = all_files.begin(); it != all_files.end(); ++it)
	  { // loop over samples
	    string sample = (*it).samplename;
	    out<< " "<< sample << ": " 
		<< (*it).Nexp_evt_cut[index][mual][f]
		<< " (eff = " << 100.*((*it).eff[index][mual][f]) 
		<< " +- " << 100.*((*it).deff[index][mual][f])
		<< " %) " << endl;

	    if(sample == "W" || sample == "Wlowpt" || sample == "QCD" || sample == "Z" || 
	       sample == "Top")
	      N_SM += (*it).Nexp_evt_cut[index][mual][f];
	    
	    if(sample == "Top")
	      out << " Total # of SM (W + QCD + Z + Top) events: " 
		   << N_SM << endl;	 
	  } // loop over samples

      } // loop over thresholds

    } // loop over tracking algorithms

  } // loop over pt & mt

  out.close();


}
