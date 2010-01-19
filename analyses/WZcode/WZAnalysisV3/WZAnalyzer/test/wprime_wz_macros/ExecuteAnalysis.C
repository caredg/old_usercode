//------------------------------------------------
// Author: Edgar Carrera
// 2010-01-13
// This macro will be used mostly for analyzing
// Wprime -> WZ -> lllnu events but it should
// work on any root-uple of events with
// variables created by the "official" CMS WZ code
// 
//------------------------------------------------

#define ExecuteAnalysis_cxx
#include <ExecuteAnalysis.h>

//--------------------------------------------------------------
void getEff(float & eff, float & deff, float Num, float Denom)
{
//--------------------------------------------------------------
  eff = Num/Denom;
  deff = TMath::Sqrt(eff * (1-eff)/Denom);

}//getEff




//--------------------------------------------------------------
string convertIntToStr(int number)
{
//--------------------------------------------------------------
   stringstream ss;
   ss << number;
   return ss.str();
}



//--------------------------------------------------------------
double deltaPhi(double phi1, double phi2)
{
//--------------------------------------------------------------

  double phi=fabs(phi1-phi2);
  return (phi <= PI) ? phi : TWOPI - phi;
}


//--------------------------------------------------------------
void Declare_Histos()
{
//--------------------------------------------------------------

  if (debugme) cout<<"Declare histos"<<endl;

  //InvMass histos
  int WZinvMin = 0;
  int WZinvMax = 1500;
  int WZbin = 150;
  hWZInvMass[0] = new TH1F("hWZInvMass_hlt","Reconstructed WZ Invariant Mass (hlt);m_{WZ} (GeV);",WZbin,WZinvMin,WZinvMax);
  hWZInvMass[1] = new TH1F("hWZInvMass_flavor","Reconstructed WZ Invariant Mass (flavor);m_{WZ} (GeV);",WZbin,WZinvMin,WZinvMax);
  hWZInvMass[2] = new TH1F("hWZInvMass_numZs","Reconstructed WZ Invariant Mass (numZs);m_{WZ} (GeV);",WZbin,WZinvMin,WZinvMax);
  
}//Declare_Histos


//Just a function to calculate DeltaR
//--------------------------------------------------------------
double deltaR(double eta1, double phi1, double eta2, double phi2)
{
//--------------------------------------------------------------

  double deta = eta1 - eta2;
  double dphi = deltaPhi(phi1, phi2);
  return sqrt(deta * deta + dphi * dphi);
}




//check if the input files are there
//----------------------------------------------------------
int Check_Files(unsigned Nfiles, vector<InputFile> & files)
{
//-----------------------------------------------------------

  if(Nfiles != files.size())
    {
      cerr << " *** Expected " << Nfiles << " files, got "
	   << files.size() << " instead..." << endl;
      return -1;
    }  
  return 0;
}





//Recruit files that are numbered from a given sample
//----------------------------------------------------------
void RecruitOrderedFiles(vector<InputFile> & files, const int& Nfiles,
                         const int& filenum_low, const int& filenum_step,
                         const string& mask1,
                         const string& mask2, const string& file_desc)
{
//-----------------------------------------------------------

    if(Check_Files(Nfiles, files)) {
        cout<<"WARNING!!!!! No Files were found for "<<
            file_desc<<"sample"<<endl;
        return;
    }

    int filenum = filenum_low;
    for(int i = 0; i != Nfiles; ++i){
        filenum = filenum_low + i*filenum_step;
        files[i].pathname = top_level_dir + mask1 + convertIntToStr(filenum) + mask2;
        files[i].description = file_desc;
    }


}// --RecruitOrderedFiles



//load the input files from the top_level_dir
//-----------------------------------------------------------
void Load_Input_Files(string file_desc, 
                    vector<InputFile> & files, 
                    float lumiPb)
{
//-----------------------------------------------------------
  if (debugme)cout<<"Loading input Files...."<<endl;

  cout << "\n Processing " << file_desc<< " files " << endl << endl;

//   if(file_desc == "QCD")
//     {
//       const int NfilesQCD = 09;
//       Nfiles = NfilesQCD;
//       if(Check_Files(Nfiles, files))
// 	return;

//       string low[NfilesQCD]= 
// 	{"100", "150", "200", "300", "400", "600", "800", "1200", "1600"};
//       string high[NfilesQCD]=
// 	{"150", "200", "300", "400", "600", "800","1200", "1600", "up"};
  
//       for(int i = 0; i != Nfiles; ++i)
// 	{
// 	  files[i].pathname = top_level_dir + string("QCD_")+low[i]+string("_")+high[i]+string("_212_Ideal_Minv_ptGlobMu.root");
// 	}
      
//     } // QCD

  //switch between different sample cases and
  //recruit the files. 
  int Nfiles = -1;
  if (!strcmp(file_desc.c_str(),"wzjj")){
          Nfiles = 20;
          const string mask1 = "WZjj_500event_";
          const string mask2 = "_outputTree.root";
          const int filenum_low = 201;
          const int filenum_step = 1;
          RecruitOrderedFiles(files,Nfiles,filenum_low,filenum_step,
                              mask1,mask2,file_desc);
  }
  else if (!strcmp(file_desc.c_str(),"wprime400")){
          Nfiles = 20;
          const string mask1 = "Wprime400_500event_";
          const string mask2 = "_outputTree.root";
          const int filenum_low = 101;
          const int filenum_step = 1;
          RecruitOrderedFiles(files,Nfiles,filenum_low,filenum_step,
                              mask1,mask2,file_desc);
  }
  else cout<<"No sample were found with the name "<<file_desc<<endl;

  
  
  
  //Loop over the files in order to get the correct tree
  //and the number of total events 
  for(int i = 0; i != Nfiles; ++i){ // loop over input files
      string pathname = files[i].pathname;
      
      TFile * file = new TFile(pathname.c_str());
      if(!(file->IsOpen()))	{
          cerr <<" *** Missing file: "<< pathname << " !!! "<<endl; 
          continue;
      }
      if (debugme)cout<<"Processing file: "<<pathname<<endl;
      
      files[i].tree = (TTree *) file->Get("WZ");
      
      TH1F* histNumEvents       = (TH1F  *)file->Get("numEvents");
      const float eventsAnalyzed  = histNumEvents->GetBinContent(1);
      delete histNumEvents;
      
      files[i].Nprod_evt = eventsAnalyzed;
      files[i].weight = lumiPb*(files[i].x_sect)/(files[i].Nprod_evt);
      
      
      cout <<"Events produced in file = "<<files[i].Nprod_evt
           <<",  # of entries = "<< files[i].tree->GetEntries() 
           << ", weight = " << files[i].weight << endl;
      
  } // loop over input files
  
  return;
}





//---------------------------------------------------------
int Load_Cross_Sections(vector<InputFile> & wzjj_files, 
                      vector<InputFile> & wprime400_files)
{
//---------------------------------------------------------
  if (debugme)cout<<"Loading cross sections...."<<endl;

  //CalHEP cross sections in pb:
  float xsec_wzjj = 0.1958;
  float xsec_wprime400 = 2.118995E-3;
  //float xsec_wprime500 = 3.7755E-3;
  //float xsec_wprime600 = 7.9664E-4;
  //float xsec_wprime700 = 5.3399E-4;
  //float xsec_wprime800 = 3.673218E-4;
  //float xsec_wprime900 = 2.58269E-4;
  //float xsec_wprime1000 = 1.84634E-4;
  //float xsec_wprime1100 = 1.32994E-4;
  //float xsec_wprime1200 = 9.702684E-5;
  


  //PYTHIA cross-section (pb)

  //qcd_files.push_back(InputFile(590000  , 450000)); // 100_150
  //qcd_files.push_back(InputFile( 83000  , 425000)); // 150_200
  //qcd_files.push_back(InputFile( 24000  , 440000)); // 200_300
  //qcd_files.push_back(InputFile(  3000  , 395000)); // 300_400
  //qcd_files.push_back(InputFile(   730  , 135000)); // 400_600
  //qcd_files.push_back(InputFile(    66  , 310000)); // 600_800
  //qcd_files.push_back(InputFile(    12  , 130000)); // 800_1200
  //qcd_files.push_back(InputFile(    0.63,  30000)); // 1200_1600
  //qcd_files.push_back(InputFile(    0.54,  25000)); // 1600_2000

  for (int j = 0; j<20;++j){
    wzjj_files.push_back(InputFile(xsec_wzjj)); // wzjj 
    wprime400_files.push_back(InputFile(xsec_wprime400)); // 400 GeV
  }
  

  return 0;
}





//Set the branch addresses that we need for the analysis
//-----------------------------------------------------------
void Set_Branch_Addresses(TTree* WZtree)
{
//-----------------------------------------------------------
  WZtree->SetBranchAddress("W_flavor",&W_flavor);
  WZtree->SetBranchAddress("Z_flavor",&Z_flavor);
//  WZtree->SetBranchAddress("triggerBitMask", &triggerBitMask);
  WZtree->SetBranchAddress("numberOfZs",&numberOfZs);
  WZtree->SetBranchAddress("WZ_invMassMinPz",&WZ_invMassMinPz);

}//Set_Branch_Addresses







//Fill Histograms
//-----------------------------------------------------------
void Fill_Histos(int index, float weight)
{
//-----------------------------------------------------------

  hWZInvMass[index]->Fill(WZ_invMassMinPz, weight);


}//Fill_Histos





//------------------------------------------------------------------------
void saveHistos(TFile * fout, string dir)
{
//------------------------------------------------------------------------

  fout->cd(); 
  fout->mkdir(dir.c_str()); 
  fout->cd(dir.c_str());

  for(int i = 0; i != Num_histo_sets; ++i){
        hWZInvMass[i]->Write();
  }

  return;

}//saveHistos






//Writing results to a txt file
//--------------------------------------------------------------------------
void printSummary(ofstream & out, const string& dir, 
                  const float& Nexp_evt, float Nexp_evt_cut[]) 
{ 
//------------------------------------------------------------------------
    if(debugme) cout<<"Writing results to a txt file"<<endl;

    out<<"$$$$$$$$$$$$$$$$$$$$$$$ Type of sample: "<<dir<<endl;
    out << " Total # of expected events = " << Nexp_evt << endl;
    
    for(int i = 0; i < Num_histo_sets; ++i){
        
        out <<"Cut # "<<i<<": expected evts = " << Nexp_evt_cut[i];
        
        //calculate efficiencies
        float eff, deff;
        if(i == 0)
            getEff(eff, deff, Nexp_evt_cut[i], Nexp_evt);
        else
            getEff(eff, deff, Nexp_evt_cut[i], Nexp_evt_cut[i-1]);
        out << ", Relative eff = "<<eff*100 << " +- " << deff*100 << "%";
        getEff(eff, deff, Nexp_evt_cut[i], Nexp_evt);
        out << ", Absolute eff = "<< eff*100 << " +- " << deff*100 << "%"
             << endl;
        
        //to do: put these results in a file
        
    } // loop over different cuts
    
    
}//printSummary






//Tabulate results after the cut has been passed
//-----------------------------------------------------------
void Tabulate_Me(int Num_surv_cut[], int& cut_index, 
                const float& weight)
{
//-----------------------------------------------------------
    if(debugme) cout<<"Tabulating results for cut_index = "
                    <<cut_index<<endl;

    //increase the number of events passing the cuts
    ++Num_surv_cut[cut_index];
    //fill the histograms
    Fill_Histos(cut_index,weight);

    //since the event has passed the cut,
    //increase the cut_index for the next cut
    if(debugme) cout<<"cut_index is now = "<<endl;
    ++cut_index;

}//Tabulate_Me





//Trigger requirements
//-----------------------------------------------------------
bool PassTriggers_Cut()
{
//-----------------------------------------------------------
    if(debugme) cout<<"Trigger requirements"<<endl;
    //implement it here
    return true;

}//--- PassTriggers_Cut()






//Check if there are valid W and Z particles in the event
//-----------------------------------------------------------
bool HasValidWandZ_Cut()
{
//-----------------------------------------------------------
    if(debugme) cout<<"Check if there are valid W and Z particles in the event"
                    <<endl;
    bool has_valid_W_and_Z = (Z_flavor && W_flavor);
    return has_valid_W_and_Z;
    
}//--- NotValidWandZ_Cut




//Check if there is more Zs than required
//-----------------------------------------------------------
bool ExeedMaxNumberOfZs_Cut(const int& max_num_Zs)
{
//-----------------------------------------------------------
    if(debugme) cout<<"Check if there is more Zs than required"<<endl;
    bool has_too_many_Zs = numberOfZs > max_num_Zs;
    return has_too_many_Zs;;

}//--- MaxNumberOfZs_Cut







//Get different types of distribution
//-----------------------------------------------------------
void Get_Distributions(vector<InputFile>& files, 
                       TFile *fout, string dir, ofstream & out)
{
//-----------------------------------------------------------
  if (debugme) cout<<"Get Distributions....."<<endl;
  
  
  Declare_Histos();
 
  int Nfiles = files.size();

  //initialize counters for expected (already weighted) 
  //number of events
  //total, and after each cut.
  float Nexp_evt = 0;
  float Nexp_evt_cut[Num_histo_sets] = {0};
  
  //loop over files
  for(int tr = 0; tr != Nfiles; ++tr){
    if(!files[tr].tree)
      continue;

    if (debugme) cout << "Processing file "<<files[tr].pathname<<endl;
    

    TTree *WZtree = files[tr].tree;    
    int nevents = WZtree->GetEntries();
    float weight = files[tr].weight;

    //get the variables that we need right here:
    Set_Branch_Addresses(WZtree);

    //counter (unweighted) events that pass each cut
    int Num_surv_cut[Num_histo_sets] = {0};
    
    
    //Loop over events:
    //The ides is to keep each cut as a separate entity 
    //so they can be better handled
    for(int i = 0; i != nevents; ++i){//event loop

      WZtree->GetEntry(i);

      //an index to indicate current cut number
      int cut_index = 0;

      //cuts
      if(!PassTriggers_Cut()) continue;
      Tabulate_Me(Num_surv_cut,cut_index,weight);
      if(!HasValidWandZ_Cut()) continue;
      Tabulate_Me(Num_surv_cut,cut_index,weight);
      if(ExeedMaxNumberOfZs_Cut(cutMaxNumZs)) continue;
      Tabulate_Me(Num_surv_cut,cut_index,weight);


    }//event loop
    
    // total # of events (before any cuts)
    Nexp_evt += nevents * weight;
    
    //Number of expected events for each cut (weighted)
    for(int ii = 0; ii < Num_histo_sets; ++ii){
      if(debugme) cout<<"Num_surv_cut["<<ii<<"] = "<<
        Num_surv_cut[ii]<<endl;
      Nexp_evt_cut[ii] += Num_surv_cut[ii] * weight;
      if(debugme) cout<<"Nexp_evt_cut["<<ii<<"] = "<<
        Nexp_evt_cut[ii]<<endl;
    }

  }//loop over files
  
  printSummary(out, dir, Nexp_evt, Nexp_evt_cut);
  saveHistos(fout, dir);

}//Get_Distributions






//-----------------------------------------------------------
void ExecuteAnalysis()
{
//-----------------------------------------------------------


  if (debugme)cout<<"Master macro to execute analysis"<<endl;

  //value of lumi to be used in the analysis
  //the weights will scale accordingly.
  float lumiPb = 100000;
 
  //name of file where to write all histograms
  TFile *fout = new TFile("Wprime_analysis.root","recreate");
 
  //containers to
  //include signal and background files
  vector<InputFile> wzjj_files;
  vector<InputFile> wprime400_files;

  //load the harcoded crosssectios
  //add background or signal containers as needed
  Load_Cross_Sections(wzjj_files, wprime400_files);
 
  //keep account of events
  string outfile("event_counts.txt");
  ofstream out(outfile.c_str());
  if(!out) { 
    cout << "Cannot open file " << outfile << endl; 
    abort();
  } 
  
  //go for the analysis now and separate 
  //nicely into background/signal-type directories
  //the results will be written under respective directories
  //Add as many as you need:
  string dir = "wzjj";
  Load_Input_Files(dir, wzjj_files, lumiPb);
  Get_Distributions(wzjj_files, fout, dir, out);
  
  dir = "wprime400";
  Load_Input_Files(dir, wprime400_files, lumiPb);
  Get_Distributions(wprime400_files, fout, dir, out);

  out.close(); 
  fout->Close();

}//ExecuteAnalysis




