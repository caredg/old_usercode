#include <ExecuteAnalysis.h>

//--------------------------------------------------------------
void getEff(float & eff, float & deff, float Num, float Denom)
{
//--------------------------------------------------------------
  eff = Num/Denom;
  deff = TMath::Sqrt(eff * (1-eff)/Denom);

}//getEff


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



//load the input files from the top_level_dir
//-----------------------------------------------------------
void Load_Input_Files(string file_desc, 
                    vector<InputFile> & files, 
                    float lumiPb)
{
//-----------------------------------------------------------
  if (debugme)cout<<"Loading input Files...."<<endl;

  int Nfiles = -1;
  cout << "\n Processing " << file_desc << " files " << endl << endl;

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

  if(file_desc == "wzjj"){
    const int Nfileswzjj = 20;
    Nfiles = Nfileswzjj;        
    if(Check_Files(Nfiles, files))
      return;
    
    //FIXME: do this in a for loop 
    string filenum[Nfileswzjj] = {"201","202","203","204","205","206","207","208","209","210","211","212","213","214","215","216","217","218", "219", "220"};
    
    
    for(int i = 0; i != Nfiles; ++i){
	  files[i].pathname = top_level_dir + string("WZjj_500event_")+filenum[i]+string("_outputTree.root");
      if (debugme)cout<<"Loading: "<<files[i].pathname<<endl;

	}
  } // wzjj


  else if(file_desc == "wprime400"){
    const int NfilesWprime = 20;
    Nfiles = NfilesWprime;
    if(Check_Files(Nfiles, files))
      return;

    //FIXME: do this in a for loop 
    string filenumWp[NfilesWprime] = {"101","102","103","104","105","106","107","108","109","110","111","112","113","114","115","116","117","118", "119", "120"};

    for(int i = 0; i != NfilesWprime; ++i){
        files[i].pathname = top_level_dir + string("Wprime400_500event_")+filenumWp[i]+string("_outputTree.root");
      if (debugme)cout<<"Loading: "<<files[i].pathname<<endl;
    }
    
  } // Wprime 400


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
//      if(! files[i].tree->GetBranch("blah")){
//        cerr << " *** Can't find wp branch in file: " << pathname<< endl;
//        files[i].tree = 0;
//        continue;
//      }
    TH1F* histNumEvents       = (TH1F  *)file->Get("numEvents");
    const float eventsAnalyzed  = histNumEvents->GetBinContent(1);
    delete histNumEvents;
        
    files[i].Nprod_evt = eventsAnalyzed;
    files[i].weight = lumiPb*(files[i].x_sect)/(files[i].Nprod_evt);


      cout <<"Events produced in file "<<files[i].Nprod_evt<<",  # of entries = "<< files[i].tree->GetEntries() << ", weight = " << files[i].weight << endl;
    
  } // loop over input files
  
  return;
}





//---------------------------------------------------------
int Load_Cross_Sections(vector<InputFile> & wzjj_files, 
                      vector<InputFile> & wprime400_files)
{
//---------------------------------------------------------
  if (debugme)cout<<"Loading cross sections...."<<endl;


  //PYTHIA cross-section (pb)

  //qcd_files.push_back(wprime::InputFile(590000  , 450000)); // 100_150
  //qcd_files.push_back(wprime::InputFile( 83000  , 425000)); // 150_200
  //qcd_files.push_back(wprime::InputFile( 24000  , 440000)); // 200_300
  //qcd_files.push_back(wprime::InputFile(  3000  , 395000)); // 300_400
  //qcd_files.push_back(wprime::InputFile(   730  , 135000)); // 400_600
  //qcd_files.push_back(wprime::InputFile(    66  , 310000)); // 600_800
  //qcd_files.push_back(wprime::InputFile(    12  , 130000)); // 800_1200
  //qcd_files.push_back(wprime::InputFile(    0.63,  30000)); // 1200_1600
  //qcd_files.push_back(wprime::InputFile(    0.54,  25000)); // 1600_2000

  for (int j = 0; j<20;++j){
    wzjj_files.push_back(InputFile(450)); // wzjj 
    wprime400_files.push_back(InputFile(1.5)); // 400 GeV
  }
  

  return 0;
}


//Set the branch addresses that we need for the analysis
//-----------------------------------------------------------
void Set_Branch_Addresses(TTree* WZtree){

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


//--------------------------------------------------------------------------
void printSummary(ofstream & out, const string& dir, 
                  const float& Nexp_evt, float Nexp_evt_cut[]) 
{ 
//------------------------------------------------------------------------

 cout << " Total # of expected events = " << Nexp_evt << endl;
 
 for(int i = 0; i < Num_histo_sets; ++i){
   
   cout <<"Cut # "<<i<<": expected evts = " << Nexp_evt_cut[i];

   //calculate efficiencies
   float eff, deff;
   if(i == 0)
     getEff(eff, deff, Nexp_evt_cut[i], Nexp_evt);
   else
     getEff(eff, deff, Nexp_evt_cut[i], Nexp_evt_cut[i-1]);
   cout << ", Relative eff = "<<eff*100 << " +- " << deff*100 << "%";
   getEff(eff, deff, Nexp_evt_cut[i], Nexp_evt);
   cout << ", Absolute eff = "<< eff*100 << " +- " << deff*100 << "%"
        << endl;
   
   //to do: put these results in a file

 } // loop over different cuts

 
}//printSummary



//Get different types of distribution
//-----------------------------------------------------------
void Get_Distributions(vector<InputFile>& files, 
                       TFile *fout, string dir, ofstream & out)
{
//-----------------------------------------------------------
  if (debugme) cout<<"Get Distributions....."<<endl;
  
  
  Declare_Histos();
 
  int Nfiles = files.size();

  //initialize counters for expected number of events
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
    
    

    for(int i = 0; i != nevents; ++i){//event loop

      WZtree->GetEntry(i);

      int index = 0;
      //do trigger here, for now pass them all
      ++Num_surv_cut[index]; 
      Fill_Histos(index,weight);

      ++index;
      //// Reject events without a valid W and Z      
      if (!Z_flavor || !W_flavor) continue;
      ++Num_surv_cut[index];
      Fill_Histos(index,weight);

      ++index;
      // Reject events with more than 1 Z bosons
      if (numberOfZs > 1) continue;
      ++Num_surv_cut[index];
      Fill_Histos(index,weight);


    }//event loop
    
    Nexp_evt += nevents * weight; // total # of events (before any cuts)
    
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
  float lumiPb = 1000;
 
  //name of file where to write all histograms
  TFile *fout = new TFile("Wprime_analysis.root","recreate");
 
  //types of files
  //include signal and background files
  vector<InputFile> wzjj_files;
  vector<InputFile> wprime400_files;

  //load the harcoded crosssectios
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
  string dir = "wzjj";
  Load_Input_Files(dir, wzjj_files, lumiPb);
  Get_Distributions(wzjj_files, fout, dir, out);
  
  dir = "wprime400";
  Load_Input_Files(dir, wprime400_files, lumiPb);
  Get_Distributions(wprime400_files, fout, dir, out);

  out.close(); 
  fout->Close();

}//ExecuteAnalysis




