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
    if(Denom){
	eff = Num/Denom;
	deff = TMath::Sqrt(eff * (1-eff)/Denom);
    }else{
	eff = deff = 0;
    }

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
double deltaEta(double eta1, double eta2)
{
//--------------------------------------------------------------
    double eta = fabs(eta1 - eta2);
    return eta;
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
  hWZInvMass[3] = new TH1F("hWZInvMass_muon","Reconstructed WZ Invariant Mass (muon);m_{WZ} (GeV);",WZbin,WZinvMin,WZinvMax);
  hWZInvMass[4] = new TH1F("hWZInvMass_elec","Reconstructed WZ Invariant Mass (elec);m_{WZ} (GeV);",WZbin,WZinvMin,WZinvMax);
  hWZInvMass[5] = new TH1F("hWZInvMass_fwdjets","Reconstructed WZ Invariant Mass (fwdjets);m_{WZ} (GeV);",WZbin,WZinvMin,WZinvMax);
  hWZInvMass[6] = new TH1F("hWZInvMass_Zmass","Reconstructed WZ Invariant Mass (Zmass);m_{WZ} (GeV);",WZbin,WZinvMin,WZinvMax);
  hWZInvMass[7] = new TH1F("hWZInvMass_Ht","Reconstructed WZ Invariant Mass (Ht);m_{WZ} (GeV);",WZbin,WZinvMin,WZinvMax);


  //Jet Multiplicity Histos
  float jmMin = 0.;
  float jmMax = 15.;
  int jmbin = int(jmMax);
  hJetMult[0] = new TH1F("hJetMult_hlt","Jet Multiplicity (hlt);Num of Jets;",jmbin,jmMin,jmMax);
  hJetMult[1] = new TH1F("hJetMult_flavor","Jet Multiplicity (flavor);Num of Jets;",jmbin,jmMin,jmMax);
  hJetMult[2] = new TH1F("hJetMult_numZs","Jet Multiplicity (numZs);Num of Jets;",jmbin,jmMin,jmMax);
  hJetMult[3] = new TH1F("hJetMult_muon","Jet Multiplicity (muon);Num of Jets;",jmbin,jmMin,jmMax);
  hJetMult[4] = new TH1F("hJetMult_elec","Jet Multiplicity (elec);Num of Jets;",jmbin,jmMin,jmMax);
  hJetMult[5] = new TH1F("hJetMult_fwdjets","Jet Multiplicity (fwdjets);Num of Jets;",jmbin,jmMin,jmMax);
  hJetMult[6] = new TH1F("hJetMult_Zmass","Jet Multiplicity (Zmass);Num of Jets;",jmbin,jmMin,jmMax);
  hJetMult[7] = new TH1F("hJetMult_Ht","Jet Multiplicity (Ht);Num of Jets;",jmbin,jmMin,jmMax);


  //DeltaEta between jets
  float jdetaMin = 0;
  float jdetaMax = 10.;
  int jdetabin = 50;
  hJetsDeltaEta[0] = new TH1F("hJetsDeltaEta_hlt","Jets Delta Eta (hlt);|#Delta#eta_{jj}| (rad);",jdetabin,jdetaMin,jdetaMax);
  hJetsDeltaEta[1] = new TH1F("hJetsDeltaEta_flavor","Jets Delta Eta (flavor);|#Delta#eta_{jj}| (rad);",jdetabin,jdetaMin,jdetaMax);
  hJetsDeltaEta[2] = new TH1F("hJetsDeltaEta_numZs","Jets Delta Eta (numZs);|#Delta#eta_{jj}| (rad);",jdetabin,jdetaMin,jdetaMax);
  hJetsDeltaEta[3] = new TH1F("hJetsDeltaEta_muon","Jets Delta Eta (muon);|#Delta#eta_{jj}| (rad);",jdetabin,jdetaMin,jdetaMax);
  hJetsDeltaEta[4] = new TH1F("hJetsDeltaEta_elec","Jets Delta Eta (elec);|#Delta#eta_{jj}| (rad);",jdetabin,jdetaMin,jdetaMax);
  hJetsDeltaEta[5] = new TH1F("hJetsDeltaEta_fwdjets","Jets Delta Eta (fwdjets);|#Delta#eta_{jj}| (rad);",jdetabin,jdetaMin,jdetaMax);
  hJetsDeltaEta[6] = new TH1F("hJetsDeltaEta_Zmass","Jets Delta Eta (Zmass);|#Delta#eta_{jj}| (rad);",jdetabin,jdetaMin,jdetaMax);
  hJetsDeltaEta[7] = new TH1F("hJetsDeltaEta_Ht","Jets Delta Eta (Ht);|#Delta#eta_{jj}| (rad);",jdetabin,jdetaMin,jdetaMax);

  //Ht Histos
  float HtMin = 0.;
  float HtMax = 1000.;
  int Htbin = 100;
  hHt[0] = new TH1F("hHt_hlt","H_{T} (hlt);Lepton Pt Sum;",Htbin,HtMin,HtMax);
  hHt[1] = new TH1F("hHt_flavor","H_{T} (flavor);Lepton Pt Sum;",Htbin,HtMin,HtMax);
  hHt[2] = new TH1F("hHt_numZs","H_{T} (numZs);Lepton Pt Sum;",Htbin,HtMin,HtMax);
  hHt[3] = new TH1F("hHt_muon","H_{T} (muon);Lepton Pt Sum;",Htbin,HtMin,HtMax);
  hHt[4] = new TH1F("hHt_elec","H_{T} (elec);Lepton Pt Sum;",Htbin,HtMin,HtMax);
  hHt[5] = new TH1F("hHt_fwdjets","H_{T} (fwdjets);Lepton Pt Sum;",Htbin,HtMin,HtMax);
  hHt[6] = new TH1F("hHt_Zmass","H_{T} (Zmass);Lepton Pt Sum;",Htbin,HtMin,HtMax);
  hHt[7] = new TH1F("hHt_Ht","H_{T} (Ht);Lepton Pt Sum;",Htbin,HtMin,HtMax);

  //Wpt Histos
  float WptMin = 0.;
  float WptMax = 200.;
  int Wptbin = 20;
  hWpt[0] = new TH1F("hWpt_hlt","p_{T}(W) (hlt);p_{T} of W;",Wptbin,WptMin,WptMax);
  hWpt[1] = new TH1F("hWpt_flavor","p_{T}(W) (flavor);p_{T} of W;",Wptbin,WptMin,WptMax);
  hWpt[2] = new TH1F("hWpt_numZs","p_{T}(W) (numZs);p_{T} of W;",Wptbin,WptMin,WptMax);
  hWpt[3] = new TH1F("hWpt_muon","p_{T}(W) (muon);p_{T} of W;",Wptbin,WptMin,WptMax);
  hWpt[4] = new TH1F("hWpt_elec","p_{T}(W) (elec);p_{T} of W;",Wptbin,WptMin,WptMax);
  hWpt[5] = new TH1F("hWpt_fwdjets","p_{T}(W) (fwdjets);p_{T} of W;",Wptbin,WptMin,WptMax);
  hWpt[6] = new TH1F("hWpt_Zmass","p_{T}(W) (Zmass);p_{T} of W;",Wptbin,WptMin,WptMax);
  hWpt[7] = new TH1F("hWpt_Wpt","p_{T}(W) (Wpt);p_{T} of W;",Wptbin,WptMin,WptMax);


  //Zpt Histos
  float ZptMin = 0.;
  float ZptMax = 200.;
  int Zptbin = 20;
  hZpt[0] = new TH1F("hZpt_hlt","p_{T}(Z) (hlt);p_{T} of Z;",Zptbin,ZptMin,ZptMax);
  hZpt[1] = new TH1F("hZpt_flavor","p_{T}(Z) (flavor);p_{T} of Z;",Zptbin,ZptMin,ZptMax);
  hZpt[2] = new TH1F("hZpt_numZs","p_{T}(Z) (numZs);p_{T} of Z;",Zptbin,ZptMin,ZptMax);
  hZpt[3] = new TH1F("hZpt_muon","p_{T}(Z) (muon);p_{T} of Z;",Zptbin,ZptMin,ZptMax);
  hZpt[4] = new TH1F("hZpt_elec","p_{T}(Z) (elec);p_{T} of Z;",Zptbin,ZptMin,ZptMax);
  hZpt[5] = new TH1F("hZpt_fwdjets","p_{T}(Z) (fwdjets);p_{T} of Z;",Zptbin,ZptMin,ZptMax);
  hZpt[6] = new TH1F("hZpt_Zmass","p_{T}(Z) (Zmass);p_{T} of Z;",Zptbin,ZptMin,ZptMax);
  hZpt[7] = new TH1F("hZpt_Zpt","p_{T}(Z) (Zpt);p_{T} of Z;",Zptbin,ZptMin,ZptMax);

  hNumEvts = new TH1F("hNumEvts","Expected # of Events",Num_histo_sets,0,Num_histo_sets);
  hEffRel = new TH1F("hEffRel","Relative Efficiency",Num_histo_sets,0,Num_histo_sets);
  hEffAbs = new TH1F("hEffAbs","Absolute Efficiency",Num_histo_sets,0,Num_histo_sets);
  

}//Declare_Histos

void Declare_Lists()
{
//--------------------------------------------------------------

  if (debugme) cout<<"Declare lists"<<endl;

  cutlist[0] = new TEntryList("cut0","HLT");
  cutlist[1] = new TEntryList("cut1","flavor");
  cutlist[2] = new TEntryList("cut2","numZs");
  cutlist[3] = new TEntryList("cut3","muon");
  cutlist[4] = new TEntryList("cut4","elec");
  cutlist[5] = new TEntryList("cut5","fwdjets");
  cutlist[6] = new TEntryList("cut6","Zmass");
  cutlist[7] = new TEntryList("cut7","Ht");

}//Declare_Lists


//See if particle(only Electron??) is in Barrell
//--------------------------------------------------------------
bool inBarrel(float eta)
{
    return TMath::Abs(eta) < maxEtaBarrel;
  
}//InBarrel

bool inEndCap(float eta)
{
    float abs_eta = TMath::Abs(eta);
    return (abs_eta > minEtaEndcap && abs_eta < maxElecEta);
  
}//InEndCap


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

    if(Nfiles > 1){
        int filenum = filenum_low;
        for(int i = 0; i != Nfiles; ++i){
            filenum = filenum_low + i*filenum_step;
            files[i].pathname = top_level_dir + mask1 + 
                convertIntToStr(filenum) + mask2;
            files[i].description = file_desc;
        }
    }
    else if (Nfiles == 1){
        files[0].pathname = top_level_dir + mask1;
        files[0].description = file_desc;
    }
    else {cout<<"Something went wrong with the files. Quiting..."
              <<endl; abort();}
    



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
          Nfiles = 1;
          const string filename = "WZjj_mergedOutputTree.root";
          RecruitOrderedFiles(files,Nfiles,0,0,
                              filename,"",file_desc);
  }
  else if (!strcmp(file_desc.c_str(),"wprime400")){
          Nfiles = 1;
          const string filename = "Wprime400_mergedOutputTree.root";
          RecruitOrderedFiles(files,Nfiles,0,0,
                              filename,"",file_desc);
  }
  else cout<<"No samples were found with the name "<<file_desc<<endl;

  
  
  
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

  //qcd_files.push_back(InputFile(590000)); // 100_150
  //qcd_files.push_back(InputFile( 83000)); // 150_200
  //qcd_files.push_back(InputFile( 24000)); // 200_300
  //qcd_files.push_back(InputFile(  3000)); // 300_400
  //qcd_files.push_back(InputFile(   730)); // 400_600
  //qcd_files.push_back(InputFile(    66)); // 600_800
  //qcd_files.push_back(InputFile(    12)); // 800_1200
  //qcd_files.push_back(InputFile(    0.63)); // 1200_1600
  //qcd_files.push_back(InputFile(    0.54)); // 1600_2000

  wzjj_files.push_back(InputFile(xsec_wzjj)); // wzjj 
  wprime400_files.push_back(InputFile(xsec_wprime400)); // 400 GeV


 

  return 0;
}





//Set the branch addresses that we need for the analysis
//-----------------------------------------------------------
void Set_Branch_Addresses(TTree* WZtree)
{
//-----------------------------------------------------------
    if(debugme) cout<<"Setting Branch Addresses"<<endl;


    WZtree->SetBranchAddress("W_flavor",&W_flavor);
    WZtree->SetBranchAddress("Z_flavor",&Z_flavor);
    WZtree->SetBranchAddress("W_pt",&W_pt);
    WZtree->SetBranchAddress("Z_pt",&Z_pt);
//  WZtree->SetBranchAddress("triggerBitMask", &triggerBitMask);
    WZtree->SetBranchAddress("numberOfZs",&numberOfZs);
    WZtree->SetBranchAddress("WZ_invMassMinPz",&WZ_invMassMinPz);
    WZtree->SetBranchAddress("jet_energy",&jet_energy);
    WZtree->SetBranchAddress("jet_et",&jet_et);
    WZtree->SetBranchAddress("jet_eta",&jet_eta);
    WZtree->SetBranchAddress("jet_mass",&jet_mass);
    WZtree->SetBranchAddress("jet_phi",&jet_phi);
    WZtree->SetBranchAddress("jet_pt",&jet_pt);
    WZtree->SetBranchAddress("jet_px",&jet_px);
    WZtree->SetBranchAddress("jet_py",&jet_py);
    WZtree->SetBranchAddress("jet_pz",&jet_pz);
  
    WZtree->SetBranchAddress("muon_eta",&muon_eta);
    WZtree->SetBranchAddress("muon_pt",&muon_pt);
    WZtree->SetBranchAddress("electron_eta",&electron_eta);
    WZtree->SetBranchAddress("electron_pt",&electron_pt);
 
    WZtree->SetBranchAddress("W_leptonIndex",&W_leptonIndex);
    WZtree->SetBranchAddress("Z_leptonIndex1",&Z_leptonIndex1);
    WZtree->SetBranchAddress("Z_leptonIndex2",&Z_leptonIndex2);
  
    WZtree->SetBranchAddress("muon_innerD0",&muon_innerD0);
    WZtree->SetBranchAddress("muon_innerD0Error",&muon_innerD0Error);
    WZtree->SetBranchAddress("muon_caloIso",&muon_caloIso);
    WZtree->SetBranchAddress("muon_trackIso",&muon_trackIso);
 
    WZtree->SetBranchAddress("electron_sigmaEtaEta", &electron_sigmaEtaEta);
    WZtree->SetBranchAddress("electron_trackIso", &electron_trackIso);
    WZtree->SetBranchAddress("electron_caloIso", &electron_caloIso);
    WZtree->SetBranchAddress("electron_deltaEtaIn", &electron_deltaEtaIn);
    WZtree->SetBranchAddress("electron_deltaPhiIn", &electron_deltaPhiIn);
    WZtree->SetBranchAddress("electron_eOverP", &electron_eOverP);
    WZtree->SetBranchAddress("electron_hOverE", &electron_hOverE);
  
}//Set_Branch_Addresses







//Fill Histograms
//-----------------------------------------------------------
void Fill_Histos(int index, float weight)
{
//-----------------------------------------------------------

    
    hWZInvMass[index]->Fill(WZ_invMassMinPz, weight);
    hJetMult[index]->Fill(jet_energy->size(),weight);
    if(jet_eta->size() > 1){
        hJetsDeltaEta[index]->Fill(deltaEta(jet_eta->at(0),
                                            jet_eta->at(1)),weight);
    }
    hHt[index]->Fill(Calc_Ht(),weight);
    hZpt[index]->Fill(Z_pt,weight);
    hWpt[index]->Fill(W_pt,weight);

}//Fill_Histos





//------------------------------------------------------------------------
void saveHistos(TFile * fout, string dir)
{
//------------------------------------------------------------------------
  if (debugme) cout<<"Save Histos....."<<endl;
  fout->cd(); 
  fout->mkdir(dir.c_str()); 
  fout->cd(dir.c_str());

  for(int i = 0; i != Num_histo_sets; ++i){
        hWZInvMass[i]->Write();
        hJetMult[i]->Write();
        hJetsDeltaEta[i]->Write();
	hHt[i]->Write();
	hWpt[i]->Write();
	hZpt[i]->Write();
	cutlist[i]->Write();
  }
  hNumEvts->Write();
  hEffRel->Write();
  hEffAbs->Write();
  
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
	hNumEvts->Fill(i,Nexp_evt_cut[i]);

        //calculate efficiencies
        float eff, deff;
        if(i == 0){
            getEff(eff, deff, Nexp_evt_cut[i], Nexp_evt);
	    hEffRel->Fill(i,eff*100);	
	}else{
            getEff(eff, deff, Nexp_evt_cut[i], Nexp_evt_cut[i-1]);
	    hEffRel->Fill(i,eff*100);
	}
        out << ", Relative eff = "<<eff*100 << " +/- " << deff*100 << "%";
        
	getEff(eff, deff, Nexp_evt_cut[i], Nexp_evt);
        hEffAbs->Fill(i,eff*100);
	out << ", Absolute eff = "<< eff*100 << " +/- " << deff*100 << "%"
             << endl;
        
        //to do: put these results in a file
        
    } // loop over different cuts
    
    
}//printSummary






//Tabulate results after the cut has been passed
//-----------------------------------------------------------
void Tabulate_Me(int Num_surv_cut[], int& cut_index, 
		 const float& weight, const int& evtnum)
{
//-----------------------------------------------------------
    if(debugme) cout<<"Tabulating results for cut_index = "
                    <<cut_index<<endl;

    //increase the number of events passing the cuts
    ++Num_surv_cut[cut_index];
    //fill the histograms
    Fill_Histos(cut_index,weight);
    
    cutlist[cut_index]->Enter(evtnum);

    //since the event has passed the cut,
    //increase the cut_index for the next cut
    ++cut_index;
    if(debugme) cout<<"cut_index is now = "<<cut_index<<endl;

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
    bool has_valid_W_and_Z = (Z_flavor && W_flavor); //Cory: We can expand this (maybe in lepton cuts)
    return has_valid_W_and_Z;
    
}//--- NotValidWandZ_Cut


//Check if there is min number of leptons
//-----------------------------------------------------------
bool PassNumLeptons_Cut()
{
//-----------------------------------------------------------
    if(debugme) cout<<"Check if there is min # of leptons"<<endl;
    if((int)(muon_pt->size() + electron_pt->size()) < minNumLeptons) return false;

    return true;
}//--- PassNumLeptons_Cut

//Check if min # of leptons have min pt
//-----------------------------------------------------------
bool PassLeptonPt_Cut()
{
//-----------------------------------------------------------
    if(debugme) cout<<"Check min pt"<<endl;
    
    int count=0;
    int size = muon_pt->size();
    for(int i=0;i<size;++i){
	if(muon_pt->at(i) > minMuonPt) ++count;
	if(count > 3) return true;
    }

    size = electron_pt->size();
    for(int i=0;i<size;++i){
	if(electron_pt->at(i) > minElecPt) ++count;
	if(count > 3) return true;
    }

    return false;
}//--- PassLeptonPt_Cut


//Check if there is more Zs than required
//-----------------------------------------------------------
bool ExeedMaxNumberOfZs_Cut(const int& max_num_Zs)
{
//-----------------------------------------------------------
    if(debugme) cout<<"Check if there is more Zs than required"<<endl;
    bool has_too_many_Zs = numberOfZs > max_num_Zs;
    return has_too_many_Zs;;

}//--- MaxNumberOfZs_Cut

//-----------------------------------------------------------
bool HasTwoEnergeticForwardJets_Cut(const float& cutMinJetE,
                                    const float& cutMinJetPt,
                                    const float& cutMinJetsDeltaEta)
{
//-----------------------------------------------------------

    
    //Here we assume that the jet collections are ordered in energy
    //and in pT, which is true if the root-uple comes from PAT-uples
    bool has_fwdjets = false;
    int njets = jet_energy->size();
    if(njets < 2) return has_fwdjets;
    
    has_fwdjets = (jet_energy->at(0) > cutMinJetE &&
                   jet_energy->at(1) > cutMinJetE) &&
        (jet_pt->at(0) > cutMinJetPt &&
         jet_pt->at(1) > cutMinJetPt) &&
        deltaEta(jet_eta->at(0),jet_eta->at(1))>cutMinJetsDeltaEta;

    return has_fwdjets;


}//-----HasTwoEnergeticForwardJets

//Check Muon Properties
//-----------------------------------------------------------
bool PassMuonCut()
{
//-----------------------------------------------------------
    if(debugme) cout<<"Check Muon Cuts"<<endl;
    
    if(abs(Z_flavor) == PDGMUON ){
	if(TMath::Abs(muon_eta->at(Z_leptonIndex1)) > maxMuonEta) return false;
	if(TMath::Abs(muon_eta->at(Z_leptonIndex2)) > maxMuonEta) return false;
	
	if(muon_pt->at(Z_leptonIndex1) < cutZMuon_pt) return false;
	if(muon_pt->at(Z_leptonIndex2) < cutZMuon_pt) return false;
    }else if (abs(W_flavor) == PDGMUON){
	if(TMath::Abs(muon_eta->at(W_leptonIndex)) > maxMuonEta) return false;
	if(muon_pt->at(W_leptonIndex) < cutWMuon_pt) return false;
	if(!PassMuonSip(W_leptonIndex)) return false;
	if(!PassMuonRelIso(W_leptonIndex)) return false;
    }
    
    return true;

}//--- Muon Properties Cut

//Check Muon Sip
//-----------------------------------------------------------
bool PassMuonSip(int index)
{
    float Sip = muon_innerD0->at(index) / muon_innerD0Error->at(index);
    return Sip < cutWmuD0;  //Cory: inner or global D0
}//--- Muon Sip Cut

//Check Muon RelIso
//-----------------------------------------------------------
bool PassMuonRelIso(int index)
{
    float relIso = muon_caloIso->at(index) + muon_trackIso->at(index);
    relIso /= muon_pt->at(index);

    return relIso < cutWmuCombRelIso;
}//--- Muon RelIso Cut


//Check Electron Properties
//-----------------------------------------------------------
bool PassElecCut()
{
//-----------------------------------------------------------
    if(debugme) cout<<"Check Electron Cuts"<<endl;
    int index;
    if(abs(Z_flavor) == PDGELEC ){
	for(int i=0; i<2; ++i){
	    if(i==0) index = Z_leptonIndex1;
	    else     index = Z_leptonIndex2;
	    
	    if(!PassElecEtaDepCut(index,PDGZ)) return false;
	}
    }else if( abs(W_flavor) == PDGELEC ){
	index = W_leptonIndex;
	if(!PassElecEtaDepCut(index,PDGW)) return false;
    }
    
    return true;
    
}//--- PassElecCut


//Check Electron Eta Dependent Cuts
//-----------------------------------------------------------
bool PassElecEtaDepCut(int index, int parent)
{
//-----------------------------------------------------------
    if(debugme) cout<<"Check Electron Eta Dep Cuts"<<endl;

    float eta = electron_eta->at(index);

    float pt = electron_pt->at(index);
    float deta = electron_deltaEtaIn->at(index);
    float dphi = electron_deltaPhiIn->at(index);
    float sigmann = electron_sigmaEtaEta->at(index);
    float EP = electron_eOverP->at(index);
    float hE = electron_hOverE->at(index);
    float calIso = electron_caloIso->at(index);
    float trkIso = electron_trackIso->at(index);
	  
    if(!PassElecPtCut(pt,eta,parent)) return false;
    if(!PassElecdEtaCut(deta,eta)) return false;
    if(!PassElecdPhiCut(dphi,eta)) return false;
    if(!PassElecSigmannCut(sigmann,eta)) return false;
    if(!PassElecEPCut(EP,eta)) return false;
    if(!PassElecHECut(hE,eta)) return false;
    if(!PassElecCalIsoCut(calIso,eta)) return false;
    if(!PassElecTrkIsoCut(trkIso,eta)) return false;
	
    return true;

}//--- PassElecPtCut

//Check Electron Pt Cut
//-----------------------------------------------------------
bool PassElecPtCut(float pt,float eta,float parent)
{
//-----------------------------------------------------------
    if(debugme) cout<<"Check Electron Eta Pt Cut"<<endl;
    
    if(inBarrel(eta)){
	if     (parent == PDGW && pt > cutWElectron_pt) return true;
	else if(parent == PDGZ && pt > cutZElectron_pt) return true;
    }else if(inEndCap(eta)){
	if     (parent == PDGW && pt > cutWElectron_pt) return true;
	else if(parent == PDGZ && pt > cutZElectron_pt) return true;
    }

   return false;

}//--- PassElecPtCut

//Check Electron dEta Cut
//-----------------------------------------------------------
bool PassElecdEtaCut(float dEta, float eta)
{
//-----------------------------------------------------------
    if(debugme) cout<<"Check Electron Eta dEta Cut"<<endl;

    if(inBarrel(eta)){
	if(dEta < cutDeltaEtaIn[0]) return true;
    }else if(inEndCap(eta)){
	if(dEta < cutDeltaEtaIn[1]) return true;
    }

   return false;

}//--- PassElecdEtaCut

//Check Electron dPhi Cut
//-----------------------------------------------------------
bool PassElecdPhiCut(float dPhi, float eta)
{
//-----------------------------------------------------------
    if(debugme) cout<<"Check Electron Eta dPhi Cut"<<endl;

    if(inBarrel(eta)){
	if(dPhi < cutDeltaPhiIn[0]) return true;
    }else if(inEndCap(eta)){
	if(dPhi < cutDeltaPhiIn[1]) return true;
    }

   return false;

}//--- PassElecdPhiCut

//Check Electron sigmann Cut
//-----------------------------------------------------------
bool PassElecSigmannCut(float sigmann, float eta)
{
//-----------------------------------------------------------
    if(debugme) cout<<"Check Electron Eta sigmann Cut"<<endl;

    if(inBarrel(eta)){
	if(sigmann < cutSigmaEtaEta[0]) return true;
    }else if(inEndCap(eta)){
	if(sigmann < cutSigmaEtaEta[1]) return true;
    }

    return false;

}//--- PassElecsigmannCut

//Check Electron EP Cut
//-----------------------------------------------------------
bool PassElecEPCut(float EP,float eta)
{
//-----------------------------------------------------------
    if(debugme) cout<<"Check Electron Eta EP Cut"<<endl;

    if(inBarrel(eta)){
	if(EP > cutEOverP[0]) return true;
    }else if(inEndCap(eta)){
	if(EP > cutEOverP[1]) return true;
    }

   return false;

}//--- PassElecEPCut

//Check Electron hE Cut
//-----------------------------------------------------------
bool PassElecHECut(float HE, float eta)
{
//-----------------------------------------------------------
    if(debugme) cout<<"Check Electron Eta hE Cut"<<endl;

    if(inBarrel(eta)){
	if(HE < cutHOverE[0]) return true;
    }else if(inEndCap(eta)){
	if(HE < cutHOverE[1]) return true;
    }

   return false;

}//--- PassElechECut

//Check Electron calIso Cut
//-----------------------------------------------------------
bool PassElecCalIsoCut(float calIso, float eta)
{
//-----------------------------------------------------------
    if(debugme) cout<<"Check Electron Eta calIso Cut"<<endl;

    if(inBarrel(eta)){
	if(calIso < cutCalRelIso[0]) return true;
    }else if(inEndCap(eta)){
	if(calIso < cutCalRelIso[1]) return true;
    }

   return false;

}//--- PassEleccalIsoCut

//Check Electron trkIso Cut
//-----------------------------------------------------------
bool PassElecTrkIsoCut(float trkIso,float eta)
{
//-----------------------------------------------------------
    if(debugme) cout<<"Check Electron Eta trkIso Cut"<<endl;

    if(inBarrel(eta)){
	if(trkIso < cutTrackRelIso[0]) return true;
    }else if(inEndCap(eta)){
	if(trkIso < cutTrackRelIso[1]) return true;
    }

   return false;

}//--- PassElectrkIsoCut


//Check Z Properties
//-----------------------------------------------------------
bool PassZCut()
{
//-----------------------------------------------------------
    if(debugme) cout<<"Check Z Cuts"<<endl;
    if(TMath::Abs(Z_mass - PDGZMASS) > ZMASSRES) return false; //GeV
    
    return true;

}//--- PassZCut

//Check Ht Properties
//-----------------------------------------------------------
bool PassHtCut()
{
//-----------------------------------------------------------
    if(debugme) cout<<"Check Ht Cuts"<<endl;
   
    return Calc_Ht() > minHt;
    
}//--- PassHtCut

//Calc Ht
//-----------------------------------------------------------
float Calc_Ht()
{
    float Ht=0;
  
    if     (W_flavor == PDGELEC) Ht += electron_pt->at(W_leptonIndex);
    else if(W_flavor == PDGMUON) Ht += muon_pt->at(W_leptonIndex);
    //else                         cout<<"W didn't decay into e,mu"<<endl;

    if     (Z_flavor == PDGELEC){
	Ht += electron_pt->at(Z_leptonIndex1);
	Ht += electron_pt->at(Z_leptonIndex2);
    }else if(Z_flavor == PDGMUON){
	Ht += muon_pt->at(Z_leptonIndex1);
	Ht += muon_pt->at(Z_leptonIndex2);
    }//else                         cout<<"Z didn't decay into e,mu"<<endl;

    return Ht;

}//--- CalcHt

//Check Ht Met Properties
//-----------------------------------------------------------
bool PassHtMetCut()
{
//-----------------------------------------------------------
    if(debugme) cout<<"Check Ht Cuts"<<endl;
    float Ht=0;
  
    if     (W_flavor == PDGELEC) Ht += electron_pt->at(W_leptonIndex);
    else if(W_flavor == PDGMUON) Ht += muon_pt->at(W_leptonIndex);
    else                         cout<<"W didn't decay into e,mu"<<endl;

    if     (Z_flavor == PDGELEC){
	Ht += electron_pt->at(Z_leptonIndex1);
	Ht += electron_pt->at(Z_leptonIndex2);
    }else if(Z_flavor == PDGMUON){
	Ht += muon_pt->at(Z_leptonIndex1);
	Ht += muon_pt->at(Z_leptonIndex2);
    }else                         cout<<"Z didn't decay into e,mu"<<endl;

    Ht += met_et;
    
    return Ht > minHtMet;
    
}//--- PassHtMetCut




//Get different types of distribution
//-----------------------------------------------------------
void Get_Distributions(vector<InputFile>& files, 
                       TFile *fout, string dir, ofstream & out)
{
//-----------------------------------------------------------
  if (debugme) cout<<"Get Distributions....."<<endl;
  
  
  Declare_Histos();
  Declare_Lists();

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

    cout << "Processing file "<<files[tr].pathname<<endl;
    

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
      if (debugme) cout<<"Processing event "<<i<<endl;
      //an index to indicate current cut number
      int cut_index = 0;

      //cuts
      if(!PassTriggers_Cut()) continue;
      Tabulate_Me(Num_surv_cut,cut_index,weight,i);
      if(!HasValidWandZ_Cut()) continue;
      Tabulate_Me(Num_surv_cut,cut_index,weight,i);
      if(ExeedMaxNumberOfZs_Cut(cutMaxNumZs)) continue;
      Tabulate_Me(Num_surv_cut,cut_index,weight,i);
      if(!PassMuonCut()) continue;
      Tabulate_Me(Num_surv_cut,cut_index,weight,i);
      if(!PassElecCut()) continue;
      Tabulate_Me(Num_surv_cut,cut_index,weight,i);
      if(!HasTwoEnergeticForwardJets_Cut(cutMinJetE,
                                         cutMinJetPt,
                                         cutMinJetsDeltaEta)) continue;
      Tabulate_Me(Num_surv_cut,cut_index,weight,i);
      
      if(!PassZCut()) continue;
      Tabulate_Me(Num_surv_cut,cut_index,weight,i);
      if(!PassHtCut()) continue;
      Tabulate_Me(Num_surv_cut,cut_index,weight,i);
      
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
  float lumiPb = 20000;
 
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




