// ROOT49 stuff
#include <TROOT.h>
#include <TClonesArray.h>
#include <TFile.h>
#include <TTree.h>
#include <TMath.h>
#include <TH1F.h>
#include <stdio.h>
#include <algorithm>
#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <assert.h>
#include <vector>

// wprime stuff
#include "UserCode/CMGWPrimeGroup/interface/wprimeEvent.h"
#include "UserCode/CMGWPrimeGroup/root_macros/input_file.h"
#include "UserCode/CMGWPrimeGroup/root_macros/loadCutsAndThresholds.h"

using std::cout; using std::endl; using std::vector; using std::string;


//debug
const bool Debugme = false;

//to control the fixed parameters:
const int myNum_trkAlgos = Num_trkAlgos; // global, tracker, tev_1st
const float myOneMuPtTrackCut = OneMuPtTrackCut; 
//const float myOneMuPtTrackCut = 300; 
const float mySumPtCut = SumPtCut; // Cone DeltaR =0.3; 
const float myCombRelCut = CombRelCut;
const unsigned mydeltaRIsoIndex = deltaRIsoIndex; //for the isolation container
const float myEtJetCut = EtJetCut;
const float myDelta_Phi = Delta_Phi;//min angle muon/jet for jet veto
const unsigned myMaxNjetsAboveThresh = MaxNjetsAboveThresh;
const float myPtTrackCut = PtTrackCut;
const float myChi2Cut = Chi2Cut;
const float myminPtMu = 100;
const float mymaxPtMu = 800;
const float mynBinPtMu = 21;




//optimization constants
//if not controlled here they are controlled in the MakeOptim.C macro.
const string optimdir = "forOptimization";
const int o_maxNumCuts = 10000;//max number of cuts
const int o_NiDR = 9;
const float o_ValDeltaRIso[o_NiDR] = {0.2,0.25,0.3,0.35,0.4,0.45,0.5,0.55,0.6};
const unsigned o_indexDeltaR[o_NiDR] = {0,1,2,3,4,5,6,7,8};
const int o_NDPhiJ = 5;
const float o_DphiJetVeto[o_NDPhiJ] = { TMath::Pi() - 0.1, 
                                       TMath::Pi() - 0.2, 
                                       TMath::Pi() - 0.3, 
                                       TMath::Pi() - 0.4,
                                       TMath::Pi() - 0.5};

//data struct for histograms
struct data{
    string algo;
    vector<TH1F*> h1f; 
  data(string ALGO, vector<TH1F*> H1F):
    algo(ALGO), h1f(H1F){}
};


//methods in this macro:
void populateOptimFile (TFile* fout,const string& dir,const string& outfile, float optNbefore[],
                        float optNafter[][myNum_trkAlgos], TH1F* hbef[],
                        const vector<data>& haft,
                        const int option, float cutVals1[],
                        const int cv1size, float cutVals2[],
                        const int cv2size);
void printSummary_Optim(TFile* fout, const string& dir,
                        float  optNbefore[], 
                        float optNafter[][myNum_trkAlgos],
                        TH1F* hbef[],
                        const vector<data>& haft,
                        const int option, float cutVals1[],
			 const int cv1size, float cutVals2[],
                        const int cv2size);
void fillHistos_MuonPt(const wprime::Muon* mu, TH1F* h,
                       const float weight,const string& algo_name);
void gatherFileBasicInfo(const wprime::InputFile& file,
                         wprime::Event*& ev, int& nevents, float& weight);
void Optimize_OnlyOneMuon(const wprime::Event* ev, 
                          wprime::Muon*& the_mu,
                          int countsAfter[], 
                          vector<TH1F*>& hvec,float o_muonTrackPt[],
                          const int o_NmuTrkPt,
                          const float& weight,
                          const string& algoName);
void Optimize_Isolation(const wprime::Muon* the_mu,
                        int countsAfter[],vector<TH1F*>& hvec,
                        float o_pTIso[],
                        const int o_NpTIso,
                        const float& weight,
                        const string& algoName);
void Optimize_JetVeto(const wprime::Event* ev, 
                      const wprime::Muon* the_mu,
                      int countsAfter[],vector<TH1F*>& hvec,
                      float o_ETJetVeto[],
                      const int o_NETjet,
                      const float& weight,
                      const string& algoName);
void Optimize_Qual(const wprime::Muon* theMu,const string& algoName,
		   int countsAfter[],vector<TH1F*>& hvec,
                   float o_pTQual[],
		   const int o_NpTQual,float o_Chi2Qual[],
                   const int o_NChi2Qual,const float& weight);

void FillCutVals(float cutVals[],const float cutMax, const float cutMin,
		 const int cutsNumSteps);
void  FillHistVectors(vector<TH1F*>& vh,const int option,
                      float cutVals1[],float cutVals2[], 
                      const int cv1size, const int cv2size, const string& algoName);
void ExecOptimization(const wprime::InputFile& file,const string& dir,
                      const bool& loopMuons,
                      float optNbefore[],
                      float optNafter[][myNum_trkAlgos], 
                      TH1F* hbef[],vector<data>& haft,
                      const int option, float cutVals1[],
                      const int cv1size, float cutVals2[],
                      const int cv2size);





// returns # of (global) muons with pt above <tracker_muon_pt>, it is algo
// dependent and it delivers the pointer to the hardest muon
//------------------------------------------------------------------------
unsigned NmuAboveThresh(float tracker_muon_pt, 
                        const wprime::Event * ev,
                        wprime::Muon*& theMu,const string& algo_name)
{
//------------------------------------------------------------------------
  unsigned N = 0;

  int nmuons = ev->mu->GetLast() + 1;
  float temp_muPt = 0;
  for(int j = 0; j != nmuons; ++j)
    { // loop over muons
      wprime::Muon * mu = (wprime::Muon *) ev->mu->At(j);
      float current_muPt = 0;
      if (!strcmp(algo_name.c_str(),"gbl")){
          if (mu->global.p.Pt() > tracker_muon_pt){
              current_muPt = mu->global.p.Pt();
              if (current_muPt > temp_muPt){
                  temp_muPt = current_muPt;
                  theMu = mu;
              }
              ++N;
          }
      }
      else if (!strcmp(algo_name.c_str(),"trk")){
          if (mu->tracker.p.Pt() > tracker_muon_pt){
             current_muPt = mu->global.p.Pt();
              if (current_muPt > temp_muPt){
                  temp_muPt = current_muPt;
                  theMu = mu;
              }
              ++N;
          }
      }
      else if (!strcmp(algo_name.c_str(),"tev")){
          if (mu->tev_1st.p.Pt() > tracker_muon_pt){
            current_muPt = mu->global.p.Pt();
              if (current_muPt > temp_muPt){
                  temp_muPt = current_muPt;
                  theMu = mu;
              }
              ++N;
          }
          
      }
      else cout<<"NumAboveThresh: WARNING!!!!!! Algo name not recognized"<<endl;

      
    } // loop over muons
  
  return N;
}//------------NumAboveThresh



//-------------------------------------------------------------------------
void populateOptimFile (TFile* fout,const string& dir,
                        const string& outfile, float optNbefore[],
                        float optNafter[][myNum_trkAlgos], TH1F* hbef[],
                        const vector<data>& haft,
                        const int option, float cutVals1[],
                        const int cv1size, float cutVals2[],
                        const int cv2size)
{
//-------------------------------------------------------------------------
  if(Debugme) cout<<"Populating files for Optimization "<<endl;
  
  ofstream out(outfile.c_str());
  if(!out) {cout << "Cannot open file " << outfile << endl; abort();}
  
  //prepare root file for histograms
  fout->cd(); 
  fout->mkdir(dir.c_str()); 
  fout->cd(dir.c_str());


  if (option == 100){
      for(int j = 0; j < myNum_trkAlgos ; ++j){
          hbef[j]->Write();
          for (int k = 0; k < cv1size; ++k){
              haft.at(j).h1f.at(k)->Write();
              out<<algo_desc_short[j]<<"\t"<<cutVals1[k]<<"\t"<<-9999<<"\t"
                 <<-9999<<"\t"<<-9999<<"\t"<<optNbefore[j]
                 <<"\t"<<optNafter[k][j]<<endl;
          }//trkPt
      }//trkAlgos
  }//option 100
  else if (option == 200 || option == 500){
      for(int j = 0; j < myNum_trkAlgos ; ++j){
          hbef[j]->Write();
          int mycounter = 0;
          for(int ii = 0; ii < o_NiDR; ++ii){
              for(int jj = 0; jj < cv1size; ++jj){
                  haft.at(j).h1f.at(mycounter)->Write();
                  out<<algo_desc_short[j]<<"\t"<<cutVals1[jj]<<"\t"
                     <<o_ValDeltaRIso[o_indexDeltaR[ii]]<<"\t"
                     <<-9999<<"\t"
                     <<-9999<<"\t"<<optNbefore[j]<<"\t"<<
                      optNafter[mycounter][j]<<endl;
                  ++mycounter;
              }//NpTIso
          }//NiDR
      }//trkAlgos
  }//option 200 or 500
  else if (option == 300 || option == 400){
      for(int j = 0; j < myNum_trkAlgos ; ++j){
          int mycounter = 0;
          for(int ll = 0; ll < o_NDPhiJ; ++ll){
              for(int kk = 0; kk < cv1size; ++kk){
                  haft.at(j).h1f.at(mycounter)->Write();
                  out<<algo_desc_short[j]<<"\t"<<cutVals1[kk]<<"\t"<<
                      o_DphiJetVeto[ll]<<"\t"<<
                      -9999<<"\t"
                     <<-9999<<"\t"<<optNbefore[j]<<"\t"<<
                      optNafter[mycounter][j]<<endl;
                  ++mycounter;
              }//NDPhiJ
          }//NETjet
      }//trkAlgos
  }//option 300 or 400
  else if (option == 600 || option == 700){
      for(int j = 0; j < myNum_trkAlgos ; ++j){
          int mycounter = 0;
          for(int ll = 0; ll < cv2size; ++ll){
              for(int kk = 0; kk < cv1size; ++kk){
                  haft.at(j).h1f.at(mycounter)->Write();
                  out<<algo_desc_short[j]<<"\t"<<cutVals1[kk]<<"\t"<<
                      cutVals2[ll]<<"\t"<<
                      -9999<<"\t"
                     <<-9999<<"\t"<<optNbefore[j]<<"\t"<<
                      optNafter[mycounter][j]<<endl;
                  ++mycounter;
              }//NDPhiJ
          }//NETjet
      }//trkAlgos
  }//option 600 or 700
}






//--------------------------------------------------------------------------
void printSummary_Optim(TFile* fout, const string& dir,
                        float  optNbefore[], 
                        float optNafter[][myNum_trkAlgos],
                        TH1F* hbef[],
                        const vector<data>& haft,
                        const int option, float cutVals1[],
			 const int cv1size, float cutVals2[],
			 const int cv2size)
{
//------------------------------------------------------------------------
  if(Debugme) cout<<"Printing summary for Optimization "<<endl;

  string outdir = optimdir;
    string cutType = "";
    string outfile = "";

    if (option == 100){cutType = cuts_desc_short[1];}//1mu
    else if (option == 200){
      cutType = cuts_desc_short[1]+cuts_desc_short[3];//1muiso
    }
    else if (option == 300){
      cutType = cuts_desc_short[1]+cuts_desc_short[4];//1mujet
    }
    else if (option == 400){
      cutType = cuts_desc_short[1]+cuts_desc_short[3]+cuts_desc_short[4];//1muisojet
    }
    else if (option == 500){
      cutType = cuts_desc_short[1]+cuts_desc_short[4]+cuts_desc_short[3];//1mujetiso
    }
    else if (option == 600){
      cutType = cuts_desc_short[1]+cuts_desc_short[3]+cuts_desc_short[4]+cuts_desc_short[5];//1mujetisoqual
    }
    else if (option == 700){
      cutType = cuts_desc_short[1]+cuts_desc_short[4]+cuts_desc_short[3]+cuts_desc_short[5];//1mujetisoqual
    }
    else {cout<<"Nothing to print for this optim study"<<endl; return;}

    outfile = outdir +"/"+ dir + "_" + cutType +".dat";
    populateOptimFile(fout,dir,outfile, optNbefore, optNafter, hbef,haft,option,
		      cutVals1,cv1size,cutVals2,cv2size);

}//--------printSummary_Optim





//Fill Histograms
//-----------------------------------------------------------
void fillHistos_MuonPt(const wprime::Muon* mu, TH1F* h,
                       const float weight,const string& algo_name)
{
//-----------------------------------------------------------

    if(!mu) {
        if(Debugme) cout<<"No muon found in the event"<<endl;
        return;
    }
    

    if (!strcmp(algo_name.c_str(),"gbl")){
        h->Fill(mu->global.p.Pt(), weight);
    }
    else if (!strcmp(algo_name.c_str(),"trk")){
        h->Fill(mu->tracker.p.Pt(), weight);
    }
    else if (!strcmp(algo_name.c_str(),"tev")){
        h->Fill(mu->tev_1st.p.Pt(), weight);
    }
    else cout<<"WARNING!!!!!! Algo name not recognized"<<endl;
  

}//fillHistos_MuonPt








//Routine to grab the information from a given file
//---------------------------------------------------------------------------
void gatherFileBasicInfo(const wprime::InputFile& file,
                         wprime::Event*& ev, int& nevents, float& weight)
{
//---------------------------------------------------------------------------

    ev = new wprime::Event();
    file.tree->SetBranchAddress("wp", &ev);
    nevents = file.tree->GetEntries();
    if(Debugme) cout<<"Number of events in the file = "<<nevents<<endl;
    weight = file.weight;

}//---------gatherFileBasicInfo()




//Optimze OnlyOneMuon cut 
//---------------------------------------------------------------------------
void Optimize_OnlyOneMuon(const wprime::Event* ev, 
                          wprime::Muon*& the_mu,
                          int countsAfter[], 
                          vector<TH1F*>& hvec,float o_muonTrackPt[],
                          const int o_NmuTrkPt,
                          const float& weight,
                          const string& algoName)
{
//---------------------------------------------------------------------------
  if(Debugme) cout<<"Running Optimize_OnlyOneMuon"<<endl;
    //Loop over the colection of thresholds
    for (int nn = 0; nn < o_NmuTrkPt; ++nn){
        float pttrkcut = o_muonTrackPt[nn];
        //if(!OnlyOneHighTrackPtMuon(ev,the_mu,pttrkcut)) continue;
        if(!OnlyOneHighTrackPtMuon(ev,pttrkcut)) continue;
        ++countsAfter[nn];
        //cout<<"countsAfter["<<nn<<"] = "<<countsAfter[nn]<<endl;
        fillHistos_MuonPt(the_mu,hvec.at(nn),weight,algoName);
    }

    return;

}//-----Optimize_OnlyOneMuon()





//---------------------------------------------------------------------------
void Optimize_Isolation(const wprime::Muon* the_mu,
                        int countsAfter[],vector<TH1F*>& hvec,
                        float o_pTIso[],
                        const int o_NpTIso,
                        const float& weight,
                        const string& algoName)

{
//---------------------------------------------------------------------------
  if(Debugme) cout<<"Running Optimize_Isolation"<<endl;
  int mycounter = 0;
  for(int ii = 0; ii < o_NiDR; ++ii){
    for(int jj = 0; jj < o_NpTIso; ++jj){
	  unsigned dRidx = o_indexDeltaR[ii];
	  float sptcut = o_pTIso[jj];
	  //if(!SumPtIsolation(the_mu,dRidx,sptcut)) {++mycounter; continue;}
      if(!CombRelIsolation(the_mu,dRidx,sptcut)) {++mycounter; continue;}
	  ++countsAfter[mycounter];
      fillHistos_MuonPt(the_mu,hvec.at(mycounter),weight,algoName);
	  ++mycounter;
    }//NDPhiJ
  }//NETjet


}//-----------Optimize_Isolation






//---------------------------------------------------------------------------
void Optimize_JetVeto(const wprime::Event* ev, 
                      const wprime::Muon* the_mu,
                      int countsAfter[],vector<TH1F*>& hvec,
                      float o_ETJetVeto[],
                      const int o_NETjet,
                      const float& weight,
                      const string& algoName)

{
//---------------------------------------------------------------------------
 if(Debugme) cout<<"Running Optimize_IsoAndJetVeto"<<endl;
  //Loop over the colection of thresholds
  int mycounter = 0;
  for(int ll = 0; ll < o_NDPhiJ; ++ll){
    for(int kk = 0; kk < o_NETjet; ++kk){
      float etjcut = o_ETJetVeto[kk]; 
      float dphicut = o_DphiJetVeto[ll];
      unsigned njets = myMaxNjetsAboveThresh;
      if (ExceedMaxNumJetsOpposedToMu(njets,etjcut,dphicut,
                                     the_mu,ev)){++mycounter;continue;}
      ++countsAfter[mycounter];
      fillHistos_MuonPt(the_mu,hvec.at(mycounter),weight,algoName);
      ++mycounter;
    }//NETjet
  }//NPhiJ
  

}//-----------Optimize_JetVeto





// //---------------------------------------------------------------------------
// void Optimize_IsoAndJetVeto(const wprime::Event* ev, 
// 			    const wprime::Muon* the_mu,
//                             int countsAfter[],
//                             const int option,float o_pTIso[],
// 			    const int o_NpTIso,float o_ETJetVeto[],
// 			    const int o_NETjet)
// {
// //---------------------------------------------------------------------------
//  if(Debugme) cout<<"Running Optimize_IsoAndJetVeto"<<endl;
//   //Loop over the colection of thresholds
//   int mycounter = 0;
//   for(int ii = 0; ii < o_NiDR; ++ii){
//     for(int jj = 0; jj < o_NpTIso; ++jj){
//       for(int ll = 0; ll < o_NDPhiJ; ++ll){
// 	for(int kk = 0; kk < o_NETjet; ++kk){
// 	  unsigned dRidx = o_indexDeltaR[ii];
// 	  float sptcut = o_pTIso[jj];
// 	  float etjcut = o_ETJetVeto[kk]; 
// 	  float dphicut = o_DphiJetVeto[ll];
// 	  unsigned njets = myMaxNjetsAboveThresh;
// 	  if (option == 200){
// 	    if(!SumPtIsolation(the_mu,dRidx,sptcut)) {++mycounter; continue;}
// 	    if (ExceedMaxNumJetsOpposedToMu(njets,etjcut,dphicut,
// 					   the_mu,ev)){++mycounter;continue;}
// 	  }//option == 200
// 	  else if (option == 300){
// 	    if (ExceedMaxNumJetsOpposedToMu(njets,etjcut,dphicut,
// 					   the_mu,ev)){++mycounter;continue;}
// 	    if(!SumPtIsolation(the_mu,dRidx,sptcut)) {++mycounter; continue;}
// 	  }//option == 300
// 	  else {cout<<"No valid option for Optim_IsoAndJetVeto, something went really wrong, please check ... quiting ...."<<endl; abort();}
	  
	 
// 	  ++countsAfter[mycounter];
// 	  ++mycounter;
// 	}//NETjet
//       }//NDPhiJ
//     }//NpTIso
//   }//NiDR
  

// }//--------Optimize_IsoAndJetVeto()





//---------------------------------------------------------------------------
void Optimize_Qual(const wprime::Muon* theMu,const string& algoName,
		   int countsAfter[],vector<TH1F*>& hvec,
                   float o_pTQual[],
		   const int o_NpTQual,float o_Chi2Qual[],
                   const int o_NChi2Qual,const float& weight)
  
{
//---------------------------------------------------------------------------
  if(Debugme) cout<<"Running Optimize_Qual"<<endl;
  //Loop over the colection of thresholds
  int mycounter = 0;
  bool fill_entry[Num_trkAlgos] = {true, true, true};
  for(int ll = 0; ll < o_NChi2Qual; ++ll){
    for(int kk = 0; kk < o_NpTQual; ++kk){
      float ptqual = o_pTQual[kk]; 
      float chi2cut = o_Chi2Qual[ll];
      float muonetacut = Muon_Eta_Cut;
      //if
      //(!CheckQuality(theMu,fill_entry,ptqual,chi2cut,muonetacut)){++mycounter;continue;}
      //IMPORTANT: Now this is not treated as a cut in the main code, so
      //this is dummy
      CheckQuality(theMu,fill_entry,ptqual,chi2cut,muonetacut);
      ++mycounter;
      ++countsAfter[mycounter];
      fillHistos_MuonPt(theMu,hvec.at(mycounter),weight,algoName);
      ++mycounter;
    }//ptQual
  }//chi2Qual
  

}//-----------Optimize_Qual






// Optimize cuts
//---------------------------------------------------------------------------
void ExecOptimization(const wprime::InputFile& file,const string& dir,
                      const bool& loopMuons,
                      float optNbefore[],
                      float optNafter[][myNum_trkAlgos], 
                      TH1F* hbef[], vector<data>& haft,
                      const int option, float cutVals1[],
                      const int cv1size, float cutVals2[],
                      const int cv2size)
{
//---------------------------------------------------------------------------
    if(Debugme) cout<<"ExecOptimization option "<<option<<endl;
    
    //gather file basic info
    wprime::Event * ev = 0;
    int nevents = 0; float weight = 0.0;
    gatherFileBasicInfo(file,ev,nevents,weight);
    
    //loop over muon algos
    for (int mual = 0; mual<myNum_trkAlgos; ++mual){//algo loop
        
        const string algoName = algo_desc_short[mual];
        //counter (unweighted) events after cuts
        int countsBefore = 0;
        int countsAfter[o_maxNumCuts] = {0};
        
        //Loop over events:
        for(int i = 0; i != nevents; ++i){ // event loop
            if(Debugme) cout<<"##########Processing event # "<<i+1<<endl;
            
            file.tree->GetEntry(i);
            int nmuons = ev->mu->GetLast() + 1;
            wprime::Muon* theMu = 0;
            
            //treat W as signal or as background, 
            //cut at 300 GeV for the hard muon.
            if (dir == "Wsignal"){
                if(file.description == "W"){
                    if(NmuAboveThresh(300.,ev,theMu,algoName)<1) continue;
                    if (ExceedMaxNumJetsOpposedToMu(myMaxNjetsAboveThresh, 
                                                       100., myDelta_Phi,
                                                       theMu,ev)) continue;
                    //if (!CheckQuality(theMu,algoName, PtTrackCut,
                    //                Chi2Cut,Muon_Eta_Cut)) continue;
                }
            }
            if (dir == "Wback"){
                if(file.description == "W"){
                    //if(NmuAboveThresh(300.,ev,theMu,algoName)>0)
                    //continue;
                    if(NmuAboveThresh(300.,ev,theMu,algoName)>=1) continue;
                    if (theMu != 0){
                        if (!ExceedMaxNumJetsOpposedToMu(myMaxNjetsAboveThresh, 
                                                        100., 
                                                        myDelta_Phi,
                                                        theMu,ev)) continue;
                        
                    }
                    //if (CheckQuality(theMu,algoName,60,
                    //               Chi2Cut,Muon_Eta_Cut)) continue;
                    
                }
            }
            //loop over muons if needed
            //the loopMuons switch manipulates the muon loop
            //functionality.  NOT IMPLEMENTED HERE YET
            for (int mi = 0; mi < nmuons; ++mi){//loop over muons
                if(Debugme) cout<<"##########Processing muon #: "<<mi+1<<endl;
                
                //get the muon
                wprime::Muon* mu = (wprime::Muon *) ev->mu->At(mi);
                theMu = mu;
                
                //optimize according to requirements sets
                if (option == 100) {//optimize 1mu cut
//                    if (!PassedHLT(ev,theMu,loopMuons)) continue;
                    if (!PassedHLT(ev)) continue;
                    //if (!IsMuonPtInRange(theMu,algoName,myminPtMu,mymaxPtMu)) continue;
                    fillHistos_MuonPt(theMu,hbef[mual],weight,algoName);
                    ++countsBefore;
                    Optimize_OnlyOneMuon(ev,theMu,countsAfter,
                                         haft.at(mual).h1f,cutVals1,
                                         cv1size,weight,algoName);

                }//100 option
                else if (option == 200 || option == 300) {
                    //optimize iso or jet veto cuts
                    //after they have passed the trigger and the hight pt single
                    //muon requirements.
//                    if (!PassedHLT(ev,theMu,loopMuons)) continue;
                    if (!PassedHLT(ev)) continue;
                    //if (!IsMuonPtInRange(theMu,algoName,myminPtMu,mymaxPtMu)) continue;
                    if (!OnlyOneHighTrackPtMuon(ev,myOneMuPtTrackCut)) continue;
                    fillHistos_MuonPt(theMu,hbef[mual],weight,algoName);
                    ++countsBefore;
                    if (option == 200) Optimize_Isolation(theMu,countsAfter,
                                                          haft.at(mual).h1f,
                                                          cutVals1,cv1size,
                                                          weight,algoName);
                    if (option == 300) Optimize_JetVeto(ev,theMu,countsAfter,
                                                        haft.at(mual).h1f,
                                                        cutVals1,cv1size,
                                                        weight,algoName);
                }//200 or 300 option
                else if(option == 400 || option == 500){
                    //if (!PassedHLT(ev,theMu,loopMuons)) continue;
                    if (!PassedHLT(ev)) continue;
                    if (!OnlyOneHighTrackPtMuon(ev,myOneMuPtTrackCut)) continue;
                    if(option == 400){
                        if (!SumPtIsolation(theMu,mydeltaRIsoIndex,
                                            mySumPtCut)) continue;
                        fillHistos_MuonPt(theMu,hbef[mual],weight,algoName);
                        ++countsBefore;
                        Optimize_JetVeto(ev,theMu,countsAfter,
                                         haft.at(mual).h1f,
                                         cutVals1,cv1size,
                                         weight,algoName);
                    }
                    if(option == 500){
                        if (ExceedMaxNumJetsOpposedToMu(myMaxNjetsAboveThresh, 
                                                       myEtJetCut, myDelta_Phi,
                                                       theMu,ev)) continue;
                        fillHistos_MuonPt(theMu,hbef[mual],weight,algoName);
                        ++countsBefore;
                        Optimize_Isolation(theMu,countsAfter,haft.at(mual).h1f,
                                           cutVals1,cv1size,weight,algoName);
                    }
                }//400 OR 500 option
                else if(option == 600 || option == 700){
//                    if (!PassedHLT(ev,theMu,loopMuons)) continue;
                    if (!PassedHLT(ev)) continue;
                    if (!OnlyOneHighTrackPtMuon(ev,myOneMuPtTrackCut)) continue;
                    if(option == 600){
                        if (!SumPtIsolation(theMu,mydeltaRIsoIndex,
                                            mySumPtCut)) continue;
                        if (ExceedMaxNumJetsOpposedToMu(myMaxNjetsAboveThresh, 
                                                       myEtJetCut, myDelta_Phi,
                                                       theMu,ev)) continue;
                        fillHistos_MuonPt(theMu,hbef[mual],weight,algoName);
                        ++countsBefore;
                        Optimize_Qual(theMu,algoName,countsAfter,haft.at(mual).h1f,
                                      cutVals1,cv1size,
                                      cutVals2,cv2size,weight);
                    }
                    if(option == 700){
                        if (ExceedMaxNumJetsOpposedToMu(myMaxNjetsAboveThresh, 
                                                       myEtJetCut, myDelta_Phi,
                                                       theMu,ev)) continue;
                        if (!SumPtIsolation(theMu,mydeltaRIsoIndex,
                                            mySumPtCut)) continue;
                        fillHistos_MuonPt(theMu,hbef[mual],weight,algoName);
                        ++countsBefore;
                        Optimize_Qual(theMu,algoName,countsAfter,haft.at(mual).h1f,
                                      cutVals1,cv1size,
                                      cutVals2,cv2size,weight);
                    }
                }//600 OR 700 options
                else {cout<<"Nothing to be optimized"<<endl; return;}
                
                if(!loopMuons) {break;}
                
            }//muon loop
            
        }//event loop
      
      
      //Number of expected events for each cut (weighted)
      //depending of the optimization study:
      if (option == 100) {
          optNbefore[mual] += countsBefore * weight;
          for(int ii = 0; ii < cv1size; ++ii){
              optNafter[ii][mual] += countsAfter[ii]*weight;
          }
      }//option == 100
      else if (option == 200 || option == 300) {
          optNbefore[mual] += countsBefore * weight;
          int mycounter = 0;
          //the order of this loop is important, has to 
          //agree with the Optimize_Isolation loop.
          for(int ii = 0; ii < o_NiDR; ++ii){
              for(int jj = 0; jj < cv1size; ++jj){
                  optNafter[mycounter][mual] += 
                      countsAfter[mycounter] * weight;
                  ++mycounter;
              }//NpTIso
          }//NiDR
      }//option == 200 or 300
      else if (option == 400 || option == 500) {
          optNbefore[mual] += countsBefore * weight;
          int mycounter = 0;
          //the order of this loop is important, has to 
          //agree with the Optimize_JetVeto loop.
          for(int ll = 0; ll < o_NDPhiJ; ++ll){
              for(int kk = 0; kk < cv1size; ++kk){
                  optNafter[mycounter][mual] += 
                      countsAfter[mycounter] * weight;
                  ++mycounter;
              }//NETjet
          }//NDPhiJ
      }//option == 400 or 500
      else if (option == 600 || option == 700) {
          optNbefore[mual] += countsBefore * weight;
          int mycounter = 0;
          //the order of this loop is important, has to 
          //agree with the Optimize_JetVeto loop.
          for(int ll = 0; ll < cv2size; ++ll){
              for(int kk = 0; kk < cv1size; ++kk){
                  optNafter[mycounter][mual] += 
                      countsAfter[mycounter] * weight;
                  ++mycounter;
              }//cv1size
          }//cv2size
      }//option == 600 or 700
      else {cout<<"Something went wrong...quiting..."<<endl; abort();}
    }//muon algo loop
    
    delete ev;
    


}//-----------GetCutsOptimization()




//---------------------------------------------------------------------------
void FillCutVals(float cutVals[],const float cutMin, const float cutMax,
		 const int cutsNumSteps)
{
//---------------------------------------------------------------------------

  const int arraysize = cutsNumSteps + 1;
  float varStep = 0;
  if(cutsNumSteps>0) varStep = double((cutMax - cutMin)/cutsNumSteps);
  for(int i = 0; i < arraysize;++i){
    cutVals[i] = cutMin + i*varStep;
    if(Debugme) cout<<"cutVals["<<i<<"] = "<<cutVals[i]<<endl;
  }

}//---------- FillCutVals






//---------------------------------------------------------------------------
void  FillHistVectors(vector<TH1F*>& vh,const int option,
                      float cutVals1[],float cutVals2[], 
                      const int cv1size, const int cv2size, const string& algoName)
{
//---------------------------------------------------------------------------

    //There is a maximum of 0.001 presicion for the naming of histograms
    //see below.  If this presicion is surpassed, there might be memory
    //leaks and problems with histograms double definitions
    char buffer[1000];
    if (option == 100){
        for (int k = 0; k < cv1size; ++k){
            sprintf(buffer,"%s_%.3f_0.000",algoName.c_str(),cutVals1[k]);
            vh.push_back(new TH1F(buffer,"", mynBinPtMu,myminPtMu,mymaxPtMu));
        }
    }//----100
    else if (option == 200 || option == 500){
        for(int ii = 0; ii < o_NiDR; ++ii){
            for(int jj = 0; jj < cv1size; ++jj){
                sprintf(buffer,"%s_%.3f_%.3f",algoName.c_str(),cutVals1[jj],o_ValDeltaRIso[o_indexDeltaR[ii]]);
                vh.push_back(new TH1F(buffer,"", mynBinPtMu,myminPtMu,mymaxPtMu));
            }
        }
    }//----opt200 or 500
    else if (option == 300 || option == 400){
        for(int ll = 0; ll < o_NDPhiJ; ++ll){
            for(int kk = 0; kk < cv1size; ++kk){
                sprintf(buffer,"%s_%.3f_%.3f",algoName.c_str(),cutVals1[kk],o_DphiJetVeto[ll]);
                vh.push_back(new TH1F(buffer,"", mynBinPtMu,myminPtMu,mymaxPtMu));
            }
        }
    }//----opt300 or 400
    else if (option == 600 || option == 700){
        for(int ll = 0; ll < cv2size; ++ll){
            for(int kk = 0; kk < cv1size; ++kk){
                sprintf(buffer,"%s_%.3f_%.3f",algoName.c_str(),cutVals1[kk],cutVals2[ll]);
                vh.push_back(new TH1F(buffer,"", mynBinPtMu,myminPtMu,mymaxPtMu));
            }
        }
    }//----opt600 or 700
    else {
        cout<<"FillHistVectors:ERROR:This option is not defined"<<endl;
    }

    return;

}//---------GetExactNumCuts



//---------------------------------------------------------------------------
void GetCutsOptimization(const wprime::InputFile& file,
                         TFile* fout,
                         const int option, const bool loopMuons,
                         const float o_min1, const float o_max1, 
                         const int o_nstep1, const float o_min2,
                         const float o_max2, const int o_nstep2)
{
//---------------------------------------------------------------------------
    if(Debugme) cout<<"$$$$$$$$$$$$$$$$$$$$GetCutsOptimization..."<<endl;
    if(Debugme) cout<<"Option = "<<option<<endl;
   
    //calculate array of cut values
    const int cv1size = o_nstep1 + 1;
    const int cv2size = o_nstep2 + 1;
    float cutVals1[cv1size];
    float cutVals2[cv2size];
    FillCutVals(cutVals1,o_min1,o_max1,o_nstep1);
    FillCutVals(cutVals2,o_min2,o_max2,o_nstep2);

    //for optimization studies
    float optNbefore[myNum_trkAlgos] = {0};
    float optNafter[o_maxNumCuts][myNum_trkAlgos];
    for (int j = 0;j<o_maxNumCuts;++j)
        for(int i = 0;i<myNum_trkAlgos;++i)
            optNafter[j][i] = 0;
    
    //declare histograms for optimization studies
    //const int maxnumcuts = cv1size*cv2size*o_NiDR*o_NDPhiJ;
    TH1F* hbef[myNum_trkAlgos];
    
    
     vector<data> haft;
     //dirty way of distinguishing histos
     vector<TH1F*> vh0f;
     vector<TH1F*> vh1f;
     vector<TH1F*> vh2f;
     //define the histos in the vectors for the three algos
     //I don't like this fix!!!
     FillHistVectors(vh0f,option,cutVals1,cutVals2,cv1size,cv2size,algo_desc_short[0]);
     FillHistVectors(vh1f,option,cutVals1,cutVals2,cv1size,cv2size,algo_desc_short[1]);
     FillHistVectors(vh2f,option,cutVals1,cutVals2,cv1size,cv2size,algo_desc_short[2]);
     
     for(int al = 0; al<myNum_trkAlgos;++al){
         char buffer1[100];
         sprintf(buffer1,"hbef_%i",al);
         hbef[al] = new TH1F(buffer1,"", mynBinPtMu,myminPtMu,mymaxPtMu);
     }
     haft.push_back(data(algo_desc_short[0],vh0f));
     haft.push_back(data(algo_desc_short[1],vh1f));
     haft.push_back(data(algo_desc_short[2],vh2f));



    //loop over background and signal files
   

     if(!file.tree)
         return;
     cout << " Processing sample " << file.description << endl;
     string dir = file.samplename;
     if (option >= 100 && option <= 700) {
         ExecOptimization(file,dir,loopMuons,
                          optNbefore,optNafter,hbef,haft,option,cutVals1,
                          cv1size,cutVals2,cv2size);    
         
     }   
     else {cout<<"Option "<<option<<" is invalid,quiting..."<<endl;abort();}

   


    //Print the results if needed according to study case
    if (option >= 100 && option <= 700){ 
        printSummary_Optim(fout,dir, optNbefore, optNafter,hbef,haft,option,cutVals1,cv1size,cutVals2,cv2size);
    }
    else  {cout<<"Nothing to print for this study"<<endl;}

    if(Debugme) cout<<"deleting histos memories"<<endl;

    //delete histograms
//    for(int al = 0; al<myNum_trkAlgos;++al){
    //      delete hbef[al];
    //  for (int cu = 0; cu < maxnumcuts;++cu){
    //      delete haft[cu][al];
    //  }
    // }//delete histograms
    
}





