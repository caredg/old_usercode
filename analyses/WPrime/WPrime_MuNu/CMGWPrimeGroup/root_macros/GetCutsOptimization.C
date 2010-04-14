// ROOT49 stuff
#include <TROOT.h>
#include <TClonesArray.h>
#include <TFile.h>
#include <TTree.h>
#include <TMath.h>
#include <TH1F.h>

// wprime stuff
#include "UserCode/CMGWPrimeGroup/interface/wprimeEvent.h"
#include "UserCode/CMGWPrimeGroup/root_macros/input_file.h"
#include "UserCode/CMGWPrimeGroup/root_macros/loadCutsAndThresholds.h"

using std::cout; using std::endl; using std::vector; using std::string;


//to control the fixed parameters:
const int myNum_trkAlgos = Num_trkAlgos; // global, tracker, tev_1st
const float myOneMuPtTrackCut = OneMuPtTrackCut; 
//const float myOneMuPtTrackCut = 300; 
const float mySumPtCut = SumPtCut; // Cone DeltaR =0.3; 
const unsigned mydeltaRIsoIndex = deltaRIsoIndex; //for the isolation container
const float myEtJetCut = EtJetCut;
const float myDelta_Phi = Delta_Phi;//min angle muon/jet for jet veto
const unsigned myMaxNjetsAboveThresh = MaxNjetsAboveThresh;
const float myPtTrackCut = PtTrackCut;
const float myChi2Cut = Chi2Cut;


//optimization constants
//if not controlled here they are controlled in the MakeOptim.C macro.
const string optimdir = "forOptimization";
const int o_maxNumCuts = 100000;//max number of cuts
const int o_NiDR = 9;
const float o_ValDeltaRIso[o_NiDR] = {0.2,0.25,0.3,0.35,0.4,0.45,0.5,0.55,0.6};
const unsigned o_indexDeltaR[o_NiDR] = {0,1,2,3,4,5,6,7,8};
const int o_NDPhiJ = 5;
const float o_DphiJetVeto[o_NDPhiJ] = { TMath::Pi() - 0.1, 
                                       TMath::Pi() - 0.2, 
                                       TMath::Pi() - 0.3, 
                                       TMath::Pi() - 0.4,
                                       TMath::Pi() - 0.5};



//methods in this macro:
void populateOptimFile (const string& outfile, float optNbefore[],
			float optNafter[][myNum_trkAlgos],
			const int option,float cutVals1[],
			 const int cv1size, float cutVals2[],
			const int cv2size);
void printSummary_Optim(const string& dir,
                        float  optNbefore[], 
                        float optNafter[][myNum_trkAlgos],
                        const int option, float cutVals1[],
			 const int cv1size,  float cutVals2[],
			const int cv2size);
void gatherFileBasicInfo(const wprime::InputFile& file,
                         wprime::Event*& ev, int& nevents, float& weight);
void Optimize_OnlyOneMuon(const wprime::Event* ev, 
                          wprime::Muon*& the_mu,
                          int countsAfter[], float o_muonTrackPt[],
			  const int o_NmuTrkPt);
void Optimize_Isolation(const wprime::Muon* the_mu,
			int countsAfter[],
			 float o_pTIso[],
			const int o_NpTIso);
void Optimize_JetVeto(const wprime::Event* ev, 
		      const wprime::Muon* the_mu,
		      int countsAfter[],
		       float o_ETJetVeto[],
		      const int o_NETjet);
void Optimize_IsoAndJetVeto(const wprime::Event* ev, 
			    const wprime::Muon* the_mu,
                            int countsAfter[],
                            const int option, float o_pTIso[],
			    const int o_NpTIso, float o_ETJetVeto[],
			    const int o_NETjet);
void Optimize_Qual(const wprime::Muon* theMu,const string& algoName,
		   int countsAfter[],float o_pTQual[],
		   const int o_NpTQual,float o_Chi2Qual[],
		   const int o_NChi2Qual);
void FillCutVals(float cutVals[],const float cutMax, const float cutMin,
		 const int cutsNumSteps);
void ExecOptimization(const wprime::InputFile& file,const string& dir,
		      const bool& loopMuons,
                         float optNbefore[],
                         float optNafter[][myNum_trkAlgos], 
                         const int option, float cutVals1[],
			 const int cv1size,  float cutVals2[],
		      const int cv2size);






//-------------------------------------------------------------------------
void populateOptimFile (const string& outfile, float optNbefore[],
			float optNafter[][myNum_trkAlgos],
			const int option, float cutVals1[],
			 const int cv1size, float cutVals2[],
			 const int cv2size)
{
//-------------------------------------------------------------------------
  if(debugme) cout<<"Populating files for Optimization "<<endl;
  
  ofstream out(outfile.c_str());
  if(!out) {cout << "Cannot open file " << outfile << endl; abort();}
  
  if (option == 100){
    for(int j = 0; j < myNum_trkAlgos ; ++j){
      for (int k = 0; k < cv1size; ++k){
	out<<algo_desc[j]<<"\t"<<cutVals1[k]<<"\t"<<-9999<<"\t"
	   <<-9999<<"\t"<<-9999<<"\t"<<optNbefore[j]<<"\t"<<optNafter[k][j]<<endl;
      }//trkPt
    }//trkAlgos
  }//option 100
  else if (option == 200 || option == 300){
    for(int j = 0; j < myNum_trkAlgos ; ++j){
      int mycounter = 0;
      for(int ii = 0; ii < o_NiDR; ++ii){
	for(int jj = 0; jj < cv1size; ++jj){
	  for(int ll = 0; ll < o_NDPhiJ; ++ll){
	    for(int kk = 0; kk < cv2size; ++kk){
	      out<<algo_desc[j]<<"\t"<<cutVals1[jj]<<"\t"<<
		o_ValDeltaRIso[o_indexDeltaR[ii]]<<"\t"
		 <<cutVals2[kk]<<"\t"<<o_DphiJetVeto[ll]<<"\t"
		 <<optNbefore[j]<<"\t"<<
		optNafter[mycounter][j]<<endl;
		++mycounter;
	    }//NDPhiJ
	  }//NETjet
	}//NpTIso
      }//NiDR
    }//trkAlgos
  }//option 200 or 300
  else if (option == 400 || option == 700){
    for(int j = 0; j < myNum_trkAlgos ; ++j){
      int mycounter = 0;
      for(int ii = 0; ii < o_NiDR; ++ii){
	for(int jj = 0; jj < cv1size; ++jj){
	  out<<algo_desc[j]<<"\t"<<cutVals1[jj]<<"\t"
	     <<o_ValDeltaRIso[o_indexDeltaR[ii]]<<"\t"
	     <<-9999<<"\t"
	     <<-9999<<"\t"<<optNbefore[j]<<"\t"<<
	    optNafter[mycounter][j]<<endl;
	  ++mycounter;
	}//NpTIso
      }//NiDR
    }//trkAlgos
  }//option 400 or 700
  else if (option == 500 || option == 600){
    for(int j = 0; j < myNum_trkAlgos ; ++j){
      int mycounter = 0;
      for(int ll = 0; ll < o_NDPhiJ; ++ll){
	for(int kk = 0; kk < cv1size; ++kk){
	  out<<algo_desc[j]<<"\t"<<cutVals1[kk]<<"\t"<<
	    o_DphiJetVeto[ll]<<"\t"<<
	    -9999<<"\t"
	   <<-9999<<"\t"<<optNbefore[j]<<"\t"<<
	    optNafter[mycounter][j]<<endl;
	  ++mycounter;
	}//NDPhiJ
      }//NETjet
    }//trkAlgos
  }//option 500 or 600
  else if (option == 800 || option == 900){
    for(int j = 0; j < myNum_trkAlgos ; ++j){
      int mycounter = 0;
      for(int ll = 0; ll < cv2size; ++ll){
	for(int kk = 0; kk < cv1size; ++kk){
	  out<<algo_desc[j]<<"\t"<<cutVals1[kk]<<"\t"<<
	    cutVals2[ll]<<"\t"<<
	    -9999<<"\t"
	   <<-9999<<"\t"<<optNbefore[j]<<"\t"<<
	    optNafter[mycounter][j]<<endl;
	  ++mycounter;
	}//NDPhiJ
      }//NETjet
    }//trkAlgos
  }//option 800 or 900
}






//--------------------------------------------------------------------------
void printSummary_Optim(const string& dir,
                        float  optNbefore[], 
                        float optNafter[][myNum_trkAlgos],
                        const int option, float cutVals1[],
			 const int cv1size, float cutVals2[],
			 const int cv2size)
{
//------------------------------------------------------------------------
  if(debugme) cout<<"Printing summary for Optimization "<<endl;

  string outdir = optimdir;
    string cutType = "";
    string outfile = "";

    if (option == 100){cutType = cuts_desc[1];}//1mu
    else if (option == 200){
      cutType = cuts_desc[1]+cuts_desc[3]+"AND"+cuts_desc[4];//1muisoANDjet
    }
    else if (option == 300){
      cutType = cuts_desc[1]+cuts_desc[4]+"AND"+cuts_desc[3];//1mujetANDiso
    }
    else if (option == 400){
      cutType = cuts_desc[1]+cuts_desc[3];//1muiso
    }
    else if (option == 500){
      cutType = cuts_desc[1]+cuts_desc[4];//1mujet
    }
    else if (option == 600){
      cutType = cuts_desc[1]+cuts_desc[3]+cuts_desc[4];//1muisojet
    }
    else if (option == 700){
      cutType = cuts_desc[1]+cuts_desc[4]+cuts_desc[3];//1mujetiso
    }
    else if (option == 800){
      cutType = cuts_desc[1]+cuts_desc[3]+cuts_desc[4]+cuts_desc[5];//1mujetisoqual
    }
    else if (option == 900){
      cutType = cuts_desc[1]+cuts_desc[4]+cuts_desc[3]+cuts_desc[5];//1mujetisoqual
    }
    else {cout<<"Nothing to print for this optim study"<<endl; return;}

    outfile = outdir +"/"+ dir + "_" + cutType +".dat";
    populateOptimFile(outfile, optNbefore, optNafter, option,
		      cutVals1,cv1size,cutVals2,cv2size);

}//--------printSummary_Optim


//Routine to grab the information from a given file
//---------------------------------------------------------------------------
void gatherFileBasicInfo(const wprime::InputFile& file,
                         wprime::Event*& ev, int& nevents, float& weight)
{
//---------------------------------------------------------------------------

    ev = new wprime::Event();
    file.tree->SetBranchAddress("wp", &ev);
    nevents = file.tree->GetEntries();
    if(debugme) cout<<"Number of events in the file = "<<nevents<<endl;
    weight = file.weight;

}//---------gatherFileBasicInfo()




//Optimze OnlyOneMuon cut 
//---------------------------------------------------------------------------
void Optimize_OnlyOneMuon(const wprime::Event* ev, 
                          wprime::Muon*& the_mu,
                          int countsAfter[], float o_muonTrackPt[],
			  const int o_NmuTrkPt)
{
//---------------------------------------------------------------------------
  if(debugme) cout<<"Running Optimize_OnlyOneMuon"<<endl;
    //Loop over the colection of thresholds
    for (int nn = 0; nn < o_NmuTrkPt; ++nn){
        float pttrkcut = o_muonTrackPt[nn];
        if(!OnlyOneHighTrackPtMuon(ev,the_mu,pttrkcut)) continue;
        ++countsAfter[nn];
    }

    return;

}//-----Optimize_OnlyOneMuon()





//---------------------------------------------------------------------------
void Optimize_Isolation(const wprime::Muon* the_mu,
			int countsAfter[],
			float o_pTIso[],
			const int o_NpTIso)

{
//---------------------------------------------------------------------------
  if(debugme) cout<<"Running Optimize_Isolation"<<endl;
  int mycounter = 0;
  for(int ii = 0; ii < o_NiDR; ++ii){
    for(int jj = 0; jj < o_NpTIso; ++jj){
	  unsigned dRidx = o_indexDeltaR[ii];
	  float sptcut = o_pTIso[jj];
	  if(!SumPtIsolation(the_mu,dRidx,sptcut)) {++mycounter; continue;}
	  ++countsAfter[mycounter];
	  ++mycounter;
    }//NDPhiJ
  }//NETjet


}//-----------Optimize_Isolation






//---------------------------------------------------------------------------
void Optimize_JetVeto(const wprime::Event* ev, 
		      const wprime::Muon* the_mu,
		      int countsAfter[],
		      float o_ETJetVeto[],
		      const int o_NETjet)

{
//---------------------------------------------------------------------------
 if(debugme) cout<<"Running Optimize_IsoAndJetVeto"<<endl;
  //Loop over the colection of thresholds
  int mycounter = 0;
  for(int ll = 0; ll < o_NDPhiJ; ++ll){
    for(int kk = 0; kk < o_NETjet; ++kk){
      float etjcut = o_ETJetVeto[kk]; 
      float dphicut = o_DphiJetVeto[ll];
      unsigned njets = myMaxNjetsAboveThresh;
      if (ExeedMaxNumJetsOpposedToMu(njets,etjcut,dphicut,
				     the_mu,ev)){++mycounter;continue;}
      ++countsAfter[mycounter];
      ++mycounter;
    }//NETjet
  }//NPhiJ
  

}//-----------Optimize_JetVeto





//---------------------------------------------------------------------------
void Optimize_IsoAndJetVeto(const wprime::Event* ev, 
			    const wprime::Muon* the_mu,
                            int countsAfter[],
                            const int option,float o_pTIso[],
			    const int o_NpTIso,float o_ETJetVeto[],
			    const int o_NETjet)
{
//---------------------------------------------------------------------------
 if(debugme) cout<<"Running Optimize_IsoAndJetVeto"<<endl;
  //Loop over the colection of thresholds
  int mycounter = 0;
  for(int ii = 0; ii < o_NiDR; ++ii){
    for(int jj = 0; jj < o_NpTIso; ++jj){
      for(int ll = 0; ll < o_NDPhiJ; ++ll){
	for(int kk = 0; kk < o_NETjet; ++kk){
	  unsigned dRidx = o_indexDeltaR[ii];
	  float sptcut = o_pTIso[jj];
	  float etjcut = o_ETJetVeto[kk]; 
	  float dphicut = o_DphiJetVeto[ll];
	  unsigned njets = myMaxNjetsAboveThresh;
	  if (option == 200){
	    if(!SumPtIsolation(the_mu,dRidx,sptcut)) {++mycounter; continue;}
	    if (ExeedMaxNumJetsOpposedToMu(njets,etjcut,dphicut,
					   the_mu,ev)){++mycounter;continue;}
	  }//option == 200
	  else if (option == 300){
	    if (ExeedMaxNumJetsOpposedToMu(njets,etjcut,dphicut,
					   the_mu,ev)){++mycounter;continue;}
	    if(!SumPtIsolation(the_mu,dRidx,sptcut)) {++mycounter; continue;}
	  }//option == 300
	  else {cout<<"No valid option for Optim_IsoAndJetVeto, something went really wrong, please check ... quiting ...."<<endl; abort();}
	  
	 
	  ++countsAfter[mycounter];
	  ++mycounter;
	}//NETjet
      }//NDPhiJ
    }//NpTIso
  }//NiDR
  

}//--------Optimize_IsoAndJetVeto()





//---------------------------------------------------------------------------
void Optimize_Qual(const wprime::Muon* theMu,const string& algoName,
		   int countsAfter[],float o_pTQual[],
		   const int o_NpTQual,float o_Chi2Qual[],
		   const int o_NChi2Qual)
  
{
//---------------------------------------------------------------------------
  if(debugme) cout<<"Running Optimize_Qual"<<endl;
  //Loop over the colection of thresholds
  int mycounter = 0;
  for(int ll = 0; ll < o_NChi2Qual; ++ll){
    for(int kk = 0; kk < o_NpTQual; ++kk){
      float ptqual = o_pTQual[kk]; 
      float chi2cut = o_Chi2Qual[ll];
      float muonetacut = Muon_Eta_Cut;
      if (!HasQuality(theMu,algoName,ptqual,
		      chi2cut,muonetacut)){++mycounter;continue;}
      ++countsAfter[mycounter];
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
                         const int option, float cutVals1[],
			 const int cv1size, float cutVals2[],
			 const int cv2size)
{
//---------------------------------------------------------------------------
  if(debugme) cout<<"ExecOptimization option "<<option<<endl;

    //gather file basic info
    wprime::Event * ev = 0;
    int nevents = 0; float weight = 0.0;
    gatherFileBasicInfo(file,ev,nevents,weight);

    //loop over muon algos
    for (int mual = 0; mual<myNum_trkAlgos; ++mual){//algo loop
      
      const string algoName = algo_desc[mual];
      //counter (unweighted) events after cuts
      int countsBefore = 0;
      int countsAfter[o_maxNumCuts] = {0};
      
      //Loop over events:
      for(int i = 0; i != nevents; ++i){ // event loop
	if(debugme) cout<<"##########Processing event # "<<i+1<<endl;
	
	file.tree->GetEntry(i);
	int nmuons = ev->mu->GetLast() + 1;
	wprime::Muon* theMu = 0;

	//treat W as signal or as background, cut at 300 GeV for the hard muon.
	if (dir == "Wsignal"){
	  if(file.description == "W"){
	    if(NmuAboveThresh(300.,ev,theMu)<1) continue;
          }
	}
	if (dir == "Wback"){
	  if(file.description == "W"){
	    if(NmuAboveThresh(300.,ev,theMu)>0) continue;
          }
	}
	//loop over muons if needed
	//the loopMuons switch manipulates the muon loop
	//functionality.  NOT IMPLEMENTED HERE YET
	for (int mi = 0; mi < nmuons; ++mi){//loop over muons
	  if(debugme) cout<<"##########Processing muon #: "<<mi+1<<endl;
	  //cut index to keep track of the order of cuts
	  //get the muon
	  wprime::Muon* mu = (wprime::Muon *) ev->mu->At(mi);
	  theMu = mu;
	  
	  //optimize according to requirements sets
	  if (option == 100) {//optimize 1mu cut
	    if (!PassedHLT(ev,theMu,loopMuons)) continue;
	    //if (!IsMuonPtInRange(theMu,algoName,minPtMu,maxPtMu)) continue;
	    ++countsBefore;
	    Optimize_OnlyOneMuon(ev,theMu,countsAfter,cutVals1,cv1size);
	  }//100 option
	  else if (option >=200 && option <= 500) {
	    //optimize iso and/or jet veto cuts and/or viceversa
	    //after they have passed the trigger and the hight pt single
	    //muon requirements.
	    if (!PassedHLT(ev,theMu,loopMuons)) continue;
	    //if (!IsMuonPtInRange(theMu,algoName,minPtMu,maxPtMu)) continue;
	    if (!OnlyOneHighTrackPtMuon(ev,theMu,myOneMuPtTrackCut)) continue;
	    ++countsBefore;
	    if (option == 200 || option == 300){
	      Optimize_IsoAndJetVeto(ev,theMu,countsAfter,option,
				     cutVals1,cv1size,cutVals2,cv2size);
	    }
	    if (option == 400) Optimize_Isolation(theMu,countsAfter,
						  cutVals1,cv1size);
	    if (option == 500) Optimize_JetVeto(ev,theMu,countsAfter,
						cutVals1,cv1size);
	  }//200 to 500 option
	  else if(option == 600 || option == 700){
	    if (!PassedHLT(ev,theMu,loopMuons)) continue;
	    if (!OnlyOneHighTrackPtMuon(ev,theMu,myOneMuPtTrackCut)) continue;
	    if(option == 600){
	      if (!SumPtIsolation(theMu,mydeltaRIsoIndex,mySumPtCut)) continue;
	      ++countsBefore;
	      Optimize_JetVeto(ev,theMu,countsAfter,cutVals1,cv1size);
	    }
	    if(option == 700){
	      if (ExeedMaxNumJetsOpposedToMu(myMaxNjetsAboveThresh, 
					     myEtJetCut, myDelta_Phi,
					     theMu,ev)) continue;
	      ++countsBefore;
	      Optimize_Isolation(theMu,countsAfter,cutVals1,cv1size);
	    }
	  }//600 OR 700 option
	  else if(option == 800 || option == 900){
	    if (!PassedHLT(ev,theMu,loopMuons)) continue;
	    if (!OnlyOneHighTrackPtMuon(ev,theMu,myOneMuPtTrackCut)) continue;
	    if(option == 800){
	      if (!SumPtIsolation(theMu,mydeltaRIsoIndex,mySumPtCut)) continue;
	      if (ExeedMaxNumJetsOpposedToMu(myMaxNjetsAboveThresh, 
					     myEtJetCut, myDelta_Phi,
					     theMu,ev)) continue;
	      ++countsBefore;
	      Optimize_Qual(theMu,algoName,countsAfter,cutVals1,cv1size,
			    cutVals2,cv2size);
	    }
	    if(option == 900){
	      if (ExeedMaxNumJetsOpposedToMu(myMaxNjetsAboveThresh, 
					     myEtJetCut, myDelta_Phi,
					     theMu,ev)) continue;
	      if (!SumPtIsolation(theMu,mydeltaRIsoIndex,mySumPtCut)) continue;
	      ++countsBefore;
	      Optimize_Qual(theMu,algoName,countsAfter,cutVals1,cv1size,
			    cutVals2,cv2size);
	    }
	  }//800 OR 900 options
	  else {cout<<"Nothing to be optimized"<<endl; return;}
	  
	  if(!loopMuons) {break;}
	  
	}//muon loop
	  
      }//event loop
      
      
      //Number of expected events for each cut (weighted)
      //depending of the optimization study:
      if (option == 100) {
	optNbefore[mual] += countsBefore * weight;
	for(int ii = 0; ii < cv1size; ++ii){
	  optNafter[ii][mual] += countsAfter[ii] * weight;
	}
      }//option == 100
      else if (option == 200 || option == 300) {
	//cout<<"@@@@@@@@@@@@@@Algo = "<<algoName<<endl;
	optNbefore[mual] += countsBefore * weight;
	int mycounter = 0;
	//the order of this loop is important, has to 
	//agree with the Optimize_IsoAndJetVeto loop.
	for(int ii = 0; ii < o_NiDR; ++ii){
	  for(int jj = 0; jj < cv1size; ++jj){
	    for(int ll = 0; ll < o_NDPhiJ; ++ll){
	      for(int kk = 0; kk < cv2size; ++kk){
		optNafter[mycounter][mual] += 
		  countsAfter[mycounter] * weight;
		++mycounter;
	      }//NNETjet
	    }//NEDPhiJ
	  }//NpTIso
	}//NiDR
      }//option == 200 or 300
      else if (option == 400 || option == 700) {
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
      }//option == 400 or 700
      else if (option == 500 || option == 600) {
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
      }//option == 500 or 600
      else if (option == 800 || option == 900) {
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
      }//option == 800 or 900
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
    if(debugme) cout<<"cutVals["<<i<<"] = "<<cutVals[i]<<endl;
  }

}//---------- FillCutVals




//---------------------------------------------------------------------------
void GetCutsOptimization(const vector<wprime::InputFile>& files,
			 const string dir, 
			 const int option, const bool loopMuons,
			 const float o_min1, const float o_max1, 
			 const int o_nstep1, const float o_min2,
			 const float o_max2, const int o_nstep2)
{
//---------------------------------------------------------------------------
    if(debugme) cout<<"$$$$$$$$$$$$$$$$$$$$GetCutsOptimization..."<<endl;

    int Nfiles = files.size();

    //calculate array of cut values
    const int cv1size = o_nstep1 + 1;
    const int cv2size = o_nstep2 + 1;
    float cutVals1[cv1size];
    float cutVals2[cv2size];
    FillCutVals(cutVals1,o_min1,o_max1,o_nstep1);
    FillCutVals(cutVals2,o_min2,o_max2,o_nstep2);

    //for optimization studies
    float optNbefore[myNum_trkAlgos] = {0};
    float optNafter[o_maxNumCuts][myNum_trkAlgos] = {0};

    //loop over background and signal files
    for(int tr = 0; tr != Nfiles; ++tr){//loop over files

        if(!files[tr].tree)
            continue;
        cout << " Processing sample " << files[tr].description << endl;
        if (option >= 100 && option <= 900) 
	  ExecOptimization(files[tr],dir,loopMuons,
			   optNbefore,optNafter,option,cutVals1,
			   cv1size,cutVals2,cv2size);       
        else {cout<<"Option "<<option<<" is invalid,quiting..."<<endl;abort();}

    }//loop over files


    //Print the results if needed according to study case
    if (option >= 100 && option <= 900) printSummary_Optim(dir, optNbefore, optNafter,option,cutVals1,cv1size,cutVals2,cv2size);
    else  {cout<<"Nothing to print for this study"<<endl;}

}





