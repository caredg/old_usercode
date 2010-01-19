// ROOT stuff
#include <TROOT.h>
#include <TClonesArray.h>
#include <TFile.h>
#include <TTree.h>
#include <TMath.h>
#include <TH1F.h>

// wprime stuff
#include "UserCode/CMGWPrimeGroup/interface/wprimeEvent.h"
#include "UserCode/CMGWPrimeGroup/root_macros/input_file.h"
#include "UserCode/CMGWPrimeGroup/root_macros/loadCuts.h"
#include "UserCode/CMGWPrimeGroup/root_macros/util.h"

// std stuff
#include <iostream>
#include <fstream>
#include <string>

using std::cout; using std::endl; using std::vector; using std::string;
const bool debugme = false;
const string algoType[Num_trkAlgos] = {"glb", "trk", "tev"};
const string varType[3] = {"iso", "jet", "pt"};




//--------------------------------------------------------------------------
void printSummary(const string& outdir, const string& dir, const string& cutType,
                  const float& Nexp_evt, float Nexp_evt_befcuts[], 
                  float Nexp_evt_cut[][Num_trkAlgos], const int& cutsNumSteps,
                  double cutVals[]) 
{ 
//--------------------------------------------------------------------------
    cout << " Total # of expected events in "<<dir<<" sample = " << Nexp_evt << endl;

    //container txt file so we can dump the results there
    //One file per sample and cut type
    string outfile = outdir +"/"+ dir + "_" + cutType +".dat";
    if (debugme) cout<<"File "<<outfile<<" has been created"<<endl;
    ofstream out(outfile.c_str());
    if(!out) { 
        cout << "Cannot open file " << outfile << endl; 
        abort();
    }

    //populate the *.dat files for optimization
    // algo cut_value total_events_before_cut total_events_after_cut
    for(int j = 0; j < Num_trkAlgos; ++j){
        for (int k = 0; k < cutsNumSteps; ++k){
            out<<algoType[j]<<"\t"<<cutVals[k]<<"\t"<<Nexp_evt_befcuts[j]<<"\t"<<Nexp_evt_cut[k][j]<<endl;
        }
    }

    out.close();
}

//---------------------------------------------------------------------------
void  Get_cut_decisions (const wprime::Muon* mu, bool before_cut[],
                         bool survived_cut[][Num_trkAlgos], 
                         const string& cutType, double cutVals[], const int& cutsNumSteps )
{
//--------------------------------------------------------------------------
    if(debugme) cout<<"$$$$$$$$$$$$$$$$$$$$$$$$$$$Get_cut_decisions"<<endl;
    if(debugme) {
        cout<<"Check booleans before loop over cut values"<<endl;
        for (int j=0;j<cutsNumSteps;++j){
            for (int k=0;k<Num_trkAlgos;++k){
                cout<<"survived_cut["<<j<<"]["<<k<<"] = "<<survived_cut[j][k]<<endl;
            }
        }
    }
    //loop over all cut values
    for (int i = 0; i < cutsNumSteps; ++i){
        if(debugme) cout<<"Cut value is = "<<cutVals[i]<<endl;
        //separate cases depending on the type of cut
        //Isolation
        if (cutType == varType[0]){
            for (int k=0;k<Num_trkAlgos;++k){
                if (!before_cut[k]) continue; 
                if (mu->SumPtIso[2] > cutVals[i]) continue;
                survived_cut[i][k] = true ;
                if(debugme) cout<<"survived_cut["<<i<<"]["<<k<<"] = "<<survived_cut[i][k]<<endl ;
            }
        }
        
    }//loop over all cut values
    
    return;
}



//--------------------------------------------------------------------------------------------
void Find_hard_muon_diffAlgos(wprime::Muon * mu, bool before_cut[])
{
//------------------------------------------------------------------------------------------
    if(debugme) cout<<"$$$$$$$$$$$$$$$$$$$$$$$$$$$Find_hard_muon_diffAlgos"<<endl;
    bool ptrange_glb = mu->global.p.Pt() >=minPtMu 
        && mu->global.p.Pt() <=maxPtMu;
    bool ptrange_trk = mu->tracker.p.Pt() >=minPtMu 
        && mu->tracker.p.Pt() <=maxPtMu;
    bool ptrange_tev = mu->tev_1st.p.Pt() >=minPtMu 
        && mu->tev_1st.p.Pt() <=maxPtMu;
    
    if(ptrange_glb) before_cut[0] = true;
    if(ptrange_trk) before_cut[1] = true;
    if(ptrange_tev) before_cut[2] = true;
    
}





 //---------------------------------------------------------------------------
void Loop_over_muons(const int& nmuons, const wprime::Event* ev, 
                     const bool& HLT, const bool& OneMuon,
                     bool before_cut[],
                     bool survived_cut[][Num_trkAlgos], 
                     const string& cutType, double cutVals[], const int& cutsNumSteps)
{
  //---------------------------------------------------------------------------
    if(debugme) cout<<"$$$$$$$$$$$$$$$$$$$$$$$$$$$Loop_over_muons"<<endl;

    // loop over muons
    for(int j = 0; j != nmuons; ++j){ 
        if(debugme) cout<<"&&&&&&&&&&&&Muon # "<<j<<endl;
        
        wprime::Muon * mu = (wprime::Muon *) ev->mu->At(j);
        
  	    // trigger decision
   	    if(!HLT) continue;
        if (debugme) cout<<"Passed HLT requirement"<<endl;
        
        // only one muon with tracker pt above OneMuPtTrackCut
  	    if(!OneMuon) continue;
        if (debugme) cout<<"Passed OneMuon requirement"<<endl;

        if (debugme){
            cout<<"Check booleans after HLT and OneMuon cuts"<<endl;
            for (int kk = 0;kk<Num_trkAlgos;++kk){
                cout<<"before_cut["<<kk<<"] = "<<before_cut[kk]<<endl;
            }
        }
        
        //find out if there is a hard muon by different algos
        //and store the results
        Find_hard_muon_diffAlgos(mu,before_cut);

        if (debugme){
            cout<<"Check booleans after basic cuts"<<endl;
            for (int jj = 0;jj<Num_trkAlgos;++jj){
                cout<<"before_cut["<<jj<<"] = "<<before_cut[jj]<<endl;
            }
        }
        
        
        //get the decisions according to the variable to be
        //optimized
        Get_cut_decisions(mu, before_cut, survived_cut,
                          cutType, cutVals, cutsNumSteps);

        if(debugme) {
        cout<<"Check booleans after loop over cut values"<<endl;
        for (int jj=0;jj<cutsNumSteps;++jj){
            for (int kk=0;kk<Num_trkAlgos;++kk){
                cout<<"survived_cut["<<jj<<"]["<<kk<<"] = "<<survived_cut[jj][kk]<<endl;
            }
        }
    }
    }//loop muons
}





//---------------------------------------------------------------------------
void GetCutsOptimization(const vector<wprime::InputFile>& files,const string& outdir,
                         const string& dir, const string& cutType,const double& cutMin,
                         const double& cutMax, const int& cutsNumSteps)
{
//---------------------------------------------------------------------------
    if(debugme) cout<<"$$$$$$$$$$$$$$$$$$$$$$$$$$$$GetCutsOptimization"<<endl;
    
    int Nfiles = files.size();
    
    //calculate array of cut values
    const double varStep = double((cutMax - cutMin)/cutsNumSteps);
    double cutVals[cutsNumSteps];
    for (int i = 0; i < cutsNumSteps; ++i){
        cutVals[i] = cutMin + i*varStep;
    }
    
    //initialize  counters
    //for expected (weighted) events
    float Nexp_evt = 0;
    float Nexp_evt_befcuts[Num_trkAlgos];
    float Nexp_evt_cut[cutsNumSteps][Num_trkAlgos];
    for (int j=0;j<cutsNumSteps;++j){
        for (int k=0;k<Num_trkAlgos;++k){
            Nexp_evt_befcuts[k] = 0;
            Nexp_evt_cut[j][k] = 0;            
        }
    }
    
    
     //loop over background and signal files
    for(int tr = 0; tr != Nfiles; ++tr){//loop over files
        
        if(!files[tr].tree)
            continue;
        cout << " Processing sample " << files[tr].description << endl;
        wprime::Event * ev = new wprime::Event();
        files[tr].tree->SetBranchAddress("wp", &ev);
        
        int nevents = files[tr].tree->GetEntries();
        if(debugme) cout<<"Number of events in the file = "<<nevents<<endl;
        float weight = files[tr].weight;
        
        //unweighted events counters
        int Num_basic_cuts[Num_trkAlgos] = {0};
        int Num_surv_cut[cutsNumSteps][Num_trkAlgos];
        for (int j=0;j<cutsNumSteps;++j){
            for (int k=0;k<Num_trkAlgos;++k){
                Num_surv_cut[j][k] = 0;
            }
        }
        
        //loop over events
        for(int i = 0; i != nevents; ++i){ // event loop
            if(debugme) cout<<"############################################Processing event # "<<i+1<<endl;
            
             files[tr].tree->GetEntry(i);

             //switches to keep track of passing status
             //before cuts that are gonna be optimized (before)
             //and after (survived) 
             bool before_cut[Num_trkAlgos] = {false};
             bool survived_cut[cutsNumSteps][Num_trkAlgos];
             for (int j = 0; j<cutsNumSteps;++j){
                 for (int k=0;k<Num_trkAlgos;++k){
                     survived_cut[j][k] = false;
                     }
             }
             int nmuons = ev->mu->GetLast() + 1;
             if (debugme) cout<<"Number of muons is = "<<nmuons<<endl;
             bool HLT = ev->HLT_Mu9;
             bool OneMuon = (NmuAboveThresh(OneMuPtTrackCut, ev) == 1);

             //loop over muons
             Loop_over_muons(nmuons,ev, HLT, OneMuon, before_cut, 
                             survived_cut, cutType, cutVals, cutsNumSteps);
             //Fill the event counters
             for (int k=0;k<Num_trkAlgos;++k){
                 if (before_cut[k]) ++Num_basic_cuts[k];
                 for (int j=0;j<cutsNumSteps;++j){
                     if (survived_cut[j][k]) ++Num_surv_cut[j][k];
                     if (debugme) cout<<"Num_basic_cuts["<<k<<"] = "<<
                         Num_basic_cuts[k]<<"\tNum_surv_cut["<<j<<"]["<<k<<"] = "
                                      <<Num_surv_cut[j][k]<<endl;
                 }
             }
              
        }//event loop

         //total # of events (before any cuts)
        Nexp_evt += nevents * weight; 
        //total num events before and after cuts
        for (int k=0;k<Num_trkAlgos;++k){
            Nexp_evt_befcuts[k] += Num_basic_cuts[k]*weight;
            for (int j=0;j<cutsNumSteps;++j){
                Nexp_evt_cut[j][k] += Num_surv_cut[j][k]*weight;
                if (debugme) cout<<"Nexp_evt_befcuts["<<k<<"] = "<<
                    Nexp_evt_befcuts[k]<<"\tNexp_evt_cut["<<j<<"]["<<k<<"] = "
                                 <<Nexp_evt_cut[j][k]<<endl;
            }
        }

        delete ev;
        
    }//loop over files
    
    printSummary(outdir,dir, cutType,Nexp_evt, Nexp_evt_befcuts, Nexp_evt_cut, cutsNumSteps, cutVals);
    
}














