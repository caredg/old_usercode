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

//methods in this macro:
void printSummary_MuonPt(ofstream & out, const string& dir, 
                         const float& Nexp_evt, 
                         float Nexp_evt_cut[][Num_trkAlgos]);
void defineHistos_MuonPt();
void defineHistos_MuonPtJetIso();
void defineHistos_MuonChargePt();
void defineHistos(const int& option);
void passMethrough(int Num_surv_cut[], int& cut_index,bool accountme[]);
void tabulateMe(int Num_surv_cut[], int& cut_index, 
                 const float& weight,const wprime::Muon* mu,
                const string& algo_name, const int& option, bool accountme[]);
void fillHistos_MuonPt(int index, float weight,const wprime::Muon* mu,
                       const string& algo_name);
void fillHistos_MuonChargePt(int index, float weight,const wprime::Muon* mu,
                        const string& algo_name);
void saveHistos_MuonPt();
void saveHistos_MuonPtJetIso();
void saveHistos_MuonChargePt();
void saveHistos(TFile * fout, string dir, const int& option);
void gatherFileBasicInfo(const wprime::InputFile& file,
                         wprime::Event*& ev, int& nevents, float& weight);

// Studies:
void GetMuonPtDistribution_JetIso(const wprime::InputFile& file);
void GetMuonPtDistribution(const wprime::InputFile& file, 
                           float Nexp_evt_cut[][Num_trkAlgos],
                           float& Nexp_evt, const int& option,
			   const bool& loopMuons);








//--------------------------------------------------------------------------
void printSummary_MuonPt(ofstream & out, const string& dir, 
                         const float& Nexp_evt, 
                         float Nexp_evt_cut[][Num_trkAlgos])
//------------------------------------------------------------------------
{
    cout<<"\n$$$$$$$$$$$$$$$$$$$$$$$ Type of sample: "<<dir<<endl;
    cout << "Total # of expected events = " << Nexp_evt << endl;

    for (int mual = 0; mual<Num_trkAlgos; ++mual){
        cout<<"++++++ Algorithm: "<<algo_desc[mual]<<endl;
        for(int i = 0; i < Num_histo_sets; ++i){

            //print in a txt file the final total number of events
            if(i == Num_histo_sets - 1)
                out << algo_desc[mual] << " " << dir << " " 
                    << Nexp_evt_cut[i][mual] << endl;

            cout <<"Cut # "<<i<<": "<<cuts_desc[i]<<" -> expected evts = " 
                << Nexp_evt_cut[i][mual];
            
            //calculate efficiencies
            float eff, deff;
            if(i == 0)
                getEff(eff, deff, Nexp_evt_cut[i][mual], Nexp_evt);
            else
                getEff(eff, deff, Nexp_evt_cut[i][mual], 
                       Nexp_evt_cut[i-1][mual]);
            cout << ", Relative eff = "<<eff*100 << " +- " << deff*100 << "%";
            getEff(eff, deff, Nexp_evt_cut[i][mual], Nexp_evt);
            cout << ", Absolute eff = "<< eff*100 << " +- " << deff*100 << "%"
                << endl;
            
        } // loop over different cuts
    }//loop over different muon algos

}//------------- printSummary_MuonPt()





//--------------------------------------------------------------
void defineHistos_MuonPt()
{
//--------------------------------------------------------------
    if (debugme) cout<<"define Muon pT histos"<<endl;

     int index = 0;
  
  hPT[index][0]= new TH1F("hPTglb_trig","Global Muon Pt with HLT",
			  nBinPtMu,minPtMu,maxPtMu);
  hPT[index][1]= new TH1F("hPTtrk_trig","Tracker Muon Pt with HLT",
			  nBinPtMu,minPtMu,maxPtMu);
  hPT[index][2]= new TH1F("hPTtev_trig","TeV-1st Muon Pt with HLT",
			  nBinPtMu,minPtMu,maxPtMu);

  
  ++index;
  
  hPT[index][0]= new TH1F("hPTglb_1mu","Global Muon Pt 1 muon",
			  nBinPtMu,minPtMu,maxPtMu);
  hPT[index][1]= new TH1F("hPTtrk_1mu","Tracker Muon Pt 1 muon",
			  nBinPtMu,minPtMu,maxPtMu);
  hPT[index][2]= new TH1F("hPTtev_1mu","TeV-1st Muon Pt 1 muon",
			  nBinPtMu,minPtMu,maxPtMu);

  //switched places with 1mu, so the "all" tag does not make sense
  //any more. Fix!
  ++index;
 
   hPT[index][0] = new TH1F("hPTglb_all","Global Muon Pt",
 			   nBinPtMu,minPtMu,maxPtMu);
   hPT[index][1] = new TH1F("hPTtrk_all","Tracker Muon Pt",
 			   nBinPtMu,minPtMu,maxPtMu);
   hPT[index][2] = new TH1F("hPTtev_all","TeV-1st Muon Pt",
 			   nBinPtMu,minPtMu,maxPtMu);

  
   
  ++index;

  hPT[index][0]= new TH1F("hPTglb_iso","Global Muon Pt iso",
			  nBinPtMu,minPtMu,maxPtMu);
  hPT[index][1]= new TH1F("hPTtrk_iso","Tracker Muon Pt iso",
			  nBinPtMu,minPtMu,maxPtMu);
  hPT[index][2]= new TH1F("hPTtev_iso","TeV-1st Muon Pt iso",
			  nBinPtMu,minPtMu,maxPtMu);
  ++index;
  
  hPT[index][0]= new TH1F("hPTglb_jveto","Global Muon Pt jet veto",
			  nBinPtMu,minPtMu,maxPtMu);
  hPT[index][1]= new TH1F("hPTtrk_jveto","Tracker Muon Pt jet veto",
			  nBinPtMu,minPtMu,maxPtMu);
  hPT[index][2]= new TH1F("hPTtev_jveto","TeV-1st Muon Pt jet veto",
			  nBinPtMu,minPtMu,maxPtMu);
  ++index;
  
  hPT[index][0]= new TH1F("hPTglb_qual","Global Muon Pt qual",
			  nBinPtMu,minPtMu,maxPtMu);
  hPT[index][1]= new TH1F("hPTtrk_qual","Tracker Muon Pt qual",
			  nBinPtMu,minPtMu,maxPtMu);
  hPT[index][2] = new TH1F("hPTtev_qual","TeV-1st Muon Pt qual",
			   nBinPtMu,minPtMu,maxPtMu);
  ++index;

}//---------defineMuonPtHistos()




//--------------------------------------------------------------
void defineHistos_MuonPtJetIso()
{
//--------------------------------------------------------------
    if (debugme) cout<<"define histos JetIso"<<endl;
    
    hNMu = new TH1F("NMu","Nb. muons in the event",10,0.,10.);
    hPtMaxMu= new TH1F("ptMaxMu","Pt mu",nBinPtMu,minPtMu,maxPtMu);
    hPtMaxMuTrackVeto= new TH1F("ptMaxMuTrackVeto","Pt mu track Veto",
                                nBinPtMu,minPtMu,maxPtMu); 
    hPtMaxMuJetVeto= new TH1F("ptMaxMuJetVeto","Pt mu JetVeto",
                              nBinPtMu,minPtMu,maxPtMu); 
    hPtMaxMuTrackVetoJetVeto= new TH1F("ptMaxMuTrackVetoJetVeto",
                                       "Pt mu TrackVeto and JetVeto ",
                                       nBinPtMu,minPtMu,maxPtMu); 
    hPtMaxMu->SetLineColor(2); 
    hPtMaxMuTrackVeto->SetLineColor(4);
    hPtMaxMuJetVeto->SetLineColor(6); 

    return;
}//------defineHistos_JetIso



//--------------------------------------------------------------
void defineHistos_MuonChargePt()
{
//--------------------------------------------------------------
    if (debugme) cout<<"define histos MuonChargePt"<<endl;
    
    hPTplus[0]= new TH1F("hPTglb_plus","(+) Global Muon Pt qual",
                         nBinPtMu,minPtMu,maxPtMu);
    hPTplus[1]= new TH1F("hPTtrk_plus","(+) Tracker Muon Pt qual",
                         nBinPtMu,minPtMu,maxPtMu);
    hPTplus[2] = new TH1F("hPTtev_plus","(+) TeV-1st Muon Pt qual",
                          nBinPtMu,minPtMu,maxPtMu);
    
    hPTminus[0]= new TH1F("hPTglb_minus","(-) Global Muon Pt qual",
                          nBinPtMu,minPtMu,maxPtMu);
    hPTminus[1]= new TH1F("hPTtrk_minus","(-) Tracker Muon Pt qual",
                          nBinPtMu,minPtMu,maxPtMu);
    hPTminus[2] = new TH1F("hPTtev_minus","(-) TeV-1st Muon Pt qual",
                           nBinPtMu,minPtMu,maxPtMu);
    return;

}//------defineHistos_MuonChargePt





//Define histograms depending on the type of study
//--------------------------------------------------------------
void defineHistos(const int& option)
{
//--------------------------------------------------------------

  if (debugme) cout<<"define histos"<<endl;
  if (debugme) cout<<"option = "<<option<<endl;

  if(option == 1) defineHistos_MuonPt();
  else if(option == 2) defineHistos_MuonPtJetIso();
  else if(option == 3) defineHistos_MuonChargePt();
  else {
      cout<<"WARNING!!! No histos were defined"<<endl;
      return;
  }

}//----defineHistos()




//Just pass through if nothing is done
//-----------------------------------------------------------
void passMethrough(int Num_surv_cut[], int& cut_index, 
		   bool accountme[])
{
//-----------------------------------------------------------
 if(debugme) cout<<"Passing through...cut index: "<<cut_index<<endl;

    //if the accountme switch is on,
    //increase the number of events passing the cuts
    //and turn the switch off so we don't over count
    if (accountme[cut_index]) {
      ++Num_surv_cut[cut_index];
      accountme[cut_index] = false;
    }
    //increase the cut number
    ++cut_index;
    

}//------------passMetrhough


//tabulate results after the cut has been passed
//TODO: Improve flexibility
//-----------------------------------------------------------
void tabulateMe(int Num_surv_cut[], int& cut_index, 
                 const float& weight,const wprime::Muon* mu,
                const string& algo_name, const int& option,
		bool accountme[])
{
//-----------------------------------------------------------
    if(debugme) cout<<"Tabulating results for cut_index = "
                    <<cut_index<<endl;
    //if the accountme switch is on,
    //increase the number of events passing the cuts
    //and turn the switch off so we don't over count
    if (accountme[cut_index]) {
      ++Num_surv_cut[cut_index];
      accountme[cut_index] = false;
    }
    //fill the histograms
    if(option == 1) fillHistos_MuonPt(cut_index,weight,mu,algo_name);
    else if(option == 3) fillHistos_MuonChargePt(cut_index,weight,mu,algo_name);
    else cout<<"WARNING!! I can't tabulate anything for this study.."<<endl;

    //since the event has passed the cut,
    //increase the cut_index for the next cut
    ++cut_index;
    if(debugme) cout<<"cut_index is now = "<<cut_index<<endl;


}//-----tabulateMe






//Fill Histograms
//-----------------------------------------------------------
void fillHistos_MuonPt(int index, float weight,const wprime::Muon* mu,
                       const string& algo_name)
{
//-----------------------------------------------------------

    if(!mu) {
        if(debugme) cout<<"No muon found in the event"<<endl;
        return;
    }
    

    if (!strcmp(algo_name.c_str(),"gbl")){
        hPT[index][0]->Fill(mu->global.p.Pt(), weight);
    }
    else if (!strcmp(algo_name.c_str(),"trk")){
        hPT[index][1]->Fill(mu->tracker.p.Pt(), weight);
    }
    else if (!strcmp(algo_name.c_str(),"tev")){
        hPT[index][2]->Fill(mu->tev_1st.p.Pt(), weight);
    }
    else cout<<"WARNING!!!!!! Algo name not recognized"<<endl;
  

}//fillHistos_MuonPt





//Fill Histograms
//-----------------------------------------------------------
void fillHistos_MuonChargePt(int index, float weight,const wprime::Muon* mu,
                       const string& algo_name)
{
//-----------------------------------------------------------

    //do nothing if there is no muon in the event
    if(!mu) {
        if(debugme) cout<<"No muon found in the event"<<endl;
        return;
    }
    
    //only fill the required histograms
    if (index < (Num_histo_sets - Num_histo_sets_chargePt)) return;
   

    if (!strcmp(algo_name.c_str(),"gbl")){
        if(mu->global.q > 0)
            hPTplus[0]->Fill(mu->global.p.Pt(), weight);
        else
            hPTminus[0]->Fill(mu->global.p.Pt(), weight);
    }
    else if (!strcmp(algo_name.c_str(),"trk")){
        if(mu->tracker.q > 0)
            hPTplus[1]->Fill(mu->tracker.p.Pt(), weight);
        else
            hPTminus[1]->Fill(mu->tracker.p.Pt(), weight);
    }
    else if (!strcmp(algo_name.c_str(),"tev")){
        if(mu->tev_1st.q > 0)
            hPTplus[2]->Fill(mu->tev_1st.p.Pt(), weight);
        else
            hPTminus[2]->Fill(mu->tev_1st.p.Pt(), weight);
    }
    else cout<<"WARNING!!!!!! Algo name not recognized"<<endl;
    
}//---------------fillHistos_MuonChargePt()





//------------------------------------------------------------------------
void saveHistos_MuonPt()
{
//------------------------------------------------------------------------
    if(debugme) cout<<"Saving MuonPt histos"<<endl;
    
    for(int i = 0; i != Num_histo_sets; ++i)
        for(int j = 0; j != Num_trkAlgos; ++j)
            hPT[i][j]->Write();

}//-----saveHistos_MuonPt






//------------------------------------------------------------------------
void saveHistos_MuonPtJetIso()
{
//------------------------------------------------------------------------
    if(debugme) cout<<"Saving MuonPtJetIso histos"<<endl;
    
    hPtMaxMu->Write(); 
    hPtMaxMuTrackVetoJetVeto->Write(); 
    hPtMaxMuJetVeto->Write();  
    hPtMaxMuTrackVeto->Write();

}//-----saveHistos_MuonPtJetIso




//------------------------------------------------------------------------
void saveHistos_MuonChargePt()
{
//------------------------------------------------------------------------
    if(debugme) cout<<"Saving MuonPtJetIso histos"<<endl;
    
    for(int j = 0; j != Num_trkAlgos; ++j)
    {
      hPTplus[j]->Write();
      hPTminus[j]->Write();
    }
    return;

}//-----saveHistos_MuonPtJetIso




//------------------------------------------------------------------------
void saveHistos(TFile * fout, string dir, const int& option)
{
//------------------------------------------------------------------------

  fout->cd(); 
  fout->mkdir(dir.c_str()); 
  fout->cd(dir.c_str());

  if(option == 1) saveHistos_MuonPt();
  else if(option == 2) saveHistos_MuonPtJetIso();
  else if(option == 3) saveHistos_MuonChargePt();
  else if(option == 4) return;
  else {
      cout<<"Something went wrong when saving...Quiting..."<<endl;abort();}
  return;

}//saveHistos




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





//Study the pt distributions for global muons
//---------------------------------------------------------------------------
void GetMuonPtDistribution_JetIso(const wprime::InputFile& file)
{
//---------------------------------------------------------------------------

    //gather file basic info
    wprime::Event * ev = 0;
    int nevents = 0; float weight = 0.0;
    gatherFileBasicInfo(file,ev,nevents,weight);

    //Loop over events:
    for(int i = 0; i != nevents; ++i){ // event loop
        if(debugme) cout<<"##########Processing event # "<<i+1<<endl;
        
        file.tree->GetEntry(i);
        //Get number of muons in the event:
        int nmuons = ev->mu->GetLast() + 1;
        hNMu->Fill(nmuons,file.weight);
        //Get the the pT distribution for the hardest
        //global muon and for hardest global muon
        //whose trackPt is greater than the PtTrackCut:
        float ptMaxMu = 0.0;
        float ptMaxMuTrackVeto  = 0.0;
        GetGlobalMuonPtMax(ev,ptMaxMu,ptMaxMuTrackVeto);
        if(ptMaxMu > 0 ){hPtMaxMu->Fill(ptMaxMu,file.weight);}
        if(ptMaxMuTrackVeto > 0 ){
            hPtMaxMuTrackVeto->Fill(ptMaxMuTrackVeto,file.weight);}
        //same pT distributions as above but now requiring the
        //jet veto (ET)
        int nJetEt = NjetAboveThresh(EtJetCut,ev);
        if(nJetEt == 0){
            if(ptMaxMu > 0 ){hPtMaxMuJetVeto->Fill(ptMaxMu,file.weight);}
            if(ptMaxMuTrackVeto > 0) {
                hPtMaxMuTrackVetoJetVeto->Fill(ptMaxMuTrackVeto, file.weight);}
        }
        
        
    }//event loop
    
    delete ev;


}//-----------GetMuonPtDistribution_JetIso







//Main Analysis study
// TODO: Maybe rename the method to a more general one because
// in principle we can plot any variable in the event with this function, 
// not only muon pT. The tabulateMe function can be re-written to
// be a templated function.
//---------------------------------------------------------------------------
void GetMuonPtDistribution(const wprime::InputFile& file, 
                           float Nexp_evt_cut[][Num_trkAlgos],
                           float& Nexp_evt, const int& option,
			   const bool& loopMuons)
{
//---------------------------------------------------------------------------

    //gather file basic info
    wprime::Event * ev = 0;
    int nevents = 0; float weight = 0.0;
    gatherFileBasicInfo(file,ev,nevents,weight);

    //loop over muon algos
    for (int mual = 0; mual<Num_trkAlgos; ++mual){//algo loop
        const string algoName = algo_desc[mual];
        //counter (unweighted) events after cuts
        int Num_surv_cut[Num_histo_sets] = {0};

        //Loop over events:
        for(int i = 0; i != nevents; ++i){ // event loop
            if(debugme) cout<<"##########Processing event #: "<<i+1<<endl;
            
            file.tree->GetEntry(i);
            int nmuons = ev->mu->GetLast() + 1;
            wprime::Muon* theMu = 0;
	    // switch to indicate whether we account the tabulation
	    // for efficiency calculation or not.  This mainly has to
	    // do with muon looping, and basically takes care of
	    // increasing the Num_surv_cut container or not.  However,
	    // we do fill the histograms for every muon regardless.
	    bool accountme[Num_histo_sets];
	    for(int mm = 0;mm<Num_histo_sets;++mm){
	      accountme[mm] = true;
	    }

	    //loop over muons if needed
	    //the loopMuons switch manipulates the muon loop
	    //functionality
	    for (int mi = 0; mi < nmuons; ++mi){//loop over muons
	       if(debugme) cout<<"##########Processing muon #: "<<mi+1<<endl;
	      //cut index to keep track of the order of cuts
	      int cut_index = 0;
	      //get the muon
	      wprime::Muon* mu = (wprime::Muon *) ev->mu->At(mi);
	      theMu = mu;
	      // apply cuts
	      //>>>>>>>>>>CUT 1
	      if (!PassedHLT(ev,theMu,loopMuons)) continue;
	      tabulateMe(Num_surv_cut,cut_index,weight,theMu,algoName,option,
			 accountme);
	      //>>>>>>>>>>CUT 2
	      //only do the onemuon cut if not looping over muons
	      if (!loopMuons){
		if (!OnlyOneHighTrackPtMuon(ev,theMu,OneMuPtTrackCut)) continue;
		tabulateMe(Num_surv_cut,cut_index,weight,theMu,algoName,option, 
			   accountme);
	      }
	      else {
		//One could put a cut on pT on the muon; for now
		//just pass through
		passMethrough(Num_surv_cut,cut_index,accountme);
	      }

	      //>>>>>>>>>>CUT 3
	      if (!IsMuonPtInRange(theMu,algoName,minPtMu,maxPtMu)) continue;
	      tabulateMe(Num_surv_cut,cut_index,weight,theMu,algoName,option,
			 accountme);
	      
	      //>>>>>>>>>>CUT 4
	      if (!SumPtIsolation(theMu,deltaRIsoIndex,SumPtCut)) continue;
	      tabulateMe(Num_surv_cut,cut_index,weight,theMu,algoName,option,
			 accountme);
	      //>>>>>>>>>>CUT 5
	      if (ExeedMaxNumJetsOpposedToMu(MaxNjetsAboveThresh, 
					     EtJetCut, Delta_Phi,
					     theMu,ev)) continue;
	      tabulateMe(Num_surv_cut,cut_index,weight,theMu,algoName,option,
			 accountme);
	      //>>>>>>>>>>CUT 6
	      if (!HasQuality(theMu,algoName,PtTrackCut,
			      Chi2Cut,Muon_Eta_Cut)) continue;
	      tabulateMe(Num_surv_cut,cut_index,weight,theMu,algoName,option,
			 accountme);
	      //break if not looping over muons because the one muon pt cut
	      //already selected the good one for us.
	      if(!loopMuons) {break;}

	    }//muon loop
        }//event loop
        
        //Number of expected events for each cut (weighted)
        for(int ii = 0; ii < Num_histo_sets; ++ii){
            Nexp_evt_cut[ii][mual] += Num_surv_cut[ii] * weight;
        }
    }//muon algo loop

    //total # of events (before any cuts)
    Nexp_evt += nevents * weight; 
    delete ev;


}//-----------GetMuonPtDistribution()






//---------------------------------------------------------------------------
void GetDistributionGeneric(const vector<wprime::InputFile>& files, 
                            TFile *fout, string dir, ofstream & out, 
                            const int option = 1, 
			    const bool loopMuons = false)
{
//---------------------------------------------------------------------------
    if(debugme) cout<<"$$$$$$$$$$$$$$$$$$$$GetDistributionGeneric"<<endl;

    //Define histograms according to the type of study
    defineHistos(option);
    
    int Nfiles = files.size();

    //initialize counters to be used in studies
    float Nexp_evt = 0;
    float Nexp_evt_cut[Num_histo_sets][Num_trkAlgos] = {0};

    //loop over background and signal files
    for(int tr = 0; tr != Nfiles; ++tr){//loop over files

        if(!files[tr].tree)
            continue;
        cout << " Processing sample " << files[tr].description << endl;
        
        if (option == 1 || option == 3) GetMuonPtDistribution(files[tr],
							      Nexp_evt_cut,
							      Nexp_evt,option,
							      loopMuons);
        else if (option == 2) GetMuonPtDistribution_JetIso(files[tr]);
        else {cout<<"Option "<<option<<" is invalid,quiting..."<<endl;abort();}

    }//loop over files


    //Print the results if needed according to study case
    if (option == 1 || option == 3) printSummary_MuonPt(out, dir, Nexp_evt, Nexp_evt_cut);
    else  {cout<<"Nothing to print for this study"<<endl;}

    //save histograms according to study case
    saveHistos(fout,dir,option);
    
}





