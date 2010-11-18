// ROOT49 stuff
#include <TROOT.h>
#include <TClonesArray.h>
#include <TFile.h>
#include <TTree.h>
#include <TMath.h>
#include <TH1F.h>
#include <TH1D.h>
#include <TH2D.h>

#include <map>

// wprime stuff
#include "UserCode/CMGWPrimeGroup/interface/wprimeEvent.h"
#include "UserCode/CMGWPrimeGroup/root_macros/input_file.h"
#include "UserCode/CMGWPrimeGroup/root_macros/loadCutsAndThresholds.h"

using std::cout; using std::endl; using std::string; using std::map;

//methods in this macro:
void printHighPtMuon(const wprime::Event * ev, const wprime::Muon* mu,
                     const bool fill_entry[][Num_trkAlgos],bool isW, TH1D* h,
                TH2D * hh);
void printSummary_MuonPt(ofstream& out, wprime::InputFile & file, float Nexp_evt);
void defineHistos_MuonPt();
void defineHistos_MuonEta();
void defineHistos_MuonPhi();
void defineHistos_MuonJetDPhi();
void defineHistos_MuonIso();
void defineHistos_TMass();
void defineHistos_MuonChargePt();
void defineHistos(int option);
void tabulateMe(int Num_surv_cut[][Num_trkAlgos][Num_flavors], int& cut_index, 
		float weight,const wprime::Event* ev, const wprime::Muon* mu,
                const bool fill_entry[][Num_trkAlgos], int option, 
                bool accountme[][Num_trkAlgos][Num_flavors],bool isW,TH1D* h,
                TH2D * hh);
void fillHistos_MuonPt(int index, float weight, const TLorentzVector **p,
                       const bool fill_entry[][Num_trkAlgos]);
void fillHistos_MuonEta(int index, float weight, const TLorentzVector **p,
                       const bool fill_entry[][Num_trkAlgos]);
void fillHistos_MuonPhi(int index, float weight, const TLorentzVector **p,
                        const bool fill_entry[][Num_trkAlgos]);

void fillHistos_MuonJetDPhi(int index, float weight, const wprime::Event * ev,
                        const TLorentzVector **p, 
			    const bool fill_entry[][Num_trkAlgos]);
void fillHistos_TMass(int index, float weight, const wprime::Event * ev,
		      const TLorentzVector **p, 
                      const bool fill_entry[][Num_trkAlgos], bool isW, TH1D* h,
                TH2D * hh);
void fillHistos_MuonIso(int index, float weight, const wprime::Muon* mu,
			const bool fill_entry[][Num_trkAlgos]);
void fillHistos_MuonChargePt(int index, float weight, const Int_t * Q, 
			     const TLorentzVector **p,
			     const bool fill_entry[][Num_trkAlgos]);
void saveHistos_MainAnalysis();
void saveHistos_MuonChargePt();
void saveHistos(TFile * fout, string dir, int option);
void gatherFileBasicInfo(const wprime::InputFile& file,
                         wprime::Event*& ev, int & nevents, float & weight);

// Studies:
void GetMuonPtDistribution(wprime::InputFile& file, float & Nexp_evt, 
			   int option, bool highestPtMuonOnly);


//--------------------------------------------------------------------------
void printSummary_MuonPt(ofstream & out, wprime::InputFile & file, float Nexp_evt)
//------------------------------------------------------------------------
{
  out << "\n Sample: " << file.samplename << endl;
  out << " Total # of expected events = " << Nexp_evt << endl;
  
  for(int f = 0; f != Num_flavors; ++f)
    { // loop over pt & mt

      out << FLAVOR_NAME[f] << endl;

      for (int mual = MuAlgo_MIN; mual <= MuAlgo_MAX; ++mual){
	out << " Algorithm: " << algo_desc_long[mual] << endl;
	float eff, deff;
	
	for(int i = 0; i != Num_selection_cuts; ++i){
	  
	  out << " Cut # " << i << ": " << cuts_desc_long[i] 
	      <<", expected # of evts = " << file.Nexp_evt_cut[i][mual][f];
	  
	  //calculate efficiencies
	  if(i == 0)
	    getEff(eff, deff, file.Nexp_evt_cut[i][mual][f]/(file.weight),
		   Nexp_evt/(file.weight));
	  else
	    getEff(eff, deff, file.Nexp_evt_cut[i][mual][f]/(file.weight),
		   file.Nexp_evt_cut[i-1][mual][f]/(file.weight));
	  out << ", Relative eff = "<<eff*100 << " +- " << deff*100 << "%";
	  getEff(eff, deff, file.Nexp_evt_cut[i][mual][f]/(file.weight), 
		 Nexp_evt/(file.weight));
	  out << ", Absolute eff = "<< eff*100 << " +- " << deff*100 << "%"
	      << endl;
	  
	  file.eff[i][mual][f] = eff;
	  file.deff[i][mual][f] = deff;

	} // loop over different cuts
	
      }//loop over different muon algos

    } // loop over pt & mt

  
}//------------- printSummary_MuonPt()


//--------------------------------------------------------------
void defineHistos_MuonPt()
{
//--------------------------------------------------------------
#if debugme
  cout << " define Muon pT histos" << endl;
#endif

  for(int cut = 0; cut != Num_selection_cuts; ++cut)
    for(int algo = MuAlgo_MIN; algo <= MuAlgo_MAX; ++algo)
      {
	string name = "hPT" + algo_desc_short[algo] + "_" + 
	  cuts_desc_short[cut];
	string title = algo_desc_long[algo] + " muon p_{T} with " + 
	  cuts_desc_long[cut];
	hPT[cut][algo] = new TH1F(name.c_str(), title.c_str(), nBinPtMu,
				  minPtMu,maxPtMu);
      }
  
}//---------defineMuonPtHistos()



//--------------------------------------------------------------
void defineHistos_TMass()
{
//--------------------------------------------------------------
#if debugme
  cout << " define TMass histos" << endl;
#endif

  for(int cut = 0; cut != Num_selection_cuts; ++cut)
    for(int algo = MuAlgo_MIN; algo <= MuAlgo_MAX; ++algo)
      {
	string name = "hTM" + algo_desc_short[algo] + "_" 
	  + cuts_desc_short[cut];
	string title = algo_desc_long[algo] + " Transv. Mass with " 
	  + cuts_desc_long[cut];
	hTM[cut][algo] = new TH1F(name.c_str(), title.c_str(), nBinTmMu,
				  minTmMu,maxTmMu);
      }
  
}//---------defineMuonPtHistos()



//--------------------------------------------------------------
void defineHistos_MuonEta()
{
//--------------------------------------------------------------
#if debugme
  cout << " define Muon ETA histos" << endl;
#endif

  for(int cut = 0; cut != Num_selection_cuts; ++cut)
    for(int algo = MuAlgo_MIN; algo <= MuAlgo_MAX; ++algo)
      {
	string name = "hETA" + algo_desc_short[algo] + "_" 
	  + cuts_desc_short[cut];
	string title = algo_desc_long[algo] + " muon #eta with " 
	  + cuts_desc_long[cut];
	hETA[cut][algo] = new TH1F(name.c_str(), title.c_str(), nBinEtaMu,
				   minEtaMu,maxEtaMu);
      }
  
}//---------defineMuonEtaHistos()



//--------------------------------------------------------------
void defineHistos_MuonPhi()
{
//--------------------------------------------------------------
#if debugme
  cout << " define Muon PHI histos" << endl;
#endif

  for(int cut = 0; cut != Num_selection_cuts; ++cut)
    for(int algo = MuAlgo_MIN; algo <= MuAlgo_MAX; ++algo)
      {
	string name = "hPHI" + algo_desc_short[algo] + "_" + 
	  cuts_desc_short[cut];
	string title = algo_desc_long[algo] + " muon #phi with " + 
	  cuts_desc_long[cut];
	hPHI[cut][algo] = new TH1F(name.c_str(), title.c_str(), nBinPhiMu,
				   minPhiMu,maxPhiMu);
      }
  
}//---------defineMuonPhiHistos()



//--------------------------------------------------------------
void defineHistos_MuonJetDPhi()
{
//--------------------------------------------------------------
#if debugme
  cout << " define Muon-Jet DPhi histos" << endl;
#endif

  for(int cut = 0; cut != Num_selection_cuts; ++cut)
    for(int algo = MuAlgo_MIN; algo <= MuAlgo_MAX; ++algo)
      {
	string name = "hMJDPHI" + algo_desc_short[algo] + "_" + 
	  cuts_desc_short[cut];
	string title = algo_desc_long[algo] + " muon-jet #Delta#phi with " + 
	  cuts_desc_long[cut];
	hMJDPHI[cut][algo] = new TH1F(name.c_str(), title.c_str(), 
				      nBinDPhiMu,minDPhiMu,maxDPhiMu);
      }
  
}//---------defineMuonJetDPhiHistos()



//--------------------------------------------------------------
void defineHistos_MuonIso()
{
//--------------------------------------------------------------
#if debugme
  cout << " define Muon ISO histos" << endl;
#endif

  for(int cut = 0; cut != Num_selection_cuts; ++cut)
    for(int algo = MuAlgo_MIN; algo <= MuAlgo_MAX; ++algo)
      {
	string name = "hISO" + algo_desc_short[algo] + "_" + 
	  cuts_desc_short[cut];
	string title = algo_desc_long[algo] + " muon isol with " + 
	  cuts_desc_long[cut];
	hISO[cut][algo] = new TH1F(name.c_str(), title.c_str(), nBinIsoMu,
				   minIsoMu,maxIsoMu);
      }
  
}//---------defineMuonPhiHistos()


//--------------------------------------------------------------
void defineHistos_MuonChargePt()
{
//--------------------------------------------------------------
#if debugme
  cout << " define histos MuonChargePt " << endl;
#endif

  for(int algo = MuAlgo_MIN; algo <= MuAlgo_MAX; ++algo)
    {
      string name = "hPT" + algo_desc_short[algo] + "_" + "plus";
      string title = " (+) " + algo_desc_long[algo] + " muon p_{T} qual ";
      hPTplus[algo] = new TH1F(name.c_str(), title.c_str(), nBinPtMu,
			       minPtMu,maxPtMu);
      name = "hPT" + algo_desc_short[algo] + "_" + "minus";
      title = " (-) " + algo_desc_long[algo] + " muon p_{T} qual ";
      hPTminus[algo] = new TH1F(name.c_str(), title.c_str(), nBinPtMu,
				minPtMu,maxPtMu);
    }
    
  return;

}//------defineHistos_MuonChargePt



// Define histograms depending on the type of study
//--------------------------------------------------------------
void defineHistos(int option)
{
//--------------------------------------------------------------

#if debugme
    cout << " define histos with option = " << option<<endl;
#endif

    if(option == 1)
      {
        defineHistos_MuonPt();
        defineHistos_MuonEta();
        defineHistos_MuonPhi();
        defineHistos_MuonJetDPhi();
        defineHistos_MuonIso();
        defineHistos_TMass();
      }
    else if(option == 2) 
      defineHistos_MuonChargePt();
    else 
    {
      cout << " WARNING!!! No histos were defined for option = " << option  
	   << endl;
      return;
    }

}//----defineHistos()


//tabulate results after the cut has been passed
//TODO: Improve flexibility
//-----------------------------------------------------------
void tabulateMe(int Num_surv_cut[][Num_trkAlgos][Num_flavors], int& cut_index, 
                float weight, const wprime::Event * ev,
                const wprime::Muon* mu,
                const bool fill_entry[][Num_trkAlgos], int option,
                bool accountme[][Num_trkAlgos][Num_flavors], bool isW,TH1D* h,
                TH2D* hh )
{
//-----------------------------------------------------------
#if debugme
  cout << " Tabulating results for cut_index = " << cut_index << endl;
#endif

  //if the accountme switch is on,
  //increase the number of events passing the cuts
  //and turn the switch off so we don't count more than once per event
  for(int f = 0; f != Num_flavors; ++f)
    for(int mual = MuAlgo_MIN; mual <= MuAlgo_MAX; ++mual)
      if (accountme[cut_index][mual][f] && fill_entry[f][mual]) 
	{
	  ++(Num_surv_cut[cut_index][mual][f]);
	  accountme[cut_index][mual][f] = false;
	}
  
  const TLorentzVector * P[Num_trkAlgos] = {
    &(mu->global.p), &(mu->tracker.p), &(mu->tpfms.p), &(mu->cocktail.p),
    &(mu->picky.p),  &(mu->tmr.p), &(mu->dyt.p)};
  //fill the histograms
  if(option == 1) {

    fillHistos_MuonPt(cut_index,weight,P, fill_entry);
    fillHistos_MuonEta(cut_index,weight,P, fill_entry);
    fillHistos_MuonPhi(cut_index,weight,P, fill_entry);
    fillHistos_MuonJetDPhi(cut_index,weight, ev, P, fill_entry);
    fillHistos_MuonIso(cut_index,weight, mu, fill_entry);
    fillHistos_TMass(cut_index,weight, ev, P, fill_entry, isW, h, hh);
  }
  else if(option == 2) 
    {
      const Int_t charge[Num_trkAlgos] = {
	mu->global.q, mu->tracker.q, mu->tpfms.q, mu->cocktail.q,
	mu->picky.q, mu->tmr.q, mu->dyt.q};
      
      fillHistos_MuonChargePt(cut_index, weight, charge, P, fill_entry);
    }
  else 
    cout << " WARNING!! I can't tabulate anything for option = " << option 
	 << endl;
  
#if debugme
  cout << " cut_index is now = " << cut_index << endl;
#endif

}//-----tabulateMe


//Fill Histograms
//-----------------------------------------------------------
void fillHistos_MuonPt(int index, float weight, const TLorentzVector **p,
                       const bool fill_entry[][Num_trkAlgos])
{
//-----------------------------------------------------------
  for(int algo = MuAlgo_MIN; algo <= MuAlgo_MAX; ++algo)
    {
      if(fill_entry[PT_INDEX][algo])
	hPT[index][algo]->Fill( (p[algo])->Pt(), weight);
    }
   
}//fillHistos_MuonPt


//-----------------------------------------------------------
void fillHistos_MuonEta(int index, float weight, const TLorentzVector **p,
			const bool fill_entry[][Num_trkAlgos])
{
//-----------------------------------------------------------
  for(int algo = MuAlgo_MIN; algo <= MuAlgo_MAX; ++algo)
    {
      if(fill_entry[PT_INDEX][algo])
	hETA[index][algo]->Fill( (p[algo])->Eta(), weight);
    }
}//fillHistos_MuonEta

//-----------------------------------------------------------
void fillHistos_MuonPhi(int index, float weight, const TLorentzVector **p,
			const bool fill_entry[][Num_trkAlgos])
{
//-----------------------------------------------------------
  for(int algo = MuAlgo_MIN; algo <= MuAlgo_MAX; ++algo)
    {
      if(fill_entry[PT_INDEX][algo])
	hPHI[index][algo]->Fill( (p[algo])->Phi(), weight);
    }

}//fillHistos_MuonPhi

//-----------------------------------------------------------
void fillHistos_MuonJetDPhi(int index, float weight, const wprime::Event * ev,
			    const TLorentzVector **p, 
			    const bool fill_entry[][Num_trkAlgos])
{
//-----------------------------------------------------------
  for(int algo = MuAlgo_MIN; algo <= MuAlgo_MAX; ++algo)
    {
      if(fill_entry[PT_INDEX][algo])
	hMJDPHI[index][algo]->Fill( XJetDPhi(*(p[algo]), ev), weight);
    }
  
}//fillHistos_MuonJetDPhi

//-----------------------------------------------------------
void fillHistos_TMass(int index, float weight, const wprime::Event * ev,
		      const TLorentzVector **p, 
                      const bool fill_entry[][Num_trkAlgos], bool isW, TH1D* h,
                TH2D * hh)
{
//-----------------------------------------------------------
    for(int algo = MuAlgo_MIN; algo <= MuAlgo_MAX; ++algo)
    {

        if(fill_entry[MT_INDEX][algo])
            hTM[index][algo]->Fill(TMass( *(p[algo]),
                                          getNewMET(ev,p, *(p[algo]),isW, h,hh))
                                   ,weight);
    }

    
}//fillHistos_TMass


//-----------------------------------------------------------
void fillHistos_MuonIso(int index, float weight, const wprime::Muon* mu,
			const bool fill_entry[][Num_trkAlgos])
{
//-----------------------------------------------------------

  if(!mu) {
#if debugme
    cout << " No muon found in the event " << endl;
#endif
    return;
  }

  for(int algo = MuAlgo_MIN; algo <= MuAlgo_MAX; ++algo){
    if(fill_entry[PT_INDEX][algo])
      hISO[index][algo]->Fill(CombRelIsolation(mu,deltaRIsoIndex),weight);
    }
  
}//fillHistos_MuonIso

//Fill Histograms
//-----------------------------------------------------------
void fillHistos_MuonChargePt(int index, float weight, const Int_t * Q, 
			     const TLorentzVector **p, 
			     const bool fill_entry[][Num_trkAlgos])
{
//-----------------------------------------------------------

  //only fill the required histograms
  if (index < (Num_selection_cuts - Num_histo_sets_chargePt)) return;

  for(int algo = MuAlgo_MIN; algo <= MuAlgo_MAX; ++algo)
    {
      if(fill_entry[PT_INDEX][algo])
	{
	  if(Q[algo] > 0)
	    hPTplus[algo]->Fill( (p[algo])->Pt(), weight);
	  else
	    hPTminus[algo]->Fill( (p[algo])->Pt(), weight);
	}
    }
    
}//---------------fillHistos_MuonChargePt()



//------------------------------------------------------------------------
void saveHistos_MainAnalysis()
{
//------------------------------------------------------------------------
#if debugme
  cout << " Saving MuonPt histos " << endl;
#endif
  for(int i = 0; i != Num_selection_cuts; ++i)
    for(int algo = MuAlgo_MIN; algo <= MuAlgo_MAX; ++algo)
      {
        hPT[i][algo]->Write();
        hETA[i][algo]->Write(); 
        hPHI[i][algo]->Write();  
        hMJDPHI[i][algo]->Write();
        hISO[i][algo]->Write();   
        hTM[i][algo]->Write();
      }
}//-----saveHistos_MainAnalysis()
  
        
//------------------------------------------------------------------------
void saveHistos_MuonChargePt()
{
//------------------------------------------------------------------------
#if debugme
  cout << " Saving MuonChargePt histos " << endl;
#endif
    
  for(int algo = MuAlgo_MIN; algo <= MuAlgo_MAX; ++algo)
    {
      hPTplus[algo]->Write();
      hPTminus[algo]->Write();
    }
  return;
}//-----saveHistos_MuonPtJetIso


//------------------------------------------------------------------------
void saveHistos(TFile * fout, string sample, int option)
{
//------------------------------------------------------------------------

  fout->cd(); 
  fout->mkdir(sample.c_str()); 
  fout->cd(sample.c_str());

  if(option == 1) saveHistos_MainAnalysis();

  else if(option == 2) saveHistos_MuonChargePt();
  else 
    {
      cout << " Something went wrong when saving with option " << option 
	   << " ...Quiting..." << endl;
      abort();
    }
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
  weight = file.weight;
  cout << " Number of events in the file = " << nevents << endl;
}//---------gatherFileBasicInfo()


//Main Analysis study
//---------------------------------------------------------------------------
void GetMuonPtDistribution(wprime::InputFile& file, float & Nexp_evt, 
			   int option, bool highestPtMuonOnly)
{
//---------------------------------------------------------------------------

  //gather file basic info
  wprime::Event * ev = 0;
  int nevents = 0; float weight = 0.0;
  gatherFileBasicInfo(file, ev, nevents, weight);
  
  // counter for (unweighted) events after cuts
  int Num_surv_cut[Num_selection_cuts][Num_trkAlgos][Num_flavors] = {{{0}}};

  // the idea is to be able to swap selection cuts without having to modify both
  // the algo_desc_short array AND the code looping over the cuts 
  // in the very same order
  selection_map cuts;
  setupCutOrder(cuts);
  bool isW = ((file.samplename=="W" || file.samplename=="Wlowpt") && dorecoil);

  //open the Z data file with info about recoil
   TH2D* hh;
   TH1D* h;
   TDirectoryFile *demo;
   TFile* File;
   if (isW){
       File = new TFile("ZMET_data.root", "READONLY");
       demo=(TDirectoryFile*)File->Get("pflow");
       hh = (TH2D*)demo->Get("hMETParalvsVBPt");
       hh->SetName("hh");
       h = (TH1D*)demo->Get("hMETPerp");
       h->SetName("h");
   }
   else{
       h = 0;
       hh = 0;
   }

  // Loop over events 
  for(int i = 0; i != nevents; ++i){ // event loop
#if debugme
    cout << " ########## Processing event #: " << i+1 << endl;
#endif
    
    file.tree->GetEntry(i);
    
    // switch to help us keep track of whether a muon has already
    // been found in current event surviving the ith-cut;
    // this will ensure that we increase Num_surv_cut maximum once per evet
    // whereas we nevertheless fill the histograms 
    // for every muon surviving the i-th cut
    bool accountme[Num_selection_cuts][Num_trkAlgos][Num_flavors];
    for(int f = 0; f != Num_flavors; ++f)
      for(int cut = 0; cut != Num_selection_cuts; ++cut)
	for(int algo = MuAlgo_MIN; algo <= MuAlgo_MAX; ++algo)
	  accountme[cut][algo][f] = true;
    
    int iMuMin = 0; int iMuMax = ev->mu->GetLast() + 1;
    
    int iMu = -1;
    if(highestPtMuonOnly)
      {
	GetTheHardestMuon(ev, iMu);
	if(iMu >=0) // there is at least one muon in the event
	  {
	    iMuMin = iMu;
	    iMuMax = iMu + 1;
	  }
      }

    wprime::Muon* theMu = 0;
    //loop over muons
    for (int mi = iMuMin; mi != iMuMax; ++mi){//loop over muons
      
#if debugme
      cout << " ########## Processing muon #: " << mi+1 << endl;
#endif
      //get the muon
      theMu = (wprime::Muon *) ev->mu->At(mi);
      
      bool fill_entry[Num_flavors][Num_trkAlgos];
      for(int f = 0; f != Num_flavors; ++f)
	for(int algo = MuAlgo_MIN; algo <= MuAlgo_MAX; ++algo)
	  fill_entry[f][algo] = true;

#if dumpHighPtMuons
      bool HighPtMuon = false;
#endif
      
      for(int cut_index = 0; cut_index != Num_selection_cuts; ++cut_index)
	{ // loop over selection cuts
	  string arg = cuts_desc_short[cut_index];
	  // call to function [as implemented in setupCutOder]
	  bool * goodPt = &fill_entry[PT_INDEX][0];
	  bool * goodMt = &fill_entry[MT_INDEX][0];
	  bool survived_cut=(cuts[arg])(ev, theMu, goodPt, goodMt,isW,h,hh);
	  if(!survived_cut)break; // skip rest of selection cuts
      tabulateMe(Num_surv_cut, cut_index, weight, ev, theMu, 
                 fill_entry, option, accountme, isW, h, hh);
	  
#if dumpHighPtMuons
	  if(file.samplename=="data" && cut_index == Num_selection_cuts-1)
	    {
	      for(int algo = MuAlgo_MIN; algo <= MuAlgo_MAX; ++algo)
		if(fill_entry[PT_INDEX][algo])
		  {
		    HighPtMuon = true;
		    break;
		  }
	      
	      //	      if(HighPtMuon || ev->run_no == 143657)
	      if(HighPtMuon)
		printHighPtMuon(ev, theMu, fill_entry);
	    }
#endif      
	  
	} // loop over selection cuts

    } //muon loop
    
  } //event loop
  
  //Number of expected events for each cut (weighted)
  for(int f = 0; f != Num_flavors; ++f)
    for(int cut = 0; cut != Num_selection_cuts; ++cut)
      for(int algo = MuAlgo_MIN; algo <= MuAlgo_MAX; ++algo)
	file.Nexp_evt_cut[cut][algo][f] += Num_surv_cut[cut][algo][f] * weight;

  //total # of events (before any cuts)
  Nexp_evt += nevents * weight; 
  
  delete ev;

}//-----------GetMuonPtDistribution()      

void printMe(const string & desc, const wprime::Muon* theMu)
{
  const wprime::Track * Tk[Num_trkAlgos] = {
    &(theMu->global), &(theMu->tracker), &(theMu->tpfms), 
    &(theMu->cocktail), &(theMu->picky),  &(theMu->tmr), &(theMu->dyt)};
  
  for(int algo = MuAlgo_MIN; algo <= MuAlgo_MAX; ++algo)
    {
      if(desc == "  pt = ")
	cout << " " <<algo_desc_short[algo] << desc << (Tk[algo])->p.Pt();
      else if(desc == " dpt = ")
	cout << " " << algo_desc_short[algo] << desc << (Tk[algo])->dpt;
      else if(desc == "layers")
	{
	  cout << " " << algo_desc_short[algo] << ": (" 
	       << (Tk[algo])->Nstrip_layer << ", " 
	       << (Tk[algo])->Npixel_layer << ") ";	  
	}
      else if(desc == "layersNoMeas")
	{
	  cout << " " << algo_desc_short[algo] << ": (" 
	       << (Tk[algo])->Nstrip_layerNoMeas << ", " 
	       << (Tk[algo])->Npixel_layerNoMeas << ") ";	  
	}
      else if(desc == "hits")
	{
	  cout << " " << algo_desc_short[algo] << ": (" 
	       << (Tk[algo])->NsiStrip_hits << ", " 
	       << (Tk[algo])->Npixel_hits 
	       << ", " << (Tk[algo])->Nmuon_hits << ") ";
	}
      else if(desc == "chi2")
	{
	  cout << " " << algo_desc_short[algo] << ": " 
	       << (Tk[algo])->chi2 << "/" << (Tk[algo])->ndof;
	}
      else
	{
	  cout << " Oops! Do not understand desc = " << desc << endl;
	  abort();
	}
    }
  cout << endl;
}

void printHighPtMuon(const wprime::Event * ev, const wprime::Muon* theMu,
                     const bool fill_entry[][Num_trkAlgos],bool isW, TH1D* h,
                TH2D * hh)
{
  cout << " Run # = " << ev->run_no << " Event # = " 
       << ev->evt_no << " LS = " << ev->LS_no << endl;

  const TLorentzVector * P[Num_trkAlgos] = {
      &(theMu->global.p), &(theMu->tracker.p), &(theMu->tpfms.p), &(theMu->cocktail.p),
      &(theMu->picky.p),  &(theMu->tmr.p), &(theMu->dyt.p)};


  cout << " default pfMET = " << ev->pfmet.Mod() << endl;
  TVector2 met = getNewMET(ev, P,theMu->global.p,isW,h,hh);
  cout << " pfMET for global = " << met.Mod() << endl;
  cout << " M_T for global = " << TMass(theMu->global.p, met) << endl;
  met = getNewMET(ev,P,theMu->cocktail.p,isW,h,hh);
  cout << " pfMET for cocktail = " << met.Mod() << endl;
  cout << " M_T for cocktail = " << TMass(theMu->cocktail.p, met) << endl;
  printMe("  pt = ", theMu);
  printMe(" dpt = ", theMu);
  cout << " # of layers: (strip, pixel) " << endl;
  printMe("layers", theMu);
  cout << " # of layers w/o measurement: (strip, pixel) " << endl;
  printMe("layersNoMeas", theMu);
  cout << " # of hits (strip, pixel, muon) " << endl;
  printMe("hits", theMu);
  cout << " Chi2/Ndof " << endl;
  printMe("chi2", theMu);
  cout << " Tracker eta =  " << theMu->tracker.p.Eta()
       << ", Tracker phi = " << theMu->tracker.p.Phi() << endl;
  cout << " # of standalone muon hits = " << theMu->Nmu_hits << endl;
  cout << " Global: " << theMu->AllGlobalMuons
       << " Tracker: " << theMu->AllTrackerMuons
       << " Standalone: " << theMu->AllStandAloneMuons
       << " Global prompt tight: " << theMu->GlobalMuonPromptTight << endl;

  cout << " Survives all cuts? ";
  for(int algo = MuAlgo_MIN; algo <= MuAlgo_MAX; ++algo)
    cout << algo_desc_short[algo] << ": PT = "<< fill_entry[PT_INDEX][algo]
	 << ": MT = " << fill_entry[MT_INDEX][algo]
	 << " ";
  cout << endl << endl;
  
}


//---------------------------------------------------------------------------
void GetDistributionGeneric(wprime::InputFile & file, 
                            TFile *fout, ofstream & out, 
                            const int option = 1, 
			    const bool highestPtMuonOnly = true)
{
//---------------------------------------------------------------------------
#if debugme
  cout << " $$$$$$$$$$$$$$$$$$$$GetDistributionGeneric " << endl;
#endif

  // Define histograms according to the type of study
  defineHistos(option);
  
  //initialize expected-event counter
  float Nexp_evt = 0;

  if(!file.tree)
    return;

  cout << " Processing sample " << file.description << endl;
  
  if (option == 1 || option == 2) 
    GetMuonPtDistribution(file, Nexp_evt,option, highestPtMuonOnly);
  else 
    {
      cout << " Option " << option << " is invalid, quiting..." << endl;
      abort();
    }
  
  // Print the results if needed according to study case
  if (option == 1 || option == 2) 
    printSummary_MuonPt(out, file, Nexp_evt);
  else  
    cout << " Nothing to print for this study " << endl;
  
  //save histograms according to study case
  saveHistos(fout, file.samplename, option);  
}





