// -*- C++ -*-
//
// Package:    ZmumuTreeMaker
// Class:      ZmumuTreeMaker
// 
/**\class ZmumuTreeMaker ZmumuTreeMaker.cc BUSeminar/ZmumuTreeMaker/src/ZmumuTreeMaker.cc

 Description: <one line class summary>

 Implementation:
     <Notes on implementation>
*/
//
// Original Author:  Edgar Fernando Carrera Jarrin
//         Created:  Tue Feb  9 11:17:17 CET 2010
// $Id$
//
//


// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "DataFormats/Common/interface/Handle.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"

#include "FWCore/Framework/interface/TriggerNames.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/MuonReco/interface/MuonFwd.h"
#include "DataFormats/METReco/interface/CaloMET.h"
#include "DataFormats/METReco/interface/MET.h"
#include "DataFormats/METReco/interface/PFMET.h"
#include "DataFormats/MuonReco/interface/MuonMETCorrectionData.h"

#include <cmath>
#include "TMath.h"
#include "HLTrigger/HLTanalyzers/interface/JetUtil.h"


#include<vector>
#include<string>
#include "TFile.h"
#include "TTree.h"
#include <stdlib.h>
#include <string.h>

//
// class decleration
//

class ZmumuTreeMaker : public edm::EDAnalyzer {
   public:
      explicit ZmumuTreeMaker(const edm::ParameterSet&);
      ~ZmumuTreeMaker();


   private:
      virtual void beginJob() ;
      virtual void analyze(const edm::Event&, const edm::EventSetup&);
      virtual void endJob() ;

      // ----------member data ---------------------------

    edm::InputTag hlTriggerResults_;  // Input tag for TriggerResults
    edm::TriggerNames triggerNames_;  // TriggerNames class
    bool _Debug;

    TFile *hfile;
    TTree *mytree;
  
    int runno;
    int evtno;
    int lumiblock;
    int trg_L1MuOpen;
    int trg_Mu9;
    int trg_DoubleMu3;

    int nmu;
     float mu_eta[200];
     float mu_phi[200];
     float mu_e[200];
     float mu_et[200];
     float mu_pt[200];
     float mu_px[200];
     float mu_py[200];
     float mu_pz[200];
    int mu_charge[200];
    
     float met;
    float met_phi;

    float tmass;
    float imass;

};

//
// constants, enums and typedefs
//

//
// static data member definitions
//

//
// constructors and destructor
//
ZmumuTreeMaker::ZmumuTreeMaker(const edm::ParameterSet& iConfig):
  hlTriggerResults_ (iConfig.getParameter<edm::InputTag> ("HLTriggerResults")),
  triggerNames_()
{
    _Debug=false;

    //std::cout<<"HL TriggerResults:"<<std::endl;
    //now do what ever initialization is needed
    LogDebug("HLTrigReport") << "HL TiggerResults: " + hlTriggerResults_.encode();


}


ZmumuTreeMaker::~ZmumuTreeMaker()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called to for each event  ------------
void
ZmumuTreeMaker::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{

    using namespace std;
    using namespace edm;

    runno = iEvent.id().run();
    evtno  = iEvent.id().event();
    lumiblock = iEvent.luminosityBlock();


    //get HLT decisions
    Handle<TriggerResults> hltresults;
    iEvent.getByLabel(hlTriggerResults_,hltresults);
    if (hltresults.isValid()) {
        int ntrigs = hltresults->size();
        if (ntrigs==0){std::cout << "%HLTInfo -- No trigger name given in TriggerResults of the input " << std::endl;}
        
        triggerNames_.init(* hltresults);
        
        trg_L1MuOpen = 0;
        trg_Mu9 = 0;
        trg_DoubleMu3 = 0;
        
        for (int itrig = 0; itrig != ntrigs; ++itrig){
            string trigName=triggerNames_.triggerName(itrig);
            bool accept = hltresults->accept(itrig);
            if (itrig == 18 && accept) trg_L1MuOpen = 1;
            if (itrig == 26 && accept) trg_Mu9 = 1;
            if (itrig == 29 && accept) trg_DoubleMu3 = 1;
        }
    }
    else { if (_Debug) std::cout << "%HLTInfo -- No Trigger Result" << std::endl;
    }
    

    //get muons
    Handle<reco::MuonCollection> Muon;
    iEvent.getByLabel("muons",Muon);
    if (Muon.isValid()) {
        reco::MuonCollection mymuons;
        mymuons = * Muon;
        std::sort(mymuons.begin(),mymuons.end(),PtGreater());
        nmu = mymuons.size();
        typedef reco::MuonCollection::const_iterator muiter;
        int imu=0;
        for (muiter i=mymuons.begin(); i!=mymuons.end(); i++) {
            mu_pt[imu] = i->pt();
            mu_phi[imu] = i->phi();
            mu_eta[imu] = i->eta();
            mu_et[imu] = i->et();
            mu_e[imu] = i->energy();
            mu_pt[imu] = i->pt();
            mu_px[imu] = i->px();
            mu_py[imu] = i->py();
            mu_pz[imu] = i->pz();
            mu_charge[imu] = i->charge();
            imu++;
        }
    }
    else {nmu = 0;}


    //get met
    Handle< edm::View<reco::CaloMET> > muCorrMEThandle;
    iEvent.getByLabel("corMetGlobalMuons", muCorrMEThandle);
    met = (muCorrMEThandle->front() ).et();
    met_phi = (muCorrMEThandle->front() ).phi();
    



    //transverse mass
    float tmass_sqr = 2*mu_et[0]*met*(1-cos(mu_phi[0]-met_phi));
    tmass = (tmass_sqr>0) ? sqrt(tmass_sqr) : 0;



    //inv mass
    float Ene = mu_e[0] + mu_e[1];
    float Px = mu_px[0] + mu_px[1];
    float Py = mu_py[0] + mu_py[1];
    float Pz = mu_pz[0] + mu_pz[1];
    float imass_sqr = Ene*Ene - Px*Px - Py*Py - Pz*Pz;
    imass = (imass_sqr>0) ? sqrt(imass_sqr) : 0;
    

    mytree->Fill();
    return;

}


// ------------ method called once each job just before starting event loop  ------------
void 
ZmumuTreeMaker::beginJob()
{
    //std::cout<<"Beginjob"<<std::endl;
    hfile = new TFile("ZtoMuMuTree.root","RECREATE");
    mytree = new TTree("mytree","");
    
    //global
    mytree->Branch("runno",&runno,"runno/I");
    mytree->Branch("evtno",&evtno,"evtno/I");
    mytree->Branch("lumiblock",&lumiblock,"lumiblock/I");
    mytree->Branch("trg_L1MuOpen",&trg_L1MuOpen,"trg_L1MuOpen/I");
    mytree->Branch("trg_Mu9",&trg_Mu9,"trg_Mu9/I");
    mytree->Branch("trg_DoubleMu3",&trg_DoubleMu3,"trg_DoubleMu3/I");
    //Muons
     mytree->Branch("nmu",&nmu,"nmu/I");
     mytree->Branch("mu_eta",mu_eta,"mu_eta[nmu]/F");
     mytree->Branch("mu_phi",mu_phi,"mu_phi[nmu]/F");
     mytree->Branch("mu_e",mu_e,"mu_e[nmu]/F");
     mytree->Branch("mu_et",mu_et,"mu_et[nmu]/F");
     mytree->Branch("mu_pt",mu_pt,"mu_pt[nmu]/F");
     mytree->Branch("mu_px",mu_px,"mu_px[nmu]/F");
     mytree->Branch("mu_py",mu_py,"mu_py[nmu]/F");
     mytree->Branch("mu_pz",mu_pz,"mu_pz[nmu]/F");
     mytree->Branch("mu_charge",mu_charge,"mu_charge[nmu]/I");
     //met
     mytree->Branch("met",&met,"met/F");
     mytree->Branch("met_phi",&met_phi,"met_phi/F");
     //masses
     mytree->Branch("tmass",&tmass,"tmass/F");
     mytree->Branch("imass",&imass,"imass/F");

    

}

// ------------ method called once each job just after ending the event loop  ------------
void 
ZmumuTreeMaker::endJob() {
    //std::cout<<"endJob"<<std::endl;
    hfile->Write();

}

//define this as a plug-in
DEFINE_FWK_MODULE(ZmumuTreeMaker);
