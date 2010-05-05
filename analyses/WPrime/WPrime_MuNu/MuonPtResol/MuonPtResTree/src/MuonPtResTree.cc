// -*- C++ -*-
//
// Package:    MuonPtResTree
// Class:      MuonPtResTree
// 
/**\class MuonPtResTree MuonPtResTree.cc MuonPtResol/MuonPtResTree/src/MuonPtResTree.cc

 Description: <one line class summary>

 Implementation:
     <Notes on implementation>
*/
//
// Original Author:  Edgar Fernando Carrera Jarrin,40 3-A01,+41227671566,
//         Created:  Tue May  4 17:31:48 CEST 2010
// $Id$
//
//

#include "MuonPtResol/MuonPtResTree/interface/MuonPtResTree.h"
// system include files
#include <memory>
#include <cmath>

#include "FWCore/Framework/interface/MakerMacros.h"

#include "DataFormats/Common/interface/TriggerResults.h"

#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/JetReco/interface/CaloJet.h"
#include "DataFormats/MuonReco/interface/MuonCocktails.h"

#include "DataFormats/METReco/interface/PFMET.h"
#include "DataFormats/METReco/interface/PFMETCollection.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TFile.h"
#include "TTree.h"

#include <iostream>
#include <fstream>
#include <map>
using std::cout; using std::endl; using std::string; using std::ifstream;
//
// class decleration
//



//
// constants, enums and typedefs
//

//
// static data member definitions
//

//
// constructors and destructor
//
MuonPtResTree::MuonPtResTree(const edm::ParameterSet& iConfig):
    muonTag_(iConfig.getParameter<edm::InputTag> ("MuonTag"))

{
   //now do what ever initialization is needed
    evt = new muresol::Event(); 
}


MuonPtResTree::~MuonPtResTree()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)
    if(evt) delete evt;
}


//
// member functions
//

// initialize histograms
void MuonPtResTree::init_histograms()
{
  string tree_title = "Muon pt resolution study";
  tree_event = fs->make<TTree>("muptres", tree_title.c_str());
  tree_event->Branch("refit", "muresol::Event", &evt, 8000, 2);
}

// initialize event info
void MuonPtResTree::init_event()
{
    gen_muons.clear();
    N_muons = N_all_muons = 0;

    evt->mu_mc->Clear(); 
    evt->mu->Clear();

    evt->HLT_L1MuOpen = evt->HLT_L1Mu = evt->HLT_Mu3 = evt->HLT_Mu5 =
        evt->HLT_Mu9 = false;
}//----init_event()



// initialize job info
void MuonPtResTree::init_run()
{
  realData = false;
  init_histograms();
}


// get the generator info, populate gen_muons, set genmu_acceptance flag
void MuonPtResTree::getGenParticles(const edm::Event & iEvent)
{
    if(realData)return;
    
    TClonesArray & mcmu = *(evt->mu_mc);
    
    int Nm = 0; 
    
    edm::Handle<reco::GenParticleCollection> genParticles;
    iEvent.getByLabel("genParticles", genParticles);
    for(size_t i = 0; i != genParticles->size(); ++i) {
        const GenParticle & p = (*genParticles)[i];
        int id = p.pdgId(); 
        if( abs(id)!=13)
            continue; // keep only muons
        
        int st = p.status(); int ch = p.charge();
        int momid = -9999;
        
        TLorentzVector p4(p.px(), p.py(), p.pz(),  p.energy());
        const Candidate * mom = p.mother();
        if(mom)
            momid = mom->pdgId();
        
        switch(abs(id))
        {
            case 13: // muon
                new(mcmu[Nm]) muresol::MCParticle(p4, ch, momid, st);
                ++Nm;
                break;
                
            default:
            {
                ; // do nothing for now
            }
            
        } // switch (abs(id))
        
        
    } // loop over genParticles
    
}



// get TeV muons
void MuonPtResTree::getTeVMuons(const edm::Event & iEvent)
{
  edm::Handle<reco::TrackToTrackMap> tevMapH_default;
  edm::Handle<reco::TrackToTrackMap> tevMapH_1stHit;
  edm::Handle<reco::TrackToTrackMap> tevMapH_picky;
  edm::Handle<reco::TrackToTrackMap> tevMapH_oddhits;

  iEvent.getByLabel("tevMuons", "default", tevMapH_default);
  tevMap_default = tevMapH_default.product();

  iEvent.getByLabel("tevMuons", "firstHit", tevMapH_1stHit);
  tevMap_1stHit = tevMapH_1stHit.product();

  iEvent.getByLabel("tevMuons", "picky", tevMapH_picky);
  tevMap_picky = tevMapH_picky.product();

  iEvent.getByLabel("tevTracksOdd1", "odd", tevMapH_oddhits);
  tevMap_oddhits = tevMapH_oddhits.product();
}





// get muons
void MuonPtResTree::getMuons(const edm::Event & iEvent)
{
  // Get the Muon Track collection from the event
  iEvent.getByLabel(muonTag_, muonCollection);
  // # of reconstructed muons in event (including standalone-only)
  N_all_muons = muonCollection->size();
}


// copy tracking info from reco::Track to muresol::Track
void MuonPtResTree::getTracking(muresol::Track & track, const reco::Track & p)
{
  TVector3 p3(p.px(), p.py(), p.pz());
  track.p.SetVectM(p3, muresol::MUON_MASS);
  track.q = p.charge();
  track.chi2 = p.chi2();
  track.dpt = p.ptError();
  track.dq_over_p = p.qoverpError();
  track.ndof = int(p.ndof());
  track.Ntot_hits = p.numberOfValidHits();
  track.Ntrk_hits = p.hitPattern().numberOfValidTrackerHits();
}




// do muon analysis
void MuonPtResTree::doMuons()
{
  TClonesArray & recomu = *(evt->mu);
  cout<<"flag1"<<endl;
  for(unsigned i = 0; i != N_all_muons; ++i) 
    { // loop over reco muons 
        MuonRef mu(muonCollection, i);
        if(!(mu->isGlobalMuon()) )
            continue; // keep only global muons
        new(recomu[N_muons]) muresol::Muon(); 
        muresol::Muon * wpmu = (muresol::Muon *) recomu[N_muons];
        wpmu->Nmu_hits = mu->standAloneMuon()->recHitsSize(); 
        getTracking(wpmu->tracker, *(mu->track()));
        getTracking(wpmu->global, *(mu->combinedMuon()) );
        ++N_muons;
//      doIsolation(mu, wpmu);

      doTeVanalysis(mu, wpmu);
      
    } // loop over reco muons


  
}

// do TeV-muon analysis
void MuonPtResTree::doTeVanalysis(reco::MuonRef mu, muresol::Muon * wpmu)
{
  if(!(mu->isGlobalMuon()) ) return; // keep only global muons

  TrackToTrackMap::const_iterator iTeV_default;
  TrackToTrackMap::const_iterator iTeV_1stHit;
  TrackToTrackMap::const_iterator iTeV_picky;
  TrackToTrackMap::const_iterator iTeV_oddhits;


      iTeV_default = tevMap_default->find(mu->globalTrack());
      iTeV_1stHit = tevMap_1stHit->find(mu->globalTrack());
      iTeV_picky = tevMap_picky->find(mu->globalTrack());
      iTeV_oddhits = tevMap_oddhits->find(mu->globalTrack());

  if(iTeV_default == tevMap_default->end() 
     || iTeV_1stHit == tevMap_1stHit->end() 
     || iTeV_picky == tevMap_picky->end() 
     || iTeV_oddhits == tevMap_oddhits->end())
    {
      cout << "-Muresol_muonreco- Warning: No Tev muons found for this event !! "
	   << endl; 
      return;
    }
  
  getTracking(wpmu->tev_1st, *(iTeV_1stHit->val) );
  getTracking(wpmu->tev_oddhits, *(iTeV_oddhits->val) );

}





// ------------ method called to for each event  ------------
void MuonPtResTree::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
    init_event();
    evt->evt_no = iEvent.id().event();
    evt->run_no = iEvent.id().run();
   
    realData = iEvent.isRealData();
    
    getGenParticles(iEvent);
    getMuons(iEvent);
    getTeVMuons(iEvent);
    doMuons();

    tree_event->Fill();
}


// ------------ method called once each job just before starting event loop  ------------
void 
MuonPtResTree::beginJob()
{

     init_run();

}

// ------------ method called once each job just after ending the event loop  ------------
void 
MuonPtResTree::endJob() {
}

//define this as a plug-in
DEFINE_FWK_MODULE(MuonPtResTree);
