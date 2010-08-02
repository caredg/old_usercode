// -*- C++ -*-
//
// Package:    MuonPtPlots
// Class:      MuonPtPlots
// 
/**\class MuonPtPlots MuonPtPlots.cc MuonPtResol/MuonPtPlots/src/MuonPtPlots.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Edgar Fernando Carrera Jarrin,40 3-A01,+41227671566,
//         Created:  Sun Jul 11 17:17:56 CEST 2010
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

#include "FWCore/ParameterSet/interface/ParameterSet.h"

//Added
#include "DataFormats/TrackReco/interface/Track.h"
#include "TrackingTools/TransientTrack/interface/TransientTrack.h"
#include "DataFormats/TrajectorySeed/interface/TrajectorySeedCollection.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/TrackReco/interface/TrackToTrackMap.h"
#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/MuonReco/interface/MuonFwd.h"
#include "DataFormats/MuonReco/interface/MuonCocktails.h"
#include "DataFormats/MuonReco/interface/MuonSelectors.h"

#include "Geometry/TrackerGeometryBuilder/interface/TrackerGeometry.h"
#include "RecoMuon/TransientTrackingRecHit/interface/MuonTransientTrackingRecHit.h"
#include "TrackingTools/TrackRefitter/interface/TrackTransformer.h"
#include "DataFormats/DetId/interface/DetId.h"
#include "TrackingTools/Records/interface/TransientRecHitRecord.h"

#include <cmath>
#include "TMath.h"
#include "TROOT.h"
#include <vector>
#include <string>
#include "TFile.h"
#include "TObject.h"
#include "TH1F.h"
#include "TTree.h"
#include <stdlib.h>
#include <iostream>
#include <fstream>
//Names
#include <sstream>


using namespace std;
//
// class declaration
//

class MuonPtPlots : public edm::EDAnalyzer {
   public:
      explicit MuonPtPlots(const edm::ParameterSet&);
      ~MuonPtPlots();
    typedef TransientTrackingRecHit::RecHitContainer RecHitContainer;
    typedef TransientTrackingRecHit::ConstRecHitContainer ConstRecHitContainer;

private:
    virtual void beginJob() ;
    virtual void analyze(const edm::Event&, const edm::EventSetup&);
    virtual void endJob() ;
    reco::Track  FindTrackInTrackMap(const reco::MuonRef& muRef,const reco::TrackToTrackMap* trkMap, bool& isValid);

    //track maps
    const reco::TrackToTrackMap * Maptev;
    const reco::TrackToTrackMap * Mapanew;
    const reco::TrackToTrackMap * Mapeve;
    const reco::TrackToTrackMap * Mapodd;

    string filename;
    float genpt;
      // ----------member data ---------------------------
    TFile *hfile;
    TTree *mytree;
    int run;
    int evt;
    int ls;
    
    int nmu;


    bool tevvalid[200];
    float tevpt[200];
    float tevipt[200];
    int tevhits[200];
    float tevres[200];

    bool anewvalid[200];
    float anewpt[200];
    float anewipt[200];
    int anewhits[200];
    float anewres[200];

    bool oddvalid[200];
    float oddpt[200];
    float oddipt[200];
    int oddhits[200];
    float oddres[200];

    bool evevalid[200];
    float evept[200];
    float eveipt[200];
    int evehits[200];
    float everes[200];

    float mres[200];

    


};



//
// constants, enums and typedefs
//
//class to host information for track 

// static data member definitions
//

//
// constructors and destructor
//
MuonPtPlots::MuonPtPlots(const edm::ParameterSet& iConfig)

{
   //now do what ever initialization is needed
    genpt = iConfig.getParameter<double>("MuonPt");
    filename = iConfig.getParameter<string>("FileName");
}


MuonPtPlots::~MuonPtPlots()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called to for each event  ------------
void
MuonPtPlots::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
   using namespace edm;

   //Handle declarations
   Handle<reco::MuonCollection> Hmu;
   Handle<reco::TrackToTrackMap> Htevm,Hanewm,Hevem,Hoddm;

   //Get the original muon collection
   iEvent.getByLabel("muons",Hmu);
   //Get the collections and the products for the refitted ones
   iEvent.getByLabel("tevMuons","firstHit", Htevm);
   Maptev = Htevm.product();
   iEvent.getByLabel("tevMuonsAltNew", "refit1", Hanewm);
   Mapanew = Hanewm.product();
   iEvent.getByLabel("tevMuonsOddHh","odd", Hoddm);
   Mapodd = Hoddm.product();
   iEvent.getByLabel("tevMuonsEveHh","eve", Hevem);
   Mapeve = Hevem.product();

   //get event information
   run = iEvent.id().run();
   evt = iEvent.id().event();
   ls = iEvent.id().luminosityBlock();

   //Loop over muon collection but take only global muons
   nmu = Hmu->size();
   for (int  it = 0; it!= nmu; ++it){
       reco::MuonRef muRef(Hmu,it);
       bool isTrackValid = false;
       if(!(muRef->isGlobalMuon())) continue;

       reco::Track Trktev = FindTrackInTrackMap(muRef,Maptev,isTrackValid);
       if(isTrackValid){
           tevvalid[it] = true;
           tevpt[it] = Trktev.pt();
           tevipt[it] = 1/tevpt[it];
           tevhits[it] = Trktev.recHitsSize();
           tevres[it] = genpt*(tevipt[it] - 1/genpt); 
       } else {tevvalid[it] = false;}


       reco::Track Trkanew = FindTrackInTrackMap(muRef,Mapanew,isTrackValid);
       if(isTrackValid){
           anewvalid[it] = true;
           anewpt[it] = Trkanew.pt();
           anewipt[it] = 1/anewpt[it];
           anewhits[it] = Trkanew.recHitsSize();
           anewres[it] = genpt*(anewipt[it] - 1/genpt); 
       } else {anewvalid[it] = false;}


       reco::Track Trkodd = FindTrackInTrackMap(muRef,Mapodd,isTrackValid);
       if(isTrackValid){
           oddvalid[it] = true;
           oddpt[it] = Trkodd.pt();
           oddipt[it] = 1/oddpt[it];
           oddhits[it] = Trkodd.recHitsSize();
           oddres[it] = genpt*(oddipt[it] - 1/genpt); 
       } else {oddvalid[it] = false;}


       reco::Track Trkeve = FindTrackInTrackMap(muRef,Mapeve,isTrackValid);
       if(isTrackValid){
           evevalid[it] = true;
           evept[it] = Trkeve.pt();
           eveipt[it] = 1/evept[it];
           evehits[it] = Trkeve.recHitsSize();
           everes[it] = genpt*(eveipt[it] - 1/genpt); 
       } else {evevalid[it] = false;}


       mres[it] = genpt*(eveipt[it] - oddipt[it]);

   }//loop over original tev first hit collection

   mytree->Fill();

}//======analyze


// ------------ method called once each job just before starting event loop  ------------
void 
MuonPtPlots::beginJob()
{

    //std::cout<<"Beginjob"<<std::endl;
    hfile = new TFile(filename.c_str(),"RECREATE");
    hfile->cd();
    mytree = new TTree("mytree","");
    mytree->Branch("run",&run,"run/I");
    mytree->Branch("evt",&evt,"evt/I");
    mytree->Branch("ls",&ls,"ls/I");

    mytree->Branch("nmu",&nmu,"nmu/I");

    mytree->Branch("tevvalid",tevvalid,"tevvalid[nmu]/O");
    mytree->Branch("tevpt",tevpt,"tevpt[nmu]/F");
    mytree->Branch("tevipt",tevipt,"tevipt[nmu]/F");
    mytree->Branch("tevhits",tevhits,"tevhits[nmu]/I");
    mytree->Branch("tevres",tevres,"tevres[nmu]/F");

    mytree->Branch("anewvalid",anewvalid,"anewvalid[nmu]/O");
    mytree->Branch("anewpt",anewpt,"anewpt[nmu]/F");
    mytree->Branch("anewipt",anewipt,"anewipt[nmu]/F");
    mytree->Branch("anewhits",anewhits,"anewhits[nmu]/I");
    mytree->Branch("anewres",anewres,"anewres[nmu]/F");

    mytree->Branch("oddvalid",oddvalid,"oddvalid[nmu]/O");
    mytree->Branch("oddpt",oddpt,"oddpt[nmu]/F");
    mytree->Branch("oddipt",oddipt,"oddipt[nmu]/F");
    mytree->Branch("oddhits",oddhits,"oddhits[nmu]/I");
    mytree->Branch("oddres",oddres,"oddres[nmu]/F");

    mytree->Branch("evevalid",evevalid,"evevalid[nmu]/O");
    mytree->Branch("evept",evept,"evept[nmu]/F");
    mytree->Branch("eveipt",eveipt,"eveipt[nmu]/F");
    mytree->Branch("evehits",evehits,"evehits[nmu]/I");
    mytree->Branch("everes",everes,"everes[nmu]/F");

    mytree->Branch("mres",mres,"mres[nmu]/F");






    
  
    
}

// ------------ method called once each job just after ending the event loop  ------------
void 
MuonPtPlots::endJob() {
    hfile->Write();
    
}


//====================Find the pt value for a track in other collection
reco::Track MuonPtPlots::FindTrackInTrackMap(const reco::MuonRef& muRef,const reco::TrackToTrackMap* trkMap, bool& isValid){

    isValid = false;
    reco::Track* trk = new reco::Track();
    reco::TrackToTrackMap::const_iterator trkit;
    trkit = trkMap->find(muRef->globalTrack());

    if (trkit == trkMap->end()){
        std::cout<<"Invalid track in this collection"<<std::endl;
        return *trk;
    }

    //std::cout<<"return track"<<std::endl;
    isValid = true;
     
    return (*(trkit->val));
    

}//=========FindPtValueInTrackMap






//define this as a plug-in
DEFINE_FWK_MODULE(MuonPtPlots);
