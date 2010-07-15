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

#include <vector>
#include <string>
#include "TFile.h"
#include "TH1F.h"
#include "TTree.h"
#include <stdlib.h>
#include <string.h>

//gen pt
float genpt = 500;
//track maps
const reco::TrackToTrackMap * tevMap;
const reco::TrackToTrackMap * tevMap_new;
const reco::TrackToTrackMap * tevMap_anew;



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
    float FindPtValueInTrackMap(const reco::MuonRef& muRef,const reco::TrackToTrackMap* trkMap);

      // ----------member data ---------------------------
     edm::InputTag theGLBMuonLabel;
    TFile *hfile;
    //TTree *mytree;
    TH1F* hresol;
    TH1F* hresol_anew;
    TH1F* hresol_new;
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
MuonPtPlots::MuonPtPlots(const edm::ParameterSet& iConfig)

{
   //now do what ever initialization is needed

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
   Handle<reco::TrackToTrackMap> Htevm,Htevnewm,Htevanewm;

   //Get the original muon collection
   iEvent.getByLabel("muons",Hmu);
   //Get the collections and the products for the refitted ones
   iEvent.getByLabel("tevMuons","firstHit", Htevm);
   tevMap = Htevm.product();
   iEvent.getByLabel("tevMuonsNew","refit1", Htevnewm);
   tevMap_new = Htevnewm.product();
   iEvent.getByLabel("tevMuonsAltNew", "refit1", Htevanewm);
   tevMap_anew = Htevanewm.product();

   
   //Loop over muon collection but take only global muons
   for (unsigned it = 0; it!= Hmu->size(); ++it){
       reco::MuonRef muRef(Hmu,it);
       if(!(muRef->isGlobalMuon())) continue;
       float tev_pt = FindPtValueInTrackMap(muRef,tevMap);;
       float tevnew_pt = FindPtValueInTrackMap(muRef,tevMap_new);
       float tevanew_pt = FindPtValueInTrackMap(muRef,tevMap_anew);

       //make resolution plot
       if (tev_pt!=9999) {
           hresol->Fill(genpt*((1/tev_pt) - (1/genpt)));
           FillXYplot(hxy,muRef,tevMap);
       }
       if (tevanew_pt!=9999) {
           hresol_anew->Fill(genpt*((1/tevanew_pt) - (1/genpt)));
           FillXYplot(hxy_anew,muRef,tevMap);
       }
       else{
           FillXYplot(hxy_anew_rej,muRef,tevMap);
       }
       if (tevnew_pt!=9999) {
           hresol_new->Fill(genpt*((1/tevnew_pt) - (1/genpt)));
           FillXYplot(hxy_new,muRef,tevMap);
       }
       else{
           FillXYplot(hxy_new_rej,muRef,tevMap);
       }
       
   }//loop over original tev first hit collection

}//======analyze


// ------------ method called once each job just before starting event loop  ------------
void 
MuonPtPlots::beginJob()
{

//std::cout<<"Beginjob"<<std::endl;
    hfile = new TFile("muonplots.root","RECREATE");
    hfile->cd();

    hresol = new TH1F("hresol",";Relative 1/p_{T} resolution;Events",100,-0.5,0.5);
    hresol_anew = new TH1F("hresol_anew",";Relative 1/p_{T} resolution;Events",100,-0.5,0.5);
    hresol_new = new TH1F("hresol_new",";Relative 1/p_{T} resolution;Events",100,-0.5,0.5);

}

// ------------ method called once each job just after ending the event loop  ------------
void 
MuonPtPlots::endJob() {
    hfile->cd();
    hresol->Write();
    hresol_anew->Write();
    hresol_new->Write();
    hfile->Close();
}


//====================Find the pt value for a track in other collection
float MuonPtPlots::FindPtValueInTrackMap(const reco::MuonRef& muRef,const reco::TrackToTrackMap* trkMap){

    reco::TrackToTrackMap::const_iterator trkit;
    trkit = trkMap->find(muRef->globalTrack());

    if (trkit == trkMap->end()){
        std::cout<<"Invalid track in this collection"<<std::endl;
        return -9999;
    }

    return (*(trkit->val)).pt();
    

}//=========FindPtValueInTrackMap



//====================Find the x and y positions of every hit in the track
//and plot
void MuonPtPlots::FillXYplot(TH2F* hh,const reco::MuonRef& muRef,const reco::TrackToTrackMap* trkMap){

    reco::TrackToTrackMap::const_iterator trkit;
    trkit = trkMap->find(muRef->globalTrack());

    if (trkit == trkMap->end()){
        std::cout<<"Invalid track in this collection"<<std::endl;
        return;
    }

    //the track
    const reco::Track glbTrack = *(trkit->val);
    //the magnetic field
    edm::ESHandle<MagneticField> theMF;
    iSetup.get<IdealMagneticFieldRecord>().get(magneticField);

    //the tracker geometry
    
    reco::TransientTrack ttrack(glbTrack,,);
    return;
    
    

}//=========FindPtValueInTrackMap


//define this as a plug-in
DEFINE_FWK_MODULE(MuonPtPlots);
