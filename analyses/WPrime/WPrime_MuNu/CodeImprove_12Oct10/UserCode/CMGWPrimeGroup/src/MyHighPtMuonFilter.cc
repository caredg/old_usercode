// -*- C++ -*-
//
// Package:    MyHighPtMuonFilter
// Class:      MyHighPtMuonFilter
// 
/**\class MyHighPtMuonFilter MyHighPtMuonFilter.cc WPrime/MyHighPtMuonFilter/src/MyHighPtMuonFilter.cc

 Description: <one line class summary>

 Implementation:
     <Notes on implementation>
*/
//
// Original Author:  Alessio GHEZZI
//         Created:  Thu Aug  7 14:48:54 CEST 2008
// $Id: MyHighPtMuonFilter.cc,v 1.1 2009/01/23 10:46:34 cleonido Exp $
//
//


// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDFilter.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "FWCore/Utilities/interface/InputTag.h"
#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/MuonReco/interface/MuonFwd.h"
#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"

//
// class declaration
//

class MyHighPtMuonFilter : public edm::EDFilter {
   public:
      explicit MyHighPtMuonFilter(const edm::ParameterSet&);
      ~MyHighPtMuonFilter();

   private:
      virtual void beginJob(const edm::EventSetup&) ;
      virtual bool filter(edm::Event&, const edm::EventSetup&);
      virtual void endJob() ;
      
      // ----------member data ---------------------------
  edm::InputTag muonTag_;
  double pt_th_;
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
MyHighPtMuonFilter::MyHighPtMuonFilter(const edm::ParameterSet& iConfig):
  muonTag_(iConfig.getParameter<edm::InputTag> ("MuonTag")),
  pt_th_(iConfig.getParameter<double>("PtTh"))
{
   //now do what ever initialization is needed

}


MyHighPtMuonFilter::~MyHighPtMuonFilter()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called on each new Event  ------------
bool
MyHighPtMuonFilter::filter(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  
   using namespace edm;
   using namespace reco;
   Handle<reco::TrackCollection> muonCollection;
   iEvent.getByLabel(muonTag_, muonCollection);
   
   for (unsigned int i=0; i<muonCollection->size(); i++) {
     TrackRef mu(muonCollection,i);
     if( mu->pt() >  pt_th_ ){return true;}
   }
   
   return false;
}

// ------------ method called once each job just before starting event loop  ------------
void 
MyHighPtMuonFilter::beginJob(const edm::EventSetup&)
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
MyHighPtMuonFilter::endJob() {
}

//define this as a plug-in
DEFINE_FWK_MODULE(MyHighPtMuonFilter);
