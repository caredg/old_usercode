// -*- C++ -*-
//
// Package:    Wprime_muonreco
// Class:      Wprime_muonreco
// 
/**\class Wprime_muonreco Wprime_muonreco.cc WPrime/MyAnalysisMuons/src/Wprime_muonreco.cc

 Description: <one line class summary>

 Implementation:
     <Notes on implementation>
*/
//
// Original Author:  Alessio Ghezzi, Christos Leonidopoulos, Silvia Goy Lopez
//
//
#include <vector>
#include <string>

#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h" 
#include "FWCore/Framework/interface/Frameworkfwd.h" 
#include "FWCore/Common/interface/TriggerNames.h"

#include "FWCore/Utilities/interface/InputTag.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "FWCore/Utilities/interface/EDMException.h"

#include "FWCore/MessageLogger/interface/MessageLogger.h"

#include "FWCore/ServiceRegistry/interface/Service.h"

#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/TrackReco/interface/TrackToTrackMap.h"
#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/MuonReco/interface/MuonFwd.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/RecoCandidate/interface/IsoDeposit.h"
#include "DataFormats/RecoCandidate/interface/IsoDepositFwd.h"
#include "DataFormats/BeamSpot/interface/BeamSpot.h"

#include "DataFormats/HepMCCandidate/interface/GenParticle.h"

#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "UserCode/CMGWPrimeGroup/interface/wprimeEvent.h"

//
// forward class declaration
//
class TFile;
class TTree;

namespace Wprime_muonreco_histo
{
 
  // cone size and SumPt for isolation
  static const unsigned nBinCone = wprime::N_CONESIZE;
  static const float minCone = wprime::MIN_CONE;
  static const float maxCone = wprime::MAX_CONE;
}

//
// class declaration
//
class Wprime_muonreco : public edm::EDAnalyzer 
{
  public:
  explicit Wprime_muonreco(const edm::ParameterSet&);
  ~Wprime_muonreco();
  
  
  private:
  bool firstEventInRun;
  virtual void beginJob() ;
  virtual void beginRun(edm::Run const &, edm::EventSetup const &);
  virtual void endRun(edm::Run const &, edm::EventSetup const &);
  virtual void analyze(const edm::Event&, const edm::EventSetup&);
  virtual void endJob() ;
  edm::InputTag pvTag_;
  edm::InputTag pvBSTag_;
  edm::InputTag muonTag_;
  edm::InputTag pfmetTag_;
  edm::InputTag HLTTag_; edm::InputTag L1Tag_;
  //      edm::InputTag isoTag_;
  edm::InputTag caloJetTag_;
  edm::InputTag pfJetTag_;
  edm::InputTag tkIsoMapTag_;
  edm::InputTag ecalIsoMapTag_;
  edm::InputTag hcalIsoMapTag_;

  struct muonTrack {
    // reference to standard tracking
    reco::MuonRef mu;
    // reference to default-TeV, 1st-Hit, picky, dyt and cocktail/optimized tracking
    reco::TrackRef TeVMuons[5];
  };

  typedef std::vector<muonTrack>::const_iterator mIt;

  struct trigEff {
    // key: (muon) trigger name, value: # of accepted events
    std::map<std::string, unsigned> trigger_count;
    // # of events processed
    unsigned Nev; 
    // # of events processed with a single reconstructed muon
    unsigned Nev_1mu; 
    // initialize
    void clear(){trigger_count.clear(); Nev = Nev_1mu = 0;}
    trigEff(){clear();}
  };

  typedef std::map<std::string, unsigned>::const_iterator tIt;

  trigEff genMuTrig; // muon-trigger efficiency wrt generator-muons
  trigEff MuTrig; // muon-trigger efficiency for all processed events

  // muon-detector acceptance
  const float detmu_acceptance;

  // (at least one) generated muon within detector acceptance 
  // (defined as |eta|< detmu_acceptance)
  bool genmu_acceptance;

  // whether this is real-data
  bool realData;

  std::vector<std::string> triggerNames;
  typedef std::vector<std::string>::const_iterator It;

  // generator-level muons
  std::vector<reco::GenParticle> gen_muons;

  // software versions used to produce HLT and RECO
  std::string HLTversion;
  std::string RECOversion;
  static const std::string INVALID_RELEASE;
  std::string sample_description;

  // TTree structures
  wprime::Event * evt;
  wprime::JobInfo * job;
  wprime::RunInfo * run;
  std::string software_version;

  //    Histograms, trees and all that
  TTree *tree_job, *tree_run, *tree_event;

  //    Root output file
  edm::Service<TFileService> fs;


  edm::Handle<reco::VertexCollection> pvCollection;
  edm::Handle<reco::VertexCollection> pvBSCollection;
  edm::Handle<reco::MuonCollection> muonCollection;
  edm::Handle<reco::IsoDepositMap> tkMapH;
  edm::Handle<reco::IsoDepositMap> ecalMapH;
  edm::Handle<reco::IsoDepositMap> hcalMapH;
  edm::Handle<reco::BeamSpot> beamSpotHandle;

  const reco::TrackToTrackMap * tevMap_default;
  const reco::TrackToTrackMap * tevMap_1stHit;
  const reco::TrackToTrackMap * tevMap_picky;
  const reco::TrackToTrackMap * tevMap_dyt;

  // # of trigger paths in HLT configuration
  unsigned N_triggers;

  // # of reconstructed muons per event
  unsigned N_muons;
  // # of all reconstructed muons per event (included standalone-only)
  unsigned N_all_muons;

// get TMR track, from 
// http://cmssw.cvs.cern.ch/cgi-bin/cmssw.cgi/CMSSW/SUSYBSMAnalysis/Zprime2muAnalysis/src/MuonCocktails.cc?revision=1.1&view=markup
  const reco::TrackRef getTMR(const reco::TrackRef& trackerTrack,
			      const reco::TrackRef& fmsTrack,
			      const double cut = 4.0) ;


  // copy tracking info from reco::Track to wprime::Track
  void getTracking(wprime::Track & track, const reco::Track & p);

  // fill in with dummy values when there is no track
  void getNullTracking(wprime::Track & track);

  // initialize run info
  void init_run(const edm::Event& iEvent);

  // initialize event info
  void init_event();

  // initialize histograms
  void init_histograms();

  // check trigger is there (at beginRun)
  void check_trigger(const edm::Event & iEvent);
  
  // print summary info over full job
  void printSummary() const;
  // print summary info for real
  void printSummary2(const trigEff & trig, const std::string & description) const;

  // get the generator info, populate gen_muons, set genmu_acceptance flag
  void getGenParticles(const edm::Event & iEvent);

  // get trigger info, update muTrig/genMuTrig
  void getTriggers(const edm::Event & iEvent);

  // get primary vertex info
  void getPVs(const edm::Event & iEvent);

  double met_x, met_y, met;

  // get ParticleFlow-MET
  void getPFMET(const edm::Event & iEvent);

  // get Jets
  void getJets(const edm::Event & iEvent);

  // get isolation
  void getIsolation(const edm::Event & iEvent);

  // get TeV muons
  void getTeVMuons(const edm::Event & iEvent);

  // get muons, update MET
  void getMuons(const edm::Event & iEvent);

  // do muon analysis
  void doMuons();
  // do TeV-muon analysis
  void doTeVanalysis(reco::MuonRef mu, wprime::Muon * wpmu);

  // do isolation
  void doIsolation(reco::MuonRef mu,  wprime::Muon * wpmu);

  void getBeamSpot(const edm::Event & iEvent);
  double correct_d0(const reco::Track & track);
};
