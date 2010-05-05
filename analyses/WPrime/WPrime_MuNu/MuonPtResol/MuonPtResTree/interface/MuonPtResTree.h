// -*- C++ -*-
//
// Package:    MuonPtResTree
// Class:      MuonPtResTree
// 
// 
//
//
#include <vector>
#include <string>

#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h" 
#include "FWCore/Framework/interface/TriggerNames.h"

#include "FWCore/ParameterSet/interface/InputTag.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "FWCore/Utilities/interface/EDMException.h"

#include "FWCore/MessageLogger/interface/MessageLogger.h"

#include "FWCore/ServiceRegistry/interface/Service.h"

#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/TrackReco/interface/TrackToTrackMap.h"
#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/MuonReco/interface/MuonFwd.h"
#include "DataFormats/RecoCandidate/interface/IsoDeposit.h"
#include "DataFormats/RecoCandidate/interface/IsoDepositFwd.h"

#include "DataFormats/HepMCCandidate/interface/GenParticle.h"

#include "PhysicsTools/UtilAlgos/interface/TFileService.h"

#include "MuonPtResol/MuonPtResTree/interface/muresolEvent.h"

//
// forward class declaration
//
class TFile;
class TTree;

namespace MuonPtResTree_histo
{
 
  // cone size and SumPt for isolation
  static const unsigned nBinCone = muresol::N_CONESIZE;
  static const float minCone = muresol::MIN_CONE;
  static const float maxCone = muresol::MAX_CONE;
}

//
// class declaration
//
class MuonPtResTree : public edm::EDAnalyzer 
{
  public:
  explicit MuonPtResTree(const edm::ParameterSet&);
  ~MuonPtResTree();
  
  
  private:
  virtual void beginJob() ;
  virtual void analyze(const edm::Event&, const edm::EventSetup&);
  virtual void endJob() ;
  edm::InputTag muonTag_;
  

  struct muonTrack {
    // reference to standard tracking
    reco::MuonRef mu;
    // reference to default-TeV, 1st-Hit, picky and cocktail/optimized tracking
    reco::TrackRef TeVMuons[4];
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
    trigEff(){trigger_count.clear(); Nev = Nev_1mu = 0;}
  };

  typedef std::map<std::string, unsigned>::const_iterator tIt;

  trigEff genMuTrig; // muon-trigger efficiency wrt generator-muons
  trigEff MuTrig; // muon-trigger efficiency for all processed events

  // muon-detector acceptance
 

  // (at least one) generated muon within detector acceptance 
  // (defined as |eta|< detmu_acceptance)
  bool genmu_acceptance;
  // whether event has been accepted by single non-isolated muon (L1 and HL) trigger
  bool muL1_acceptance;
  bool muHLT_acceptance;

  // "golden" single-muon trigger names
  std::string muHLT_20x;
  std::string muHLT_21x;
  std::string muL1;

  // whether this is real-data
  bool realData;

  edm::TriggerNames triggerNames;

  // generator-level muons
  std::vector<reco::GenParticle> gen_muons;

  // software versions used to produce HLT and RECO
  std::string HLTversion;
  std::string RECOversion;
  static const std::string INVALID_RELEASE;
  std::string sample_description;
  // # of produced events before filtering
  int Nprod_evt;

  static bool is31x(const std::string & release_string);
  static bool is22x(const std::string & release_string);
  static bool is21x(const std::string & release_string);
  static bool is20x(const std::string & release_string);

  // TTree structures
  muresol::Event * evt;
  muresol::JobInfo * job;
  std::string software_version;

  //    Histograms, trees and all that
  TTree *tree_job, *tree_event;

  //    Root output file
  edm::Service<TFileService> fs;


  edm::Handle<reco::MuonCollection> muonCollection;
  edm::Handle<reco::IsoDepositMap> tkMapH;
  edm::Handle<reco::IsoDepositMap> ecalMapH;
  edm::Handle<reco::IsoDepositMap> hcalMapH;

  const reco::TrackToTrackMap * tevMap_default;
  const reco::TrackToTrackMap * tevMap_1stHit;
  const reco::TrackToTrackMap * tevMap_picky;
  const reco::TrackToTrackMap * tevMap_oddhits;

  // # of trigger paths in HLT configuration
  unsigned N_triggers;

  // # of reconstructed muons per event
  unsigned N_muons;
  // # of all reconstructed muons per event (included standalone-only)
  unsigned N_all_muons;

  // copy tracking info from reco::Track to muresol::Track
  void getTracking(muresol::Track & track, const reco::Track & p);

  // initialize run info
  void init_run();

  // initialize event info
  void init_event();

  // initialize histograms
  void init_histograms();

  // initialize trigger structure
  void init_trigger(const edm::Handle<edm::TriggerResults> & hltresults);
  
  // print summary info over full job
  void printSummary() const;
  // print summary info for real
  void printSummary2(const trigEff & trig, const std::string & description) const;

  // get the generator info, populate gen_muons, set genmu_acceptance flag
  void getGenParticles(const edm::Event & iEvent);

  // get trigger info, update muTrig/genMuTrig, set muL1/HLT_acceptance flag
  void getTriggers(const edm::Event & iEvent);

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
  void doTeVanalysis(reco::MuonRef mu, muresol::Muon * wpmu);

  // do isolation
  void doIsolation(reco::MuonRef mu,  muresol::Muon * wpmu);

};
