#include "UserCode/CMGWPrimeGroup/interface/Wprime_muonreco.h"
//

// system include files
#include <memory>
#include <cmath>

#include "FWCore/Framework/interface/MakerMacros.h"

#include "DataFormats/Common/interface/TriggerResults.h"

#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/JetReco/interface/CaloJet.h"
#include "DataFormats/MuonReco/interface/MuonCocktails.h"
#include "DataFormats/MuonReco/interface/MuonSelectors.h"

#include "DataFormats/METReco/interface/PFMET.h"
#include "DataFormats/METReco/interface/PFMETCollection.h"
#include "DataFormats/JetReco/interface/PFJetCollection.h"

#include "HLTrigger/HLTcore/interface/HLTConfigProvider.h"

#include "TFile.h"
#include "TTree.h"

#include <iostream>
#include <fstream>
#include <map>

using std::cout; using std::endl; using std::string; using std::ifstream;

using namespace Wprime_muonreco_histo;

/// Constructor
Wprime_muonreco::Wprime_muonreco(const edm::ParameterSet& iConfig):
  pvTag_(iConfig.getParameter<edm::InputTag> ("pvTag")),
  pvBSTag_(iConfig.getParameter<edm::InputTag> ("pvBSTag")),
  muonTag_(iConfig.getParameter<edm::InputTag> ("MuonTag")),
  pfmetTag_(iConfig.getParameter<edm::InputTag> ("pfMetTag")),
  HLTTag_(iConfig.getParameter<edm::InputTag>( "HLTriggerResults" ) ),
  //  isoTag_(iConfig.getParameter<edm::InputTag> ("IsolationTag")),
  caloJetTag_(iConfig.getParameter<edm::InputTag> ("caloJetTag")),
  pfJetTag_(iConfig.getParameter<edm::InputTag> ("pfJetTag")),
  tkIsoMapTag_(iConfig.getParameter<edm::InputTag> ("TkIsoMapTag")),
  ecalIsoMapTag_(iConfig.getParameter<edm::InputTag> ("EcalIsoMapTag")),
  hcalIsoMapTag_(iConfig.getParameter<edm::InputTag> ("HcalIsoMapTag")),
  detmu_acceptance(iConfig.getParameter<double>("Detmu_acceptance")),
  sample_description(iConfig.getParameter<string>("description"))
{
   //now do what ever initialization is needed
  tree_job = tree_run = tree_event = 0;

  evt = new wprime::Event(); job = new wprime::JobInfo();
  run = new wprime::RunInfo();
  job->sample = sample_description;
  software_version = "V00-00-00";
  firstEventInRun = false;
}


/// Destructor
Wprime_muonreco::~Wprime_muonreco()
{
  if(evt) delete evt; if(job) delete job; if(run) delete run;
}

// get the generator info, populate gen_muons, set genmu_acceptance flag
void Wprime_muonreco::getGenParticles(const edm::Event & iEvent)
{
  if(realData)return;

  TClonesArray & mcmu = *(evt->mu_mc);
  TClonesArray & mcneu = *(evt->neu_mc);
  TClonesArray & mcw = *(evt->w_mc);
  TClonesArray & mcwp = *(evt->wp_mc);

  int Nm = 0; int Nn = 0; int Nw = 0; int Nwp = 0;

  edm::Handle<reco::GenParticleCollection> genParticles;
  iEvent.getByLabel("genParticles", genParticles);
  for(size_t i = 0; i != genParticles->size(); ++i) {
    const GenParticle & p = (*genParticles)[i];
    int id = p.pdgId(); 
    if( abs(id)!=13 && abs(id) != 14 && abs(id) != 24 && abs(id) != 34)
      continue; // keep only muons, neutrino, W or W'

    int st = p.status(); int ch = p.charge();
    int momid = -9999;

    TLorentzVector p4(p.px(), p.py(), p.pz(),  p.energy());
    const Candidate * mom = p.mother();
    if(mom)
      momid = mom->pdgId();

    switch(abs(id))
      {
      case 13: // muon
	new(mcmu[Nm]) wprime::MCParticle(p4, ch, momid, st);
	++Nm;
	
	if(st == 1)
	  {
	    gen_muons.push_back(p);
	    if(abs(p.eta())  < detmu_acceptance)
	      genmu_acceptance = true;
	    
	  }

	break;
	
      case 14: // muon-neutrino
	new(mcneu[Nn]) wprime::MCParticle(p4, ch, momid, st);
	++Nn;

	break;
	
      case 24: // W
	new(mcw[Nw]) wprime::MCParticle(p4, ch, momid, st);
	++Nw;
      
	break;
	
      case 34: // W'
	new(mcwp[Nwp]) wprime::MCParticle(p4, ch, momid, st);
	++Nwp;
	break;

      default:
	{
	  ; // do nothing for now
	}

      } // switch (abs(id))


  } // loop over genParticles

}

// get trigger info, update MuTrig/genMuTrig
void Wprime_muonreco::getTriggers(const edm::Event & iEvent)
{
  edm::Handle<edm::TriggerResults> hltresults;
  iEvent.getByLabel(HLTTag_,hltresults);

  It trig = triggerNames.begin();
  for (unsigned itrig = 0; itrig != N_triggers; ++itrig)
    { // loop over triggers
      string trigName = *trig; ++trig;
      
      // skip loop for non-muon triggers....
      // the trigger class name (i.e. Muon) should be release-dependent ???
      if(trigName.find("Mu") == string::npos)continue; // in 21x and later
      //      if(trigName.find("HLT1Muon") == string::npos)continue; // in 20x
      
      bool accept = hltresults->accept(itrig);

      // this is where we keep statistics (to be printed out at end of job)
      if(accept)
	{ // extract event counts for trigger efficiencies
	  tIt it;
	  if(genmu_acceptance)
	    {
	      it = genMuTrig.trigger_count.find(trigName);
	      if(it != genMuTrig.trigger_count.end())
		(genMuTrig.trigger_count[trigName])++; 
	      else
		genMuTrig.trigger_count[trigName] = 1;
	    }
	  it = MuTrig.trigger_count.find(trigName);
	  if(it != MuTrig.trigger_count.end())
	    (MuTrig.trigger_count[trigName])++;
	  else
	    MuTrig.trigger_count[trigName] = 1;
	} // extract event counts for trigger efficiencies

      if(trigName == "HLT_L1MuOpen")
	evt->HLT_L1MuOpen = accept? 1 : 0;
      else if (trigName == "HLT_L1Mu")
	evt->HLT_L1Mu = accept? 1 : 0;
      else if (trigName == "HLT_Mu3")
	evt->HLT_Mu3 = accept? 1 : 0;
      else if (trigName == "HLT_Mu5")
	evt->HLT_Mu5 = accept? 1 : 0;
      else if (trigName == "HLT_Mu7")
	evt->HLT_Mu7 = accept? 1 : 0;
      else if (trigName == "HLT_Mu9")
	evt->HLT_Mu9 = accept? 1 : 0;
      else if (trigName == "HLT_Mu11")
	evt->HLT_Mu11 = accept? 1 : 0;
      else if (trigName == "HLT_L2Mu5")
	evt->HLT_L2Mu5 = accept? 1 : 0;
      else if (trigName == "HLT_L2Mu9")
	    evt->HLT_L2Mu9 = accept? 1 : 0;
      else if (trigName == "HLT_L2Mu11")
	evt->HLT_L2Mu11 = accept? 1 : 0;
      else if (trigName == "HLT_L2Mu15")
	evt->HLT_L2Mu15 = accept? 1 : 0;
      else if (trigName == "HLT_L2Mu25")
	evt->HLT_L2Mu25 = accept? 1 : 0;
    } // loop over triggers

}

// get Particle-Flow MET
void Wprime_muonreco::getPFMET(const edm::Event & iEvent)
{
  // Get the MET collection from the event
  edm::Handle<reco::PFMETCollection> pfmetCollection;
  iEvent.getByLabel(pfmetTag_, pfmetCollection);
  
  PFMETCollection::const_iterator pfMET = pfmetCollection->begin();
  met_x += pfMET->px();
  met_y += pfMET->py();
  evt->pfmet.Set(met_x, met_y);
}

// get Jets
void Wprime_muonreco::getJets(const edm::Event & iEvent)
{
  // Get the Jet collections from the event
  edm::Handle<reco::CaloJetCollection> jetCaloCollection;
  iEvent.getByLabel(caloJetTag_, jetCaloCollection);
  edm::Handle<reco::PFJetCollection> jetPFCollection;
  iEvent.getByLabel(pfJetTag_, jetPFCollection);

  TClonesArray & cjet = *(evt->calojet);
  TClonesArray & pfjet = *(evt->pfjet);

  //get calo jets
  int j = 0;

  for (CaloJetCollection::const_iterator jet_ = jetCaloCollection->begin(); 
       jet_ !=jetCaloCollection->end(); ++jet_) { // loop over jets
      
      new(cjet[j]) TLorentzVector(jet_->px(), jet_->py(), jet_->pz(), 
                                 jet_->energy());
      ++j;
  } // loop over calojets
  

  //get pf jets
  j = 0;
  for (PFJetCollection::const_iterator jet_ = jetPFCollection->begin(); 
       jet_ !=jetPFCollection->end(); ++jet_) { // loop over jets
      
      new(pfjet[j]) TLorentzVector(jet_->px(), jet_->py(), jet_->pz(), 
                                 jet_->energy());
      ++j;
  } // loop over pfjets

}  

// get Isolation
void Wprime_muonreco::getIsolation(const edm::Event & iEvent)
{
  // take iso deposits for tracks 
  // (contains (eta,phi, pt) of tracks in R<X (1.0) around each muon)
  iEvent.getByLabel(tkIsoMapTag_, tkMapH);

  // take iso deposits for ECAL 
  // (contains (eta,phi, pt) of ecal in R<X (1.0) around each muon)
  iEvent.getByLabel(ecalIsoMapTag_, ecalMapH);

  // take iso deposits for HCAL 
  // (contains (eta,phi, pt) of hcal in R<X (1.0) around each muon)
  iEvent.getByLabel(hcalIsoMapTag_, hcalMapH);
}

// get TeV muons
void Wprime_muonreco::getTeVMuons(const edm::Event & iEvent)
{
  edm::Handle<reco::TrackToTrackMap> tevMapH_default;
  edm::Handle<reco::TrackToTrackMap> tevMapH_1stHit;
  edm::Handle<reco::TrackToTrackMap> tevMapH_picky;
  edm::Handle<reco::TrackToTrackMap> tevMapH_dyt;

  iEvent.getByLabel("tevMuons", "default", tevMapH_default);
  tevMap_default = tevMapH_default.product();

  iEvent.getByLabel("tevMuons", "firstHit", tevMapH_1stHit);
  tevMap_1stHit = tevMapH_1stHit.product();

  iEvent.getByLabel("tevMuons", "picky", tevMapH_picky);
  tevMap_picky = tevMapH_picky.product();

  iEvent.getByLabel("tevMuons", "dyt", tevMapH_dyt);
  tevMap_dyt = tevMapH_dyt.product();
}

// get primary vertex info
void Wprime_muonreco::getPVs(const edm::Event & iEvent)
{
  iEvent.getByLabel(pvTag_, pvCollection);
  iEvent.getByLabel(pvBSTag_, pvBSCollection);
  evt->Npv = pvCollection->size();
  evt->NpvBS = pvBSCollection->size();
}

// get muons
void Wprime_muonreco::getMuons(const edm::Event & iEvent)
{
  // Get the Muon Track collection from the event
  iEvent.getByLabel(muonTag_, muonCollection);

  if(RECOversion == INVALID_RELEASE)
    { // need a better way of doing this - in beginRun???
      const edm::Provenance& prov=iEvent.getProvenance(muonCollection.id());
      RECOversion = prov.releaseVersion();
      cout << " CMSSW release used to produce muon RECO : "
	   << RECOversion << endl;
      
      job->RECOversion = RECOversion;
    }

  // # of reconstructed muons in event (including standalone-only)
  N_all_muons = muonCollection->size();
}

// copy tracking info from reco::Track to wprime::Track
void Wprime_muonreco::getTracking(wprime::Track & track, const reco::Track & p)
{
  TVector3 p3(p.px(), p.py(), p.pz());
  track.p.SetVectM(p3, wprime::MUON_MASS);
  track.q = p.charge();
  track.chi2 = p.chi2();
  track.d0 = p.d0();
  track.dd0 = p.d0Error();
  track.dpt = p.ptError();
  track.dq_over_p = p.qoverpError();
  track.ndof = int(p.ndof());
  track.Nstrip_layer = p.hitPattern().stripLayersWithMeasurement();
  track.Npixel_layer = p.hitPattern().pixelLayersWithMeasurement();
  track.Nstrip_layerNoMeas = p.hitPattern().stripLayersWithoutMeasurement();
  track.Npixel_layerNoMeas = p.hitPattern().pixelLayersWithoutMeasurement();
  track.NsiStrip_hits = p.hitPattern().numberOfValidStripHits();
  track.Npixel_hits = p.hitPattern().numberOfValidPixelHits();
  track.Nmuon_hits = p.hitPattern().numberOfValidMuonHits();
  track.Ntot_hits = p.numberOfValidHits();
  track.Ntrk_hits = p.hitPattern().numberOfValidTrackerHits();
}

// fill in with dummy values when there is no track
void Wprime_muonreco::getNullTracking(wprime::Track & track)
{
  float wrong = -99999;
  TVector3 p3(wrong, wrong, wrong);
  track.p.SetVectM(p3, wrong);
  track.q = wrong;
  track.chi2 = track.d0 = track.dd0 = track.dpt = track.dq_over_p = wrong;
  track.ndof = track.Nstrip_layer = track.Npixel_layer = 
    track.Nstrip_layerNoMeas = track.Npixel_layerNoMeas = 
    track.NsiStrip_hits = track.Npixel_hits = track.Nmuon_hits = 
    track.Ntot_hits = track.Ntrk_hits = -1;
}

// do muon analysis
void Wprime_muonreco::doMuons()
{
  TClonesArray & recomu = *(evt->mu);

  for(unsigned i = 0; i != N_all_muons; ++i) 
    { // loop over reco muons 
      
      reco::MuonRef mu(muonCollection, i);
      if(!(mu->isGlobalMuon()) )
	continue; // keep only global muons
      
      new(recomu[N_muons]) wprime::Muon(); 
      wprime::Muon * wpmu = (wprime::Muon *) recomu[N_muons];
      wpmu->Nmu_hits = mu->standAloneMuon()->recHitsSize(); 
      wpmu->Nmatches = mu->numberOfMatches();

      getTracking(wpmu->tracker, *(mu->track()));
      getTracking(wpmu->global, *(mu->combinedMuon()) );

      wpmu->GlobalMuonPromptTight = muon::isGoodMuon(*mu, muon::GlobalMuonPromptTight);
      wpmu->TMLastStationLoose = muon::isGoodMuon(*mu, muon::TMLastStationLoose);
      wpmu->TMLastStationTight = muon::isGoodMuon(*mu, muon::TMLastStationTight);
      wpmu->TMLastStationAngTight = muon::isGoodMuon(*mu, muon::TMLastStationAngTight);
      wpmu->AllGlobalMuons = muon::isGoodMuon(*mu, muon::AllGlobalMuons);
      wpmu->AllStandAloneMuons = muon::isGoodMuon(*mu, muon::AllStandAloneMuons);
      wpmu->AllTrackerMuons = muon::isGoodMuon(*mu, muon::AllTrackerMuons);


      ++N_muons;
      
      doIsolation(mu, wpmu);

      doTeVanalysis(mu, wpmu);
      
    } // loop over reco muons

  if(1 == N_muons)
    {
      ++(MuTrig.Nev_1mu);
      if(genmu_acceptance)
	++(genMuTrig.Nev_1mu);
    
    }
  
}

// get TMR track, from 
// http://cmssw.cvs.cern.ch/cgi-bin/cmssw.cgi/CMSSW/SUSYBSMAnalysis/Zprime2muAnalysis/src/MuonCocktails.cc?revision=1.1&view=markup
const reco::TrackRef Wprime_muonreco::getTMR(const reco::TrackRef& trackerTrack,
					     const reco::TrackRef& fmsTrack,
					     const double cut) 
{
  double probTK  = 0;
  double probFMS = 0;
  
  if (trackerTrack.isNonnull() && trackerTrack->numberOfValidHits())
    probTK = muon::trackProbability(trackerTrack);
  if (fmsTrack.isNonnull() && fmsTrack->numberOfValidHits())
    probFMS = muon::trackProbability(fmsTrack);
  
  bool TKok  = probTK > 0;
  bool FMSok = probFMS > 0;

  if (TKok && FMSok) {
    if (probFMS - probTK > cut)
      return trackerTrack;
    else
      return fmsTrack;
  }
  else if (FMSok)
    return fmsTrack;
  else if (TKok)
    return trackerTrack;
  else
    return reco::TrackRef();
}


// do TeV-muon analysis
void Wprime_muonreco::doTeVanalysis(reco::MuonRef mu, wprime::Muon * wpmu)
{
  if(!(mu->isGlobalMuon()) ) return; // keep only global muons

  TrackToTrackMap::const_iterator iTeV_default;
  TrackToTrackMap::const_iterator iTeV_1stHit;
  TrackToTrackMap::const_iterator iTeV_picky;
  TrackToTrackMap::const_iterator iTeV_dyt;

  iTeV_default = tevMap_default->find(mu->globalTrack());
  iTeV_1stHit = tevMap_1stHit->find(mu->globalTrack());
  iTeV_picky = tevMap_picky->find(mu->globalTrack());
  iTeV_dyt = tevMap_dyt->find(mu->globalTrack());

  bool TeVfailed = false;

  if(iTeV_1stHit == tevMap_1stHit->end())
    {
      getNullTracking(wpmu->tpfms);
      TeVfailed = true;
    }
  else
    getTracking(wpmu->tpfms, *(iTeV_1stHit->val) );

  if(iTeV_picky == tevMap_picky->end())
    {
      getNullTracking(wpmu->picky);
      TeVfailed = true;
    }
  else
    getTracking(wpmu->picky, *(iTeV_picky->val) );

  if(iTeV_dyt == tevMap_dyt->end())
    {
      getNullTracking(wpmu->dyt);
      TeVfailed = true;
    }
  else
    getTracking(wpmu->dyt, *(iTeV_dyt->val) );


  if(TeVfailed)
    {
      getNullTracking(wpmu->cocktail);
    }
  else
    {
      reco::TrackRef cocktail = 
	muon::tevOptimized(mu->combinedMuon(), mu->track(), 
			   *tevMap_default, *tevMap_1stHit, 
			   *tevMap_picky);
      getTracking(wpmu->cocktail, *cocktail);
    }

  if(iTeV_1stHit == tevMap_1stHit->end())
    getNullTracking(wpmu->tmr);
  else
    {
      reco::TrackRef tmr = getTMR(mu->track(), iTeV_1stHit->val);
      getTracking(wpmu->tmr, *tmr);
    }

}



// do isolation
void Wprime_muonreco::doIsolation(reco::MuonRef mu,  wprime::Muon * wpmu)
{
  // General Isolation
  const reco::IsoDeposit tkDep((*tkMapH)[mu]);
  const reco::IsoDeposit ecalDep((*ecalMapH)[mu]);
  const reco::IsoDeposit hcalDep((*hcalMapH)[mu]);

  for(unsigned i = 0; i != nBinCone; ++i) 
    { // loop over cone sizes
      float coneSize = minCone+i*(maxCone-minCone)/nBinCone;
      double ptsum_cone=tkDep.depositWithin(coneSize); 
      wpmu->SumPtIso[i] = ptsum_cone;

      double ntk_cone=tkDep.depositAndCountWithin(coneSize).second;
      wpmu->NtrkIso[i] = ntk_cone;

      double ecetsum_cone=ecalDep.depositWithin(coneSize); 
      wpmu->ECALIso[i] = ecetsum_cone;

      double hcetsum_cone=hcalDep.depositWithin(coneSize); 
      wpmu->HCALIso[i] = hcetsum_cone;
    } // loop over cone sizes

}

// ------------ method called to for each event  ------------
void
Wprime_muonreco::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  if(firstEventInRun)
    {
      init_run(iEvent);
      firstEventInRun = false;
    }

  init_event();
  evt->evt_no = iEvent.id().event();
  evt->run_no = iEvent.id().run();
  evt->LS_no  = iEvent.id().luminosityBlock();

  getGenParticles(iEvent);

  ++(MuTrig.Nev); // # of processed events
  if(genmu_acceptance)
    // # of processed events with gen-muon within det.acceptance
    ++(genMuTrig.Nev);

  ++(run->Nproc_evt);

  getTriggers(iEvent);

  getPVs(iEvent);

  getPFMET(iEvent);

  getJets(iEvent);

  getIsolation(iEvent);

  getMuons(iEvent);
  getTeVMuons(iEvent);
  doMuons();

  tree_event->Fill();
}

// initialize histograms
void Wprime_muonreco::init_histograms()
{
  string common_desc = " - UserCode/CMGWPrimeGroup version " + software_version;
  string tree_title = "Job info" + common_desc;
  tree_job = fs->make<TTree>("jobinfo", tree_title.c_str());
  tree_title = "Run info" + common_desc;
  tree_run = fs->make<TTree>("runinfo", tree_title.c_str());
  tree_title = "Wprime kinematic info per event" + common_desc;
  tree_event = fs->make<TTree>("wprime", tree_title.c_str());

  tree_job->Branch("job", "wprime::JobInfo", &job, 8000, 2);
  tree_run->Branch("run", "wprime::RunInfo", &run, 8000, 2);
  tree_event->Branch("wp", "wprime::Event", &evt, 8000, 2);
}

const string Wprime_muonreco::INVALID_RELEASE = "invalid release number";

// initialize run info
void Wprime_muonreco::init_run(const edm::Event& iEvent)
{
  realData = iEvent.isRealData();
  check_trigger(iEvent);
}

// initialize event info
void Wprime_muonreco::init_event()
{
  gen_muons.clear();
  N_muons = N_all_muons = 0;
  genmu_acceptance = false;
  met_x = met_y = met = 0.0;

  evt->mu_mc->Clear(); evt->neu_mc->Clear(); 
  evt->w_mc->Clear(); evt->wp_mc->Clear(); 
  evt->mu->Clear(); 
  evt->calojet->Clear();evt->pfjet->Clear();

  evt->reset_triggers();
 
}


// ------------ method called once each job just before starting event loop  ----
void Wprime_muonreco::beginJob()
{
  const char * _path = getenv("CMSSW_BASE");
  const char * _filename = "/src/UserCode/CMGWPrimeGroup/VERSION";

  if(_path == 0)
    {
      cerr << " *** Wprime_muonreco::beginJob *** " << endl;
      cerr << " *** Error! can't locate working directory ! ***\n";
      cerr << " (maybe the UserCode/CMGWPrimeGroup directory is missing?)\n"
	   << endl;
    }
  else
    {
      char verFile[1024];
      sprintf(verFile, "%s%s", _path, _filename);
      cout << " Found version file: " << verFile << endl;

      ifstream cin; cin.open(verFile);
      while(cin >> software_version)
	{
	  ; // write version # into string
	}
      cout << " Extracted version: " << software_version << endl;
    }

  HLTversion = RECOversion = INVALID_RELEASE;

  genMuTrig.clear(); MuTrig.clear();
  init_histograms();
 
}

void Wprime_muonreco::endRun(edm::Run const &, edm::EventSetup const &)
{
  cout << " Wprime_muonreco::endRun: completed run " << run->run_no 
       << ", processed " << run->Nproc_evt << " events " << endl;
  tree_run->Fill();
}

void Wprime_muonreco::beginRun(edm::Run const & iRun, 
			       edm::EventSetup const & iSetup)
{
  firstEventInRun = true;

  HLTConfigProvider hltConfig;
  bool changed(true);
  if(! hltConfig.init(iRun, iSetup, HLTTag_.process(), changed))
    {
      cerr << "  -- Failed to init HLT Config" << endl; 
      abort(); // what is the proper way of handling this??? skip run?
    }
 
  triggerNames = hltConfig.triggerNames();
  N_triggers = hltConfig.size();
  run->HLTmenu = hltConfig.tableName();

  run->run_no = iRun.run(); 
  run->Nproc_evt = 0;

  cout << "\n Run # " << run->run_no << " recorderd with HLT menu " 
       << run->HLTmenu << endl;
}

// ------------ method called once each job just after ending the event loop  ------------
void Wprime_muonreco::endJob() {

  printSummary();
  tree_job->Fill();
}

// print summary info over full job
void Wprime_muonreco::printSummary() const
{
  cout << " Processed " << MuTrig.Nev << " events" << endl;
  cout << " " << genMuTrig.Nev << " of those contained at least one MC muon"
       << " within |eta| < " << detmu_acceptance  << endl;
  
  std::ostringstream description; 
  string disclaimer = " WITHOUT confirming that all trigger paths have been enabled for all runs";
  description << " Trigger efficiencies for events\n w/ gen-muon within det. acceptance (" << genMuTrig.Nev << " events processed)\n";
  description<< disclaimer << endl;
  if(!realData)
    printSummary2(genMuTrig, description.str());

  description.str(""); // interesting way to clear the description content...
  description << " Trigger efficiencies for all events (" << MuTrig.Nev 
		<< " events processed)\n";
  description<< disclaimer << endl;
  printSummary2(MuTrig, description.str());
}

// print summary info for real
void Wprime_muonreco::printSummary2(const trigEff & trig, 
				    const string & description) const
{
  string longline("=============================================================");
  cout << longline << endl;
  cout<< description << endl << endl;
  
  for(tIt it = trig.trigger_count.begin(); it != trig.trigger_count.end();
      ++it)
    {
      float eff = 1.*it->second/(trig.Nev);
      float d_eff = sqrt(eff*(1-eff)/trig.Nev);
      cout << " TriggerName " << it->first;
      cout << ", Accepted events: " << it->second
	   << " (" << 100*eff << " +- " << 100*d_eff << " %) " << endl;
    }
  cout << endl;
}

// check trigger is there (at beginRun)
void Wprime_muonreco::check_trigger(const edm::Event & iEvent)
{
  edm::Handle<edm::TriggerResults> hltresults;
  iEvent.getByLabel(HLTTag_,hltresults);
  
  if (! hltresults.isValid() ) 
    { 
      cerr << "  -- No HLTRESULTS" << endl; 
      abort(); // what is the proper way of handling this??? skip run?
    }
    
  const edm::Provenance & prov = iEvent.getProvenance(hltresults.id());
  HLTversion = prov.releaseVersion();
  cout << " CMSSW release used to produce trigger decisions: "
       << HLTversion << endl;
  
  run->HLTversion = HLTversion;
}


//define this as a plug-in
DEFINE_FWK_MODULE(Wprime_muonreco);
