// -*- C++ -*-
//
// Package:    ZTemplateForSampling
// Class:      ZTemplateForSampling
// 
//
// Original Author:  Antonio Villela 
// Modified (after 22X) by Maria Cepeda 
// Adpated to W' Analysis by Edgar Carrera


// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Utilities/interface/InputTag.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "DataFormats/FWLite/interface/Handle.h"
#include "DataFormats/FWLite/interface/Event.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"


#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/Candidate/interface/CompositeCandidate.h"

#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/MuonReco/interface/MuonFwd.h"
#include "DataFormats/MuonReco/interface/MuonCocktails.h"
#include "DataFormats/MuonReco/interface/MuonSelectors.h"

#include "DataFormats/METReco/interface/CaloMET.h"
#include "DataFormats/METReco/interface/PFMET.h"
#include "DataFormats/METReco/interface/PFMETCollection.h"

#include "DataFormats/Candidate/interface/Particle.h"

#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "DataFormats/BeamSpot/interface/BeamSpot.h"

#include "DataFormats/Common/interface/View.h"
#include "DataFormats/Common/interface/AssociationVector.h"
#include "DataFormats/Common/interface/AssociationMap.h"
#include "DataFormats/Common/interface/RefToBase.h"

#include "FWCore/Common/interface/TriggerNames.h"
#include "DataFormats/Common/interface/TriggerResults.h"

#include "TH1D.h"
#include "TH2D.h"
#include <iostream>
#include <fstream>
#include <vector>
#include <algorithm>


//
// class declaration
//

class ZTemplateForSampling : public edm::EDAnalyzer {
   public:
      explicit ZTemplateForSampling(const edm::ParameterSet&);
      ~ZTemplateForSampling();

  bool GoodEWKMuon(reco::MuonRef Mu);
    //reco::Track GetCockTailMu(const reco::MuonRef& mu);

   private:
      virtual void beginJob() ;
      virtual void analyze(const edm::Event&, const edm::EventSetup&);
      virtual void endJob() ;



   TH1D *hZmass, *hVBPt, *hVBPhi, *hVBEta, *hVBPt2;
   TH1D *hMET, *hMETParal, *hMETPerp, *hMETPhi, *hMETPhiMinusVBPhi;
   TH2D *hMETvsVBPt, *hMETPhivsVBPhi,*hMETParalvsVBPt, *hMETPerpvsVBPt; 

   TH1D *hptMu, *hProbe_ptMu, *hIsoAll_ptMu, *hIso, *hProbe_Iso;

   Long_t nEvts;
  
  edm::InputTag trigTag_;
  edm::InputTag muonTag_;
  edm::InputTag metTag_;
  edm::InputTag jetTag_;
    std::vector<std::string> expressions;

    //   const reco::TrackToTrackMap * tevMap_default;
    //const reco::TrackToTrackMap * tevMap_1stHit;
    //const reco::TrackToTrackMap * tevMap_picky;


    double wMass_;
  double zMass_;

  double minZmass_;
  double maxZmass_;
  
  double ptCut_;
  double etaCut_;
  bool isRelativeIso_;
  bool isCombinedIso_;
  double isoCut03_;
  double mtMin_;
  double mtMax_;
  double metMin_;
  double metMax_;
  double acopCut_;

  double dxyCut_;
  double normalizedChi2Cut_;
  int trackerHitsCut_;
  int pixelHitsCut_;
  int muonHitsCut_;
  bool isAlsoTrackerMuon_;
  int nMatchesCut_;




  int selectByCharge_;

   
};

ZTemplateForSampling::ZTemplateForSampling(const edm::ParameterSet& iConfig):
      // Input collections
      trigTag_(iConfig.getUntrackedParameter<edm::InputTag> ("TrigTag", edm::InputTag("TriggerResults::HLT"))),
      muonTag_(iConfig.getUntrackedParameter<edm::InputTag> ("MuonTag", edm::InputTag("muons"))),
      metTag_(iConfig.getUntrackedParameter<edm::InputTag> ("METTag")),
      expressions(iConfig.getParameter<std::vector<std::string> >("MuonTrig")),
      //Especific Parameters for Template: 
      wMass_(iConfig.getUntrackedParameter<double> ("WMass", 80.403)),
      zMass_(iConfig.getUntrackedParameter<double> ("ZMass", 91.1876)),

      //Z Selection
      minZmass_(iConfig.getUntrackedParameter<double> ("MinZMass",70.)),
      maxZmass_(iConfig.getUntrackedParameter<double> ("MaxZMass",110.)),

      // Main cuts for W 
      ptCut_(iConfig.getUntrackedParameter<double>("PtCut", 20.)),
      etaCut_(iConfig.getUntrackedParameter<double>("EtaCut", 2.1)),
      isRelativeIso_(iConfig.getUntrackedParameter<bool>("IsRelativeIso", true)),
      isCombinedIso_(iConfig.getUntrackedParameter<bool>("IsCombinedIso", true)),
      isoCut03_(iConfig.getUntrackedParameter<double>("IsoCut03", 0.15)),
      mtMin_(iConfig.getUntrackedParameter<double>("MtMin", -999999.)),
      mtMax_(iConfig.getUntrackedParameter<double>("MtMax", 999999.)),
      metMin_(iConfig.getUntrackedParameter<double>("MetMin", -999999.)),
      metMax_(iConfig.getUntrackedParameter<double>("MetMax", 999999.)),
      acopCut_(iConfig.getUntrackedParameter<double>("AcopCut", 999.)),

      // Muon quality cuts
      dxyCut_(iConfig.getUntrackedParameter<double>("DxyCut", 0.2)),   // dxy < 0.2 cm 
      normalizedChi2Cut_(iConfig.getUntrackedParameter<double>("NormalizedChi2Cut", 10.)), // chi2/ndof (of global fit) <10.0
      trackerHitsCut_(iConfig.getUntrackedParameter<int>("TrackerHitsCut", 11)),  // Tracker Hits >10 
      pixelHitsCut_(iConfig.getUntrackedParameter<int>("PixelHitsCut", 1)), // Pixel Hits >0
      muonHitsCut_(iConfig.getUntrackedParameter<int>("MuonHitsCut", 1)),  // Valid Muon Hits >0 
      isAlsoTrackerMuon_(iConfig.getUntrackedParameter<bool>("IsAlsoTrackerMuon", true)),
      nMatchesCut_(iConfig.getUntrackedParameter<int>("NMatchesCut", 2)), // At least 2 Chambers with matches 

      // W+/W- Selection
      selectByCharge_(iConfig.getUntrackedParameter<int>("SelectByCharge", 0))
{
   LogDebug("ZTemplateForSampling")<<" ZTemplateForSamplingAnalyzer constructor called";
   
}


ZTemplateForSampling::~ZTemplateForSampling()
{

}

bool ZTemplateForSampling::GoodEWKMuon(reco::MuonRef Mu){
     bool goodMuon=true;
     if (!Mu->isGlobalMuon()||!Mu->isTrackerMuon()){goodMuon=false;}
     else{
     reco::BeamSpot vertexBeamSpot;     edm::Handle<reco::BeamSpot> recoBeamSpotHandle;

            reco::TrackRef gm = Mu->globalTrack();

            // Quality cuts
            double dxy = gm->dxy(vertexBeamSpot.position());
            double normalizedChi2 = gm->normalizedChi2();
            int trackerHits = gm->hitPattern().numberOfValidTrackerHits();
            int pixelHits = gm->hitPattern().numberOfValidPixelHits();
            int muonHits = gm->hitPattern().numberOfValidMuonHits();
            int nMatches = Mu->numberOfMatches();

            goodMuon = (fabs(dxy)<=dxyCut_)&&(normalizedChi2<=normalizedChi2Cut_)&&(trackerHits>=trackerHitsCut_)&&(pixelHits>=pixelHitsCut_)&&(muonHits>=muonHitsCut_)&&(nMatches>=nMatchesCut_);

     }

     return goodMuon;
}



// ------------ method called to for each event  ------------
void
ZTemplateForSampling::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
   using namespace edm;
   TH1::SetDefaultSumw2(true);
   bool trigger_fired = false;

   // Trigger
   Handle<TriggerResults> triggerResults;
   if (!iEvent.getByLabel(trigTag_, triggerResults)) {
       LogError("") << ">>> TRIGGER collection does not exist !!!";
       return;
   }
   const edm::TriggerNames & triggerNames = iEvent.triggerNames(*triggerResults);
   const int sizetrignames = triggerNames.size();
   //loop over triggers of interest
   for (unsigned int i = 0; i < expressions.size(); ++i){
       string mytrig = expressions[i];
       int itrig1 = triggerNames.triggerIndex(mytrig);
       if (itrig1 == sizetrignames) continue;
       if (triggerResults->accept(itrig1)) {
           trigger_fired = true;
           LogTrace("") << ">>> Trigger bit: " 
                        << trigger_fired << " (" << mytrig << ")";
           break;
       }
   
   }//loop over triggers of interest
   if (!trigger_fired) return;

   //Prepare Muon collections and maps
   Handle<reco::MuonCollection> muonCollection;
   iEvent.getByLabel("muons", muonCollection);
   //Handle<reco::TrackToTrackMap> tevMapH_default;
   //Handle<reco::TrackToTrackMap> tevMapH_1stHit;
   //Handle<reco::TrackToTrackMap> tevMapH_picky;
   //iEvent.getByLabel("tevMuons", "default", tevMapH_default);
   //tevMap_default = tevMapH_default.product();
   //iEvent.getByLabel("tevMuons", "firstHit", tevMapH_1stHit);
   //tevMap_1stHit = tevMapH_1stHit.product();
   //iEvent.getByLabel("tevMuons", "picky", tevMapH_picky);
   //tevMap_picky = tevMapH_picky.product();
  
  
   nEvts++;
   LogDebug("ZTemplateForSampling")<< ">>> Event number: " << nEvts;


   // GET MET:
      // MET
      Handle<View<MET> > metCollection;
      if (!iEvent.getByLabel(metTag_, metCollection)) {
            LogError("") << ">>> MET collection does not exist !!!";
            return;
      }
      //const MET& Met = metCollection->at(0);
      edm::Ptr<reco::MET> met(metCollection,0);

     double MET_Ex=met->px();
     double MET_Ey=met->py();
     double MET_Et=met->pt();


     
     // Do Z selection
     bool eventSelected = true, eventSelected_foriso=true;


     //reco::Candidate::const_iterator zmmCand, zbegin = zMuMu.begin(), zend = zMuMu.end();
     //for(zmmCand=zbe; zmmCand!=zend; ++zmmCand){

     double pt1=0, pt2=0; unsigned int i1=0,i2=1, nGLBmuons=0;

     // Get the two highest pt muons :-).
        for (unsigned int i=0; i<muonCollection->size(); i++){
                MuonRef mu(muonCollection,i);
        LogTrace("") << "mu pt,eta,iso= " << mu->pt() << " , " << mu->eta()<<"  global? "<<mu->isGlobalMuon();
                bool quality=GoodEWKMuon(mu);
                if (!quality)continue;
                nGLBmuons++;
                //reco::Track cktmu = GetCockTailMu(mu);
                if (mu->pt()>pt1){pt1=mu->pt(); i1=i;}
                //if (cktmu.pt()>pt1){pt1=cktmu.pt(); i1=i;}
        }

        for (unsigned int j=0; j<muonCollection->size(); j++){
                if (j==i1)continue;
                MuonRef mu(muonCollection,j);
                bool quality=GoodEWKMuon(mu);
                if (!quality)continue;
                //reco::Track cktmu = GetCockTailMu(mu);
                if (mu->pt()>pt2){pt2=mu->pt(); i2=j;}
                //if (cktmu.pt()>pt2){pt2=cktmu.pt(); i2=j;}
        }
        LogTrace("")<<"global muons: "<<nGLBmuons;
        if (nGLBmuons<2) return;

        // Build Z
        MuonRef muon; MuonRef neutrino;
        reco::MuonRef muon1(muonCollection,i1);
        reco::MuonRef muon2(muonCollection,i2);
        //reco::Track cktmuon1 = GetCockTailMu(muon1);
        //reco::Track cktmuon2 = GetCockTailMu(muon2);
        

        double SumPt1 = muon1->isolationR03().sumPt; double isovar1=SumPt1;
        double Cal1   = muon1->isolationR03().emEt + muon1->isolationR03().hadEt; if(isCombinedIso_)isovar1+=Cal1;
        if (isRelativeIso_) isovar1 /= muon1->pt();
        bool iso1 = (isovar1<=isoCut03_);

        double SumPt2 = muon2->isolationR03().sumPt; double isovar2=SumPt2;
        double Cal2   = muon2->isolationR03().emEt + muon2->isolationR03().hadEt; if(isCombinedIso_)isovar2+=Cal2;
        if (isRelativeIso_) isovar2 /= muon2->pt();
        bool iso2 = (isovar2<=isoCut03_);

        LogTrace("") << "muon 1 pt,eta,iso= " << muon1->pt() << " , " << muon1->eta() << " , " << iso1;
        LogTrace("") << "muon 2 pt,eta,iso= " << muon2->pt() << " , " << muon2->eta() << " , " << iso2;

        // Cuts for Z?
       const math::XYZTLorentzVector myZ (muon1->px() + muon2->px(), muon1->py() + muon2->py(), muon1->pz() + muon2->pz(), muon1->p() + muon2->p());

       if ( muon1->pt()<ptCut_ || fabs(muon1->eta())>etaCut_ || muon2->pt()<ptCut_ || fabs(muon2->eta())>etaCut_ ) {
		eventSelected_foriso=false; 
 		if (!iso1 || !iso2){ eventSelected=false;}
       } else {LogTrace("")<<"Muons for Z passed the kinematic cuts..";}

      if((myZ.mass() < minZmass_)||(myZ.mass() > maxZmass_)) {eventSelected=false; eventSelected_foriso=false;}
      else{LogTrace("")<<myZ.mass()<<"--->Z passed the mass cuts :-)";}



     LogTrace("")<< "after z loop";
     if(eventSelected){
     // Sum back the muons (in case muon corrections used)
     MET_Ex += muon1->px();
     MET_Ex += muon2->px();
     MET_Ey += muon1->py();
     MET_Ey += muon2->py();
      
     double MET_Et = sqrt(MET_Ex*MET_Ex + MET_Ey*MET_Ey);
     math::XYZTLorentzVector myMET(MET_Ex,MET_Ey,0.,MET_Et);

     // Change coordinates to the Z in the x^{prime} axis
     double cosZ = myZ.px()/myZ.pt();
     double sinZ = myZ.py()/myZ.pt();
     double MET_Paral = myMET.px()*cosZ + myMET.py()*sinZ;
     double MET_Perp = -myMET.px()*sinZ + myMET.py()*cosZ;

     double scaleZtoW = wMass_/zMass_;


     hZmass->Fill(myZ.mass());
     hVBPt->Fill(myZ.pt()*scaleZtoW);
     hVBPt2->Fill(myZ.pt()*scaleZtoW);
     hVBPhi->Fill(myZ.phi());
     hVBEta->Fill(myZ.eta());

     hMET->Fill(MET_Et*scaleZtoW);
     hMETParal->Fill(MET_Paral*scaleZtoW);
     hMETPerp->Fill(MET_Perp);
     hMETPhi->Fill(myMET.phi());

     hMETvsVBPt->Fill(myZ.pt()*scaleZtoW,MET_Et*scaleZtoW);
     hMETPhiMinusVBPhi->Fill(myMET.phi() - myZ.phi());
     hMETPhivsVBPhi->Fill(myZ.phi(),myMET.phi());
     hMETParalvsVBPt->Fill(myZ.pt()*scaleZtoW,MET_Paral*scaleZtoW);
     hMETPerpvsVBPt->Fill(myZ.pt()*scaleZtoW,MET_Perp);

     } // Ends of check over event selection


    // For Isolation T&P: we need to drop the cuts on isolation for the muons in the Z Candidate selection.
     if(eventSelected_foriso){
        hptMu->Fill(muon1->pt()); hptMu->Fill(muon2->pt());
        hIso->Fill(isovar1);
        hIso->Fill(isovar2);
        if(iso1){
                hIsoAll_ptMu->Fill(muon1->pt());
                hProbe_Iso->Fill(isovar1);
                if(iso2){ hProbe_ptMu->Fill(muon2->pt());}
        }
        if(iso2){
                hIsoAll_ptMu->Fill(muon2->pt());
                hProbe_Iso->Fill(isovar2);
                if(iso1){hProbe_ptMu->Fill(muon1->pt());}
        }

      } // Ends of check over event selection 

} // End of analyzer


void 
ZTemplateForSampling::beginJob()
{
   nEvts=0;  

   //const int nbins = 14;
//  double pt_bin[15]= {0, 2, 4, 6, 8, 10, 13, 16, 20, 25, 30, 40,
//  60,100,200};
   const int nbins = 15;
  double pt_bin[16]= {0, 2, 4, 6, 8, 10, 13, 16, 20, 25, 30, 40, 60,100,150,400};

  /*
  pt_bin[0]=0; 
  for (int i=1; i<52; i++){
            if(i<20)           pt_bin[i]=pt_bin[i-1]+0.5;
            if(i>=20 && i<30)  pt_bin[i]=pt_bin[i-1]+1;
            if(i>=30 && i<40)  pt_bin[i]=pt_bin[i-1]+2;
            if(i>=40)          pt_bin[i]=pt_bin[i-1]+5;
      printf(" %4.d --> %4.f \n",i,pt_bin[i]);
  }
    */  
   

   edm::Service<TFileService> fs;
   hZmass = fs->make<TH1D>("hZmass","hZmass",100, 70., 110.);
   hVBPt2  = fs->make<TH1D>("hVBPt2","hVBPt",400,0,400);   
   hVBPt = fs->make<TH1D>("hVBPt","hVBPt",nbins,pt_bin);
   hVBPhi = fs->make<TH1D>("hVBPhi","hVBPhi",101, -M_PI, M_PI);
   hVBEta = fs->make<TH1D>("hVBEta","hVBEta",101, -5.0, 5.0);
   hMET = fs->make<TH1D>("hMET","hMET",4000, 0., 400.);
   hMETParal = fs->make<TH1D>("hMETParal","hMETParal",476, -75., 400.);
   hMETPerp = fs->make<TH1D>("hMETPerp","hMETPerp",201, -100., 100.);
   hMETPhi = fs->make<TH1D>("hMETPhi","hMETPhi",101, -M_PI, M_PI);
   hMETPhiMinusVBPhi = fs->make<TH1D>("hMETPhiMinusVBPhi","hMETPhiMinusVBPhi",101, -M_PI, M_PI);
   hMETvsVBPt = fs->make<TH2D>("hMETvsVBPt","hMETvsVBPt",nbins,pt_bin,400, 0., 400.);
   hMETPhivsVBPhi  = fs->make<TH2D>("hMETPhivsVBPhi","hMETPhivsVBPhi",101, -M_PI,M_PI,101, -M_PI,M_PI);
   hMETParalvsVBPt  = fs->make<TH2D>("hMETParalvsVBPt","hMETParalvsVBPt",nbins,pt_bin,476, -75., 400.);
   hMETPerpvsVBPt  = fs->make<TH2D>("hMETPerpvsVBPt","hMETPerpvsVBPt",nbins,pt_bin,201, -100., 100.);
   /*

   edm::Service<TFileService> fs;
   hZmass = fs->make<TH1D>("hZmass","hZmass",100, 60., 120.);
   hVBPt = fs->make<TH1D>("hVBPt","hVBPt",nbins,pt_bin);
   hVBPt2  = fs->make<TH1D>("hVBPt2","hVBPt",100,0,100);   
   hVBPhi = fs->make<TH1D>("hVBPhi","hVBPhi",100, -3.141593, 3.141593);
   hVBEta = fs->make<TH1D>("hVBEta","hVBEta",100, -5.0, 5.0);
   hMET = fs->make<TH1D>("hMET","hMET",100, 0., 40.);
   hMETParal = fs->make<TH1D>("hMETParal","hMETParal",100, -20., 60.);
   hMETPerp = fs->make<TH1D>("hMETPerp","hMETPerp",100, -40., 40.);
   hMETPhi = fs->make<TH1D>("hMETPhi","hMETPhi",100, -3.141593, 3.141593);
   hMETPhiMinusVBPhi = fs->make<TH1D>("hMETPhiMinusVBPhi","hMETPhiMinusVBPhi",100, -3.141593, 3.141593);


   hMETvsVBPt = fs->make<TH2D>("hMETvsVBPt","hMETvsVBPt",nbins,pt_bin,50, 0., 80.);
   hMETPhivsVBPhi  = fs->make<TH2D>("hMETPhivsVBPhi","hMETPhivsVBPhi",100, -3.141593, 3.141593, 100, -3.141593, 3.141593);
   hMETParalvsVBPt  = fs->make<TH2D>("hMETParalvsVBPt","hMETParalvsVBPt",nbins,pt_bin,50, -20., 60.);
   hMETPerpvsVBPt  = fs->make<TH2D>("hMETPerpvsVBPt","hMETPerpvsVBPt",nbins,pt_bin,50, -40., 40.);
   */
// Isolation plots
   
   hptMu = fs->make<TH1D>("hptMu","All Muons",100, 0., 100.);
   hIsoAll_ptMu = fs->make<TH1D>("hISOAll_ptMu","Pt of Isolated Muons (no t&P)", 100, 0., 100.);
   hProbe_ptMu = fs->make<TH1D>("hProbe_ptMu","T&P Probes",100, 0., 100.);
   hIso= fs->make<TH1D>("hIso","Isolation",1000,0,1);
   hProbe_Iso= fs->make<TH1D>("hProbe_Iso","Isolation of Probes",1000,0,1);
}

void 
ZTemplateForSampling::endJob() {
}


// reco::Track ZTemplateForSampling::GetCockTailMu(const reco::MuonRef& mu) {

//     reco::TrackRef cocktail = 
//         muon::tevOptimized(mu->combinedMuon(), mu->track(), 
//                            *tevMap_default, *tevMap_1stHit, 
//                            *tevMap_picky);
     
//     return (*cocktail);

// }


//define this as a plug-in
DEFINE_FWK_MODULE(ZTemplateForSampling);
