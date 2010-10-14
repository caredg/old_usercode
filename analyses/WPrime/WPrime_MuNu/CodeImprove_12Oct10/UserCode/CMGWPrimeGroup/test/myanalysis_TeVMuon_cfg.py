import FWCore.ParameterSet.Config as cms

process = cms.Process("Demo2")
process.source = cms.Source("PoolSource",
    # replace input file with one you want to use
 #fileNames = cms.untracked.vstring('file:/home/cleonido/wprime/Summer09MC/W/F866EB68-E39C-DE11-8F01-00145EDD7879.root')
# fileNames = cms.untracked.vstring('file:/tmp/WprimeMuSkim_46_3_gAm.root')
                            fileNames = cms.untracked.vstring('file:/home/work_jarvis/data_2010B_fix.root')
)

# Number of events to process
process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(-1)
#    input = cms.untracked.int32(100)
)

# configuration for Wprime muon reconstruction
process.StdMu = cms.EDAnalyzer("Wprime_muonreco",
    # input tags defined here
    pvTag = cms.InputTag("offlinePrimaryVertices"),
    pvBSTag = cms.InputTag("offlinePrimaryVerticesWithBS"),
    MuonTag = cms.InputTag("muons"),
    pfMetTag = cms.InputTag("pfMet"),
    pfJetTag = cms.InputTag("ak5PFJets"),
    caloJetTag = cms.InputTag("ak5CaloJets"),
    TkIsoMapTag = cms.InputTag("muIsoDepositTk"),
    EcalIsoMapTag = cms.InputTag("muIsoDepositCalByAssociatorTowers","ecal"),
    HcalIsoMapTag = cms.InputTag("muIsoDepositCalByAssociatorTowers","hcal"),

    HLTriggerResults = cms.InputTag("TriggerResults","","HLT"),
     #	InputTag HLTriggerResults = "TriggerResults::HLT2"

   # sample description
    description = cms.string('Single Muon 7 TeV data skim with pt > 10 GeV'),

    # muon-detector eta acceptance
    Detmu_acceptance = cms.double(2.4)

)

process.TFileService = cms.Service("TFileService",
   fileName = cms.string('HighPtMuon_7TeV_Run2010.root')
#   fileName = cms.string('Wmu_31x_IDEAL_7TeV_1.root')
#   fileName = cms.string('wprime_1TeV.root')
)

process.p = cms.Path(process.StdMu)

