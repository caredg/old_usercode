import FWCore.ParameterSet.Config as cms

process = cms.Process("Demo2")

#The global tag is needed if the option extractL1Prescales is set to True
#Just switch to False if this is not needed or if there is some problem
process.load('Configuration/StandardSequences/FrontierConditions_GlobalTag_cff')
process.GlobalTag.globaltag = 'GR10_H_V9::All'
process.GlobalTag.connect   = 'frontier://FrontierProd/CMS_COND_31X_GLOBALTAG'
process.GlobalTag.pfnPrefix = cms.untracked.string('frontier://FrontierProd/')


process.source = cms.Source("PoolSource",
    # replace input file with one you want to use
 #fileNames = cms.untracked.vstring('file:/home/cleonido/wprime/Summer09MC/W/F866EB68-E39C-DE11-8F01-00145EDD7879.root')
# fileNames = cms.untracked.vstring('file:/tmp/WprimeMuSkim_46_3_gAm.root')
                            fileNames = cms.untracked.vstring(
                               # 'file:/localdata/data_repo/wprime_munu/temp/WprimeMuSkim_89_1_Uqm.root',
                                                              'file:/localdata/data_repo/wprime_munu/temp/WprimeMuSkim_9_1_XYU.root')
                            #fileNames = cms.untracked.vstring('file:/localdata/data_repo/wprime_munu/temp/wprime_RECO_1000GeV_2_1_ZXU.root')
#                            fileNames = cms.untracked.vstring('file:/home/work_jarvis/data_2010B_fix.root')
)

# Number of events to process
process.maxEvents = cms.untracked.PSet(
    #input = cms.untracked.int32(-1)
    input = cms.untracked.int32(100)
)

# configuration for Wprime muon reconstruction
process.StdMu = cms.EDAnalyzer("Wprime_muonreco",
    # input tags defined here
    pvTag = cms.InputTag("offlinePrimaryVertices"),
    pvBSTag = cms.InputTag("offlinePrimaryVerticesWithBS"),
    MuonTag = cms.InputTag("muons"),
    tevMuonLabel = cms.string("tevMuons"),
    pfMetTag = cms.InputTag("pfMet"),
    pfJetTag = cms.InputTag("ak5PFJets"),
    caloJetTag = cms.InputTag("ak5CaloJets"),
    TkIsoMapTag = cms.InputTag("muIsoDepositTk"),
    EcalIsoMapTag = cms.InputTag("muIsoDepositCalByAssociatorTowers","ecal"),
    HcalIsoMapTag = cms.InputTag("muIsoDepositCalByAssociatorTowers","hcal"),

    HLTriggerResults = cms.InputTag("TriggerResults","","HLT"),
    #	InputTag HLTriggerResults = "TriggerResults::HLT2"

    #this needs a global tag above.
    #By default, the HLT prescales are extracted
    extractL1Prescales = cms.bool(True),
    triggerConditions = cms.vstring('HLT_L1MuOpen',
                                    'HLT_L1MuOpen_v*',
                                    'HLT_L1Mu',
                                    'HLT_L1Mu_v*',
                                    'HLT_L2Mu5',
                                    'HLT_L2Mu5_v*',
                                    'HLT_L2Mu9'
                                    'HLT_L2Mu9_v*'
                                    'HLT_L2Mu11'
                                    'HLT_L2Mu11_v*'
                                    'HLT_L2Mu15'
                                    'HLT_L2Mu15_v*'
                                    'HLT_L2Mu25'
                                    'HLT_L2Mu25_v*'
                                    'HLT_Mu3',
                                    'HLT_Mu3_v*',
                                    'HLT_Mu5',
                                    'HLT_Mu5_v*',
                                    'HLT_Mu7',
                                    'HLT_Mu7_v*',
                                    'HLT_Mu9',
                                    'HLT_Mu9_v*',
                                    'HLT_Mu11',
                                    'HLT_Mu11_v*',
                                    'HLT_Mu13',
                                    'HLT_Mu13_v*',
                                    'HLT_Mu15',
                                    'HLT_Mu15_v*',
                                    'HLT_Mu17',
                                    'HLT_Mu17_v*',
                                    'HLT_Mu19',
                                    'HLT_Mu19_v*',
                                    'HLT_Mu21',
                                    'HLT_Mu21_v*',
                                    'HLT_IsoMu9',
                                    'HLT_IsoMu9_v*',
                                    'HLT_IsoMu11',
                                    'HLT_IsoMu11_v*',
                                    'HLT_IsoMu13',
                                    'HLT_IsoMu13_v*',
                                    'HLT_IsoMu15',
                                    'HLT_IsoMu15_v*',
                                    'HLT_IsoMu17',
                                    'HLT_IsoMu17_v*',
                                    'HLT_IsoMu19',
                                    'HLT_IsoMu19_v*',
                                    'HLT_IsoMu21',
                                    'HLT_IsoMu21_v*',
                                    ),   


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

