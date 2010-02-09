import FWCore.ParameterSet.Config as cms

process = cms.Process("Demo")

process.load("FWCore.MessageService.MessageLogger_cfi")

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

process.source = cms.Source("PoolSource",
    # replace 'myfile.root' with the source file you want to use
    fileNames = cms.untracked.vstring(
    #zmumu
      # 'file:9A42A6DC-63B1-DE11-933F-000423D99F1E.root',
        #wmu
       #         'file:D0F6EDD7-72B1-DE11-AC81-001D09F29619.root'

'/store/relval/CMSSW_3_1_4/RelValZMM/GEN-SIM-RECO/STARTUP31X_V2-v1/0006/9A42A6DC-63B1-DE11-933F-000423D99F1E.root',
       '/store/relval/CMSSW_3_1_4/RelValZMM/GEN-SIM-RECO/STARTUP31X_V2-v1/0005/BE87B5F5-DCB0-DE11-A83A-001D09F23174.root',
#       '/store/relval/CMSSW_3_1_4/RelValZMM/GEN-SIM-RECO/STARTUP31X_V2-v1/0005/72574DE2-D9B0-DE11-B6B9-001D09F2512C.root',
#       '/store/relval/CMSSW_3_1_4/RelValZMM/GEN-SIM-RECO/STARTUP31X_V2-v1/0005/2A6D41CE-DBB0-DE11-A14F-0019B9F70468.root'

#'/store/relval/CMSSW_3_1_4/RelValWM/GEN-SIM-RECO/STARTUP31X_V2-v1/0006/D0F6EDD7-72B1-DE11-AC81-001D09F29619.root',
#       '/store/relval/CMSSW_3_1_4/RelValWM/GEN-SIM-RECO/STARTUP31X_V2-v1/0006/C004562D-3FB1-DE11-8FA3-001D09F24DDF.root',
#       '/store/relval/CMSSW_3_1_4/RelValWM/GEN-SIM-RECO/STARTUP31X_V2-v1/0006/9C49D64F-3AB1-DE11-9344-001D09F24FEC.root',
#       '/store/relval/CMSSW_3_1_4/RelValWM/GEN-SIM-RECO/STARTUP31X_V2-v1/0006/507D35D3-3DB1-DE11-8C24-001D09F2A465.root'

#Wjets 80-120
#'/store/relval/CMSSW_3_1_4/RelValWjet_Pt_80_120/GEN-SIM-RECO/MC_31X_V3-v1/0006/5AF66D47-72B1-DE11-BF58-001D09F2B2CF.root',
#       '/store/relval/CMSSW_3_1_4/RelValWjet_Pt_80_120/GEN-SIM-RECO/MC_31X_V3-v1/0005/DE107E07-5FB0-DE11-9F09-000423D990CC.root',
#       '/store/relval/CMSSW_3_1_4/RelValWjet_Pt_80_120/GEN-SIM-RECO/MC_31X_V3-v1/0005/7A7D8AE4-5DB0-DE11-80AB-001D09F29524.root',
#       '/store/relval/CMSSW_3_1_4/RelValWjet_Pt_80_120/GEN-SIM-RECO/MC_31X_V3-v1/0005/624E8518-5DB0-DE11-9742-001D09F29619.root',
#       '/store/relval/CMSSW_3_1_4/RelValWjet_Pt_80_120/GEN-SIM-RECO/MC_31X_V3-v1/0005/342261A8-5FB0-DE11-A171-000423D998BA.root'


#qcd 80-120
#'/store/relval/CMSSW_3_1_4/RelValQCD_Pt_80_120/GEN-SIM-RECO/START_31X_V4A-v1/0000/F89B1E02-B1B9-DE11-A7AA-001D09F2447F.root',
#       '/store/relval/CMSSW_3_1_4/RelValQCD_Pt_80_120/GEN-SIM-RECO/START_31X_V4A-v1/0000/6C7257CF-ABB9-DE11-9EC5-001D09F2960F.root',
#       '/store/relval/CMSSW_3_1_4/RelValQCD_Pt_80_120/GEN-SIM-RECO/START_31X_V4A-v1/0000/30CE260E-C4B9-DE11-9AC8-001D09F28F0C.root',
#       '/store/relval/CMSSW_3_1_4/RelValQCD_Pt_80_120/GEN-SIM-RECO/START_31X_V4A-v1/0000/127D5ED5-A7B9-DE11-890E-001D09F2437B.root'

#ttbar
#'/store/relval/CMSSW_3_1_4/RelValTTbar/GEN-SIM-RECO/STARTUP31X_V2-v1/0006/AC0641BB-73B1-DE11-A138-001D09F291D2.root',
#       '/store/relval/CMSSW_3_1_4/RelValTTbar/GEN-SIM-RECO/STARTUP31X_V2-v1/0005/EE8AA50B-CCB0-DE11-843B-000423D94E70.root',
#       '/store/relval/CMSSW_3_1_4/RelValTTbar/GEN-SIM-RECO/STARTUP31X_V2-v1/0005/E6E11BD7-C7B0-DE11-A4A8-001D09F291D2.root',
#       '/store/relval/CMSSW_3_1_4/RelValTTbar/GEN-SIM-RECO/STARTUP31X_V2-v1/0005/CC21C5DB-C6B0-DE11-994E-000423D98868.root',
#       '/store/relval/CMSSW_3_1_4/RelValTTbar/GEN-SIM-RECO/STARTUP31X_V2-v1/0005/B2EED459-C7B0-DE11-9CC1-001D09F29597.root',
#       '/store/relval/CMSSW_3_1_4/RelValTTbar/GEN-SIM-RECO/STARTUP31X_V2-v1/0005/6668B3E0-C5B0-DE11-9A9F-001D09F28F0C.root',
#       '/store/relval/CMSSW_3_1_4/RelValTTbar/GEN-SIM-RECO/STARTUP31X_V2-v1/0005/2A1AC0DD-C8B0-DE11-8689-001D09F276CF.root',
#       '/store/relval/CMSSW_3_1_4/RelValTTbar/GEN-SIM-RECO/STARTUP31X_V2-v1/0005/0E1BD752-C6B0-DE11-BEAB-001D09F28EC1.root' 
        
    )
)

process.maxEvents = cms.untracked.PSet(
#     input = cms.untracked.int32(-1)
   input = cms.untracked.int32(400)
)


process.demo = cms.EDAnalyzer('ZmumuTreeMaker',
                              HLTriggerResults = cms.InputTag("TriggerResults","","HLT")
)


process.p = cms.Path(process.demo)
