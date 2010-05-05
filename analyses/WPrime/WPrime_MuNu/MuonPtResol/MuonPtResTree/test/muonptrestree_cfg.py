import FWCore.ParameterSet.Config as cms

process = cms.Process("MUPTRES")

process.load("FWCore.MessageService.MessageLogger_cfi")

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )
process.MessageLogger.cerr.FwkReport.reportEvery = 1

process.source = cms.Source("PoolSource",
    # replace 'myfile.root' with the source file you want to use
    fileNames = cms.untracked.vstring(
        'file:mytrajectories.root'
    )
)

process.MuPtRes = cms.EDAnalyzer('MuonPtResTree',
                                 
                                 # input tags defined here
                                 MuonTag = cms.InputTag("muons"),
                                 
                                 )

process.TFileService = cms.Service("TFileService",
                                   fileName = cms.string('mures_rootuple.root')
                                   )

process.p = cms.Path(process.MuPtRes)
