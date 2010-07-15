import FWCore.ParameterSet.Config as cms

process = cms.Process("MYPLOTS")

process.load("FWCore.MessageService.MessageLogger_cfi")

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )
process.MessageLogger.cerr.FwkReport.reportEvery = 1

process.source = cms.Source("PoolSource",
    # replace 'myfile.root' with the source file you want to use
    fileNames = cms.untracked.vstring(
        'file:newMuonTracksMURES.root'
    ),
    duplicateCheckMode = cms.untracked.string("noDuplicateCheck")

)

process.myplots = cms.EDAnalyzer('MuonPtPlots'
)


process.p = cms.Path(process.myplots)
