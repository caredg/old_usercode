import FWCore.ParameterSet.Config as cms
import copy

process = cms.Process("ZToMuMuHistos")

process.load("FWCore.MessageLogger.MessageLogger_cfi")
#process.load("UserCode.MCepeda.GoodData_cfi")

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(

    #'/store/data/Run2010B/Mu/RECO/PromptReco-v2/000/146/511/E0E62F3F-8FC7-DF11-B66F-0019B9F709A4.root'
'file:/localdata/data_repo/wprime_munu/temp/WprimeMuSkim_89_1_Uqm.root' 

)
)


process.options = cms.untracked.PSet(
    wantSummary = cms.untracked.bool(True)
)

process.MessageLogger = cms.Service("MessageLogger",
    debugModules = cms.untracked.vstring('met_corrected'),
    categories = cms.untracked.vstring('FwkJob','FwkReport','FwkSummary','Root_NoDictionary'),
    cout = cms.untracked.PSet(
        threshold = cms.untracked.string('ERROR')
    ),
    destinations = cms.untracked.vstring('cout')
)

process.AdaptorConfig = cms.Service("AdaptorConfig")


process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(-1)
)

print "Number of files to process is %s" % (len(process.source.fileNames))

process.met_corrected = cms.EDAnalyzer('ZTemplateForSampling',
    METTag = cms.untracked.InputTag('corMetGlobalMuons'),
    MuonTrig = cms.vstring('HLT_Mu9',
                           'HLT_Mu11',
                           'HLT_Mu13',
                           'HLT_Mu15',
                           'HLT_Mu17',
                           'HLT_Mu19',
                           'HLT_Mu21',
                           'HLT_Mu13_v1',
                           'HLT_Mu15_v1',
                           'HLT_Mu17_v1',
                           'HLT_Mu19_v1',
                           'HLT_Mu21_v1'
                           ),
    EtaCut = cms.untracked.double(1.8)                                   
                                       )
process.pflow = process.met_corrected.clone(
    METTag = cms.untracked.InputTag('pfMet'),
)

process.tcmet = process.met_corrected.clone(
    METTag = cms.untracked.InputTag('tcMet'),
)


#CAREFUL: MET is assumed corrected for muons!

process.TFileService = cms.Service(
    "TFileService",
    fileName = cms.string("ZMET_data.root")
)


process.endPath = cms.EndPath( 
    process.pflow*
    process.met_corrected*
    process.tcmet
)

