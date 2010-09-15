
import FWCore.ParameterSet.Config as cms

process = cms.Process('MURESO')

# import of standard configurations
process.load('Configuration.StandardSequences.Services_cff')
process.load('SimGeneral.HepPDTESSource.pythiapdt_cfi')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.StandardSequences.MixingNoPileUp_cff')
process.load('Configuration.StandardSequences.GeometryDB_cff')
process.load('Configuration.StandardSequences.MagneticField_38T_cff')
process.load('Configuration.StandardSequences.Reconstruction_cff')
process.load('Configuration.StandardSequences.EndOfProcess_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
process.load('Configuration.EventContent.EventContent_cff')

process.configurationMetadata = cms.untracked.PSet(
    version = cms.untracked.string('$Revision: 1.5 $'),
    annotation = cms.untracked.string('recoTeVMuon nevts:10'),
    name = cms.untracked.string('PyReleaseValidation')
)
process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(-1)
)
process.options = cms.untracked.PSet(

)
# Input source
process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring('/castor/cern.ch/user/e/ecarrera/analyses/wprime_munu/singleMuonEvts/382/SingleMuon_RECO_1000_START38_V10_20000evts.root')
           
)

# Output definition
process.output = cms.OutputModule("PoolOutputModule",
    splitLevel = cms.untracked.int32(0),
    outputCommands = process.RECOSIMEventContent.outputCommands,
    fileName = cms.untracked.string('file:myresolution.root'),
    dataset = cms.untracked.PSet(
        dataTier = cms.untracked.string('GEN-SIM-RECO'),
        filterName = cms.untracked.string('')
    )
)

# Additional output definition

# Other statements
process.GlobalTag.globaltag = 'START38_V10::All'


#Trying to access localPosition for a RecHit that was
#read from disk, but since CMSSW_2_1_X local positions are transient.
#If you want to get coarse position/error estimation from disk,
#please set: ComputeCoarseLocalPositionFromDisk = True 
#to the TransientTrackingRecHitBuilder you are using
#from RecoTracker/TransientTrackingRecHit/python/TTRHBuilders_cff.py
process.ttrhbwr.ComputeCoarseLocalPositionFromDisk = True

#check the refitter without any change
process.tevMuonsNew = process.tevMuons.clone()
process.tevMuonsNew.RefitterParameters.HitsToKeep = cms.untracked.string("all")
process.tevMuonsNew.RefitterParameters.NewSeed = cms.untracked.string("none")
process.tevMuonsNew.RefitterParameters.RefitDirection = 'RinsideOut'
#refit with odd hits
process.tevMuonsOdd = process.tevMuons.clone()
process.tevMuonsOdd.RefitterParameters.HitsToKeep = cms.untracked.string("odd")
process.tevMuonsOdd.RefitterParameters.NewSeed = cms.untracked.string("pairORtriplet")
process.tevMuonsOdd.RefitterParameters.RefitDirection = 'RinsideOut'
#refit with even hits
process.tevMuonsEven = process.tevMuons.clone()
process.tevMuonsEven.RefitterParameters.HitsToKeep = cms.untracked.string("even")
process.tevMuonsEven.RefitterParameters.NewSeed = cms.untracked.string("pairORtriplet")
process.tevMuonsEven.RefitterParameters.RefitDirection = 'RinsideOut'

#save everything, no need to use pre-defined stuff
process.output.outputCommands = cms.untracked.vstring('keep *')

# Path and EndPath definitions
#process.reconstruction_step = cms.Path(process.reconstruction)
#process.reconstruction_step = cms.Path(process.tevMuons)
process.reconstruction_step = cms.Path(process.tevMuonsNew+process.tevMuonsOdd+process.tevMuonsEven)
process.endjob_step = cms.Path(process.endOfProcess)
process.out_step = cms.EndPath(process.output)

# Schedule definition
process.schedule = cms.Schedule(process.reconstruction_step,process.endjob_step,process.out_step)


