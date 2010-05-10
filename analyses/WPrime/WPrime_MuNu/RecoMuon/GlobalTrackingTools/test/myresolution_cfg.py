import FWCore.ParameterSet.Config as cms

process = cms.Process("MURES")

refitJitterScale = 100.0
outfilename ="newMuonTracksMURES.root" 


from RecoMuon.TrackingTools.MuonServiceProxy_cff import MuonServiceProxy
from RecoMuon.TrackingTools.MuonTrackLoader_cff import MuonUpdatorAtVertexAnyDirection
from RecoMuon.GlobalTrackingTools.GlobalMuonRefitter_cff import GlobalMuonRefitter
from RecoTracker.TrackProducer.TrackRefitterP5_cfi import TrackRefitterP5

#process.load("Alignment.CommonAlignmentProducer.GlobalPosition_Fake_cff")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
process.GlobalTag.globaltag = 'MC_3XY_V9A::All'

process.load("Configuration.EventContent.EventContent_cff")
process.load("Configuration.StandardSequences.Generator_cff")
process.load("Configuration.StandardSequences.Geometry_cff")
process.load("Configuration.StandardSequences.MagneticField_cff")
#process.load("Configuration.StandardSequences.MixingNoPileUp_cff")
process.load("Configuration.StandardSequences.Reconstruction_cff")

process.load("Geometry.CMSCommonData.cmsIdealGeometryXML_cff")
process.load("Geometry.CommonDetUnit.globalTrackingGeometry_cfi")
process.load("Geometry.MuonNumbering.muonNumberingInitialization_cfi")
process.load("Configuration.StandardSequences.Geometry_cff")
process.load("Configuration.StandardSequences.MagneticField_cff")

process.load("RecoMuon.TrackingTools.MuonServiceProxy_cff")
process.load("RecoMuon.TrackingTools.MuonTrackLoader_cff")


process.load("FWCore.MessageService.MessageLogger_cfi")
#process.load("RecoMuon.GlobalMuonProducer.tevMuons_cfi")
process.load("RecoLocalTracker.Configuration.RecoLocalTracker_cff")

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(500) )

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
    'file:/localdata/ecarrera/analyses/wprime_to_munu/trackPt_resolution/singleMu_data/singleMuPt100/0AB2FE6B-2DE4-DE11-9E2C-0018F3D0969C.root'
    )
)


process.generalTracks1= cms.EDProducer("TevMuonProducer",
    process.MuonTrackLoaderForGLB,
    #    InputTag MuonCollectionLabel = standAloneMuons:UpdatedAtVtx
    MuonServiceProxy,
    RefitIndex = cms.vint32(1),
    Refits = cms.vstring('refit1'),
    MuonCollectionLabel = cms.InputTag("generalTracks"),
    RefitterParameters = cms.PSet(
        GlobalMuonRefitter,
    	StabStab = cms.untracked.int32(0),
    	StabScale= cms.untracked.double(refitJitterScale),
    )
)
process.generalTracksEve1= cms.EDProducer("TevMuonProducer",
    process.MuonTrackLoaderForGLB,
    #    InputTag MuonCollectionLabel = standAloneMuons:UpdatedAtVtx
    MuonServiceProxy,
    RefitIndex = cms.vint32(1),
    Refits = cms.vstring('eve'),
    MuonCollectionLabel = cms.InputTag("generalTracks"),
    RefitterParameters = cms.PSet(
        GlobalMuonRefitter,
    	StabStab = cms.untracked.int32(12),
    	StabScale= cms.untracked.double(refitJitterScale),
    )
)
process.generalTracksOdd1= cms.EDProducer("TevMuonProducer",
    process.MuonTrackLoaderForGLB,
    #    InputTag MuonCollectionLabel = standAloneMuons:UpdatedAtVtx
    process.MuonServiceProxy,
    RefitIndex = cms.vint32(1),
    Refits = cms.vstring('odd'), 
	MuonCollectionLabel = cms.InputTag("generalTracks"),
    RefitterParameters = cms.PSet(
        GlobalMuonRefitter,
    	StabStab = cms.untracked.int32(13),
    	StabScale= cms.untracked.double(refitJitterScale),
    )
)
process.generalTracks1.MuonCollectionLabel = cms.InputTag("generalTracks")
process.generalTracksEve1.MuonCollectionLabel = cms.InputTag("generalTracks")
process.generalTracksOdd1.MuonCollectionLabel = cms.InputTag("generalTracks")




process.load("RecoTracker.TransientTrackingRecHit.TTRHBuilders_cff")
##
process.WithTrackAngle1= cms.ESProducer("TkTransientTrackingRecHitBuilderESProducer",
	StripCPE = cms.string('StripCPEfromTrackAngle'),
	ComponentName = cms.string('WithTrackAngle1'),
	PixelCPE = cms.string('PixelCPEGeneric'),
	Matcher = cms.string('StandardMatcher'),
	ComputeCoarseLocalPositionFromDisk = cms.bool(True)
)

process.generalTracks1.RefitterParameters.TrackerRecHitBuilder = cms.string('WithTrackAngle1')
process.generalTracksEve1.RefitterParameters.TrackerRecHitBuilder = cms.string('WithTrackAngle1')
process.generalTracksOdd1.RefitterParameters.TrackerRecHitBuilder = cms.string('WithTrackAngle1')

#the new refits:
process.tevMuonsNew = process.generalTracks1.clone()
process.tevMuonsNew.MuonCollectionLabel = cms.InputTag("globalMuons")
process.tevMuonsNew.RefitIndex = cms.vint32(2)
process.tevMuonsNew.RefitterParameters.StabStab = 0
#nojitter
process.tevMuonsNewNoJit = process.generalTracks1.clone()
process.tevMuonsNewNoJit.MuonCollectionLabel = cms.InputTag("globalMuons")
process.tevMuonsNewNoJit.RefitIndex = cms.vint32(2)
process.tevMuonsNewNoJit.RefitterParameters.StabStab = 0
process.tevMuonsNewNoJit.RefitterParameters.StabScale = 1
#with odd/even hits in the tracker and odd/even hits in the first
#muon chamber with hits
process.tevMuonsOddHh = process.generalTracksOdd1.clone()
process.tevMuonsEveHh = process.generalTracksEve1.clone()
process.tevMuonsOddHh.MuonCollectionLabel = cms.InputTag("globalMuons")
process.tevMuonsEveHh.MuonCollectionLabel = cms.InputTag("globalMuons")
process.tevMuonsOddHh.RefitIndex = cms.vint32(2)
process.tevMuonsEveHh.RefitIndex = cms.vint32(2)
process.tevMuonsOddHh.RefitterParameters.StabStab = 13
process.tevMuonsEveHh.RefitterParameters.StabStab = 12
#with odd/even tracker layers and odd/even hits in first muon chamber
#with hits
process.tevMuonsOddLh= process.generalTracksOdd1.clone()
process.tevMuonsEveLh = process.generalTracksEve1.clone()
process.tevMuonsOddLh.MuonCollectionLabel = cms.InputTag("globalMuons")
process.tevMuonsEveLh.MuonCollectionLabel = cms.InputTag("globalMuons")
process.tevMuonsOddLh.RefitIndex = cms.vint32(2)
process.tevMuonsEveLh.RefitIndex = cms.vint32(2)
process.tevMuonsOddLh.RefitterParameters.StabStab = 131
process.tevMuonsEveLh.RefitterParameters.StabStab = 121







#process.preFit = cms.Path(process.siPixelRecHits*process.siStripMatchedRecHits*process.recopixelvertexing*process.ckftracks
##	process.reconstruction
#	+process.generalTracks1*process.generalTracksEve1*process.generalTracksOdd1
#	*process.generalTracksEve2*process.generalTracksOdd2
#	*process.globalTracks1*process.globalTracksEve1*process.globalTracksOdd1*process.tevTracksEve1*process.tevTracksOdd1)


process.reFit = cms.Path(process.tevMuonsNew*process.tevMuonsNewNoJit+process.tevMuonsOddHh*process.tevMuonsEveHh+process.tevMuonsOddLh*process.tevMuonsEveLh)



process.o1 = cms.OutputModule("PoolOutputModule",
  fileName = cms.untracked.string("myfile.root"),
)
process.o1.fileName = outfilename 
process.outpath = cms.EndPath(process.o1)

process.source.inputCommands = cms.untracked.vstring(
	 "keep *")

process.schedule=cms.Schedule(process.reFit,process.outpath)




process.MessageLogger.cerr.FwkReport.reportEvery = 1

#activate debugging
#process.load('logger_cfi')
#process.MessageLogger._categoryCanTalk('RecoMuon')
#process.MessageLogger._categoryCanTalk('TrackTransformer')
#process.MessageLogger._categoryCanTalk('TrackFitters')
#process.MessageLogger._moduleCanTalk('tevTracksEve1')

