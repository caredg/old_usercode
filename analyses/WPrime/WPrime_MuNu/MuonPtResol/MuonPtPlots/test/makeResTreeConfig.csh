#!/bin/tcsh

if ( $#argv < 3 ) then
  echo "Usage: source gun_some_muons.csh <PTMUON> <GT> <NEVT>"
  exit 1
endif

set PTMUON = $1
set GT = $2
set NEVT = $3

/bin/rm -f muonrestree_${PTMUON}_${GT}_cfg.py

cat >  muonrestree_${PTMUON}_${GT}_cfg.py <<EOF

import FWCore.ParameterSet.Config as cms

process = cms.Process("MYTREE")

process.load("FWCore.MessageService.MessageLogger_cfi")

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(${NEVT}) )
process.MessageLogger.cerr.FwkReport.reportEvery = 1

process.source = cms.Source("PoolSource",
    # replace 'myfile.root' with the source file you want to use
    fileNames = cms.untracked.vstring(
        'file:newMuonTracksMURES_${PTMUON}_${GT}.root'
    ),
    duplicateCheckMode = cms.untracked.string("noDuplicateCheck")

)

process.myplots = cms.EDAnalyzer('MuonPtPlots',
                                 MuonPt = cms.double(${PTMUON}),
                                 FileName = cms.string("murestree_${PTMUON}_${GT}.root")
)


process.p = cms.Path(process.myplots)

EOF

cmsenv
cmsRun muonrestree_${PTMUON}_${GT}_cfg.py >&! log_murestree_${PTMUON}_${GT}.log
