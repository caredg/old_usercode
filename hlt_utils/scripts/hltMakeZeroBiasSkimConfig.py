#!/usr/bin/env python
############################################################################
#
# Edgar Carrera
# ecarrera@cern.ch
#
# This script takes as input the path for the RAW dataset of a given run
# that resides in CASTOR
# and outputs a configuration file that can be used to skim
# on the L1 "enabled" seeds that were used in that run
# requiring also that the ZeroBias tech bit 4 was fired.
# The purpose of such skim is to use it for timing studies.
#
# It also creates a *.csh file that can be used to quickly submit a LSF batch job
#
# One can use this script as:
# python hltMakeZeroBiasSkim.py -p /castor/cern.ch/cms/store/data/Run2010B/MinimumBias/RAW/v1/000/147/219/ -t 33500 -z 34.82 -o /castor/cern.ch/user/e/ecarrera/hlt_debug_skims/forHilton/run147219_L1SkimZB.root
# for example.
###########################################################################

"""
   usage: %prog [options]
   -p, --path = PATH: Full CASTOR path of the RAW dataset for the run that will be skimmed.
   -o, --out  = OUTDIR: Full CASTOR name for the output file.
   -z, --zbrate = ZBRATE: Rate in Hz for the ZB trigger, bit 4 in the run
   -t, --trate = TRATE: Total L1 rate in Hz
"""


import os,sys
import subprocess
import string, re
import fileinput
import commands
from time import gmtime, localtime, strftime


# _____________________OPTIONS_______________________________________________

############################################################################
# Code taken from http://code.activestate.com/recipes/278844/
############################################################################
import optparse
USAGE = re.compile(r'(?s)\s*usage: (.*?)(\n[ \t]*\n|$)')
def nonzero(self): # will become the nonzero method of optparse.Values
    "True if options were given"
    for v in self.__dict__.itervalues():
        if v is not None: return True
    return False

optparse.Values.__nonzero__ = nonzero # dynamically fix optparse.Values

class ParsingError(Exception): pass

optionstring=""

def exit(msg=""):
    raise SystemExit(msg or optionstring.replace("%prog",sys.argv[0]))

def parse(docstring, arglist=None):
    global optionstring
    optionstring = docstring
    match = USAGE.search(optionstring)
    if not match: raise ParsingError("Cannot find the option string")
    optlines = match.group(1).splitlines()
    try:
        p = optparse.OptionParser(optlines[0])
        for line in optlines[1:]:
            opt, help=line.split(':')[:2]
            short,long=opt.split(',')[:2]
            if '=' in opt:
                action='store'
                long=long.split('=')[0]
            else:
                action='store_true'
            p.add_option(short.strip(),long.strip(),
                         action = action, help = help.strip())
    except (IndexError,ValueError):
        raise ParsingError("Cannot parse the option string correctly")
    return p.parse_args(arglist)


#_________________________________________________________________________




#######################################################
def dump_L1_info(path):
#######################################################
    #note: we should get this information in a more efficient way
    #get rid of trash
    os.system("rm -f auxfile.py")
    os.system("rm -f l1dump.txt")
    #pick the first file from the dataset to use it as source of info
    str_nsls = "nsls "+path
    mypipe = subprocess.Popen(str_nsls,shell=True,stdout=subprocess.PIPE)
    files = mypipe.communicate()[0].split()
    pfile = path+"/"+files[0]
    #get the global tag from the Provenance
    print "Extracting GlobalTag from the Provenance, this might take some time due to CASTOR I/O, be patient or pre-stage the file..."
    str_gtag = "edmProvDump rfio:"+pfile+"|grep globaltag|head -1| awk '{print $5}'"
    mypipe = subprocess.Popen(str_gtag,shell=True,stdout=subprocess.PIPE)
    gtag = mypipe.communicate()[0].split("\'")[1]
    global GTAG
    GTAG = gtag
    #make auxiliary python config file to print
    #L1 mask for the run
    auxfile = open("auxfile.py","a")
    auxfile.write("useGlobalTag = '%s'\n" % gtag)
    auxfile.write(
"""
import FWCore.ParameterSet.Config as cms

process = cms.Process("L1SKIM")
# initialize  MessageLogger
process.load("FWCore.MessageService.MessageLogger_cfi")
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
process.load("L1Trigger.GlobalTriggerAnalyzer.l1GtTrigReport_cfi")


process.MessageLogger.categories.append('L1GtTrigReport')
process.l1GtTrigReport.PrintVerbosity = 1
process.MessageLogger.cerr.FwkReport.reportEvery = 1 




process.GlobalTag.connect   = 'frontier://FrontierProd/CMS_COND_31X_GLOBALTAG'
process.GlobalTag.pfnPrefix = cms.untracked.string('frontier://FrontierProd/')
process.GlobalTag.globaltag = useGlobalTag 



process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
                                "rfio:dummy.root"
                                )
)
process.maxEvents = cms.untracked.PSet(
   input = cms.untracked.int32(5)
)

process.e = cms.EndPath(process.l1GtTrigReport)
"""

        )
    auxfile.write("process.source.fileNames = cms.untracked.vstring('rfio:%s')" % pfile)
    auxfile.close()
    #execute the auxfile to get the dump in a file
    print "Extracting L1 info, this might take some time due to CASTOR I/O, be patient or pre-stage the file..."
    os.system("cmsRun auxfile.py >& l1dump.txt")
    #!!!!!!!!!!!!!!Temporary fixes for seeds that have different names
    #in the output file, don't know why
    #fix the Jet* seeds.  They need the "U" for uncorrected
    os.system("sed -i -e 's#Jet\([0-9]*\)#Jet\\1U#g' l1dump.txt")
    #get rid of the QE8 in the L1_Mu* seeds
    os.system("sed -i -e 's#Mu\([0-9]*\)QE8#Mu\\1#g' l1dump.txt")
    #get rid of the _Ext extension
    os.system("sed -i -e 's#_Ext##g' l1dump.txt")
    #clean up
    os.system("rm -f auxfile.py")

    
    

#######################################################
def parse_L1algo_info():
#######################################################

    algolist = []
    file = open("l1dump.txt","r")
    parseMe = False
    #loop over lines in the file
    for line in file.readlines():
        #split the line 
        sline = line.split()
        #activate the parsing to be safe
        if (line.find("Algorithm Key") != -1):
            parseMe = True
        #switch off parsing when done
        if (line.find("Technical Trigger Key") !=-1):
                parseMe = False
        #make sure to find an algo seed
        if (parseMe):
            if (line.find("L1_")!=-1):
                #check the mask
                if (sline[2] == "0"):
                    algolist.append(sline[0])
       

    return algolist
                

#######################################################
def parse_L1tt_info():
#######################################################

    techlist = []
    file = open("l1dump.txt","r")
    parseMe = False
    #loop over lines in the file
    for line in file.readlines():
        #split the line 
        sline = line.split()
        #activate the parsing to be safe
        if (line.find("Technical Trigger Key") !=-1):
            parseMe = True
        #switch off parsing when done with the tt bits
        if (len(sline) == 0):
                parseMe = False
        if (parseMe):
            #make sure to find an tt seed but skip tt4 (ZB)
            if (sline[0].isdigit() and sline[0]!= "4"):
                if (sline[2] == "0"):
                    techlist.append(sline[0])
        

    return techlist


#######################################################
def get_L1bits_list(dicOpt):
#######################################################
   
    #get path to RAW files
    path = dicOpt['path']
    #dump L1 information in a file
    dump_L1_info(path)
    #parse the resulting file
    algolist = parse_L1algo_info()
    techlist = parse_L1tt_info()

    #clean up
    os.system("rm -f l1dump.txt")
    
    return algolist, techlist


#######################################################
def print_config_file(dicOpt,algolist,techlist):
#######################################################
    path = dicOpt['path']
    out = dicOpt['out']
    zbrate = dicOpt['zbrate']
    trate = dicOpt['trate']

    #calculate prescale for ZB path
    prezb = int(float(trate)/float(zbrate))
    print "the prescale for ZB is %i" % prezb
    
    #clear config
    os.system("rm -f skimL1Seeds.py")
    #open a config file to append
    cfg = open("skimL1Seeds.py","a")
    
    #--------make the first part of the file, which is the same
    #--------for any run
    cfg.write("useGlobalTag = '%s'\n" % GTAG)
    cfg.write(
"""
import FWCore.ParameterSet.Config as cms

process = cms.Process("TSKIM")

# initialize  MessageLogger
process.load("FWCore.MessageService.MessageLogger_cfi")
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')


#Deal with the global tag
process.GlobalTag.connect   = 'frontier://FrontierProd/CMS_COND_31X_GLOBALTAG'
process.GlobalTag.pfnPrefix = cms.untracked.string('frontier://FrontierProd/')
process.GlobalTag.globaltag = useGlobalTag


#Input the source file but feed it at the end of the file
process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
                                'rfio:dummy.root'
                                )
)




#Unpack the L1 information a la HLT menu
process.hltGtDigis = cms.EDProducer( "L1GlobalTriggerRawToDigi",
    DaqGtInputTag = cms.InputTag( "source" ),
    DaqGtFedId = cms.untracked.int32( 813 ),
    ActiveBoardsMask = cms.uint32( 0xffff ),
    UnpackBxInEvent = cms.int32( 5 ),
    Verbosity = cms.untracked.int32( 0 )
)

process.hltGctDigis = cms.EDProducer( "GctRawToDigi",
    inputLabel = cms.InputTag( "source" ),
    gctFedId = cms.untracked.int32( 745 ),
    hltMode = cms.bool( True ),
    numberOfGctSamplesToUnpack = cms.uint32( 1 ),
    numberOfRctSamplesToUnpack = cms.uint32( 1 ),
    unpackSharedRegions = cms.bool( False ),
    unpackerVersion = cms.uint32( 0 )
)

process.hltL1GtObjectMap = cms.EDProducer( "L1GlobalTrigger",
    GmtInputTag = cms.InputTag( "hltGtDigis" ),
    GctInputTag = cms.InputTag( "hltGctDigis" ),
    CastorInputTag = cms.InputTag( "castorL1Digis" ),
    ProduceL1GtDaqRecord = cms.bool( False ),
    ProduceL1GtEvmRecord = cms.bool( False ),
    ProduceL1GtObjectMapRecord = cms.bool( True ),
    WritePsbL1GtDaqRecord = cms.bool( False ),
    ReadTechnicalTriggerRecords = cms.bool( True ),
    EmulateBxInEvent = cms.int32( 1 ),
    AlternativeNrBxBoardDaq = cms.uint32( 0 ),
    AlternativeNrBxBoardEvm = cms.uint32( 0 ),
    BstLengthBytes = cms.int32( -1 ),
    TechnicalTriggersInputTags = cms.VInputTag( 'simBscDigis' ),
    RecordLength = cms.vint32( 3, 0 )
)

process.hltL1extraParticles = cms.EDProducer( "L1ExtraParticlesProd",
    produceMuonParticles = cms.bool( True ),
    muonSource = cms.InputTag( "hltGtDigis" ),
    produceCaloParticles = cms.bool( True ),
    isolatedEmSource = cms.InputTag( 'hltGctDigis','isoEm' ),
    nonIsolatedEmSource = cms.InputTag( 'hltGctDigis','nonIsoEm' ),
    centralJetSource = cms.InputTag( 'hltGctDigis','cenJets' ),
    forwardJetSource = cms.InputTag( 'hltGctDigis','forJets' ),
    tauJetSource = cms.InputTag( 'hltGctDigis','tauJets' ),
    etTotalSource = cms.InputTag( "hltGctDigis" ),
    etHadSource = cms.InputTag( "hltGctDigis" ),
    etMissSource = cms.InputTag( "hltGctDigis" ),
    htMissSource = cms.InputTag( "hltGctDigis" ),
    hfRingEtSumsSource = cms.InputTag( "hltGctDigis" ),
    hfRingBitCountsSource = cms.InputTag( "hltGctDigis" ),
    centralBxOnly = cms.bool( True ),
    ignoreHtMiss = cms.bool( False )
)


#construct the seeds for the L1 bits needed
#first create a general process, one for tech and one for algo
#then clone them as needed
process.l1tTT = cms.EDFilter( "HLTLevel1GTSeed",
    L1UseL1TriggerObjectMaps = cms.bool( True ),
    L1NrBxInEvent = cms.int32( 3 ),
    L1TechTriggerSeeding = cms.bool( True ),
    L1UseAliasesForSeeding = cms.bool( True ),
    L1SeedsLogicalExpression = cms.string( "0" ),
    L1GtReadoutRecordTag = cms.InputTag( "hltGtDigis" ),
    L1GtObjectMapTag = cms.InputTag( "hltL1GtObjectMap" ),
    L1CollectionsTag = cms.InputTag( "hltL1extraParticles" ),
    L1MuonCollectionTag = cms.InputTag( "hltL1extraParticles" )
)

process.l1tALGO = cms.EDFilter( "HLTLevel1GTSeed",
    L1UseL1TriggerObjectMaps = cms.bool( True ),
    L1NrBxInEvent = cms.int32( 3 ),
    L1TechTriggerSeeding = cms.bool( False ),
    L1UseAliasesForSeeding = cms.bool( True ),
    L1SeedsLogicalExpression = cms.string( "L1SingleEG2" ),
    L1GtReadoutRecordTag = cms.InputTag( "hltGtDigis" ),
    L1GtObjectMapTag = cms.InputTag( "hltL1GtObjectMap" ),
    L1CollectionsTag = cms.InputTag( "hltL1extraParticles" ),
    L1MuonCollectionTag = cms.InputTag( "hltL1extraParticles" )
)

#unmask all bits, I don't know why this is necessary, it used to work
#without it in Run2010A

process.L1GtTriggerMaskAlgoTrigTrivialProducer = cms.ESProducer( "L1GtTriggerMaskAlgoTrigTrivialProducer",
  appendToDataLabel = cms.string( "" ),
  TriggerMask = cms.vuint32( 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 )
)
process.L1GtTriggerMaskTechTrigTrivialProducer = cms.ESProducer( "L1GtTriggerMaskTechTrigTrivialProducer",
  appendToDataLabel = cms.string( "" ),
  TriggerMask = cms.vuint32( 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 )
)

"""
        )
   
    #----------clone the L1 seeders; first the tt 4
    cfg.write("\n#clone the L1 seeders as needed for the L1 bits of interest")
    cfg.write(
"""
process.l1tTT4 = process.l1tTT.clone(
    L1SeedsLogicalExpression = cms.string('4')
    )
"""
        )
    #---------clone the algo seeders
    for seed in algolist:
        #modules can't have "_"
        theseed = ""
        tmplist = seed.split("L1_")[1].split("_")
        for name in tmplist:
            theseed += name
        cfg.write("process.l1t%s = process.l1tALGO.clone(\n" % theseed)
        cfg.write("L1SeedsLogicalExpression = cms.string('%s')\n)\n" % seed)
    #---------clone the tt seeders
    for tseed in techlist:
        cfg.write("process.l1tTT%s = process.l1tTT.clone(\n" % tseed)
        cfg.write("L1SeedsLogicalExpression = cms.string('%s')\n)\n" % tseed)

    #---------smart prescaler to skim according to the L1 rate in a given run
    cfg.write(
"""
#smart prescaler to skim according to the L1 rate in a given run
process.hltPreHLTIMINGSKIMSmart = cms.EDFilter( "TriggerResultsFilter",
\ttriggerConditions = cms.vstring(
"""
    )
    cfg.write("\t\t'HLT_TTT4 / %i',\n" % prezb)
    for seed in algolist:
        cfg.write("\t\t'HLT_T%s',\n" % seed.split("L1_")[1].rstrip())
    for tseed in techlist:
        cfg.write("\t\t'HLT_TTT%s',\n" % tseed)
    cfg.write(
"""
 ),
    \thltResults = cms.InputTag( "TriggerResults" ),
    \tl1tResults = cms.InputTag( "" ),
    \tl1tIgnoreMask = cms.bool( False ),
    \tl1techIgnorePrescales = cms.bool(False),
    \tdaqPartitions = cms.uint32( 1 ),
    \tthrow = cms.bool( True )
    )
"""
        )

    #--------Write the output file
    cfg.write("process.hltOutputHLTIMINGSKIM = cms.OutputModule('PoolOutputModule',\n")
    cfg.write("\t\t\t\tfileName = cms.untracked.string('rfio:%s')," % out) 
    cfg.write(
"""
                                     outputCommands = cms.untracked.vstring(
                                                     'drop *','keep FEDRawDataCollection_source_*_*'
                                                     ),
                                     SelectEvents = cms.untracked.PSet(
                                                     SelectEvents = cms.vstring(
                                                     'HLT_TTT4',
"""
        )
    for seed in algolist:
        cfg.write("\t\t'HLT_T%s',\n" % seed.split("L1_")[1].rstrip())
    for tseed in techlist:
        cfg.write("\t\t'HLT_TTT%s',\n" % tseed)
    cfg.write("\t\t)\n\t),\n)")

    #---------sequences
    cfg.write(
"""
#Sequences
process.HLTL1UnpackerSequence = cms.Sequence( process.hltGtDigis + process.hltGctDigis + process.hltL1GtObjectMap + process.hltL1extraParticles )
"""
        )
    #---------paths
    cfg.write(
"""
\n#Paths (TTT4 always the same)
process.HLT_TTT4 = cms.Path(process.HLTL1UnpackerSequence + process.l1tTT4)
"""
        )
    for seed in algolist:
        #modules can't have "_"
        theseed = ""
        tmplist = seed.split("L1_")[1].split("_")
        for name in tmplist:
            theseed += name
        cfg.write("process.HLT_T%s = cms.Path(process.HLTL1UnpackerSequence + process.l1tTT4 + process.l1t%s)\n" % (seed.split("L1_")[1].rstrip(),theseed))
    for tseed in techlist:
        cfg.write("process.HLT_TTT%s = cms.Path(process.HLTL1UnpackerSequence + process.l1tTT4 + process.l1tTT%s)\n" % (tseed,tseed))

    #--------endpaths, messagelogger config, etc
    cfg.write(
"""
#EndPaths
process.HLTIMINGSKIMOutput = cms.EndPath(process.hltPreHLTIMINGSKIMSmart + process.hltOutputHLTIMINGSKIM)

#configure message logger
process.MessageLogger.cerr.FwkReport.reportEvery = 1 
process.MessageLogger.cerr_stats.threshold = 'INFO' 
#activate summary
process.options = cms.untracked.PSet(
    wantSummary = cms.untracked.bool( True )
)

#Set maximum number of events to process
process.maxEvents = cms.untracked.PSet(
   input = cms.untracked.int32(-1)
)


#append reports
if 'hltTrigReport' in process.__dict__:
    process.hltTrigReport.HLTriggerResults       = cms.InputTag( 'TriggerResults','',process.name_() )
if 'hltPreHLTIMINGSKIMSmart' in process.__dict__:
    process.hltPreHLTIMINGSKIMSmart.TriggerResultsTag = cms.InputTag( 'TriggerResults','',process.name_() )

process.MessageLogger.categories.append('L1GtTrigReport')
process.MessageLogger.categories.append('HLTrigReport')
"""
        )

    #------------all the files in dataset
    #get all the MinBias RAW files for this run
    str_nsls = "nsls "+path
    mypipe = subprocess.Popen(str_nsls,shell=True,stdout=subprocess.PIPE)
    files = mypipe.communicate()[0].split()
    #print files
    cfg.write("\nprocess.source.fileNames = (\n")
    for k in files:
        kfile = "\t'rfio:"+path+"/"+k+"',\n"
        cfg.write(kfile)
    cfg.write(")\n")
    cfg.write(
"""
#in case no ordering is needed
#process.source.noEventSort = cms.untracked.bool(False)
#to select only certain lumi sections
#process.source.lumisToProcess = cms.untracked. VLuminosityBlockRange('132601:10-132601:100')\n
#in case you need to limit the file size to ~4.0Gb to avoid problems
#Several files will be created if this limit is surpassed
#The value is in Kb.
process.hltOutputHLTIMINGSKIM.maxSize = cms.untracked.int32(4000000)
"""
        )
    
    cfg.close()

    print "\nfile skimL1Seeds.py has been created"



#######################################################
def print_batch_script():
#######################################################
    #clean beforehand
    os.system("rm -f skimL1Seeds.csh")
    #get the current working directory
    wdir = os.getcwd()
    #write batch script file
    fbatch = open("skimL1Seeds.csh","a")
    fbatch.write("#!/bin/csh\n")
    fbatch.write("cd %s" % wdir)
    fbatch.write(
"""
cmsenv
cmsRun skimL1Seeds.py
"""
        )
    fbatch.close()
    print "batch file skimL1Seeds.csh has been created"
    #hint on how to submit the job
    print "Now you can sumbit a batch job as:"
    print "bsub -q 8nh -J skimL1ZB < skimL1Seeds.csh"
    

#######################################################
def get_default_options(option):
#######################################################
    dicOpt = {}

    dicOpt['path']= str(option.path)
    dicOpt['out']= str(option.out)
    dicOpt['zbrate']= option.zbrate
    dicOpt['trate']= option.trate


    return dicOpt





#######################################################
if __name__ =='__main__':
#######################################################

    #import optionparser
    option,args = parse(__doc__)
    if not args and not option:
        exit()
    
    
    #set default options
    dicOpt = get_default_options(option)

    #print configuration
    for k in dicOpt:
        print str(k)+" = "+str(dicOpt[k])
    
    #get unmasked L1 bit list
    algolist,techlist = get_L1bits_list(dicOpt)

    #print configuration file
    print_config_file(dicOpt,algolist,techlist)

    #prepare script for batch submission
    print_batch_script()
