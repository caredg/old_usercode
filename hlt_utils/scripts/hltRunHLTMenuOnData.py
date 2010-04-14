#!/usr/bin/env python
############################################################################
#
# Author: Edgar Carrera
# ecarrera@cern.ch
#
# Based on various scripts by:
#    Bryan Dhames,
#    Connor Henderson, and
#    Maurizio Pierini.
#        
# 
# Oct 5, 2009
# Severely modified Dec 11, 2009: Now use getHLT.py from the release area.
# 
############################################################################

"""
   usage: %prog [options]
   -k, --hltkey = HLTKEY: HLT menu configuration. Ex: 'orcoff:/cdaq/cosmic/commissioning2009/CRAFT/v2.4/HLT/V2' or '/online/blah/blah/HLT/V1'
   -n, --events = EVENTS: number of events to be run. Default is -1 (run over all events)
   -l, --list   = LIST: file with list of files (dataset).  If not specified a dummy file is inserted. (NOT IMPLEMENTED)
   -L, --l1over  =  L1OVER: Enable L1 menu override, using the given payload from the database. Ex L1_mymenu_v5
   -g, --gtag = GLOBAL : Global tag. e.g., GR_P1_V1 (leave the ::All part out).  If not specified, the tag in the implementation is taken.
   -i, --id = ID : Id in file name.  Default is GRunData.
   -t, --timing : Switches on the timing module. Default is off.
   -r, --run = RUNNUMBER : Run number
   -s, --singlerun : Use --singleRun in error_stream_collector.pl script if errstream studies is activated.
   -e, --era = ERA : This activates errstream studies. Era for the error_stream_collector.pl script. Ex: Data
   -m, --dontmess : Do not make any additional modification on top of what the getHLT.py script already does
"""


import os,sys
import string, re
import fileinput
import commands
from time import gmtime, localtime, strftime

#flag to debug
DEBUG = False
printConfig = True

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
def run_errstream_collector(singleRun, era, irun):
#######################################################
    if (DEBUG):
        print "Runing run_errstream_collector"
        
   
    if os.path.exists("/afs/cern.ch/user/c/chenders/public/hlt/error_stream_collector/error_stream_collector.pl"):
        os.system("cp /afs/cern.ch/user/c/chenders/public/hlt/error_stream_collector/error_stream_collector.pl .")
    else:
        print "No error_stream_collector.pl script could be found at /afs/cern.ch/user/c/chenders/public/hlt/error_stream_collector/"
        exit(1)

    os.system("./error_stream_collector.pl "+singleRun+" --era="+era+" --doCopy "+irun)
    #os.system("./error_stream_collector.pl "+singleRun+" --era= "+era+" "+irun)
    outdir = os.popen("echo /tmp/$USER").readline().rstrip()


    return outdir
             
    
    


#######################################################
def run_on_errorstream(configfile, dicOpt):
#######################################################

    if (DEBUG):
        print "running run_on_errorstream"
        
    singleRun = dicOpt['singlerun']
    era = dicOpt['era']
    irun = dicOpt['run']
    #run the script error_stream_collector.pl and return
    #the location of the copied files
    errDir = run_errstream_collector(singleRun,era,irun)

    #From Maurizio:
    inputdir = errDir
    print "inputdir = "+inputdir
    # check the input directory
    if not os.path.isdir(inputdir):
        print "ERROR: input directory does not exist"
        sys.exit(1)
    # set the output directory
    outputdir = "errstream"
    if os.path.isdir(outputdir) == 0: os.system("mkdir "+outputdir)
    else: print "WARNING: Directory "+outputdir+" already exists."
    # set the subdirectories:
    #1) root
    if os.path.isdir(outputdir+"/root") == 0: os.system("mkdir "+outputdir+"/root")
    else: print "WARNING: Directory "+outputdir+"/root already exists."
    #2) cfg
    if os.path.isdir(outputdir+"/cfg") == 0: os.system("mkdir "+outputdir+"/cfg")
    else: print "WARNING: Directory "+outputdir+"/cfg already exists."
    #3) log
    if os.path.isdir(outputdir+"/log") == 0: os.system("mkdir "+outputdir+"/log")
    else: print "WARNING: Directory "+outputdir+"/log already exists."
    # create input files list
    os.system("ls "+inputdir+"/*.dat > myFilesList.tmp")
    list = open("myFilesList.tmp")
    #IF REQUESTED convert the HLT configuration ONLINE->OFFLINE
    #HLT = "none"
    #for iIN in range(3, len(sys.argv)):
    #if sys.argv[iIN].startswith('-HLT='): HLT = sys.argv[iIN].split('-HLT=')[1]
    #continue
#if HLT != "none": os.system("/afs/cern.ch/user/m/mpierini/public/TMI/ErrorStream/getHLT.py "+HLT+" GRunData")
# loop over the input files and produce the output file
    print "\nLooping over files.... be patient" 
    for line in list:
        filename = line.split("/")[3].split(".dat")[0]
        print "Working on file = "+filename
        #copy the tenmplate file
        os.system("cp /afs/cern.ch/user/m/mpierini/public/TMI/ErrorStream/converterError_cfg.py template_cfg.py")
        #create the conversion _py.cfg file 
        template  = open("template_cfg.py")
        mycfgfile = open(outputdir+"/cfg/"+filename+"_cfg.py","w")
        for templateline in template: mycfgfile.write(templateline.replace("file:in", "file:"+inputdir+"/"+filename).replace("file:out", "file:"+outputdir+"/root/"+filename))
        template.close()
        mycfgfile.close()
         #run the .dat->.root conversion _py.cfg file
        os.system("cmsRun "+outputdir+"/cfg/"+filename+"_cfg.py >&"+outputdir+"/log/convert_"+filename+".log")
        
        # run the requested HLT menu on the ErrorStream .root file
        # and log the output 
        template2 = open(configfile)
        mycfgfile2 = open(outputdir+"/cfg/"+filename+"_HLT_cfg.py","w")
        for templateline in template2:
            if templateline.find(".root") != -1:
                mycfgfile2.write(templateline.replace("/tmp/InputCollection","").replace("file:","file:"+outputdir+"/root/").replace(".root",filename+".root"))
            else:
                mycfgfile2.write(templateline.replace("process.hltTriggerType +", " "))
                continue
            continue
        template2.close()
        mycfgfile2.close()
        os.system("cmsRun "+outputdir+"/cfg/"+filename+"_HLT_cfg.py >& "+outputdir+"/log/convert_"+filename+"_HLT.log")
        os.system("mv output*"+filename+"* "+outputdir+"/root/.")
        continue
    list.close()
    if not (DEBUG):
        os.system("rm myFilesList.tmp")
        os.system("rm template_cfg.py")

    
    return 0



#######################################################
def modify_config_file(configfile,dicOpt):
#######################################################

    outName = configfile
    numevents = dicOpt['events']
    poolout = dicOpt['filename']
    timing = dicOpt['timing']
    myID = dicOpt['id']
    
    #Modify Message logger to report every 50 events and WARN
    mlogger_str = "sed -i.original -e '0, /reportEvery = cms.untracked.int32( 1 )/s//reportEvery = cms.untracked.int32( 50 )/' -e '0, /threshold = cms.untracked.string( \"ERROR\" )/s//threshold = cms.untracked.string( \"WARNING\" )/' -e 's#\(^ *limit = cms.untracked.int32( \)0 )#\\1100000 )#' "+outName
    os.system(mlogger_str)

    #Modify Event Number content if requested
    if numevents:
        nevents_str = "sed -i.nume 's#\(^ *input = cms.untracked.int32(\) .* )#\\1 "+numevents+ " )#' "+outName
        os.system(nevents_str)

    #add output file and check if timing is needed    
    myfile = open(outName,'a')

    myfile.write("""
#Output
process.output = cms.OutputModule("PoolOutputModule",
    splitLevel = cms.untracked.int32(0),
    outputCommands = cms.untracked.vstring( 'drop *', 'keep HLTPerformanceInfo_*_*_*'),
    fileName = cms.untracked.string('HLT.root'),
    dataset = cms.untracked.PSet(
        dataTier = cms.untracked.string('RECO'),
        filterName = cms.untracked.string('')
    )
)
\n""")
    
    #add timing module if requested
    if timing:
        myfile.write("""
#timer
process.PathTimerService = cms.Service( "PathTimerService" )
process.hltTimer = cms.EDProducer( "PathTimerInserter" )
process.HLTOutput_edgar = cms.EndPath( process.hltTimer + process.output)
\n""")
           

    myfile.close()
    
  
    return 0



#######################################################
def get_basic_config(dicOpt):
#######################################################

    dbName = dicOpt['hltkey']
    myID = dicOpt['id']
    myl1ov = dicOpt['l1over']
    mygtag = dicOpt['gtag']
    outName = "OnLine_HLT_"+myID+".py"

    #check if require global tag
    if (mygtag):
        getgtag = "-- globaltag "+mygtag+" "
    else:
        getgtag = ""


    getHLT = "$CMSSW_RELEASE_BASE/src/HLTrigger/Configuration/test/getHLT.py"
    #The options might change from release to release:
    dogetHLT = getHLT+" --data --force "+getgtag+dbName+" "+myID
    os.system(dogetHLT)
    #print dogetHLT


    return outName


#######################################################
def get_default_options(option):
#######################################################
    dicOpt = {}

    #hlt menu
    dicOpt['hltkey']= str(option.hltkey)

    # number of events
    if not option.events:
        dicOpt['events'] = "-1"
    else:
        dicOpt['events'] = str(option.events)
    
    # list of files    
    if not option.list:
            dicOpt['filename'] = "dummy_input.root"
    else:
        dicOpt['filename'] = "dummy_input.root"

    # override L1 menu?
    if not option.l1over:
            dicOpt['l1over'] = None
    else:
        dicOpt['l1over'] = "--l1"

    # global tag
    if not option.gtag:
        dicOpt['gtag'] = None
    else:
        dicOpt['gtag'] = str(option.gtag)

    #Id in file name
    if not option.id:
        dicOpt['id'] = "GRunData"
    else:
        dicOpt['id'] = str(option.id)

    #Timing
    if not option.timing:
        dicOpt['timing'] = False
    else:
        dicOpt['timing'] = option.timing

    #run number
    if not option.run:
        dicOpt['run'] = None
    else:
        dicOpt['run'] = str(option.run)

    #Single run
    if not option.singlerun:
        dicOpt['singlerun'] = ""
    else:
        dicOpt['singlerun'] = "--singleRun"

    #Era
    if not option.era:
        dicOpt['era'] = None
    else:
        dicOpt['era'] = str(option.era)

    #Modification to offline
    if not option.dontmess:
        dicOpt['dontmess'] = False
    else:
        dicOpt['dontmess'] = option.dontmess

     

    return dicOpt




#######################################################
if __name__ =='__main__':
#######################################################

    #import optionparser
    option,args = parse(__doc__)
    if not args and not option:
        exit()

    #safety check for HLT menu config    
    if not option.hltkey:
        print " you need to provide at least the HLT configuration key string"
        optionparse.exit()


    #set default options
    dicOpt = get_default_options(option)

    if (printConfig):
        for k in dicOpt:
            print str(k)+" = "+str(dicOpt[k])
    
    #get the basic configuration file
    configfile = get_basic_config(dicOpt)
    #modify cfg file to run offline if requested
    # exit if no modification is requiered
    if not dicOpt['dontmess']:
        modify_config_file(configfile,dicOpt)
    else:
        print "The config dump has not been modified"
    #run on error stream files if requested
    if dicOpt['era'] and dicOpt['run']:
        run_on_errorstream(configfile,dicOpt)
    else:
        if dicOpt['run']:
            print "NOTE: You need --era to activate running on errorstream"

    
    
    
    
    
