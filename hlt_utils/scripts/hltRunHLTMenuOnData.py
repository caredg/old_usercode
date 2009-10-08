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
#
# 
############################################################################

"""
   usage: %prog [options]
   -k, --hltkey = HLTKEY: HLT menu configuration. Ex: /cdaq/cosmic/commissioning2009/CRAFT/v2.4/HLT/V2
   -n, --events = EVENTS: number of events to be run. Default is -1 (run over all events)
   -l, --list   = LIST: file with list of files (dataset).  If not specified a dummy file is inserted. (NOT IMPLEMENTED)
   -g, --gtag = GLOBAL : Global tag. e.g., GR_P1_V1 (leave the ::All part out).  If not specified, the tag in the implementation is taken. (NOT IMPLEMENTED)
   -d, --devel : Turns off the orcoff switch for extracting HLT menu config. Default is on.
   -u, --usecase = USERCASE  : If GEN-HLT is specified, a stripped file is generated for the GEN-HLT workflow.  Default is ONLINE for minimal modifications.
   -i, --id = ID : Id in file name.  Default is GRunData.
   -t, --timing : Switches on the timing module. Default is off.
   -r, --run = RUNNUMBER : Run number
   -s, --singlerun : Use --singleRun in error_stream_collector.pl script if errstream studies is activated.
   -e, --era = ERA : This activates errstream studies. Era for the error_stream_collector.pl script. Ex: MWGR_40_2009
   -m, --online : Do not make the modification to run offline.  Default is to make the modifications to run offline.
"""


import os,sys
import string, re
import fileinput
import commands
from time import gmtime, localtime, strftime

#flag to debug
DEBUG = True


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
def run_on_errorstream(dicOpt):
#######################################################
    if (DEBUG):
        print "running run_on_errorstream"
        
    if dicOpt['era'] and dicOpt['run']:
        singleRun = dicOpt['singlerun']
        era = dicOpt['era']
        irun = dicOpt['run']
        #run the script error_stream_collector.pl and return
        #the location of the copied files
        errDir = run_errstream_collector(singleRun,era,irun)
    else:
        print "You need --era to activate running on errorstream"
        return 0

    #From Maurizio:
    inputdir = errDir
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
    for line in list:
        filename = line.split("/")[3].split(".dat")[0]
        #copy the tenmplate file
        os.system("cp /afs/cern.ch/user/m/mpierini/public/TMI/ErrorStream/converterError_cfg.py template_cfg.py")
        #create the conversion _py.cfg file 
        template  = open("template_cfg.py")
        mycfgfile = open(outputdir+"/cfg/"+filename+"_cfg.py","w")
        for templateline in template: mycfgfile.write(templateline.replace("file:in", "file:"+inputdir+"/"+filename).replace("file:out", "file:"+outputdir+"/root/"+filename))
        template.close()
        mycfgfile.close()
        # run the .dat->.root conversion _py.cfg file
        os.system("cmsRun "+outputdir+"/cfg/"+filename+"_cfg.py >&"+outputdir+"/log/convert_"+filename+".log")
        # run the requested HLT menu on the ErrorStream .root file
        # and log the output 
        template2 = open("OnLine_HLT_GRunData.py")
        mycfgfile2 = open(outputdir+"/cfg/"+filename+"_HLT_cfg.py","w")
        for templateline in template2:
            if templateline.find(".root") != -1:
                mycfgfile2.write(templateline.replace("RawInput_GRunData","").replace("file:","file:"+outputdir+"/root/").replace(".root",filename+".root"))
            else:
                mycfgfile2.write(templateline.replace("process.hltTriggerType +", " "))
                continue
            continue
        template2.close()
        mycfgfile2.close()
        os.system("cmsRun "+outputdir+"/cfg/"+filename+"_HLT_cfg.py >& "+outputdir+"/log/convert_"+filename+"_HLT.log")
        continue
    list.close()
    #os.system("rm OnLine_HLT_GRunData.py")
    os.system("rm myFilesList.tmp")
    os.system("rm template_cfg.py")
    
    return 0



#######################################################
def modify_config_file(configfile,dicOpt):
#######################################################

    # exit if no modification is requested
    nomodif = dicOpt['online']
    if nomodif:
        return 0
    
    outName = configfile
    numevents = dicOpt['events']
    poolout = dicOpt['filename']
    gtag = dicOpt['gtag']
    useCase = dicOpt['usecase']
    timing = dicOpt['timing']
    myID = dicOpt['id']
    
    #
    # in situ replacements:
    #    - Convert ShmStreamConsumer to PoolOutputModule
    #
    addAnaly = True
    for line in fileinput.input(outName,inplace=1):
        if line.find("HLTAnalyzerEndpath")>= 0:
            addAnaly = False
        if line.find("ShmStreamConsumer") >= 0:
            line = line.replace('ShmStreamConsumer','PoolOutputModule')
            outputName = line[line.find(".")+1:line.find("=")-1]
            print line[:-1]
            print "    fileName = cms.untracked.string('file:HLT_" + outputName + ".root'),"
            print "    basketSize = cms.untracked.int32(4096),"
        else:
            print line[:-1]


    #
    # Overwrite ProcessName, PoolSource and include maxEvents
    #
    os.system("cat >> "+outName+" <<EOI\nprocess.setName_('HLT"+myID+"')\nEOI\n")
    os.system("cat >> "+outName+" <<EOI\nprocess.source = cms.Source(\"PoolSource\",\nEOI\n")
    if myID.find("Data") >= 0:
        os.system("cat >> "+outName+" <<EOI\n    fileNames = cms.untracked.vstring('file:RawInput_"+myID+".root')\nEOI\n")
    else:
        os.system("cat >> "+outName+" <<EOI\n    fileNames = cms.untracked.vstring('file:RelVal_DigiL1Raw_"+myID+".root')\nEOI\n")
    os.system("cat >> "+outName+" <<EOI\n)\nEOI\n")

    os.system("cat >> "+outName+" <<EOI\nprocess.maxEvents = cms.untracked.PSet(\nEOI\n")
    os.system("cat >> "+outName+" <<EOI\n    input = cms.untracked.int32("+numevents+")\nEOI\n")
    os.system("cat >> "+outName+" <<EOI\n)\nEOI\n")

    #
    # Overwrite GlobalTag
    #
    print "WARNING: Global Tags are harcoded for now check if they are correct!!!"
    os.system("cat >> "+outName+" <<EOI\nprocess.GlobalTag.connect = 'frontier://FrontierProd/CMS_COND_31X_GLOBALTAG'\nEOI\n")

    if myID.find("Data") >= 0:
        os.system("cat >> "+outName+" <<EOI\n# Use this for running on CRAFT 09 data\nEOI\n")
        os.system("cat >> "+outName+" <<EOI\nprocess.GlobalTag.globaltag = 'GR09_P_V3::All'\nEOI\n")
    else:
       
        if myID.find("8E29") >= 0:
            os.system("cat >> "+outName+" <<EOI\nprocess.GlobalTag.globaltag = 'STARTUP31X_V2::All'\nEOI\n")
        elif myID.find("GRun") >= 0:            
            os.system("cat >> "+outName+" <<EOI\nprocess.GlobalTag.globaltag = 'STARTUP31X_V2::All'\nEOI\n")
        elif myID.find("1E31") >= 0:
            os.system("cat >> "+outName+" <<EOI\nprocess.GlobalTag.globaltag = 'MC_31X_V1::All'\nEOI\n")
        elif myID.find("HIon") >= 0:
            os.system("cat >> "+outName+" <<EOI\nprocess.GlobalTag.globaltag = 'MC_31X_V1::All'\nEOI\n")
        else:
            os.system("cat >> "+outName+" <<EOI\nprocess.GlobalTag.globaltag = 'MC_31X_V1::All'\nEOI\n")


    # 
    # Running on Monte Carlo: change RAW input from "source" to "rawDataCollector"
    #
    prefix = "process."
    if useCase == "GEN-HLT":
        prefix = ""

    #
    # Special case: Running on 0T data
    #
    if myID.find("GRun") >= 0:
        os.system("cat >> "+outName+" <<EOI\n\n# Uncomment these lines to run on 0T data\nEOI\n")
        os.system("cat >> "+outName+" <<EOI\n# process.GlobalTag.globaltag = 'GR09_P_V3::All'\nEOI\n")
        os.system("cat >> "+outName+" <<EOI\n# process.localUniform = cms.ESProducer( \"UniformMagneticFieldESProducer\",\nEOI\n")
        os.system("cat >> "+outName+" <<EOI\n#   ZFieldInTesla = cms.double( 0.0 ),\nEOI\n")
        os.system("cat >> "+outName+" <<EOI\n#   appendToDataLabel = cms.string( \"\" )\nEOI\n")
        os.system("cat >> "+outName+" <<EOI\n# )\nEOI\n")
        os.system("cat >> "+outName+" <<EOI\n# process.es_prefer_localUniform = cms.ESPrefer( \"UniformMagneticFieldESProducer\", \"localUniform\" )\nEOI\n")
        os.system("cat >> "+outName+" <<EOI\n# "+prefix+"FastSteppingHelixPropagatorAny.SetVBFPointer = True\nEOI\n")
        os.system("cat >> "+outName+" <<EOI\n# "+prefix+"FastSteppingHelixPropagatorOpposite.SetVBFPointer = True\nEOI\n")
        os.system("cat >> "+outName+" <<EOI\n# "+prefix+"SteppingHelixPropagatorAlong.SetVBFPointer = True\nEOI\n")
        os.system("cat >> "+outName+" <<EOI\n# "+prefix+"SteppingHelixPropagatorAny.SetVBFPointer = True\nEOI\n")
        os.system("cat >> "+outName+" <<EOI\n# "+prefix+"SteppingHelixPropagatorOpposite.SetVBFPointer = True\nEOI\n")
        os.system("cat >> "+outName+" <<EOI\n# "+prefix+"VolumeBasedMagneticFieldESProducer.label = 'VolumeBasedMagneticField'\nEOI\n")
        os.system("cat >> "+outName+" <<EOI\n# End of 0T changes\nEOI\n")
        
        
        
    if myID.find("MC") >= 0:
        os.system("cat >> "+outName+" <<EOI\n\n# Data vs. Monte Carlo changes\nEOI\n")
        os.system("cat >> "+outName+" <<EOI\n"+prefix+"hltGtDigis.DaqGtInputTag = 'rawDataCollector'\nEOI\n")
        os.system("cat >> "+outName+" <<EOI\n"+prefix+"hltGctDigis.inputLabel = 'rawDataCollector'\nEOI\n")
        os.system("cat >> "+outName+" <<EOI\n"+prefix+"hltSiPixelDigis.InputLabel = 'rawDataCollector'\nEOI\n")
        os.system("cat >> "+outName+" <<EOI\n"+prefix+"hltSiStripRawToClustersFacility.ProductLabel = 'rawDataCollector'\nEOI\n")
        os.system("cat >> "+outName+" <<EOI\n"+prefix+"hltEcalRawToRecHitFacility.sourceTag = 'rawDataCollector'\nEOI\n")
        os.system("cat >> "+outName+" <<EOI\n"+prefix+"hltESRawToRecHitFacility.sourceTag = 'rawDataCollector'\nEOI\n")
        os.system("cat >> "+outName+" <<EOI\n"+prefix+"hltHcalDigis.InputLabel = 'rawDataCollector'\nEOI\n")
        os.system("cat >> "+outName+" <<EOI\n"+prefix+"hltMuonCSCDigis.InputObjects = 'rawDataCollector'\nEOI\n")
        os.system("cat >> "+outName+" <<EOI\n"+prefix+"hltMuonDTDigis.inputLabel = cms.untracked.InputTag( 'rawDataCollector' )\nEOI\n")
        os.system("cat >> "+outName+" <<EOI\n"+prefix+"hltMuonRPCDigis.InputLabel = 'rawDataCollector'\nEOI\n")
        if myID.find("GRun") >= 0:
            os.system("cat >> "+outName+" <<EOI\n"+prefix+"hltPixelFEDSizeFilter.rawData = 'rawDataCollector'\nEOI\n")
            os.system("cat >> "+outName+" <<EOI\n"+prefix+"hltEcalCalibrationRaw.rawInputLabel = 'rawDataCollector'\nEOI\n")
        os.system("cat >> "+outName+" <<EOI\n"+prefix+"hltAlCaHcalFEDSelector.rawInputLabel = 'rawDataCollector'\nEOI\n")
        if useCase != "GEN-HLT":
            os.system("cat >> "+outName+" <<EOI\n"+prefix+"DQML1Scalers.fedRawData = 'rawDataCollector'\nEOI\n")
            os.system("cat >> "+outName+" <<EOI\n"+prefix+"DQMHLTScalers.triggerResults = cms.InputTag( 'TriggerResults','','HLT"+myID+"' )\nEOI\n")
        
    #
    # Include the HLTAnalyzerEndpath
    #
    if myID.find("GRun") >= 0  and addAnaly:
        os.system("cat >> "+outName+" <<EOI\n\n# Attach HLTAnalyzerEndpath\nEOI\n")
        os.system("cat >> "+outName+" <<EOI\n"+prefix+"hltL1GtTrigReport = cms.EDAnalyzer( \"L1GtTrigReport\",\nEOI\n")
        os.system("cat >> "+outName+" <<EOI\n    UseL1GlobalTriggerRecord = cms.bool( False ),\nEOI\n")
        os.system("cat >> "+outName+" <<EOI\n    L1GtRecordInputTag = cms.InputTag( \"hltGtDigis\" )\nEOI\n")
        os.system("cat >> "+outName+" <<EOI\n)\nEOI\n")
        os.system("cat >> "+outName+" <<EOI\n"+prefix+"hltTrigReport = cms.EDAnalyzer( \"HLTrigReport\",\nEOI\n")
        os.system("cat >> "+outName+" <<EOI\n    HLTriggerResults = cms.InputTag( 'TriggerResults','','HLT"+myID+"' )\nEOI\n")
        os.system("cat >> "+outName+" <<EOI\n)\nEOI\n")
    
        os.system("cat >> "+outName+" <<EOI\n"+prefix+"HLTAnalyzerEndpath = cms.EndPath( "+prefix+"hltL1GtTrigReport + "+prefix+"hltTrigReport )\nEOI\n")
        if useCase == "GEN-HLT":
            os.system("cat >> "+outName+" <<EOI\n"+prefix+"HLTSchedule.append( HLTAnalyzerEndpath ) \nEOI\n")
    

    #
    # The following is stolen from cmsDriver's ConfigBuilder.py - addCustomise
    #
    if useCase != "GEN-HLT":        
        # let python search for that package and do syntax checking at the same time
        packageName = 'HLTrigger/Configuration/customL1THLT_Options.py'.replace(".py","").replace(".","/")
        package = __import__(packageName)
        
        # now ask the package for its definition and pick .py instead of .pyc
        customiseFile = package.__file__.rstrip("c")
        
        
        final_snippet='\n\n# Automatic addition of the customisation function\n'
        for line in file(customiseFile,'r'):
            if "import FWCore.ParameterSet.Config" in line:
                continue
            final_snippet += line
            
        final_snippet += '\n\n# End of customisation function definition'
        final_snippet += "\n\nprocess = customise(process)\n"
            
        os.system("cat >> "+outName+" <<EOI\n"+final_snippet+"EOI\n")

  
    return 0



#######################################################
def get_basic_config(dicOpt):
#######################################################

    dbName = dicOpt['hltkey']
    myorcoff = dicOpt['devel']
    useCase = dicOpt['usecase']
    myID = dicOpt['id']
    
    
    if useCase == "ONLINE":
        outName = "OnLine_HLT_"+myID+".py"
    else:    
        outName = "HLT_"+myID+"_cff.py"

    if os.path.exists(outName):
        print outName, "already exists - abort!"
        exit()
    else:
        # Initialize everything
        essources = "  " 
        esmodules = "  "
        modules   = "  "
        services  = "  "
        paths     = "  "
        psets     = "  "
        edsources = "  "

        # Required for online-compliant menus
        edsources = "--noedsources"
        
        if useCase == "GEN-HLT":
            essources = "--essources "
            essources += "-SiStripQualityFakeESSource,"
            essources += "-GlobalTag,"
            essources += "-HepPDTESSource,"
            essources += "-XMLIdealGeometryESSource,"
            essources += "-eegeom,"
            essources += "-es_hardcode,"
            essources += "-magfield"
            
            esmodules = "--esmodules "
            esmodules += "-CSCGeometryESModule,"
            esmodules += "-CaloGeometryBuilder,"
            esmodules += "-CaloTowerHardcodeGeometryEP,"
            esmodules += "-DTGeometryESModule,"
            esmodules += "-EcalBarrelGeometryEP,"
            esmodules += "-EcalElectronicsMappingBuilder,"
            esmodules += "-EcalEndcapGeometryEP,"
            esmodules += "-EcalLaserCorrectionService,"    
            esmodules += "-EcalPreshowerGeometryEP,"
            esmodules += "-HcalHardcodeGeometryEP,"
            esmodules += "-HcalTopologyIdealEP,"
            esmodules += "-MuonNumberingInitialization,"
            esmodules += "-RPCGeometryESModule,"
            esmodules += "-SiStripGainESProducer,"
            esmodules += "-SiStripRecHitMatcherESProducer,"
            esmodules += "-StripCPEfromTrackAngleESProducer,"
            esmodules += "-TrackerDigiGeometryESModule,"
            esmodules += "-TrackerGeometricDetESModule,"
            esmodules += "-VolumeBasedMagneticFieldESProducer,"    
            esmodules += "-ZdcHardcodeGeometryEP,"
            esmodules += "-hcal_db_producer,"
            esmodules += "-l1GtTriggerMenuXml,"
            esmodules += "-sistripconn"
            
            services  = "--services -MessageLogger"
            
            paths     = "--paths -HLTOutput,-AlCaOutput,-ESOutput,-MONOutput"
            
            psets     = "--psets -maxEvents,-options"
            

            myGet = "edmConfigFromDB --cff --configName " + dbName + " " + essources + " " + esmodules + " " + modules + " " + services + " " + paths + " " + psets + " > " + outName
            os.system(myGet)
            

            # Inline replacements: Remove prescalers for removed endpaths
            inPrescaleService = False
            inPreServESOutput = False
            inPreServMONOutput = False
            preCounter = -1 
            for line in fileinput.input(outName,inplace=1):
                if inPrescaleService and line.find("pathName") >= 0:
                    inPreServESOutput = False
                    inPreServMONOutput = False
                    if line.find("ESOutput") >= 0:
                        inPreServESOutput = True
                        preCounter = 0
                    elif line.find("MONOutput") >= 0:
                        inPreServMONOutput = True
                        preCounter = 0
                if line.find("PrescaleService") >= 0:
                    inPrescaleService = True
                    print line[:-1]
                elif inPrescaleService:
                    if inPreServESOutput or inPreServMONOutput:
                        preCounter += 1
                        if preCounter > 3:
                            inPreServESOutput = False
                            inPreServMONOutput = False
                            preCounter = -1
                            print line[:-1]
                        else:
                            print line[:-1]
                            if line.find(")") == 0:
                                inPrescaleService = False
                else:
                    print line[:-1]
            
        else:
            myGet = "edmConfigFromDB "+ myorcoff +" --configName " + dbName + " " + edsources + " " + essources + " " + esmodules + " " + modules + " " + services + " " + paths + " " + psets + " > " + outName
            os.system(myGet)


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

    # global tag
    if not option.gtag:
        dicOpt['gtag'] = None
    else:
        dicOpt['gtag'] = str(option.gtag)

    # switch for orcoff
    if not option.devel:
        dicOpt['devel'] = "--orcoff"
    else:
        dicOpt['devel'] = ""

    # usecase option
    if not option.usecase:
        dicOpt['usecase'] = "ONLINE"
    else:
        dicOpt['usecase'] = str(option.usecase)

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
    if not option.online:
        dicOpt['online'] = False
    else:
        dicOpt['online'] = option.online

     

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
    if (DEBUG):
        for k in dicOpt:
            print str(k)+" = "+str(dicOpt[k])
    
    #get the basic configuration file
    configfile = get_basic_config(dicOpt)
    #modify cfg file to run offline is requested
    modify_config_file(configfile,dicOpt)
    #run on error stream files if requested
    run_on_errorstream(dicOpt)
    
    
    
