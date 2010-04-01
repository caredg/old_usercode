#!/usr/bin/python

##############################################################
# Creator: Edgar Carrera
# Date: Sept 25, 2009
# Script to analyze events within a region
# of interest in the plot of Total module time
# per event, which is part of
# the output of the hltTimeSummary macro.
#
#
#
# 
#############################################################
import os, string, re,sys
from time import gmtime, localtime, strftime
#if not built-in in the current python version, try to import Set:
try:
    set()
except NameError:
    from sets import Set as set



#just a debug flag
DEBUG = False



###############################################################
def usage():
###############################################################

    if (DEBUG):
        print "This is the usage function"
        
    print '\n'
    print 'Usage: '+sys.argv[0]+' <tmin in msec> <tmax in msec> <optional HLT cfg file>'
    print 'e.g.:  '+sys.argv[0]+' 40 50 runHLT_cfg.py'
    print 'NOTE: 0<tmin<tmax and both need not be integers but 0 is not allowed (fixme later)'
    print 'the optional file is any config py file which is going to be'
    print 'expanded to just look at the events of interest'





###############################################################
def extend_HLTcfg_file(hltcfg_file,peakEvtsFile):
###############################################################
    
    if (DEBUG):
        print "Extending config file "+hltcfg_file

    #extend the cfg file
    fin = open(hltcfg_file,"a")
    fin.write("\ntemplist = open(\""+peakEvtsFile+"\",\"r\").readlines()\n")
    fin.write("mylist = []\n")
    fin.write("for k in templist:\n")
    fin.write("\tmylist.append(k.rstrip())\n")
    fin.write("\n")
    fin.write("process.source.eventsToProcess = cms.untracked.VEventRange(mylist)\n")
    fin.close()




###############################################################
def write_dat_file(peakEvtsFile,setPeakEvts):
###############################################################

    if (DEBUG):
        print "Creating *.dat file with peak events"

    fout = open(peakEvtsFile,"w")
    for line in setPeakEvts:
        fout.write(line+"\n")

    fout.close()
    return 0




###############################################################
def parse_hltTimingOutput(inputfile):
###############################################################

    if (DEBUG):
        print "Parsing hltTimingSummary "+ inputfile +" file."

    #nicely parse the summary file using lines that
    #start with 'Run:' and create a set with format Run:EvtNum
    fin = open(inputfile)
    evtList = []
    for line in fin.readlines():
        aux = line.find('Run:')
        if aux == 0:
            evtList.append(line.split()[1].split(",")[0]+":"+line.split()[3])

    if (len(evtList) == 0):
        print "No events were found in the "+inputfile+" file. Are you sure you are using the latest version of hltTimingSummary macro?"

    evtSet = set(evtList)

    if (DEBUG):
        print len(evtSet)

    return evtSet




    
###############################################################
def run_hltTimingSummary(tmin, outfile):
###############################################################

    print "Running hltTimingSummary macro for file "+outfile
    print "Be patient ...."

    #Hardwire the input file and some parameters.
    inrootfile = 'hlt.root' #default input in the hltTimingSummary script
    exclusionfile = 'exclude.txt'
    totaltimehist_maxX = str(200)
    totaltimehist_bin = str(2)
    pdfverbosity = str(0)

    #check if input file exists
    if  not(os.path.isfile(inrootfile)):
        print inrootfile+" does not exist. Please check."
        sys.exit(1)

    #check if the file that contains the hltGetRaw module for exclusion
    #exists and issue a warning if not
    if not(os.path.isfile(exclusionfile)):
        print "WARNING! exclude.txt does not exist"
        
    hltTS_command = "$CMSSW_BASE/test/$SCRAM_ARCH/hltTimingSummary -i "+ inrootfile +" -o "+ outfile + " -f -c -s -e "+ exclusionfile + " -l " + tmin + " -p "+ pdfverbosity
              

    if (DEBUG):
        print hltTS_command

    os.system(hltTS_command)

    return 0






###############################################################
def get_eventNumbers(tmin,tmax,hltcfg_file):
###############################################################

    if (DEBUG):
        print "Getting event numbers in the peak"

    #give specific names to outputfiles    
    outfile1 = 'peakoutfile_tt'+ tmin
    outfile2 = 'peakoutfile_tt'+ tmax
    outsummary1 = outfile1+"-summary.txt"
    outsummary2 = outfile2+"-summary.txt"
    #give name to name of final dat file with all
    #the events of interest
    peakEvtsFile = "peakevents_tt_"+tmin+"_"+tmax+".dat"

    # run hltTimingSummary macro
    if os.path.isfile(outsummary1):
        print "File "+outsummary1+" exists. Using the already created file..."
    else:
        run_hltTimingSummary(tmin, outfile1)

    if os.path.isfile(outsummary2):
        print "File "+outsummary2+" exists. Using the already created file..."
    else:
        run_hltTimingSummary(tmax, outfile2)

    #make sure the outputfiles were created
    if  not(os.path.isfile(outsummary1) and os.path.isfile(outsummary2)):
        print "File "+outsummary1+" and/or "+outsummary2+" do not exist."
        sys.exit(1)

    #parse the output files, get set of events
    setEvents1 = parse_hltTimingOutput(outsummary1)
    setEvents2 = parse_hltTimingOutput(outsummary2)

    #look at the events in set1 but not in set2
    #which defines the region we are interested in
    setPeakEvts = setEvents1.difference(setEvents2)

    #print number of events
    print "Num of events from file "+outsummary1+": "+str(len(setEvents1))
    print "Num of events from file "+outsummary2+": "+str(len(setEvents2))
    print "Num of events in the peak ["+tmin+":"+tmax+"]: "+str(len(setPeakEvts))

    #write events to a *.dat file
    write_dat_file(peakEvtsFile,setPeakEvts)
    
    #take a runHLT config file and make an extension to it
    #to use just those events that we just parsed and picked
    if (os.path.isfile(hltcfg_file) and hltcfg_file != "dummy"):
        extend_HLTcfg_file(hltcfg_file,peakEvtsFile)
        print "File "+hltcfg_file+" was expanded to process only"
        print "events of interest.  Now you can run something like:"
        print "cmsRun "+hltcfg_file+" >& mylog.log"
    else:
        print "No HLT config file to modify.  No action taken"

    return 0

    
    


###############################################################
def main():
###############################################################

    #TO DO: work on more sophisticated parsing of parameters

    #check the number of parameters
    numarg = len(sys.argv)
    if numarg < 3:
        usage()
        return 1

    tmin = sys.argv[1]
    tmax = sys.argv[2]
    if (numarg == 4):
        hltcfg_file = sys.argv[3]
    else:
        hltcfg_file = "dummy"

    if (DEBUG):
        print tmin
        print tmax

    #check for crazy requirements
    if (float(tmin)<0 or float(tmax)<0 or float(tmax)<float(tmin)):
        usage()
        sys.exit(1)

    #run the script
    get_eventNumbers(tmin,tmax,hltcfg_file)    
    





###############################################################
if __name__ =='__main__':
###############################################################    
    sys.exit(main())



