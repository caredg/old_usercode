#!/usr/bin/env python
############################################################################
#
# Edgar Carrera
# ecarrera@cern.ch
#
# Sep 07, 2010
# This script takes as input a root file with a timing measurement,
# and a HLT menu python configuration file or a txt file with a list of paths.
# It outputs a list of paths and the total mean time for the menu
# without that given path.
# A CMSSW area are needs to be set and the hltTimingSummary macro needs
# to be compatible with the options that this script calls.
# Also, a exclude.txt file is assumed to exist, containing the modules
# to exclude for I/O, otherwise the timing results will not be correct.
############################################################################

"""
   usage: %prog [options]
   -f, --hfile = HFILE: (Use this option XOR option -l) HLT menu python configuration file.  As extracted by hltGetConfiguration, for example.
   -l, --lfile = LFILE: (Use this option XOR option -f) Text file with a list of paths to exclude.
   -m, --mfile = MFILE: Root file containing the CPU timing measurement.
"""


import os,sys
import string, re
import fileinput
import commands
from time import gmtime, localtime, strftime
from ROOT import TFile, TH1F


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
def get_time_without_paths(hltpaths,dicOpt):
#######################################################

    themfile = dicOpt['mfile']
    results = {}
    #loop over paths
    for path in hltpaths:
        #printout information
        print "Working on path "+path+", be patient, this might take a while"
        #run the timing script
        time_str = "$CMSSW_BASE/test/$SCRAM_ARCH/hltTimingSummary -i "+themfile+" -o outfile -f -s -c -p 0 -e exclude.txt -x "+path
        os.system(time_str)
        #get the the total time info through PyROOT
        rootf = TFile("outfile.root","READ")
        h1 = rootf.Get("totalTime")
        totaltime = h1.GetMean()
        results[path] = totaltime

    return results


#######################################################
def parse_hltconfig(dicOpt):
#######################################################

    allpaths = []
    #if python file is given
    if dicOpt['hfile'] is not None:
        hltfile = dicOpt['hfile']
        thefile = open(hltfile,"r")
        for line in thefile.readlines():
            if (line.find('process.HLT_') != -1 or line.find('process.AlCa_') != -1 or line.find('process.HLTriggerFinalPath') != -1 or line.find('process.DQM_') != -1):
                path = line.split("=")[0].split(".")[1].rstrip()
                allpaths.append(path)
    #if text file is given
    else:
        hltfile = dicOpt['lfile']
        thefile = open(hltfile,"r")
        for line in thefile.readlines():
            if (line.find('HLT_') != -1 or line.find('AlCa_') != -1 or line.find('HLTriggerFinalPath') != -1 or line.find('DQM_') != -1):
                path = line.rstrip()
                allpaths.append(path)


    return allpaths



#######################################################
def get_default_options(option):
#######################################################
    dicOpt = {}

    dicOpt['mfile']= str(option.mfile)

    if not option.hfile:
        dicOpt['hfile']= None
    else:
        dicOpt['hfile']= str(option.hfile)
        
    if not option.lfile:
        dicOpt['lfile']= None
    else:
        dicOpt['lfile']= str(option.lfile)

    return dicOpt




#######################################################
if __name__ =='__main__':
#######################################################

    #import optionparser
    option,args = parse(__doc__)
    if not args and not option:
        exit()

    #safety checks 
    if not option.hfile and not option.lfile:
        print "provide the HLT menu python configuration file or a text file with a list of paths to be excluded"
        exit()
        
    if option.hfile and option.lfile:
        print "provide the HLT menu python configuration file or a text file with a list of paths to be excluded, but not both at the same time"
        exit()
        
    if not option.mfile:
        print "provide the name of the root file with the timing measurement"
        exit()


    #set default options
    dicOpt = get_default_options(option)

    #print configuration
    for k in dicOpt:
        print str(k)+" = "+str(dicOpt[k])


    #parse the python configuration file or txt file for
    #the list of HLT paths
    hltpaths = parse_hltconfig(dicOpt)
    
    #run the hltTimingSummary macro to get the numbers needed
    timingresults = get_time_without_paths(hltpaths,dicOpt)

    #print results (order first)
    K = timingresults.keys()
    K.sort()
    outf = open("timingresults.txt","w")
    print "RESULTS:"
    for path in K:
        finalres =  "%s: %0.3f\n" % (path,timingresults[path])
        print finalres.rstrip()
        outf.write(finalres)
    
    
