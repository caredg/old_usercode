#!/usr/bin/python

##############################################################
# Creator: Edgar Carrera
# Date: Oct 1, 2009
# 
#
#
# 
#############################################################
import os, string, re,sys
from time import gmtime, localtime, strftime


#just a debug flag
DEBUG = False



###############################################################
def usage():
###############################################################

    if (DEBUG):
        print "This is the usage function"
        
    print '\n'
    print 'Usage: '+sys.argv[0]+' <text summary log file>'
    print 'e.g.:  '+sys.argv[0]+' full.log\n'
    print "This script parses the HLT-Report and TrigReport"
    print "sections of the a log file obtained after running"
    print "HLT win CMSSW.  It basically searches for"
    print "those paths that didn't finish running and"
    print "the modules that were not executed.  The first"
    print "module appearing is always the one that caused"
    print "the failure of the path"


###############################################################
def parse_TrigReport(infile):
###############################################################

    fin = open(infile,"r")
    allPaths = {}
    activateParsing = False
    for line in fin.readlines():
        linestuff = line.split()
        aux = line.find('TrigReport')
        if (aux == 0):
            saux = line.find('Modules in Path')
            saux2 = line.find('Module Summary')
            if (saux != -1):
                activateParsing = True
                hltPath = linestuff[5]
                allPaths[hltPath] = []
            if (saux2 != -1):
                activateParsing = False
            if (activateParsing):
                if (linestuff[3]==linestuff[5]):
                    allPaths[hltPath].append(linestuff[7])

    print "\nPATHS THAT DID NOT FINISH RUNNING "
    print "FIRST MODULE INDICATES THE ONE CAUSING THE FAILURE"
    print "THE REST OF THE MODULES WERE NEVER RUN IN THE PATH"

    mycounter = 1
    for k in allPaths:
        if allPaths[k]:
            print "\n"+str(mycounter) +". %%%%%%%%%%% PATH "+k+":"
            for m in allPaths[k]:
                print m
            mycounter+=1
                

    return 0




###############################################################
def parse_HLTReport(infile):
###############################################################

    fin = open(infile,"r")
    quitePaths = []
    for line in fin.readlines():
        aux = line.find('HLT-Report')
        linestuff = line.split()
        if (aux == 0 and linestuff[1].isdigit()):
            if linestuff[3] == "0":
                quitePaths.append(linestuff[5])

    print "\nPaths that didn't fire: "
    for k in quitePaths:
        print k

    return 0
            



###############################################################
def main():
###############################################################


    #check the number of parameters
    numarg = len(sys.argv)
    if numarg < 2:
        usage()
        return 1

    infile = sys.argv[1]
    
    #check if input file exists
    if  not(os.path.isfile(infile)):
        print infile+" does not exist. Please check."
        sys.exit(1)
    
    #run the script
    parse_HLTReport(infile)
    parse_TrigReport(infile)
    



###############################################################
if __name__ =='__main__':
###############################################################    
    sys.exit(main())

