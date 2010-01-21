#!/usr/bin/python

##############################################################
# Author: Edgar Carrera
# Date: Oct 14, 2009
# 
#
#
# 
#############################################################
import os, string, re,sys
from time import gmtime, localtime, strftime


#just a debug flag
DEBUG = True



###############################################################
def usage():
###############################################################

    if (DEBUG):
        print "This is the usage function"
        
    print '\n'
    print 'Usage: '+sys.argv[0]+' <file with a list of log files>'
    print 'e.g.:  '+sys.argv[0]+' full.log\n'
    print "This script parses the HLT-Report and TrigReport"
    print "sections of the a group of log files obtained after running"
    print "HLT with CMSSW.  It basically searches for"
    print "those paths that didn't finish running and"
    print "the modules that were not executed.  The first"
    print "module appearing is always the one that caused"
    print "the failure of the path."


###############################################################
def parse_TrigReport(inflist):
###############################################################

    lin = open(inflist,"r")
    allPaths = {}
    for line in lin.readlines():
        if (DEBUG): print "File: "+ line.rstrip()
        fin = open(line.rstrip(),"r")
        activateParsing = False
        firstTime = False
        for fline in fin.readlines():
            flinestuff = fline.split()
            aux = fline.find('TrigReport')
            if (aux == 0):
                saux = fline.find('Modules in Path')
                saux2 = fline.find('Module Summary')
                if (saux != -1):
                    activateParsing = True
                    myIdx = 0
                    if not flinestuff[5] in allPaths:
                        hltPath = flinestuff[5]
                        allPaths[hltPath] = []
                        allPaths[hltPath].append([])#for modules
                        allPaths[hltPath].append([])#for # pass evts
                        firstTime = True
                    else:
                        hltPath = flinestuff[5]
                if (saux2 != -1):
                    activateParsing = False
                if (activateParsing):
                    if (flinestuff[4].isdigit()):
                        if (firstTime):
                            allPaths[hltPath][0].append(flinestuff[7])
                            allPaths[hltPath][1].append(int(flinestuff[4]))
                        else:
                            allPaths[hltPath][1][myIdx]+=int(flinestuff[4])
                            myIdx+=1
                            
                            
    print "\nPATHS THAT DID NOT FINISH RUNNING "
    print "FIRST MODULE INDICATES THE ONE CAUSING THE FAILURE"
    print "THE REST OF THE MODULES WERE NEVER RUN IN THE PATH"


    #order results
    Ks = allPaths.keys()
    Ks.sort()

    #print report
    mycounter = 1
    for k in Ks:
        printPath = True
        for m in range(len(allPaths[k][0])):
            if allPaths[k][1][m] == 0:
                if printPath:
                    print "\n"+str(mycounter) +". %%%%%%%%%%% PATH "+k+":"
                    mycounter+=1
                    printPath = False
                print allPaths[k][0][m]

                

    return 0




###############################################################
def parse_HLTReport(inflist):
###############################################################

    lin = open(inflist,"r")
    quitePaths = []
    allPaths = {}
    for line in lin.readlines():
        fin = open(line.rstrip(),"r")
        for fline in fin.readlines():
#            aux = fline.find('HLT-Report')
            aux = fline.find('TrigReport')
            if (aux == 0):
                saux = fline.find('Modules in Path') 
                if( saux != -1):
                    break
            flinestuff = fline.split()
            if (aux == 0 and flinestuff[1].isdigit()):
                if not flinestuff[7] in allPaths:
                    allPaths[flinestuff[7]] = int(flinestuff[4])
                else:
                    allPaths[flinestuff[7]]+=int(flinestuff[4])

    # sort keys for loop
    Ks = allPaths.keys()
    Ks.sort()

    # print results
    print "\nPrinting all Paths and # of events passsing..."
    for k in Ks:
        print k+": "+str(allPaths[k])
        if not allPaths[k] or allPaths[k] == 0:
            quitePaths.append(k)
            
    
    # prints paths that didn't fire        
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

    inflist = sys.argv[1]
    
    #check if input file exists
    if  not(os.path.isfile(inflist)):
        print inflist+" does not exist. Please check."
        sys.exit(1)
    
    #run the script
    parse_HLTReport(inflist)
    parse_TrigReport(inflist)
    



###############################################################
if __name__ =='__main__':
###############################################################    
    sys.exit(main())

