#!/usr/bin/env python
############################################################################
#
# Edgar Carrera
# ecarrera@cern.ch
#
# This script compares the running modules within paths of two
# timing measurements (files).
#  It looks first at the file with more paths and starts with this, then
#  looks for the same path in the second file.
#  If found, it goes module by module checkin if they match,
#  and if they do, it compares their mean running time.
#  The script also prints those cases where the discrepancy
#  is large either by looking at their ratio or their difference. These
#  thresholds are set in the header of the script.
#  The naming of the legends follow the names of the files.
#  One can run the script as
#  ./hltCompareRunningModules <file1.root> <file2.root> -b >! list.txt
# The option -b will put the script running in batch mode.
# The list.txt file will contain the cases with large discrepancies.
###########################################################################

"""
   usage: %prog <file1> <file2>
 
"""


import os,sys
import subprocess
import string, re
import fileinput
import commands
from time import gmtime, localtime, strftime
from ROOT import *

#just a debug flag
DEBUG = False
#threshold for the ratio of the running modules
THRATIO = 10.0
#threshold for the difference of the running modules, used
#when one instance is zero
THDIFF = 200

gROOT.SetStyle("Plain")
gStyle.SetOptStat(0)

###############################################################
def usage():
###############################################################

    if (DEBUG):
        print "This is the usage function"
        
    print '\n'
    print 'Usage: '+sys.argv[0]+' <file1> <file2>'
    print 'e.g.:  '+sys.argv[0]+' outputTiming1.root outputTiming2.root\n'


###############################################################
def get_histos_info(infile):
###############################################################
    repo = {}
    f = TFile(infile,"READ")
    dir = f.GetListOfKeys()
    for k in dir:
        h = k.ReadObj()
        allnames = h.GetName()
        themean = h.GetMean()
        ishist = allnames.split("_").count("moduleInPathScaledTime")
        #get only the running modules in path
        if ishist > 0:
            fullname = allnames.split("_")
            fullname.remove("moduleInPathScaledTime")
            #form path name
            pathname = ""
            for j in range(len(fullname)-1):
                pathname = pathname + fullname[j]
                if j is not (len(fullname)-2):
                    pathname = pathname + "_"
            #get module name
            modname = fullname[len(fullname)-1]
            #store information paths and mods
            if not pathname in repo:
                        repo[pathname] = {} #for modules and their means
            repo[pathname][modname] = themean


    #print to check
    if (DEBUG):
        mycounter = 1
        for p in repo:
            printPath = True
            for m in repo[p]:
                if printPath:
                    print "\n"+str(mycounter) +". %%%%%%%%%%% PATH "+p+":"
                    mycounter+=1
                    printPath = False
                #print m +": "+str(repo[p][m])

    return repo;
            



###############################################################
def plot_results(repo1,repo2,file1,file2):
###############################################################

    #decide which container to start with
    #start with the one with more paths
    len1 = len(repo1)
    len2 = len(repo2)
    if len1 > len2:
        container1 = repo1
        container2 = repo2
        fname1 = file1.strip(".root")
        fname2 = file2.strip(".root")
    else:
        container1 = repo2
        container2 = repo1
        fname1 = file2.strip(".root")
        fname2 = file1.strip(".root")

    #start plotting path by path
    mycountp = 0
    for p in container1:
        if DEBUG:
            if mycountp > 10:
                break
            
        #continue if the path is not in both containers
        isPath = container2.get(p)
        if isPath is None:
            print "No plot was made for path %s" % p
            continue
        
        #decide which path container to use
        #start with the one with more modules
        nbins1 = len(container1[p])
        nbins2 = len(container2[p])
        if nbins1 > nbins2:
            modules1 = container1[p]
            modules2 = container2[p]
            nbins = nbins1
            ffname1 = fname1
            ffname2 = fname2
        else:
            modules1 = container2[p]
            modules2 = container1[p]
            nbins = nbins2
            ffname1 = fname2
            ffname2 = fname1

        #define histograms
        hname1 = ffname1 + p
        hname2 = ffname2 + p
        h1 = TH1F(hname1,p,nbins,0,nbins)
        h2 = TH1F(hname2,p,nbins,0,nbins)
        
        #fill histos
        maxmean = 0
        for m in modules1:
            h1.Fill(m,modules1[m])
            #check for max val of mean
            if modules1[m] > maxmean:
                maxmean = modules1[m]
            #check if the other container has the module
            theModMean = modules2.get(m)
            if theModMean is None:
                h2.Fill(m,0)
                ratiomod = 1
                print "Module %s in path %s was no found in one of the files" % m,p
            else:
                #check for the max val of mean
                if theModMean > maxmean:
                    maxmean = theModMean
                h2.Fill(m,theModMean)
                #check if the discrepancy is large and print
                if theModMean > 0 and modules1[m] > 0:
                    ratiomod = theModMean/modules1[m]
                    diffmod = abs(modules1[m] - theModMean)
                    if ratiomod > 0 and (ratiomod >= THRATIO or ratiomod <= (1/THRATIO)):
                        print "Path %s, module %s: mean1 = %s ms, mean2 = %s ms, ratio = %s, diff = %s ms" % (p,m,modules1[m],theModMean,str(ratiomod),str(diffmod))
                else:
                    ratiomod = 0.
                    diffmod = abs(modules1[m] - theModMean)
                    if diffmod > 0 and diffmod > THDIFF:
                        print "Path %s, module %s: mean1 = %s ms, mean2 = %s ms, ratio = %s, diff = %s ms" % (p,m,modules1[m],theModMean,str(ratiomod),str(diffmod))
                
        #draw the two histos and save the plot
        can = TCanvas("can","",1500,600);
        can.cd();
        h1.Draw();
        h1.SetLineWidth(2)
        h1.SetMaximum(maxmean+10)
        h1.GetYaxis().SetTitle("msec")
        h1.GetXaxis().SetLabelSize(0.03)
        h2.SetLineColor(2);
        h2.Draw("same")
        h2.SetLineWidth(2)
        can.SetBottomMargin(0.5)
        can.SetLogy()
        #lg = TLegend(0.59, 0.67, 0.89, 0.89);
        lg = TLegend(0.79, 0.89, 1.0, 0.99);
        lg.AddEntry(h1,ffname1,"L");
        lg.AddEntry(h2,ffname2,"L");
        lg.SetTextSize(0.03);
        lg.SetBorderSize(0);
        lg.SetFillColor(0);
        lg.Draw("same");
        canname = p + ".gif"
        can.Print(canname)
        del can
        
        mycountp+=1
        
###############################################################
def main():
###############################################################
    #check the number of parameter
    numarg = len(sys.argv)
    if numarg < 2:
        usage()
        return 1

    infile1 = sys.argv[1]
    infile2 = sys.argv[2]

    #check if input files exist
    if  not(os.path.isfile(infile1)):
        print infile1+" does not exist. Please check."
        sys.exit(1)
    if  not(os.path.isfile(infile2)):
        print infile2+" does not exist. Please check."
        sys.exit(1)


    #get the histos info in containers 
    repo1 = get_histos_info(infile1)
    repo2 = get_histos_info(infile2)

    #plot results
    plot_results(repo1,repo2,infile1,infile2)

#######################################################
if __name__ =='__main__':
#######################################################
    sys.exit(main())
