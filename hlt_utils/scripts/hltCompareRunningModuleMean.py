#!/usr/bin/env python
############################################################################
#
# Edgar Carrera
# ecarrera@cern.ch
#
# This script compares the module running time for each module for a given 
# list of paths between two measurements (usually made on data taken with
# different conditions)
# One assumes here that the same HLT menu was used for both kind of data 
# File1 is usually the reference while File2 is the HPU file or a measurement
# with other conditions

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
DEBUG = True


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
def get_paths_and_modules_repo(infile, pathlist):
###############################################################
    #prepare the two container we will need later 	
    repo = {}
    f = TFile(infile,"READ")
    dir = f.GetListOfKeys()
    for k in dir:
        h = k.ReadObj()
	allnames = h.GetName()
        ishist = allnames.split("_").count("moduleInPathScaledTimeSummary")
	#get the modules for the paths of interest	
	if ishist > 0:
            fullname = allnames.split("_")
            fullname.remove("moduleInPathScaledTimeSummary")
            #form path name
            pathname = ""
            for j in range(len(fullname)):
                pathname = pathname + fullname[j]
                if j is not (len(fullname)-1):
                    pathname = pathname + "_"
            if not pathname in pathlist:
		continue
	    if not pathname in repo:
                repo[pathname] = []
	    else:
		print "repeated path, please check ....quiting"
		sys.exit(1)
            #get the module names in order for this path
	    #zero the repos
	    nbins = h.GetXaxis().GetNbins()
	    for b in range(nbins):
		modname = h.GetXaxis().GetBinLabel(b+1)
		repo[pathname].append(modname)
    #print repo	
    return repo
	    
###############################################################
def get_histos_info(infile, modAndPathList):
###############################################################
    repo = {} #for entries
    repom = {} #for mean values
    f = TFile(infile,"READ")
    dir = f.GetListOfKeys()
    for k in dir:
        h = k.ReadObj()
        allnames = h.GetName()
	nentries = h.GetEntries()
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
            #get the pathlist
	    pathlist = []
	    for path in modAndPathList:
	 	pathlist.append(path)
	    #skip if histo is not related to our paths
	    if not pathname in pathlist:
		continue	
	    
            #get module name
            modname = fullname[len(fullname)-1]
            #store information paths and mods

            if not pathname in repo:
            	repo[pathname] = {}
            	repom[pathname] = {}

            repo[pathname][modname] = nentries
            repom[pathname][modname] = themean
#    print repo
#    print repom	
    return repo,repom            



###############################################################
def calculate_rejection_info(repo,modAndPathList):
###############################################################
	repoj = {}
	for path in repo:
		if DEBUG: print "%%%%%%%%%%%% "+path
		if not path in repoj:
			repoj[path] = {}
		for modIdx in range(len(modAndPathList[path])):
			modname = modAndPathList[path][modIdx]
			if modIdx is not (len(modAndPathList[path])-1):
				if repo[path][modname] > 0:
					rejFact = 1-(repo[path][modAndPathList[path][modIdx+1]]/repo[path][modname])
					if DEBUG:
						print modname+": "+str(repo[path][modname])+", "+str(modAndPathList[path][modIdx+1])+": "+str(repo[path][modAndPathList[path][modIdx+1]])+", rejection: "+str(rejFact)

				else:
					rejFact = 0.
			else: rejFact = 0.
			repoj[path][modname] = rejFact
	#print repoj
	return repoj	


###############################################################
def plot_results(repo1,repo2,modAndPathList):
###############################################################

    #Always start with the reference container, which here we assume
    #is number one
    plen = len(repo1)
    container1 = repo1
    container2 = repo2
    
    #start plotting path by path
    mycountp = 0
    hp = TH1F("hp","Ratio of Path's running time",plen,0,plen)
    for p in container1:
       
        #continue if the path is not in both containers
        isPath = container2.get(p)
        if isPath is None:
            print "No plot was made for path..please check %s" % p
            continue
        
	#use the modules in the reference path
        nbins1 = len(container1[p])
        nbins2 = len(container2[p])
        modules1 = container1[p]
        modules2 = container2[p]
        nbins = nbins1
        

        #define histograms
	type_of_study = "means"
	ffname1 = "means_ref"
	ffname2 = "means_HPU"
        hname1 = ffname1+"_" + p
        hname2 = ffname2+"_" + p
        h1 = TH1F(hname1,p,nbins,0,nbins)
        h2 = TH1F(hname2,p,nbins,0,nbins)
        
        #fill histos
        maxmean = 0
	for modIdx in range(len(modAndPathList[p])):
	    m = modAndPathList[p][modIdx]
            h1.Fill(m,modules1[m])
            #check for max val of mean (it actually works for any quantity)
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
              
        #fill the total running path time ratio histo
        ratiopaths = h2.Integral()/h1.Integral()
        #print " "
        #print p
        #print ratiopaths
        hp.Fill(p,ratiopaths)
        
        #draw the two histos and save the plot
        can = TCanvas("can","",1800,800);
        can.cd();
        h1.Draw();
        h1.SetLineWidth(2)
	#check if the numbers in the histograms are fractions or whole numbers
	#in order to adjust the maximum
	if maxmean <= 1:
		h1.SetMaximum(1.1)
	        h1.GetYaxis().SetTitle("rejection fraction")
	else:
	        h1.SetMaximum(maxmean+10)
	        h1.GetYaxis().SetTitle("msec")		

        h1.GetXaxis().SetLabelSize(0.03)
        h2.SetLineColor(2);
        h2.Draw("same")
        h2.SetLineWidth(2)
        can.SetBottomMargin(0.55)
        can.SetLogy()
        can.SetGridx()
        can.SetGridy()
        #lg = TLegend(0.59, 0.67, 0.89, 0.89);
        lg = TLegend(0.79, 0.89, 1.0, 0.99);
        lg.AddEntry(h1,ffname1,"L");
        lg.AddEntry(h2,ffname2,"L");
        lg.SetTextSize(0.025);
        lg.SetBorderSize(0);
        lg.SetFillColor(0);
        lg.Draw("same");
        canname = type_of_study+"_"+p + ".gif"
        can.Print(canname)
        del can
        
        mycountp+=1

    #plot the path summary
    can2 = TCanvas("can2","",1500,600)
    can2.cd()
    hp.Draw()
    hp.SetLineWidth(2)
    can2.SetBottomMargin(0.5)
    can2.SetGridx()
    can2.SetGridy()
    hp.GetXaxis().SetLabelSize(0.03)
    can2.Print("pathRunningTimeRatio.gif")
    
        
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

	
    #Give a list of paths:
   # pathlist = ['HLT_Ele80_CaloIdVT_TrkIdT_v3','HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v10','HLT_Photon36_CaloIdL_IsoVL_Photon22_CaloIdL_IsoVL_v7','HLT_Photon90EBOnly_CaloIdVL_IsoL_TriPFJet25_v5','HLT_HT350_Ele5_CaloIdVL_CaloIsoVL_TrkIdVL_TrkIsoVL_PFMHT45_v11','HLT_HT350_MHT110_v3','HLT_Jet370_NoJetID_v10','HLT_Jet370_v10','HLT_QuadJet80_v5','HLT_MET200_v7','HLT_PixelTracks_Multiplicity100_v8','HLT_Mu17_TkMu8_v4','HLT_Mu8_Ele17_CaloIdT_CaloIsoVL_v8','HLT_Mu10_R029_MR200_v6','HLT_Dimuon0_Jpsi_Muon_v11','HLT_Mu100_eta2p1_v5','HLT_IsoMu30_eta2p1_v7', 'HLT_DoubleIsoPFTau45_Trk5_eta2p1_v8']
  
    #slowest paths in HPU
    #pathlist = ['HLT_CentralJet46_CentralJet38_DiBTagIP3D_v7', 'HLT_CentralJet60_CentralJet53_DiBTagIP3D_v6', 'HLT_QuadJet45_IsoPFTau45_v13', 'HLT_R020_MR300_CentralJet40_BTagIP_v4', 'HLT_R030_MR200_CentralJet40_BTagIP_v4', 'HLT_R014_MR400_CentralJet40_BTagIP_v4', 'DST_HT350_RunPF_v1', 'HLT_R014_MR450_CentralJet40_BTagIP_v4', 'HLT_R020_MR350_CentralJet40_BTagIP_v4', 'HLT_DoubleJet60_ForwardBackward_v10']

    #slowest paths in HPU without hltPixelTracks
    pathlist = ['HLT_CentralJet46_CentralJet38_DiBTagIP3D_v7', 'HLT_CentralJet60_CentralJet53_DiBTagIP3D_v6', 'HLT_QuadJet45_IsoPFTau45_v13', 'HLT_R020_MR300_CentralJet40_BTagIP_v4', 'HLT_R030_MR200_CentralJet40_BTagIP_v4', 'HLT_R014_MR400_CentralJet40_BTagIP_v4', 'DST_HT350_RunPF_v1', 'HLT_R014_MR450_CentralJet40_BTagIP_v4', 'HLT_R020_MR350_CentralJet40_BTagIP_v4', 'HLT_TkIso10Mu5_Ele8_CaloIdT_CaloIsoVVL_TrkIdVL_Mass8_HT150_v5']

    #get ordered lists for modules in the paths
    modAndPathList = get_paths_and_modules_repo(infile1, pathlist)
    #fill out the repos with the right information
    repo_entries1, repo_means1 = get_histos_info(infile1,modAndPathList)
    repo_entries2, repo_means2  = get_histos_info(infile2,modAndPathList)
    #repo_rejection1 = calculate_rejection_info(repo_entries1, modAndPathList)
    #print repo_rejection1
    #repo_rejection2 = calculate_rejection_info(repo_entries2, modAndPathList)
    #print repo_rejection2
    #plot results
    plot_results(repo_means1,repo_means2,modAndPathList)

#######################################################
if __name__ =='__main__':
#######################################################
    sys.exit(main())
