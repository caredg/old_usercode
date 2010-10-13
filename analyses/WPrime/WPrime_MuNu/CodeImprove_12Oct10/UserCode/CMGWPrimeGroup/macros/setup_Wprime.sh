#!/bin/bash
if [ -f /etc/bashrc ]; then
      . /etc/bashrc
else
echo -e " No /etc/bashrc file! So confused..."
echo -e " Reflect, repent, and try with a better script next time"
echo -e " Exiting...."
exit
fi

# Definitions
RELEASE_VERSION=CMSSW_3_8_5
WORKING_AREA=V150_385
# end definitions

export RELEASE_VERSION WORKING_AREA

echo -e  "******************************************"
echo -en " Setting up new working area: "
echo -e $WORKING_AREA
echo -e  "******************************************"

#setup new working area
mkdir $WORKING_AREA
cd $WORKING_AREA
scramv1 p CMSSW $RELEASE_VERSION
cd $RELEASE_VERSION/src

echo -e  "\n**************************"
echo -e  " Checking out the code..."
echo -e  "**************************"
cvs -Q co -r OCT13-10_V3 -d UserCode UserCode/ecarrera/analyses/WPrime/WPrime_MuNu/CodeImprove_12Oct10/UserCode
# cvs -Q co UserCode/CMGWPrimeGroup

echo -e "\n************************"
echo -e " Done. Now will compile"
echo -e "************************"
sleep 2
source /afs/cern.ch/cms/sw/cmsset_default.sh
# Not sure how to make script aware of cmsenv alias...
#cmsenv
scramv1 runtime -sh
scramv1 b # (lots of output here, but nothing to worry about)
#cd ../../../

echo -e "\n*************************************************"
echo -e " Done compiling UserCode/CMGWPrimeGroup"
echo -e " Now making symbolic link to example config file"
echo -e "*************************************************"
echo -e " \n"
ln -s UserCode/CMGWPrimeGroup/test/myanalysis_TeVMuon_cfg.py .
