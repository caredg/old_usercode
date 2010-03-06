#!/usr/bin/env python
############################################################################
#
# Author: Edgar Carrera
# ecarrera@cern.ch
# Script to search for core files in filter unit machines.  it is based
# on the timestamp of the corefiles and a given range of time.
############################################################################

"""
   usage: %prog [options]
   -s, --startdate = STARTDATE: Required start date of run/range in this format: 'Feb 29, 18:30'. 
   -e, --enddate = ENDDATE: Optional end date of run/range in this format: 'Mar 1, 06:00' Default is now.
"""


import os,sys
import datetime
import string, re
import fileinput
import commands
import time

#flag to debug
debugme = False
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
def run_parsing(dicOpt):
#######################################################

    runstart = dicOpt['startdate']
    runend = dicOpt['enddate']
    
    #this year
    currentYear = str(datetime.datetime.now().year)

    #time mask to be used
    maskme = "%Y %b %d %H:%M"

    #parse the startdate option and get epoch seconds
    sday = runstart.split(",")[0]
    shour = runstart.split(",")[1].lstrip()
    stime = currentYear+" "+sday+" "+shour
    sepoch = int(time.mktime(time.strptime(stime, maskme)))

    #parse the enddate option and get epoch in seconds
    if runend == "":
        eepoch = time.time()
    else:    
        eday = runend.split(",")[0]
        ehour = runend.split(",")[1].lstrip()
        etime = currentYear+" "+eday+" "+ehour
        eepoch = int(time.mktime(time.strptime(etime, maskme)))
    
    #open the log file
    myfile = open("error.log","r")

    #initialize bufu name
    mybufu = "N/A"

    #loop over lines in error log
    for line in myfile:
        hasbufu = line.find("bufu-")+1
        mysline = line.split()
        mssize = len(mysline)
        #if line has bufu in it and core info
        #just get the bufu name and the core file
        if hasbufu and mssize>1:
            #protect against many bufus in one line
            lastbufu = mysline[mssize - 1].find("bufu-")+1
            if lastbufu:
                mybufu = mysline[mssize-1]
            elif not lastbufu:
                mybufu = mysline[0]
                mycore = mysline[9]
                mystime = mysline[6]+" "+mysline[7]+" "+mysline[8]
                ctime = currentYear+" "+mystime
                cepoch = int(time.mktime(time.strptime(ctime, maskme)))
                #make a decision to print or not
                if cepoch>=sepoch and cepoch<=eepoch:
                    print mybufu+mycore
                    if debugme:
                        print mybufu+" "+mycore+"\t"+mystime
        #if line has bufu but nothing else store the
        #name until you find a core that fits in the range
        elif hasbufu and mssize == 1:
            mybufu = mysline[0].rstrip()
        elif not hasbufu and mssize >=9:
            mycore = mysline[8]
            mystime = mysline[5]+" "+mysline[6]+" "+mysline[7]
            ctime = currentYear+" "+mystime
            cepoch = int(time.mktime(time.strptime(ctime, maskme)))
            #make a decision to print or not
            if cepoch>=sepoch and cepoch<=eepoch:
                print mybufu+mycore
                if debugme:
                    print mybufu+" "+mycore+"\t"+mystime

                
#######################################################
def get_default_options(option):
#######################################################
    dicOpt = {}
    
    if not option.startdate:
        print "you need a run/range start date"
        exit()
    else:
        dicOpt['startdate'] = str(option.startdate)

    if not option.enddate:
        dicOpt['enddate'] = ""
    else:
        dicOpt['enddate'] = str(option.enddate)

    return dicOpt




#######################################################
def get_default_options(option):
#######################################################
    dicOpt = {}
    
    if not option.startdate:
        print "you need a run/range start date"
        exit()
    else:
        dicOpt['startdate'] = str(option.startdate)

    if not option.enddate:
        dicOpt['enddate'] = ""
    else:
        dicOpt['enddate'] = str(option.enddate)

    return dicOpt





#######################################################
if __name__ =='__main__':
#######################################################

    #import optionparser
    option,args = parse(__doc__)
    if not args and not option:
        exit()

    if not option.startdate:
        print "\nNo filtering will be done, you need at least the start date."
        exit()

    #set default options
    dicOpt = get_default_options(option)

    if (printConfig):
        for k in dicOpt:
            print str(k)+" = "+str(dicOpt[k])    

    run_parsing(dicOpt)

    
