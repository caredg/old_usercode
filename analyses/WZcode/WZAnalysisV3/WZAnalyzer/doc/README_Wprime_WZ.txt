----------------------------------------------------------
THIS DOCUMENT:

- Instructions on how to run the "ExecuteAnalysis" macro
- Description of how the W'-> WZ (ExecuteAnalysis) macro
  works

WARNING!!! These instructions will be updated constantly as
the work progresses.  Not everything is implemented yet.
Please contribute if you feel like it!!

Edgar Carrera
Jan 13, 2010: version original
----------------------------------------------------------




-------------------------------------
REQUIREMENTS, ASSUMPTIONS, and NOTES:
-------------------------------------


- Work in an area where ROOT can be called.  It can very well be 
your laptop.  No CMSSW libraries required.
- Have the necessary signal and background files (in *.root format) 
in a directory that can be directly accessed.  
In the below example, the absolute path to them
is defined in the *.h file as
/localdata/ecarrera/analyses/wprime_to_wz/root_uples/
- The ROOT Trees in the root-tuples are named "WZ"
- Each root-uple has a histogram called numEvents, which give information
of the total number of events that were processed.  This is according to
the instruction at:
https://twiki.cern.ch/twiki/bin/viewauth/CMS/WZWorkPageSummer09
- These macros are inspired in the original macro "WZHistogramMaker.C"
written by Jeff Klukas that can be found at:
http://cmssw.cvs.cern.ch/cgi-bin/cmssw.cgi/UserCode/Klukas/WZAnalysisV3/WZAnalyzer/test/WZHistogramMaker.C?view=log
- $> represents the shell prompt



----------------------------
HOW TO RUN THE MACROS
----------------------------

1. Download or copy the macros to your favorite area.  In the case of the
example: wprime_wz_macros.  NOTE: You can browse the code at: http://cmssw.cvs.cern.ch/cgi-bin/cmssw.cgi/UserCode/ecarrera/analyses/WZcode/WZAnalysisV3/WZAnalyzer/test/wprime_wz_macros/


$> cvs co -d wprime_wz_macros UserCode/ecarrera/analyses/WZcode/WZAnalysisV3/WZAnalyzer/test/wprime_wz_macros

$> cd wprime_wz_macros

NOTE: You can browse the code at: http://cmssw.cvs.cern.ch/cgi-bin/cmssw.cgi/UserCode/ecarrera/analyses/WZcode/WZAnalysisV3/WZAnalyzer/test/wprime_wz_macros/


2. Run the ROOT macro in a compiled (+) form:

$wprime_wz_macros> root -l ExecuteAnalysis.C+


3. After it is finished, you will find:

a. The printout that includes the number of expected (weighted) events after
each cut, as well as absolute and relative efficiencies for all the samples
that were input. (For now the cuts
are just some demo cuts).

b. A ROOT file called Wprime_analysis.root that contains the signal and
background histograms (after each cut) nicely separated into directories so
the plotting can be done much more easily.



----------------------------------
DESCRIPTION OF HOW THE MACRO WORKS
----------------------------------

Objective: The sole objective of this macro is to run over background and
signal files that were generated and produced with the procedures described at:
https://twiki.cern.ch/twiki/bin/view/Sandbox/WprimeToWZGenerationInstructions
https://twiki.cern.ch/twiki/bin/viewauth/CMS/WZWorkPageSummer09
In the example above, all the files (signal and background) live at:
/localdata/ecarrera/analyses/wprime_to_wz/root_uples/
The idea is to be able to run the whole analysis quickly over the samples
available, as well as being able to expand and perfect the code in a clean
way; all this while keeping the simplicity of the macro.


The macro works with two files, ExecuteAnalysis.C and ExecuteAnalysis.h.

---++ ExecuteAnalysis.h:

This file defines several variables used in the analysis.  In particular,
it creates a "struct" to store the information about the different signal
and background files.  

Also, here is where the analysis variables that we are going to need 
(those present in the root-uples) are declared.  For example:

// +++++++++++++++++++Variables to store Branch Addresses:
Int_t           W_flavor;
Int_t           Z_flavor;
Int_t           numberOfZs;
Float_t         WZ_invMassMinPz;

You can quickly add a new variable here if needed.

The "Num_histo_sets" variable will be used to match, generally, the number of
cuts applied in the analysis.  That way we can see different distributions
(like the WZ invariant mass) after each cut is applied.


---++ ExecuteAnalysis.C:

Within this file, the method "ExecuteAnalysis()" drives the whole analysis.
In general, the macro is very well commented, so its understanding is
easier.

Basically, we define a *.root file called "Wprime_analysis.root" where we can
store all the histograms produced and a text file called "event_counts.txt"
where we will put all the numerical results.

Containers are created for the signal and background files.  Ex:
vector<InputFile> wzjj_files;
vector<InputFile> wprime400_files;

We load the cross sections for the different processes using the method
"Load_Cross_Sections()".  Note that this method allows you to load cross
sections for samples generated in bins, like it will be the case for QCD
multijet background.  The argument load should be expanded as we
incorporate more signal and background samples.

We then iniciate the analysis by loading the input files using the method
Load_Input_Files() which acquires files sample by sample.  Here, the
weights are set for the events according to the integrated luminosity that
was passed to the method.

Once we have loaded the files, we use the method Get_Distributions() to
loop over each file corresponding to each sample and the events contained
in there, applying the required cuts and store the number of events that
pass each cut.  We fill whatever histograms we have declared and defined in
the "ExecuteAnalysis.h" file and in the method Declare_Histos(),
respectively. This can be eventually expanded easily.
Finally we print the summary using the printSummary() method and save the
histograms using the saveHistos() method, which nicely separate samples in directories.


