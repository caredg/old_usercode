#include "UserCode/CMGWPrimeGroup/root_macros/input_file.h"
#include "UserCode/CMGWPrimeGroup/interface/wprimeEvent.h"

#include <string>
#include <vector>
#include <iostream>
#include <fstream>

#include <TTree.h>
#include <TFile.h>

using std::string; using std::vector;
using std::cout; using std::endl; using std::cerr;

void parseLine(const string & new_line, wprime::InputFile * in_file, 
	       float & lumi_ipb);

string top_level_dir = "";

void loadInputFiles(vector<wprime::InputFile> & files, float & lumiPb)
{
  ifstream topdir_file("UserCode/CMGWPrimeGroup/config/top_directory.txt");
  getline(topdir_file, top_level_dir);
  if(top_level_dir.empty())
    {
      cerr << " *** Failed to load top level directory! " << endl;
      return;
    }

  cout << "\n Reading root-tuples from directory " << top_level_dir << endl;

  ifstream in("UserCode/CMGWPrimeGroup/config/samples_cross_sections_v8.txt");
  string new_line; wprime::InputFile * new_file = 0;
  while(getline(in, new_line))
    {
      // if DONE, we are done!
      if(new_line == "DONE")break;

      if(new_line.find("samplename = ") != string::npos)
	// new file found! create structure to put in info
	new_file = new wprime::InputFile();
      
      if(!new_line.empty())
	parseLine(new_line, new_file, lumiPb);
      else
	{
	  if(new_file->samplename != "data")
	    new_file->weight = lumiPb*(new_file->x_sect)
	      /(new_file->Nprod_evt);
	  else
	    new_file->weight = 1;
	  cout << "; Weight to be applied on " << new_file->samplename
	       << " sample = " << new_file->weight << endl << endl;
	  // all info should now be logged in; check!
	  new_file->checkFile();
	  // if we made it here, everything looks good: add to vector of MC samples
	  files.push_back(*new_file);
	  // release memory
	  delete new_file;
	}
    }
}

void parseLine(const string & new_line, wprime::InputFile * in_file, 
	       float & lumi_ipb)
{
  unsigned int i = 0;
  i = new_line.find("samplename = ");
  if(i != string::npos)
    {
      in_file->samplename = new_line.substr(13, new_line.length() - 13);
      return;
    }

  i = new_line.find("description = ");
  if(i != string::npos)
    {
      in_file->description = new_line.substr(14, new_line.length() - 14);

      // this is where we extract the integrated luminosity
      // from the description of the data sample; expect something like
      // "blah-blah-blah random description, XYZ pb^-1" ==> XYZ is what we want
      // NB: the data sample must appear first in samples_cross_sections.txt !
      if(lumi_ipb < 0)
	{
	  unsigned int begin = in_file->description.find(",") + 2;
	  unsigned int end = in_file->description.find("pb^-1");
	  string integ_lumi = in_file->description.substr(begin, end-begin);
	  lumi_ipb = atof(integ_lumi.c_str());

	  cout << "\n Will apply weights to MC samples to get distributions for "
	       << lumi_ipb << " pb^-1" << endl << endl;
	  if(lumi_ipb <= 0)
	    {
	      cout << " *** Oops, read in \"" << integ_lumi 
		   << "\" which I failed to parse! Please seek help... " 
		   << endl;
	      abort();
	    }
	}
      return;
    }

  i = new_line.find("pathname = ");
  if(i != string::npos)
    {
      string pathname = top_level_dir + new_line.substr(11, new_line.length() - 11);
      in_file->pathname = pathname;

      TFile * file = new TFile(pathname.c_str());
      if(!(file->IsOpen()))
	{
	  cerr <<" *** Missing file: "<< pathname << " !!! "<<endl; 
	  abort();
	}
      in_file->tree = (TTree *) file->Get("StdMu/wprime");
      if(! in_file->tree->GetBranch("wp"))
	{
	  cerr << " *** Can't find wp branch in file: " << pathname<< endl;
	  in_file->tree = 0;
	  return;
	}
      
      cout << " Opened file: " << in_file->samplename;
      cout << " (" << in_file->description << ") " << endl;
      cout << " # of entries = " << in_file->tree->GetEntries();

      return;
    }

  i = new_line.find("x-section = ");
  if(i != string::npos)
    {
      in_file->x_sect = atof(new_line.substr(12, new_line.length() - 12).c_str());
      return;
    }

  i = new_line.find("Nprod_evt = ");
  if(i != string::npos)
    {
      in_file->Nprod_evt = atoi(new_line.substr(12, new_line.length() - 12).c_str());
      return;
    }
  
}
