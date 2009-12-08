// Author: Edgar Carrera
// Root C++ script to loop over the information
// provided by the 
// DataFormats/HLTReco/test/hltTimingSummary.cpp
// macro
// It needs to be compiled.  You can do it by executing something like:
// rm hlttiming_histos.root hlt_timing ;g++ `root-config --cflags --glibs` -o hlt_timing hlt_timing.cpp

#include <vector>
#include <set>
#include <string>
#include <assert.h>
#include <stdio.h>
#include <algorithm>
#include <iostream>
#include <stdlib.h>
#include <TROOT.h>
#include <TKey.h>
#include <TObject.h>
#include <TClass.h>
#include <TCanvas.h>
#include <TFile.h>
#include <TH1.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TF2.h>
#include <TPaletteAxis.h>
#include <math.h>

using namespace std;


struct HLTData {
    
    double hlt_value;
    string hlt_name;
    string hlt_subname;
    HLTData(): hlt_name(""),hlt_subname(""),hlt_value(-1.0){}
    HLTData(string HLT_NAME,string HLT_SUBNAME,double HLT_VALUE):
        hlt_name(HLT_NAME),hlt_subname(HLT_SUBNAME),hlt_value(HLT_VALUE){}
    bool operator() (const HLTData &d1, const HLTData &d2) {
        return (d1.hlt_value > d2.hlt_value);
    }
    
};

TFile* f = new TFile("outfile.root","READ");
TH1F* hne = (TH1F*)f->Get("totalTime");
int totalNevents = hne->GetEntries();
struct HLTData mysort; //to nicely sort the data
bool debugme = false;

//Function to split a string into two pieces
//-----------------------------------------------------------------------------
void splitstring(string str, string seperater, 
                 string &first, string &second) 
{
//-----------------------------------------------------------------------------
     int i = (int)str.find(seperater); //find seperator
     if(i != -1)
     {
          int y = 0;
          if(!str.empty())
          {
               while(y != i)
               {
                    first += str[y++]; //creating first string
               }
               y = y+(int)seperater.length(); //jumping forward seperater length
               while(y != str.length())
               {
                    second += str[y++]; //creating second string
               }
               
          }
     }
     else
     {
          first = str;
          second = "NULL"; //if seperator is not there then second string == null 
     }
}//split string function






//Sort a vector decreasently by the second value in data structure
//--------------------------------------------------------------------------------
void sort_hlt_vector(vector<HLTData>& hlt_vector)
{
//--------------------------------------------------------------------------------
    if (debugme) cout<<"Sorting vector"<<endl;
    sort(hlt_vector.begin(),hlt_vector.end(),mysort);
    return;


}//sort_hlt_vector







//Save histogram
//------------------------------------------------------------------------------------
void save_histogram(TH1F* h)
{
//-------------------------------------------------------------------------------------
    if (debugme) cout<<"Saving histogram"<<endl;
    TFile myfile("hlttiming_histos.root","UPDATE");
    h->Write();
    myfile.Close();

}//Save histogram





//Save histogram
//------------------------------------------------------------------------------------
void save_histogram(TH2F* h)
{
//-------------------------------------------------------------------------------------
    if (debugme) cout<<"Saving histogram"<<endl;
    TFile myfile("hlttiming_histos.root","UPDATE");
    h->Write();
    myfile.Close();

}//Save histogram






//Loop over keys in file and return vector of histograms found
//by a string wildcard
//---------------------------------------------------------------------------
vector<TH1*> get_vector_of_histos (vector<string> wildcard )
{
//---------------------------------------------------------------------------
    if (debugme) cout<<"Getting vector from file"<<endl;
    vector<TH1*> histos_vector;
    TIter nextkey(f->GetListOfKeys());
    TKey *key;
    TH1* h;
    while ((key=(TKey*)nextkey())) {
        TObject* obj = key->ReadObj();
        if (obj->IsA()->InheritsFrom("TH1")) {
            h = (TH1*)obj;
            TString histName = h->GetName();
            for(int i=0;i<wildcard.size();++i){
                if (!histName.Contains(wildcard[i].c_str())) continue;
                histos_vector.push_back(h);
            }
        }
    }//---while key

    return histos_vector;

}//get_vector_of_histos







//fill a histogram with hlt_vector data
//---------------------------------------------------------------------------------------
TH1F* fill_hlt_histogram(const vector<HLTData>& hlt_vector)
{
//---------------------------------------------------------------------------------------   
    if (debugme) cout<<"Filling histograms"<<endl;
    Int_t _size = 0;
    for (vector<HLTData>::const_iterator it=hlt_vector.begin(); 
         it!=hlt_vector.end(); ++it){
        if (!it->hlt_value>0) continue;
        ++_size;
    }
    
    TH1F* hlt_hist = new TH1F("hlt_hist","",_size,0,_size);
    for (vector<HLTData>::const_iterator it=hlt_vector.begin(); 
         it!=hlt_vector.end(); ++it){
        if (!it->hlt_value>0) continue;
        hlt_hist->Fill(it->hlt_name.c_str(),it->hlt_value);
    }

    return hlt_hist;
}//fill a histo


//fill a histogram with hlt_vector data
//---------------------------------------------------------------------------------------
TH2F* fill2_hlt_histogram(const vector<HLTData>& hlt_vector)
{
//---------------------------------------------------------------------------------------   
    if (debugme) cout<<"Filling histograms"<<endl;
    Int_t _size = 0;
    set<string> _sizemodules;
    for (vector<HLTData>::const_iterator it=hlt_vector.begin(); 
         it!=hlt_vector.end(); ++it){
        if (!it->hlt_value>0) continue;
        _sizemodules.insert(it->hlt_name);
    }
    _size = _sizemodules.size();
    TH2F* hlt_hist = new TH2F("hlt_hist2","",_size,0,_size,100,0,100);
    for (vector<HLTData>::const_iterator it=hlt_vector.begin(); 
         it!=hlt_vector.end(); ++it){
        if (!it->hlt_value>0) continue;
        //cout<<it->hlt_name<<"\t"<<it->hlt_value<<endl;
        hlt_hist->Fill(it->hlt_name.c_str(),it->hlt_value,1);
    }

    return hlt_hist;
}//fill a histo



//fill a histogram with hlt_vector data
//---------------------------------------------------------------------------------------
TH2F* fill22_hlt_histogram(const vector<HLTData>& hlt_vector)
{
//---------------------------------------------------------------------------------------   
    if (debugme) cout<<"Filling histograms"<<endl;
    Int_t _sizemod = 0;
    Int_t _sizep = 0;
    set<string> _sizemodules;
    set<string> _sizepaths;
    for (vector<HLTData>::const_iterator it=hlt_vector.begin(); 
         it!=hlt_vector.end(); ++it){
        if (!it->hlt_value>0) continue;
        _sizemodules.insert(it->hlt_subname);
        _sizepaths.insert(it->hlt_name);
    }
    
    _sizemod = _sizemodules.size();
    _sizep = _sizepaths.size();

    TH2F* hlt_hist = new TH2F("hlt_hist22","",_sizep,0,_sizep,_sizemod,0,_sizemod);
    for (vector<HLTData>::const_iterator it=hlt_vector.begin(); 
         it!=hlt_vector.end(); ++it){
        if (!it->hlt_value>0) continue;
        //cout<<it->hlt_name<<"\t"<<it->hlt_value<<endl;
        hlt_hist->Fill(it->hlt_name.c_str(),it->hlt_subname.c_str(),it->hlt_value);
    }

    return hlt_hist;
}//fill a histo





//Check path times
//---------------------------------------------------------------------------
void check_pathtimes (Double_t _threshold, string incOrnormalTime)
{
//---------------------------------------------------------------------------
    Double_t threshold = _threshold;//in miliseconds
    cout<<"\n$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$Total number of events: "<<totalNevents<<endl;
    cout<<"Paths with "<<incOrnormalTime<<" more than "<<threshold<<" msec per event:\n"<<endl;
    vector<HLTData> pathtime_vector; 
    vector<string> mywildcards;
    mywildcards.push_back(incOrnormalTime);
    TH1* h;
    vector<TH1*> pathtime_histvect = get_vector_of_histos(mywildcards);
    for(vector<TH1*>::const_iterator vit = pathtime_histvect.begin(),
            vitend = pathtime_histvect.end();vit!=vitend;++vit){
        h = (*vit);
        //if (debugme) cout<<h<<endl;
        string longName = h->GetName();
        string shortName = longName.erase(0,9);
        Int_t nbins = h->GetNbinsX();
        Double_t binwidth =  h->GetBinWidth(1);
        //cout<<"binwidth = "<<binwidth<<endl;
        assert(binwidth>0);
        Int_t skipbins = Int_t(ceil(threshold/binwidth))+1;
        //cout<<"skipbins = "<<skipbins<<endl;
        Double_t highEntries = h->Integral(skipbins,nbins);
        //Double_t path_discriminator = 100*highEntries/histEntries;
        Double_t path_discriminator = highEntries;
      if(path_discriminator>0){
          cout<<"pathName = "<<longName<<"\tnumEvents = "<<path_discriminator<<endl;
      }
      pathtime_vector.push_back(HLTData(shortName,"",path_discriminator));
    }
    
    sort_hlt_vector(pathtime_vector);
    TH1F* pathtime_hist = fill_hlt_histogram(pathtime_vector);
    pathtime_hist->SetTitle(Form("Num events (out of %i) per path that took > %0.0f ms to run",threshold,totalNevents));
    save_histogram(pathtime_hist);
    
    
    return;
}//checkhisto




//Check average path time plot
//---------------------------------------------------------------------------
void check_average_pathtime (Double_t _threshold)
{
//---------------------------------------------------------------------------
    Double_t threshold = _threshold;//in miliseconds
    cout<<"\n$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$Total number of events: "<<totalNevents<<endl;
    cout<<"Paths with average time more than "<<threshold<<" msec:\n"<<endl;
    vector<HLTData> averagepath_vector; 
    vector<string> mywildcards;
    mywildcards.push_back("pathTimeSummary");
    TH1* h;
    vector<TH1*> pathtime_histvect = get_vector_of_histos(mywildcards);
    for(vector<TH1*>::const_iterator vit = pathtime_histvect.begin(),
            vitend = pathtime_histvect.end();vit!=vitend;++vit){
        h = (*vit);
        Int_t nbins = h->GetNbinsX();
        Double_t binwidth =  h->GetBinWidth(1);
        for (int bin=1 ; bin<=nbins ; ++bin){
            Double_t avpathtime = h->GetBinContent(bin);
            string pathlabel = h->GetXaxis()->GetBinLabel(bin);
            if (avpathtime>threshold){
                cout<<"PathName = "<<pathlabel<<"\tAverageTime = "<<avpathtime<<" msec"<<endl;
                averagepath_vector.push_back(HLTData(pathlabel,"",avpathtime));
            }
        }
    }
    
    //to do, plot averagepath_vector
    
    return;
}//checkhisto




//Check path times
//---------------------------------------------------------------------------
// void check_module_running_time (Double_t _threshold)
// {
// //---------------------------------------------------------------------------
//     if (debugme) cout<<"Checking module times"<<endl;
//     Double_t threshold = _threshold;//in miliseconds
//     vector<HLTData> moduletime_vector; 
//     TH1* h;
//     vector<TH1*> moduletime_histvect = get_vector_of_histos("moduleInPathScaledTimeSummary_");
//      for(vector<TH1*>::const_iterator vit = moduletime_histvect.begin(),
//              vitend = moduletime_histvect.end();vit!=vitend;++vit){
//          h = (*vit);
//          if (h->GetMaximum()<_threshold) continue;
//          string pathName = h->GetName();
//          pathName = pathName.erase(0,30);
//          Int_t nbins = h->GetNbinsX();
//          Double_t binwidth =  h->GetBinWidth(1);
         
//          for (int bin=1 ; bin<=nbins ; ++bin){
//              Double_t moduletime = h->GetBinContent(bin);
//              if (moduletime<_threshold) continue;
//              string modulelabel = h->GetXaxis()->GetBinLabel(bin);
//              moduletime_vector.push_back(HLTData(modulelabel,pathName,
//                                                  moduletime));
//          }
//      }

//      sort_hlt_vector(moduletime_vector);
//      TH2F* moduletime_hist = fill2_hlt_histogram(moduletime_vector);
//      moduletime_hist->SetTitle(Form("Average running time per module with time_threshold>%0.0f ms.  Third dimension shows number of ocurrences.",threshold));
//      save_histogram(moduletime_hist);
    
    
//     return;
// }//checkhisto



//Check path times
//---------------------------------------------------------------------------
void check_module_per_path_running_time (Double_t _threshold)
{
//---------------------------------------------------------------------------

    Double_t threshold = _threshold;//in miliseconds
    cout<<"\n$$$$$$$$$$$$$$$$$$$$$$$$Total number of events: "<<totalNevents<<endl;
    cout<<"Modules (in corresponding paths) that took more than "<<threshold<<" msec to run per event:"<<endl;
    vector<HLTData> modulepath_time_vector;
    vector<string> mywildcards;
    mywildcards.push_back("moduleInPathScaledTime_HLT_");
    mywildcards.push_back("moduleInPathScaledTime_AlCa_");
    TH1* h;
    vector<TH1*> modulepath_time_histvect = get_vector_of_histos(mywildcards);
     for(vector<TH1*>::const_iterator vit = modulepath_time_histvect.begin(),
             vitend = modulepath_time_histvect.end();vit!=vitend;++vit){
         h = (*vit);
         string longName = h->GetName();
         string shortName = longName.erase(0,27);
         string pathName = "";
         string modName = "";
         splitstring(shortName,"_hlt",pathName,modName);
         modName = "hlt"+modName;
         //cout<<h->GetName()<<"\t\t"<<pathName<<"\t"<<modName<<endl;
         Int_t nbins = h->GetNbinsX();
         Double_t binwidth =  h->GetBinWidth(1);
         assert(binwidth>0);
         Int_t skipbins = Int_t(ceil(threshold/binwidth))+1;
         Double_t highEntries = h->Integral(skipbins,nbins);
         //Double_t path_discriminator = 100*highEntries/histEntries;
         Double_t time_discriminator = highEntries;
         if(time_discriminator>0){
             cout<<"shortPathName = "<<pathName<<"\tmodName = "<<modName<<"\tnumEvents = "<<time_discriminator<<endl;
         }
         modulepath_time_vector.push_back(HLTData(pathName,modName,time_discriminator));
         
         
         
     }
     
     sort_hlt_vector(modulepath_time_vector);
     TH2F* modulepath_time_hist = fill22_hlt_histogram(modulepath_time_vector);
     modulepath_time_hist->SetTitle(Form("Modules (in corresponding paths) that take >%0.0f ms to run per event.  Third dimension shows number of events out of %i.",threshold,totalNevents));
     save_histogram(modulepath_time_hist);
    
    return;
}//checkhisto



//The main() function
//-----------------------------------------------------------------------------
int main()
{
//-----------------------------------------------------------------------------

    //look at pathtimes (could be incPathTime_ or just pathTime_)
    //for the moment, if you need plots, don't do incPathTime and pathTime
    //at the same time.  Printing txt file is fine.

    //check_pathtimes(40,"pathTime_");
    check_average_pathtime(10);
    check_average_pathtime(40);
    check_pathtimes(100,"pathTime_");
    check_pathtimes(5,"incPathTime_");
    check_module_per_path_running_time(40);
    check_module_per_path_running_time(100);
    //check_pathtimes(10,"incPathTime_");
    //check_pathtimes(150);
    //check_pathtimes(200);
    //check_module_running_time(10);
    //check_module_running_time(40);
    //check_module_per_path_running_time(20);
    //check_module_per_path_running_time(30);



    
    
    return 0;

}//main()



