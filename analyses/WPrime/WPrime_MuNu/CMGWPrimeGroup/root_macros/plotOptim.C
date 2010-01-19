#include <stdio.h>
#include <algorithm>
#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <assert.h>
#include <vector>
#include <string>
#include <math.h>
#include <TAxis.h>
#include <TCanvas.h>
#include <TGraphErrors.h>
#include <TPaveStats.h>
#include <TStyle.h>
#include <TLatex.h>
#include <TLegend.h>
#include <TROOT.h>
#include <TFile.h>
#include <TMultiGraph.h>
#include <TLegend.h>

using namespace std;


struct data{
    string algo;
    float cut_val, ntotevt, npassevt;
    data(string ALGO, float CUT_VAL, float NTOTEVT, float NPASSEVT):
        algo(ALGO), cut_val(CUT_VAL), ntotevt(NTOTEVT), npassevt(NPASSEVT){}
};

struct datasignif{
    string algo; float thevalue;
    datasignif(string ALGO, float THEVALUE):
        algo(ALGO), thevalue(THEVALUE){}

};

//location of the files containing the *.dat files with the information
//for optimization
const string dirfiles = "forOptimization/";

//switch for debugging
const bool debugme = false;

//number of algorithms
const int Num_trkAlgos = 3;
const string algoType[Num_trkAlgos] = {"glb", "trk", "tev"};


//------------------------------------------------------------
void parse_file(const string& filename, vector<data>& mydata)
{
//------------------------------------------------------------
    if(debugme) cout<<"Parsing files"<<endl;

    vector< string> typesOfBack;
    typesOfBack.push_back("QCD");
    typesOfBack.push_back("Z");
    typesOfBack.push_back("W");
    typesOfBack.push_back("Top");

    //background is just one of the processess:
    if (filename.find("_")<1000){
        ifstream myfile(filename.c_str());
        string _algo;
        float _cut_val,_ntotevt,_npassevt;
        while(myfile>>_algo>>_cut_val>>_ntotevt>>_npassevt){
            mydata.push_back(data(_algo,_cut_val,_ntotevt,_npassevt));
        }
    }
    //background is the sum of all backgrounds
    else{
        for (unsigned int j = 0; j<typesOfBack.size(); ++j){
            //here filename is the var
            string _filename = dirfiles+typesOfBack.at(j)+"_"+filename+".dat";
            if(debugme) cout<<"Parsing file "<<_filename<<endl;
            ifstream myfiles(_filename.c_str());
            string talgo; float tcutval,tntot,tnpass;
            if (j == 0){
                while(myfiles>>talgo>>tcutval>>tntot>>tnpass){
                    mydata.push_back(data(talgo,tcutval,tntot,tnpass));
                }
            }
            else{
                int mycounter = 0;
                while(myfiles>>talgo>>tcutval>>tntot>>tnpass){
                    mydata.at(mycounter).ntotevt+=tntot;
                    mydata.at(mycounter).npassevt+=tnpass;
                    ++mycounter;
                }
            }
        }
    }

    
}//-------

//------------------------------------------------------------
int get_algo_dim(const vector<datasignif>& signif, 
                          const string& algo)
//------------------------------------------------------------
{
    if(debugme) cout<<"\t Getting dimensions for each algo"<<endl;
    int temp_dim = 0;
    for (unsigned int j = 0;j<signif.size();++j)
        if (signif.at(j).algo == algo) ++temp_dim;

    return temp_dim;
    

}//get_dim_for_each_algo



//------------------------------------------------------------
void fill_arrays_and_graph(TGraphErrors* myg[], 
                           const vector<datasignif>& signif,
                           const vector<datasignif>& cuts)
//------------------------------------------------------------
{
    if(debugme) cout<<"Filling arrays and graphing"<<endl;
     //get the number of points for each muon algorithm
    //current algos: glb, trk, tev
  
    //loop over algos
    for (int na = 0;na<Num_trkAlgos;++na){
        string thealgo = algoType[na];
        const int mydim = get_algo_dim(signif,thealgo);
        if (debugme) cout<<"mydim = "<<mydim<<endl;
                
        //define arrays for plots
        float significance[mydim];
        float cut_range[mydim];
        float signif_err[mydim];
        float cut_range_err[mydim];
        
        //Fill arrays
        int mycounter = 0;
        for (unsigned int i=0;i<signif.size();++i){
            if(signif.at(i).algo == thealgo){
                if(debugme) cout<<"thealgo = "<<thealgo<<endl;
                //make sure we don't go out of limits
                assert(mycounter <= mydim);
                significance[mycounter] = signif.at(i).thevalue;
                cut_range[mycounter] = cuts.at(i).thevalue;
                signif_err[mycounter] = 0;
                cut_range_err[mycounter] = 0;
                if (debugme) cout<<"significance["<<mycounter<<"] = "
                                 <<significance[mycounter]
                                 <<"\t cut_range["<<mycounter
                                 <<"] = "<<cut_range[mycounter]<<endl;
                
                ++mycounter;
            }
        }
        
        
        
        //set colors and markers for curves depending on the algo
        int colorme = na*3+1;
        int markerme = na+22;
        
        myg[na] = new TGraphErrors(mydim,cut_range,significance,
                               cut_range_err,signif_err);
        
        
    
        myg[na]->SetMarkerColor(colorme);
        myg[na]->SetMarkerStyle(markerme);
        myg[na]->SetMarkerSize(1.3);
        myg[na]->SetLineColor(colorme);
        myg[na]->SetLineWidth(2);

    }//loop over algos

    
    

}//fill_arrays_for_plot





//--------------------------------------------------------------------
void calc_signif_and_cutrange(string file_sig, string file_bkg,
                              vector<datasignif>& signif,
                              vector<datasignif>& cuts)
{
//---------------------------------------------------------------------
    if(debugme) cout<<"Getting vectors according to input files"<<endl;

    vector<data> mysignal;
    vector<data> mybackground;
    vector<float> sigeff;
    vector<float> bkg;
    
    parse_file(file_sig,mysignal);
    parse_file(file_bkg,mybackground);
    
    assert(mysignal.size() == mybackground.size());
    const int mynumpoints = mysignal.size();

    //fill significance and cuts vector
    for(int j = 0;j<mynumpoints;++j){
        //fill cuts vector
        cuts.push_back(datasignif(mysignal.at(j).algo,mysignal.at(j).cut_val));
        //calculate signal efficiency
        float sgl = mysignal.at(j).npassevt;
        float sigeff = 0;
        if(mysignal.at(j).ntotevt!=0) 
            sigeff = mysignal.at(j).npassevt/mysignal.at(j).ntotevt;
        //calculate significance
        //float my_a = 5.; //number of sigmas
        float backg = mybackground.at(j).npassevt;
        //float _signif = sigeff/((my_a/2)+sqrt(backg));
        float _signif = sgl/sqrt(backg+sgl);
        string _algo = mysignal.at(j).algo;
        signif.push_back(datasignif(_algo,_signif));
    }

 }//---------------------------




//--------------------------------------------------------------------
void plotOptim(const string& var, const string& signame,
                       const string& backname)
{
//--------------------------------------------------------------------

    if(debugme) cout<<"Starting plotOptim function"<<endl;
    
    //gROOT->SetStyle("Plain");
//    gROOT->ProcessLine(".L tdrstyle.C");
//    gROOT->ProcessLine("setTDRStyle()");
  
    //gStyle->SetPadLeftMargin(0.13);    

    vector<datasignif> signif;
    vector<datasignif> cuts;
    string file1 = dirfiles+signame+"_"+var+".dat";
    string file2 = "";
    //deal with the background being the sum of all of them
    //the file name carries the variable type in that case
    if(backname == "all") file2 = var; 
    //or just a single type of background
    else file2 = dirfiles+backname+"_"+var+".dat";
    //get vectors according to input files
    calc_signif_and_cutrange(file1,file2,signif,cuts);
   
   
    //prepare TGraphsErrors for the thre algos we have
    TGraphErrors*  myg[Num_trkAlgos];
    fill_arrays_and_graph(myg,signif,cuts);

    //fill a multigraph
    if(debugme) cout<<"Form convenient information strings"<<endl;
    string mytitle = "Optimization "+var+", sig: "+signame+", bkg: "+backname;
    //string myYtitle = "#frac{#epsilon_{S}}{a/2 + #sqrt{B}}";
    //string myYtitle = "#frac{S}{#sqrt{B}}";
    string myYtitle = "S/ #sqrt{S+B}";
    string myXtitle = var+" [GeV]";
   
    if(debugme) cout<<"Creating multigraph"<<endl;
    TMultiGraph* mg = new TMultiGraph();
    mg->Add(myg[0]);
    mg->Add(myg[1]);
    mg->Add(myg[2]);
    mg->Draw("APC");
    mg->SetTitle(mytitle.c_str());
    mg->GetXaxis()->SetTitle(myXtitle.c_str());
    mg->GetYaxis()->SetTitle(myYtitle.c_str());
    mg->GetYaxis()->SetTitleOffset(1.4);

     if(debugme) cout<<"Creating legend"<<endl;
     TLegend* leg = new TLegend(0.8,0.6,0.6,0.8);
     for (int j = 0; j<Num_trkAlgos;++j){
         leg->AddEntry(myg[j],algoType[j].c_str(),"lp");
     }
     leg->SetTextSize(0.04);
     leg->SetLineColor(0);
     leg->SetLineStyle(0);
     leg->SetLineWidth(0);
     leg->SetFillColor(0);
     leg->Draw("same");

    //string gifname = "plots/"+var+"_"+signame+"_"+backname+".gif";
    //  c1.Print(gifnameweb.c_str());
    //c1->Print(gifname.c_str());

}//-----------------------------------
  



