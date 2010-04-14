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
#include <TMath.h>
#include <TROOT.h>
#include <TFile.h>
#include <TMultiGraph.h>
#include <TLegend.h>

using namespace std;


struct data{
    string algo;
    float cut_val1,cut_val2,cut_val3,cut_val4, ntotevt, npassevt;
  data(string ALGO, float CUT_VAL1 ,float CUT_VAL2,float CUT_VAL3,
       float CUT_VAL4,float NTOTEVT, float NPASSEVT):
    algo(ALGO), cut_val1(CUT_VAL1),cut_val2(CUT_VAL2),
    cut_val3(CUT_VAL3),cut_val4(CUT_VAL4), 
    ntotevt(NTOTEVT), npassevt(NPASSEVT){}
};

struct datasignif{
  string algo; 
  float signifval,dsignif,theval1, theval2,theval3,theval4;
  datasignif(string ALGO, float SIGNIFVAL, float DSIGNIF,
	     float THEVAL1,float THEVAL2,
	     float THEVAL3,float THEVAL4):
    algo(ALGO), signifval(SIGNIFVAL), dsignif(DSIGNIF),
    theval1(THEVAL1),theval2(THEVAL2),
    theval3(THEVAL3),theval4(THEVAL4){}
  
};

//location of the files containing the *.dat files with the information
//for optimization
const string dirfiles = "forOptimization/";
//fix this!! we should be able to read this values directly
const float o_ValDeltaRIso[9] = {0.2,0.25,0.3,0.35,0.4,0.45,0.5,0.55,0.6};
const float o_DphiJetVeto[5] = { TMath::Pi() - 0.1, 
                                       TMath::Pi() - 0.2, 
                                       TMath::Pi() - 0.3, 
                                       TMath::Pi() - 0.4,
                                       TMath::Pi() - 0.5};
//switch for debugging
const bool debugme = false;

//number of algorithms
//const int Num_trkAlgos = 3;
//const string algoType[Num_trkAlgos] = {"glb", "trk", "tev"};
//const int Num_trkAlgos = 1;
//const string algoType[Num_trkAlgos] = {"glb"};

// Calculate efficiencies
//------------------------------------------------------------------------
void getEff(float & eff, float & deff,float Num,float Denom)
{
//------------------------------------------------------------------------
  eff = Num/Denom;
  deff = TMath::Sqrt(eff * (1-eff)/Denom);
}//---------------getEff


//------------------------------------------------------------------------
void punzi(float & signif, float& dsignif, float sigeff, float dsigeff,
	   float backg, float a)
{
//------------------------------------------------------------------------
  if(backg>=0){
    signif = sigeff/((a/2)+sqrt(backg));
    dsignif = (signif/sigeff)*(dsigeff-(signif/2));
  }

}//-----punzi()



//------------------------------------------------------------
void parse_file(const string& filename, const string& backname,
		vector<data>& mydata)
{
//------------------------------------------------------------
    if(debugme) cout<<"Parsing files"<<endl;

    //this only affects the "all" options for the backgrounds.
    vector< string> typesOfBack;
    if (backname == "all"){
      typesOfBack.push_back("QCD");
      typesOfBack.push_back("Z");
      //typesOfBack.push_back("W");
      typesOfBack.push_back("Wback");
      typesOfBack.push_back("Top");
    }
    if (backname == "W+all"){
      typesOfBack.push_back("QCD");
      typesOfBack.push_back("Z");
      typesOfBack.push_back("W");
      typesOfBack.push_back("Top");
    }
    else if(backname == "QCD+Top"){
       typesOfBack.push_back("QCD");
       typesOfBack.push_back("Top");
    }
    //background is just one of the processess:
    //string.find returns the max possible number if not found.
    if (filename.find("_")!=string::npos){
        ifstream myfile(filename.c_str());
        string _algo;
        float _cv1,_cv2,_cv3,_cv4,_ntotevt,_npassevt;
        while(myfile>>_algo>>_cv1>>_cv2>>_cv3>>_cv4>>_ntotevt>>_npassevt){
            mydata.push_back(data(_algo,_cv1,_cv2,_cv3,_cv4,_ntotevt,_npassevt));
        }
    }
    //background is the sum of all backgrounds
    else{
        for (unsigned int j = 0; j<typesOfBack.size(); ++j){
            //here filename is the var
            string _filename = dirfiles+typesOfBack.at(j)+"_"+filename+".dat";
            if(debugme) cout<<"Parsing file "<<_filename<<endl;
            ifstream myfiles(_filename.c_str());
            string talgo; 
	    float tcv1,tcv2,tcv3,tcv4,tntot,tnpass;
            if (j == 0){
	      while(myfiles>>talgo>>tcv1>>tcv2>>tcv3>>tcv4>>tntot>>tnpass){
                    mydata.push_back(data(talgo,tcv1,tcv2,tcv3,tcv4,
					  tntot,tnpass));
                }
            }
            else{
                int mycounter = 0;
                while(myfiles>>talgo>>tcv1>>tcv2>>tcv3
		      >>tcv4>>tntot>>tnpass){
                    mydata.at(mycounter).ntotevt+=tntot;
                    mydata.at(mycounter).npassevt+=tnpass;
                    ++mycounter;
                }
            }
        }
    }

    
}//-------



//------------------------------------------------------------
void fill_arrays_and_graph(TGraphErrors* myg[],
			   const int myPar1dim, const int myPar2dim,
			   const vector<datasignif>& signif)
//------------------------------------------------------------
{
    if(debugme) cout<<"Filling arrays and graphing"<<endl;
    //currently, don't care about algos, just do it
    //for gbl

    //just one algo for now
    string thealgo = "gbl";
    //dimension of X-axis for graphs
    const int myXdim = myPar1dim;
    //dimension of # of curves due to par2
    const int myNcurves = myPar2dim;
    if (debugme) cout<<"myXdim = "<<myXdim<<endl;
    if (debugme) cout<<"myNcurves = "<<myNcurves<<endl;


    int myjump = 0;

    for (int p2=0;p2<myNcurves;++p2){
      if(debugme) cout<<"define arrays for plots"<<endl;
      float significance[myXdim];
      float cut_range[myXdim];
      float signif_err[myXdim];
      float cut_range_err[myXdim];
      

      if(debugme) cout<<"Fill arrays"<<endl;
      int mycounter = 0;
      for (int i=myjump;i<(myjump+myXdim);++i){
	if(signif.at(i).algo != thealgo) return;
	if(debugme) cout<<"thealgo = "<<thealgo<<endl;
	significance[mycounter] = signif.at(i).signifval;
	cut_range[mycounter] = signif.at(i).theval1;
	signif_err[mycounter] = signif.at(i).dsignif;
	//signif_err[mycounter] = 0;
	cut_range_err[mycounter] = 0;
	if (debugme) cout<<"significance["<<mycounter<<"] = "
			 <<significance[mycounter]<<
		       "\t signif_err["<<mycounter<<"] = "
			 <<signif_err[mycounter]
			 <<"\t cut_range["<<mycounter
			 <<"] = "<<cut_range[mycounter]<<endl;
	
	++mycounter;
      }
      
      
      
      //set colors and markers for curves
      int colorme = p2+1;
      int markerme = p2+22;
      
      myg[p2] = new TGraphErrors(myXdim,cut_range,significance,
				 cut_range_err,signif_err);
      
      
      
      myg[p2]->SetMarkerColor(colorme);
      myg[p2]->SetMarkerStyle(markerme);
      myg[p2]->SetMarkerSize(1.3);
      myg[p2]->SetLineColor(colorme);
      //myg[p2]->SetLineWidth(2);

      myjump+=myXdim;

    }//loop over par2vals    

}//fill_arrays_for_plot



//--------------------------------------------------------------------
void calc_signif_and_cutrange(string file_sig,string file_bkg,
			      const string& backname,
                              vector<datasignif>& signif,
			      const string& myestimator)
{
//---------------------------------------------------------------------
    if(debugme) cout<<"Getting vectors according to input files"<<endl;

    vector<data> mysignal;
    vector<data> mybackground;
    vector<float> sigeff;
    vector<float> bkg;
    
    parse_file(file_sig,backname,mysignal);
    parse_file(file_bkg,backname,mybackground);
    
    assert(mysignal.size() == mybackground.size());
    const int mynumpoints = mysignal.size();

    //fill significance vector
    for(int j = 0;j<mynumpoints;++j){
        //calculate signal efficiency
        float sgl = mysignal.at(j).npassevt;
        float sigeff = 0;
	float dsigeff = 0;
        if(mysignal.at(j).ntotevt!=0) {
	  getEff(sigeff,dsigeff,mysignal.at(j).npassevt,
		 mysignal.at(j).ntotevt);
	}
        //calculate significance
        float my_a = 5.; //number of sigmas
        float backg = mybackground.at(j).npassevt;
	float _signif = 0;
	float _dsignif = 0;
	//usual significance-like estimator
	if((myestimator == "sig1") && ((backg+sgl)>0)){ 
	  _signif = sgl/sqrt(backg+sgl);
	}
	//another significance-like estimator
	else if ((myestimator == "sig2") && (backg>0)){
	  _signif = sgl/sqrt(backg);
	}
	//the punzi estimator
	//http://arxiv.org/pdf/physics/0308063v2
	else if(myestimator == "punzi"){
	  punzi(_signif,_dsignif,sigeff,dsigeff,backg,my_a);
	}
	else if(myestimator == "eff"){
	  _signif = sigeff;
	  _dsignif = dsigeff;
	}
	if(debugme) cout<<"cut_val1 = "<<mysignal.at(j).cut_val1<<"\tsgl = "<<sgl<<"\t backg = "<<backg<<"\t signif = "<<_signif<<endl;
        signif.push_back(datasignif(mysignal.at(j).algo,
				    _signif,_dsignif,
				    mysignal.at(j).cut_val1,
				    mysignal.at(j).cut_val2,
				    mysignal.at(j).cut_val3,
				    mysignal.at(j).cut_val4));
    }

 }//---------------------------




//--------------------------------------------------------------------
void plotCutOpt_TwoPars(const string& var, const string& signame,
			const string& backname, TGraphErrors* myg[],
			const int myPar1dim, const int myPar2dim,
			const vector<datasignif>& signif, 
			const string& myestimator)
//--------------------------------------------------------------------
{

  fill_arrays_and_graph(myg,myPar1dim,myPar2dim,signif);

  //fill a multigraph
  if(debugme) cout<<"Form convenient information strings"<<endl;
  string mytitle = "Optimization "+var+", sig: "+signame+", bkg: "+backname;
  string myYtitle = "";
  if (myestimator == "punzi"){
    myYtitle = "#frac{#epsilon_{S}}{a/2 + #sqrt{B}}";
  }
  else if (myestimator == "sig1"){
    myYtitle = "S/ #sqrt{S+B}";
  }
  else if (myestimator == "sig2"){
    myYtitle = "#frac{S}{#sqrt{B}}";
  }
  else if (myestimator == "eff"){
    myYtitle = "signal eff";
    mytitle = "Optimization "+var+", sig: "+signame;
  }
  
  string myXtitle = "";
  if (var == "1mu") myXtitle = "muon p^{trk}_{T} (GeV)";
  if (var == "1muiso" || var == "1mujetiso") myXtitle = "#Sigma pT (GeV)";
  if (var == "1mujet" || var == "1muisojet") myXtitle = "jet ET (GeV)";
  if (var == "1mujetisoqual" || var == "1muisojetqual") myXtitle = "muon p^{trk}_{T} (GeV)";
  

  if(debugme) cout<<"Creating multigraph"<<endl;
  TMultiGraph* mg = new TMultiGraph();
  for (int j=0;j<myPar2dim;++j){
    mg->Add(myg[j]);
  }

  TCanvas* c1 = new TCanvas("c1","",1200,800);
  c1->cd();


  if(debugme) cout<<"Drawing..."<<endl;
  mg->Draw("APC");
  mg->SetTitle(mytitle.c_str());
  mg->GetXaxis()->SetTitle(myXtitle.c_str());
  mg->GetYaxis()->SetTitle(myYtitle.c_str());
  mg->GetYaxis()->SetTitleOffset(1.5);
  mg->GetYaxis()->SetTitleSize(0.03);
  mg->GetYaxis()->SetLabelSize(0.025);
  mg->GetXaxis()->SetLabelSize(0.025);

  
  if(debugme) cout<<"Creating legend"<<endl;

  if(var != "1mu"){
    TLegend* leg = new TLegend(0.8,0.3,0.6,0.7);
    for (int j = 0; j<myPar2dim;++j){
      char buffer[50];
      if (var == "1muiso" || var == "1mujetiso") sprintf(buffer,"%.3f",o_ValDeltaRIso[j]);
      else if (var == "1mujet" || var == "1muisojet") sprintf(buffer,"%.3f",o_DphiJetVeto[j]);
      else if (var == "1mujetisoqual" || var == "1muisojetqual") sprintf(buffer,"%.3f",signif.at(j*myPar1dim).theval2);
      else { cout<<"Something went wrong while plotting... quit...."<<endl;abort();}
      leg->AddEntry(myg[j],buffer,"lp");
    }
    leg->SetTextSize(0.04);
    leg->SetLineColor(0);
    leg->SetLineStyle(0);
    leg->SetLineWidth(0);
    leg->SetFillColor(0);
    leg->Draw("same");
  }
  string gifname = "~/www/analyses/wprime/optim/"+myestimator+"_"+var+"_"+signame+"_"+backname+".gif";
  //c1->Print(gifnameweb.c_str());
  c1->Print(gifname.c_str());
  
  
}//------------plotTwoPars




//--------------------------------------------------------------------
void get_dimensions(const vector<datasignif>& signif,int& mydimPar1,
		  int& mydimPar2)
{
//--------------------------------------------------------------------

  int tempalgosize = 0;
  int temp2ndsize = 0;
  
  //get size of the algo in 0th column in dat file
  string firstalgo = signif.at(0).algo;
  for (unsigned int j = 0;j<signif.size();++j){
    if (signif.at(j).algo == firstalgo) ++tempalgosize;
    else {break;}
  }
  //get the size of 2nd column in dat file
  //and calculate dimensions
  float first2ndpar = signif.at(0).theval2;
  for (unsigned int j = 0;j<signif.size();++j){
    if (signif.at(j).theval2 == first2ndpar) ++temp2ndsize;
    else {break;}
  }
  //get the dimensions of par1 and par2(if available) by
  //checking the sizes calculated earlier
  if (temp2ndsize > tempalgosize){
    mydimPar1 = tempalgosize;
    //do nothing for mydimPar2 because its 1 already
  }
  else{
    mydimPar1 = temp2ndsize;
    mydimPar2 = int(tempalgosize/temp2ndsize);
  }

  

}//---------get_dimension


//--------------------------------------------------------------------
void plotOptim(const string& var, const string& signame,
	       const string& backname,
	       const string myestimator="punzi")
{
//--------------------------------------------------------------------

    if(debugme) cout<<"Starting plotOptim function..."<<endl;
    

    //gROOT->SetStyle("Plain");
//    gROOT->ProcessLine(".L tdrstyle.C");
//    gROOT->ProcessLine("setTDRStyle()");
  
    //gStyle->SetPadLeftMargin(0.13);    


    vector<datasignif> signif;

    string file1 = dirfiles+signame+"_"+var+".dat";
    string file2 = "";
    //deal with the background being the sum of all of them
    //the file name carries the variable type in that case
    if(backname == "all" || backname.find("+")!=string::npos) file2 = var; 
    //or just a single type of background
    else file2 = dirfiles+backname+"_"+var+".dat";
    //get vectors according to input files
    calc_signif_and_cutrange(file1,file2,backname,signif,myestimator);
   
   
    //prepare TGraphsErrors for different options
    //for now just concentrate on one algo type
    int mydimPar1 = 1;
    int mydimPar2 = 1;
    get_dimensions(signif,mydimPar1,mydimPar2);
    const int myArrsize = mydimPar2;
    TGraphErrors* myg[myArrsize];
    plotCutOpt_TwoPars(var, signame, backname,myg,
		       mydimPar1,mydimPar2,signif,myestimator);
   
    



}//-----------------------------------
  



