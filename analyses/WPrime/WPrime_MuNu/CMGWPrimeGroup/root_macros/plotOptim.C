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
#include <TPad.h>
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
#include <TKey.h>
#include <TH1F.h>
#include <TString.h>

using namespace std;



//data struct for histograms in corresponding 
//sample folder
struct datahist{
    string dir;
    vector<TH1F*> h1f; 
    datahist(string DIR, vector<TH1F*> H1F):
        dir(DIR), h1f(H1F){}
};


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
//const string dirfiles = "forOptimization_1mujet_0To120_6step_WsigJetVeto60/";
//fix this!! we should be able to read this values directly
const float o_ValDeltaRIso[9] = {0.2,0.25,0.3,0.35,0.4,0.45,0.5,0.55,0.6};
const float o_DphiJetVeto[5] = { TMath::Pi() - 0.1, 
                                 TMath::Pi() - 0.2, 
                                 TMath::Pi() - 0.3, 
                                 TMath::Pi() - 0.4,
                                 TMath::Pi() - 0.5};
//punzi estimator value for the number of sigmas
const float punzi_a_value = 3;
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
    //background is the sum of all backgrounds or combined backgrounds
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
                           const vector<datasignif>& signif,
                           const string& algo,const vector<string>& myAlgoName )
//------------------------------------------------------------
{
    if(debugme) cout<<"Filling arrays and graphing"<<endl;
    //currently, don't care about algos, just do it
    //for gbl

    //just one algo for now
    const string thealgo = algo;
    if (debugme) cout<<"fill_arrays: algo = "<<thealgo<<endl;
    //dimension of X-axis for graphs
    const int myXdim = myPar1dim;
    //dimension of # of curves due to par2
    const int myNcurves = myPar2dim;
    if (debugme) cout<<"myXdim = "<<myXdim<<endl;
    if (debugme) cout<<"myNcurves = "<<myNcurves<<endl;
    //determine where to jump for the algo
    int indexAlgoInUse = -1;
    for(unsigned int j = 0;j<myAlgoName.size();++j){
        if(algo == myAlgoName.at(j)) indexAlgoInUse = j;
    }
    //jump to the correct algorithm
    int firstcolumnsize = myPar1dim*myPar2dim;
    int jumpAlgo = indexAlgoInUse*firstcolumnsize;
    //to jump to the next set of par2's
    int jumpPar2 = 0;


    //loop over par2 values
    for (int p2=0;p2<myNcurves;++p2){
        if(debugme) cout<<"define arrays for plots"<<endl;
        float significance[myXdim];
        float cut_range[myXdim];
        float signif_err[myXdim];
        float cut_range_err[myXdim];
      
      
        if(debugme) cout<<"Fill arrays"<<endl;
        int mycounter = 0;
        for (int i=jumpAlgo+jumpPar2;i<(jumpAlgo+jumpPar2+myXdim);++i){
            if(signif.at(i).algo != thealgo) continue;
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
    

        jumpPar2+=myXdim;
        
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
        float my_a = punzi_a_value; //number of sigmas
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
                        const string& myestimator,const string& algo,
                        const vector<string>& myAlgoName)
//--------------------------------------------------------------------
{

    fill_arrays_and_graph(myg,myPar1dim,myPar2dim,signif,algo,myAlgoName);

  //fill a multigraph
  if(debugme) cout<<"Form convenient information strings"<<endl;
  string mytitle = "Optimization "+var+", sig: "+signame+", bkg: "+
      backname+", algo:"+algo;
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
  //if (var == "1muiso" || var == "1mujetiso") myXtitle = "#Sigma pT (GeV)";
  if (var == "1muiso" || var == "1mujetiso") myXtitle = "(#Sigma pT+ Et_{ECAL}+ Et_{HCAL})/p^{trk}_{T} (GeV)";
  if (var == "1mujet" || var == "1muisojet") myXtitle = "jet ET (GeV)";
  if (var == "1mujetisoqual" || var == "1muisojetqual") myXtitle = "#chi^{2}/nof";
  

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
  //c1->Print(gifname.c_str());
  
  
}//------------plotTwoPars




//--------------------------------------------------------------------
void get_dimensions(const vector<datasignif>& signif,
                    vector<string>& myAlgoName,
                    int& mydimPar1,int& mydimPar2)
{
//--------------------------------------------------------------------

    //sizes of repetitive elements in the first and third columns
    int firstcolsize = 0;
    int thirdcolsize = 0;
    
    //to get the name and dim of algos
    string firstalgo = signif.at(0).algo;
    string tempalgo = firstalgo;
    myAlgoName.push_back(signif.at(0).algo);
    //get size of the algo in 0th column in dat file
    for (unsigned int j = 0;j<signif.size();++j){
        if (signif.at(j).algo == firstalgo) ++firstcolsize;
        else {break;}
    }
    //get the algo names
    for (unsigned int j = 0;j<signif.size();++j){
        if (signif.at(j).algo != tempalgo) {
            myAlgoName.push_back(signif.at(j).algo);
            tempalgo = signif.at(j).algo;
        }
    }

  //get the size of 3rd column in dat file
  //and calculate dimensions
  float first2ndpar = signif.at(0).theval2;
  for (unsigned int j = 0;j<signif.size();++j){
    if (signif.at(j).theval2 == first2ndpar) ++thirdcolsize;
    else {break;}
  }
  //get the dimensions of algom, par1 and par2(if available) by
  //checking the sizes calculated earlier
  if (thirdcolsize >= firstcolsize){
      //the above condition happens when the 3rd column are all -9999, i.e,
      //there is only one parameter.
    mydimPar1 = firstcolsize;
    //do nothing for mydimPar2 because its 1 already
  }
  else{
      mydimPar1 = thirdcolsize;
      mydimPar2 = int(firstcolsize/thirdcolsize);
  }


  

}//---------get_dimension







//---------------------------------------------------------------------------
void plotEffRatios(TH1F* effh[],const int& effhdim,
                   TH1F* myrefhist,
                   const string& algon)
{
//---------------------------------------------------------------------------
    if(debugme) cout<<"processing plotEffRatios...."<<endl;
   
    gStyle->SetOptStat(00000);
    gStyle->SetErrorX(0);
    
    

    TCanvas* c1 = new TCanvas("c1","",800,1200);
    c1->Divide(1,2);
    //plot ref histogram
    c1->cd(1);
    gPad->SetLogy();
    myrefhist->Draw();
    myrefhist->SetTitle("Total Distribution with reference cut");
    myrefhist->GetXaxis()->SetTitle("Muon p_{T} (GeV/c)");
    

    //now plot the efficiency histos
    c1->cd(2);
    int firsthist = 0; //to track the first plotted hist
    //legends
    TLegend * lg = new TLegend(0.59, 0.67, 0.89, 0.89);
    lg->SetTextSize(0.03);
    lg->SetBorderSize(0);
    lg->SetFillColor(0);
    TLine *line1 = new TLine(effh[0]->fXaxis.GetXmin(), 1, 
                             effh[0]->fXaxis.GetXmax(), 1); 
    lg->AddEntry(line1,myrefhist->GetName(),"L");

    //loop for drawing
    for (int j = 0; j<effhdim;++j){
        string hname = effh[j]->GetName();
        if(debugme) cout<<"hname = "<<hname<<endl;
        if(hname.find(algon) != string::npos){
            if(firsthist == 0){
                effh[j]->SetTitle("Relative Distribution change with respect to the reference");
                effh[j]->GetXaxis()->SetTitle("Muon p_{T} (GeV/c)");
                effh[j]->SetLineWidth(2);
                effh[j]->SetMinimum(0.8);
                effh[j]->SetMaximum(1.2);
                effh[j]->SetLineColor(j+2);
                effh[j]->Draw();
                ++firsthist;
            }
            else{
                effh[j]->SetLineWidth(2);
                effh[j]->SetLineColor(j+2);
                effh[j]->Draw("same");
            }

            //legends
            lg->AddEntry(effh[j],hname.c_str(),"L");
            
            
        }//find(algon)
    }//effhdim
    
    lg->Draw();
    
    line1->Draw();

      
    


}//-------plotEffRatios






//---------------------------------------------------------------------------
void makeEffRatioHists(const vector<TH1F*> totalhists,
                       int& refhidx,
                      TH1F* effh[],const string& algon, 
                      const float& par1, const float& par2)
{
//---------------------------------------------------------------------------
    if(debugme) cout<<"processing makeEffRatioHist"<<endl;
    char refhistname[1000];
    sprintf(refhistname,"%s_%.3f_%.3f",algon.c_str(),par1,par2);
    int refhistidx = -1;
    //initialize the eff histos
    int ed =0; //eff dim counter
    for(unsigned int j = 0;j<totalhists.size();++j){
        string hname = totalhists.at(j)->GetName();
        string newhname = "eff_"+hname;
        //cout<<"newhname = "<<newhname<<endl;
        if(hname != refhistname) {
            effh[ed] = (TH1F*)totalhists.at(j)->Clone();
            effh[ed]->SetName(newhname.c_str());
            effh[ed]->Reset();
            if(debugme) cout<<"effh["<<ed<<"] = "<<effh[ed]->GetName()<<endl;
            ++ed;
        }
        else {
            refhistidx = j; 
            if(true) cout<<"reference hist index is = "<<refhistidx<<endl;
        }
    }//init eff histos
    
    //check that we found a reference histogram
    if (refhistidx <0) {
        cout<<"No reference histogram found. Something went wrong, please check"<<endl;
    }
    

    //fill the eff histograms
    ed = 0; //reset the eff dim counter
    for (int k = 0;k<int(totalhists.size());++k){
        if (k == refhistidx){
            continue;
        }
        effh[ed]->Divide(totalhists.at(k),
                        totalhists.at(refhistidx), 1, 1, "B");
        ++ed;
    }

    //pass the reference histogram for drawing purposes
    refhidx = refhistidx;
    


}//----------makeEffRatioHist





//Get a vector of histograms which are the sum of the requested signal
//and backgrounds samples.  This can be accomodated to the user's needs
//---------------------------------------------------------------------------
vector<TH1F*> getTotalHists (vector<datahist>& hSamples,const string& signame,
                    const string& backname)
{
//---------------------------------------------------------------------------
    if(debugme) cout<<"Processing getTotalHIsts...."<<endl;


    //in case there are various backgrounds requested, form a vector
    //with them according to the needs.
    vector< string> typesOfBack;
    if (backname == "all"){
        typesOfBack.push_back("QCD");
        typesOfBack.push_back("Z");
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
    
    //get the signal (call it total as we will dump everything there) and
    //the backgrounds that will be added.
    vector<TH1F*> totalvec;
    vector<vector<TH1F*> > backmvec;
    for (unsigned int samp = 0;samp<hSamples.size();++samp){
        //get signal sample
        if (hSamples.at(samp).dir == signame) {
            cout<<"Signal is: "
                <<hSamples.at(samp).dir<<endl;
            totalvec.assign(hSamples.at(samp).h1f.begin(),hSamples.at(samp).h1f.end());
        }
        
        //if the background is just one sample,get it
        if (hSamples.at(samp).dir == backname){
            backmvec.push_back(hSamples.at(samp).h1f);
        }
        //if there are more samples in the background get them according
        //to the typesOfBack vector
        else{
            vector<string>::iterator it = find(typesOfBack.begin(),
                                               typesOfBack.end(),hSamples.at(samp).dir);
            if (it!=typesOfBack.end()) backmvec.push_back(hSamples.at(samp).h1f);
            else{cout<<"WARNING: sample:"<<hSamples.at(samp).dir
                     <<" was not considered in the total background vector"<<endl;}
        }
    }//sample loop


    //make sure the sizes didn't get screwed up, at least check the first one
    assert(totalvec.size() == backmvec.at(0).size());
    //add the backgrounds to the signal to form the total
    for (unsigned int cu = 0; cu<totalvec.size();++cu){
        totalvec.at(cu)->Sumw2();
        for(unsigned int mv = 0; mv<backmvec.size();++mv){
            backmvec.at(mv).at(cu)->Sumw2();
            totalvec.at(cu)->Add(backmvec.at(mv).at(cu));
        }
        
        //rebin
        //totalvec.at(cu)->Rebin(2);

        // set to zero all entries with (# of ev) < 2.5 
        //to avoid extra large error bars
        const float epsilon = 2.;
        for(int y = 1; y <= totalvec.at(cu)->GetNbinsX(); ++y){
            if(totalvec.at(cu)->GetBinContent(y) < epsilon){
                totalvec.at(cu)->SetBinContent(y, 0);
            }
        }//set to zero bins with small epsilon entries
    }//loop over totalvec


    

    return totalvec;


}//------------getHistTotal





//store directories and histograms within them in a container
//---------------------------------------------------------------------------
void get_histos (const string& filename,const string& wildcard, 
                 vector<datahist>& hSamples,const float& fixpar)
{
//---------------------------------------------------------------------------
    if (debugme) cout<<"Starting get_histos"<<endl;
    //transform the auxiliary float to a string to check.  If this is found
    //in the histogram name, i.e if its differen than -1,
    // the parameter will be fixed. This is to have a handle
    //on the number of curves in the plot, otherwise it becomes really busy
    char fixparstr[100];
    sprintf(fixparstr,"%.3f",fixpar);
    TFile* f = new TFile(filename.c_str(),"read");
    TIter nextkey(f->GetListOfKeys());
    TKey *key;
    while ((key=(TKey*)nextkey())) {
        TObject* obj = key->ReadObj();
        string dirname = obj->GetName();
        vector<TH1F*> histos_vector;
        if(debugme) cout<<"dirname = "<<dirname<<endl;
        TDirectory* tdir = (TDirectory*) f->GetDirectory(dirname.c_str());
        TIter hkeynext(tdir->GetListOfKeys());
        TKey *hkey;
        TH1F* h;
        while ((hkey = (TKey*)hkeynext())){
            TObject* objh = hkey->ReadObj();
            if (objh->IsA()->InheritsFrom("TH1F")){
                h = (TH1F*)objh;
                TString histName = h->GetName();
                if(debugme) cout<<"Hist name = "<<h->GetName()<<endl;
                if (!histName.Contains(wildcard.c_str())) continue;
                if (fixpar!= -1. && !histName.Contains(fixparstr)) continue;
                histos_vector.push_back(h);
            }//---if
        }//----while hist

        
        hSamples.push_back(datahist(dirname,histos_vector));

        
        
    }//---while dir
    
    return;

}//get_vector_of_histos




// Function to check the ratio of muon pT distributions when 
// compared to a default set of cuts.  The default parameters are set
// in par1 and par2.  For example, for isolation, par1 is sumPt and
// par2 is deltaR.  If, for some cut there is only one parameter,
// use the wildcard 99999 for the second par.  The precision of the 
// parameters is restricted to 0.001
//--------------------------------------------------------------------
void plotOptim(const string& var, const string& signame,
               const string& backname,const string& algo,
               const float& par1, const float& par2, 
               const float& fixpar = -1.)
{
//--------------------------------------------------------------------

    if(debugme) cout<<"Starting plot hist eff ratios function..."<<endl;
    
    //store all the histograms for different samples
    //but just for the algorithm and conditions required
    vector<datahist> hSamples;
    string filename = dirfiles+"h_"+var+".root";
    get_histos(filename,algo,hSamples,fixpar);
    //add the histograms according to users needs
    vector<TH1F*> totalhists = getTotalHists(hSamples,signame,backname);
    //effh to store the efficiency plots
    //dimension is all minus the reference plot
    const unsigned int effhdim = totalhists.size()-1;
    TH1F* effh[effhdim];
    //store the reference histogram for drawing
    int refhidx = -1;
    makeEffRatioHists(totalhists,refhidx,effh,algo,par1,par2);
    plotEffRatios(effh,effhdim,totalhists.at(refhidx),algo);

}




// Plot significance-like plots using different estimators
//--------------------------------------------------------------------
void plotOptim(const string& var, const string& signame,
               const string& backname,const string& algo,
               const string& myestimator)
{
//--------------------------------------------------------------------

    if(debugme) cout<<"Starting plot with estimators function..."<<endl;
    

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
   
   
    //prepare dimensions(names) of algos and pars
    vector<string> myAlgoName;
    int mydimPar1 = 1;//number of parameters 1
    int mydimPar2 = 1;//number of parameters 2
    get_dimensions(signif,myAlgoName,mydimPar1,mydimPar2);

    //prepare TGraphsErrors for different options
    //because we need different curves for different par2's:
    //loop over algos if requested
    const int myArrsize = mydimPar2;
    TGraphErrors* myg[myArrsize];
    plotCutOpt_TwoPars(var,signame, backname,myg,
                       mydimPar1,mydimPar2,signif,myestimator,algo,myAlgoName);
    



}//---------------------plotWithEstimators
  







