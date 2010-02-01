#include <stdio.h>
#include <math.h>
#include <string>
#include <iostream>
#include "TFile.h"
#include "TTree.h"
#include "TH1F.h"
#include "TCanvas.h"
#include "TRandom.h"
#include "TMath.h"

void DrawandSave(TFile* fin, string title);
void PlotEff(TFile* fin, string title);

void
MakePlots(){  
    const int Nplots = 5;

    TFile *fin = TFile::Open("Wprime_analysis.root");
    string title[Nplots]; 
    string efftitle[3];

    title[0] = "hHt_muon";       
    title[1] = "hWZInvMass_elec";
    title[2] = "hWZInvMass_muon";
    title[3] = "hWpt_muon";      
    title[4] = "hZpt_muon";      

    efftitle[0] = "hEffAbs";
    efftitle[1] = "hEffRel";
    efftitle[2] = "hNumEvts";

    for(int i=0;i<Nplots;++i) DrawandSave(fin,title[i]);
    
    for(int i=0;i<3;++i) PlotEff(fin,efftitle[i]);
}


void
DrawandSave(TFile* fin, string title){
    string dir1 = "wprime400";
    string dir2 = "wzjj";

    string h1_name = dir1 + "/" + title;
    string h2_name = dir2 + "/" + title;
    string filename = title + ".eps";

    float sum,max1,max2,max;
    int first, last;

    TCanvas c1;
    TH1F *h1,*h2;
    TAxis *axis;

    h1 = (TH1F*) fin->Get(h1_name.c_str());
    axis = h1->GetXaxis();
    first = axis->GetFirst();
    last  = axis->GetLast();
    sum = h1->Integral(first+1,last); 
    h1->Scale(1./sum);
    axis->SetRange(first+1,last);

    h2 = (TH1F*) fin->Get(h2_name.c_str());
    axis = h2->GetXaxis();
    first = axis->GetFirst();
    last  = axis->GetLast();
    sum = h2->Integral(first+1,last); 
    h2->Scale(1./sum);
    axis->SetRange(first+1,last);

    max1 = h1->GetMaximum();
    max2 = h2->GetMaximum();

    max = TMath::Max(max1,max2);
    h1->SetMaximum(max);

    h1->SetLineColor(kRed);  h1->Draw();
    h2->SetLineColor(kBlue); h2->Draw("same");

    c1.SaveAs(filename.c_str());
}

void
MakePlots2(){
    
}

void
PlotEff(TFile* fin, string title){
    string dir1 = "wprime400";
    string dir2 = "wzjj";

    string h1_name = dir1 + "/" + title;
    string h2_name = dir2 + "/" + title;
    string filename = title + ".eps";

    TCanvas c1;
    TH1F *h1,*h2;
        
    h1 = (TH1F*) fin->Get(h1_name.c_str());
    h2 = (TH1F*) fin->Get(h2_name.c_str());

    h1->SetStats(kFALSE);
    h2->SetStats(kFALSE);

    h1->SetLineColor(kRed);  h1->Draw();
    h2->SetLineColor(kBlue); h2->Draw("same");

    c1.SaveAs(filename.c_str());
}
