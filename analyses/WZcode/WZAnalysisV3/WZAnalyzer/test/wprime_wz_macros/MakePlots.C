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

void
MakePlots(){  
   
    TFile *fin = TFile::Open("Wprime_analysis.root");
    string title; 
    ///////////
    title = "hHt_muon";        DrawandSave(fin,title);
    title = "hWZInvMass_elec"; DrawandSave(fin,title);
    title = "hWZInvMass_muon"; DrawandSave(fin,title);
    title = "hWpt_muon";       DrawandSave(fin,title);
    title = "hZpt_muon";       DrawandSave(fin,title);
    
    /*
    h1 = (TH1F*) fin->Get("wprime400/hHt_muon");
    h2 = (TH1F*) fin->Get("wzjj/hHt_muon");

    sum = h1->Integral(); h1->Scale(1./sum);
    sum = h2->Integral(); h2->Scale(1./sum);
			      
    h1->SetLineColor(kRed);  h1->Draw();
    h2->SetLineColor(kBlue); h2->Draw("same");

    c1->SaveAs("hHt_muon.eps");

    ////////////
    
    h1 = (TH1F*) fin->Get("wprime400/hWZInvMass_muon");
    h2 = (TH1F*) fin->Get("wzjj/hWZInvMass_muon");

    sum = h1->Integral(); h1->Scale(1./sum);
    sum = h2->Integral(); h2->Scale(1./sum);
			      
    h1->SetLineColor(kRed);  h1->Draw();
    h2->SetLineColor(kBlue); h2->Draw("same");

    c1->SaveAs("hWZInvMass_muon.eps");

    ////////////

    h1 = (TH1F*) fin->Get("wprime400/hWZInvMass_elec");
    h2 = (TH1F*) fin->Get("wzjj/hWZInvMass_elec");

    sum = h1->Integral(); h1->Scale(1./sum);
    sum = h2->Integral(); h2->Scale(1./sum);
			      
    h1->SetLineColor(kRed);  h1->Draw();
    h2->SetLineColor(kBlue); h2->Draw("same");

    c1->SaveAs("hWZInvMass_elec.eps");

 ////////////

    h1 = (TH1F*) fin->Get("wprime400/hWpt_muon");
    h2 = (TH1F*) fin->Get("wzjj/hWpt_muon");

    sum = h1->Integral(); h1->Scale(1./sum);
    sum = h2->Integral(); h2->Scale(1./sum);
			      
    h1->SetLineColor(kRed);  h1->Draw();
    h2->SetLineColor(kBlue); h2->Draw("same");

    c1->SaveAs("hWpt_muon.eps");

 ////////////

    h1 = (TH1F*) fin->Get("wprime400/hZpt_muon");
    h2 = (TH1F*) fin->Get("wzjj/hZpt_muon");

    sum = h1->Integral(); h1->Scale(1./sum);
    sum = h2->Integral(); h2->Scale(1./sum);
			      
    h1->SetLineColor(kRed);  h1->Draw();
    h2->SetLineColor(kBlue); h2->Draw("same");

    c1->SaveAs("hZpt_muon.eps");
    */
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
