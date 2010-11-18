#include <TFile.h>
#include <TH1F.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TStyle.h>

#include <string>
#include <iostream>

#include "wprime_histo_constants.h"

using std::string; using std::cout; using std::endl;

// select tracking algorithm
// from 0 to Num_trkAlgos-1 (see wprime_histo_constants.h)
// NB: this is currently common for Muon-pt and transverse mass distributions
const unsigned tracking_option = 3;

// algorithm description
string algo; 

// histogram title;
string desc;


void doPlots(TFile * _file0, int option);

float Lumi_ipb = -1;

void plotMuPtandMT2()
{
  assert(tracking_option >=0 && tracking_option <= Num_trkAlgos-1);

  string input_file = "Wprime_analysis.root";
  TFile *_file0 = TFile::Open(input_file.c_str());
  if(!_file0 || _file0->IsZombie())
    {
      cout << " *** Ooops! Cannot find input file " << input_file << endl;
      return;
    }

  string histo_lumi = "lumi_ipb";
  TH1F * lumi_ipb = (TH1F* )_file0->Get(histo_lumi.c_str());
  Lumi_ipb = lumi_ipb->GetBinContent(1);

  gStyle->SetOptStat(00000);

  algo = algo_desc_short[tracking_option];
  char temp[1024]; sprintf(temp, " muons  (STARTUP, %4.2f pb^{-1})", Lumi_ipb);
  desc = algo_desc_long[tracking_option] + temp;
  doPlots(_file0, 1); // muon-pt
  doPlots(_file0, 2); // transverse mass

}

bool badHisto(TH1F * h, string s)
{
  if(!h)
    {
      cout << " *** Oops! Didn't find " << s << " histogram... " << endl;
      return true;
    }

  return false;
}

// option = 1 (muon-pt), option = 2 (transverse mass)
void doPlots(TFile * _file0, int option)
{

  //  gStyle->Reset();
  gStyle->SetFillColor(1);


  TCanvas * c1 = new TCanvas();
  c1->SetLogy();

  string prefix = "";
  if(option == 1)
    prefix = "hPT";
  else if(option == 2)
    prefix = "hTM";

  string histo = prefix + algo + "_" + final_histo_desc;
  string histoW = "W/" + histo;
  TH1F * w = (TH1F* ) _file0->Get(histoW.c_str());
  if(badHisto(w, "W"))
    return;

  string histoWlpt = "Wlowpt/" + histo;
  TH1F * wlpt = (TH1F* ) _file0->Get(histoWlpt.c_str());
  if(badHisto(wlpt, "Wlowpt"))
      return;

  string histoQ = "QCD/" + histo;
  TH1F * qcd = (TH1F* )_file0->Get(histoQ.c_str());
  if(badHisto(qcd, "QCD"))
    return;

  string histoZ = "Z/" + histo;
  TH1F * z = (TH1F* )_file0->Get(histoZ.c_str());
  if(badHisto(z, "Z"))
    return;
  
  string histoT = "Top/" + histo;
  TH1F * top = (TH1F* )_file0->Get(histoT.c_str());
  if(badHisto(top, "Top"))
    return;

  string histo10 = "wprime1.0/" + histo;
  TH1F * wp10 = (TH1F* )_file0->Get(histo10.c_str());
  if(badHisto(wp10, "wprime1.0"))
    return;

  string histo15 = "wprime1.5/" + histo;
  TH1F * wp15 = (TH1F* )_file0->Get(histo15.c_str());
  if(badHisto(wp15, "wprime1.5"))
    return;

  string histo20 = "wprime2.0/" + histo;
  TH1F * wp20 = (TH1F* )_file0->Get(histo20.c_str());
  if(badHisto(wp20, "wprime2.0"))
    return;

  string hdata = "data/" + histo;
  TH1F * data = (TH1F*) _file0->Get(hdata.c_str());
  if(badHisto(data, "data"))
    return;

  string hname = "tot_bgd";
  if(option == 1)
    hname += "_mupt";
  else if(option == 2)
    hname += "_TM";
  // this is the total background distribution (W + QCD + top + Z/DY)
  TH1F * bgd = new TH1F(hname.c_str(), desc.c_str(), 
			top->GetNbinsX(), top->GetXaxis()->GetXmin(), 
			top->GetXaxis()->GetXmax());
  bgd->Add(z);
  bgd->Add(qcd);  
  bgd->Add(top);
  bgd->Add(w);
  bgd->Add(wlpt);

  string desc = "";
  if(option == 1)
    desc = " p_{T} distribution";
  else if(option == 2)
    desc = " muon + (ckt-corrected) pfMET M_{T} distribution";
  string new_title = algo_desc_long[tracking_option] + desc;
  data->SetTitle(new_title.c_str());
  //data->SetMarkerStyle(4);
  //data->SetMarkerSize(1.3);
  data->SetMarkerStyle(8);
  //data->SetMarkerSize(1.);
  if(option == 1)
    {
      data->GetXaxis()->SetTitle("Muon p_{T} (GeV/c)");
      data->GetXaxis()->SetRangeUser(0, 500);
    }
  else if(option == 2)
    {
      data->GetXaxis()->SetTitle("M_{T} (GeV/c^{2})");
      data->GetXaxis()->SetRangeUser(0, 800);
    }

  if(data->GetMinimum() < 0.00001)data->SetMinimum(0.00001);
  data->Draw("e");
  bgd->Draw("same");
 


  const int fill_style_sig = 3001;
  const int fill_style_bgd = 3001;

  bgd->SetLineColor(kAzure+1);
  wp10->SetLineColor(kRed);
  wp15->SetLineColor(kRed+1);
  wp20->SetLineColor(kRed+2);
  bgd->SetFillColor(kAzure+1);
  wp10->SetFillColor(kRed);
  wp15->SetFillColor(kRed+1);
  wp20->SetFillColor(kRed+2);
  bgd->SetFillStyle(fill_style_bgd);
  wp10->SetFillStyle(fill_style_sig);
  wp15->SetFillStyle(fill_style_sig);
  wp20->SetFillStyle(fill_style_sig);
 
   
  wp10->Draw("same");
  wp15->Draw("same");
  wp20->Draw("same");

  char temp2[1024]; sprintf(temp2, " data (%4.2f pb^{-1})", Lumi_ipb);

  TLegend * lg = new TLegend(0.59, 0.67, 0.89, 0.89);
  lg->SetTextSize(0.03);
  lg->SetBorderSize(0);
  lg->SetFillColor(0);
  lg->AddEntry(bgd, "Total bgd", "F");
  lg->AddEntry(wp10, "W ' (1.0 TeV)", "F");
  lg->AddEntry(wp15, "W ' (1.5 TeV)", "F");
  lg->AddEntry(wp20, "W ' (2.0 TeV)", "F");
  lg->AddEntry(data, temp2, "LP");
  lg->Draw();
  
  string desc2 = "";
  if(option == 1)
    desc2 = "MuPt";
  else
    desc2 = "TM";
  char temp3[1024]; sprintf(temp3, "_%4.2fipb", Lumi_ipb);
  string file = desc2 + temp3 + ".gif";

  c1->SaveAs(file.c_str());
  //  c1->SaveAs(file2.c_str());
  //  delete c1;

}
