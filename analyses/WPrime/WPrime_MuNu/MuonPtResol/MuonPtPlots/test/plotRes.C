#define plotResHeader_cxx
#include "plotResHeader.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TH1.h>
#include <TF1.h>
#include <iostream>
#include <TProfile.h>
#include <TMath.h>
#include <TLegend.h>

TString plotdirname = "plots";

void plotRes(TString fileName,float ptmuon, TString gt)
{
    gROOT->ProcessLine(".L CMSStyle.C");
    gROOT->ProcessLine("CMSstyle()");
    
    

    plotResHeader t(fileName,ptmuon,gt);
  t.Loop();
}



void plotResHeader::Loop()
{

    //define histograms
    float maxptbin = thept+(thept/2.);
    float minptbin = thept-(thept/2.);
    int nptbins = int((maxptbin-minptbin)/5.);
  
    int nresbins = 50;
    float maxresbin = 0.5;
    float minresbin = -0.5;
    
    int nhitbins = 60;
    float minhitbin = 0;
    float maxhitbin = 60;
    
    int ndresbins = 20;
    float mindresbin = 0;
    float maxdresbin = 5;

    int nreshitbins = 50;
    float maxreshitbin = 0.5;
    float minreshitbin = 0;

    //individual pt distributions and resolution plots
    TH1F* d_pt_tev = new TH1F("d_pt_tev","Muon p_{T} distribution;(GeV);Events",nptbins,minptbin,maxptbin);
    TH1F* r_gen_tev = new TH1F("r_gen_tev","Relative muon (1/p_{T}) resolution;relative (1/p_{T}) resolution;Events",nresbins,minresbin,maxresbin);
    TH1F* d_pt_anew = new TH1F("d_pt_anew","Muon p_{T};(GeV);Events",nptbins,minptbin,maxptbin);
    TH1F* r_gen_anew = new TH1F("r_gen_anew","Relative muon (1/p_{T}) resolution;relative (1/p_{T}) resolution;Events",nresbins,minresbin,maxresbin);
    TH1F* d_pt_eve = new TH1F("d_pt_eve","Muon p_{T};(GeV);Events",nptbins,minptbin,maxptbin);
    TH1F* r_gen_eve = new TH1F("r_gen_eve","Relative muon (1/p_{T}) resolution;relative (1/p_{T}) resolution;Events",nresbins,minresbin,maxresbin);
    TH1F* d_pt_odd = new TH1F("d_pt_odd","Muon p_{T};(GeV);Events",nptbins,minptbin,maxptbin);
    TH1F* r_gen_odd = new TH1F("r_gen_odd","Relative muon (1/p_{T}) resolution;relative (1/p_{T}) resolution;Events",nresbins,minresbin,maxresbin);

    //resolution with the method
    TH1F* r_eve_odd = new TH1F("r_eve_odd","Relative muon (1/p_{T}) resolution;relative (1/p_{T}) resolution;Events",nresbins,minresbin,maxresbin);
 
   //resolution vs original number of hits
    TH2F* rh_tev = new TH2F("rh_tev","Relative muon (1/p_{T}) resolution vs number of hits;num hits;relative (1/p_{T}) resolution",nhitbins,minhitbin,maxhitbin,nreshitbins,minreshitbin,maxreshitbin);
    TH2F* rh_eve_odd = new TH2F("rh_eve_odd","Relative muon (1/p_{T}) resolution vs number of hits;num hits;relative (1/p_{T}) resolution",nhitbins,minhitbin,maxhitbin,nreshitbins,minreshitbin,maxreshitbin);
    //ratio of new/original resolutions vs original hits
    TH2F* Drh = new TH2F("Drh","Ratio of resolution and true resolution;num hits;(res)/(true res)",nhitbins,minhitbin,maxhitbin,ndresbins,mindresbin,maxdresbin);
    

    //original resolution for those that succeed in getting a new
    //resolution and for those that don't
    TH1F* r_gen_tev_succ = new TH1F("r_gen_tev_succ","Relative muon (1/p_{T}) resolution;relative (1/p_{T}) resolution;Events",nresbins,minresbin,maxresbin);
    TH1F* r_gen_tev_fail = new TH1F("r_gen_tev_fail","Relative muon (1/p_{T}) resolution;relative (1/p_{T}) resolution;Events",nresbins,minresbin,maxresbin);
    
    
    

   if (fChain == 0) return;

   Long64_t nentries = fChain->GetEntriesFast();

   Long64_t nbytes = 0, nb = 0;
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
      // if (Cut(ientry) < 0) continue;
    
      //discard events with more than 2 muons (this is weird I don't know
      //why it happens) and loop over muons
      if (nmu>2) continue;
      for (int imu = 0; imu<nmu; ++imu){

          //individual pt distributions and resolution plots
          if (tevvalid[imu]){
              d_pt_tev->Fill(tevpt[imu]);
              r_gen_tev->Fill(tevres[imu]);
          }
          if (anewvalid[imu]){
              d_pt_anew->Fill(anewpt[imu]);
              r_gen_anew->Fill(anewres[imu]);
          }
          if (oddvalid[imu]){
              d_pt_odd->Fill(oddpt[imu]);
              r_gen_odd->Fill(oddres[imu]);
          }
          if (evevalid[imu]){
              d_pt_eve->Fill(evept[imu]);
              r_gen_eve->Fill(everes[imu]);
          }

          //resolution with the new method
          if (evevalid[imu] && oddvalid[imu]){
              r_eve_odd->Fill(mres[imu]);
          }

          //plots resolution vs number of hits
          if(tevvalid[imu]){
              rh_tev->Fill(fabs(tevhits[imu]),tevres[imu]);
          }
          if (evevalid[imu] && oddvalid[imu]){
              rh_eve_odd->Fill(fabs(tevhits[imu]),mres[imu]);
          }

          //plot the ratio of two resolutions vs number of hits
          if(tevvalid[imu] && evevalid[imu] && oddvalid[imu]){
              float theratio = fabs(mres[imu])/fabs(tevres[imu]);
              Drh->Fill(tevhits[imu],theratio);
          }

          //compare resolution for those that succeed getting new seeds and
          //for those that don't
          if (tevvalid[imu]){
              if (evevalid[imu] && oddvalid[imu]){
                  r_gen_tev_succ->Fill(tevres[imu]);
              }
              if ((!evevalid[imu]) || (!oddvalid[imu])){
                  r_gen_tev_fail->Fill(tevres[imu]);
              }
          }


      }//loop over muons
   }//loop over events

   //plotting
   float mymean = 0;
   float mysig = 0;
   TString canvname = "";
   TString fileplot = "";

   //tev
   TCanvas* c_d_pt_tev = new TCanvas("c_d_pt_tev","c_d_pt_tev");
   c_d_pt_tev->cd();
   d_pt_tev->Draw();
   c_d_pt_tev->Update();
   canvname = c_d_pt_tev->GetName();
   fileplot = plotdirname+"/"+canvname+"_"+strpt+"_"+thegt+".eps";
   c_d_pt_tev->Print(fileplot);
   TCanvas* c_r_gen_tev = new TCanvas("c_r_gen_tev","c_r_gen_tev");
   c_r_gen_tev->cd();
   mymean = r_gen_tev->GetMean();
   mysig = r_gen_tev->GetRMS();
   r_gen_tev->Fit("gaus","R","",mymean-mysig,mymean+mysig);
   c_r_gen_tev->Update();
   canvname = c_r_gen_tev->GetName();
   fileplot = plotdirname+"/"+canvname+"_"+strpt+"_"+thegt+".eps";
   c_r_gen_tev->Print(fileplot);

   //anew
   TCanvas* c_d_pt_anew = new TCanvas("c_d_pt_anew","c_d_pt_anew");
   c_d_pt_anew->cd();
   d_pt_anew->Draw();
   c_d_pt_anew->Update();
   canvname = c_d_pt_anew->GetName();
   fileplot = plotdirname+"/"+canvname+"_"+strpt+"_"+thegt+".eps";
   c_d_pt_anew->Print(fileplot);
   TCanvas* c_r_gen_anew = new TCanvas("c_r_gen_anew","c_r_gen_anew");
   c_r_gen_anew->cd();
   mymean = r_gen_anew->GetMean();
   mysig = r_gen_anew->GetRMS();
   r_gen_anew->Fit("gaus","R","",mymean-mysig,mymean+mysig);
   c_r_gen_anew->Update();
   canvname = c_r_gen_anew->GetName();
   fileplot = plotdirname+"/"+canvname+"_"+strpt+"_"+thegt+".eps";
   c_r_gen_anew->Print(fileplot);

   //odd
   TCanvas* c_d_pt_odd = new TCanvas("c_d_pt_odd","c_d_pt_odd");
   c_d_pt_odd->cd();
   d_pt_odd->Draw();
   c_d_pt_odd->Update();
   canvname = c_d_pt_odd->GetName();
   fileplot = plotdirname+"/"+canvname+"_"+strpt+"_"+thegt+".eps";
   c_d_pt_odd->Print(fileplot);
   TCanvas* c_r_gen_odd = new TCanvas("c_r_gen_odd","c_r_gen_odd");
   c_r_gen_odd->cd();
   mymean = r_gen_odd->GetMean();
   mysig = r_gen_odd->GetRMS();
   r_gen_odd->Fit("gaus","R","",mymean-mysig,mymean+mysig);
   c_r_gen_odd->Update();
   canvname = c_r_gen_odd->GetName();
   fileplot = plotdirname+"/"+canvname+"_"+strpt+"_"+thegt+".eps";
   c_r_gen_odd->Print(fileplot);

   //even
   TCanvas* c_d_pt_eve = new TCanvas("c_d_pt_eve","c_d_pt_eve");
   c_d_pt_eve->cd();
   d_pt_eve->Draw();
   c_d_pt_eve->Update();
   canvname = c_d_pt_eve->GetName();
   fileplot = plotdirname+"/"+canvname+"_"+strpt+"_"+thegt+".eps";
   c_d_pt_eve->Print(fileplot);
   TCanvas* c_r_gen_eve = new TCanvas("c_r_gen_eve","c_r_gen_eve");
   c_r_gen_eve->cd();
   mymean = r_gen_eve->GetMean();
   mysig = r_gen_eve->GetRMS();
   r_gen_eve->Fit("gaus","R","",mymean-mysig,mymean+mysig);
   c_r_gen_eve->Update();
   canvname = c_r_gen_eve->GetName();
   fileplot = plotdirname+"/"+canvname+"_"+strpt+"_"+thegt+".eps";
   c_r_gen_eve->Print(fileplot);
   
   //resolution new method
   TCanvas* c_r_eve_odd = new TCanvas("c_r_eve_odd","c_r_eve_odd");
   c_r_eve_odd->cd();
   mymean = r_eve_odd->GetMean();
   mysig = r_eve_odd->GetRMS();
   r_eve_odd->Fit("gaus","R","",mymean-mysig,mymean+mysig);
   c_r_eve_odd->Update();
   canvname = c_r_eve_odd->GetName();
   fileplot = plotdirname+"/"+canvname+"_"+strpt+"_"+thegt+".eps";
   c_r_eve_odd->Print(fileplot);

   //plots resolution vs number of hits
   TCanvas* c_rh_tev = new TCanvas("c_rh_tev","c_rh_tev");
   c_rh_tev->cd();
   TProfile* P_rh_tev = (TProfile*) rh_tev->ProfileX("P_rh_tev",1,-1,"");
   P_rh_tev->SetLineColor(kRed);
   P_rh_tev->SetMarkerColor(kRed);
   TF1* F_rh_tev = new TF1("F_rh_tev","pol1",13,40);
   F_rh_tev->SetLineColor(kRed);
   P_rh_tev->Draw();
   //P_rh_tev->Fit(F_rh_tev,"","",16,40);
   canvname = c_rh_tev->GetName();
   fileplot = plotdirname+"/"+canvname+"_"+strpt+"_"+thegt+".eps";
   c_rh_tev->Print(fileplot);

   TCanvas* c_rh_eve_odd = new TCanvas("c_rh_eve_odd","c_rh_eve_odd");
   c_rh_eve_odd->cd();
   TProfile* P_rh_eve_odd = (TProfile*) rh_eve_odd->ProfileX("P_rh_eve_odd",1,-1,"");
   P_rh_eve_odd->SetLineColor(kBlack);
   P_rh_eve_odd->SetMarkerColor(kBlack);
   TF1* F_rh_eve_odd = new TF1("F_rh_eve_odd","pol1",13,40);
   //F_rh_eve_odd->SetParLimits(2,0,99999);
   F_rh_eve_odd->SetLineColor(kBlack);
   P_rh_eve_odd->Draw();
   //P_rh_eve_odd->Fit(F_rh_eve_odd,"","",16,40);
   canvname = c_rh_eve_odd->GetName();
   fileplot = plotdirname+"/"+canvname+"_"+strpt+"_"+thegt+".eps";
   c_rh_eve_odd->Print(fileplot);

   //plot the two previous plots in one single canvas
   TCanvas* c_all_rh = new TCanvas("c_all_rh","c_all_rh");
   c_all_rh->cd();
   P_rh_eve_odd->Draw();
   P_rh_tev->Draw("same");
   P_rh_eve_odd->SetStats(0);
   P_rh_tev->SetStats(0);
    TLegend* leg_all_rh = new TLegend(0.7,0.7,0.5,0.8);
    leg_all_rh->SetTextSize(0.04);
    leg_all_rh->AddEntry(P_rh_eve_odd,"Resolution","lp");
    leg_all_rh->AddEntry(P_rh_tev,"True Resolution","lp");
    leg_all_rh->SetLineColor(0);
    leg_all_rh->SetLineStyle(0);
    leg_all_rh->SetLineWidth(0);
    leg_all_rh->SetFillColor(0);
    leg_all_rh->Draw("same");
    canvname = c_all_rh->GetName();
    fileplot = plotdirname+"/"+canvname+"_"+strpt+"_"+thegt+".eps";
    c_all_rh->Print(fileplot);

    //the ratio of new/true (resolution) vs number of this
    TCanvas* c_Drh = new TCanvas("c_Drh","c_Drh");
    c_Drh->cd();
    TProfile* P_Drh = (TProfile*) Drh->ProfileX("P_Drh",1,-1,"");
    P_Drh->Draw();
    TF1* F_Drh = new TF1("F_Drh","pol0",13,40);
    P_Drh->Fit(F_Drh,"","",18,40);
    canvname = c_Drh->GetName();
    fileplot = plotdirname+"/"+canvname+"_"+strpt+"_"+thegt+".eps";
    c_Drh->Print(fileplot);


    //original resolution for those that succeed in getting a new
    //resolution and for those that don't
    TCanvas* c_r_gen_tev_succ = new TCanvas("c_r_gen_tev_succ","c_r_gen_tev_succ");
    c_r_gen_tev_succ->cd();
    mymean = r_gen_tev_succ->GetMean();
    mysig = r_gen_tev_succ->GetRMS();
    r_gen_tev_succ->Fit("gaus","R","",mymean-mysig,mymean+mysig);
    c_r_gen_tev_succ->Update();
    canvname = c_r_gen_tev_succ->GetName();
    fileplot = plotdirname+"/"+canvname+"_"+strpt+"_"+thegt+".eps";
    c_r_gen_tev_succ->Print(fileplot);

    TCanvas* c_r_gen_tev_fail = new TCanvas("c_r_gen_tev_fail","c_r_gen_tev_fail");
    c_r_gen_tev_fail->cd();
    mymean = r_gen_tev_fail->GetMean();
    mysig = r_gen_tev_fail->GetRMS();
    r_gen_tev_fail->Fit("gaus","R","",mymean-mysig,mymean+mysig);
    c_r_gen_tev_fail->Update();
    canvname = c_r_gen_tev_fail->GetName();
    fileplot = plotdirname+"/"+canvname+"_"+strpt+"_"+thegt+".eps";
    c_r_gen_tev_fail->Print(fileplot);

}


