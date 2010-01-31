{
    TFile *fin = TFile::Open("Wprime_analysis.root");

    TCanvas c1;
    TH1F *h1, *h2;

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

}
