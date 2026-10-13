void compare_merge(){
  TFile *file1 = new TFile("merge_r1.root"); double pot_1 = 4.42457e19 ;
  TFile *file2 = new TFile("merge_r3.root"); double pot_3 = 6.27737e+18;
  //TFile *file3 = new TFile("merge.root");

  TH1F *h10 = (TH1F*)file1->Get("hmc_obsch_3");
  TH1F *h20 = (TH1F*)file2->Get("hmc_obsch_3");
  TH1F *h30 = (TH1F*)file1->Get("hmc_obsch_4");
  TH1F *h40 = (TH1F*)file2->Get("hmc_obsch_4");
  //  TH1F *h30 = (TH1F*)file3->Get("hmc_obsch_3");

  TH1F *h11 = (TH1F*)file1->Get("hdata_obsch_3");
  TH1F *h21 = (TH1F*)file2->Get("hdata_obsch_3");
  TH1F *h31 = (TH1F*)file1->Get("hdata_obsch_4");
  TH1F *h41 = (TH1F*)file2->Get("hdata_obsch_4");
  //  TH1F *h31 = (TH1F*)file3->Get("hdata_obsch_3");

  h10->Add(h30);
  h20->Add(h40);

  h11->Add(h31);
  h21->Add(h41);
  
  h10->Scale(5e19/pot_1);
  h20->Scale(5e19/pot_3);

  h11->Scale(5e19/pot_1);
  h21->Scale(5e19/pot_3);
  // h10->Add(h20);

  // h30->Draw();
  h10->SetLineColor(1);
  h20->SetLineColor(2);
  h10->Draw("");
  h10->SetTitle("numuCC");
  h10->GetYaxis()->SetRangeUser(0,1600);
  h20->Draw("same");

  TGraphErrors *g1 = new TGraphErrors();
  TGraphErrors *g2 = new TGraphErrors();
  for (Int_t i=0;i!=h11->GetNbinsX();i++){
    g1->SetPoint(i, h11->GetBinCenter(i+1), h11->GetBinContent(i+1));
    g2->SetPoint(i, h21->GetBinCenter(i+1), h21->GetBinContent(i+1));
    
    g1->SetPointError(i,0, sqrt(h11->GetBinContent(i+1)/5e19*pot_1) * 5e19/pot_1);
    g2->SetPointError(i,0, sqrt(h21->GetBinContent(i+1)/5e19*pot_3) * 5e19/pot_3);
    
    //h11->Draw("same");
    //h11->SetLineColor(4);
    //h21->Draw("same");
    //h21->SetLineColor(6);
  }

  g1->Draw("*same");
  g2->Draw("*same");
  g1->SetMarkerStyle(20);
  g2->SetMarkerStyle(20);
  g1->SetLineColor(1);
  g2->SetLineColor(2);
  g1->SetMarkerColor(1);
  g2->SetMarkerColor(2);
}
