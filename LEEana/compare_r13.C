void compare_r13(){
  /*
  TFile *file1 = new TFile("./processed_checkout_rootfiles/run1_data_bnb_merge.root"); double pot_1 = 4.47375e+19; // 4.42457e19 (with good cut)...
  TFile *file3 = new TFile("./processed_checkout_rootfiles/checkout_data_bnb_run3_1e19.root"); double pot_3 = 9e+18;// 7.05836e18 (without good cut);

   // TFile *file1 = new TFile("./processed_checkout_rootfiles/checkout_data_extbnb_run1.root"); double pot_1 = 7.37133e+19;  // 7.73973e+19
   // TFile *file3 = new TFile("./processed_checkout_rootfiles/checkout_data_extbnb_run3.root"); double pot_3 = 2.95674e+20; // 3.15411e20 ...

  TH1F *h10 = new TH1F("h10","h10",25,0,2500);
  TH1F *h20 = new TH1F("h20","h20",25,0,2500);
  TGraphErrors *g1 = new TGraphErrors();
  TGraphErrors *g2 = new TGraphErrors();
  {
    TTree *T_eval = (TTree*)file1->Get("wcpselection/T_eval");
    TTree *T_BDTvars = (TTree*)file1->Get("wcpselection/T_BDTvars");
    TTree *T_KINEvars = (TTree*)file1->Get("wcpselection/T_KINEvars");
    T_eval->AddFriend(T_BDTvars,"T_BDTvars");
    T_eval->AddFriend(T_KINEvars,"T_KINEvars");

    T_eval->Project("h10","T_KINEvars.kine_reco_Enu","T_BDTvars.numu_score>0.9");
    //T_eval->Project("h10","T_KINEvars.kine_reco_Enu","T_BDTvars.numu_cc_flag>=0");
    for (Int_t i=0;i!=h10->GetNbinsX();i++){
      g1->SetPoint(i, h10->GetBinCenter(i+1), h10->GetBinContent(i+1) * 5e19/pot_1);
      g1->SetPointError(i, 0, sqrt(h10->GetBinContent(i+1)) * 5e19/pot_1);
    }
  }


  {
    TTree *T_eval = (TTree*)file3->Get("wcpselection/T_eval");
    TTree *T_BDTvars = (TTree*)file3->Get("wcpselection/T_BDTvars");
    TTree *T_KINEvars = (TTree*)file3->Get("wcpselection/T_KINEvars");
    T_eval->AddFriend(T_BDTvars,"T_BDTvars");
    T_eval->AddFriend(T_KINEvars,"T_KINEvars");
    T_eval->Project("h20","T_KINEvars.kine_reco_Enu","T_BDTvars.numu_score>0.9");
    //T_eval->Project("h20","T_KINEvars.kine_reco_Enu","T_BDTvars.numu_cc_flag>=0");
    for (Int_t i=0;i!=h20->GetNbinsX();i++){
      g2->SetPoint(i, h20->GetBinCenter(i+1), h20->GetBinContent(i+1) * 5e19/pot_3);
      g2->SetPointError(i, 0, sqrt(h20->GetBinContent(i+1)) * 5e19/pot_3);
    }
  }

  g1->Draw("A*");
  g2->Draw("*same");
  g1->SetMarkerColor(1);
  g1->SetMarkerStyle(20);
  g1->SetLineColor(1);
  g2->SetMarkerColor(2);
  g2->SetMarkerStyle(21);
  g2->SetLineColor(2);

  */

  TFile *file1 = new TFile("merge_r1.root");
  TFile *file3 = new TFile("merge_r3.root");
 
  TCanvas *c1 = new TCanvas("c1","c1",1200,600);
  c1->Divide(2,1);
  c1->cd(1);

  TH1F *he1 = (TH1F*)file1->Get("histo_1");
  TH1F *he3 = (TH1F*)file3->Get("histo_1");
  he1->Scale(5e19/4.427e19);
  he3->Scale(5e19/9.0e18);

  he1->Draw();
  he3->Draw("same");
  he3->SetLineWidth(2);
  he1->SetLineWidth(2);
  he3->SetLineColor(2);

  he1->GetYaxis()->SetLabelSize(0.07);
  he1->GetYaxis()->SetTitleSize(0.07);
  he1->GetXaxis()->SetLabelSize(0.07);
  he1->GetXaxis()->SetTitleSize(0.07);
  he1->GetXaxis()->SetNdivisions(506);
  he1->GetYaxis()->SetNdivisions(506);

  he1->SetTitle("#nu_{e}CC FC @ 5e19 POT");
  he1->SetXTitle("E_{#nu}^{reco} (GeV)");

  
  c1->cd(2);
  TH1F *hm1 = (TH1F*)file1->Get("histo_3");
  TH1F *hm3 = (TH1F*)file3->Get("histo_3");
  hm1->Scale(5e19/4.427e19);
  hm3->Scale(5e19/9.0e18);

  hm1->Draw();
  hm3->Draw("same");
  hm3->SetLineWidth(2);
  hm1->SetLineWidth(2);
  hm3->SetLineColor(2);

  hm1->GetYaxis()->SetLabelSize(0.07);
  hm1->GetYaxis()->SetTitleSize(0.07);
  hm1->GetXaxis()->SetLabelSize(0.07);
  hm1->GetXaxis()->SetTitleSize(0.07);
  hm1->GetXaxis()->SetNdivisions(506);
  hm1->GetYaxis()->SetNdivisions(506);

  hm1->SetTitle("#nu_{#mu}CC FC @ 5e19 POT");
  hm1->SetXTitle("E_{#nu}^{reco} (GeV)");

  TLegend *le1 = new TLegend(0.6,0.6,0.89,0.89);
  le1->AddEntry(hm1,"Run 1","l");
  le1->AddEntry(hm3,"Run 3","l");
  le1->Draw();
  
  
}
