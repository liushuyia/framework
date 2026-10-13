void compare_nueCC(){
  TFile *file1 = new TFile("merge_r1.root");
  TFile *file2 = new TFile("./hist_rootfiles/DetVar/cov_LYDown.root");
  TFile *file3 = new TFile("./hist_rootfiles/XsFlux/cov_1.root");
  TH1F *h1 = (TH1F*)file1->Get("histo_1");
  h1->Scale(6./4.42723);
  TH1F *h2 = (TH1F*)file2->Get("pred_covch_1");
  h2->Scale(6./5.);
  TH1F *h3 = (TH1F*)file3->Get("pred_covch_1");
  h1->Draw("");
  h1->SetLineColor(1);
  h2->Draw("same");
  h2->SetLineColor(2);
  h3->SetLineColor(4);
  h3->Draw("same");

  TLegend *le1 = new TLegend(0.6,0.6,0.89,0.89);
  le1->AddEntry(h1,"CV","l");
  le1->AddEntry(h2,"Det (no BG)","l");
  le1->AddEntry(h3,"Xs+Flux","l");
  le1->Draw();
  
}
