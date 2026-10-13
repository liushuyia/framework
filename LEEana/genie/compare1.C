void compare1(){
  TFile *file1 = new TFile("nue_g2t1.root");
  TH1F *h_enu_ccqe1 = (TH1F*)file1->Get("h_enu_ccqe");
  TH1F *h_enu_ccmec1 = (TH1F*)file1->Get("h_enu_ccmec");
  TH1F *h_enu_ccres1 = (TH1F*)file1->Get("h_enu_ccres");
  TFile *file2 = new TFile("nue_g3.root");
  TH1F *h_enu_ccqe2 = (TH1F*)file2->Get("h_enu_ccqe");
  TH1F *h_enu_ccmec2 = (TH1F*)file2->Get("h_enu_ccmec");
  TH1F *h_enu_ccres2 = (TH1F*)file2->Get("h_enu_ccres");

  h_enu_ccqe1->Add(h_enu_ccmec1);
  h_enu_ccqe1->Add(h_enu_ccres1);

  h_enu_ccqe2->Add(h_enu_ccmec2);
  h_enu_ccqe2->Add(h_enu_ccres2);

  h_enu_ccqe1->Divide(h_enu_ccqe2);
  
  h_enu_ccqe1->Draw();
  // h_enu_ccmec1->Draw("same");
  // h_enu_ccres1->Draw("same");
  // h_enu_ccmec1->SetLineColor(2);
  // h_enu_ccres1->SetLineColor(4);
  
  h_enu_ccqe2->Draw("same");
  // h_enu_ccmec2->Draw("same");
  // h_enu_ccres2->Draw("same");
  // h_enu_ccmec2->SetLineColor(2);
  // h_enu_ccres2->SetLineColor(4);
  h_enu_ccqe2->SetMarkerStyle(20);
  // h_enu_ccmec2->SetMarkerStyle(20);
  // h_enu_ccres2->SetMarkerStyle(20);
  
}
