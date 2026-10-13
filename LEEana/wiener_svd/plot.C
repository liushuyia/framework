void plot(){
  TFile *file = new TFile("merge_xs.root");
  TH1F *h1 = (TH1F*)file->Get("hdata_obsch_1");
  TH1F *h2 = (TH1F*)file->Get("hmc_obsch_1");
  TH1F *h3 = (TH1F*)file->Get("hdata_obsch_2");
  TH1F *h4 = (TH1F*)file->Get("hmc_obsch_2");

  TCanvas *c1 = new TCanvas("c1","c1",1200,600);
  c1->Divide(2,1);
  c1->cd(1);
  h1->Draw("E");
  h1->SetLineColor(1);
  h2->Draw("histsame");
  h2->SetLineColor(2);
  h1->SetTitle("#nu_{#mu}CC Fully Contained");
  h1->SetXTitle("E^{rec}_{#nu} (MeV)");
  TLegend *le1 = new TLegend(0.6,0.6,0.89,0.89);
  le1->SetHeader("5.3e19 POT");
  le1->AddEntry(h1,"Data (stat. only)","l");
  le1->AddEntry(h2,"Prediction","l");
  le1->Draw();

  h1->SetLineWidth(2);
  h2->SetLineWidth(2);
  
  h1->GetXaxis()->SetTitleSize(0.06);
  h1->GetXaxis()->SetLabelSize(0.06);
  h1->GetYaxis()->SetTitleSize(0.06);
  h1->GetYaxis()->SetLabelSize(0.06);
  
  c1->cd(2);
  h3->Draw("E");
  h3->SetLineColor(1);
  h4->Draw("histsame");
  h4->SetLineColor(2);
  h3->SetTitle("#nu_{#mu}CC Partially Contained");
  h3->SetXTitle("E^{rec}_{#nu} (MeV)");

  h3->GetXaxis()->SetTitleSize(0.06);
  h3->GetXaxis()->SetLabelSize(0.06);
  h3->GetYaxis()->SetTitleSize(0.06);
  h3->GetYaxis()->SetLabelSize(0.06);

  h3->SetLineWidth(2);
  h4->SetLineWidth(2);
}
