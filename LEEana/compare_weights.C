void compare_weights(){
  gStyle->SetOptStat(0);
  

  
  TFile *file1 = new TFile("merge_all_new_weights.root");
  TFile *file2 = new TFile("merge_all.root");

  TH1F *he1 = (TH1F*)file1->Get("histo_1");
  TH1F *hm1 = (TH1F*)file1->Get("histo_3");

  TH1F *he2 = (TH1F*)file2->Get("histo_1");
  TH1F *hm2 = (TH1F*)file2->Get("histo_3");

  TCanvas *c1 = new TCanvas("c1","c1",1200,600);
  c1->Divide(2,1);
  c1->cd(1);
  he1->Draw();
  he2->SetLineColor(1);
  he2->Draw("same");
  he1->SetXTitle("E_{#nu}^{reco} (MeV)");
  he1->SetYTitle("# Events");
  he1->SetTitle("#nu_{e}CC FC");

  he1->GetYaxis()->SetLabelSize(0.07);
  he1->GetYaxis()->SetTitleSize(0.07);
  he1->GetXaxis()->SetLabelSize(0.07);
  he1->GetXaxis()->SetTitleSize(0.07);
  he1->GetXaxis()->SetNdivisions(506);
  he1->GetYaxis()->SetNdivisions(506);

  
  c1->cd(2);
  hm1->Draw();
  hm2->SetLineColor(1);
  hm2->Draw("same");
  hm1->SetXTitle("E_{#nu}^{reco} (MeV)");
  hm1->SetYTitle("# Events");
  hm1->SetTitle("#nu_{#mu}CC FC");
  he1->SetLineWidth(2);
  he2->SetLineWidth(2);
  hm1->SetLineWidth(2);
  hm2->SetLineWidth(2);
  he1->SetLineColor(2);
  hm1->SetLineColor(2);

  hm1->GetYaxis()->SetLabelSize(0.07);
  hm1->GetYaxis()->SetTitleSize(0.07);
  hm1->GetXaxis()->SetLabelSize(0.07);
  hm1->GetXaxis()->SetTitleSize(0.07);
  hm1->GetXaxis()->SetNdivisions(506);
  hm1->GetYaxis()->SetNdivisions(506);

  
  TLegend *le1 = new TLegend(0.6,0.6,0.89,0.89);
  le1->SetHeader("Open Data");
  le1->AddEntry(hm1,"new weights","l");
  le1->AddEntry(hm2,"old weights","l");
  le1->Draw();
  
  //std::cout << he1->Integral(0,6) << " " << he2->Integral(0,6) << std::endl;
}
