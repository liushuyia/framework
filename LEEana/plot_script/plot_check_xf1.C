void plot_check_xf1(int ch = 1){
  TFile *file1 = new TFile("merge.root");
  TH1F **h1 = new TH1F*[9];
  TH1F **h2 = new TH1F*[9];
  TFile *file2 = new TFile(Form("./hist_rootfiles/XsFlux/cov_%d.root",ch));
  for (Int_t i=0;i!=9;i++){
    h1[i] = (TH1F*)file1->Get(Form("histo_%d",i+1));
    h2[i] = (TH1F*)file2->Get(Form("pred_covch_%d",i+1));
  }

  TCanvas *c1 = new TCanvas("c1","c1",1200,900);
  c1->Divide(4,3);
  for (Int_t i=0;i!=9;i++){
    c1->cd(i+1);
    h1[i]->Draw();
    h2[i]->Draw("same");
    h2[i]->SetLineColor(2);
  }
}
