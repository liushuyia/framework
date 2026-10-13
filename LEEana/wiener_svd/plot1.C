void plot1(){
  TFile *file = new TFile("wiener.root");
  TH2D *hR = (TH2D*)file->Get("hR");
  TH1D *hmeas = (TH1D*)file->Get("hmeas");
  TH2D *hcov_tot = (TH2D*)file->Get("hcov_tot");
  //hmeas->Draw();
  //hR->Draw("COLZ");
  hcov_tot->Draw("COLZ");
}
