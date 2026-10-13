void plot_correlation(){
  auto file = TFile::Open("output.root");
  auto unfcov = (TH2D*)file->Get("unfcov");
  auto rho = (TH2D*)unfcov->Clone("rho"); // correlation coefficient
  rho->SetTitle("Correlation Coefficient");


  for (int i=0; i<rho->GetNbinsX(); i++) {
    for (int j=0; j<rho->GetNbinsY(); j++) {
      rho->SetBinContent(i+1, j+1, rho->GetBinContent(i+1,j+1)/ std::sqrt(unfcov->GetBinContent(i+1,i+1) * unfcov->GetBinContent(j+1,j+1)));
    }
  }

  auto c2 = new TCanvas("c2","c2",1200,600);
  c2->Divide(2,1);

  c2->cd(1);
  unfcov->GetXaxis()->SetTitle("E^{#nu}_{true} Bin");
  unfcov->GetYaxis()->SetTitle("E^{#nu}_{true} Bin");
  unfcov->GetXaxis()->SetTitleSize(0.05);
  unfcov->GetYaxis()->SetTitleSize(0.05);
  unfcov->Draw("colz");

  c2->cd(2);
  rho->GetXaxis()->SetTitle("E^{#nu}_{true} Bin");
  rho->GetYaxis()->SetTitle("E^{#nu}_{true} Bin");
  rho->GetXaxis()->SetTitleSize(0.05);
  rho->GetYaxis()->SetTitleSize(0.05);
  rho->Draw("colz");

}
