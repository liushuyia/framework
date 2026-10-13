void plot_check_xs(){
  TFile *file = new TFile("./hist_rootfiles/XsFlux/cov_xs.root");
  TH1F *pred_covch_1 = (TH1F*)file->Get("pred_covch_1");
  TH1F *pred_covch_1_signal = (TH1F*)file->Get("pred_covch_1_signal");
  TH2F *pred_covch_1_R = (TH2F*)file->Get("pred_covch_1_R");
  TVectorD *vec_mean_17 = (TVectorD*)file->Get("vec_mean_17");
  TVectorD *vec_signal_17 = (TVectorD*)file->Get("vec_signal_17");
  TMatrixD *mat_R_17 = (TMatrixD*)file->Get("mat_R_17");
  
  TCanvas *c1 = new TCanvas("c1","c1",1200,800);
  c1->Divide(2,2);
  c1->cd(1);
  
  pred_covch_1->Draw();
  TGraph *g1 = new TGraph();
  for (Int_t i=0;i!=pred_covch_1->GetNbinsX()+1;i++){
    g1->SetPoint(i,i+1, pred_covch_1->GetBinContent(i+1));
  }
  g1->Draw("A*");
  g1->SetMarkerStyle(20);

  TGraph *g2 = new TGraph();
  for (Int_t i=0;i!=pred_covch_1_R->GetNbinsX()+1;i++){
    double sum = 0;
    for (Int_t j=0;j!=pred_covch_1_R->GetNbinsY();j++){
      sum += pred_covch_1_R->GetBinContent(i+1,j+1);
    }
    g2->SetPoint(i,i+1,sum);
  }
  g2->Draw("*same");
  g2->SetMarkerStyle(20);
  g2->SetMarkerColor(2);

  // TGraph *g3 = new TGraph();
  // for (Int_t i=0;i!=pred_covch_1_R->GetNbinsX()+1;i++){
  //   g3->SetPoint(i,i+1, (*vec_mean_17)[i]);
  // }
  // g3->Draw("*same");
  // g3->SetMarkerStyle(20);
  // g3->SetMarkerColor(4);
  
  c1->cd(2);
  pred_covch_1_signal->Draw("");
  vec_signal_17->Draw("same");
  //std::cout << pred_covch_1_signal->GetBinContent(0) << " " << pred_covch_1_signal->GetBinContent(15) << std::endl;
  c1->cd(3);
  mat_R_17->Draw("COLZ");

  TVectorD pred = (*mat_R_17) * (*vec_signal_17);
  TGraph *g4 = new TGraph();
  for (Int_t i=0;i!=26;i++){
    g4->SetPoint(i,i+1,pred(i));
  }
  c1->cd(1);
  g4->Draw("*same");
  g4->SetMarkerColor(4);
  g4->SetMarkerStyle(20);

  // TFile *file = new TFile("merge_xs.root");
  // TH1F *histo_1 = (TH1F*)file->Get("histo_1");
  // TH1F *histo_2 = (TH1F*)file->Get("histo_2");

  // TMatrixD* mat_R = (TMatrixD*)file->Get("mat_R");
  // TVectorD* vec_signal = (TVectorD*)file->Get("vec_signal");
  // TVectorD result = (*mat_R)*(*vec_signal);
  // TGraph *g1 = new TGraph();
  // for (Int_t i=0;i!=52;i++){
  //   g1->SetPoint(i,i+1,result(i));
  // }
  // g1->Draw("A*");
  // g1->SetMarkerColor(1);
  // g1->SetMarkerStyle(20);
  // TGraph *g2 = new TGraph();
  // for (Int_t i=0;i!=26;i++){
  //   g2->SetPoint(i,i+1,histo_1->GetBinContent(i+1));
  // }
  // for (Int_t i=0;i!=26;i++){
  //   g2->SetPoint(26+i,i+27,histo_2->GetBinContent(i+1));
  // }
  // g2->Draw("*same");
  // g2->SetMarkerColor(2);
  // g2->SetMarkerStyle(20);
}
