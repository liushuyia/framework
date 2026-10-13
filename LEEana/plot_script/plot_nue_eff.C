void plot_nue_eff(){
  TFile *file = new TFile("./checkout_rootfiles_correct_bdt/nue_test.root");
  TTree *T_eval = (TTree*)file->Get("wcpselection/T_eval");
  TTree *T_BDTvars = (TTree*)file->Get("wcpselection/T_BDTvars");
  TTree *T_PFeval = (TTree*)file->Get("wcpselection/T_PFeval");
  TTree *T_KINEvars = (TTree*)file->Get("wcpselection/T_KINEvars");

  T_PFeval->AddFriend(T_eval);
  T_PFeval->AddFriend(T_BDTvars);

  TH2F *h1 = new TH2F("h1","h1",14,0,1.400,20,-1,1);
  TH2F *h2 = new TH2F("h2","h2",14,0,1.400,20,-1,1);
  T_PFeval->Project("h1","truth_showerMomentum[2]/sqrt(pow(truth_showerMomentum[0],2)+pow(truth_showerMomentum[1],2)+pow(truth_showerMomentum[2],2)):truth_showerMomentum[3]","(T_eval.truth_vtxInside==1)");
  T_PFeval->Project("h2","truth_showerMomentum[2]/sqrt(pow(truth_showerMomentum[0],2)+pow(truth_showerMomentum[1],2)+pow(truth_showerMomentum[2],2)):truth_showerMomentum[3]","(T_eval.truth_vtxInside==1 && T_BDTvars.nue_score>7.0 && T_BDTvars.numu_cc_flag>=0)");

  TCanvas *c1 = new TCanvas("c1","c1",1200,900);
  c1->Divide(2,2);
  c1->cd(1);
  h1->Draw("COLZ");
  h1->SetTitle("Before Cuts");
  h1->SetXTitle("E_{shower} (GeV)");
  h1->SetYTitle("cos(#theta_{shower})");
  c1->cd(2);
  TH2F *h3 = (TH2F*)h2->Clone("h3");
  h3->Draw("COLZ");
  h3->SetTitle("After Cuts");
  h3->SetXTitle("E_{shower} (GeV)");
  h3->SetYTitle("cos(#theta_{shower})");
  c1->cd(3);
  h2->Divide(h1);
  h2->Draw("COLZ");
  h2->SetTitle("Efficiency");
  h2->SetXTitle("E_{shower} (GeV)");
  h2->SetYTitle("cos(#theta_{shower})");
  c1->cd(4);
  TGraph *g1 = new TGraph();
  TGraph *g2 = new TGraph();
  TGraph *g3 = new TGraph();
  TGraph *g4 = new TGraph();
  for (Int_t i=0;i!=20;i++){
    g1->SetPoint(i,h2->GetYaxis()->GetBinCenter(i+1),h2->GetBinContent(2,i+1));
    g2->SetPoint(i,h2->GetYaxis()->GetBinCenter(i+1),h2->GetBinContent(3,i+1));
    g3->SetPoint(i,h2->GetYaxis()->GetBinCenter(i+1),h2->GetBinContent(4,i+1));
    g4->SetPoint(i,h2->GetYaxis()->GetBinCenter(i+1),h2->GetBinContent(5,i+1));
  }

  
  g4->Draw("AL*");
  g3->Draw("L*same");
  g2->Draw("L*same");
  g1->Draw("L*same");

  g4->GetYaxis()->SetRangeUser(0.0,0.5);
  g4->SetLineColor(1); g4->SetMarkerStyle(20); g4->SetMarkerColor(1);
  g3->SetLineColor(2); g3->SetMarkerStyle(20); g3->SetMarkerColor(2);
  g2->SetLineColor(4); g2->SetMarkerStyle(20); g2->SetMarkerColor(4);
  g1->SetLineColor(6); g1->SetMarkerStyle(20); g1->SetMarkerColor(6);

  g4->GetXaxis()->SetTitle("cos(#theta_{shower})");
  g4->GetYaxis()->SetTitle("Efficiency");
  TLegend *le1 = new TLegend(0.6,0.6,0.89,0.89);
  le1->AddEntry(g4,"400-500 MeV","lp");
  le1->AddEntry(g3,"300-400 MeV","lp");
  le1->AddEntry(g2,"200-300 MeV","lp");
  le1->AddEntry(g1,"100-200 MeV","lp");
  le1->Draw();
  
  
}
