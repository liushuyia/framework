void plot_eff_nueCC_full(){
  TFile *file1 =  new TFile("./processed_checkout_rootfiles/checkout_prodgenie_bnb_intrinsic_nue_overlay_run123_all.root");
 
  TH1F *h10 = new TH1F("h10","h10",25,0,2500);
  TH1F *h11 = new TH1F("h11","h11",25,0,2500);
  TH1F *h12 = new TH1F("h12","h12",25,0,2500);

  TH1F *h30 = new TH1F("h30","h30",25,0,2500);
  TH1F *h31 = new TH1F("h31","h31",25,0,2500);
  TH1F *h32 = new TH1F("h32","h32",25,0,2500);

  TH1F *h40 = new TH1F("h40","h40",25,0,2500);
  TH1F *h41 = new TH1F("h41","h41",25,0,2500);
  TH1F *h42 = new TH1F("h42","h42",25,0,2500);

  TTree *T_eval = (TTree*)file1->Get("wcpselection/T_eval");
  TTree *T_BDTvars = (TTree*)file1->Get("wcpselection/T_BDTvars");
  TTree *T_KINEvars = (TTree*)file1->Get("wcpselection/T_KINEvars");
  TTree *T_PFeval = (TTree*)file1->Get("wcpselection/T_PFeval");
  T_eval->AddFriend(T_BDTvars,"T_BDTvars");
  T_eval->AddFriend(T_KINEvars,"T_KINEvars");
  T_eval->AddFriend(T_PFeval,"T_PFeval");
  
 
  
  T_eval->Project("h10","truth_nuEnergy","weight_cv*weight_spline*(abs(truth_nuPdg)==12 && truth_isCC==1 && truth_vtxInside==1)");
  T_eval->Project("h11","truth_nuEnergy","weight_cv*weight_spline*(abs(truth_nuPdg)==12 && truth_isCC==1 && truth_vtxInside==1 && T_BDTvars.numu_cc_flag>=0 && T_BDTvars.nue_score>7.0 && match_isFC==1)");
  T_eval->Project("h12","truth_nuEnergy","weight_cv*weight_spline*(abs(truth_nuPdg)==12 && truth_isCC==1 && truth_vtxInside==1 && T_BDTvars.numu_cc_flag>=0 && T_BDTvars.nue_score>7.0 )");


  T_eval->Project("h30","T_PFeval.truth_showerMomentum[3]*1000.","weight_cv*weight_spline*(abs(truth_nuPdg)==12 && truth_isCC==1 && truth_vtxInside==1)");
  T_eval->Project("h31","T_PFeval.truth_showerMomentum[3]*1000.","weight_cv*weight_spline*(abs(truth_nuPdg)==12 && truth_isCC==1 && truth_vtxInside==1 && T_BDTvars.numu_cc_flag>=0 && T_BDTvars.nue_score>7.0 && match_isFC==1)");
  T_eval->Project("h32","T_PFeval.truth_showerMomentum[3]*1000.","weight_cv*weight_spline*(abs(truth_nuPdg)==12 && truth_isCC==1 && truth_vtxInside==1 && T_BDTvars.numu_cc_flag>=0 && T_BDTvars.nue_score>7.0 )");

  T_eval->Project("h40","truth_nuEnergy-T_PFeval.truth_showerMomentum[3]*1000.","weight_cv*weight_spline*(abs(truth_nuPdg)==12 && truth_isCC==1 && truth_vtxInside==1)");
  T_eval->Project("h41","truth_nuEnergy-T_PFeval.truth_showerMomentum[3]*1000.","weight_cv*weight_spline*(abs(truth_nuPdg)==12 && truth_isCC==1 && truth_vtxInside==1 && T_BDTvars.numu_cc_flag>=0 && T_BDTvars.nue_score>7.0 && match_isFC==1)");
  T_eval->Project("h42","truth_nuEnergy-T_PFeval.truth_showerMomentum[3]*1000.","weight_cv*weight_spline*(abs(truth_nuPdg)==12 && truth_isCC==1 && truth_vtxInside==1 && T_BDTvars.numu_cc_flag>=0 && T_BDTvars.nue_score>7.0 )");


  
  TGraph *g1 = new TGraph();
  TGraph *g2 = new TGraph();

  TGraph *g3 = new TGraph();
  TGraph *g4 = new TGraph();

  TGraph *g5 = new TGraph();
  TGraph *g6 = new TGraph();
  

  for (Int_t i=0;i!=h10->GetNbinsX();i++){
    Double_t x = h10->GetBinCenter(i+1);
    Double_t y = h11->GetBinContent(i+1)/h10->GetBinContent(i+1);
    if (std::isnan(y)) y = 0;
    g1->SetPoint(i,x,y);
    y = h12->GetBinContent(i+1)/h10->GetBinContent(i+1);
    if (std::isnan(y)) y = 0;
    g2->SetPoint(i,x,y);


    x = h30->GetBinCenter(i+1);
    y = h31->GetBinContent(i+1)/h30->GetBinContent(i+1);
    if (std::isnan(y)) y = 0;
    g3->SetPoint(i,x,y);
    y = h32->GetBinContent(i+1)/h30->GetBinContent(i+1);
    if (std::isnan(y)) y = 0;
    g4->SetPoint(i,x,y);

    x = h40->GetBinCenter(i+1);
    y = h41->GetBinContent(i+1)/h40->GetBinContent(i+1);
    if (std::isnan(y)) y = 0;
    g5->SetPoint(i,x,y);
    y = h42->GetBinContent(i+1)/h40->GetBinContent(i+1);
    if (std::isnan(y)) y = 0;
    g6->SetPoint(i,x,y);
    
  }

  TCanvas *c1 = new TCanvas("c1","c1",1200,400);
  c1->Divide(3,1);
  
  c1->cd(1);
  g1->Draw("AL*");
  g2->Draw("*Lsame");
  g1->SetMarkerColor(1);
  g2->SetMarkerColor(2);
  g1->SetMarkerStyle(20);
  g2->SetMarkerStyle(20);
  g1->SetLineColor(1);
  g2->SetLineColor(2);
  g1->GetYaxis()->SetRangeUser(0,1);

  g1->GetXaxis()->SetTitle("E^{#nu}_{true} (MeV)");
  g1->SetTitle("BNB #nu_{e}CC Efficiency");
  TLegend *le1 = new TLegend(0.6,0.6,0.89,0.89);
  le1->AddEntry(g1,"FC","lp");
  le1->AddEntry(g2,"FC+PC","lp");
 
  le1->Draw();
  
  c1->cd(2);
  
  g3->Draw("AL*");
  g4->Draw("*Lsame");
  g3->SetMarkerColor(1);
  g4->SetMarkerColor(2);
  g3->SetMarkerStyle(20);
  g4->SetMarkerStyle(20);
  g3->SetLineColor(1);
  g4->SetLineColor(2);
  g3->GetYaxis()->SetRangeUser(0,1);
  g3->GetXaxis()->SetTitle("E^{e}_{true} (MeV)");

  c1->cd(3);
  
  g5->Draw("AL*");
  g6->Draw("*Lsame");
  g5->SetMarkerColor(1);
  g6->SetMarkerColor(2);
  g5->SetMarkerStyle(20);
  g6->SetMarkerStyle(20);
  g5->SetLineColor(1);
  g6->SetLineColor(2);
  g5->GetYaxis()->SetRangeUser(0,1);
  g5->GetXaxis()->SetTitle("#nu (MeV)");
  
  

  // h10->Draw();
  // h30->Scale(h10->GetSum()/h30->GetSum());
  // h30->SetLineColor(2);
  // h30->Draw("same");
  
}
