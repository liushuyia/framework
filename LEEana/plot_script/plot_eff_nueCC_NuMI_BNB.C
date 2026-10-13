void plot_eff_nueCC_NuMI_BNB(){
  TFile *file1 =  new TFile("./processed_checkout_rootfiles/checkout_prodgenie_bnb_intrinsic_nue_overlay_run123_all.root");
  TFile *file3 = new TFile("./processed_checkout_rootfiles/checkout_prodgenie_numi_intrinsic_nue_overlay_run1.root");
 
  TH1F *h10 = new TH1F("h10","h10",8,0,6000);
  TH1F *h11 = new TH1F("h11","h11",8,0,6000);
  TH1F *h12 = new TH1F("h12","h12",8,0,6000);

  TH1F *h30 = new TH1F("h30","h30",8,0,6000);
  TH1F *h31 = new TH1F("h31","h31",8,0,6000);
  TH1F *h32 = new TH1F("h32","h32",8,0,6000);

 

  TTree *T_eval = (TTree*)file1->Get("wcpselection/T_eval");
  TTree *T_BDTvars = (TTree*)file1->Get("wcpselection/T_BDTvars");
  TTree *T_KINEvars = (TTree*)file1->Get("wcpselection/T_KINEvars");
  TTree *T_PFeval = (TTree*)file1->Get("wcpselection/T_PFeval");
  T_eval->AddFriend(T_BDTvars,"T_BDTvars");
  T_eval->AddFriend(T_KINEvars,"T_KINEvars");
  T_eval->AddFriend(T_PFeval,"T_PFeval");

  
  TTree *T_eval3 = (TTree*)file3->Get("wcpselection/T_eval");
  TTree *T_BDTvars3 = (TTree*)file3->Get("wcpselection/T_BDTvars");
  TTree *T_KINEvars3 = (TTree*)file3->Get("wcpselection/T_KINEvars");
  TTree *T_PFeval3 = (TTree*)file3->Get("wcpselection/T_PFeval");
  T_eval3->AddFriend(T_BDTvars3,"T_BDTvars");
  T_eval3->AddFriend(T_KINEvars3,"T_KINEvars");
  T_eval3->AddFriend(T_PFeval3,"T_PFeval");
  
 
  
  T_eval->Project("h10","truth_nuEnergy","weight_cv*weight_spline*(abs(truth_nuPdg)==12 && truth_isCC==1 && truth_vtxInside==1)");
  T_eval->Project("h11","truth_nuEnergy","weight_cv*weight_spline*(abs(truth_nuPdg)==12 && truth_isCC==1 && truth_vtxInside==1 && T_BDTvars.numu_cc_flag>=0 && T_BDTvars.nue_score>7.0 && match_isFC==1)");
  T_eval->Project("h12","truth_nuEnergy","weight_cv*weight_spline*(abs(truth_nuPdg)==12 && truth_isCC==1 && truth_vtxInside==1 && T_BDTvars.numu_cc_flag>=0 && T_BDTvars.nue_score>7.0 )");


  T_eval3->Project("h30","truth_nuEnergy","weight_cv*weight_spline*(abs(truth_nuPdg)==12 && truth_isCC==1 && truth_vtxInside==1)");
  T_eval3->Project("h31","truth_nuEnergy","weight_cv*weight_spline*(abs(truth_nuPdg)==12 && truth_isCC==1 && truth_vtxInside==1 && T_BDTvars.numu_cc_flag>=0 && T_BDTvars.nue_score>7.0 && match_isFC==1)");
  T_eval3->Project("h32","truth_nuEnergy","weight_cv*weight_spline*(abs(truth_nuPdg)==12 && truth_isCC==1 && truth_vtxInside==1 && T_BDTvars.numu_cc_flag>=0 && T_BDTvars.nue_score>7.0 )");

 


  
  TGraph *g1 = new TGraph();
  TGraph *g2 = new TGraph();

  TGraph *g3 = new TGraph();
  TGraph *g4 = new TGraph();

  

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

  
    
  }

  TCanvas *c1 = new TCanvas("c1","c1",800,600);
 
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
  g1->SetTitle(" #nu_{e}CC Efficiency");
 
  
  
  
  g3->Draw("L*same");
  g4->Draw("*Lsame");
  g3->SetMarkerColor(4);
  g4->SetMarkerColor(6);
  g3->SetMarkerStyle(21);
  g4->SetMarkerStyle(21);
  g3->SetLineColor(4);
  g4->SetLineColor(6);
  
   TLegend *le1 = new TLegend(0.6,0.6,0.89,0.89);
  le1->AddEntry(g1,"BNB FC","lp");
  le1->AddEntry(g2,"BNB FC+PC","lp");
  le1->AddEntry(g3,"NuMI FC","lp");
  le1->AddEntry(g4,"NuMI FC+PC","lp");
  le1->Draw();
 
  

  // h10->Draw();
  // h30->Scale(h10->GetSum()/h30->GetSum());
  // h30->SetLineColor(2);
  // h30->Draw("same");
  
}
