void plot_kdar(){
  TFile *file = new TFile("processed_checkout_rootfiles/checkout_data_numi_run1_morestat.root");
  TTree *T_eval = (TTree*)file->Get("wcpselection/T_eval");
  TTree *T_BDTvars = (TTree*)file->Get("wcpselection/T_BDTvars");
  TTree *T_KINEvars = (TTree*)file->Get("wcpselection/T_KINEvars");
  TTree *T_PFeval = (TTree*)file->Get("wcpselection/T_PFeval");
  T_eval->AddFriend(T_BDTvars,"T_BDTvars");
  T_eval->AddFriend(T_KINEvars,"T_KINEvars");
  T_eval->AddFriend(T_PFeval,"T_PFeval");

 

  TFile *file1 = new TFile("processed_checkout_rootfiles/checkout_prodgenie_numi_overlay_run1.root");
  TTree *T_eval1 = (TTree*)file1->Get("wcpselection/T_eval");
  TTree *T_BDTvars1 = (TTree*)file1->Get("wcpselection/T_BDTvars");
  TTree *T_KINEvars1 = (TTree*)file1->Get("wcpselection/T_KINEvars");
  T_eval1->AddFriend(T_BDTvars1,"T_BDTvars");
  T_eval1->AddFriend(T_KINEvars1,"T_KINEvars");

   TCanvas *c1 = new TCanvas("c1","c1",1200,600);
   c1->Divide(2,1);

   c1->cd(1);
  TH1F *h1 = new TH1F("h1","h1",50,150,300);
  T_eval->Project("h1","T_KINEvars.kine_reco_Enu","T_BDTvars.numu_score>0.0&&match_isFC==1");
  h1->Draw();
  h1->SetXTitle("E^{rec}_{#nu} (MeV)");
  h1->SetTitle("KDAR #nu (3 MeV per bin)");
  
   c1->cd(2);
  TH1F *h2 = new TH1F("h2","h2",50,0,300);
  TH1F *h3 = new TH1F("h3","h3",50,150,300);
  TH1F *h4 = new TH1F("h4","h4",50,150,300);
  //T_eval1->Project("h2","truth_nuEnergy"," weight_cv*weight_spline*(T_BDTvars.numu_cc_flag >=0 &&match_isFC==1 && truth_vtxInside==1&& truth_isCC==1 && abs(truth_nuPdg)==14)");
  //T_eval1->Project("h3","truth_nuEnergy"," weight_cv*weight_spline*(T_BDTvars.numu_score >=0 &&match_isFC==1 && truth_vtxInside==1&& truth_isCC==1 && abs(truth_nuPdg)==14)");
  T_eval1->Project("h2","T_KINEvars.kine_reco_Enu"," weight_cv*weight_spline*(T_BDTvars.numu_score >=0 &&match_isFC==1 && truth_vtxInside==1&& truth_isCC==1 && abs(truth_nuPdg)==14 && truth_nuEnergy>232 && truth_nuEnergy < 240)");
  // T_eval1->Project("h3","truth_nuEnergy"," weight_cv*weight_spline*(T_BDTvars.numu_cc_flag >=0 &&match_isFC==1 && truth_vtxInside==1&& truth_isCC==1 && abs(truth_nuPdg)==14)");
  // T_eval1->Project("h4","truth_nuEnergy"," weight_cv*weight_spline*(T_BDTvars.numu_score >=0.9 &&match_isFC==1 && truth_vtxInside==1&& truth_isCC==1 && abs(truth_nuPdg)==14)");
  h2->Draw();
  h2->SetLineColor(1);
  h3->Draw("same");
  h3->SetLineColor(2);
  h4->Draw("same");
  h4->SetLineColor(4);
  //  h2->GetYaxis()->SetRangeUser(-20,600);
  h2->SetTitle("Truth #nu energy (MeV)");
  //h2->SetTitle("Reco #nu energy (MeV)");

  TLegend *le1 = new TLegend(0.6,0.6,0.89,0.89);
  le1->AddEntry(h2,"all");
  le1->AddEntry(h3,"Generic #nu Det.");
  le1->AddEntry(h4,"Selected");
  le1->Draw();

  // TH1F *h100 = new TH1F("h100","h100",150,0,1500);
  // TH1F *h101 = new TH1F("h101","h101",150,0,1500);
  // TH1F *h102 = new TH1F("h102","h102",150,0,1500);
  // // T_eval1->Project("h100","truth_energyInside","weight_cv*weight_spline*(truth_vtxInside==1&& truth_isCC==1 && abs(truth_nuPdg)==14)");
  // //T_eval1->Project("h101","truth_energyInside","weight_cv*weight_spline*(T_BDTvars.numu_cc_flag >=0.0 && truth_vtxInside==1 && truth_isCC==1 && abs(truth_nuPdg)==14)");
  // T_eval1->Project("h100","truth_nuEnergy","weight_cv*weight_spline*(truth_vtxInside==1&& truth_isCC==1 && abs(truth_nuPdg)==14)");
  // T_eval1->Project("h101","truth_nuEnergy","weight_cv*weight_spline*(T_BDTvars.numu_cc_flag >=0.0 && truth_vtxInside==1 && truth_isCC==1 && abs(truth_nuPdg)==14)");
  // T_eval1->Project("h102","truth_nuEnergy","weight_cv*weight_spline*(match_isFC==1 && T_BDTvars.numu_score >=0 && T_BDTvars.numu_cc_flag >=0.0 && truth_vtxInside==1 && truth_isCC==1 && abs(truth_nuPdg)==14)");
  // h101->Divide(h100);
  // h102->Divide(h100);
  
  // h101->Draw();
  // h102->SetLineColor(2);
  // h102->Draw("same");
  // h101->GetYaxis()->SetRangeUser(0,1.);
  // le1->Draw();
  
  // T_eval->Draw("reco_muonMomentum[1]:reco_muonMomentum[2]","T_BDTvars.numu_score >0.0 &&match_isFC==1 && T_KINEvars.kine_reco_Enu>224 && T_KINEvars.kine_reco_Enu<232 && reco_muonMomentum[3]>0.","*");
}
