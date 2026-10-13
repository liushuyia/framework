void plot_eff_NCpi0(){
  // TFile *file1 =  new TFile("./processed_checkout_rootfiles/checkout_prodgenie_bnb_nu_overlay_run1.root");
  // TFile *file1 =  new TFile("checkout_prodgenie_bnb_nu_overlay_run1.root");
  TFile *file1 =  new TFile("../../checkout_rootfiles_correct_bdt/checkout_prodgenie_bnb_nu_overlay_run1.root");
  //TFile *file1 =  new TFile("./processed_checkout_rootfiles/DetVar/WCPcheckout_prodgenie_bnb_nu_overlay_WCP_DetVar_CV_run3b.root");
  // TFile *file3 =  new TFile("./processed_checkout_rootfiles/checkout_prodgenie_bnb_nu_overlay_run3.root");
  // TFile *file3 =  new TFile("checkout_prodgenie_bnb_nu_overlay_run3.root");
  TFile *file3 =  new TFile("../../checkout_rootfiles_correct_bdt/checkout_prodgenie_bnb_nu_overlay_run3.root");
  //TFile *file3 =  new TFile("./processed_checkout_rootfiles/DetVar/WCPcheckout_prodgenie_bnb_nu_overlay_WCP_DetVar_CV_run3b.root");
  //TFile *file3 =  new TFile("./processed_checkout_rootfiles/DetVar/WCPcheckout_prodgenie_bnb_nu_overlay_WCP_DetVar_LYDown_run3b.root");

  // TH1F *h10 = new TH1F("h10","h10",12,0,1200);
  // TH1F *h11 = new TH1F("h11","h11",12,0,1200);
  // TH1F *h12 = new TH1F("h12","h12",12,0,1200);
  // TH1F *h30 = new TH1F("h30","h30",12,0,1200);
  // TH1F *h31 = new TH1F("h31","h31",12,0,1200);
  // TH1F *h32 = new TH1F("h32","h32",12,0,1200);

  TH1F *h10 = new TH1F("h10","h10",8,0,800);
  TH1F *h11 = new TH1F("h11","h11",8,0,800);
  TH1F *h12 = new TH1F("h12","h12",8,0,800);
  TH1F *h30 = new TH1F("h30","h30",8,0,800);
  TH1F *h31 = new TH1F("h31","h31",8,0,800);
  TH1F *h32 = new TH1F("h32","h32",8,0,800);


  TTree *T_eval = (TTree*)file1->Get("wcpselection/T_eval");
  TTree *T_PFeval = (TTree*)file1->Get("wcpselection/T_PFeval");
  TTree *T_BDTvars = (TTree*)file1->Get("wcpselection/T_BDTvars");
  TTree *T_KINEvars = (TTree*)file1->Get("wcpselection/T_KINEvars");
  T_eval->AddFriend(T_BDTvars,"T_BDTvars");
  T_eval->AddFriend(T_KINEvars,"T_KINEvars");
  T_eval->AddFriend(T_PFeval,"T_PFeval");
  T_eval->SetAlias("alpha","fabs(T_PFeval.truth_pio_energy_1 - T_PFeval.truth_pio_energy_2)/(T_PFeval.truth_pio_energy_1 + T_PFeval.truth_pio_energy_2)");
  T_eval->SetAlias("truth_pi0Energy","135.0 * (sqrt(2./(1-alpha*alpha)/(1-cos(T_PFeval.truth_pio_angle/180.*3.1415926)))-1)");
  T_eval->SetAlias("flag_numuCC","T_BDTvars.numu_cc_flag>=0 && T_BDTvars.numu_score>0.9");
  T_eval->SetAlias("flag_cc_pi0","(T_KINEvars.kine_pio_flag==1 && T_KINEvars.kine_pio_vtx_dis < 9 ) && T_KINEvars.kine_pio_energy_1 > 40 && T_KINEvars.kine_pio_energy_2 > 25 && T_KINEvars.kine_pio_dis_1 < 110 && T_KINEvars.kine_pio_dis_2 < 120 && T_KINEvars.kine_pio_angle > 0 && T_KINEvars.kine_pio_angle < 174 && T_KINEvars.kine_pio_mass > 22 && T_KINEvars.kine_pio_mass < 300");
  T_eval->SetAlias("flag_nueCC","T_BDTvars.numu_cc_flag >=0 && T_BDTvars.nue_score > 7.0");
  T_eval->SetAlias("flag_NC","(!T_BDTvars.cosmict_flag) && T_BDTvars.numu_score < 0");
  T_eval->SetAlias("flag_pi0","(T_KINEvars.kine_pio_flag==1 && T_KINEvars.kine_pio_vtx_dis < 9 || T_KINEvars.kine_pio_flag==2) && T_KINEvars.kine_pio_energy_1 > 40 && T_KINEvars.kine_pio_energy_2 > 25 && T_KINEvars.kine_pio_dis_1 < 110 && T_KINEvars.kine_pio_dis_2 < 120 && T_KINEvars.kine_pio_angle > 0 && T_KINEvars.kine_pio_angle < 174 && T_KINEvars.kine_pio_mass > 22 && T_KINEvars.kine_pio_mass < 300");


  TTree *T_eval3 = (TTree*)file3->Get("wcpselection/T_eval");
  TTree *T_PFeval3 = (TTree*)file3->Get("wcpselection/T_PFeval");
  TTree *T_BDTvars3 = (TTree*)file3->Get("wcpselection/T_BDTvars");
  TTree *T_KINEvars3 = (TTree*)file3->Get("wcpselection/T_KINEvars");
  T_eval3->AddFriend(T_BDTvars3,"T_BDTvars");
  T_eval3->AddFriend(T_KINEvars3,"T_KINEvars");
  T_eval3->AddFriend(T_PFeval3,"T_PFeval");
  T_eval3->SetAlias("alpha","fabs(truth_pio_energy_1 - truth_pio_energy_2)/(truth_pio_energy_1 + truth_pio_energy_2)");
  T_eval3->SetAlias("truth_pi0Energy","135.0 * (sqrt(2./(1-alpha*alpha)/(1-cos(truth_pio_angle/180.*3.1415926)))-1)");
  T_eval3->SetAlias("flag_numuCC","T_BDTvars.numu_cc_flag>=0 && T_BDTvars.numu_score>0.9");
  T_eval3->SetAlias("flag_cc_pi0","(T_KINEvars.kine_pio_flag==1 && T_KINEvars.kine_pio_vtx_dis < 9 ) && T_KINEvars.kine_pio_energy_1 > 40 && T_KINEvars.kine_pio_energy_2 > 25 && T_KINEvars.kine_pio_dis_1 < 110 && T_KINEvars.kine_pio_dis_2 < 120 && T_KINEvars.kine_pio_angle > 0 && T_KINEvars.kine_pio_angle < 174 && T_KINEvars.kine_pio_mass > 22 && T_KINEvars.kine_pio_mass < 300");
  T_eval3->SetAlias("flag_nueCC","T_BDTvars.numu_cc_flag >=0 && T_BDTvars.nue_score > 7.0");
  T_eval3->SetAlias("flag_NC","(!T_BDTvars.cosmict_flag) && T_BDTvars.numu_score < 0");
  T_eval3->SetAlias("flag_pi0","(T_KINEvars.kine_pio_flag==1 && T_KINEvars.kine_pio_vtx_dis < 9 || T_KINEvars.kine_pio_flag==2) && T_KINEvars.kine_pio_energy_1 > 40 && T_KINEvars.kine_pio_energy_2 > 25 && T_KINEvars.kine_pio_dis_1 < 110 && T_KINEvars.kine_pio_dis_2 < 120 && T_KINEvars.kine_pio_angle > 0 && T_KINEvars.kine_pio_angle < 174 && T_KINEvars.kine_pio_mass > 22 && T_KINEvars.kine_pio_mass < 300");


  
  // NCpi0
  T_eval->Project("h10","truth_pi0Energy","weight_cv*weight_spline*(abs(truth_nuPdg)==14 && truth_isCC==1 && truth_vtxInside==1 && truth_NprimPio>0)");
  T_eval->Project("h12","truth_pi0Energy","weight_cv*weight_spline*(abs(truth_nuPdg)==14 && truth_isCC==1 && truth_vtxInside==1 && truth_NprimPio>0 && flag_NC && flag_pi0 && (!flag_nueCC) )");

  T_eval3->Project("h30","truth_pi0Energy","weight_cv*weight_spline*(abs(truth_nuPdg)==14 && truth_isCC==1 && truth_vtxInside==1 && truth_NprimPio>0)");
  T_eval3->Project("h32","truth_pi0Energy","weight_cv*weight_spline*(abs(truth_nuPdg)==14 && truth_isCC==1 && truth_vtxInside==1 && truth_NprimPio>0 && flag_NC && flag_pi0 && (!flag_nueCC) )");


   TGraph *g1 = new TGraph();
   TGraph *g2 = new TGraph();
 
   for (Int_t i=0;i!=h10->GetNbinsX();i++){
     Double_t x = h10->GetBinCenter(i+1);
     Double_t y = h12->GetBinContent(i+1)/h10->GetBinContent(i+1);
     if (std::isnan(y)) y = 0;
     g1->SetPoint(i,x,y);
 
     y = h32->GetBinContent(i+1)/h30->GetBinContent(i+1);
     if (std::isnan(y)) y = 0;
     g2->SetPoint(i,x,y);
     
   }
   g1->Draw("AL*");
   // g2->Draw("*Lsame");
   g1->SetMarkerColor(1);
   g2->SetMarkerColor(2);
   g1->SetMarkerStyle(20);
   g2->SetMarkerStyle(20);
   g1->SetLineColor(1);
   g2->SetLineColor(2);
   g1->GetYaxis()->SetRangeUser(0,1);
 
   // g3->Draw("*Lsame");
   // g4->Draw("*Lsame");
   // g3->SetMarkerStyle(21);
   // g4->SetMarkerStyle(21);
   // g3->SetLineColor(4);
   // g4->SetLineColor(6);
   // g3->SetMarkerColor(4);
   // g4->SetMarkerColor(6);
   // 
   // g3->SetLineStyle(2);
   // g4->SetLineStyle(2);
   
   
   g1->GetXaxis()->SetTitle("E^{#pi0}_{true} (MeV)");
   g1->SetTitle("NC#pi^{0} Efficiency");
   TLegend *le1 = new TLegend(0.6,0.6,0.89,0.89);
   // le1->AddEntry(g1,"FC","lp");
   // le1->AddEntry(g2,"FC+PC","lp");
   // le1->AddEntry(g1,"run 1 FC","lp");
   // le1->AddEntry(g2,"run 1 FC+PC","lp");
   // le1->AddEntry(g3,"run 3 FC","lp");
   // le1->AddEntry(g4,"run 3 FC+PC","lp");
   // le1->Draw();
  
  
}
