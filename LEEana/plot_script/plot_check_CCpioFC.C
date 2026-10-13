void plot_check_CCpioFC(int run = 1){
  TFile *file1, *file2, *file3, *file4, *file5, *file6, *file7;
  Double_t pot_1, pot_2, pot_3, pot_4, pot_5, pot_6, pot_data;

  if (run==1){
    file1 = new TFile("./processed_checkout_rootfiles/checkout_prodgenie_bnb_intrinsic_nue_overlay_run1.root");
    pot_1 = 5.63425e+22;
    file2 = new TFile("./processed_checkout_rootfiles/checkout_prodgenie_bnb_intrinsic_nue_overlay_LowE_run1.root");
    pot_2 = 2.68057e+23;
    file3 = new TFile("./processed_checkout_rootfiles/checkout_prodgenie_bnb_nu_overlay_run1.root");
    pot_3 = 6.29927e+20;
    file4 = new TFile("./processed_checkout_rootfiles/checkout_prodgenie_bnb_nu_overlay_LowE_run1.root");
    pot_4 = 7.7364e+20;
    file5 = new TFile("./processed_checkout_rootfiles/checkout_prodgenie_dirt_overlay_run1.root");
    pot_5 = 2.3922e+20;
    file6 = new TFile("./processed_checkout_rootfiles/checkout_data_extbnb_run1.root");
    pot_6 = 7.37133e+19;
    file7 = new TFile("./processed_checkout_rootfiles/run1_data_bnb_merge.root");
    pot_data = 4.47375e+19;
    

  }else if (run==3){
    file1 = new TFile("./processed_checkout_rootfiles/checkout_prodgenie_bnb_intrinsic_nue_overlay_run3.root");
    pot_1 = 4.14518e+22;
    file2 = new TFile("./processed_checkout_rootfiles/checkout_prodgenie_bnb_intrinsic_nue_overlay_LowE_run3.root");
    pot_2 = 3.11223e+23;
    file3 = new TFile("./processed_checkout_rootfiles/prev_checkout_run3_bnb_nu_POT1.22E21.root");
    pot_3 = 6.04258e+20;
    file4 = new TFile("./processed_checkout_rootfiles/checkout_prodgenie_bnb_nu_overlay_LowE_run3.root");
    pot_4 = 7.35109e+20;
    file5 = new TFile("./processed_checkout_rootfiles/checkout_prodgenie_dirt_overlay_run3.root");
    pot_5 = 1.79278e+20;
    file6 = new TFile("./processed_checkout_rootfiles/checkout_data_extbnb_run3.root");
    pot_6 = 2.95674e+20;
    file7 = new TFile("./processed_checkout_rootfiles/checkout_data_bnb_run3_1e19.root");
    pot_data = 7.05836e18;
  }
  
  // TFile *file3 = new TFile("processed_files/checkout_run1_bnb_nu_POT1.2E21.root"); Double_t pot_3 = 5.926e+20;
  // TFile *file4 = new TFile("processed_files/checkout_run1_bnb_nu_LowE_POT1.6E21.root"); Double_t pot_4 = 8.238e+20;
  // TFile *file5 = new TFile("processed_files/checkout_run1_dirt_POT2.4E20.root");
  // Double_t pot_5 = 2.442e+20;
  // TFile *file6 = new TFile("processed_files/checkout_run1_extbnb_C1_gt10_wcp_v00_14_00_POT1.2E20.root");
  // Double_t pot_6 = 6.000e+19;
  
  // Double_t pot_data = 4.5e19;

  TH1F *h10 = new TH1F("h10","h10",20,0,2000);
  TH1F *h11 = new TH1F("h11","h11",20,0,2000);

  TH1F *h20 = new TH1F("h20","h20",20,0,2000);
  TH1F *h21 = new TH1F("h21","h21",20,0,2000);

  TH1F *h30 = new TH1F("h30","h30",20,0,2000);
  TH1F *h31 = new TH1F("h31","h31",20,0,2000);

  TH1F *h40 = new TH1F("h40","h40",20,0,2000);
  TH1F *h41 = new TH1F("h41","h41",20,0,2000);

  TH1F *h50 = new TH1F("h50","h50",20,0,2000);
  TH1F *h51 = new TH1F("h51","h51",20,0,2000);
  
  TH1F *h60 = new TH1F("h60","h60",20,0,2000);
  TH1F *h61 = new TH1F("h61","h61",20,0,2000);
  

  TString pio_energy = "135 * (sqrt(2./(1-pow((T_KINEvars.kine_pio_energy_1 - T_KINEvars.kine_pio_energy_2)/(T_KINEvars.kine_pio_energy_1 + T_KINEvars.kine_pio_energy_2),2))/(1-cos(T_KINEvars.kine_pio_angle/180.*3.1415926)))-1)";
  
  {
    TTree *T_eval = (TTree*)file3->Get("wcpselection/T_eval");
    TTree *T_BDTvars = (TTree*)file3->Get("wcpselection/T_BDTvars");
    TTree *T_KINEvars = (TTree*)file3->Get("wcpselection/T_KINEvars");
    T_eval->AddFriend(T_BDTvars,"T_BDTvars");
    T_eval->AddFriend(T_KINEvars,"T_KINEvars");

    T_eval->Project("h10",pio_energy,"weight_cv*weight_spline*(T_BDTvars.numu_cc_flag>=0 &&match_isFC==1 && (truth_nuEnergy <=400 ) && T_BDTvars.numu_score>0.9  && T_BDTvars.nue_score <=7.0 && (T_KINEvars.kine_pio_flag==1 && T_KINEvars.kine_pio_vtx_dis < 9  && T_KINEvars.kine_pio_energy_1 > 40 && T_KINEvars.kine_pio_energy_2 > 25 && T_KINEvars.kine_pio_dis_1 < 110 && T_KINEvars.kine_pio_dis_2 < 120 && T_KINEvars.kine_pio_angle > 0 && T_KINEvars.kine_pio_angle < 174  && T_KINEvars.kine_pio_mass > 22 && T_KINEvars.kine_pio_mass < 300))");
    T_eval->Project("h11",pio_energy,"weight_cv*weight_cv*weight_spline*weight_spline*(T_BDTvars.numu_cc_flag>=0 &&match_isFC==1 && (truth_nuEnergy <=400 ) && T_BDTvars.numu_score>0.9  && T_BDTvars.nue_score <=7.0&& (T_KINEvars.kine_pio_flag==1 && T_KINEvars.kine_pio_vtx_dis < 9  && T_KINEvars.kine_pio_energy_1 > 40 && T_KINEvars.kine_pio_energy_2 > 25 && T_KINEvars.kine_pio_dis_1 < 110 && T_KINEvars.kine_pio_dis_2 < 120 && T_KINEvars.kine_pio_angle > 0 && T_KINEvars.kine_pio_angle < 174  && T_KINEvars.kine_pio_mass > 22 && T_KINEvars.kine_pio_mass < 300))");
    h10->Scale(pot_data/(pot_3+pot_4));
    h11->Scale(pow(pot_data/(pot_3+pot_4),2));

    T_eval->Project("h20",pio_energy,"weight_cv*weight_spline*(T_BDTvars.numu_cc_flag>=0 && match_isFC==1 && !(truth_nuEnergy <=400 ) && T_BDTvars.numu_score>0.9  && T_BDTvars.nue_score <=7.0&& (T_KINEvars.kine_pio_flag==1 && T_KINEvars.kine_pio_vtx_dis < 9  && T_KINEvars.kine_pio_energy_1 > 40 && T_KINEvars.kine_pio_energy_2 > 25 && T_KINEvars.kine_pio_dis_1 < 110 && T_KINEvars.kine_pio_dis_2 < 120 && T_KINEvars.kine_pio_angle > 0 && T_KINEvars.kine_pio_angle < 174  && T_KINEvars.kine_pio_mass > 22 && T_KINEvars.kine_pio_mass < 300))");
    T_eval->Project("h21",pio_energy,"weight_cv*weight_cv*weight_spline*weight_spline*(T_BDTvars.numu_cc_flag>=0  &&match_isFC==1 && !(truth_nuEnergy <=400 )&& T_BDTvars.numu_score>0.9  && T_BDTvars.nue_score <=7.0&& (T_KINEvars.kine_pio_flag==1 && T_KINEvars.kine_pio_vtx_dis < 9  && T_KINEvars.kine_pio_energy_1 > 40 && T_KINEvars.kine_pio_energy_2 > 25 && T_KINEvars.kine_pio_dis_1 < 110 && T_KINEvars.kine_pio_dis_2 < 120 && T_KINEvars.kine_pio_angle > 0 && T_KINEvars.kine_pio_angle < 174  && T_KINEvars.kine_pio_mass > 22 && T_KINEvars.kine_pio_mass < 300))");
    h20->Scale(pot_data/(pot_3));
    h21->Scale(pow(pot_data/(pot_3),2));
  }

  {
    TTree *T_eval = (TTree*)file4->Get("wcpselection/T_eval");
    TTree *T_BDTvars = (TTree*)file4->Get("wcpselection/T_BDTvars");
    TTree *T_KINEvars = (TTree*)file4->Get("wcpselection/T_KINEvars");
    T_eval->AddFriend(T_BDTvars,"T_BDTvars");
    T_eval->AddFriend(T_KINEvars,"T_KINEvars");

    T_eval->Project("h30",pio_energy,"weight_cv*weight_spline*(T_BDTvars.numu_cc_flag>=0 &&match_isFC==1 && (truth_nuEnergy <=400 ) && T_BDTvars.numu_score>0.9  && T_BDTvars.nue_score <=7.0&& (T_KINEvars.kine_pio_flag==1 && T_KINEvars.kine_pio_vtx_dis < 9  && T_KINEvars.kine_pio_energy_1 > 40 && T_KINEvars.kine_pio_energy_2 > 25 && T_KINEvars.kine_pio_dis_1 < 110 && T_KINEvars.kine_pio_dis_2 < 120 && T_KINEvars.kine_pio_angle > 0 && T_KINEvars.kine_pio_angle < 174  && T_KINEvars.kine_pio_mass > 22 && T_KINEvars.kine_pio_mass < 300))");
    T_eval->Project("h31",pio_energy,"weight_cv*weight_cv*weight_spline*weight_spline*(T_BDTvars.numu_cc_flag>=0 &&match_isFC==1 && (truth_nuEnergy <=400 ) && T_BDTvars.numu_score>0.9 && T_BDTvars.nue_score <=7.0&& (T_KINEvars.kine_pio_flag==1 && T_KINEvars.kine_pio_vtx_dis < 9  && T_KINEvars.kine_pio_energy_1 > 40 && T_KINEvars.kine_pio_energy_2 > 25 && T_KINEvars.kine_pio_dis_1 < 110 && T_KINEvars.kine_pio_dis_2 < 120 && T_KINEvars.kine_pio_angle > 0 && T_KINEvars.kine_pio_angle < 174  && T_KINEvars.kine_pio_mass > 22 && T_KINEvars.kine_pio_mass < 300))");
    h30->Scale(pot_data/(pot_3+pot_4));
    h31->Scale(pow(pot_data/(pot_3+pot_4),2));
  }

  {
     TTree *T_eval = (TTree*)file5->Get("wcpselection/T_eval");
    TTree *T_BDTvars = (TTree*)file5->Get("wcpselection/T_BDTvars");
    TTree *T_KINEvars = (TTree*)file5->Get("wcpselection/T_KINEvars");
    T_eval->AddFriend(T_BDTvars,"T_BDTvars");
    T_eval->AddFriend(T_KINEvars,"T_KINEvars");

    T_eval->Project("h40",pio_energy,"weight_cv*weight_spline*(T_BDTvars.numu_cc_flag>=0 &&match_isFC==1  && T_BDTvars.numu_score>0.9  && T_BDTvars.nue_score <=7.0&& (T_KINEvars.kine_pio_flag==1 && T_KINEvars.kine_pio_vtx_dis < 9  && T_KINEvars.kine_pio_energy_1 > 40 && T_KINEvars.kine_pio_energy_2 > 25 && T_KINEvars.kine_pio_dis_1 < 110 && T_KINEvars.kine_pio_dis_2 < 120 && T_KINEvars.kine_pio_angle > 0 && T_KINEvars.kine_pio_angle < 174  && T_KINEvars.kine_pio_mass > 22 && T_KINEvars.kine_pio_mass < 300))");
    T_eval->Project("h41",pio_energy,"weight_cv*weight_cv*weight_spline*weight_spline*(T_BDTvars.numu_cc_flag>=0 &&match_isFC==1 && T_BDTvars.numu_score>0.9 && T_BDTvars.nue_score <=7.0&& (T_KINEvars.kine_pio_flag==1 && T_KINEvars.kine_pio_vtx_dis < 9  && T_KINEvars.kine_pio_energy_1 > 40 && T_KINEvars.kine_pio_energy_2 > 25 && T_KINEvars.kine_pio_dis_1 < 110 && T_KINEvars.kine_pio_dis_2 < 120 && T_KINEvars.kine_pio_angle > 0 && T_KINEvars.kine_pio_angle < 174  && T_KINEvars.kine_pio_mass > 22 && T_KINEvars.kine_pio_mass < 300))");
    h40->Scale(pot_data/(pot_5));
    h41->Scale(pow(pot_data/(pot_5),2));
  }

  {
    TTree *T_eval = (TTree*)file6->Get("wcpselection/T_eval");
    TTree *T_BDTvars = (TTree*)file6->Get("wcpselection/T_BDTvars");
    TTree *T_KINEvars = (TTree*)file6->Get("wcpselection/T_KINEvars");
    T_eval->AddFriend(T_BDTvars,"T_BDTvars");
    T_eval->AddFriend(T_KINEvars,"T_KINEvars");

    T_eval->Project("h50",pio_energy,"(T_BDTvars.numu_cc_flag>=0 &&match_isFC==1  && T_BDTvars.numu_score>0.9  && T_BDTvars.nue_score <=7.0&& (T_KINEvars.kine_pio_flag==1 && T_KINEvars.kine_pio_vtx_dis < 9  && T_KINEvars.kine_pio_energy_1 > 40 && T_KINEvars.kine_pio_energy_2 > 25 && T_KINEvars.kine_pio_dis_1 < 110 && T_KINEvars.kine_pio_dis_2 < 120 && T_KINEvars.kine_pio_angle > 0 && T_KINEvars.kine_pio_angle < 174  && T_KINEvars.kine_pio_mass > 22 && T_KINEvars.kine_pio_mass < 300))");
    T_eval->Project("h51",pio_energy,"(T_BDTvars.numu_cc_flag>=0 &&match_isFC==1 && T_BDTvars.numu_score>0.9 && T_BDTvars.nue_score <=7.0&& (T_KINEvars.kine_pio_flag==1 && T_KINEvars.kine_pio_vtx_dis < 9  && T_KINEvars.kine_pio_energy_1 > 40 && T_KINEvars.kine_pio_energy_2 > 25 && T_KINEvars.kine_pio_dis_1 < 110 && T_KINEvars.kine_pio_dis_2 < 120 && T_KINEvars.kine_pio_angle > 0 && T_KINEvars.kine_pio_angle < 174  && T_KINEvars.kine_pio_mass > 22 && T_KINEvars.kine_pio_mass < 300))");
    h50->Scale(pot_data/(pot_6));
    h51->Scale(pow(pot_data/(pot_6),2));
  }

  {
    TTree *T_eval = (TTree*)file7->Get("wcpselection/T_eval");
    TTree *T_BDTvars = (TTree*)file7->Get("wcpselection/T_BDTvars");
    TTree *T_KINEvars = (TTree*)file7->Get("wcpselection/T_KINEvars");
    T_eval->AddFriend(T_BDTvars,"T_BDTvars");
    T_eval->AddFriend(T_KINEvars,"T_KINEvars");

    T_eval->Project("h60",pio_energy,"(T_BDTvars.numu_cc_flag>=0 &&match_isFC==1  && T_BDTvars.numu_score>0.9  && T_BDTvars.nue_score <=7.0&& (T_KINEvars.kine_pio_flag==1 && T_KINEvars.kine_pio_vtx_dis < 9  && T_KINEvars.kine_pio_energy_1 > 40 && T_KINEvars.kine_pio_energy_2 > 25 && T_KINEvars.kine_pio_dis_1 < 110 && T_KINEvars.kine_pio_dis_2 < 120 && T_KINEvars.kine_pio_angle > 0 && T_KINEvars.kine_pio_angle < 174  && T_KINEvars.kine_pio_mass > 22 && T_KINEvars.kine_pio_mass < 300))");
    T_eval->Project("h61",pio_energy,"(T_BDTvars.numu_cc_flag>=0 &&match_isFC==1 && T_BDTvars.numu_score>0.9 && T_BDTvars.nue_score <=7.0&& (T_KINEvars.kine_pio_flag==1 && T_KINEvars.kine_pio_vtx_dis < 9  && T_KINEvars.kine_pio_energy_1 > 40 && T_KINEvars.kine_pio_energy_2 > 25 && T_KINEvars.kine_pio_dis_1 < 110 && T_KINEvars.kine_pio_dis_2 < 120 && T_KINEvars.kine_pio_angle > 0 && T_KINEvars.kine_pio_angle < 174  && T_KINEvars.kine_pio_mass > 22 && T_KINEvars.kine_pio_mass < 300))");
    // h60->Scale(pot_data/(pot_6));
    // h61->Scale(pow(pot_data/(pot_6),2));
  }

  h10->Add(h20); h11->Add(h21);
  h10->Add(h30); h11->Add(h31);
  h10->Add(h40); h11->Add(h41);
  h10->Add(h50); h11->Add(h51);
  

  h60->Draw();
  h10->Draw("same");
  h60->SetLineColor(2);

  for (Int_t i=0;i!=h10->GetNbinsX()+1;i++){
    std::cout << i << " " << h10->GetBinContent(i+1) << " " << sqrt(h11->GetBinContent(i+1)) << " : " << h60->GetBinContent(i+1) << std::endl;
  }
  
}
