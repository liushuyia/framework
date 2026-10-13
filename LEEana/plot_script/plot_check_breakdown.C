void plot_check_breakdown(){
  TFile *file = new TFile("processed_checkout_rootfiles/checkout_prodgenie_bnb_nu_overlay_run1.root");
  TTree *T_eval = (TTree*)file->Get("wcpselection/T_eval");
  TTree *T_BDTvars = (TTree*)file->Get("wcpselection/T_BDTvars");
  TTree *T_KINEvars = (TTree*)file->Get("wcpselection/T_KINEvars");
  TTree *T_PFeval = (TTree*)file->Get("wcpselection/T_PFeval");
  
  TH1F *h10 = new TH1F("h10","h10",25,0,2500);
  TH1F *h11 = new TH1F("h11","h11",25,0,2500);
  TH1F *h12 = new TH1F("h12","h12",25,0,2500);
  TH1F *h13 = new TH1F("h13","h13",25,0,2500);
  TH1F *h14 = new TH1F("h14","h14",25,0,2500);
  TH1F *h15 = new TH1F("h15","h15",25,0,2500);
  TH1F *h16 = new TH1F("h16","h16",25,0,2500);
  TH1F *h17 = new TH1F("h17","h17",25,0,2500);

  T_eval->AddFriend(T_BDTvars,"T_BDTvars");
  T_eval->AddFriend(T_KINEvars,"T_KINEvars");
  T_eval->AddFriend(T_PFeval,"T_PFeval");
  
  T_eval->Project("h10","T_KINEvars.kine_reco_Enu","weight_cv*weight_spline*(T_BDTvars.numu_score>=0.9 &&match_isFC!=1 && T_BDTvars.nue_score<=7 && !(T_KINEvars.kine_pio_flag==1 && T_KINEvars.kine_pio_vtx_dis < 9  && T_KINEvars.kine_pio_energy_1 > 40 && T_KINEvars.kine_pio_energy_2 > 25 && T_KINEvars.kine_pio_dis_1 < 110 && T_KINEvars.kine_pio_dis_2 < 120 && T_KINEvars.kine_pio_angle > 0 && T_KINEvars.kine_pio_angle < 174  && T_KINEvars.kine_pio_mass > 22 && T_KINEvars.kine_pio_mass < 300))");
  std::cout << "total: " << h10->GetSum() << std::endl;

  // no match
  T_eval->Project("h11","T_KINEvars.kine_reco_Enu","weight_cv*weight_spline*(T_BDTvars.numu_score>=0.9 &&match_isFC!=1 && T_BDTvars.nue_score<=7 && !(T_KINEvars.kine_pio_flag==1 && T_KINEvars.kine_pio_vtx_dis < 9  && T_KINEvars.kine_pio_energy_1 > 40 && T_KINEvars.kine_pio_energy_2 > 25 && T_KINEvars.kine_pio_dis_1 < 110 && T_KINEvars.kine_pio_dis_2 < 120 && T_KINEvars.kine_pio_angle > 0 && T_KINEvars.kine_pio_angle < 174  && T_KINEvars.kine_pio_mass > 22 && T_KINEvars.kine_pio_mass < 300) && match_completeness_energy/truth_energyInside >=0.1 && truth_vtxInside==0)");
  std::cout << "outFV: " << h11->GetSum() << std::endl;

  // nueCCinFV
  T_eval->Project("h12","T_KINEvars.kine_reco_Enu","weight_cv*weight_spline*(T_BDTvars.numu_score>=0.9 &&match_isFC!=1 && T_BDTvars.nue_score<=7 && !(T_KINEvars.kine_pio_flag==1 && T_KINEvars.kine_pio_vtx_dis < 9  && T_KINEvars.kine_pio_energy_1 > 40 && T_KINEvars.kine_pio_energy_2 > 25 && T_KINEvars.kine_pio_dis_1 < 110 && T_KINEvars.kine_pio_dis_2 < 120 && T_KINEvars.kine_pio_angle > 0 && T_KINEvars.kine_pio_angle < 174  && T_KINEvars.kine_pio_mass > 22 && T_KINEvars.kine_pio_mass < 300) && match_completeness_energy/truth_energyInside >=0.1 && truth_vtxInside==1 && abs(truth_nuPdg)==12 && truth_isCC==1)");
  std::cout << "nueCCinFV: " << h12->GetSum() << std::endl;

  // numuCCinFV
   T_eval->Project("h13","T_KINEvars.kine_reco_Enu","weight_cv*weight_spline*(T_BDTvars.numu_score>=0.9 &&match_isFC!=1 && T_BDTvars.nue_score<=7 && !(T_KINEvars.kine_pio_flag==1 && T_KINEvars.kine_pio_vtx_dis < 9  && T_KINEvars.kine_pio_energy_1 > 40 && T_KINEvars.kine_pio_energy_2 > 25 && T_KINEvars.kine_pio_dis_1 < 110 && T_KINEvars.kine_pio_dis_2 < 120 && T_KINEvars.kine_pio_angle > 0 && T_KINEvars.kine_pio_angle < 174  && T_KINEvars.kine_pio_mass > 22 && T_KINEvars.kine_pio_mass < 300) && match_completeness_energy/truth_energyInside >=0.1 && truth_vtxInside==1 && abs(truth_nuPdg)==14 && truth_isCC==1 && T_PFeval.truth_NprimPio==0)");
  std::cout << "numuCCinFV: " << h13->GetSum() << std::endl;

  //CCpi0inFV
   T_eval->Project("h14","T_KINEvars.kine_reco_Enu","weight_cv*weight_spline*(T_BDTvars.numu_score>=0.9 &&match_isFC!=1 && T_BDTvars.nue_score<=7 && !(T_KINEvars.kine_pio_flag==1 && T_KINEvars.kine_pio_vtx_dis < 9  && T_KINEvars.kine_pio_energy_1 > 40 && T_KINEvars.kine_pio_energy_2 > 25 && T_KINEvars.kine_pio_dis_1 < 110 && T_KINEvars.kine_pio_dis_2 < 120 && T_KINEvars.kine_pio_angle > 0 && T_KINEvars.kine_pio_angle < 174  && T_KINEvars.kine_pio_mass > 22 && T_KINEvars.kine_pio_mass < 300) && match_completeness_energy/truth_energyInside >=0.1 && truth_vtxInside==1 && abs(truth_nuPdg)==14 && truth_isCC==1 && T_PFeval.truth_NprimPio>0)");
  std::cout << "CCpi0inFV: " << h14->GetSum() << std::endl;

  
  //NCinFV
    T_eval->Project("h15","T_KINEvars.kine_reco_Enu","weight_cv*weight_spline*(T_BDTvars.numu_score>=0.9 &&match_isFC!=1 && T_BDTvars.nue_score<=7 && !(T_KINEvars.kine_pio_flag==1 && T_KINEvars.kine_pio_vtx_dis < 9  && T_KINEvars.kine_pio_energy_1 > 40 && T_KINEvars.kine_pio_energy_2 > 25 && T_KINEvars.kine_pio_dis_1 < 110 && T_KINEvars.kine_pio_dis_2 < 120 && T_KINEvars.kine_pio_angle > 0 && T_KINEvars.kine_pio_angle < 174  && T_KINEvars.kine_pio_mass > 22 && T_KINEvars.kine_pio_mass < 300) && match_completeness_energy/truth_energyInside >=0.1 && truth_vtxInside==1 && truth_isCC==0 && T_PFeval.truth_NprimPio==0)");
  std::cout << "NCinFV: " << h15->GetSum() << std::endl;
  

  //NCpi0inFV
  T_eval->Project("h16","T_KINEvars.kine_reco_Enu","weight_cv*weight_spline*(T_BDTvars.numu_score>=0.9 &&match_isFC!=1 && T_BDTvars.nue_score<=7 && !(T_KINEvars.kine_pio_flag==1 && T_KINEvars.kine_pio_vtx_dis < 9  && T_KINEvars.kine_pio_energy_1 > 40 && T_KINEvars.kine_pio_energy_2 > 25 && T_KINEvars.kine_pio_dis_1 < 110 && T_KINEvars.kine_pio_dis_2 < 120 && T_KINEvars.kine_pio_angle > 0 && T_KINEvars.kine_pio_angle < 174  && T_KINEvars.kine_pio_mass > 22 && T_KINEvars.kine_pio_mass < 300) && match_completeness_energy/truth_energyInside >=0.1 && truth_vtxInside==1 && truth_isCC==0 && T_PFeval.truth_NprimPio>0)");
  std::cout << "NCpi0inFV: " << h16->GetSum() << std::endl;
  

  //badmatch
  T_eval->Project("h17","T_KINEvars.kine_reco_Enu","weight_cv*weight_spline*(T_BDTvars.numu_score>=0.9 &&match_isFC!=1 && T_BDTvars.nue_score<=7 && !(T_KINEvars.kine_pio_flag==1 && T_KINEvars.kine_pio_vtx_dis < 9  && T_KINEvars.kine_pio_energy_1 > 40 && T_KINEvars.kine_pio_energy_2 > 25 && T_KINEvars.kine_pio_dis_1 < 110 && T_KINEvars.kine_pio_dis_2 < 120 && T_KINEvars.kine_pio_angle > 0 && T_KINEvars.kine_pio_angle < 174  && T_KINEvars.kine_pio_mass > 22 && T_KINEvars.kine_pio_mass < 300) && match_completeness_energy <0.1 * truth_energyInside )");
  std::cout << "badmatch: " << h17->GetSum() << std::endl;

  
}


