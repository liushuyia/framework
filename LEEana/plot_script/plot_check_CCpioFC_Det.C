void plot_check_CCpioFC_Det(int flag=1){
  TFile *file3;

   if (flag==1){
    file3 = new TFile("./hist_rootfiles/DetVar/WCP_checkout_prodgenie_bnb_nu_overlay_WCP_DetVar_LYDown_run3b.root");
    std::cout << "LYDown" << std::endl;
  }else if (flag==2){
    file3 = new TFile("./hist_rootfiles/DetVar/WCP_checkout_prodgenie_bnb_nu_overlay_WCP_DetVar_LYRayleigh_run3b.root");
    std::cout << "LYRayleigh" << std::endl;
  }else if (flag==3){
    file3 = new TFile("./hist_rootfiles/DetVar/WCP_checkout_prodgenie_bnb_nu_overlay_WCP_DetVar_Recomb2_run3b.root");
    std::cout << "Recomb2" << std::endl;
  }else if (flag==4){
    file3 = new TFile("./hist_rootfiles/DetVar/WCP_checkout_prodgenie_bnb_nu_overlay_WCP_DetVar_SCE_run3b.root");
    std::cout << "SCE" << std::endl;
  }else if (flag==5){
    file3 = new TFile("./hist_rootfiles/DetVar/WCP_checkout_prodgenie_bnb_nu_overlay_WCP_DetVar_WireModdEdX_run3b.root");
    std::cout << "WireModdEdX" << std::endl;
  }else if (flag==6){
    file3 = new TFile("./hist_rootfiles/DetVar/WCP_checkout_prodgenie_bnb_nu_overlay_WCP_DetVar_WireModThetaXZ_run3b.root");
    std::cout << "WireModThetaXZ" << std::endl;
  }else if (flag==7){
    file3 = new TFile("./hist_rootfiles/DetVar/WCP_checkout_prodgenie_bnb_nu_overlay_WCP_DetVar_WireModThetaYZ_run3b.root");
    std::cout << "WireModThetaYZ" << std::endl;
  }else if (flag==8){
    file3 = new TFile("./hist_rootfiles/DetVar/WCP_checkout_prodgenie_bnb_nu_overlay_WCP_DetVar_WireModX_run3b.root");
    std::cout << "WireModX" << std::endl;
  }else if (flag==9){
    file3 = new TFile("./hist_rootfiles/DetVar/WCP_checkout_prodgenie_bnb_nu_overlay_WCP_DetVar_WireModYZ_run3b.root");
    std::cout << "WireModYZ" << std::endl;
  }

  Double_t pot_3 = 0;//8.8408e+19;
  TTree *T_pot3 = (TTree*)file3->Get("wcpselection/T_pot_cv");
  Double_t pot_tor875;
  T_pot3->SetBranchAddress("pot_tor875",&pot_tor875);
  for (Int_t i=0;i!=T_pot3->GetEntries();i++){
    T_pot3->GetEntry(i);
    pot_3 += pot_tor875;
  }
  std::cout << pot_3 << std::endl;
 
  
  Double_t pot_data = 5.0e19;

  TH1F *h10 = new TH1F("h10","h10",20,0,2000);
  TH1F *h11 = new TH1F("h11","h11",20,0,2000);

  TH1F *h20 = new TH1F("h20","h20",20,0,2000);
  TH1F *h21 = new TH1F("h21","h21",20,0,2000);

  TH1F *h30 = new TH1F("h30","h30",20,0,2000);
  TH1F *h31 = new TH1F("h31","h31",20,0,2000);

  TH1F *h40 = new TH1F("h40","h40",20,0,2000);
  TH1F *h41 = new TH1F("h41","h41",20,0,2000);


  TString pio_energy = "135 * (sqrt(2./(1-pow((T_KINEvars.kine_pio_energy_1 - T_KINEvars.kine_pio_energy_2)/(T_KINEvars.kine_pio_energy_1 + T_KINEvars.kine_pio_energy_2),2))/(1-cos(T_KINEvars.kine_pio_angle/180.*3.1415926)))-1)";
  
  {
    TTree *T_eval = (TTree*)file3->Get("wcpselection/T_eval_cv");
    TTree *T_BDTvars = (TTree*)file3->Get("wcpselection/T_BDTvars_cv");
    TTree *T_KINEvars = (TTree*)file3->Get("wcpselection/T_KINEvars_cv");
    T_eval->AddFriend(T_BDTvars,"T_BDTvars");
    T_eval->AddFriend(T_KINEvars,"T_KINEvars");

    TTree *T_eval_det = (TTree*)file3->Get("wcpselection/T_eval_det");
    TTree *T_BDTvars_det = (TTree*)file3->Get("wcpselection/T_BDTvars_det");
    TTree *T_KINEvars_det = (TTree*)file3->Get("wcpselection/T_KINEvars_det");
    T_eval_det->AddFriend(T_BDTvars_det,"T_BDTvars");
    T_eval_det->AddFriend(T_KINEvars_det,"T_KINEvars");

    T_eval->Project("h10",pio_energy,"weight_cv*weight_spline*(T_BDTvars.numu_cc_flag>=0 &&match_isFC==1 && (truth_nuEnergy <=400 ) && T_BDTvars.numu_score>0.9  && (!T_BDTvars.nue_score<7) && (T_KINEvars.kine_pio_flag==1 && T_KINEvars.kine_pio_vtx_dis < 9  && T_KINEvars.kine_pio_energy_1 > 40 && T_KINEvars.kine_pio_energy_2 > 25 && T_KINEvars.kine_pio_dis_1 < 110 && T_KINEvars.kine_pio_dis_2 < 120 && T_KINEvars.kine_pio_angle > 0 && T_KINEvars.kine_pio_angle < 174  && T_KINEvars.kine_pio_mass > 22 && T_KINEvars.kine_pio_mass < 300))");
    T_eval->Project("h11",pio_energy,"weight_cv*weight_cv*weight_spline*weight_spline*(T_BDTvars.numu_cc_flag>=0 &&match_isFC==1 && (truth_nuEnergy <=400 ) && T_BDTvars.numu_score>0.9  && (!T_BDTvars.nue_score<7) && (T_KINEvars.kine_pio_flag==1 && T_KINEvars.kine_pio_vtx_dis < 9  && T_KINEvars.kine_pio_energy_1 > 40 && T_KINEvars.kine_pio_energy_2 > 25 && T_KINEvars.kine_pio_dis_1 < 110 && T_KINEvars.kine_pio_dis_2 < 120 && T_KINEvars.kine_pio_angle > 0 && T_KINEvars.kine_pio_angle < 174  && T_KINEvars.kine_pio_mass > 22 && T_KINEvars.kine_pio_mass < 300))");
    h10->Scale(pot_data/(pot_3));
    h11->Scale(pow(pot_data/(pot_3),2));

    T_eval->Project("h20",pio_energy,"weight_cv*weight_spline*(T_BDTvars.numu_cc_flag>=0 && match_isFC==1 && !(truth_nuEnergy <=400 ) && T_BDTvars.numu_score>0.9  && (!T_BDTvars.nue_score<7)&& (T_KINEvars.kine_pio_flag==1 && T_KINEvars.kine_pio_vtx_dis < 9  && T_KINEvars.kine_pio_energy_1 > 40 && T_KINEvars.kine_pio_energy_2 > 25 && T_KINEvars.kine_pio_dis_1 < 110 && T_KINEvars.kine_pio_dis_2 < 120 && T_KINEvars.kine_pio_angle > 0 && T_KINEvars.kine_pio_angle < 174  && T_KINEvars.kine_pio_mass > 22 && T_KINEvars.kine_pio_mass < 300))");
    T_eval->Project("h21",pio_energy,"weight_cv*weight_cv*weight_spline*weight_spline*(T_BDTvars.numu_cc_flag>=0  &&match_isFC==1 && !(truth_nuEnergy <=400 )&& T_BDTvars.numu_score>0.9  && (!T_BDTvars.nue_score<7) && (T_KINEvars.kine_pio_flag==1 && T_KINEvars.kine_pio_vtx_dis < 9  && T_KINEvars.kine_pio_energy_1 > 40 && T_KINEvars.kine_pio_energy_2 > 25 && T_KINEvars.kine_pio_dis_1 < 110 && T_KINEvars.kine_pio_dis_2 < 120 && T_KINEvars.kine_pio_angle > 0 && T_KINEvars.kine_pio_angle < 174  && T_KINEvars.kine_pio_mass > 22 && T_KINEvars.kine_pio_mass < 300))");
    h20->Scale(pot_data/(pot_3));
    h21->Scale(pow(pot_data/(pot_3),2));


    T_eval_det->Project("h30",pio_energy,"weight_cv*weight_spline*(T_BDTvars.numu_cc_flag>=0 &&match_isFC==1 && (truth_nuEnergy <=400 ) && T_BDTvars.numu_score>0.9  && (!T_BDTvars.nue_score<7) && (T_KINEvars.kine_pio_flag==1 && T_KINEvars.kine_pio_vtx_dis < 9  && T_KINEvars.kine_pio_energy_1 > 40 && T_KINEvars.kine_pio_energy_2 > 25 && T_KINEvars.kine_pio_dis_1 < 110 && T_KINEvars.kine_pio_dis_2 < 120 && T_KINEvars.kine_pio_angle > 0 && T_KINEvars.kine_pio_angle < 174  && T_KINEvars.kine_pio_mass > 22 && T_KINEvars.kine_pio_mass < 300))");
    T_eval_det->Project("h40",pio_energy,"weight_cv*weight_spline*(T_BDTvars.numu_cc_flag>=0 && match_isFC==1 && !(truth_nuEnergy <=400 ) && T_BDTvars.numu_score>0.9  && (!T_BDTvars.nue_score<7)&& (T_KINEvars.kine_pio_flag==1 && T_KINEvars.kine_pio_vtx_dis < 9  && T_KINEvars.kine_pio_energy_1 > 40 && T_KINEvars.kine_pio_energy_2 > 25 && T_KINEvars.kine_pio_dis_1 < 110 && T_KINEvars.kine_pio_dis_2 < 120 && T_KINEvars.kine_pio_angle > 0 && T_KINEvars.kine_pio_angle < 174  && T_KINEvars.kine_pio_mass > 22 && T_KINEvars.kine_pio_mass < 300))");
    h30->Scale(pot_data/(pot_3));
    h40->Scale(pot_data/(pot_3));
  }

 

  h10->Add(h20); h11->Add(h21);
  h30->Add(h40);
  
  
  h10->Draw();
  h30->Draw("same");
  h30->SetLineColor(2);

  for (Int_t i=0;i!=h10->GetNbinsX()+1;i++){
    //    std::cout << i << " " << h10->GetBinContent(i+1) << " " << sqrt(h11->GetBinContent(i+1)) << std::endl;
  }
  
}
