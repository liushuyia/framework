void plot_check_nueCCFC_Det(Int_t flag=1){
  
  TFile *file1; TFile *file3;

  if (flag == 1){
    file1 = new TFile("./hist_rootfiles/DetVar/WCP_checkout_prodgenie_bnb_intrinsic_nue_overlay_WCP_DetVar_LYDown_run3b.root");
    file3 = new TFile("./hist_rootfiles/DetVar/WCP_checkout_prodgenie_bnb_nu_overlay_WCP_DetVar_LYDown_run3b.root");
    std::cout << "LYDown" << std::endl;
  }else if (flag==2){
    file1 = new TFile("./hist_rootfiles/DetVar/WCP_checkout_prodgenie_bnb_intrinsic_nue_overlay_WCP_DetVar_LYRayleigh_run3b.root");
    file3 = new TFile("./hist_rootfiles/DetVar/WCP_checkout_prodgenie_bnb_nu_overlay_WCP_DetVar_LYRayleigh_run3b.root");
    std::cout << "LYRayleigh" << std::endl;
  }else if (flag==3){
    file1 = new TFile("./hist_rootfiles/DetVar/WCP_checkout_prodgenie_bnb_intrinsic_nue_overlay_WCP_DetVar_Recomb2_run3b.root");
    file3 = new TFile("./hist_rootfiles/DetVar/WCP_checkout_prodgenie_bnb_nu_overlay_WCP_DetVar_Recomb2_run3b.root");
    std::cout << "Recomb2" << std::endl;
  }else if (flag==4){
    file1 = new TFile("./hist_rootfiles/DetVar/WCP_checkout_prodgenie_bnb_intrinsic_nue_overlay_WCP_DetVar_SCE_run3b.root");
    file3 = new TFile("./hist_rootfiles/DetVar/WCP_checkout_prodgenie_bnb_nu_overlay_WCP_DetVar_SCE_run3b.root");
    std::cout << "SCE" << std::endl;
  }else if (flag==5){
    file1 = new TFile("./hist_rootfiles/DetVar/WCP_checkout_prodgenie_bnb_intrinsic_nue_overlay_WCP_DetVar_WireModdEdX_run3b.root");
    file3 = new TFile("./hist_rootfiles/DetVar/WCP_checkout_prodgenie_bnb_nu_overlay_WCP_DetVar_WireModdEdX_run3b.root");
    std::cout << "WireModdEdX" << std::endl;
  }else if (flag==6){
    file1 = new TFile("./hist_rootfiles/DetVar/WCP_checkout_prodgenie_bnb_intrinsic_nue_overlay_WCP_DetVar_WireModThetaXZ_run3b.root");
    file3 = new TFile("./hist_rootfiles/DetVar/WCP_checkout_prodgenie_bnb_nu_overlay_WCP_DetVar_WireModThetaXZ_run3b.root");
    std::cout << "WireModThetaXZ" << std::endl;
  }else if (flag==7){
    file1 = new TFile("./hist_rootfiles/DetVar/WCP_checkout_prodgenie_bnb_intrinsic_nue_overlay_WCP_DetVar_WireModThetaYZ_run3b.root");
    file3 = new TFile("./hist_rootfiles/DetVar/WCP_checkout_prodgenie_bnb_nu_overlay_WCP_DetVar_WireModThetaYZ_run3b.root");
    std::cout << "WireModThetaYZ" << std::endl;
  }else if (flag==8){
    file1 = new TFile("./hist_rootfiles/DetVar/WCP_checkout_prodgenie_bnb_intrinsic_nue_overlay_WCP_DetVar_WireModX_run3b.root");
    file3 = new TFile("./hist_rootfiles/DetVar/WCP_checkout_prodgenie_bnb_nu_overlay_WCP_DetVar_WireModX_run3b.root");
    std::cout << "WireModX" << std::endl;
  }else if (flag==9){
    file1 = new TFile("./hist_rootfiles/DetVar/WCP_checkout_prodgenie_bnb_intrinsic_nue_overlay_WCP_DetVar_WireModYZ_run3b.root");
    file3 = new TFile("./hist_rootfiles/DetVar/WCP_checkout_prodgenie_bnb_nu_overlay_WCP_DetVar_WireModYZ_run3b.root");
    std::cout << "WireModYZ" << std::endl;
  }


  Double_t pot_1 = 0;//3.80432e+22;
  Double_t pot_3 = 0;//8.8408e+19;
  TTree *T_pot1 = (TTree*)file1->Get("wcpselection/T_pot_cv");
  TTree *T_pot3 = (TTree*)file3->Get("wcpselection/T_pot_cv");
  Double_t pot_tor875;
  T_pot1->SetBranchAddress("pot_tor875",&pot_tor875);
  T_pot3->SetBranchAddress("pot_tor875",&pot_tor875);
  for (Int_t i=0;i!=T_pot1->GetEntries();i++){
    T_pot1->GetEntry(i);
    pot_1 += pot_tor875;
  }
   for (Int_t i=0;i!=T_pot3->GetEntries();i++){
    T_pot3->GetEntry(i);
    pot_3 += pot_tor875;
  }
   std::cout << pot_1 << " " << pot_3 << std::endl;
  
  
  Double_t pot_data = 5.0e19;

  TH1F *h10 = new TH1F("h10","h10",25,0,2500);
  TH1F *h11 = new TH1F("h11","h11",25,0,2500);

  TH1F *h20 = new TH1F("h20","h20",25,0,2500);
  TH1F *h21 = new TH1F("h21","h21",25,0,2500);

  TH1F *h30 = new TH1F("h30","h30",25,0,2500);
  TH1F *h31 = new TH1F("h31","h31",25,0,2500);

  
  TH1F *h40 = new TH1F("h40","h40",25,0,2500);
  TH1F *h41 = new TH1F("h41","h41",25,0,2500);

  TH1F *h50 = new TH1F("h50","h50",25,0,2500);
  TH1F *h51 = new TH1F("h51","h51",25,0,2500);

  TH1F *h60 = new TH1F("h60","h60",25,0,2500);
  TH1F *h61 = new TH1F("h61","h61",25,0,2500);



  {
    TTree *T_eval = (TTree*)file1->Get("wcpselection/T_eval_cv");
    TTree *T_BDTvars = (TTree*)file1->Get("wcpselection/T_BDTvars_cv");
    TTree *T_KINEvars = (TTree*)file1->Get("wcpselection/T_KINEvars_cv");
    T_eval->AddFriend(T_BDTvars);
    T_eval->AddFriend(T_KINEvars);

    TTree *T_eval_det = (TTree*)file1->Get("wcpselection/T_eval_det");
    TTree *T_BDTvars_det = (TTree*)file1->Get("wcpselection/T_BDTvars_det");
    TTree *T_KINEvars_det = (TTree*)file1->Get("wcpselection/T_KINEvars_det");
    T_eval_det->AddFriend(T_BDTvars_det);
    T_eval_det->AddFriend(T_KINEvars_det);

    //    std::cout << T_eval_cv->GetEntries() << " " << T_BDTvars_cv->GetEntries() << " " << T_KINEvars_cv->GetEntries() << std::endl;

    // T_eval->Scan("run:event:T_KINEvars_cv.kine_reco_Enu:weight_cv*weight_spline","weight_cv*weight_spline*(T_BDTvars_cv.numu_cc_flag>=0 && T_BDTvars_cv.nue_score>7 && match_isFC==1 && (truth_nuEnergy <=400 &&truth_vtxX > -1 && truth_vtxX <= 254.3 &&  truth_vtxY >-115.0 && truth_vtxY<=117.0 && truth_vtxZ > 0.6 && truth_vtxZ <=1036.4))");
    
    T_eval->Project("h10","T_KINEvars_cv.kine_reco_Enu","weight_cv*weight_spline*(T_BDTvars_cv.numu_cc_flag>=0 && T_BDTvars_cv.nue_score>7 && match_isFC==1 && (truth_nuEnergy <=400 &&truth_vtxX > -1 && truth_vtxX <= 254.3 &&  truth_vtxY >-115.0 && truth_vtxY<=117.0 && truth_vtxZ > 0.6 && truth_vtxZ <=1036.4))");
    T_eval->Project("h11","T_KINEvars_cv.kine_reco_Enu","weight_cv*weight_cv*weight_spline*weight_spline*(1+0*weight_lee)*(1+0*weight_lee)*(T_BDTvars_cv.numu_cc_flag>=0 && T_BDTvars_cv.nue_score>7 &&match_isFC==1 && (truth_nuEnergy <=400 &&truth_vtxX > -1 && truth_vtxX <= 254.3 &&  truth_vtxY >-115.0 && truth_vtxY<=117.0 && truth_vtxZ > 0.6 && truth_vtxZ <=1036.4))");

    h10->Scale(pot_data/(pot_1));
    h11->Scale(pow(pot_data/(pot_1),2));

    //    std::cout << h10->GetEntries() << " " << h10->GetSum()  << std::endl;
    
    
    T_eval->Project("h20","T_KINEvars_cv.kine_reco_Enu","weight_cv*weight_spline*(1+0*weight_lee)*(T_BDTvars_cv.numu_cc_flag>=0 && T_BDTvars_cv.nue_score>7 &&match_isFC==1 && !(truth_nuEnergy <=400 &&truth_vtxX > -1 && truth_vtxX <= 254.3 &&  truth_vtxY >-115.0 && truth_vtxY<=117.0 && truth_vtxZ > 0.6 && truth_vtxZ <=1036.4))");
    T_eval->Project("h21","T_KINEvars_cv.kine_reco_Enu","weight_cv*weight_cv*weight_spline*weight_spline*(1+0*weight_lee)*(1+0*weight_lee)*(T_BDTvars_cv.numu_cc_flag>=0 && T_BDTvars_cv.nue_score>7 &&match_isFC==1 && !(truth_nuEnergy <=400 &&truth_vtxX > -1 && truth_vtxX <= 254.3 &&  truth_vtxY >-115.0 && truth_vtxY<=117.0 && truth_vtxZ > 0.6 && truth_vtxZ <=1036.4))");
    h20->Scale(pot_data/(pot_1));
    h21->Scale(pow(pot_data/(pot_1),2));

    T_eval_det->Project("h30","T_KINEvars_det.kine_reco_Enu","weight_cv*weight_spline*(T_BDTvars_det.numu_cc_flag>=0 && T_BDTvars_det.nue_score>7 && match_isFC==1 && (truth_nuEnergy <=400 &&truth_vtxX > -1 && truth_vtxX <= 254.3 &&  truth_vtxY >-115.0 && truth_vtxY<=117.0 && truth_vtxZ > 0.6 && truth_vtxZ <=1036.4))");
    T_eval_det->Project("h31","T_KINEvars_det.kine_reco_Enu","weight_cv*weight_spline*(1+0*weight_lee)*(T_BDTvars_det.numu_cc_flag>=0 && T_BDTvars_det.nue_score>7 &&match_isFC==1 && !(truth_nuEnergy <=400 &&truth_vtxX > -1 && truth_vtxX <= 254.3 &&  truth_vtxY >-115.0 && truth_vtxY<=117.0 && truth_vtxZ > 0.6 && truth_vtxZ <=1036.4))");
    h30->Scale(pot_data/(pot_1));
    h31->Scale(pot_data/(pot_1));
    
  }

 

  {
    TTree *T_eval = (TTree*)file3->Get("wcpselection/T_eval_cv");
    TTree *T_BDTvars = (TTree*)file3->Get("wcpselection/T_BDTvars_cv");
    TTree *T_KINEvars = (TTree*)file3->Get("wcpselection/T_KINEvars_cv");
    T_eval->AddFriend(T_BDTvars);
    T_eval->AddFriend(T_KINEvars);

    TTree *T_eval_det = (TTree*)file3->Get("wcpselection/T_eval_det");
    TTree *T_BDTvars_det = (TTree*)file3->Get("wcpselection/T_BDTvars_det");
    TTree *T_KINEvars_det = (TTree*)file3->Get("wcpselection/T_KINEvars_det");
    T_eval_det->AddFriend(T_BDTvars_det);
    T_eval_det->AddFriend(T_KINEvars_det);

    // T_eval->Project("h40","T_KINEvars_cv.kine_reco_Enu","weight_cv*weight_spline*(numu_cc_flag>=0 && nue_score>7 &&match_isFC==1 && (truth_nuEnergy <=400 ) && !(truth_isCC==1 && abs(truth_nuPdg)==12 && truth_vtxX > -1 && truth_vtxX <= 254.3 &&  truth_vtxY >-115.0 && truth_vtxY<=117.0 && truth_vtxZ > 0.6 && truth_vtxZ <=1036.4))");
    // T_eval->Project("h41","T_KINEvars_cv.kine_reco_Enu","weight_cv*weight_cv*weight_spline*weight_spline*(numu_cc_flag>=0 && nue_score>7 &&match_isFC==1 && (truth_nuEnergy <=400 )&& !(truth_isCC==1 && abs(truth_nuPdg)==12 && truth_vtxX > -1 && truth_vtxX <= 254.3 &&  truth_vtxY >-115.0 && truth_vtxY<=117.0 && truth_vtxZ > 0.6 && truth_vtxZ <=1036.4))");
    // h40->Scale(pot_data/(pot_3));
    // h41->Scale(pow(pot_data/(pot_3),2));

    // T_eval->Project("h50","T_KINEvars_cv.kine_reco_Enu","weight_cv*weight_spline*(numu_cc_flag>=0 && nue_score>7 &&match_isFC==1 && !(truth_nuEnergy <=400 )&& !(truth_isCC==1 && abs(truth_nuPdg)==12 && truth_vtxX > -1 && truth_vtxX <= 254.3 &&  truth_vtxY >-115.0 && truth_vtxY<=117.0 && truth_vtxZ > 0.6 && truth_vtxZ <=1036.4))");
    // T_eval->Project("h51","T_KINEvars_cv.kine_reco_Enu","weight_cv*weight_cv*weight_spline*weight_spline*(numu_cc_flag>=0 && nue_score>7 &&match_isFC==1 && !(truth_nuEnergy <=400 )&& !(truth_isCC==1 && abs(truth_nuPdg)==12 && truth_vtxX > -1 && truth_vtxX <= 254.3 &&  truth_vtxY >-115.0 && truth_vtxY<=117.0 && truth_vtxZ > 0.6 && truth_vtxZ <=1036.4))");
    // h50->Scale(pot_data/(pot_3));
    // h51->Scale(pow(pot_data/(pot_3),2));


    //  T_eval_det->Project("h60","T_KINEvars_det.kine_reco_Enu","weight_cv*weight_spline*(numu_cc_flag>=0 && nue_score>7 &&match_isFC==1 && (truth_nuEnergy <=400 ) && !(truth_isCC==1 && abs(truth_nuPdg)==12 && truth_vtxX > -1 && truth_vtxX <= 254.3 &&  truth_vtxY >-115.0 && truth_vtxY<=117.0 && truth_vtxZ > 0.6 && truth_vtxZ <=1036.4))");
    //  T_eval_det->Project("h61","T_KINEvars_det.kine_reco_Enu","weight_cv*weight_spline*(numu_cc_flag>=0 && nue_score>7 &&match_isFC==1 && !(truth_nuEnergy <=400 )&& !(truth_isCC==1 && abs(truth_nuPdg)==12 && truth_vtxX > -1 && truth_vtxX <= 254.3 &&  truth_vtxY >-115.0 && truth_vtxY<=117.0 && truth_vtxZ > 0.6 && truth_vtxZ <=1036.4))");
    //  h60->Scale(pot_data/pot_3);
    //  h61->Scale(pot_data/pot_3);
  }

  

  // std::cout << h10->GetBinContent(3) << " " << h11->GetBinContent(3) << std::endl;
  // std::cout << h20->GetBinContent(3) << " " << h21->GetBinContent(3) << std::endl;
  // //std::cout << h30->GetBinContent(3) << " " << h31->GetBinContent(3) << std::endl;
  // std::cout << h40->GetBinContent(3) << " " << h41->GetBinContent(3) << std::endl;
  // std::cout << h50->GetBinContent(3) << " " << h51->GetBinContent(3) << std::endl;
  //  std::cout << h60->GetBinContent(3) << " " << h61->GetBinContent(3) << std::endl;
  //std::cout << h70->GetBinContent(3) << " " << h71->GetBinContent(3) << std::endl;
  //  std::cout << h80->GetBinContent(3) << " " << h81->GetBinContent(3) << std::endl;

  

  h10->Add(h20);  h11->Add(h21);
  //  h10->Add(h30);  h11->Add(h31);
  h10->Add(h40);  h11->Add(h41);
  h10->Add(h50);  h11->Add(h51);

  h30->Add(h31);
  h30->Add(h60);
  h30->Add(h61);
  //h10->Add(h70);  h11->Add(h71);
  //h10->Add(h80);  h11->Add(h81);

  h10->Draw();
  h30->Draw("same");
  h30->SetLineColor(2);

  for (Int_t i=0;i!=h10->GetNbinsX()+1;i++){
    std::cout << i << " " << h10->GetBinContent(i+1) << " " << sqrt(h11->GetBinContent(i+1)) << std::endl;
  }

  
}
