void plot_check_far_sideband(int run=1){
  TFile *file1, *file2, *file3, *file4, *file5, *file6, *file7;
  Double_t pot_1, pot_2, pot_3, pot_4, pot_5, pot_6, pot_data;

  if (run ==1){
    //file1 = new TFile("./processed_checkout_rootfiles/checkout_prodgenie_bnb_intrinsic_nue_overlay_run1.root");
    // pot_1 = 5.63425e+22;
    //    file2 = new TFile("./processed_checkout_rootfiles/checkout_prodgenie_bnb_intrinsic_nue_overlay_LowE_run1.root");
    //pot_2 = 2.68057e+23;
    file3 = new TFile("./processed_checkout_rootfiles/checkout_prodgenie_bnb_nu_overlay_run1.root");
    pot_3 = 6.40718e+20;
    //file4 = new TFile("./processed_checkout_rootfiles/checkout_prodgenie_bnb_nu_overlay_LowE_run1.root");
    //pot_4 = 7.7364e+20;
    file5 = new TFile("./processed_checkout_rootfiles/old_weights/checkout_prodgenie_dirt_overlay_run1.root");
    pot_5 = 2.3922e+20;
    file6 = new TFile("./processed_checkout_rootfiles/checkout_data_extbnb_run1.root");
    pot_6 = 2.20176e+20;
    file7 = new TFile("./processed_checkout_rootfiles/checkout_data_bnb_run1_5e19.root");
    pot_data = 4.42723e+19;
  }

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

  TH1F *h70 = new TH1F("h70","h70",25,0,2500);
  TH1F *h71 = new TH1F("h71","h71",25,0,2500);

  {
    TTree *T_eval = (TTree*)file3->Get("wcpselection/T_eval");
    TTree *T_BDTvars = (TTree*)file3->Get("wcpselection/T_BDTvars");
    TTree *T_KINEvars = (TTree*)file3->Get("wcpselection/T_KINEvars");
    T_eval->AddFriend(T_BDTvars,"T_BDTvars");
    T_eval->AddFriend(T_KINEvars,"T_KINEvars");

    T_eval->Project("h10","T_KINEvars.kine_reco_Enu","weight_cv*weight_spline*(T_BDTvars.nue_score<=0 || T_KINEvars.kine_reco_Enu >=800)");
    h10->Scale(pot_data/(pot_3));
    

  }

  {
    TTree *T_eval = (TTree*)file5->Get("wcpselection/T_eval");
    TTree *T_BDTvars = (TTree*)file5->Get("wcpselection/T_BDTvars");
    TTree *T_KINEvars = (TTree*)file5->Get("wcpselection/T_KINEvars");
    T_eval->AddFriend(T_BDTvars,"T_BDTvars");
    T_eval->AddFriend(T_KINEvars,"T_KINEvars");

    T_eval->Project("h20","T_KINEvars.kine_reco_Enu","weight_cv*weight_spline*(T_BDTvars.nue_score<=0 || T_KINEvars.kine_reco_Enu >=800)");
    h20->Scale(pot_data/(pot_5));
    

  }

  {
    TTree *T_eval = (TTree*)file6->Get("wcpselection/T_eval");
    TTree *T_BDTvars = (TTree*)file6->Get("wcpselection/T_BDTvars");
    TTree *T_KINEvars = (TTree*)file6->Get("wcpselection/T_KINEvars");
    T_eval->AddFriend(T_BDTvars,"T_BDTvars");
    T_eval->AddFriend(T_KINEvars,"T_KINEvars");

    T_eval->Project("h30","T_KINEvars.kine_reco_Enu","(T_BDTvars.nue_score<=0 || T_KINEvars.kine_reco_Enu >=800)");
    h30->Scale(pot_data/(pot_6));
    

  }


  {
    TTree *T_eval = (TTree*)file7->Get("wcpselection/T_eval");
    TTree *T_BDTvars = (TTree*)file7->Get("wcpselection/T_BDTvars");
    TTree *T_KINEvars = (TTree*)file7->Get("wcpselection/T_KINEvars");
    T_eval->AddFriend(T_BDTvars,"T_BDTvars");
    T_eval->AddFriend(T_KINEvars,"T_KINEvars");

    T_eval->Project("h40","T_KINEvars.kine_reco_Enu","(T_BDTvars.nue_score<=0 || T_KINEvars.kine_reco_Enu >=800)");
    // h40->Scale(pot_data/(pot_data));
    

  }
  THStack *hs = new THStack("hs","");
  hs->Add(h10);
  h10->SetFillColor(kGreen);
  hs->Add(h20);
  h20->SetFillColor(kRed);
  hs->Add(h30);
  h30->SetFillColor(kBlue);
  //  h40->Draw("");
  hs->Draw("hist");
  

  h40->Draw("same");
  h40->SetLineColor(6);
  h40->SetLineWidth(2);
    }
