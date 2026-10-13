void plot_check_numuCCFC_xf(int run =1){
  TFile *file1;
  TString option;
  if (run==1){
    file1 = new TFile("./processed_checkout_rootfiles/prodgenie_bnb_nu_overlay_run1/expskin_FluxUnisim.root");
    option = "expskin_FluxUnisim";
  }else if(run==2){
    file1 = new TFile("processed_checkout_rootfiles/prodgenie_bnb_nu_overlay_run1/horncurrent_FluxUnisim.root");
    option = "horncurrent_FluxUnisim";
  }else if(run==3){
    file1 = new TFile("processed_checkout_rootfiles/prodgenie_bnb_nu_overlay_run1/kminus_PrimaryHadronNormalization.root");
    option = "kminus_PrimaryHadronNormalization";
  }else if(run==4){
    file1 = new TFile("processed_checkout_rootfiles/prodgenie_bnb_nu_overlay_run1/kplus_PrimaryHadronFeynmanScaling.root");
    option = "kplus_PrimaryHadronFeynmanScaling";
  }else if(run==5){
    file1 = new TFile("processed_checkout_rootfiles/prodgenie_bnb_nu_overlay_run1/kzero_PrimaryHadronSanfordWang.root");
    option = "kzero_PrimaryHadronSanfordWang";
  }else if(run==6){
    file1 = new TFile("processed_checkout_rootfiles/prodgenie_bnb_nu_overlay_run1/nucleoninexsec_FluxUnisim.root");
    option = "nucleoninexsec_FluxUnisim";
  }else if(run==7){
    file1 = new TFile("processed_checkout_rootfiles/prodgenie_bnb_nu_overlay_run1/nucleonqexsec_FluxUnisim.root");
    option = "nucleonqexsec_FluxUnisim";
  }else if(run==8){
    file1 = new TFile("processed_checkout_rootfiles/prodgenie_bnb_nu_overlay_run1/nucleontotxsec_FluxUnisim.root");
    option = "nucleontotxsec_FluxUnisim";
  }else if(run==9){
    file1 = new TFile("processed_checkout_rootfiles/prodgenie_bnb_nu_overlay_run1/piminus_PrimaryHadronSWCentralSplineVariation.root");
    option = "piminus_PrimaryHadronSWCentralSplineVariation";
  }else if(run==10){
    file1 = new TFile("processed_checkout_rootfiles/prodgenie_bnb_nu_overlay_run1/pioninexsec_FluxUnisim.root");
    option = "pioninexsec_FluxUnisim";
  }else if(run==11){
    file1 = new TFile("processed_checkout_rootfiles/prodgenie_bnb_nu_overlay_run1/pionqexsec_FluxUnisim.root");
    option = "pionqexsec_FluxUnisim";
  }else if(run==12){
    file1 = new TFile("processed_checkout_rootfiles/prodgenie_bnb_nu_overlay_run1/piontotxsec_FluxUnisim.root");
    option = "piontotxsec_FluxUnisim";
  }else if(run==13){
    file1 = new TFile("processed_checkout_rootfiles/prodgenie_bnb_nu_overlay_run1/piplus_PrimaryHadronSWCentralSplineVariation.root");
    option = "piplus_PrimaryHadronSWCentralSplineVariation";
  }else if(run==14){
    file1 = new TFile("processed_checkout_rootfiles/prodgenie_bnb_nu_overlay_run1/UBGenieFluxSmallUni.root");
    option = "All_UBGenie";
  }

   Double_t pot_3 = 0;
  TTree *T_pot3 = (TTree*)file1->Get("wcpselection/T_pot");
  Double_t pot_tor875;
  T_pot3->SetBranchAddress("pot_tor875",&pot_tor875);
  for (Int_t i=0;i!=T_pot3->GetEntries();i++){
    T_pot3->GetEntry(i);
    pot_3 += pot_tor875;
  }
  std::cout << pot_3 << std::endl;

  double pot_data = 6e19;

  TH1F *h10 = new TH1F("h10","h10",25,0,2500);
  TH1F *h11 = new TH1F("h11","h11",25,0,2500);

  TH1F *h20 = new TH1F("h20","h20",25,0,2500);
  TH1F *h21 = new TH1F("h21","h21",25,0,2500);

  TString cut = "T_eval.weight_cv*T_eval.weight_spline*(T_BDTvars.numu_cc_flag>=0 &&match_isFC==1  && T_BDTvars.numu_score>0.9 && T_BDTvars.nue_score<=7 && !(T_KINEvars.kine_pio_flag==1 && T_KINEvars.kine_pio_vtx_dis < 9  && T_KINEvars.kine_pio_energy_1 > 40 && T_KINEvars.kine_pio_energy_2 > 25 && T_KINEvars.kine_pio_dis_1 < 110 && T_KINEvars.kine_pio_dis_2 < 120 && T_KINEvars.kine_pio_angle > 0 && T_KINEvars.kine_pio_angle < 174  && T_KINEvars.kine_pio_mass > 22 && T_KINEvars.kine_pio_mass < 300) && T_eval.weight_change==0)";
  {
    TTree *T_eval = (TTree*)file1->Get("wcpselection/T_eval");
    TTree *T_BDTvars = (TTree*)file1->Get("wcpselection/T_BDTvars");
    TTree *T_KINEvars = (TTree*)file1->Get("wcpselection/T_KINEvars");
    TTree *T_weight = (TTree*)file1->Get("wcpselection/T_weight");    
    T_eval->AddFriend(T_BDTvars,"T_BDTvars");
    T_eval->AddFriend(T_KINEvars,"T_KINEvars");
    T_eval->AddFriend(T_weight,"T_weight");
    T_eval->Project("h10","T_KINEvars.kine_reco_Enu",cut);
    
    h10->Scale(pot_data/(pot_3));
    
  }

  h10->Draw();

  for (Int_t i=0;i!=h10->GetNbinsX();i++){
    std::cout << i << " " << h10->GetBinContent(i+1) << std::endl;
  }

  TMatrixD *cov = new TMatrixD(26,26);
  TMatrixD *cov1 = new TMatrixD(26,26);
  {
    TTree *T_eval = (TTree*)file1->Get("wcpselection/T_eval");
    TTree *T_BDTvars = (TTree*)file1->Get("wcpselection/T_BDTvars");
    TTree *T_KINEvars = (TTree*)file1->Get("wcpselection/T_KINEvars");
    TTree *T_weight = (TTree*)file1->Get("wcpselection/T_weight");
    T_eval->AddFriend(T_BDTvars,"T_BDTvars");
    T_eval->AddFriend(T_KINEvars,"T_KINEvars");
    T_eval->AddFriend(T_weight,"T_weight");

    h11->Reset();
    for (Int_t i=0;i!=10;i++){
      std::cout << i << std::endl;
      TString cut1 = "T_weight." + option + Form("[%d]*",i) + cut;     
      if (run==14)
	cut1 = "1./T_eval.weight_cv * " + cut1;
      //std::cout << cut1 << std::endl;
      T_eval->Project("h11","T_KINEvars.kine_reco_Enu",cut1);
      h11->Scale(pot_data/(pot_3));    

      h11->Add(h10,-1);
      for (Int_t ii=0;ii!=26;ii++){
      	for (Int_t j=0;j!=26;j++){
      	  (*cov)(ii,j) += h11->GetBinContent(ii+1) * h11->GetBinContent(j+1);
      	}
      }
    }
  }

  (*cov) *= 1./10.;

  h11->SetLineColor(2);
  h11->Draw("same");

  TFile *file2 = new TFile(Form("./hist_rootfiles/XsFlux/cov_%d.root",run));
  TMatrixD *cov2 = (TMatrixD*)file2->Get(Form("cov_xf_mat_%d",run));
  for (Int_t i=0;i!=26;i++){
    for (Int_t j=0;j!=26;j++){
      (*cov1)(i,j) = (*cov2)(26+26+i,26+26+j);
    }
  }
  
  TCanvas *c1 = new TCanvas("c1","c1",1200,800);
  c1->Divide(2,2);
  
  c1->cd(1);
  cov->Draw("COLZ");
  //cov->SetName("this one");
  c1->cd(2);
  cov1->Draw("COLZ");
  //cov1->SetName("program");
  c1->cd(3);
  h10->Draw();

  TH1F *h50 = (TH1F*)file2->Get("pred_covch_3");
  h50->Draw("same");
  h50->SetLineColor(2);
  
  for (Int_t i=0;i!=26;i++){
    std::cout << (*cov)(i,i) << " " << (*cov1)(i,i) << std::endl;
  }
  
}
