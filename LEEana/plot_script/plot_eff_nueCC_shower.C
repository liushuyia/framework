double delta_binomial(double a, double b){ // error of a/b
  // return a/b * sqrt(1./a+1./b);
  double p = a/b;
  // cout << a << " " << b << " " << sqrt(b*p*(1-p))/b << endl;
  return sqrt(b*p*(1-p)) / b;
}

void plot_eff_nueCC_shower(){
  TFile *file1 =  new TFile("./processed_checkout_rootfiles/checkout_prodgenie_bnb_intrinsic_nue_overlay_run123_all.root");
  // TFile *file1 =  new TFile("./processed_checkout_rootfiles/checkout_prodgenie_bnb_intrinsic_nue_overlay_run1.root");
  // TFile *file1 =  new TFile("../../checkout_rootfiles_correct_bdt/checkout_prodgenie_bnb_intrinsic_nue_overlay_run1.root");
  // TFile *file3 =  new TFile("./processed_checkout_rootfiles/checkout_prodgenie_bnb_nu_overlay_run3.root");
  // TFile *file3 =  new TFile("../../checkout_rootfiles_correct_bdt/checkout_prodgenie_bnb_intrinsic_nue_overlay_run3.root");


  TH2F *h10 = new TH2F("h10","h10",14,0,1.400, 10,-1,1); // Eshower vs. costh
  TH2F *h11 = new TH2F("h11","h11",14,0,1.400, 10,-1,1);
  TH2F *h12 = new TH2F("h12","h12",14,0,1.400, 10,-1,1);
  TH2F *h30 = new TH2F("h30","h30",14,0,1.400, 10,-1,1);
  TH2F *h31 = new TH2F("h31","h31",14,0,1.400, 10,-1,1);
  TH2F *h32 = new TH2F("h32","h32",14,0,1.400, 10,-1,1);

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
  T_eval->SetAlias("truth_showerCosth","truth_showerMomentum[2]/sqrt(pow(truth_showerMomentum[2],2) + pow(truth_showerMomentum[0],2))");


  // TTree *T_eval3 = (TTree*)file3->Get("wcpselection/T_eval");
  // TTree *T_PFeval3 = (TTree*)file3->Get("wcpselection/T_PFeval");
  // TTree *T_BDTvars3 = (TTree*)file3->Get("wcpselection/T_BDTvars");
  // TTree *T_KINEvars3 = (TTree*)file3->Get("wcpselection/T_KINEvars");
  // T_eval3->AddFriend(T_BDTvars3,"T_BDTvars");
  // T_eval3->AddFriend(T_KINEvars3,"T_KINEvars");
  // T_eval3->AddFriend(T_PFeval3,"T_PFeval");
  // T_eval3->SetAlias("alpha","fabs(truth_pio_energy_1 - truth_pio_energy_2)/(truth_pio_energy_1 + truth_pio_energy_2)");
  // T_eval3->SetAlias("truth_pi0Energy","135.0 * (sqrt(2./(1-alpha*alpha)/(1-cos(truth_pio_angle/180.*3.1415926)))-1)");
  // T_eval3->SetAlias("flag_numuCC","T_BDTvars.numu_cc_flag>=0 && T_BDTvars.numu_score>0.9");
  // T_eval3->SetAlias("flag_cc_pi0","(T_KINEvars.kine_pio_flag==1 && T_KINEvars.kine_pio_vtx_dis < 9 ) && T_KINEvars.kine_pio_energy_1 > 40 && T_KINEvars.kine_pio_energy_2 > 25 && T_KINEvars.kine_pio_dis_1 < 110 && T_KINEvars.kine_pio_dis_2 < 120 && T_KINEvars.kine_pio_angle > 0 && T_KINEvars.kine_pio_angle < 174 && T_KINEvars.kine_pio_mass > 22 && T_KINEvars.kine_pio_mass < 300");
  // T_eval3->SetAlias("flag_nueCC","T_BDTvars.numu_cc_flag >=0 && T_BDTvars.nue_score > 7.0");
  // T_eval3->SetAlias("flag_NC","(!T_BDTvars.cosmict_flag) && T_BDTvars.numu_score < 0");
  // T_eval3->SetAlias("flag_pi0","(T_KINEvars.kine_pio_flag==1 && T_KINEvars.kine_pio_vtx_dis < 9 || T_KINEvars.kine_pio_flag==2) && T_KINEvars.kine_pio_energy_1 > 40 && T_KINEvars.kine_pio_energy_2 > 25 && T_KINEvars.kine_pio_dis_1 < 110 && T_KINEvars.kine_pio_dis_2 < 120 && T_KINEvars.kine_pio_angle > 0 && T_KINEvars.kine_pio_angle < 174 && T_KINEvars.kine_pio_mass > 22 && T_KINEvars.kine_pio_mass < 300");

  
  T_eval->Project("h10","truth_showerCosth:truth_showerMomentum[3]","weight_cv*weight_spline*(abs(truth_nuPdg)==12 && truth_isCC==1 && truth_vtxInside==1)");
  T_eval->Project("h11","truth_showerCosth:truth_showerMomentum[3]","weight_cv*weight_spline*(abs(truth_nuPdg)==12 && truth_isCC==1 && truth_vtxInside==1 && flag_nueCC && match_isFC==0)");
  T_eval->Project("h12","truth_showerCosth:truth_showerMomentum[3]","weight_cv*weight_spline*(abs(truth_nuPdg)==12 && truth_isCC==1 && truth_vtxInside==1 && flag_nueCC && match_isFC==0)");
  h12->Divide(h10); // efficiency 

  h10->SetTitle("Before selection");
  h10->GetXaxis()->SetTitle("E_{shower} (GeV)");
  h10->GetYaxis()->SetTitle("cos#theta_{shower}");
  h11->SetTitle("After selection");
  h11->GetXaxis()->SetTitle("E_{shower} (GeV)");
  h11->GetYaxis()->SetTitle("cos#theta_{shower}");
  h12->SetTitle("Efficiency");
  h12->GetXaxis()->SetTitle("E_{shower} (GeV)");
  h12->GetYaxis()->SetTitle("cos#theta_{shower}");

  auto c1 = new TCanvas("c1","c1",800,800);
  c1->Divide(2,2);
  c1->cd(1); h10->Draw("colz");
  c1->cd(2); h11->Draw("colz");
  c1->cd(3); h12->Draw("colz");

  c1->cd(4);
  auto gh1 = new TGraphErrors();
  auto gh2 = new TGraphErrors();
  auto gh3 = new TGraphErrors();
  auto gh4 = new TGraphErrors();

  for(int i=0; i<10; i++) {

    double a = h11->GetBinContent(2,i+1);
    double b = h10->GetBinContent(2,i+1);
    gh1->SetPoint(i,h10->GetYaxis()->GetBinCenter(i+1), a/b);
    gh1->SetPointError(i,0, delta_binomial(a,b));

    a = h11->GetBinContent(3,i+1);
    b = h10->GetBinContent(3,i+1);
    gh2->SetPoint(i,h10->GetYaxis()->GetBinCenter(i+1), a/b);
    gh2->SetPointError(i,0, delta_binomial(a,b));

    a = h11->GetBinContent(4,i+1);
    b = h10->GetBinContent(4,i+1);
    gh3->SetPoint(i,h10->GetYaxis()->GetBinCenter(i+1), a/b);
    gh3->SetPointError(i,0, delta_binomial(a,b));

    a = h11->GetBinContent(5,i+1);
    b = h10->GetBinContent(5,i+1);
    gh4->SetPoint(i,h10->GetYaxis()->GetBinCenter(i+1), a/b);
    gh4->SetPointError(i,0, delta_binomial(a,b));
  }

  gh1->SetLineColor(6);
  gh1->SetMarkerColor(6);
  gh1->SetMarkerStyle(20);
  // gh1->SetMarkerSize(2);

  gh2->SetLineColor(4);
  gh2->SetMarkerColor(4);
  gh2->SetMarkerStyle(20);
  // gh2->SetMarkerSize(2);

  gh3->SetLineColor(2);
  gh3->SetMarkerColor(2);
  gh3->SetMarkerStyle(20);
  // gh3->SetMarkerSize(2);

  gh4->SetLineColor(1);
  gh4->SetMarkerColor(1);
  gh4->SetMarkerStyle(20);
  // gh4->SetMarkerSize(2);

  gh1->GetYaxis()->SetRangeUser(0,0.5);
  gh1->GetXaxis()->SetTitle("cos#theta_{shower}");
  gh1->GetYaxis()->SetTitle("Efficiency");
  gh2->GetXaxis()->SetTitle("cos#theta_{shower}");
  gh2->GetYaxis()->SetTitle("Efficiency");
  gh3->GetXaxis()->SetTitle("cos#theta_{shower}");
  gh3->GetYaxis()->SetTitle("Efficiency");
  gh4->GetXaxis()->SetTitle("cos#theta_{shower}");
  gh4->GetYaxis()->SetTitle("Efficiency");

  gh1->Draw("APL");
  gh2->Draw("same PL");
  gh3->Draw("same PL");
  gh4->Draw("same PL");

  
  // auto c2 = new TCanvas("c2","c2",800,800);
  // gh1->Draw("APL");
  // auto c3 = new TCanvas("c3","c3",800,800);
  // gh2->Draw("APL");
  // auto c4 = new TCanvas("c4","c4",800,800);
  // gh3->Draw("APL");
  // auto c5 = new TCanvas("c5","c5",800,800);
  // gh4->Draw("APL");

 
  
}
