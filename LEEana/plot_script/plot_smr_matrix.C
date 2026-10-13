void plot_smr_matrix(){

  gStyle->SetPalette(kRainBow);

  TFile *file1 =  new TFile("../processed_checkout_rootfiles/checkout_prodgenie_bnb_nu_overlay_run123_all.root");

  TTree *T_eval = (TTree*)file1->Get("wcpselection/T_eval");
  TTree *T_BDTvars = (TTree*)file1->Get("wcpselection/T_BDTvars");
  TTree *T_KINEvars = (TTree*)file1->Get("wcpselection/T_KINEvars");
  TTree *T_PFeval = (TTree*)file1->Get("wcpselection/T_PFeval");
  T_eval->AddFriend(T_BDTvars,"T_BDTvars");
  T_eval->AddFriend(T_KINEvars,"T_KINEvars");
  T_eval->AddFriend(T_PFeval, "T_PFeval");
  T_eval->SetAlias("Etrue","T_eval.truth_nuEnergy");
  T_eval->SetAlias("Evis","T_eval.match_energy");
  T_eval->SetAlias("Ereco","T_KINEvars.kine_reco_Enu");

  T_eval->SetAlias("Emu","T_PFeval.truth_muonMomentum[3]");
  T_eval->SetAlias("Emureco","T_PFeval.reco_muonMomentum[3]");

  T_eval->SetAlias("nu","T_eval.truth_nuEnergy*1e-3 - T_PFeval.truth_muonMomentum[3]");
  T_eval->SetAlias("Ehad","T_PFeval.reco_muonMomentum[3]>0 ? T_KINEvars.kine_reco_Enu*1e-3 - T_PFeval.reco_muonMomentum[3] : -1");

  // canvas 1
  auto c1 = new TCanvas("c1","c1",1100,500);
  c1->Divide(2,1);
  int nbinsx = 10;
  double xbins[] = {0.2, 0.540, 0.705, 0.805, 0.920, 1.050, 1.200, 1.375, 1.570, 2.050, 4.000}; // GeV
  int nbinsy = 23;
  double ybins[24];
  for(int i=0; i<24; i++) {
    ybins[i] = 0.2 + i*0.1; 
  } 

  auto h1FC = new TH2F("h1FC","",nbinsx,xbins, nbinsy,ybins);
  auto h1PC = new TH2F("h1PC","",nbinsx,xbins, nbinsy,ybins);
  h1FC->GetXaxis()->SetTitle("True #nu Energy (GeV)");
  h1FC->GetYaxis()->SetTitle("Reco #nu Energy (GeV)");
  h1PC->GetXaxis()->SetTitle("True #nu Energy (GeV)");
  h1PC->GetYaxis()->SetTitle("Reco #nu Energy (GeV)");

  c1->cd(1); gPad->SetLogz();
  T_eval->Project("h1FC", "Ereco*1e-3:Etrue*1e-3","weight_cv*weight_spline*(abs(truth_nuPdg)==14 && truth_isCC==1 && truth_vtxInside==1 && T_BDTvars.numu_cc_flag>=0 && T_BDTvars.numu_score>0.9 && match_isFC>0)");
  h1FC->Draw("colz");
  auto l1 = new TLine(0.2,0.2,2.5,2.5);
  l1->SetLineStyle(9);
  l1->SetLineWidth(4);
  l1->Draw("same");

  c1->cd(2); gPad->SetLogz();
  T_eval->Project("h1PC", "Ereco*1e-3:Etrue*1e-3","weight_cv*weight_spline*(abs(truth_nuPdg)==14 && truth_isCC==1 && truth_vtxInside==1 && T_BDTvars.numu_cc_flag>=0 && T_BDTvars.numu_score>0.9 && match_isFC==0)");
  h1PC->Draw("colz");
  l1->Draw("same");

  // canvas 2
  auto c2 = new TCanvas("c2","c2",1100,500);
  c2->Divide(2,1);
  int nbinsx2 = 11;
  double xbins2[] = {0.106, 0.226, 0.296, 0.386, 0.505, 0.577, 0.659, 0.753, 0.861, 0.984, 1.285, 2.506};
  int nbinsy2 = 19;
  double ybins2[20];
  for(int i=0; i<20; i++) {
    ybins2[i] = 0.1+i*0.1; 
  } 

  auto h2FC = new TH2F("h2FC","",nbinsx2,xbins2, nbinsy2,ybins2);
  auto h2PC = new TH2F("h2PC","",nbinsx2,xbins2, nbinsy2,ybins2);
  h2FC->GetXaxis()->SetTitle("True #mu Energy (GeV)");
  h2FC->GetYaxis()->SetTitle("Reco #mu Energy (GeV)");
  h2PC->GetXaxis()->SetTitle("True #mu Energy (GeV)");
  h2PC->GetYaxis()->SetTitle("Reco #mu Energy (GeV)");

  c2->cd(1); gPad->SetLogz();
  T_eval->Project("h2FC", "Emureco:Emu","weight_cv*weight_spline*(abs(truth_nuPdg)==14 && truth_isCC==1 && truth_vtxInside==1 && T_BDTvars.numu_cc_flag>=0 && T_BDTvars.numu_score>0.9 && match_isFC>0)");
  h2FC->Draw("colz");
  auto l2 = new TLine(0.1,0.1,2.0,2.0);
  l2->SetLineStyle(9);
  l2->SetLineWidth(4);
  l2->Draw("same");

  c2->cd(2); gPad->SetLogz();
  T_eval->Project("h2PC", "Emureco:Emu","weight_cv*weight_spline*(abs(truth_nuPdg)==14 && truth_isCC==1 && truth_vtxInside==1 && T_BDTvars.numu_cc_flag>=0 && T_BDTvars.numu_score>0.9 && match_isFC==0)");
  h2PC->Draw("colz");
  l2->Draw("same");

  // canvas 3
  auto c3 = new TCanvas("c3","c3",1100,500);
  c3->Divide(2,1);
  int nbinsx3 = 8;
  double xbins3[] = {0.03, 0.15, 0.275, 0.411, 0.502, 0.614, 0.75, 1.12, 2.5}; // GeV
  int nbinsy3 = 25;
  double ybins3[26];
  for(int i=0; i<26; i++) {
    ybins3[i] = i*0.1; 
  } 

  auto h3FC = new TH2F("h3FC","",nbinsx3,xbins3, nbinsy3,ybins3);
  auto h3PC = new TH2F("h3PC","",nbinsx3,xbins3, nbinsy3,ybins3);
  h3FC->GetXaxis()->SetTitle("True #nu = E_{#nu} - E_{#mu} (GeV)");
  h3FC->GetYaxis()->SetTitle("Reco E_{had} (GeV)");
  h3PC->GetXaxis()->SetTitle("True #nu = E_{#nu} - E_{#mu} (GeV)");
  h3PC->GetYaxis()->SetTitle("Reco E_{had} (GeV)");

  c3->cd(1); gPad->SetLogz();
  T_eval->Project("h3FC", "Ehad:nu","weight_cv*weight_spline*(abs(truth_nuPdg)==14 && truth_isCC==1 && truth_vtxInside==1 && T_BDTvars.numu_cc_flag>=0 && T_BDTvars.numu_score>0.9 && match_isFC>0)");
  h3FC->Draw("colz");
  auto l3 = new TLine(0.0,0.0,2.5,2.5);
  l3->SetLineStyle(9);
  l3->SetLineWidth(4);
  l3->Draw("same");

  c3->cd(2); gPad->SetLogz();
  T_eval->Project("h3PC", "Ehad:nu","weight_cv*weight_spline*(abs(truth_nuPdg)==14 && truth_isCC==1 && truth_vtxInside==1 && T_BDTvars.numu_cc_flag>=0 && T_BDTvars.numu_score>0.9 && match_isFC==0)");
  h3PC->Draw("colz");
  l3->Draw("same");

}
