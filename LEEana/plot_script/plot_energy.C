
void plot_energy(){
  TFile *file1 =  new TFile("processed_checkout_rootfiles/checkout_prodgenie_bnb_nu_overlay_run123_all.root");

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
  T_eval->SetAlias("nu","T_PFeval.truth_muonMomentum[3]>0 ? T_eval.truth_nuEnergy*1e-3 - T_PFeval.truth_muonMomentum[3] : -1");

  auto c1 = new TCanvas("c1","c1",1300,500);
  c1->Divide(3,1);

  c1->cd(1);
  T_eval->Draw("Etrue*1e-3 >> h1(100,0,4)","weight_cv*weight_spline*(abs(truth_nuPdg)==14 && truth_isCC==1 && truth_vtxInside==1 && T_BDTvars.numu_cc_flag>=0 && T_BDTvars.numu_score>0.9 && match_isFC>0)");
  T_eval->Draw("Etrue*1e-3 >> h1PC(100,0,4)","weight_cv*weight_spline*(abs(truth_nuPdg)==14 && truth_isCC==1 && truth_vtxInside==1 && T_BDTvars.numu_cc_flag>=0 && T_BDTvars.numu_score>0.9 && match_isFC==0)");
  auto h1 = (TH2F*)gROOT->FindObject("h1");
  h1->GetXaxis()->SetTitle("E_{#nu}^{true} (GeV)");
  h1->SetTitle("");
  h1->SetLineColor(4);
  auto h1PC = (TH2F*)gROOT->FindObject("h1PC");
  h1PC->GetXaxis()->SetTitle("E_{#nu}^{true} (GeV)");
  h1PC->SetTitle("");
  h1PC->SetLineColor(2);
  h1PC->Draw("hist");
  h1->Draw("hist same");

  auto lg = new TLegend(0.6,0.6,0.8,0.8);
  lg->AddEntry(h1,"FC","le")->SetTextColor(4);
  lg->AddEntry(h1PC,"PC","le")->SetTextColor(2);
  lg->SetBorderSize(0);
  lg->Draw();

  vector<float> xbins1 = {0.2, 0.540, 0.705, 0.805, 0.920, 1.050, 1.200, 1.375, 1.570, 2.050, 4.000}; // GeV
  int i=1;
  for(auto x: xbins1){
    auto l = new TLine(x,0,x,8000);
    l->Draw("same");
    if(x<4.000) {TLatex t; t.DrawLatex(x,8000,to_string(i).c_str()); i++;}
  }

  c1->cd(2);
  T_eval->Draw("Emu >> h2(100,0,2.5)","weight_cv*weight_spline*(abs(truth_nuPdg)==14 && truth_isCC==1 && truth_vtxInside==1 && T_BDTvars.numu_cc_flag>=0 && T_BDTvars.numu_score>0.9 && match_isFC>0)");
  T_eval->Draw("Emu >> h2PC(100,0,2.5)","weight_cv*weight_spline*(abs(truth_nuPdg)==14 && truth_isCC==1 && truth_vtxInside==1 && T_BDTvars.numu_cc_flag>=0 && T_BDTvars.numu_score>0.9 && match_isFC==0)");
  auto h2 = (TH2F*)gROOT->FindObject("h2");
  h2->GetXaxis()->SetTitle("E_{#mu}^{true} (GeV)");
  h2->SetTitle("");
  h2->SetLineColor(4);
  auto h2PC = (TH2F*)gROOT->FindObject("h2PC");
  h2PC->GetXaxis()->SetTitle("E_{#mu}^{true} (GeV)");
  h2PC->SetTitle("");
  h2PC->SetLineColor(2);
  h2PC->Draw("hist");
  h2->Draw("hist same");

  xbins1.clear();
  xbins1 = {0.106, 0.226, 0.296, 0.386, 0.505, 0.577, 0.659, 0.753, 0.861, 0.984, 1.285, 2.506}; // GeV
  i=1;
  for(auto x: xbins1){
    auto l = new TLine(x,0,x,6000);
    l->Draw("same");
    if (x<2.506) {TLatex t; t.DrawLatex(x,6000,to_string(i).c_str()); i++;}
  }


  c1->cd(3);
  T_eval->Draw("nu >> h3(100,0,2.5)","weight_cv*weight_spline*(T_PFeval.truth_muonMomentum[3]>0 && abs(truth_nuPdg)==14 && truth_isCC==1 && truth_vtxInside==1 && T_BDTvars.numu_cc_flag>=0 && T_BDTvars.numu_score>0.9 && match_isFC>0)");
  T_eval->Draw("nu >> h3PC(100,0,2.5)","weight_cv*weight_spline*(T_PFeval.truth_muonMomentum[3]>0 && abs(truth_nuPdg)==14 && truth_isCC==1 && truth_vtxInside==1 && T_BDTvars.numu_cc_flag>=0 && T_BDTvars.numu_score>0.9 && match_isFC==0)");
  auto h3 = (TH2F*)gROOT->FindObject("h3");
  h3->GetXaxis()->SetTitle("#nu (E_{#nu}^{true} - E_{#mu}^{true}, GeV)");
  h3->SetTitle("");
  h3->SetLineColor(4);
  auto h3PC = (TH2F*)gROOT->FindObject("h3PC");
  h3PC->GetXaxis()->SetTitle("#nu (E_{#nu}^{true} - E_{#mu}^{true}, GeV)");
  h3PC->SetTitle("");
  h3PC->SetLineColor(2);
  h3PC->Draw("hist");
  h3->Draw("hist same");

  xbins1.clear();
  xbins1 = {0.03, 0.15, 0.275, 0.411, 0.502, 0.614, 0.75, 1.12, 2.5}; // GeV
  i=1;
  for(auto x: xbins1){
    auto l = new TLine(x,0,x,10000);
    l->Draw("same");
    if (x<2.5) {TLatex t; t.DrawLatex(x,10000,to_string(i).c_str()); i++;}
  }


}
