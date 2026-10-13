TGraphAsymmErrors* get_Enu_rez(TH2F* h2, int nfirst=3){ // nfirst: fitst bin for the calculation
  TGraphAsymmErrors* gEnuRes = new TGraphAsymmErrors();
  for(int i=1; i<=h2->GetNbinsX(); i++)
  {
      if(i<nfirst) continue;
      TH1D* htemp = h2->ProjectionY("htemp", i, i, "");
      double mean = htemp->GetBinCenter(htemp->GetMaximumBin());
      double xq0 = htemp->Integral(1,htemp->GetMaximumBin())/htemp->Integral();
      double par[3];
      par[1]=mean;
      double xq = xq0;
      // std::cout<<"Mean quantile: "<<xq<<std::endl;
      xq = xq0*(1-0.6827);
      htemp->GetQuantiles(1,&par[0],&xq);
      // std::cout<<"Lower quantile: "<<xq<<std::endl;
      // std::cout<<"Lower rms: "<<par[1]-par[0]<<std::endl;
      // std::cout << h2->GetXaxis()->GetBinCenter(i) << " " << par[1]-par[0] << " " ;
      xq = xq0 + (1-xq0)*0.6827;
      htemp->GetQuantiles(1,&par[2],&xq);
      // std::cout<<"higher quantile: "<<xq<<std::endl;
      // std::cout<<"higher rms: "<<par[2]-par[1]<<std::endl;
      // std::cout << par[2]-par[1] << std::endl;
      int N=gEnuRes->GetN();
      gEnuRes->SetPoint(N, h2->GetXaxis()->GetBinCenter(i), mean);
      gEnuRes->SetPointError(N, 0.5*h2->GetXaxis()->GetBinWidth(1),  0.5*h2->GetXaxis()->GetBinWidth(1), par[1]-par[0], par[2]-par[1] );
  }
  gEnuRes->SetMarkerStyle(20);
  return gEnuRes;

}

void plot_pat_recog_nu_rez(){
  TFile *file1 =  new TFile("processed_checkout_rootfiles/checkout_prodgenie_bnb_nu_overlay_run123_all.root");
  // TFile *file1 =  new TFile("processed_checkout_rootfiles/checkout_prodgenie_bnb_intrinsic_nue_overlay_run123_all.root");

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

  auto c1 = new TCanvas("c1","c1",1800,600);
  c1->Divide(3,1);

  // fully-contained numu (or nue) CC (inside active vol) tagged with BDT
  c1->cd(1);
  T_eval->Draw("Evis*1e-3:Etrue*1e-3 >> h1(60,0,3,60,0,3)","weight_cv*weight_spline*(abs(truth_nuPdg)==14 && truth_isCC==1 && truth_vtxInside==1 && T_BDTvars.numu_cc_flag>=0 && T_BDTvars.numu_score>0.9 && match_isFC>0)","colz");
  // T_eval->Draw("Evis*1e-3:Etrue*1e-3 >> h1(60,0,3,60,0,3)","weight_cv*weight_spline*(abs(truth_nuPdg)==12 && truth_isCC==1 && truth_vtxInside==1 && T_BDTvars.numu_cc_flag>=0 && T_BDTvars.nue_score>7.0 && match_isFC>0)","colz");
  auto h1 = (TH2F*)gROOT->FindObject("h1");
  h1->GetXaxis()->SetTitle("E_{#nu}^{true} (GeV)");
  h1->GetYaxis()->SetTitle("E_{#nu}^{vis} (GeV)");
  h1->SetTitle("");
  gPad->SetLogz();
  auto l1 = new TLine(0.1,0.1,3,3);
  l1->SetLineColor(1);
  l1->SetLineStyle(2);
  l1->SetLineWidth(2);
  l1->Draw();


  c1->cd(2);
  T_eval->Draw("Ereco*1e-3:Etrue*1e-3 >> h2(60,0,3,60,0,3)","weight_cv*weight_spline*(abs(truth_nuPdg)==14 && truth_isCC==1 && truth_vtxInside==1 && T_BDTvars.numu_cc_flag>=0 && T_BDTvars.numu_score>0.9 && match_isFC>0)","colz");
  // T_eval->Draw("Ereco*1e-3:Etrue*1e-3 >> h2(60,0,3,60,0,3)","weight_cv*weight_spline*(abs(truth_nuPdg)==12 && truth_isCC==1 && truth_vtxInside==1 && T_BDTvars.numu_cc_flag>=0 && T_BDTvars.nue_score>7.0 && match_isFC>0)","colz");
  auto h2 = (TH2F*)gROOT->FindObject("h2");
  h2->GetXaxis()->SetTitle("E_{#nu}^{true} (GeV)");
  h2->GetYaxis()->SetTitle("E_{#nu}^{rec} (GeV)");
  h2->SetTitle("");
  gPad->SetLogz();
  auto l2 = new TLine(0.1,0.1,3,3);
  l2->SetLineColor(1);
  l2->SetLineStyle(2);
  l2->SetLineWidth(2);
  l2->Draw();

  c1->cd(3);
  T_eval->Draw("Ereco/Etrue:Etrue*1e-3 >> h3(30,0,3,40,0,2)","weight_cv*weight_spline*(abs(truth_nuPdg)==14 && truth_isCC==1 && truth_vtxInside==1 && T_BDTvars.numu_cc_flag>=0 && T_BDTvars.numu_score>0.9 && match_isFC>0)","colz");
  // T_eval->Draw("Ereco/Etrue:Etrue*1e-3 >> h3(30,0,3,40,0,2)","weight_cv*weight_spline*(abs(truth_nuPdg)==12 && truth_isCC==1 && truth_vtxInside==1 && T_BDTvars.numu_cc_flag>=0 && T_BDTvars.nue_score>7.0 && match_isFC>0)","colz");
  auto h3 = (TH2F*)gROOT->FindObject("h3");
  h3->GetXaxis()->SetTitle("E_{#nu}^{true} (GeV)");
  h3->GetYaxis()->SetTitle("E_{#nu}^{rec} / E_{#nu}^{true}");
  h3->SetTitle("");
  gPad->SetLogz();
  auto gEnuRes = get_Enu_rez(h3);
  gEnuRes->Draw("P same");

}
