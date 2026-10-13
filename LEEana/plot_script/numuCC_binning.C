vector<double> get_xbins(double xstart, double xend, double resolution) {
  vector<double> xbins;
  double xcurr = xstart;
  while (xcurr < xend) {
    xbins.push_back(xcurr);
    xcurr *= (1+resolution)/(1-resolution);
  }
  xbins.push_back(xend);
  return xbins;
}

// derive the bin width from the end of the histogram
vector<double> get_xbins2(double xstart, double xend, double resolution) {
  vector<double> xbins;
  double xcurr = xend;
  while (xcurr > xstart) {
    xbins.push_back(xcurr);
    xcurr *= (1-resolution)/(1+resolution);
  }
  xbins.push_back(xstart);
  std::reverse(xbins.begin(), xbins.end());
  return xbins;
}

vector<double> rebin_xbins(TH1F* h, int nrebin=10) {
  vector<double> ret; // bin edges to return

  int nbins = h->GetNbinsX();
  if (nbins < nrebin) {
    cout << "WARNING: cannot rebin to a more bins! Keeps the current bin width." << endl;
    return ret;
  }
  double ncounts = h->Integral();
  double ncounts_perbin = ncounts / (1.0*nrebin);


  double xstart = h->GetBinLowEdge(1);
  double xend = h->GetBinLowEdge(nbins+1);
  ret.push_back(xstart);

  double accum=0; // accumulted counts in each bin
  for (int i=0; i<nbins; i++) {
    double content = h->GetBinContent(i+1);
    accum += content;
    // cout << "bin: " << i << " low: " << h->GetBinLowEdge(i+1) << " high: " << h->GetBinLowEdge(i+1) + h->GetBinWidth(i+1) << " content: " << content << endl;
    if (accum > ncounts_perbin) {
      ret.push_back(h->GetBinLowEdge(i+1) + h->GetBinWidth(i+1));
      accum = 0;
    }
  }
  if (accum > 0) {
    ret.push_back(xend);
  }
  return ret;
}

void numuCC_binning(){

  TFile *file1 =  new TFile("./processed_checkout_rootfiles/checkout_prodgenie_bnb_nu_overlay_run1.root");

  // for Enu
  TH1F *h10 = new TH1F("h10","h10",25,110,2500);
  TH1F *h11 = new TH1F("h11","h11",25,110,2500);
  TH1F *h12 = new TH1F("h12","h12",25,110,2500);
  // for Emu
  TH1F *h20 = new TH1F("h20","h20",25,105.7,2500);
  TH1F *h21 = new TH1F("h21","h21",25,105.7,2500);
  TH1F *h22 = new TH1F("h22","h22",25,105.7,2500);
  // for Etransfer
  TH1F *h30 = new TH1F("h30","h30",25,0,2500);
  TH1F *h31 = new TH1F("h31","h31",25,0,2500);
  TH1F *h32 = new TH1F("h32","h32",25,0,2500);

  TTree *T_eval = (TTree*)file1->Get("wcpselection/T_eval");
  TTree *T_BDTvars = (TTree*)file1->Get("wcpselection/T_BDTvars");
  TTree *T_KINEvars = (TTree*)file1->Get("wcpselection/T_KINEvars");
  TTree *T_PFeval = (TTree*)file1->Get("wcpselection/T_PFeval");
  T_eval->AddFriend(T_BDTvars,"T_BDTvars");
  T_eval->AddFriend(T_KINEvars,"T_KINEvars");
  T_eval->AddFriend(T_PFeval, "T_PFeval");
  T_eval->SetAlias("Ehad","truth_nuEnergy - 1000.0*T_PFeval.truth_muonMomentum[3]");


  double resolution = 0.2; // 20%
  double percentage = 1.0/3;

  // auto c0 = new TCanvas("c0","c0",600,450); 
  // c0->cd();

  // // T_eval->Draw("1000.0*T_PFeval.truth_muonMomentum[3]:Ehad >> h(100,0,2500,100,0,2500)","weight_cv*weight_spline*(abs(truth_nuPdg)==14 && truth_isCC==1 && truth_vtxInside==1 && T_BDTvars.numu_cc_flag>=0 && T_BDTvars.numu_score>0.9 && Ehad>0 && T_PFeval.truth_muonMomentum[3]>0)","colz");
  // // T_eval->Draw("truth_nuEnergy:1000.0*T_PFeval.truth_muonMomentum[3] >> h(100,0,2500,100,0,2500)","weight_cv*weight_spline*(abs(truth_nuPdg)==14 && truth_isCC==1 && truth_vtxInside==1 && T_BDTvars.numu_cc_flag>=0 && T_BDTvars.numu_score>0.9 && match_isFC==0 && Ehad>0 && T_PFeval.truth_muonMomentum[3]>0)","colz");
  // // T_eval->Draw("Ehad >> h(100,0,2500)","weight_cv*weight_spline*(abs(truth_nuPdg)==14 && truth_isCC==1 && truth_vtxInside==1 && T_BDTvars.numu_cc_flag>=0 && T_BDTvars.numu_score>0.9 && Ehad>0 && T_PFeval.truth_muonMomentum[3]>0 && match_isFC==1)");
  // // T_eval->Draw("1000.0*T_PFeval.truth_muonMomentum[3] >> h(100,0,2500)","weight_cv*weight_spline*(abs(truth_nuPdg)==14 && truth_isCC==1 && truth_vtxInside==1 && T_BDTvars.numu_cc_flag>=0 && T_BDTvars.numu_score>0.9 && Ehad>0 && T_PFeval.truth_muonMomentum[3]>0)");
  // // T_eval->Draw("truth_nuEnergy >> h(300,0,3000)","weight_cv*weight_spline*(abs(truth_nuPdg)==14 && truth_isCC==1 && truth_vtxInside==1 && T_BDTvars.numu_cc_flag>=0 && T_BDTvars.numu_score>0.9)");
  // // T_eval->Draw("T_KINEvars.kine_reco_Enu","weight_cv*weight_spline*(abs(truth_nuPdg)==14 && truth_isCC==1 && truth_vtxInside==1 && truth_nuEnergy>2250 && truth_nuEnergy<4000 && T_BDTvars.numu_cc_flag>=0 && T_BDTvars.numu_score>0.9 && match_isFC==1)");

  // resolution = 0.2; // 20%
  // percentage = 1.0/3;
  // vector<double> xbins0 = get_xbins2(200, 4000, resolution * percentage); 
  // // vector<double> xvec, yvec;
  // // for (int i=1; i<xbins0.size(); i++) {
  // //   xvec.push_back( 0.5*(xbins0.at(i) + xbins0.at(i-1)) );
  // //   yvec.push_back( 0.5*(xbins0.at(i) - xbins0.at(i-1))/(0.5*(xbins0.at(i) + xbins0.at(i-1))) );
  // // }
  // // for(int i=0; i<xvec.size(); i++) {
  // //   cout << xvec.at(i) << " " << yvec.at(i) << endl;
  // // }

  // // int nbins = 7;
  // // double xbins[] = {200, 300, 450, 675, 1000, 1500, 2250, 4000}; // 20% resolution
  // // int nbins = 10;
  // // double xbins[] = {200, 540, 705, 805, 920, 1050, 1200, 1375, 1570, 2050, 4000};
  // // TH1F* h = new TH1F("h","",nbins, xbins);
  // int nbins = xbins0.size() - 1;
  // cout << "percentage of resolution: " << percentage << " nbins: " << nbins << endl;
  // TH1F* h = new TH1F("h","",nbins, xbins0.data());
  // // TH1F* h = new TH1F("h","",50,0,5000);
  // T_eval->Draw("truth_nuEnergy >> h","weight_cv*weight_spline*(abs(truth_nuPdg)==14 && truth_isCC==1 && truth_vtxInside==1)");
  // // T_eval->Draw("truth_nuEnergy >> h(50,0,5000)","weight_cv*weight_spline*(abs(truth_nuPdg)==14 && truth_isCC==1 && truth_vtxInside==1 && match_isFC>=0 && T_BDTvars.numu_cc_flag>=0 && T_BDTvars.numu_score>0.9)","hist");
  // // auto h = (TH1F*)gROOT->FindObject("h");
  // // h->GetXaxis()->SetRangeUser(0,4000);
  // // h->GetYaxis()->SetRangeUser(0,17000);
  // h->SetLineColor(4);
  // h->SetLineStyle(2);
  // h->GetXaxis()->SetTitle("E^{#nu} [MeV]");
  // h->SetTitle("");
  // h->Draw("hist");

  // // Rebin the true spectrum with mostly equal number of events
  // vector<double> xbins_new = rebin_xbins(h, 13); // we want 10 bins, but the last bin
  //                                                // always has low counts, need to combine
  //                                                // the last two bins
  // if (xbins_new.size()>1){
  //   cout << "xbins_new nbins: " << xbins_new.size()-1 << endl;
  //   for (auto x: xbins_new) {
  //     cout << x << endl;
  //   }
  //   // // replot the true spectrum
  //   h->Delete();
  //   h = new TH1F("h","",xbins_new.size()-1, xbins_new.data());
  //   T_eval->Draw("truth_nuEnergy >> h","weight_cv*weight_spline*(abs(truth_nuPdg)==14 && truth_isCC==1 && truth_vtxInside==1)");
  //   h->SetLineColor(4);
  //   h->SetLineStyle(2);
  //   h->GetXaxis()->SetTitle("E^{#nu} [MeV]");
  //   h->SetTitle("");
  //   h->GetYaxis()->SetRangeUser(0,50000);
  //   h->Draw("hist");
  // }

  // // TH1F* h2 = new TH1F("h2","",nbins, xbins);
  // TH1F* h2 = new TH1F("h2","",25,0,2500);
  // T_eval->Draw("T_KINEvars.kine_reco_Enu >> h2","weight_cv*weight_spline*(abs(truth_nuPdg)==14 && truth_isCC==1 && truth_vtxInside==1 && T_BDTvars.numu_cc_flag>=0 && T_BDTvars.numu_score>0.9)", "hist same");
  // // T_eval->Draw("T_KINEvars.kine_reco_Enu >> h2(50,0,5000)","weight_cv*weight_spline*(abs(truth_nuPdg)==14 && truth_isCC==1 && truth_vtxInside==1 && match_isFC>=0 && T_BDTvars.numu_cc_flag>=0 && T_BDTvars.numu_score>0.9)","hist same");
  // // T_eval->Draw("truth_nuEnergy >> h2(50,0,5000)","weight_cv*weight_spline*(abs(truth_nuPdg)==14 && truth_isCC==1 && truth_vtxInside==1 && match_isFC>=0)","hist same");
  // // auto h2 = (TH1F*)gROOT->FindObject("h2");
  // // h2->Draw("hist same");

  // auto lg = new TLegend(0.6,0.6,0.8,0.8, "FC+PC");
  // lg->AddEntry(h, "E^{#nu}_{true}", "l")->SetTextColor(4);
  // lg->AddEntry(h2, "E^{#nu}_{reco}", "l")->SetTextColor(1);
  // lg->Draw();

  // // reco muon energy
  // auto c02 = new TCanvas("c02","c02",600,450); 
  // c02->cd();

  // resolution = 0.2; // 20%
  // percentage = 1.0/3;
  // vector<double> xbins_mu = get_xbins2(106, 2506, resolution * percentage);
  // int nbins_mu = xbins_mu.size() -1 ;
  // TH1F* hmu = new TH1F("hmu","",nbins_mu, xbins_mu.data());

  // // int nbins_mu = 5;
  // // double xbins_mu[] = {106, 177, 295, 491, 818, 1365}; // 25% resolution
  // // TH1F* hmu = new TH1F("hmu","",nbins_mu, xbins_mu);
  // // TH1F* hmu = new TH1F("hmu","",24,106,2506);
  // T_eval->Draw("1000.0*T_PFeval.truth_muonMomentum[3] >> hmu","weight_cv*weight_spline*(abs(truth_nuPdg)==14 && truth_isCC==1 && truth_vtxInside==1)");
  // hmu->SetLineColor(4);
  // hmu->SetLineStyle(2);
  // hmu->GetXaxis()->SetTitle("E^{#mu} [MeV]");
  // hmu->GetYaxis()->SetRangeUser(0,30000);
  // hmu->SetTitle("");
  // hmu->Draw("hist");
  // // cout << "hmu > 2506 MeV: " << hmu->GetBinContent(25) / hmu->Integral() << endl; // 1.3% for > 2506 MeV

  // // Rebin the true spectrum with mostly equal number of events
  // vector<double> xbins_mu_new = rebin_xbins(hmu, 14); // we want 10 bins, but the last bin
  //                                                // always has low counts, need to combine
  //                                                // the last two bins
  // if (xbins_mu_new.size()>1){
  //   cout << "xbins_mu_new nbins: " << xbins_mu_new.size()-1 << endl;
  //   for (auto x: xbins_mu_new) {
  //     cout << x << endl;
  //   }
  //   // // replot the true spectrum
  //   hmu->Delete();
  //   hmu = new TH1F("hmu","",xbins_mu_new.size()-1, xbins_mu_new.data());
  //   T_eval->Draw("1000.0*T_PFeval.truth_muonMomentum[3] >> hmu","weight_cv*weight_spline*(abs(truth_nuPdg)==14 && truth_isCC==1 && truth_vtxInside==1)");
  //   hmu->SetLineColor(4);
  //   hmu->SetLineStyle(2);
  //   hmu->GetXaxis()->SetTitle("E^{#mu} [MeV]");
  //   hmu->SetTitle("");
  //   hmu->GetYaxis()->SetRangeUser(0,30000);
  //   hmu->Draw("hist");
  // }


  // // TH1F* hmu2 = new TH1F("hmu2","",nbins_mu, xbins_mu);
  // TH1F* hmu2 = new TH1F("hmu2","",17,106,1806);
  // T_eval->Draw("1000.0*T_PFeval.reco_muonMomentum[3] >> hmu2","weight_cv*weight_spline*(abs(truth_nuPdg)==14 && truth_isCC==1 && truth_vtxInside==1 && T_BDTvars.numu_cc_flag>=0 && T_BDTvars.numu_score>0.9 && T_PFeval.reco_muonMomentum[3]>0)", "hist same");
  // // hmu2->Draw("hist same");

  // auto lgmu = new TLegend(0.6,0.6,0.8,0.8, "FC+PC");
  // lgmu->AddEntry(hmu, "E^{#mu}_{true}", "l")->SetTextColor(4);
  // lgmu->AddEntry(hmu2, "E^{#mu}_{reco}", "l")->SetTextColor(1);
  // lgmu->Draw();


  // hadronic energy
  auto c03 = new TCanvas("c03","c03",600,450); 
  c03->cd();

  resolution = 0.3; // 30%
  percentage = 1.0/3;
  vector<double> xbins_had = get_xbins2(30, 2500, resolution * percentage);
  int nbins_had = xbins_had.size() -1 ;
  cout << "percentage of resolution: " << percentage << " nbins: " << nbins_had << endl;
  TH1F* hhad = new TH1F("hhad","",nbins_had, xbins_had.data());

  // TH1F* hhad = new TH1F("hhad","",250,0,2500);
  
  T_eval->Draw("truth_nuEnergy - 1000.0*T_PFeval.truth_muonMomentum[3] >> hhad","weight_cv*weight_spline*(abs(truth_nuPdg)==14 && truth_isCC==1 && truth_vtxInside==1)");
  hhad->SetLineColor(4);
  hhad->SetLineStyle(2);
  hhad->GetXaxis()->SetTitle("E^{trans} [MeV]");
  hhad->GetYaxis()->SetRangeUser(0,35000);
  hhad->SetTitle("");
  hhad->Draw("hist");
  // cout << "hhad > 2500 MeV: " << hhad->GetBinContent(26) / hhad->Integral() << endl;
  // ~0% < 30MeV
  // 0.7% > 2500MeV

  // Rebin the true spectrum with mostly equal number of events
  vector<double> xbins_had_new = rebin_xbins(hhad, 17); // we want 10 bins, but the last bin
                                                 // always has low counts, need to combine
                                                 // the last two bins
  if (xbins_had_new.size()>1){
    cout << "xbins_had_new nbins: " << xbins_had_new.size()-1 << endl;
    for (auto x: xbins_had_new) {
      cout << x << endl;
    }
    // // replot the true spectrum
    hhad->Delete();
    hhad = new TH1F("hhad","",xbins_had_new.size()-1, xbins_had_new.data());
    T_eval->Draw("truth_nuEnergy - 1000.0*T_PFeval.truth_muonMomentum[3] >> hhad","weight_cv*weight_spline*(abs(truth_nuPdg)==14 && truth_isCC==1 && truth_vtxInside==1)");
    hhad->SetLineColor(4);
    hhad->SetLineStyle(2);
    hhad->GetXaxis()->SetTitle("E^{had} [MeV]");
    hhad->SetTitle("");
    hhad->GetYaxis()->SetRangeUser(0,35000);
    hhad->Draw("hist");
  }


  // TH1F* hhad2 = new TH1F("hhad2","",nbins_had, xbins_had);
  TH1F* hhad2 = new TH1F("hhad2","",20,0,2000);
  T_eval->Draw("T_KINEvars.kine_reco_Enu - 1000.0*T_PFeval.reco_muonMomentum[3] >> hhad2","weight_cv*weight_spline*(abs(truth_nuPdg)==14 && truth_isCC==1 && truth_vtxInside==1 && T_BDTvars.numu_cc_flag>=0 && T_BDTvars.numu_score>0.9 && T_PFeval.reco_muonMomentum[3]>0)", "hist same");
  // hmu2->Draw("hist same");

  auto lghad = new TLegend(0.6,0.6,0.8,0.8, "FC+PC");
  lghad->AddEntry(hhad, "E^{trans}_{true}", "l")->SetTextColor(4);
  lghad->AddEntry(hhad2, "E^{had}_{reco}", "l")->SetTextColor(1);
  lghad->Draw();



  // // efficiency vs energies
  // auto c1 = new TCanvas("c1","c1",1200,400); 
  // c1->Divide(3,1);
  
  // T_eval->Project("h10","truth_nuEnergy","weight_cv*weight_spline*(abs(truth_nuPdg)==14 && truth_isCC==1 && truth_vtxInside==1)");
  // T_eval->Project("h11","truth_nuEnergy","weight_cv*weight_spline*(abs(truth_nuPdg)==14 && truth_isCC==1 && truth_vtxInside==1 && T_BDTvars.numu_cc_flag>=0 && T_BDTvars.numu_score>0.9 && match_isFC==1)");
  // T_eval->Project("h12","truth_nuEnergy","weight_cv*weight_spline*(abs(truth_nuPdg)==14 && truth_isCC==1 && truth_vtxInside==1 && T_BDTvars.numu_cc_flag>=0 && T_BDTvars.numu_score>0.9 )");


  // T_eval->Project("h20","1000.0*T_PFeval.truth_muonMomentum[3]","weight_cv*weight_spline*(abs(truth_nuPdg)==14 && truth_isCC==1 && truth_vtxInside==1)");
  // T_eval->Project("h21","1000.0*T_PFeval.truth_muonMomentum[3]","weight_cv*weight_spline*(abs(truth_nuPdg)==14 && truth_isCC==1 && truth_vtxInside==1 && T_BDTvars.numu_cc_flag>=0 && T_BDTvars.numu_score>0.9 && match_isFC==1)");
  // T_eval->Project("h22","1000.0*T_PFeval.truth_muonMomentum[3]","weight_cv*weight_spline*(abs(truth_nuPdg)==14 && truth_isCC==1 && truth_vtxInside==1 && T_BDTvars.numu_cc_flag>=0 && T_BDTvars.numu_score>0.9 )");

  // T_eval->Project("h30","truth_nuEnergy - 1000.0*T_PFeval.truth_muonMomentum[3]","weight_cv*weight_spline*(abs(truth_nuPdg)==14 && truth_isCC==1 && truth_vtxInside==1)");
  // T_eval->Project("h31","truth_nuEnergy - 1000.0*T_PFeval.truth_muonMomentum[3]","weight_cv*weight_spline*(abs(truth_nuPdg)==14 && truth_isCC==1 && truth_vtxInside==1 && T_BDTvars.numu_cc_flag>=0 && T_BDTvars.numu_score>0.9 && match_isFC==1)");
  // T_eval->Project("h32","truth_nuEnergy - 1000.0*T_PFeval.truth_muonMomentum[3]","weight_cv*weight_spline*(abs(truth_nuPdg)==14 && truth_isCC==1 && truth_vtxInside==1 && T_BDTvars.numu_cc_flag>=0 && T_BDTvars.numu_score>0.9 )");
  

  // c1->cd(1);
  // TGraph *g1 = new TGraph();
  // TGraph *g2 = new TGraph();
  // for (Int_t i=0;i!=h10->GetNbinsX();i++){
  //   Double_t x = h10->GetBinCenter(i+1);
  //   Double_t y = h11->GetBinContent(i+1)/h10->GetBinContent(i+1);
  //   if (std::isnan(y)) y = 0;
  //   g1->SetPoint(i,x,y);
  //   y = h12->GetBinContent(i+1)/h10->GetBinContent(i+1);
  //   if (std::isnan(y)) y = 0;
  //   g2->SetPoint(i,x,y);
    
  // }
  // g1->Draw("AL*");
  // g2->Draw("*Lsame");
  // g1->SetMarkerColor(1);
  // g2->SetMarkerColor(2);
  // g1->SetMarkerStyle(20);
  // g2->SetMarkerStyle(20);
  // g1->SetLineColor(1);
  // g2->SetLineColor(2);
  // g1->GetYaxis()->SetRangeUser(0,1);
  // g1->GetXaxis()->SetRangeUser(100, 2500);

  // g1->GetXaxis()->SetTitle("E^{#nu}_{true} (MeV)");
  // // g1->GetXaxis()->SetTitle("E^{#mu}_{true} (MeV)");
  // // g1->GetXaxis()->SetTitle("E^{trans}_{true} (MeV)");
  // g1->SetTitle("Efficiency");
  // TLegend *le1 = new TLegend(0.6,0.6,0.89,0.89);
  // le1->AddEntry(g1,"FC","lp");
  // le1->AddEntry(g2,"FC+PC","lp");
  // le1->Draw();
  

  // c1->cd(2);
  // TGraph *g21 = new TGraph();
  // TGraph *g22 = new TGraph();
  // for (Int_t i=0;i!=h20->GetNbinsX();i++){
  //   Double_t x = h20->GetBinCenter(i+1);
  //   Double_t y = h21->GetBinContent(i+1)/h20->GetBinContent(i+1);
  //   if (std::isnan(y)) y = 0;
  //   g21->SetPoint(i,x,y);
  //   y = h22->GetBinContent(i+1)/h20->GetBinContent(i+1);
  //   if (std::isnan(y)) y = 0;
  //   g22->SetPoint(i,x,y);
    
  // }
  // g21->Draw("AL*");
  // g22->Draw("*Lsame");
  // g21->SetMarkerColor(1);
  // g22->SetMarkerColor(2);
  // g21->SetMarkerStyle(20);
  // g22->SetMarkerStyle(20);
  // g21->SetLineColor(1);
  // g22->SetLineColor(2);
  // g21->GetYaxis()->SetRangeUser(0,1);
  // g21->GetXaxis()->SetRangeUser(100, 2500);

  // g21->GetXaxis()->SetTitle("E^{#mu}_{true} (MeV)");
  // g21->SetTitle("Efficiency");
  // TLegend *le2 = new TLegend(0.6,0.6,0.89,0.89);
  // le2->AddEntry(g21,"FC","lp");
  // le2->AddEntry(g22,"FC+PC","lp");
  // le2->Draw();
 
  // c1->cd(3);
  // TGraph *g31 = new TGraph();
  // TGraph *g32 = new TGraph();
  // for (Int_t i=0;i!=h30->GetNbinsX();i++){
  //   Double_t x = h30->GetBinCenter(i+1);
  //   Double_t y = h31->GetBinContent(i+1)/h30->GetBinContent(i+1);
  //   if (std::isnan(y)) y = 0;
  //   g31->SetPoint(i,x,y);
  //   y = h32->GetBinContent(i+1)/h30->GetBinContent(i+1);
  //   if (std::isnan(y)) y = 0;
  //   g32->SetPoint(i,x,y);
    
  // }
  // g31->Draw("AL*");
  // g32->Draw("*Lsame");
  // g31->SetMarkerColor(1);
  // g32->SetMarkerColor(2);
  // g31->SetMarkerStyle(20);
  // g32->SetMarkerStyle(20);
  // g31->SetLineColor(1);
  // g32->SetLineColor(2);
  // g31->GetYaxis()->SetRangeUser(0,1);

  // // g31->GetXaxis()->SetTitle("E^{#nu-Ar trans}_{true} (MeV)");
  // g31->GetXaxis()->SetTitle("E^{trans}_{true} (MeV)");
  // g31->SetTitle("Efficiency");
  // TLegend *le3 = new TLegend(0.6,0.6,0.89,0.89);
  // le3->AddEntry(g31,"FC","lp");
  // le3->AddEntry(g32,"FC+PC","lp");
  // le3->Draw();
}
