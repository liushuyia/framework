// mean value in a bin
double get_mean_x(TGraph* gh, double low, double high, double threshold=0.11){ // 110MeV threshold for CCQE
  int N = 20;
  double sum_yx=0, sum_y=0;
  for (int i=0; i<=N; i++) {
    double x = low + (high-low)/(1.0*N) * i;
    double y = gh->Eval(x);
    if(x>threshold) {
      sum_yx += y*x;
      sum_y += y;
    }
  }
  // cout << "low: " << low << " high: " << high << " mean: " << sum_yx/sum_y << endl;
  return sum_yx/sum_y; // mean in x
}

// flux weighted xsec
double get_mean_sigma(TGraph* gh, TGraph* gx, double low, double high, double threshold=0.11){
  int N = 20;
  double sum_fxsec=0, sum_f=0;
  for (int i=0; i<=N; i++) {
    double x = low + (high-low)/(1.0*N) * i;
    double F = gh->Eval(x);
    double xsec = gx->Eval(x);
    if (x>threshold){
      sum_fxsec += F*xsec;
      sum_f += F;
    }
  }
  return sum_fxsec/sum_f;
}

// GoF
double calc_GoF(TMatrixD matrix_pred, TMatrixD matrix_data, TMatrixD cov){
  TMatrixD md = matrix_pred - matrix_data;
  TMatrixD mdT = matrix_pred - matrix_data;
  mdT.T();
  cov.Invert();
  auto mret = md * cov * mdT;
  cout << mret.GetNrows() << " x " << mret.GetNcols() << " chi2/NDF= " << mret(0,0) << "/" << md.GetNcols() << endl;
  return mret(0,0);
}

// main function
void plot_xs_tot(int opt=2){

  std::vector<double> rt_data2ubmc, err_rt_data2ubmc;  // ratio data-to-uBTune

  auto c1 = new TCanvas("c1","c1",800,600);

  // unfolded result
  auto uBxsec = TFile::Open("output.root");
  auto unfold = (TH1D*)uBxsec->Get("unf");
  auto absError = (TH1D*)uBxsec->Get("absError");

  // real binning
  // int nbins = 11;
  // double xbins1[] = {0,0.5,0.6,0.7,0.8,0.9,1,1.1,1.2,1.5,2.1,5.0}; // GeV
  // int nbins = 10;
  // double xbins1[] = {0.2, 0.540, 0.705, 0.805, 0.920, 1.050, 1.200, 1.375, 1.570, 2.050, 4.000}; // GeV
  int nbins = 7;
  double xbins1[] = {0.2, 0.54, 0.705, 0.920, 1.2, 1.57, 2.05, 4.0}; // GeV

  TH1D* unfold1 = new TH1D("unfold_real","",nbins, xbins1);
  TVectorD unfold1_vd(nbins);
  unfold1->GetXaxis()->SetTitle("E_{#nu} [GeV]");
  unfold1->GetYaxis()->SetTitle("#sigma_{CC} (#times 10^{-36} cm^{2} / Ar)");
  for(int i=0; i<nbins; i++) {
    unfold1->SetBinContent(i+1, unfold->GetBinContent(i+1));
    unfold1->SetBinError(i+1, absError->GetBinContent(i+1));
    unfold1_vd(i) = unfold->GetBinContent(i+1);
  }
  unfold1->SetMaximum(1);
  // unfold1->Draw("E");

  // set assymetric error bar in x-axis
  // auto flux = TFile::Open("../flux_info/gh_averaged_numu_flux.root");
  // auto gh_flux = (TGraph*)flux->Get("gh_averaged_numu_flux");
  auto flux = TFile::Open("../flux_info/numi_flux_graphs.root");
  auto gh_flux = (TGraph*)flux->Get("nue");
  // data points
  std::vector<double> x_v, y_v;
  std::vector<double> exl_v, exh_v;
  std::vector<double> eyl_v, eyh_v;
  for (int i=0; i<nbins; i++){
    double low = unfold1->GetBinLowEdge(i+1);
    double width = unfold1->GetBinWidth(i+1);
    double y = unfold1->GetBinContent(i+1);
    double ey = unfold1->GetBinError(i+1);
    double x = get_mean_x(gh_flux, low, low+width);

    x_v.push_back(x); 
    y_v.push_back(y); 
    exl_v.push_back(x-low); 
    exh_v.push_back(low+width-x); 
    eyl_v.push_back(ey); 
    eyh_v.push_back(ey); 
  }


  // 0. xsec model
  TFile *file1 = new TFile("../genie/xsec_graphs_v2_tune1.root");
  TFile *file2 = new TFile("../genie/xsec_graphs_v3_10a_02_11a.root");
  // double mc_weight[10] = {1.21993, 1.19518, 1.19197, 1.17844, 1.16532, 1.15571, 1.15111, 1.13284, 1.11594, 1.0744};
  // double mc_weight[10] = {1.5178, 1.26847, 1.21674, 1.18273, 1.15804, 1.14223, 1.13518, 1.11799, 1.10339, 1.0671}; // no spline
 
  TGraph *g1 = (TGraph*)file1->Get("nu_e_Ar40/tot_cc");
  TGraph *g2 = (TGraph*)file2->Get("nu_e_Ar40/tot_cc");

  // scale by 1/100
  int N1 = g1->GetN();
  double* X1 = g1->GetX();
  double* Y1 = g1->GetY();
  for (int i=0; i<N1; i++) {
    Y1[i] *= 0.01;
  }
  int N2 = g2->GetN();
  double* X2 = g2->GetX();
  double* Y2 = g2->GetY();
  for (int i=0; i<N2; i++) {
    Y2[i] *= 0.01;
  }
  g1 = new TGraph(N1,X1,Y1);
  g2 = new TGraph(N2,X2,Y2);

  // if (opt==0){
  //   cout << "original GENIE" << endl;
  //   g1->Draw("Lsame");
  //   g1->GetXaxis()->SetRangeUser(0,3);
  //   g1->GetYaxis()->SetRangeUser(0,120);
  //   g2->Draw("Lsame");
  //   g1->SetLineColor(1);
  //   g2->SetLineColor(1);
  //   g1->SetLineWidth(2);
  //   g2->SetLineWidth(2);
  //   g1->SetLineStyle(2);
  // }

  // 1. xsec model + real binning
  // flux weighted
  std::vector<double> erg, xsec1, xsec2;
  for (int i=0; i<nbins; i++){
    double low = unfold1->GetBinLowEdge(i+1);
    double width = unfold1->GetBinWidth(i+1);
    double flux_weighted_xsec1 = get_mean_sigma(gh_flux, g1, low, low+width);
    double flux_weighted_xsec2 = get_mean_sigma(gh_flux, g2, low, low+width);
    double mean_e = get_mean_x(gh_flux, low, low+width);
    erg.push_back(mean_e);
    xsec1.push_back(flux_weighted_xsec1);
    xsec2.push_back(flux_weighted_xsec2);
    // xsec2.push_back(flux_weighted_xsec2 * mc_weight[i] * 1.015);
  }

  auto g1_flux_weighted = new TGraph(erg.size(), erg.data(), xsec1.data());
  auto g2_flux_weighted = new TGraph(erg.size(), erg.data(), xsec2.data());

 
  // if (opt ==1){
  //   cout << "weighted " << endl;
  //   g1_flux_weighted->Draw("PLsame");
  //   g1_flux_weighted->GetXaxis()->SetRangeUser(0,3);
  //   g1_flux_weighted->GetYaxis()->SetRangeUser(0,120);
  //   g2_flux_weighted->Draw("PLsame");
  //   g1_flux_weighted->SetLineColor(1);
  //   g2_flux_weighted->SetLineColor(1);
  //   g1_flux_weighted->SetLineWidth(2);
  //   g2_flux_weighted->SetLineWidth(2);
  //   g1_flux_weighted->SetLineStyle(2);
  //   g1_flux_weighted->SetMarkerColor(1);
  //   g2_flux_weighted->SetMarkerColor(1);
  //   g1_flux_weighted->SetMarkerSize(1);
  //   g2_flux_weighted->SetMarkerSize(1);
  // }

  // plot MC truth from CV
  auto fmerge = TFile::Open("merge_xs.root");
  TVectorD * vec_signal = (TVectorD*)fmerge->Get("vec_signal");
  int nbin_true = vec_signal->GetNoElements();
  TH1D *htrue_signal = new TH1D("htrue_signal","htrue_signal",nbins, xbins1);
  for(int i=0; i<nbins; i++) {
    htrue_signal->SetBinContent(i+1, (*vec_signal)(i));
  }
  // htrue_signal->SetLineStyle(2);
  // htrue_signal->Draw("same");
 

  // 2. xsec model + real binning + Ac smearing  
  double* x1 = g1_flux_weighted->GetX();
  double* y1 = g1_flux_weighted->GetY();
  double* x2 = g2_flux_weighted->GetX();
  double* y2 = g2_flux_weighted->GetY();
  vector<double> x3, y3;
  auto smear = (TH2D*)uBxsec->Get("smear"); // Ac from unfolding
  int nrows = smear->GetNbinsX();
  int ncols = smear->GetNbinsY();
  TMatrixD m(nrows, ncols);
  TVectorD v1(nbins), v2(nbins), v3(nbins);
  // if (opt==2){

    for (int i=0; i<nrows; i++) {
      v1(i) = y1[i];
      v2(i) = y2[i];
      v3(i) = htrue_signal->GetBinContent(i+1);
      for (int j=0; j<ncols; j++) {
        m(i,j) = smear->GetBinContent(i+1, j+1);
        // if (i==j) m(i,j) = 1;
        // else m(i,j) = 0;
      }
    }

    // with Ac smearing
    TVectorD Ac_v1 = m * v1;
    TVectorD Ac_v2 = m * v2;
    TVectorD Ac_v3 = m * v3;
    // TVectorD Ac_v1 =  v1;
    // TVectorD Ac_v2 =  v2;
    // TVectorD Ac_v3 =  v3;
    for (int i=0; i<nbins; i++) {
      y1[i] = Ac_v1(i);
      y2[i] = Ac_v2(i);
      x3.push_back(x1[i]);
      y3.push_back(Ac_v3(i));
    }
    
    // without Ac smearing
    // for (int i=0; i<nbins; i++) {
    //   x3.push_back(x1[i]);
    //   y3.push_back(v3(i));
    //   cout << x1[i] << " " << v3(i) << endl;
    // }
    // TVector Ac_add = m*v3 - v3; // additional uncertianty to add to data points
    // TVector Ac_add = m*unfold1_vd - unfold1_vd; // additional uncertianty to add to data points
    TVector Ac_add1 = m*v1 - v1; // additional uncertianty to add to data points
    TVector Ac_add2 = m*v2 - v2; // additional uncertianty to add to data points
    TVector Ac_add3 = m*v3 - v3; // additional uncertianty to add to data points
    TVector Ac_add4 = m*unfold1_vd - unfold1_vd; // additional uncertianty to add to data points
    // auto gr = new TGraphAsymmErrors(x_v.size(),x_v.data(),y_v.data(),exl_v.data(),exh_v.data(),eyl_v.data(),eyh_v.data());

    auto g1_flux_weighted_smeared = new TGraph(nbins, x1, y1); // GENIE  v2
    auto g2_flux_weighted_smeared = new TGraph(nbins, x2, y2); // v3
    auto g3_flux_weighted_smeared = new TGraph(nbins, x3.data(), y3.data()); // MC truth
   
    for (int i=0; i<nbins; i++) {
      rt_data2ubmc.push_back( y_v.at(i)/ y3.at(i));
      err_rt_data2ubmc.push_back( eyl_v.at(i)/ y3.at(i));
    }
 
    cout << "weighted smeared" << endl;
    g1_flux_weighted_smeared->SetTitle("");
    g1_flux_weighted_smeared->GetXaxis()->SetTitle("E_{#nu} [GeV]");
    // g1_flux_weighted_smeared->GetYaxis()->SetTitle("10^{-36} cm^{2} / Ar");
    g1_flux_weighted_smeared->GetYaxis()->SetTitle("#sigma_{CC} (#times 10^{-36} cm^{2} / Ar)");
    g1_flux_weighted_smeared->GetXaxis()->SetRangeUser(0,4);
    g1_flux_weighted_smeared->GetXaxis()->SetLimits(0,4);
    g1_flux_weighted_smeared->GetYaxis()->SetRangeUser(0,1);
    g1_flux_weighted_smeared->Draw("APL");

    g2_flux_weighted_smeared->Draw("PLsame");
    g1_flux_weighted_smeared->SetLineColor(1);
    g2_flux_weighted_smeared->SetLineColor(1);
    g1_flux_weighted_smeared->SetLineWidth(2);
    g2_flux_weighted_smeared->SetLineWidth(2);
    g1_flux_weighted_smeared->SetLineStyle(2);
    g1_flux_weighted_smeared->SetMarkerColor(1);
    g2_flux_weighted_smeared->SetMarkerColor(1);
    g1_flux_weighted_smeared->SetMarkerSize(1);
    g2_flux_weighted_smeared->SetMarkerSize(1);
    g3_flux_weighted_smeared->SetLineWidth(2);
    g3_flux_weighted_smeared->SetLineStyle(1);
    g3_flux_weighted_smeared->SetLineColor(8);
    g3_flux_weighted_smeared->Draw("Lsame");

  // }


  // GOF (model+real bin+Ac vs. data)
  TMatrixD matrix_pred1(1,nbins);
  TMatrixD matrix_pred2(1,nbins);
  TMatrixD matrix_pred3(1,nbins);
  TMatrixD matrix_data(1,nbins);
  cout << "nbins: " << nbins << endl;
  for(int i=0; i<nbins; i++){
    matrix_pred1(0,i) = y1[i]; // GENIE v2
    matrix_pred2(0,i) = y2[i]; // GENIE v3

    matrix_pred3(0,i) = y3[i]; // MC truth

    matrix_data(0,i) = y_v.at(i);
  }

  auto cov = (TH2D*)uBxsec->Get("unfcov");
  auto bias = (TH1D*)uBxsec->Get("bias");
  TMatrixD unfcov(nbins, nbins); // fill cov matrix
  for (int i=0; i<nbins; i++) {
    for (int j=0; j<nbins; j++) {
      unfcov(i,j) = cov->GetBinContent(i+1, j+1);
      // in the case that Ac smearing is not added,
      // use Ac-I diagnal term as additional uncer.
      if (i==j) {
        double add = ( Ac_add1(i) * Ac_add1(i) + Ac_add2(i) * Ac_add2(i) +Ac_add3(i) * Ac_add3(i) + Ac_add4(i) * Ac_add4(i) ) * 0.25;
        // unfcov(i,i) += add;

        // cout << "unfold1 # " << i+1 << " error: " << unfold1->GetBinError(i+1) << " add: " << add << endl;
        // unfold1->SetBinError(i+1, sqrt( unfcov(i,i) ) );
        // unfold1->SetBinError(i+1, 0);
      //   cout << "cov diagonal # " << i << " " <<  sqrt(cov->GetBinContent(i+1, i+1)) << " bias: " << Ac_add(i) << endl;
      }
    }
  }


  for(int i=0; i<nbins; i++) {
    eyl_v[i] = sqrt(unfcov(i,i));
    eyh_v[i] = sqrt(unfcov(i,i));
  }
  auto gr = new TGraphAsymmErrors(x_v.size(),x_v.data(),y_v.data(),exl_v.data(),exh_v.data(),eyl_v.data(),eyh_v.data());
  gr->GetXaxis()->SetTitle("E_{#nu} [GeV]");
  gr->GetYaxis()->SetTitle("#sigma_{CC} (#times 10^{-36} cm^{2} / Ar)");
  gr->SetMarkerColor(4);
  gr->SetMarkerStyle(21);
  gr->SetMarkerSize(1);
  gr->SetLineWidth(2);
  gr->SetLineColor(4);
  gr->SetMaximum(1);
  gr->SetTitle("");
  gr->Draw("Psame");



  double chi2_m1 = calc_GoF(matrix_pred1, matrix_data, unfcov);
  double chi2_m2 = calc_GoF(matrix_pred2, matrix_data, unfcov);
  double chi2_m3 = calc_GoF(matrix_pred3, matrix_data, unfcov);

  auto lg = new TLegend(0.11,0.7,0.5,0.89);
  lg->SetBorderSize(0);
  lg->AddEntry(gr, "Data","lpe")->SetTextColor(4);
  lg->AddEntry(g1_flux_weighted_smeared, Form("GENIE v2, #chi^{2}/dof=%.1f/%d",chi2_m1,nbins),"lp");
  lg->AddEntry(g2_flux_weighted_smeared, Form("GENIE v3, #chi^{2}/dof=%.1f/%d",chi2_m2,nbins),"lp");
  lg->AddEntry(g3_flux_weighted_smeared, Form("#muB Tuned, #chi^{2}/dof=%.1f/%d",chi2_m3,nbins), "l");
  lg->Draw();


  // ratio to ubTune Model
  auto c3 = new TCanvas("c3","c3",800,600);
  c3->cd();
  vector<double> xv2;
  for(int i=0; i<nbins; i++) {
    xv2.push_back(x_v.at(i) *1000.0);
    exl_v.at(i) *= 1000.0;
    exh_v.at(i) *= 1000.0;
  }

  auto gh_rt_data2ubmc = new TGraphAsymmErrors(nbins, xv2.data(), rt_data2ubmc.data(), exl_v.data(),exh_v.data(), err_rt_data2ubmc.data(), err_rt_data2ubmc.data());


  gh_rt_data2ubmc->SetTitle("");
  gh_rt_data2ubmc->GetXaxis()->SetTitle("E_{#nu} (MeV)");
  gh_rt_data2ubmc->GetYaxis()->SetTitle("Ratio to #muB Tuned model");
  gh_rt_data2ubmc->SetLineColor(4);
  gh_rt_data2ubmc->SetLineWidth(4);
  gh_rt_data2ubmc->SetMarkerColor(4);
  gh_rt_data2ubmc->SetMarkerStyle(20);
  gh_rt_data2ubmc->SetMarkerSize(2);
  gh_rt_data2ubmc->Draw("AP");


  // FC reco
  auto hdata_obsch_1 = (TH1F*)fmerge->Get("hdata_obsch_1");
  auto histo_3 = (TH1F*)fmerge->Get("histo_3"); // bkgd
  auto histo_5 = (TH1F*)fmerge->Get("histo_5"); // ext
  auto histo_1 = (TH1F*)fmerge->Get("histo_1"); // real nu

  hdata_obsch_1->Add(histo_3,-1);
  hdata_obsch_1->Add(histo_5,-1);
  hdata_obsch_1->Divide(histo_1);
  hdata_obsch_1->SetLineColor(1);
  hdata_obsch_1->SetMarkerColor(1);
  hdata_obsch_1->SetMarkerStyle(20);
  hdata_obsch_1->SetLineWidth(2);
  hdata_obsch_1->Draw("same");
  // PC reco
  auto hdata_obsch_2 = (TH1F*)fmerge->Get("hdata_obsch_2");
  auto histo_4 = (TH1F*)fmerge->Get("histo_4"); // bkgd
  auto histo_6 = (TH1F*)fmerge->Get("histo_6"); // ext
  auto histo_2 = (TH1F*)fmerge->Get("histo_2"); // real nu

  hdata_obsch_2->Add(histo_4,-1);
  hdata_obsch_2->Add(histo_6,-1);
  hdata_obsch_2->Divide(histo_2);
  hdata_obsch_2->SetLineColor(2);
  hdata_obsch_2->SetMarkerColor(2);
  hdata_obsch_2->SetMarkerStyle(20);
  hdata_obsch_2->SetLineWidth(2);
  hdata_obsch_2->Draw("same");


  auto lg1 = new TLegend(0.11,0.7,0.5,0.89);
  lg1->SetBorderSize(0);
  lg1->AddEntry(gh_rt_data2ubmc, "Unfolded Xs (NuMI)","lpe")->SetTextColor(4);
  lg1->AddEntry(hdata_obsch_1,"FC numuCC Data/MC in E^{reco}_{#nu}","lpe");
  lg1->AddEntry(hdata_obsch_2,"PC numuCC Data/MC in E^{reco}_{#nu}","lpe");
  lg1->Draw();


  auto ofile = new TFile("ccinc-numi-nue.root","recreate");
  gh_rt_data2ubmc->SetName("ccinc_numi_nue");  
  gh_rt_data2ubmc->Write();
  ofile->Close();


}

