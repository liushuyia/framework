vector< vector<double> > get_bias_res(TH2F* hh, double threshold=0){
  // auto g1 = TGraphErrors();
  vector<double> gx,gy,gex,gey;

  int nx= hh->GetNbinsX();
  int ny= hh->GetNbinsY();
  for(int ix=0; ix<nx; ix++){
    double ymin = hh->GetYaxis()->GetXmin();
    double ymax = hh->GetYaxis()->GetXmax();
    TH1F hc("hc","hc",ny, ymin, ymax);
    for(int iy=0; iy<ny; iy++){
      double c = hh->GetBinContent(ix+1, iy+1);
      if (c>0) hc.SetBinContent(iy+1, c);
    }

    double mean=0, rms=0;
    double par[3]; 
    if (hc.GetEntries()>0) {
      // mean = hc.GetSum()/ hc.GetNbinsX(); 
      // rms = hc.GetRMS();
      // maxbin = hc.GetMaximumBin();
      // mean = hc.GetBinCenter(maxbin);

      // cout << "ix: " << ix << " hc.entries: " << hc.GetEntries() << endl;
      
      double xq = 0.5;
      hc.GetQuantiles(1,&par[1],&xq);
      xq = 0.5 + 0.34;
      hc.GetQuantiles(1,&par[0],&xq);
      xq = 0.5 - 0.34;
      hc.GetQuantiles(1,&par[2],&xq);
      par[2] = sqrt((pow(par[0]-par[1],2)+pow(par[2]-par[1],2))/2.);
      
      double trueE = hh->GetXaxis()->GetBinCenter(ix+1);
      mean = (par[1] - trueE ) / trueE;
      rms = par[2] / trueE;

      // g1.SetPoint(ix, hh->GetXaxis()->GetBinCenter(ix+1), mean);
      // g1.SetPointError(ix, 0, rms);
      // cout << ix << " " << hh->GetXaxis()->GetBinCenter(ix+1) << " " << mean << " " << rms << endl;
      gx.push_back(hh->GetXaxis()->GetBinCenter(ix+1));
      gy.push_back(mean);
      gex.push_back(0);
      gey.push_back(rms);
    }
  }
  // return g1;
  vector<vector<double> > ret;
  ret.push_back(gx);
  ret.push_back(gy);
  ret.push_back(gex);
  ret.push_back(gey);
  return ret;
}

vector< vector<double> > get_bias_maxbin(TH2F* hh, double threshold=0){
  // auto g1 = TGraphErrors();
  vector<double> gx,gy,gex,gey;

  int nx= hh->GetNbinsX();
  int ny= hh->GetNbinsY();
  for(int ix=0; ix<nx; ix++){
    double ymin = hh->GetYaxis()->GetXmin();
    double ymax = hh->GetYaxis()->GetXmax();
    TH1F hc("hc","hc",ny, ymin, ymax);
    for(int iy=0; iy<ny; iy++){
      double c = hh->GetBinContent(ix+1, iy+1);
      if (c>0) hc.SetBinContent(iy+1, c);
    }

    double mean=0, rms=0;
    if (hc.GetEntries()>0) {
      // mean = hc.GetSum()/ hc.GetNbinsX(); 
      // rms = hc.GetRMS();
      int maxbin = hc.GetMaximumBin();
      double xmean = hc.GetBinCenter(maxbin);
      
      // cout << "ix: " << ix << " hc.entries: " << hc.GetEntries() << endl;
      
      double trueE = hh->GetXaxis()->GetBinCenter(ix+1);
      mean = (xmean - trueE ) / trueE;

      // g1.SetPoint(ix, hh->GetXaxis()->GetBinCenter(ix+1), mean);
      // g1.SetPointError(ix, 0, rms);
      // cout << ix << " " << hh->GetXaxis()->GetBinCenter(ix+1) << " " << mean << " " << rms << endl;
      gx.push_back(hh->GetXaxis()->GetBinCenter(ix+1));
      gy.push_back(mean);
      gex.push_back(0);
      gey.push_back(0);
    }
  }
  // return g1;
  vector<vector<double> > ret;
  ret.push_back(gx);
  ret.push_back(gy);
  ret.push_back(gex);
  ret.push_back(gey);
  return ret;
}


void plot_mu_energy_resolution(){
  // TFile *file1 =  new TFile("./processed_checkout_rootfiles/checkout_prodgenie_bnb_nu_overlay_run1.root");
  TFile *file1 =  new TFile("checkout_prodgenie_bnb_nu_overlay/checkout_prodgenie_bnb_nu_overlay_run123.root");

  TTree *T_eval = (TTree*)file1->Get("wcpselection/T_eval");
  TTree *T_BDTvars = (TTree*)file1->Get("wcpselection/T_BDTvars");
  TTree *T_KINEvars = (TTree*)file1->Get("wcpselection/T_KINEvars");
  TTree *T_PFeval = (TTree*)file1->Get("wcpselection/T_PFeval");
  T_eval->AddFriend(T_BDTvars,"T_BDTvars");
  T_eval->AddFriend(T_KINEvars,"T_KINEvars");
  T_eval->AddFriend(T_PFeval, "T_PFeval");
  T_eval->SetAlias("truth_muEnergy","1000.0*T_PFeval.truth_muonMomentum[3]");
  T_eval->SetAlias("truth_hadEnergy","truth_nuEnergy - 1000.0*T_PFeval.truth_muonMomentum[3]");
  T_eval->SetAlias("reco_muEnergy","1000.0*T_PFeval.reco_muonMomentum[3]");
  T_eval->SetAlias("reco_hadEnergy","T_KINEvars.kine_reco_Enu - 1000.0*T_PFeval.reco_muonMomentum[3]");
  T_eval->SetAlias("reco_nuEnergy","T_KINEvars.kine_reco_Enu");

  auto c1 = new TCanvas("c1","c1",1200,900);
  c1->Divide(4,3);

  // double xbins1[] = {0.2, 0.540, 0.705, 0.805, 0.920, 1.050, 1.200, 1.375, 1.570, 2.050, 4.000}; // GeV, neutrino energy, 10 bins
  double xbins1[] = {0.106, 0.226, 0.296, 0.386, 0.505, 0.577, 0.659, 0.753, 0.861, 0.984, 1.285, 2.506}; // GeV, muon energy, 11bins
  // double xbins1[] = {0.03, 0.1, 0.15, 0.225, 0.275, 0.336, 0.411, 0.502, 0.614, 0.75, 1.12, 2.5}; // GeV, nu = enu - emu, 11bins
  
  vector<TH1F*> hfc(25);
  vector<TH1F*> hpc(25);
  for(int i=0; i<11; i++) {
    c1->cd(i+1);
    double etrue1 = xbins1[i] * 1000.0;
    double etrue2 = xbins1[i+1] * 1000.0;
    // double reco1 = i*100.0;
    // double reco2 = (i+1)*100.0;

    hfc.at(i) = new TH1F(Form("hfc%d",i),Form("hfc%d",i),25,0,2.5);
    hfc.at(i)->SetLineColor(1);
    hpc.at(i) = new TH1F(Form("hpc%d",i),Form("hpc%d",i),25,0,2.5);
    hpc.at(i)->SetLineColor(2);

    T_eval->Project(Form("hfc%d",i), "1e-3*reco_muEnergy", Form("weight_cv*weight_spline*(abs(truth_nuPdg)==14 && truth_isCC==1 && truth_vtxInside==1 && T_BDTvars.numu_cc_flag>=0 && T_BDTvars.numu_score>0.9 && match_isFC==1 && truth_muEnergy>%f && truth_muEnergy<%f)", etrue1, etrue2) );
    T_eval->Project(Form("hpc%d",i), "1e-3*reco_muEnergy", Form("weight_cv*weight_spline*(abs(truth_nuPdg)==14 && truth_isCC==1 && truth_vtxInside==1 && T_BDTvars.numu_cc_flag>=0 && T_BDTvars.numu_score>0.9 && match_isFC==0 && truth_muEnergy>%f && truth_muEnergy<%f)", etrue1, etrue2) );

    // hfc.at(i)->SetTitle(Form("True energy bin: %.3f - %.3f", 1e-3*etrue1, 1e-3*etrue2));
    hfc.at(i)->SetTitle("");
    hfc.at(i)->GetXaxis()->SetTitle("Reco Muon Energy (GeV)");
    hfc.at(i)->GetYaxis()->SetTitle("Entries");
    hfc.at(i)->SetTitleSize(0.1);
    hfc.at(i)->GetXaxis()->SetTitleSize(0.06);
    hfc.at(i)->Draw("hist");
    hpc.at(i)->Draw("hist same");


    double hmax = hfc.at(i)->GetMaximum()>hpc.at(i)->GetMaximum() ? hfc.at(i)->GetMaximum()*1.05 : hpc.at(i)->GetMaximum()*1.05;
    hfc.at(i)->GetYaxis()->SetRangeUser(0, hmax);
    hpc.at(i)->GetYaxis()->SetRangeUser(0, hmax);

    // auto l1 = new TLine(etrue1*1e-3, 0, etrue1*1e-3, hmax);
    // auto l2 = new TLine(etrue2*1e-3, 0, etrue2*1e-3, hmax);
    // l1->SetLineColor(4);
    // l2->SetLineColor(4);
    // l1->SetLineStyle(2);
    // l2->SetLineStyle(2);
    // l1->SetLineWidth(2);
    // l2->SetLineWidth(2);
    // l1->Draw("same"); 
    // l2->Draw("same"); 
    auto b1 = new TBox(etrue1*1e-3, 0, etrue2*1e-3, hmax);
    b1->SetFillStyle(3004);
    b1->SetFillColor(4);
    b1->Draw("same");
  }

  // auto l1 = new TLine(0,1,0,2);
  // l1->SetLineColor(4);
  // l1->SetLineStyle(2);
  auto b1 = new TBox(0,1,0,2);
  b1->SetFillStyle(3004);
  b1->SetFillColor(4);

  c1->cd(12);
  auto lg = new TLegend(0.1,0.3,0.8,0.8);
  lg->AddEntry(hfc.at(0), "FC","l");
  lg->AddEntry(hpc.at(0), "PC","l")->SetTextColor(2);
  lg->AddEntry(b1, "E^{true}_{#mu} Bins","f")->SetTextColor(4);
  lg->Draw();

}
