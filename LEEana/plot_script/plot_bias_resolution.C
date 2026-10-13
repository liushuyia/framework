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


void plot_bias_resolution(){
  // TFile *file1 =  new TFile("./processed_checkout_rootfiles/checkout_prodgenie_bnb_nu_overlay_run1.root");
  TFile *file1 =  new TFile("./checkout_prodgenie_bnb_nu_overlay/checkout_prodgenie_bnb_nu_overlay_run123.root");

  // TH1F *h10 = new TH1F("h10","h10",25,0,2500);
  // TH1F *h11 = new TH1F("h11","h11",25,0,2500);
  // TH1F *h12 = new TH1F("h12","h12",25,0,2500);

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

  auto c1 = new TCanvas("c1","c1",1800,600);
  c1->Divide(3,1);

  c1->cd(1);
  T_eval->Draw("reco_nuEnergy*1e-3 >> h1(30,0.110,3)","weight_cv*weight_spline*(abs(truth_nuPdg)==14 && truth_isCC==1 && truth_vtxInside==1 && T_BDTvars.numu_cc_flag>=0 && T_BDTvars.numu_score>0.9 && match_isFC>0)","hist");
  auto h1 = (TH2F*)gROOT->FindObject("h1");
  h1->GetXaxis()->SetTitle("Reco E_{#nu} [GeV]");
  h1->SetTitle("");

  c1->cd(2);
  T_eval->Draw("reco_muEnergy*1e-3 >> h2(30,0.106,3)","weight_cv*weight_spline*(abs(truth_nuPdg)==14 && truth_isCC==1 && truth_vtxInside==1 && T_BDTvars.numu_cc_flag>=0 && T_BDTvars.numu_score>0.9 && match_isFC>0)","hist");
  auto h2 = (TH2F*)gROOT->FindObject("h2");
  h2->GetXaxis()->SetTitle("Reco E_{#mu} [GeV]");
  h2->SetTitle("");

  c1->cd(3);
  T_eval->Draw("reco_hadEnergy*1e-3 >> h3(30,0,3)","weight_cv*weight_spline*(abs(truth_nuPdg)==14 && truth_isCC==1 && truth_vtxInside==1 && T_BDTvars.numu_cc_flag>=0 && T_BDTvars.numu_score>0.9 && match_isFC>0 && reco_muEnergy>0)","hist");
  auto h3 = (TH2F*)gROOT->FindObject("h3");
  h3->GetXaxis()->SetTitle("Reco E_{had} [GeV]");
  h3->SetTitle("");

  auto c2 = new TCanvas("c2","c2",1800,600);
  c2->Divide(3,1);

  c2->cd(1);
  T_eval->Draw("reco_nuEnergy*1e-3:truth_nuEnergy*1e-3 >> h21(60,0.110,3,60,0.110,3)","weight_cv*weight_spline*(abs(truth_nuPdg)==14 && truth_isCC==1 && truth_vtxInside==1 && T_BDTvars.numu_cc_flag>=0 && T_BDTvars.numu_score>0.9 && match_isFC>0)","colz");
  auto h21 = (TH2F*)gROOT->FindObject("h21");
  h21->GetXaxis()->SetTitle("True E_{#nu} [GeV]");
  h21->GetYaxis()->SetTitle("Reco E_{#nu} [GeV]");
  h21->SetTitle("");
  gPad->SetLogz();
  auto l1 = new TLine(0,0,3,3);
  l1->SetLineColor(1);
  l1->SetLineStyle(2);
  l1->SetLineWidth(2);
  l1->Draw();

//
  auto c2x = new TCanvas("c2x","c2x",600,450);
  // T_eval->Draw("reco_nuEnergy/truth_nuEnergy:truth_nuEnergy*1e-3 >> h21x(30,0.110,3, 50,0,2)","weight_cv*weight_spline*(abs(truth_nuPdg)==14 && truth_isCC==1 && truth_vtxInside==1 && T_BDTvars.numu_cc_flag>=0 && T_BDTvars.numu_score>0.9 && match_isFC>0)","colz");
  T_eval->Draw("reco_nuEnergy/truth_nuEnergy:truth_nuEnergy*1e-3 >> h21x(30,0,3, 50,0,2)","weight_cv*weight_spline*(abs(truth_nuPdg)==14 && truth_isCC==1 && truth_vtxInside==1 && T_BDTvars.numu_cc_flag>=0 && T_BDTvars.numu_score>0.9 && match_isFC>0)","colz");
  auto h21x = (TH2F*)gROOT->FindObject("h21x");
  h21x->GetXaxis()->SetTitle("True E_{#nu} [GeV]");
  h21x->GetYaxis()->SetTitle("Reco E_{#nu} / True E_{#nu}");
  h21x->SetTitle("");
  // gPad->SetLogz();
  // auto l1 = new TLine(0,0,3,3);
  // l1->SetLineColor(1);
  // l1->SetLineStyle(2);
  // l1->SetLineWidth(2);
  // l1->Draw();

  TGraphAsymmErrors* gEnuRes = new TGraphAsymmErrors();
  for(int i=1; i<=h21x->GetNbinsX(); i++)
  {
      if(i<3) continue;
      TH1D* htemp = h21x->ProjectionY("htemp", i, i, "");
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
      std::cout << h21x->GetXaxis()->GetBinCenter(i) << " " << par[1]-par[0] << " " ;
      xq = xq0 + (1-xq0)*0.6827;
      htemp->GetQuantiles(1,&par[2],&xq);
      // std::cout<<"higher quantile: "<<xq<<std::endl;
      // std::cout<<"higher rms: "<<par[2]-par[1]<<std::endl;
      std::cout << par[2]-par[1] << std::endl;
      gEnuRes->SetPoint(i-2, h21x->GetXaxis()->GetBinCenter(i), mean);
      gEnuRes->SetPointError(i-2, 0.5*h21x->GetXaxis()->GetBinWidth(1),  0.5*h21x->GetXaxis()->GetBinWidth(1), par[1]-par[0], par[2]-par[1] );
  }
  h21x->Draw("colz");
  gEnuRes->Draw("P same");
   
//

  c2->cd(2);
  T_eval->Draw("reco_muEnergy*1e-3:truth_muEnergy*1e-3 >> h22(60,0.106,3,60,0.106,3)","weight_cv*weight_spline*(abs(truth_nuPdg)==14 && truth_isCC==1 && truth_vtxInside==1 && T_BDTvars.numu_cc_flag>=0 && T_BDTvars.numu_score>0.9 && match_isFC>0)","colz");
  auto h22 = (TH2F*)gROOT->FindObject("h22");
  h22->GetXaxis()->SetTitle("True E_{#mu} [GeV]");
  h22->GetYaxis()->SetTitle("Reco E_{#mu} [GeV]");
  h22->SetTitle("");
  gPad->SetLogz();
  l1->Draw();

  c2->cd(3);
  T_eval->Draw("reco_hadEnergy*1e-3:truth_hadEnergy*1e-3 >> h23(30,0,3,30,0,3)","weight_cv*weight_spline*(abs(truth_nuPdg)==14 && truth_isCC==1 && truth_vtxInside==1 && T_BDTvars.numu_cc_flag>=0 && T_BDTvars.numu_score>0.9 && match_isFC>0 && reco_muEnergy>0)","colz");
  auto h23 = (TH2F*)gROOT->FindObject("h23");
  h23->GetXaxis()->SetTitle("True E_{trans} [GeV]");
  h23->GetYaxis()->SetTitle("Reco E_{had}  [GeV]");
  h23->SetTitle("");
  gPad->SetLogz();
  l1->Draw();

  auto c3 = new TCanvas("c3","c3",1800,600);
  c3->Divide(3,1);

  c3->cd(1);
  auto vv1 = get_bias_res(h21);
  auto g1 = new TGraphErrors(vv1.at(0).size(), vv1.at(0).data(), vv1.at(1).data(), vv1.at(2).data(), vv1.at(3).data());
  g1->SetTitle("");
  g1->GetXaxis()->SetTitle("True E_{#nu} [GeV]");
  g1->SetMarkerColor(1);
  g1->SetMarkerStyle(20);
  g1->SetMarkerSize(1);
  g1->SetMaximum(1);
  g1->SetMinimum(-1);
  g1->GetXaxis()->SetRangeUser(0.110, 3);
  g1->Draw("AP");
  auto vvv1 = get_bias_maxbin(h21);
  auto g12 = new TGraphErrors(vvv1.at(0).size(), vvv1.at(0).data(), vvv1.at(1).data(), vvv1.at(2).data(), vvv1.at(3).data());
  g12->SetMarkerColor(2);
  g12->SetMarkerStyle(20);
  g12->SetMarkerSize(1);
  g12->SetLineColor(2);
  g12->Draw("PLsame");


  c3->cd(2);
  auto vv2 = get_bias_res(h22);
  auto g2 = new TGraphErrors(vv2.at(0).size(), vv2.at(0).data(), vv2.at(1).data(), vv2.at(2).data(), vv2.at(3).data());
  g2->SetTitle("");
  g2->GetXaxis()->SetTitle("True E_{#mu} [GeV]");
  g2->SetMarkerColor(1);
  g2->SetMarkerStyle(20);
  g2->SetMarkerSize(1);
  g2->SetMaximum(1);
  g2->SetMinimum(-1);
  g2->GetXaxis()->SetRangeUser(0.106, 3);
  g2->Draw("AP");
  auto vvv2 = get_bias_maxbin(h22);
  auto g22 = new TGraphErrors(vvv2.at(0).size(), vvv2.at(0).data(), vvv2.at(1).data(), vvv2.at(2).data(), vvv2.at(3).data());
  g22->SetMarkerColor(2);
  g22->SetMarkerStyle(20);
  g22->SetMarkerSize(1);
  g22->SetLineColor(2);
  g22->Draw("PLsame");


  c3->cd(3);
  auto vv3 = get_bias_res(h23);
  auto g3 = new TGraphErrors(vv3.at(0).size(), vv3.at(0).data(), vv3.at(1).data(), vv3.at(2).data(), vv3.at(3).data());
  g3->SetTitle("");
  g3->GetXaxis()->SetTitle("True E_{trans} [GeV]");
  g3->SetMarkerColor(1);
  g3->SetMarkerStyle(20);
  g3->SetMarkerSize(1);
  g3->SetMaximum(1);
  g3->SetMinimum(-1);
  g3->GetXaxis()->SetRangeUser(0., 3);
  g3->Draw("AP");
  auto vvv3 = get_bias_maxbin(h23);
  auto g32 = new TGraphErrors(vvv3.at(0).size(), vvv3.at(0).data(), vvv3.at(1).data(), vvv3.at(2).data(), vvv3.at(3).data());
  g32->SetMarkerColor(2);
  g32->SetMarkerStyle(20);
  g32->SetMarkerSize(1);
  g32->SetLineColor(2);
  g32->Draw("PLsame");


}
