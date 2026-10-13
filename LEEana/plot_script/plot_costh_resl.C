

void plot_costh_resl(){
  gStyle->SetErrorX(10);
  TFile *file1 =  new TFile("./checkout_prodgenie_bnb_nu_overlay/checkout_prodgenie_bnb_nu_overlay_run123.root");

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

  T_eval->SetAlias("reco_mu_costh", "T_PFeval.reco_muonMomentum[2] / sqrt( pow(T_PFeval.reco_muonMomentum[0],2) + pow(T_PFeval.reco_muonMomentum[1],2) + pow(T_PFeval.reco_muonMomentum[2],2) )");
  T_eval->SetAlias("truth_mu_costh", "T_PFeval.truth_muonMomentum[2] / sqrt( pow(T_PFeval.truth_muonMomentum[0],2) + pow(T_PFeval.truth_muonMomentum[1],2) + pow(T_PFeval.truth_muonMomentum[2],2) )");

  auto c2 = new TCanvas("c2","c2",1200,600);
  c2->Divide(2,1);

  c2->cd(1);
  T_eval->Draw("reco_mu_costh:truth_mu_costh >> h22(15,-1,1,15,-1,1)","weight_cv*weight_spline*(abs(truth_nuPdg)==14 && truth_isCC==1 && truth_vtxInside==1 && T_BDTvars.numu_cc_flag>=0 && T_BDTvars.numu_score>0.9 && T_PFeval.reco_muonMomentum[3]>0 && match_isFC>0)","colz");
  auto h22 = (TH2F*)gROOT->FindObject("h22");
  h22->GetXaxis()->SetTitle("True cos#theta_{#mu}");
  h22->GetYaxis()->SetTitle("Reco cos#theta_{#mu}");
  h22->GetYaxis()->SetTitleOffset(1.0);
  h22->SetTitle("FC");
  gPad->SetLogz();

  c2->cd(2);
  T_eval->Draw("reco_mu_costh:truth_mu_costh >> h23(15,-1,1,15,-1,1)","weight_cv*weight_spline*(abs(truth_nuPdg)==14 && truth_isCC==1 && truth_vtxInside==1 && T_BDTvars.numu_cc_flag>=0 && T_BDTvars.numu_score>0.9 && T_PFeval.reco_muonMomentum[3]>0 && match_isFC==0)","colz");
  auto h23 = (TH2F*)gROOT->FindObject("h23");
  h23->GetXaxis()->SetTitle("True cos#theta_{#mu}");
  h23->GetYaxis()->SetTitle("Reco cos#theta_{#mu}");
  h23->GetYaxis()->SetTitleOffset(1.0);
  h23->SetTitle("PC");
  gPad->SetLogz();

  auto c3 = new TCanvas("c3","c3",1200,600);
  c3->Divide(2,1);

  c3->cd(1);
  T_eval->Draw("reco_mu_costh - truth_mu_costh:truth_mu_costh >> h32(20,-1,1,50,-0.5,0.5)","weight_cv*weight_spline*(abs(truth_nuPdg)==14 && truth_isCC==1 && truth_vtxInside==1 && T_BDTvars.numu_cc_flag>=0 && T_BDTvars.numu_score>0.9 && T_PFeval.reco_muonMomentum[3]>0 && match_isFC>0)","colz");
  auto h32 = (TH2F*)gROOT->FindObject("h32");
  h32->GetXaxis()->SetTitle("True cos#theta_{#mu}");
  h32->GetYaxis()->SetTitle("(Reco - True) cos#theta_{#mu}");
  h32->GetYaxis()->SetTitleOffset(1.25);
  h32->SetTitle("FC");
  gPad->SetLogz();

  TGraphAsymmErrors* gEnuRes = new TGraphAsymmErrors();
  for(int i=1; i<=h32->GetNbinsX(); i++)
  {
      // if(i<3) continue;
      TH1D* htemp = h32->ProjectionY("htemp", i, i, "");
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
      xq = xq0 + (1-xq0)*0.6827;
      htemp->GetQuantiles(1,&par[2],&xq);
      // std::cout<<"higher quantile: "<<xq<<std::endl;
      // std::cout<<"higher rms: "<<par[2]-par[1]<<std::endl;
      // std::cout << par[2]-par[1] << std::endl;
      cout << h32->GetXaxis()->GetBinCenter(i) << endl;
      int N=gEnuRes->GetN();
      gEnuRes->SetPoint(N, h32->GetXaxis()->GetBinCenter(i), mean);
      gEnuRes->SetPointError(N, 0.5*h32->GetXaxis()->GetBinWidth(1),  0.5*h32->GetXaxis()->GetBinWidth(1), par[1]-par[0], par[2]-par[1] );
  }
  gEnuRes->SetLineWidth(3);
  gEnuRes->SetMarkerStyle(20);
  gEnuRes->Draw("P same");


  c3->cd(2);
  T_eval->Draw("reco_mu_costh - truth_mu_costh:truth_mu_costh >> h33(20,-1,1,50,-0.5,0.5)","weight_cv*weight_spline*(abs(truth_nuPdg)==14 && truth_isCC==1 && truth_vtxInside==1 && T_BDTvars.numu_cc_flag>=0 && T_BDTvars.numu_score>0.9 && T_PFeval.reco_muonMomentum[3]>0 && match_isFC==0)","colz");
  auto h33 = (TH2F*)gROOT->FindObject("h33");
  h33->GetXaxis()->SetTitle("True cos#theta_{#mu}");
  h33->GetYaxis()->SetTitle("(Reco - True) cos#theta_{#mu}");
  h33->GetYaxis()->SetTitleOffset(1.25);
  h33->SetTitle("PC");
  gPad->SetLogz();

  TGraphAsymmErrors* gEnuRes2 = new TGraphAsymmErrors();
  for(int i=1; i<=h33->GetNbinsX(); i++)
  {
      TH1D* htemp = h33->ProjectionY("htemp", i, i, "");
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
      xq = xq0 + (1-xq0)*0.6827;
      htemp->GetQuantiles(1,&par[2],&xq);
      // std::cout<<"higher quantile: "<<xq<<std::endl;
      // std::cout<<"higher rms: "<<par[2]-par[1]<<std::endl;
      // std::cout << par[2]-par[1] << std::endl;
      int N=gEnuRes2->GetN();
      gEnuRes2->SetPoint(N, h33->GetXaxis()->GetBinCenter(i), mean);
      gEnuRes2->SetPointError(N, 0.5*h33->GetXaxis()->GetBinWidth(1),  0.5*h33->GetXaxis()->GetBinWidth(1), par[1]-par[0], par[2]-par[1] );
  }
  gEnuRes2->SetLineWidth(3);
  gEnuRes2->SetMarkerStyle(20);
  gEnuRes2->Draw("P same");


}
