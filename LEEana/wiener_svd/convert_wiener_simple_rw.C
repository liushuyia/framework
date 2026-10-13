void convert_wiener_simple_rw(){
  gStyle->SetOptStat(0);
  
  TFile *file1 = new TFile("merge_xs.root");
  TMatrixD *mat_collapse = (TMatrixD*)file1->Get("mat_collapse"); // collapse
  TMatrixD *cov_mat_add = (TMatrixD*)file1->Get("cov_mat_add"); // additional uncertainties
  TVectorD *vec_signal = (TVectorD*)file1->Get("vec_signal"); // xs in pbar
  TMatrixD *mat_R = (TMatrixD*)file1->Get("mat_R"); // smearing matrix

  TH1F *hdata_obsch_1 = (TH1F*)file1->Get("hdata_obsch_1"); // BNB data FC
  TH1F *hdata_obsch_2 = (TH1F*)file1->Get("hdata_obsch_2"); // BNB data PC
  TH1F *hmc_obsch_1 = (TH1F*)file1->Get("hmc_obsch_1");  // Total pred FC
  TH1F *hmc_obsch_2 = (TH1F*)file1->Get("hmc_obsch_2"); // Total pred PC
  TH1F *histo_1 = (TH1F*)file1->Get("histo_1"); // MC signal FC
  TH1F *histo_2 = (TH1F*)file1->Get("histo_2"); // MC signal PC
  TH1F *histo_3 = (TH1F*)file1->Get("histo_3"); // MC bkg FC
  TH1F *histo_4 = (TH1F*)file1->Get("histo_4"); // MC bkg PC
  TH1F *histo_5 = (TH1F*)file1->Get("histo_5"); // EXT FC
  TH1F *histo_6 = (TH1F*)file1->Get("histo_6"); // EXT PC

  int nbin_true = vec_signal->GetNoElements();

  std::vector<double> cv_data_obsch_1;  std::vector<double> err_data_obsch_1;
  std::vector<double> cv_data_obsch_2;  std::vector<double> err_data_obsch_2;
  std::vector<double> cv_mc_obsch_1;    std::vector<double> err_mc_obsch_1;
  std::vector<double> cv_mc_obsch_2;    std::vector<double> err_mc_obsch_2;
  std::vector<double> cv_histo_1;       std::vector<double> err_histo_1;
  std::vector<double> cv_histo_2;       std::vector<double> err_histo_2;
  std::vector<double> cv_histo_3;       std::vector<double> err_histo_3;
  std::vector<double> cv_histo_4;       std::vector<double> err_histo_4;
  std::vector<double> cv_histo_5;       std::vector<double> err_histo_5;
  std::vector<double> cv_histo_6;       std::vector<double> err_histo_6;

  for(int i=0;i!=nbin_true;i++) cv_data_obsch_1.push_back(hdata_obsch_1->GetBinContent(i+1));
  for(int i=0;i!=nbin_true;i++) cv_data_obsch_2.push_back(hdata_obsch_2->GetBinContent(i+1));
  for(int i=0;i!=nbin_true;i++) cv_mc_obsch_1.push_back(hmc_obsch_1->GetBinContent(i+1));
  for(int i=0;i!=nbin_true;i++) cv_mc_obsch_2.push_back(hmc_obsch_2->GetBinContent(i+1));
  for(int i=0;i!=nbin_true;i++) cv_histo_1.push_back(histo_1->GetBinContent(i+1));
  for(int i=0;i!=nbin_true;i++) cv_histo_2.push_back(histo_2->GetBinContent(i+1));
  for(int i=0;i!=nbin_true;i++) cv_histo_3.push_back(histo_3->GetBinContent(i+1));
  for(int i=0;i!=nbin_true;i++) cv_histo_4.push_back(histo_4->GetBinContent(i+1));
  for(int i=0;i!=nbin_true;i++) cv_histo_5.push_back(histo_5->GetBinContent(i+1));
  for(int i=0;i!=nbin_true;i++) cv_histo_6.push_back(histo_6->GetBinContent(i+1));

  // true signal ...
  TH1D *htrue_signal = new TH1D("htrue_signal","htrue_signal",nbin_true, 0.5, nbin_true+0.5);
  for(int i=0;i!=nbin_true;i++){
    htrue_signal->SetBinContent(i+1,(*vec_signal)(i));
  }
  Int_t nbin_meas = hdata_obsch_1->GetNbinsX() + 1 + hdata_obsch_2->GetNbinsX() + 1;

  // Cross section Numerator for data and prediction
  TH1D *hmeas = new TH1D("hmeas","hmeas",nbin_meas,0.5, nbin_meas+0.5);
  TH1D *hpred = new TH1D("hpred","hpred",nbin_meas,0.5, nbin_meas+0.5);
  for (Int_t i=0;i!=hdata_obsch_1->GetNbinsX()+1;i++){
    hmeas->SetBinContent(i+1,                                   hdata_obsch_1->GetBinContent(i+1) - histo_3->GetBinContent(i+1) - histo_5->GetBinContent(i+1) );
    hmeas->SetBinContent(i+1 + hdata_obsch_1->GetNbinsX() + 1 , hdata_obsch_2->GetBinContent(i+1) - histo_4->GetBinContent(i+1) - histo_6->GetBinContent(i+1) );

    //std::cout << hmc_obsch_2->GetBinContent(i+1) -  histo_4->GetBinContent(i+1) - histo_6->GetBinContent(i+1) << " " << histo_2->GetBinContent(i+1) << std::endl;

    hpred->SetBinContent(i+1,                                   hmc_obsch_1->GetBinContent(i+1) );
    hpred->SetBinContent(i+1 + hdata_obsch_1->GetNbinsX() + 1 , hmc_obsch_2->GetBinContent(i+1) );
  }
  // response matrix
  
  TMatrixD mat_collapse_T(mat_collapse->GetNcols(), mat_collapse->GetNrows()); 
  mat_collapse_T.Transpose(*mat_collapse);
  TMatrixD mat_R_collapse = (mat_collapse_T)*(*mat_R);
  
  ////mat_R_collapse.Draw("COLZ");

  TH2D* hR = new TH2D("hR","hR",nbin_meas,0.5,nbin_meas+0.5, nbin_true, 0.5, nbin_true+0.5);
  for (Int_t i=0;i!=nbin_true;i++){
    for (Int_t j=0;j!=nbin_meas;j++){
      hR->SetBinContent(j+1,i+1,mat_R_collapse(j,i));
    }
  }
  ////  hR->Draw("COLZ");
  
  // covariance matrix
  TH2D *hcov_tot = new TH2D("hcov_tot","hcov_tot",nbin_meas,0.5,nbin_meas+0.5,nbin_meas,0.5,nbin_meas+0.5);
  // MC statistical uncertainties ...
  TH2D *hcov_mcstat = new TH2D("hcov_mcstat","hcov_mcstat",nbin_meas,0.5,nbin_meas+0.5,nbin_meas,0.5,nbin_meas+0.5);
  ifstream infile("mcstat/xs_tot.log");
  double temp, temp1,err2;
  infile >> temp >> temp;
  for (Int_t i=0;i!=nbin_meas;i++){
    infile >> temp >> temp >> temp >> err2 >> temp;
    hcov_mcstat->SetBinContent(i+1,i+1,err2);
  }
  ////hcov_mcstat->Draw("COLZ");

  // statistical uncertainties
  TH2D *hcov_stat = new TH2D("hcov_stat","hcov_stat",nbin_meas,0.5,nbin_meas+0.5,nbin_meas,0.5,nbin_meas+0.5);
  for (Int_t i=0;i!=hdata_obsch_1->GetNbinsX()+1;i++){
    double meas = hdata_obsch_1->GetBinContent(i+1);
    double pred = hmc_obsch_1->GetBinContent(i+1);
    double content;
    if (pred !=0){
      if (meas == 0){ content = pred/2.;
      }else{ content = 3./(1./meas+2./pred);
      }
    }else{ content = 0;
    }
    hcov_stat->SetBinContent(i+1,i+1,content);
  }
  for (Int_t i=0;i!=hdata_obsch_2->GetNbinsX()+1;i++){
    double meas = hdata_obsch_2->GetBinContent(i+1);
    double pred = hmc_obsch_2->GetBinContent(i+1);
    double content;
    if (pred !=0){
      if (meas == 0){ content = pred/2.;
      }else{ content = 3./(1./meas+2./pred);
      }
    }else{ content = 0;
    }
    hcov_stat->SetBinContent(hdata_obsch_1->GetNbinsX()+1+i+1,hdata_obsch_1->GetNbinsX()+1+i+1,content);
  }
  //hcov_stat->Draw("COLZ");
  
  // additional uncertainty
  TH2D *hcov_add = new TH2D("hcov_add","hcov_add",nbin_meas,0.5,nbin_meas+0.5,nbin_meas,0.5,nbin_meas+0.5);
  TMatrixD mat_add = mat_collapse_T * (*cov_mat_add) * (*mat_collapse);
  ////  mat_add.Draw("COLZ");
  for (Int_t i=0;i!=nbin_meas;i++){
    hcov_add->SetBinContent(i+1,i+1,mat_add(i,i));
  }
  //hcov_add->Draw("COLZ");
  
  // detector systematics
  TH2D *hcov_det = new TH2D("hcov_det","hcov_det",nbin_meas,0.5,nbin_meas+0.5,nbin_meas,0.5,nbin_meas+0.5);
  int nold = mat_collapse->GetNrows();
  std::cout << nold << std::endl;
  TVectorD vec_nominal(nold);
  int noffset = 0;
  for (Int_t i=0;i!=histo_1->GetNbinsX()+1;i++){
    vec_nominal(noffset + i) = histo_1->GetBinContent(i+1);
  }
  noffset += histo_1->GetNbinsX()+1;
  for (Int_t i=0;i!=histo_2->GetNbinsX()+1;i++){
    vec_nominal(noffset + i) = histo_2->GetBinContent(i+1);
  }
  noffset += histo_2->GetNbinsX()+1;
  for (Int_t i=0;i!=histo_3->GetNbinsX()+1;i++){
    vec_nominal(noffset + i) = histo_3->GetBinContent(i+1);
  }
  noffset += histo_3->GetNbinsX()+1;
  for (Int_t i=0;i!=histo_4->GetNbinsX()+1;i++){
    vec_nominal(noffset + i) = histo_4->GetBinContent(i+1);
  }
  noffset += histo_4->GetNbinsX()+1;
  for (Int_t i=0;i!=histo_5->GetNbinsX()+1;i++){
    vec_nominal(noffset + i) = histo_5->GetBinContent(i+1);
  }
  noffset += histo_5->GetNbinsX()+1;
  for (Int_t i=0;i!=histo_6->GetNbinsX()+1;i++){
    vec_nominal(noffset + i) = histo_6->GetBinContent(i+1);
  }
  noffset += histo_6->GetNbinsX()+1;
  //vec_nominal.Draw();
  std::map<int, TString> map_det_str;
  map_det_str[1] = "./DetVar/cov_LYDown.root";
  map_det_str[2] = "./DetVar/cov_LYRayleigh.root";
  map_det_str[3] = "./DetVar/cov_Recomb2.root";
  map_det_str[4] = "./DetVar/cov_SCE.root";
  // map_det_str[5] = "./DetVar/cov_WMdEdx.root";
  map_det_str[6] = "./DetVar/cov_WMThetaXZ.root";
  map_det_str[7] = "./DetVar/cov_WMThetaYZ.root";
  map_det_str[8] = "./DetVar/cov_WMX.root";
  map_det_str[9] = "./DetVar/cov_WMYZ.root";
  map_det_str[10] = "./DetVar/cov_LYatt.root";
  TMatrixD frac_det(nold,nold);
  for (auto it = map_det_str.begin(); it!=map_det_str.end(); it++){
    int idx = it->first;
    TFile tmp_file(it->second);
    TMatrixD *tmp_matrix = (TMatrixD*)tmp_file.Get(Form("frac_cov_det_mat_%d",idx));
    frac_det += (*tmp_matrix);
  }

  for (Int_t i=0;i!=nold;i++){
    for (Int_t j=0;j!=nold;j++){
      frac_det(i,j) *= vec_nominal(i) * vec_nominal(j);
    }
  }
  ////frac_det.Draw("COLZ");
  TMatrixD mat_det = mat_collapse_T * frac_det * (*mat_collapse);
  ////mat_det.Draw("COLZ");
  for (Int_t i=0;i!=nbin_meas;i++){
    for (Int_t j=0;j!=nbin_meas;j++){
      hcov_det->SetBinContent(i+1,j+1,mat_det(i,j));
    }
  }
  //hcov_det->Draw("COLZ");
  
  // Flux systematics
  TH2D *hcov_flux = new TH2D("hcov_flux","hcov_flux",nbin_meas,0.5,nbin_meas+0.5,nbin_meas,0.5,nbin_meas+0.5);
  TMatrixD frac_flux(nold,nold);
  for(Int_t i=1;i!=17;i++){
    TFile tmp_file(Form("./XsFlux/cov_%d.root",i));
    TMatrixD *tmp_matrix = (TMatrixD*)tmp_file.Get(Form("frac_cov_xf_mat_%d",i));
    frac_flux += *tmp_matrix;
  }
  for (Int_t i=0;i!=nold;i++){
    for (Int_t j=0;j!=nold;j++){
      frac_flux(i,j) *= vec_nominal(i) * vec_nominal(j);
    }
  }
  ////frac_flux.Draw("COLZ");
  TMatrixD mat_flux = mat_collapse_T * frac_flux * (*mat_collapse);
  ////mat_flux.Draw("COLZ");
  for (Int_t i=0;i!=nbin_meas;i++){
    for (Int_t j=0;j!=nbin_meas;j++){
      hcov_flux->SetBinContent(i+1,j+1,mat_flux(i,j));
    }
  }
  //hcov_flux->Draw("COLZ");

  // Xs systematics
  TH2D *hcov_xs = new TH2D("hcov_xs","hcov_xs",nbin_meas,0.5,nbin_meas+0.5,nbin_meas,0.5,nbin_meas+0.5);
  {
    TFile tmp_file("./XsFlux/cov_xs.root");   // cross section mode
    //TFile tmp_file("./XsFlux/cov_17.root"); // nominal 
    TMatrixD *frac_xs = (TMatrixD*)tmp_file.Get("frac_cov_xf_mat_17");
    
    for (Int_t i=0;i!=nold;i++){
      for (Int_t j=0;j!=nold;j++){
        (*frac_xs)(i,j) *= vec_nominal(i) * vec_nominal(j);
      }
    }
    ////frac_flux.Draw("COLZ");
    TMatrixD mat_xs = mat_collapse_T * (*frac_xs) * (*mat_collapse);
    for (Int_t i=0;i!=nbin_meas;i++){
      for (Int_t j=0;j!=nbin_meas;j++){
        hcov_xs->SetBinContent(i+1,j+1,mat_xs(i,j));
      }
    }
    //hcov_xs->Draw("COLZ");
  }
  // rw uncor sys
 TH2D *hcov_rw = new TH2D("hcov_rw","hcov_rw",nbin_meas,0.5,nbin_meas+0.5,nbin_meas,0.5,nbin_meas+0.5);
  {
    TFile tmp_file("./XsFlux/cov_rw.root");   // cross section mode
    TMatrixD *frac_rw = (TMatrixD*)tmp_file.Get("frac_cov_xf_mat_18");

    for (Int_t i=0;i!=nold;i++){
      for (Int_t j=0;j!=nold;j++){
        (*frac_rw)(i,j) *= vec_nominal(i) * vec_nominal(j);
      }
    }
    ////frac_flux.Draw("COLZ");
    TMatrixD mat_rw = mat_collapse_T * (*frac_rw) * (*mat_collapse);
    for (Int_t i=0;i!=nbin_meas;i++){
      for (Int_t j=0;j!=nbin_meas;j++){
        hcov_rw->SetBinContent(i+1,j+1,mat_rw(i,j));
      }
    }
    //hcov_xs->Draw("COLZ");
  } 

 // rw cor sys
 TH2D *hcov_rw_cor = new TH2D("hcov_rw_cor","hcov_rw_cor",nbin_meas,0.5,nbin_meas+0.5,nbin_meas,0.5,nbin_meas+0.5);
  {
    TFile tmp_file("./XsFlux/cov_rw_cor.root");   // cross section mode
    TMatrixD *frac_rw_cor = (TMatrixD*)tmp_file.Get("frac_cov_xf_mat_19");

    for (Int_t i=0;i!=nold;i++){
      for (Int_t j=0;j!=nold;j++){
        (*frac_rw_cor)(i,j) *= vec_nominal(i) * vec_nominal(j);
      }
    }
    ////frac_flux.Draw("COLZ");
    TMatrixD mat_rw_cor = mat_collapse_T * (*frac_rw_cor) * (*mat_collapse);
    for (Int_t i=0;i!=nbin_meas;i++){
      for (Int_t j=0;j!=nbin_meas;j++){
        hcov_rw_cor->SetBinContent(i+1,j+1,mat_rw_cor(i,j));
      }
    }
    //hcov_xs->Draw("COLZ");
  } 

  vec_nominal.Draw();

  hcov_tot->Add(hcov_stat);
  hcov_tot->Add(hcov_mcstat);
  hcov_tot->Add(hcov_add);
  hcov_tot->Add(hcov_flux);
  hcov_tot->Add(hcov_det);
  hcov_tot->Add(hcov_xs);
  hcov_tot->Add(hcov_rw);
  hcov_tot->Add(hcov_rw_cor);

  //hcov_tot->Draw("COLZ");

  int reco_bin_1 = hdata_obsch_1->GetNbinsX();
  int reco_bin_2 = hdata_obsch_1->GetNbinsX() + 1;
  int reco_bin_3 = reco_bin_2 + hdata_obsch_2->GetNbinsX();
  int reco_bin_4 = reco_bin_2 + hdata_obsch_2->GetNbinsX() + 1;

  ofstream myfile;
  myfile.open ("values.txt");

  std::cout << "----------------------------------------------------------------------------" << std::endl;
  for (Int_t i=0;i!=nbin_meas;i++){
    if(i==0){
      myfile << "                  Reconstructed space fractional uncertainties (diagonal terms / prediction)" << std::endl;
      myfile << "# BNB  MCBKG   EXT    NUM     MCSIG  err_tot  err_stat  err_mcstat  err_add  err_flux  err_det  err_xs  err_rw  err_rw_cor" << std::endl;
      std::cout << "                  Reconstructed space fractional uncertainties (diagonal terms / prediction)" << std::endl;
      std::cout << "# BNB  MCBKG   EXT    NUM     MCSIG  err_tot  err_stat  err_mcstat  err_add  err_flux  err_det  err_xs  err_rw  err_rw_cor" << std::endl;
      std::cout << "----------------------------------------------------------------------------" << std::endl;
    }
    if(i==reco_bin_1) std::cout << "----------------------------------------------------------------------------" << std::endl;
    if(i==reco_bin_2) std::cout << "----------------------------------------------------------------------------" << std::endl;
    if(i==reco_bin_3) std::cout << "----------------------------------------------------------------------------" << std::endl;
    if(i==reco_bin_4) std::cout << "----------------------------------------------------------------------------" << std::endl;     

    //if (hpred->GetBinContent(i+1)==0) { cout << i << endl; continue; }

    if(i<reco_bin_2){
      std::cout << i << " " << hdata_obsch_1->GetBinContent(i+1) << " " << histo_3->GetBinContent(i+1) << " " << histo_5->GetBinContent(i+1) << " " << hmeas->GetBinContent(i+1) << " " << hpred->GetBinContent(i+1) << " "
                                                                                                                 << sqrt(hcov_tot->GetBinContent(i+1,i+1))/hpred->GetBinContent(i+1) << " "
                                                                                                                 << sqrt(hcov_stat->GetBinContent(i+1,i+1))/hpred->GetBinContent(i+1) << " "
                                                                                                                 << sqrt(hcov_mcstat->GetBinContent(i+1,i+1))/hpred->GetBinContent(i+1) << " "
                                                                                                                 << sqrt(hcov_add->GetBinContent(i+1,i+1))/hpred->GetBinContent(i+1) << " "
                                                                                                                 << sqrt(hcov_flux->GetBinContent(i+1,i+1))/hpred->GetBinContent(i+1) << " "
                                                                                                                 << sqrt(hcov_det->GetBinContent(i+1,i+1))/hpred->GetBinContent(i+1) << " "
                                                                                                                 << sqrt(hcov_xs->GetBinContent(i+1,i+1))/hpred->GetBinContent(i+1) << " "
														<< sqrt(hcov_rw->GetBinContent(i+1,i+1))/hpred->GetBinContent(i+1) << " "
														<< sqrt(hcov_rw_cor->GetBinContent(i+1,i+1))/hpred->GetBinContent(i+1) << " "
                                                                                                                 << std::endl;
      myfile << i << " " << hdata_obsch_1->GetBinContent(i+1) << " " << histo_3->GetBinContent(i+1) << " " << histo_5->GetBinContent(i+1) << " " << hmeas->GetBinContent(i+1) << " " << hpred->GetBinContent(i+1) << " "
                                                                                                                 << sqrt(hcov_tot->GetBinContent(i+1,i+1))/hpred->GetBinContent(i+1) << " "
                                                                                                                 << sqrt(hcov_stat->GetBinContent(i+1,i+1))/hpred->GetBinContent(i+1) << " "
                                                                                                                 << sqrt(hcov_mcstat->GetBinContent(i+1,i+1))/hpred->GetBinContent(i+1) << " "
                                                                                                                 << sqrt(hcov_add->GetBinContent(i+1,i+1))/hpred->GetBinContent(i+1) << " "
                                                                                                                 << sqrt(hcov_flux->GetBinContent(i+1,i+1))/hpred->GetBinContent(i+1) << " "
                                                                                                                 << sqrt(hcov_det->GetBinContent(i+1,i+1))/hpred->GetBinContent(i+1) << " "
                                                                                                                 << sqrt(hcov_xs->GetBinContent(i+1,i+1))/hpred->GetBinContent(i+1) << " "
														<< sqrt(hcov_rw->GetBinContent(i+1,i+1))/hpred->GetBinContent(i+1) << " "
														<< sqrt(hcov_rw_cor->GetBinContent(i+1,i+1))/hpred->GetBinContent(i+1) << " "
                                                                                                                 << std::endl; 
    }
    if(i>=reco_bin_2 && i<=reco_bin_4){
      int j = i-reco_bin_1;
      std::cout << i << " " << hdata_obsch_2->GetBinContent(j) << " " << histo_4->GetBinContent(j) << " " << histo_6->GetBinContent(j) << " " << hmeas->GetBinContent(i+1) << " " << hpred->GetBinContent(i+1) << " "
                                                                                                                 << sqrt(hcov_tot->GetBinContent(i+1,i+1))/hpred->GetBinContent(i+1) << " "
                                                                                                                 << sqrt(hcov_stat->GetBinContent(i+1,i+1))/hpred->GetBinContent(i+1) << " "
                                                                                                                 << sqrt(hcov_mcstat->GetBinContent(i+1,i+1))/hpred->GetBinContent(i+1) << " "
                                                                                                                 << sqrt(hcov_add->GetBinContent(i+1,i+1))/hpred->GetBinContent(i+1) << " "
                                                                                                                 << sqrt(hcov_flux->GetBinContent(i+1,i+1))/hpred->GetBinContent(i+1) << " "
                                                                                                                 << sqrt(hcov_det->GetBinContent(i+1,i+1))/hpred->GetBinContent(i+1) << " "
                                                                                                                 << sqrt(hcov_xs->GetBinContent(i+1,i+1))/hpred->GetBinContent(i+1) << " "
														<< sqrt(hcov_rw->GetBinContent(i+1,i+1))/hpred->GetBinContent(i+1) << " "
														<< sqrt(hcov_rw_cor->GetBinContent(i+1,i+1))/hpred->GetBinContent(i+1) << " "

                                                                                                                 << std::endl;
      myfile << i << " " << hdata_obsch_2->GetBinContent(j) << " " << histo_4->GetBinContent(j) << " " << histo_6->GetBinContent(j) << " " << hmeas->GetBinContent(i+1) << " " << hpred->GetBinContent(i+1) << " "
                                                                                                                 << sqrt(hcov_tot->GetBinContent(i+1,i+1))/hpred->GetBinContent(i+1) << " "
                                                                                                                 << sqrt(hcov_stat->GetBinContent(i+1,i+1))/hpred->GetBinContent(i+1) << " "
                                                                                                                 << sqrt(hcov_mcstat->GetBinContent(i+1,i+1))/hpred->GetBinContent(i+1) << " "
                                                                                                                 << sqrt(hcov_add->GetBinContent(i+1,i+1))/hpred->GetBinContent(i+1) << " "
                                                                                                                 << sqrt(hcov_flux->GetBinContent(i+1,i+1))/hpred->GetBinContent(i+1) << " "
                                                                                                                 << sqrt(hcov_det->GetBinContent(i+1,i+1))/hpred->GetBinContent(i+1) << " "
                                                                                                                 << sqrt(hcov_xs->GetBinContent(i+1,i+1))/hpred->GetBinContent(i+1) << " "
														<< sqrt(hcov_rw->GetBinContent(i+1,i+1))/hpred->GetBinContent(i+1) << " "
														<< sqrt(hcov_rw_cor->GetBinContent(i+1,i+1))/hpred->GetBinContent(i+1) << " "

                                                                                                                 << std::endl;
     
    }
  }
  std::cout << "----------------------------------------------------------------------------" << std::endl;
  myfile.close();


  // stack plots ...
  TH1F *h10 = new TH1F("h10","h10",nbin_meas,0.5,nbin_meas+0.5); // stat
  TH1F *h20 = new TH1F("h20","h10",nbin_meas,0.5,nbin_meas+0.5); // mcstat
  TH1F *h30 = new TH1F("h30","h10",nbin_meas,0.5,nbin_meas+0.5); // add
  TH1F *h40 = new TH1F("h40","h10",nbin_meas,0.5,nbin_meas+0.5); // det
  TH1F *h50 = new TH1F("h50","h10",nbin_meas,0.5,nbin_meas+0.5); // flux
  TH1F *h60 = new TH1F("h60","h10",nbin_meas,0.5,nbin_meas+0.5); // xs
  TH1F *h61 = new TH1F("h61","h10",nbin_meas,0.5,nbin_meas+0.5); // rw
  TH1F *h62 = new TH1F("h62","h10",nbin_meas,0.5,nbin_meas+0.5); // rw_cor
  TH1F *h70 = new TH1F("h70","h10",nbin_meas,0.5,nbin_meas+0.5); // total
  
  for (Int_t i=0;i!=nbin_meas;i++){
    // std::cout << i << " " << sqrt(hcov_tot->GetBinContent(i+1,i+1))/hpred->GetBinContent(i+1) << " "
    // 	      << sqrt(hcov_stat->GetBinContent(i+1,i+1))/hpred->GetBinContent(i+1) << " "
    // 	      << sqrt(hcov_mcstat->GetBinContent(i+1,i+1))/hpred->GetBinContent(i+1) << " "
    // 	      << sqrt(hcov_add->GetBinContent(i+1,i+1))/hpred->GetBinContent(i+1) << " "
    // 	      << sqrt(hcov_flux->GetBinContent(i+1,i+1))/hpred->GetBinContent(i+1) << " "
    // 	      << sqrt(hcov_det->GetBinContent(i+1,i+1))/hpred->GetBinContent(i+1) << " "
    // 	      << sqrt(hcov_xs->GetBinContent(i+1,i+1))/hpred->GetBinContent(i+1) << " "
    // 	      << std::endl;

    if (hpred->GetBinContent(i+1)==0) continue;

    h10->SetBinContent(i+1,sqrt(hcov_stat->GetBinContent(i+1,i+1))/hpred->GetBinContent(i+1));
    h20->SetBinContent(i+1,sqrt(hcov_mcstat->GetBinContent(i+1,i+1))/hpred->GetBinContent(i+1));
    h30->SetBinContent(i+1,sqrt(hcov_add->GetBinContent(i+1,i+1))/hpred->GetBinContent(i+1));
    h40->SetBinContent(i+1,sqrt(hcov_flux->GetBinContent(i+1,i+1))/hpred->GetBinContent(i+1));
    h50->SetBinContent(i+1,sqrt(hcov_det->GetBinContent(i+1,i+1))/hpred->GetBinContent(i+1));
    h60->SetBinContent(i+1,sqrt(hcov_xs->GetBinContent(i+1,i+1))/hpred->GetBinContent(i+1));
    h61->SetBinContent(i+1,sqrt(hcov_rw->GetBinContent(i+1,i+1))/hpred->GetBinContent(i+1));
    h62->SetBinContent(i+1,sqrt(hcov_rw_cor->GetBinContent(i+1,i+1))/hpred->GetBinContent(i+1));
    h70->SetBinContent(i+1,sqrt(hcov_tot->GetBinContent(i+1,i+1))/hpred->GetBinContent(i+1));
  }

  h70->Draw();
  h70->SetXTitle("Bin no");
  h70->SetTitle("Relative Uncertainties");
  h70->SetLineColor(1);
  h70->SetLineWidth(2);
  
  h10->Draw("same");  //stat
  h10->SetLineColor(9);
  h20->Draw("same");  // mcstat
  h20->SetLineColor(8);
  h30->Draw("same"); //dirt
  h30->SetLineColor(3);
  h40->Draw("same"); //flux
  h40->SetLineColor(2);
  h50->Draw("same"); //det
  h50->SetLineColor(6);
  
  h60->Draw("same"); // xs
  h60->SetLineColor(4);

  h61->Draw("same"); // rw
  h61->SetLineColor(kYellow+1);

  h62->Draw("same"); // rw_cor
  h62->SetLineColor(kYellow-7);

  h10->SetLineWidth(2);
  h20->SetLineWidth(2);
  h30->SetLineWidth(2);
  h40->SetLineWidth(2);
  h50->SetLineWidth(2);
  h60->SetLineWidth(2);
  h61->SetLineWidth(2);
  h62->SetLineWidth(2);

  h70->GetYaxis()->SetRangeUser(0,1.);

  TLegend *le1 = new TLegend(0.6,0.6,0.89,0.89);
  le1->AddEntry(h70,"Total","l");
  le1->AddEntry(h10,"Stat.","l");
  le1->AddEntry(h20,"MC stat.","l");
  le1->AddEntry(h30,"Dirt","l");
  le1->AddEntry(h40,"Flux","l");
  le1->AddEntry(h50,"Det.","l");
  le1->AddEntry(h60,"Xs","l");
  le1->AddEntry(h61,"RW","l");
  le1->AddEntry(h62,"RW cor.","l");
  le1->Draw();
  
  
  

 TFile *file = new TFile("wiener.root","RECREATE");
 htrue_signal->SetDirectory(file);
 hmeas->SetDirectory(file);
 hR->SetDirectory(file);
 hcov_stat->SetDirectory(file);
 hcov_mcstat->SetDirectory(file);
 hcov_add->SetDirectory(file);
 hcov_det->SetDirectory(file);
 hcov_flux->SetDirectory(file);
 hcov_xs->SetDirectory(file);
 hcov_rw->SetDirectory(file);
 hcov_rw_cor->SetDirectory(file);
 hcov_tot->SetDirectory(file);
 file->Write();
 file->Close();
  
}
