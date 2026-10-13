void convert_wiener(int n_diff_bins_2D=9, int n_diff_bins_3D=4){

  bool just_stat_uncertainty      = false;
  bool use_fakedata               = false;
  bool new_noise                  = true;
  int add_noise                   = 0;

  bool draw = false;
  bool data_wgu = false;
  bool systematics_wgu = false;

  gStyle->SetOptStat(0);
  
  std::string filename1     =     "merge_xs.root";
  if (data_wgu) { filename1 = "wgu/merge_xs.root"; }
  TFile *file1 = new TFile(filename1.c_str());
  TMatrixD *mat_collapse = (TMatrixD*)file1->Get("mat_collapse"); // collapse				// 108 x 36  (36 = 18x2)
  TMatrixD *cov_mat_add  = (TMatrixD*)file1->Get("cov_mat_add"); // additional uncertainties
  TVectorD *vec_signal   = (TVectorD*)file1->Get("vec_signal"); // xs in pbar
  TMatrixD *mat_R        = (TMatrixD*)file1->Get("mat_R"); // smearing matrix				// 108 x 11  (108 = 36x3 = 18x2x3)

  //std::cout << "mat_R      row,col    = " << mat_R->GetNrows() << ", " << mat_R->GetNcols() << std::endl;				//(6 x n_E_reco x n_theta_reco) x (n_E_true x n_theta_true)
  //std::cout << "collapse   row,col    = " << mat_collapse->GetNrows() << ", " << mat_collapse->GetNcols() << std::endl;		//(6 x n_E_reco x n_theta_reco) x (2 x n_E_reco x n_theta_reco)
  //std::cout << "additional row,col    = " << cov_mat_add->GetNrows() << ", " << cov_mat_add->GetNcols() << std::endl;			//(6 x n_E_reco x n_theta_reco) x (6 x n_E_reco x n_theta_reco)
  //std::cout << "signal vec row        = " << vec_signal->GetNrows() << std::endl;							//n_E_true x n_theta_true

  int n_Ebins_reco = mat_R->GetNrows()/n_diff_bins_2D/n_diff_bins_3D/3/2;	// N_ebin * N_theta * (sig,bkg,ext) * (FC,PC) * (Pmu,Ehad)
  int nbins_histo = n_Ebins_reco*n_diff_bins_2D*n_diff_bins_3D;			// 16*9*4
std::cout << "nbins_histo = " << nbins_histo << std::endl;
  TH1F *hdata_obsch_1 = new TH1F("hdata_obsch_1_merge","hdata_obsch_1_merge", nbins_histo, 0, nbins_histo);
  TH1F *hdata_obsch_2 = new TH1F("hdata_obsch_2_merge","hdata_obsch_2_merge", nbins_histo, 0, nbins_histo);
  TH1F *hmc_obsch_1   = new TH1F("hmc_obsch_1_merge",  "hmc_obsch_1_merge",   nbins_histo, 0, nbins_histo);
  TH1F *hmc_obsch_2   = new TH1F("hmc_obsch_2_merge",  "hmc_obsch_2_merge",   nbins_histo, 0, nbins_histo);
  TH1F *histo_1       = new TH1F("histo_1_merge",      "histo_1_merge",       nbins_histo, 0, nbins_histo);
  TH1F *histo_2       = new TH1F("histo_2_merge",      "histo_2_merge",       nbins_histo, 0, nbins_histo);
  TH1F *histo_3       = new TH1F("histo_3_merge",      "histo_3_merge",       nbins_histo, 0, nbins_histo);
  TH1F *histo_4       = new TH1F("histo_4_merge",      "histo_4_merge",       nbins_histo, 0, nbins_histo);
  TH1F *histo_5       = new TH1F("histo_5_merge",      "histo_5_merge",       nbins_histo, 0, nbins_histo);
  TH1F *histo_6       = new TH1F("histo_6_merge",      "histo_6_merge",       nbins_histo, 0, nbins_histo);

  for (int k=0;k<n_diff_bins_3D;k++) {
    for (int i=0;i<n_diff_bins_2D;i++) {
      int FC_index = 1+i+k*n_diff_bins_2D;
      int PC_index = 1+i+k*n_diff_bins_2D+n_diff_bins_2D*n_diff_bins_3D;
      std::string FC_data_str = "hdata_obsch_" + std::to_string(FC_index);					// data numuCC FC
      std::string PC_data_str = "hdata_obsch_" + std::to_string(PC_index);					// data numuCC PC
std::cout << FC_data_str << ",  " << PC_data_str << std::endl;
      if (use_fakedata) {
        FC_data_str = "hmc_obsch_" + std::to_string(FC_index);							// fakedata numuCC FC
        PC_data_str = "hmc_obsch_" + std::to_string(PC_index);							// fakedata numuCC PC
      }
      std::string FC_mc_str   = "hmc_obsch_"   + std::to_string(FC_index);					// pred numuCC FC
      std::string PC_mc_str   = "hmc_obsch_"   + std::to_string(PC_index);					// pred numuCC PC
      std::string histo_str_1 = "histo_"       + std::to_string(FC_index + 0*2*n_diff_bins_2D*n_diff_bins_3D);	// sig numuCC FC
      std::string histo_str_2 = "histo_"       + std::to_string(PC_index + 0*2*n_diff_bins_2D*n_diff_bins_3D);	// sig numuCC PC
      std::string histo_str_3 = "histo_"       + std::to_string(FC_index + 1*2*n_diff_bins_2D*n_diff_bins_3D);	// bkg numuCC FC
      std::string histo_str_4 = "histo_"       + std::to_string(PC_index + 1*2*n_diff_bins_2D*n_diff_bins_3D);	// bkg numuCC PC
      std::string histo_str_5 = "histo_"       + std::to_string(FC_index + 2*2*n_diff_bins_2D*n_diff_bins_3D);	// EXT FC
      std::string histo_str_6 = "histo_"       + std::to_string(PC_index + 2*2*n_diff_bins_2D*n_diff_bins_3D);	// EXT PC
      TH1F *hdata_obsch_1_tmp = (TH1F*)file1->Get(FC_data_str.c_str());						// data numuCC FC
      TH1F *hdata_obsch_2_tmp = (TH1F*)file1->Get(PC_data_str.c_str());						// data numuCC PC
      TH1F *hmc_obsch_1_tmp   = (TH1F*)file1->Get(FC_mc_str.c_str());						// pred numuCC FC
      TH1F *hmc_obsch_2_tmp   = (TH1F*)file1->Get(PC_mc_str.c_str());						// pred numuCC PC
      TH1F *histo_1_tmp       = (TH1F*)file1->Get(histo_str_1.c_str());						// sig numuCC FC
      TH1F *histo_2_tmp       = (TH1F*)file1->Get(histo_str_2.c_str());						// sig numuCC PC
      TH1F *histo_3_tmp       = (TH1F*)file1->Get(histo_str_3.c_str());						// bkg numuCC FC
      TH1F *histo_4_tmp       = (TH1F*)file1->Get(histo_str_4.c_str());						// bkg numuCC PC
      TH1F *histo_5_tmp       = (TH1F*)file1->Get(histo_str_5.c_str());						// EXT FC
      TH1F *histo_6_tmp       = (TH1F*)file1->Get(histo_str_6.c_str());						// EXT PC

      for (int j=0;j<n_Ebins_reco;j++) {
        int index = (k*n_diff_bins_2D+i)*n_Ebins_reco+j+1;
        histo_1->SetBinContent(      index, histo_1_tmp->GetBinContent(j+1));
        histo_2->SetBinContent(      index, histo_2_tmp->GetBinContent(j+1));
        histo_3->SetBinContent(      index, histo_3_tmp->GetBinContent(j+1));
        histo_4->SetBinContent(      index, histo_4_tmp->GetBinContent(j+1));
        histo_5->SetBinContent(      index, histo_5_tmp->GetBinContent(j+1));
        histo_6->SetBinContent(      index, histo_6_tmp->GetBinContent(j+1));
        hdata_obsch_1->SetBinContent(index, hdata_obsch_1_tmp->GetBinContent(j+1));
        hdata_obsch_2->SetBinContent(index, hdata_obsch_2_tmp->GetBinContent(j+1));
std::cout << "index, j+1 = " << index<< ",  " << (j+1) << std::endl;
        hmc_obsch_1->SetBinContent(  index, hmc_obsch_1_tmp->GetBinContent(j+1));
        hmc_obsch_2->SetBinContent(  index, hmc_obsch_2_tmp->GetBinContent(j+1));
      }
    }
  }

  int nbin_true = vec_signal->GetNoElements();
  int nbin_meas = hdata_obsch_1->GetNbinsX() + hdata_obsch_2->GetNbinsX();

  std::cout << "nbin_meas = " << nbin_meas << std::endl;
  // response matrix
  
  TMatrixD mat_collapse_T(mat_collapse->GetNcols(), mat_collapse->GetNrows()); 
  mat_collapse_T.Transpose(*mat_collapse);
  TMatrixD mat_R_collapse = (mat_collapse_T)*(*mat_R);
  
  //mat_R_collapse.Draw("COLZ");

  TH2D* hR = new TH2D("hR","hR",nbin_meas,0.5,nbin_meas+0.5, nbin_true, 0.5, nbin_true+0.5);
  for (Int_t i=0;i!=nbin_true;i++){
    for (Int_t j=0;j!=nbin_meas;j++){
      hR->SetBinContent(j+1,i+1,mat_R_collapse(j,i));
    }
  }
  //  hR->Draw("COLZ");
  
  // covariance matrix
  TH2D *hcov_tot    = new TH2D("hcov_tot",   "hcov_tot",   nbin_meas,0.5,nbin_meas+0.5,nbin_meas,0.5,nbin_meas+0.5);
  // MC statistical uncertainties ...
  TH2D *hcov_mcstat = new TH2D("hcov_mcstat","hcov_mcstat",nbin_meas,0.5,nbin_meas+0.5,nbin_meas,0.5,nbin_meas+0.5);
  ifstream infile("mcstat/xs_tot.log");
  double temp, temp1,err2;
  infile >> temp >> temp;
  for (Int_t i=0;i!=nbin_meas;i++){
    infile >> temp >> temp >> temp >> err2 >> temp;
    hcov_mcstat->SetBinContent(i+1,i+1,err2);
  }
  //hcov_mcstat->Draw("COLZ");

  // statistical uncertainties
  TH2D *hcov_stat = new TH2D("hcov_stat","hcov_stat",nbin_meas,0.5,nbin_meas+0.5,nbin_meas,0.5,nbin_meas+0.5);
  for (int k=0;k<n_diff_bins_3D;k++) {
    for (int i=0;i<n_diff_bins_2D;i++) {
      for (int j=0;j<n_Ebins_reco;j++) {
        int FC_index = (i+k*n_diff_bins_2D                              ) * n_Ebins_reco + j+1;	//FC
        int PC_index = (i+k*n_diff_bins_2D+n_diff_bins_2D*n_diff_bins_3D) * n_Ebins_reco + j+1;	//PC

        double meas = hdata_obsch_1->GetBinContent(FC_index);
        double pred = hmc_obsch_1->GetBinContent(FC_index);
        double content;
        if (pred !=0){
          if (meas == 0){
            content = pred/2.;
          }else{
	    content = 3./(1./meas+2./pred);
          }
        }else{
          content = 0;
        }
        hcov_stat->SetBinContent(FC_index,FC_index,content);
        //std::cout << "i,j,FC,PC, FC_content = " << i << ",  " << j << ",  " << FC_index << ",     " << PC_index << ",     " << content << std::endl;

        meas = hdata_obsch_2->GetBinContent(FC_index);
        pred = hmc_obsch_2->GetBinContent(FC_index);
        if (pred !=0){
          if (meas == 0){
	    content = pred/2.;
          }else{
	    content = 3./(1./meas+2./pred);
          }
        }else{
          content = 0;
        }
        hcov_stat->SetBinContent(PC_index,PC_index,content);
      }
    }
  }

  // additional uncertainty
  TH2D *hcov_add = new TH2D("hcov_add","hcov_add",nbin_meas,0.5,nbin_meas+0.5,nbin_meas,0.5,nbin_meas+0.5);
  TMatrixD mat_add = mat_collapse_T * (*cov_mat_add) * (*mat_collapse);
  //  mat_add.Draw("COLZ");
  for (Int_t i=0;i!=nbin_meas;i++){
    hcov_add->SetBinContent(i+1,i+1,mat_add(i,i));
  }
  //hcov_add->Draw("COLZ");
  
  // detector systematics
  TH2D *hcov_det = new TH2D("hcov_det","hcov_det",nbin_meas,0.5,nbin_meas+0.5,nbin_meas,0.5,nbin_meas+0.5);
  int nold = mat_collapse->GetNrows();
  TVectorD vec_nominal(nold);

  for (int k=0;k<n_diff_bins_3D;k++) {
    for (int i=0;i<n_diff_bins_2D;i++) {
      for (int j=0;j<n_Ebins_reco;j++) {

        int index0 = (k*n_diff_bins_2D+i)*n_Ebins_reco + j;
        int index1 = index0 + n_Ebins_reco*n_diff_bins_2D*n_diff_bins_3D*0;
        int index2 = index0 + n_Ebins_reco*n_diff_bins_2D*n_diff_bins_3D*1;
        int index3 = index0 + n_Ebins_reco*n_diff_bins_2D*n_diff_bins_3D*2;
        int index4 = index0 + n_Ebins_reco*n_diff_bins_2D*n_diff_bins_3D*3;
        int index5 = index0 + n_Ebins_reco*n_diff_bins_2D*n_diff_bins_3D*4;
        int index6 = index0 + n_Ebins_reco*n_diff_bins_2D*n_diff_bins_3D*5;
        vec_nominal(index1) = histo_1->GetBinContent(index0+1);
        vec_nominal(index2) = histo_2->GetBinContent(index0+1);
        vec_nominal(index3) = histo_3->GetBinContent(index0+1);
        vec_nominal(index4) = histo_4->GetBinContent(index0+1);
        vec_nominal(index5) = histo_5->GetBinContent(index0+1);
        vec_nominal(index6) = histo_6->GetBinContent(index0+1);
      }
    }
  }
  //vec_nominal.Draw();

  std::map<int, TString> map_det_str;
  std::string detvar_prefix = "./DetVar/";
  if (systematics_wgu) { detvar_prefix = "./wgu/DetVar/"; }
  map_det_str[1]  = detvar_prefix+"cov_LYDown.root";
  map_det_str[2]  = detvar_prefix+"cov_LYRayleigh.root";
  map_det_str[3]  = detvar_prefix+"cov_Recomb2.root";
  map_det_str[4]  = detvar_prefix+"cov_SCE.root";
  // map_det_str[5]  = detvar_prefix+"cov_WMdEdx.root";
  map_det_str[6]  = detvar_prefix+"cov_WMThetaXZ.root";
  map_det_str[7]  = detvar_prefix+"cov_WMThetaYZ.root";
  map_det_str[8]  = detvar_prefix+"cov_WMX.root";
  map_det_str[9]  = detvar_prefix+"cov_WMYZ.root";
  map_det_str[10] = detvar_prefix+"cov_LYatt.root";
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
  TMatrixD mat_det = mat_collapse_T * frac_det * (*mat_collapse);
  for (Int_t i=0;i!=nbin_meas;i++){
    for (Int_t j=0;j!=nbin_meas;j++){
      hcov_det->SetBinContent(i+1,j+1,mat_det(i,j));
    }
  }
  
  // Flux systematics
  TH2D *hcov_flux = new TH2D("hcov_flux","hcov_flux",nbin_meas,0.5,nbin_meas+0.5,nbin_meas,0.5,nbin_meas+0.5);
  TMatrixD frac_flux(nold,nold);
  std::string flux_str            =     "./XsFlux/cov_%d.root";
  if (systematics_wgu) { flux_str = "./wgu/XsFlux/cov_%d.root"; }
  for(Int_t i=1;i!=17;i++){
    TFile tmp_file(Form(flux_str.c_str(),i));
    TMatrixD *tmp_matrix = (TMatrixD*)tmp_file.Get(Form("frac_cov_xf_mat_%d",i));
    frac_flux += *tmp_matrix;
  }
  for (Int_t i=0;i!=nold;i++){
    for (Int_t j=0;j!=nold;j++){
      frac_flux(i,j) *= vec_nominal(i) * vec_nominal(j);
    }
  }
  //frac_flux.Draw("COLZ");
  TMatrixD mat_flux = mat_collapse_T * frac_flux * (*mat_collapse);
  //mat_flux.Draw("COLZ");
  for (Int_t i=0;i!=nbin_meas;i++){
    for (Int_t j=0;j!=nbin_meas;j++){
      hcov_flux->SetBinContent(i+1,j+1,mat_flux(i,j));
    }
  }
  // hcov_flux->Draw("COLZ");

  // Xs systematics
  TH2D *hcov_xs      = new TH2D("hcov_xs","hcov_xs",nbin_meas,0.5,nbin_meas+0.5,nbin_meas,0.5,nbin_meas+0.5);

  {
    TFile tmp_file("./XsFlux/cov_xs.root");
    //TFile tmp_file("./XsFlux/cov_17.root");
    TMatrixD *frac_xs = (TMatrixD*)tmp_file.Get("frac_cov_xf_mat_17");
    
    for (Int_t i=0;i!=nold;i++){
      for (Int_t j=0;j!=nold;j++){
	(*frac_xs)(i,j) *= vec_nominal(i) * vec_nominal(j);
      }
    }
    //frac_flux.Draw("COLZ");
    TMatrixD mat_xs = mat_collapse_T * (*frac_xs) * (*mat_collapse);

    for (Int_t i=0;i!=nbin_meas;i++){
      for (Int_t j=0;j!=nbin_meas;j++){
	hcov_xs->SetBinContent(i+1,j+1,mat_xs(i,j));
      }
    }
    //    hcov_xs->Draw("COLZ");
  }
  //vec_nominal.Draw();


  //before storing the measured and predicted values, apply a gaussian noise to them for debugging purposes
  TH1D *hnoise           = new TH1D("hnoise",          "hnoise"          ,nbin_meas,                              0.5,  nbin_meas                               +0.5);
  TH1D *hnoise_collapsed = new TH1D("hnoise_collapsed","hnoise_collapsed",nbin_meas/n_diff_bins_2D/n_diff_bins_3D,0.5, (nbin_meas/n_diff_bins_2D/n_diff_bins_3D)+0.5);

  if (new_noise) {
    int seed = std::random_device{}();
    std::cout << "seed = " << seed << std::endl;
    std::mt19937 generator(seed);	// https://en.cppreference.com/w/cpp/numeric/random/mersenne_twister_engine
    for (int k=0;k<n_diff_bins_3D;k++) {
      for (int i=0;i<n_diff_bins_2D;i++) {
        for (int j=0;j<2*n_Ebins_reco;j++) {
          int index = 2*(k*n_diff_bins_2D+i)*n_Ebins_reco + j+1;

          double sigma = std::sqrt(hcov_stat->GetBinContent(index,index));
          std::normal_distribution<double> dist(0, sigma);
          double noise = dist(generator);

          hnoise->SetBinContent(index, noise);
          hnoise_collapsed->SetBinContent(j+1, hnoise_collapsed->GetBinContent(j+1)+noise);
          //std::cout << "i,j,index, sigma, noise, hnoise_collapsed = " << i << ",  " << j << ",   " << index << ",     " << sigma << ",     " << noise << ",     " << hnoise_collapsed->GetBinContent(j+1) << std::endl;
        }
      }
    }
  } else {
    //read in niose histograms
    TFile *noise_file = new TFile("/uboone/data/users/lcoopert/LEE/LEEana_xs_Emuon/wiener_svd/noise_file.root");
    hnoise           = (TH1D*)noise_file->Get("hnoise"); 
    hnoise_collapsed = (TH1D*)noise_file->Get("hnoise_collapsed"); 
  }

  // true signal ...
  TH1D *htrue_signal = new TH1D("htrue_signal","htrue_signal",nbin_true, 0.5, nbin_true+0.5);
  for(int i=0;i!=nbin_true;i++){
    htrue_signal->SetBinContent(i+1,(*vec_signal)(i));
  }
  // measurement ...
  TH1D *hmeas = new TH1D("hmeas","hmeas",nbin_meas,0.5, nbin_meas+0.5);
  TH1D *hpred = new TH1D("hpred","hpred",nbin_meas,0.5, nbin_meas+0.5);
  for (int k=0;k<n_diff_bins_3D;k++) {
    for (int i=0;i<n_diff_bins_2D;i++) {
      for (int j=0;j<n_Ebins_reco;j++) {
        int FC_index = (i+k*n_diff_bins_2D                              ) * n_Ebins_reco + j+1;	//FC
        int PC_index = (i+k*n_diff_bins_2D+n_diff_bins_2D*n_diff_bins_3D) * n_Ebins_reco + j+1;	//PC
        hmeas->SetBinContent(FC_index, hdata_obsch_1->GetBinContent(FC_index) - histo_3->GetBinContent(FC_index) - histo_5->GetBinContent(FC_index) + add_noise*hnoise->GetBinContent(FC_index));
        hmeas->SetBinContent(PC_index, hdata_obsch_2->GetBinContent(FC_index) - histo_4->GetBinContent(FC_index) - histo_6->GetBinContent(FC_index) + add_noise*hnoise->GetBinContent(PC_index));

        hpred->SetBinContent(FC_index, hmc_obsch_1->GetBinContent(FC_index) );
        hpred->SetBinContent(PC_index, hmc_obsch_2->GetBinContent(FC_index) );
      }
    }
  }


  hcov_tot->Add(hcov_stat);
  if (!just_stat_uncertainty) {
    hcov_tot->Add(hcov_mcstat);
    hcov_tot->Add(hcov_add);
    hcov_tot->Add(hcov_flux);
    hcov_tot->Add(hcov_det);
    hcov_tot->Add(hcov_xs);
  }

  for (Int_t i=0;i!=nbin_meas;i++){
    std::cout << i << " " 
              << sqrt(hcov_tot->GetBinContent(i+1,i+1))/hpred->GetBinContent(i+1) << " "
	      << sqrt(hcov_stat->GetBinContent(i+1,i+1))/hpred->GetBinContent(i+1) << " "
	      << sqrt(hcov_mcstat->GetBinContent(i+1,i+1))/hpred->GetBinContent(i+1) << " "
	      << sqrt(hcov_add->GetBinContent(i+1,i+1))/hpred->GetBinContent(i+1) << " "
	      << sqrt(hcov_flux->GetBinContent(i+1,i+1))/hpred->GetBinContent(i+1) << " "
	      << sqrt(hcov_det->GetBinContent(i+1,i+1))/hpred->GetBinContent(i+1) << " "
	      << sqrt(hcov_xs->GetBinContent(i+1,i+1))/hpred->GetBinContent(i+1) << " " 
	      << std::endl;
  }



  // stack plots ...

  TH1F *h10 = new TH1F("h10","h10",nbin_meas,0.5,nbin_meas+0.5); // stat
  TH1F *h20 = new TH1F("h20","h10",nbin_meas,0.5,nbin_meas+0.5); // mcstat
  TH1F *h30 = new TH1F("h30","h10",nbin_meas,0.5,nbin_meas+0.5); // add
  TH1F *h40 = new TH1F("h40","h10",nbin_meas,0.5,nbin_meas+0.5); // det
  TH1F *h50 = new TH1F("h50","h10",nbin_meas,0.5,nbin_meas+0.5); // flux
  TH1F *h60 = new TH1F("h60","h10",nbin_meas,0.5,nbin_meas+0.5); // xs
  TH1F *h70 = new TH1F("h70","h10",nbin_meas,0.5,nbin_meas+0.5); // total

  TH1F *h10_rel = new TH1F("h10_rel","h10_rel",nbin_meas,0.5,nbin_meas+0.5); // stat
  TH1F *h20_rel = new TH1F("h20_rel","h10_rel",nbin_meas,0.5,nbin_meas+0.5); // mcstat
  TH1F *h30_rel = new TH1F("h30_rel","h10_rel",nbin_meas,0.5,nbin_meas+0.5); // add
  TH1F *h40_rel = new TH1F("h40_rel","h10_rel",nbin_meas,0.5,nbin_meas+0.5); // det
  TH1F *h50_rel = new TH1F("h50_rel","h10_rel",nbin_meas,0.5,nbin_meas+0.5); // flux
  TH1F *h60_rel = new TH1F("h60_rel","h10_rel",nbin_meas,0.5,nbin_meas+0.5); // xs
  TH1F *h70_rel = new TH1F("h70_rel","h10_rel",nbin_meas,0.5,nbin_meas+0.5); // total



  for (Int_t i=0;i!=nbin_meas;i++){
    if (hpred->GetBinContent(i+1)!=0){
      h10->SetBinContent(i+1,sqrt(hcov_stat->GetBinContent(i+1,i+1))/hpred->GetBinContent(i+1));
      h20->SetBinContent(i+1,sqrt(hcov_mcstat->GetBinContent(i+1,i+1))/hpred->GetBinContent(i+1));
      h30->SetBinContent(i+1,sqrt(hcov_add->GetBinContent(i+1,i+1))/hpred->GetBinContent(i+1));
      h40->SetBinContent(i+1,sqrt(hcov_flux->GetBinContent(i+1,i+1))/hpred->GetBinContent(i+1));
      h50->SetBinContent(i+1,sqrt(hcov_det->GetBinContent(i+1,i+1))/hpred->GetBinContent(i+1));
      h60->SetBinContent(i+1,sqrt(hcov_xs->GetBinContent(i+1,i+1))/hpred->GetBinContent(i+1));
      h70->SetBinContent(i+1,sqrt(hcov_tot->GetBinContent(i+1,i+1))/hpred->GetBinContent(i+1));

      h10_rel->SetBinContent(i+1, hcov_stat->GetBinContent(i+1,i+1)   / hcov_tot->GetBinContent(i+1,i+1));
      h20_rel->SetBinContent(i+1, hcov_mcstat->GetBinContent(i+1,i+1) / hcov_tot->GetBinContent(i+1,i+1));
      h30_rel->SetBinContent(i+1, hcov_add->GetBinContent(i+1,i+1)    / hcov_tot->GetBinContent(i+1,i+1));
      h40_rel->SetBinContent(i+1, hcov_flux->GetBinContent(i+1,i+1)   / hcov_tot->GetBinContent(i+1,i+1));
      h50_rel->SetBinContent(i+1, hcov_det->GetBinContent(i+1,i+1)    / hcov_tot->GetBinContent(i+1,i+1));
      h60_rel->SetBinContent(i+1, hcov_xs->GetBinContent(i+1,i+1)     / hcov_tot->GetBinContent(i+1,i+1));
      h70_rel->SetBinContent(i+1, hcov_tot->GetBinContent(i+1,i+1)    / hcov_tot->GetBinContent(i+1,i+1));
    }
  }

  //DetVar covariance and correlation matrices
  int nbins_cov = hcov_det->GetNbinsX();
  TH2D* detvar_covariance_matrix  = new TH2D("detvar_covariance_matrix", "Detector Systematic Covariance Matrix", nbins_cov,0.5,0.5+nbins_cov,nbins_cov,0.5,0.5+nbins_cov);
  TH2D* detvar_correlation_matrix = new TH2D("detvar_correlation_matrix","Detector Systematic Correlation Matrix",nbins_cov,0.5,0.5+nbins_cov,nbins_cov,0.5,0.5+nbins_cov);
  for (int i=0;i<nbins_cov;i++) {
    for (int j=0;j<nbins_cov;j++) {
      double cov_ij = sqrt(hcov_det->GetBinContent(i+1,j+1));
      double cov_ii = sqrt(hcov_det->GetBinContent(i+1,i+1));
      double cov_jj = sqrt(hcov_det->GetBinContent(j+1,j+1));
      double val = cov_ij/std::sqrt(cov_ii*cov_jj);      
      if (cov_ij==0) { val = 0; }
      if (i==j)      { val = 1; }
      if (cov_ij > 4000) { cov_ij = 4000; }
      detvar_covariance_matrix->SetBinContent(i+1,j+1,cov_ij);
      detvar_correlation_matrix->SetBinContent(i+1,j+1,val);
    }
  }

  h10->SetTitle("Stat");
  h20->SetTitle("MCstat");
  h30->SetTitle("Dirt");
  h40->SetTitle("Flux");
  h50->SetTitle("Det");
  h60->SetTitle("XS");
  h70->SetTitle("Relative Uncertainties");
  h70->GetYaxis()->SetTitle("Diagonal Uncertainty / Bin Content");
  h70->SetXTitle("Reco Bin Index");
  h70->GetXaxis()->SetTitleSize(.05);
  h10->SetLineColor(kCyan);
  h20->SetLineColor(kYellow);
  h30->SetLineColor(kGreen);
  h40->SetLineColor(kRed);
  h50->SetLineColor(kMagenta);
  h60->SetLineColor(kBlue);
  h70->SetLineColor(1);
  h10->SetLineWidth(2);
  h20->SetLineWidth(2);
  h30->SetLineWidth(2);
  h40->SetLineWidth(2);
  h50->SetLineWidth(2);
  h60->SetLineWidth(2);
  h70->SetLineWidth(2);
  h70->GetYaxis()->SetRangeUser(0,2.5);

  TLegend *legend1 = new TLegend(0.6,0.6,0.89,0.89);
  legend1->AddEntry(h70,"Total","l");
  legend1->AddEntry(h10,"Stat.","l");
  legend1->AddEntry(h20,"MC stat.","l");
  legend1->AddEntry(h30,"Dirt","l");
  legend1->AddEntry(h40,"Flux","l");
  legend1->AddEntry(h50,"Det.","l");
  legend1->AddEntry(h60,"Xs","l");

  //budget plot
  h10_rel->SetTitle("Stat");
  h20_rel->SetTitle("MCstat");
  h30_rel->SetTitle("Dirt");
  h40_rel->SetTitle("Flux");
  h50_rel->SetTitle("Det");
  h60_rel->SetTitle("XS");
  h70_rel->SetTitle("Relative Uncertainties");
  h70_rel->SetXTitle("Reco Bin Index");
  h10_rel->SetLineColor(kCyan);
  h20_rel->SetLineColor(kYellow);
  h30_rel->SetLineColor(kGreen);
  h40_rel->SetLineColor(kRed);
  h50_rel->SetLineColor(kMagenta);
  h60_rel->SetLineColor(kBlue);
  h70_rel->SetLineColor(1);
  h10_rel->SetFillColor(kCyan);
  h20_rel->SetFillColor(kYellow);
  h30_rel->SetFillColor(kGreen);
  h40_rel->SetFillColor(kRed);
  h50_rel->SetFillColor(kMagenta);
  h60_rel->SetFillColor(kBlue);
  h70_rel->SetFillColor(1);

  THStack *hstack = new THStack("hs","");
  hstack->Add(h30_rel);
  hstack->Add(h20_rel);
  hstack->Add(h10_rel);
  hstack->Add(h60_rel);
  hstack->Add(h40_rel);
  hstack->Add(h50_rel);

  TLegend *legend2 = new TLegend(0.6,0.6,0.89,0.89);
  legend2->AddEntry(h30_rel,"Dirt","f");
  legend2->AddEntry(h20_rel,"MC stat.","f");
  legend2->AddEntry(h10_rel,"Stat.","f");
  legend2->AddEntry(h60_rel,"Xs","f");
  legend2->AddEntry(h40_rel,"Flux","f");
  legend2->AddEntry(h50_rel,"Det.","f");
 
  
  if (draw) {
    TCanvas* c1 = new TCanvas("c1","c1");
    c1->cd();
    gPad->SetBottomMargin(0.15);
    h70->Draw();
    h10->Draw("same"); // stat
    h20->Draw("same"); // mcstat
    h30->Draw("same"); // dirt
    h40->Draw("same"); // flux
    h50->Draw("same"); // det
    h60->Draw("same"); // xs
    legend1->Draw("same");

    TCanvas* c2 = new TCanvas("c2","c2");
    c2->cd();
    gPad->SetBottomMargin(0.15);
    //hstack->GetXaxis()->SetTitleSize(.05);
    //hstack->GetXaxis()->SetTitle("Reco Bin Index");
    hstack->Draw();
    legend2->Draw("same");

    TCanvas* c_det_cov = new TCanvas("c_det_cov","c_det_cov");
    c_det_cov->cd();
    gPad->SetBottomMargin(0.15);
    detvar_covariance_matrix->GetYaxis()->SetTitleOffset(.82);
    detvar_covariance_matrix->GetYaxis()->SetTitleSize(.05);
    detvar_covariance_matrix->GetXaxis()->SetTitleSize(.05);
    detvar_covariance_matrix->GetXaxis()->SetTitle("Reco Bin Index");
    detvar_covariance_matrix->GetYaxis()->SetTitle("Reco Bin Index");
    detvar_covariance_matrix->Draw("colz");

    TCanvas* c_det_corr = new TCanvas("c_det_corr","c_det_corr");
    c_det_corr->cd();
    gPad->SetBottomMargin(0.15);
    detvar_correlation_matrix->GetYaxis()->SetTitleOffset(.82);
    detvar_correlation_matrix->GetYaxis()->SetTitleSize(.05);
    detvar_correlation_matrix->GetXaxis()->SetTitleSize(.05);
    detvar_correlation_matrix->GetXaxis()->SetTitle("Reco Bin Index");
    detvar_correlation_matrix->GetYaxis()->SetTitle("Reco Bin Index");
    detvar_correlation_matrix->Draw("colz");
  }



  if (new_noise) {
    TFile *noise_file = new TFile("noise_file.root","RECREATE");
    hnoise->SetDirectory(noise_file);
    hnoise_collapsed->SetDirectory(noise_file);
    noise_file->Write();
    noise_file->Close();
  }

 TFile *file = new TFile("wiener.root","RECREATE");
 file->cd();
 hdata_obsch_1->SetDirectory(file);
 hdata_obsch_2->SetDirectory(file);
 hmc_obsch_1->SetDirectory(file);
 hmc_obsch_2->SetDirectory(file);
 htrue_signal->SetDirectory(file);
 hmeas->SetDirectory(file);
 hpred->SetDirectory(file);
 hR->SetDirectory(file);
 hcov_stat->SetDirectory(file);
 hcov_mcstat->SetDirectory(file);
 hcov_add->SetDirectory(file);
 hcov_det->SetDirectory(file);
 hcov_flux->SetDirectory(file);
 hcov_xs->SetDirectory(file);
 hcov_tot->SetDirectory(file);

 if (!draw) {
   h10->SetDirectory(file);
   h20->SetDirectory(file);
   h30->SetDirectory(file);
   h40->SetDirectory(file);
   h50->SetDirectory(file);
   h60->SetDirectory(file);
   h70->SetDirectory(file);
   legend1->Write("legend1");

   h10_rel->SetDirectory(file);
   h20_rel->SetDirectory(file);
   h30_rel->SetDirectory(file);
   h40_rel->SetDirectory(file);
   h50_rel->SetDirectory(file);
   h60_rel->SetDirectory(file);
   h70_rel->SetDirectory(file);
   hstack->Write("hstack");
   legend2->Write("legend2");

   detvar_covariance_matrix->SetDirectory(file);
   detvar_correlation_matrix->SetDirectory(file);
 }

 file->Write();
 file->Close();

}
