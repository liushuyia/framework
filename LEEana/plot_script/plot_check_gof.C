void plot_check_gof(int lee){
  //TFile *file1 = new TFile("merge_all.root");
  TFile *file1 = new TFile("merge_all_new_weights.root");
  TMatrixD *cov_mat_add = (TMatrixD*)file1->Get("cov_mat_add");
  TMatrixD *mat_collapse_save = (TMatrixD*)file1->Get("mat_collapse");
  TH1F **pred = new TH1F*[16];
  for (Int_t i=0;i!=16;i++){
    pred[i] = (TH1F*)file1->Get(Form("histo_%d",i+1));
  }
  TH1F **hdata = new TH1F*[7];
  TH1F **hmc = new TH1F*[7];
  for (Int_t i=0;i!=7;i++){
    hdata[i] = (TH1F*)file1->Get(Form("hdata_obsch_%d",i+1));
    hmc[i] = (TH1F*)file1->Get(Form("hmc_obsch_%d",i+1));
  }

  TMatrixD *mat_collapse = new TMatrixD(26*4+11*3+26*2+26*4+11*3,26*4+33);
  *mat_collapse = *mat_collapse_save;

  if (lee == 0){ // observe LEE
    for (Int_t i=0;i!=52;i++){
      //     std::cout << (*mat_collapse)(26*4+11*3+i,i) << std::endl;
      (*mat_collapse)(26*4+11*3+i,i) = 0;
    }
  }
  
  TMatrixD *mat_collapse_T = new TMatrixD(mat_collapse->GetNcols(),mat_collapse->GetNrows());
  mat_collapse_T->Transpose(*mat_collapse);


  double fac_mc_stat = 1;
  bool flag_xf = true;
  bool flag_det = true;
  bool flag_add = true;
  bool flag_mcstat = true;

  int det_choice = 1; //1 for current, 2 for 10stat, 3 for no random
  int xf_choice = 1; // 1 for new weights, 2 for old weights
  
  // std::cout << mat_collapse->GetNrows() << " " << mat_collapse->GetNcols() << std::endl;
  // std::cout << mat_collapse_T->GetNrows() << " " << mat_collapse_T->GetNcols() << std::endl;
  // std::cout << cov_mat_add->GetNcols() << std::endl;
  
  //check input ...
  TMatrixD *pred_vec = new TMatrixD(cov_mat_add->GetNcols(),1);
  int num = 0;
  for (Int_t i=0;i!=16;i++){
    for (Int_t j=0;j!=pred[i]->GetNbinsX()+1;j++){
      (*pred_vec)(num,0) = pred[i]->GetBinContent(j+1);
      num++;
    }
  }
  TMatrixD pred_obs_vec = (*mat_collapse_T) * (*pred_vec);
  // std::cout << pred_vec->GetNrows() << " " << pred_vec->GetNcols() << std::endl;
  
  TMatrixD data_vec( pred_obs_vec.GetNrows(),1);
  num = 0;
  for (Int_t i=0;i!=7;i++){
    for (Int_t j=0;j!=hdata[i]->GetNbinsX()+1;j++){
      data_vec(num,0) = hdata[i]->GetBinContent(j+1);
      num ++;
    }
  }
  
  // std::cout << pred_obs_vec.GetNrows() << " " << pred_obs_vec.GetNcols() << std::endl;
  
  //  for (Int_t i=0;i!=26;i++){
  //  std::cout << i << " " << hmc[0]->GetBinContent(i+1) << " " << pred_obs_vec(i,0) << " " << (*data_vec)(i,0) << std::endl;
  // }

  // load MC stat covariance matrix ... (at observed space ... )
  TMatrixD *cov_mat_mcstat = new TMatrixD(pred_obs_vec.GetNrows(),pred_obs_vec.GetNrows());

  if (lee == 0){
    ifstream infile("mc_stat/0.log");
    double temp,err2;
    infile >> temp >> temp;
    
    for (Int_t i=0;i!=pred_obs_vec.GetNrows();i++){
      infile >> temp >> temp >> temp >> err2 >> temp;
      (*cov_mat_mcstat)(i,i) = err2;
    }
  }else{
    ifstream infile("mc_stat/33.log");
    double temp,err2;
    infile >> temp >> temp;
    
    for (Int_t i=0;i!=pred_obs_vec.GetNrows();i++){
      infile >> temp >> temp >> temp >> err2 >> temp;
      (*cov_mat_mcstat)(i,i) = err2;
    }
  }

  cov_mat_mcstat->Draw("COLZ");
  

  TMatrixD *cov_mat_det = new TMatrixD(cov_mat_add->GetNrows(), cov_mat_add->GetNcols());

  //Det variation ...
  if(det_choice ==1){
    TString name[9]={"hist_rootfiles/DetVar/both/cov_LYDown.root",  //0
		     "hist_rootfiles/DetVar/both/cov_LYRayleigh.root", // 1
		     "hist_rootfiles/DetVar/both/cov_Recomb2.root", // 2
		     "hist_rootfiles/DetVar/both/cov_SCE.root", // 3
		     "hist_rootfiles/DetVar/both/cov_WMdEdx.root", // 4   not useful
		     "hist_rootfiles/DetVar/both/cov_WMThetaXZ.root", // 5
		     "hist_rootfiles/DetVar/both/cov_WMThetaYZ.root", // 6
		     "hist_rootfiles/DetVar/both/cov_WMX.root",  //7 
		     "hist_rootfiles/DetVar/both/cov_WMYZ.root"  //8
    };
      for (Int_t i=0;i!=9;i++){
    if(i==4) continue;
    TFile file1(name[i]);
    TMatrixD *frac_cov_det_mat = (TMatrixD*)file1.Get(Form("frac_cov_det_mat_%d",i+1));
    for (Int_t j=0;j!=frac_cov_det_mat->GetNrows();j++){
      for (Int_t k=0;k!=frac_cov_det_mat->GetNcols();k++){
	(*frac_cov_det_mat)(j,k) *= (*pred_vec)(j,0) * (*pred_vec)(k,0);
      }
    }
    (*cov_mat_det) += (*frac_cov_det_mat);
  }

  }else if (det_choice == 2){
     TString name[9]={"hist_rootfiles/DetVar/10stat/cov_LYDown.root",  //0
  		   "hist_rootfiles/DetVar/10stat/cov_LYRayleigh.root", // 1
  		   "hist_rootfiles/DetVar/10stat/cov_Recomb2.root", // 2
  		   "hist_rootfiles/DetVar/10stat/cov_SCE.root", // 3
  		   "hist_rootfiles/DetVar/10stat/cov_WMdEdx.root", // 4   not useful
  		   "hist_rootfiles/DetVar/10stat/cov_WMThetaXZ.root", // 5
  		   "hist_rootfiles/DetVar/10stat/cov_WMThetaYZ.root", // 6
  		   "hist_rootfiles/DetVar/10stat/cov_WMX.root",  //7 
  		   "hist_rootfiles/DetVar/10stat/cov_WMYZ.root"  //8
     };
       for (Int_t i=0;i!=9;i++){
    if(i==4) continue;
    TFile file1(name[i]);
    TMatrixD *frac_cov_det_mat = (TMatrixD*)file1.Get(Form("frac_cov_det_mat_%d",i+1));
    for (Int_t j=0;j!=frac_cov_det_mat->GetNrows();j++){
      for (Int_t k=0;k!=frac_cov_det_mat->GetNcols();k++){
	(*frac_cov_det_mat)(j,k) *= (*pred_vec)(j,0) * (*pred_vec)(k,0);
      }
    }
    (*cov_mat_det) += (*frac_cov_det_mat);
  }

  }else if (det_choice ==3){
    TString name[9]={"hist_rootfiles/DetVar/no_random/cov_LYDown.root",  //0
		     "hist_rootfiles/DetVar/no_random/cov_LYRayleigh.root", // 1
		     "hist_rootfiles/DetVar/no_random/cov_Recomb2.root", // 2
		     "hist_rootfiles/DetVar/no_random/cov_SCE.root", // 3
		     "hist_rootfiles/DetVar/no_random/cov_WMdEdx.root", // 4   not useful
		     "hist_rootfiles/DetVar/no_random/cov_WMThetaXZ.root", // 5
		     "hist_rootfiles/DetVar/no_random/cov_WMThetaYZ.root", // 6
		     "hist_rootfiles/DetVar/no_random/cov_WMX.root",  //7 
		     "hist_rootfiles/DetVar/no_random/cov_WMYZ.root"  //8
    };
      for (Int_t i=0;i!=9;i++){
    if(i==4) continue;
    TFile file1(name[i]);
    TMatrixD *frac_cov_det_mat = (TMatrixD*)file1.Get(Form("frac_cov_det_mat_%d",i+1));
    for (Int_t j=0;j!=frac_cov_det_mat->GetNrows();j++){
      for (Int_t k=0;k!=frac_cov_det_mat->GetNcols();k++){
	(*frac_cov_det_mat)(j,k) *= (*pred_vec)(j,0) * (*pred_vec)(k,0);
      }
    }
    (*cov_mat_det) += (*frac_cov_det_mat);
  }

  }
  
 

  cov_mat_det->Draw("COLZ");

  TMatrixD *cov_mat_xf = new TMatrixD(cov_mat_add->GetNrows(), cov_mat_add->GetNcols());

  if (xf_choice==1){
    // Xs and Flux covariance matrix ...
    for (Int_t i=0;i!=17;i++){
      TFile file1(Form("hist_rootfiles/XsFlux/cov_%d.root",i+1));
      TMatrixD *frac_cov_xf_mat = (TMatrixD*)file1.Get(Form("frac_cov_xf_mat_%d",i+1));
      for (Int_t j=0;j!=frac_cov_xf_mat->GetNrows();j++){
	for (Int_t k=0;k!=frac_cov_xf_mat->GetNcols();k++){
	  (*frac_cov_xf_mat)(j,k) *= (*pred_vec)(j,0) * (*pred_vec)(k,0);
	}
      }
      (*cov_mat_xf) += (*frac_cov_xf_mat);
    }
  }else{
    for (Int_t i=0;i!=14;i++){
      TFile file1(Form("hist_rootfiles/XsFlux/old_weights/cov_%d.root",i+1));
      TMatrixD *frac_cov_xf_mat = (TMatrixD*)file1.Get(Form("frac_cov_xf_mat_%d",i+1));
      for (Int_t j=0;j!=frac_cov_xf_mat->GetNrows();j++){
	for (Int_t k=0;k!=frac_cov_xf_mat->GetNcols();k++){
	  (*frac_cov_xf_mat)(j,k) *= (*pred_vec)(j,0) * (*pred_vec)(k,0);
	}
      }
      (*cov_mat_xf) += (*frac_cov_xf_mat);
    }
  }
    
  cov_mat_xf->Draw("COLZ");

  TMatrixD *cov_mat_tot_before = new TMatrixD(cov_mat_add->GetNrows(), cov_mat_add->GetNcols());
  
  if (flag_xf) (*cov_mat_tot_before) += (*cov_mat_xf);
  if (flag_det) (*cov_mat_tot_before) += (*cov_mat_det);
  if (flag_add) {
    (*cov_mat_tot_before) += (*cov_mat_add); // additional covariance matrix ...
    // for (Int_t i=0;i!=cov_mat_add->GetNcols();i++){
    //   std::cout << i << " " << (*cov_mat_add)(i,i) << std::endl;
    // }
  }

  TMatrixD cov_mat_tot = (*mat_collapse_T) * (*cov_mat_tot_before) * (*mat_collapse);
  if (flag_mcstat) cov_mat_tot += (*cov_mat_mcstat); // mc statistical uncertainties ...
  
  cov_mat_tot.Draw("COLZ");


  // add stat term ...
  for    (Int_t i=0;i!=cov_mat_tot.GetNcols();i++){
    double tmp = cov_mat_tot(i,i);
    //    if (i>=104 && i<104+11){
    //  std::cout << i << " " <<cov_mat_tot(i,i) << std::endl;
    // }
    
    if ( cov_mat_tot(i,i)==0 && pred_obs_vec(i,0) ==0){
    // if ( data_vec(i,0)==0 && pred_obs_vec(i,0) ==0){
      for (Int_t j=0;j!=cov_mat_tot.GetNrows();j++){
	cov_mat_tot(j,i) = 0;
	cov_mat_tot(i,j) = 0;
      }
      cov_mat_tot(i,i) = 1e9; // very large uncertainties ...
      std::cout << "Bin: " << i << " not good" << std::endl;
    }else{
      if (pred_obs_vec(i,0)>=0.461 || data_vec(i,0)==0){
	cov_mat_tot(i,i) += pred_obs_vec(i,0); // add stat term ...
      }else{
	cov_mat_tot(i,i) += pow(pred_obs_vec(i,0) - data_vec(i,0),2)/(2.*(pred_obs_vec(i,0) - data_vec(i,0) + data_vec(i,0) * log(data_vec(i,0)/pred_obs_vec(i,0))));
	//cov_mat_tot(i,i) += pred_obs_vec(i,0); // add stat term ...
      }
    }

    // if (i>=104 && i<104+11){
    //   std::cout << i << " " <<cov_mat_tot(i,i) - tmp << std::endl;
    // }
  }

  
  {
    // goodness of numu CC FC, pick a single channel ...
    int offset = 0;
    int nbin = 52+52+33;
    
    TMatrixD transform(cov_mat_tot.GetNrows(),nbin);
    for (Int_t i = 0; i!= nbin;i++){
      transform(offset+i,i)= 1; //  numu CC FC
    }
    TMatrixD transform_T(nbin,cov_mat_tot.GetNrows());
    transform_T.Transpose(transform);
    TMatrixD cov_temp = transform_T * (cov_mat_tot) * transform;

    
    TMatrixD diff_vec(nbin,1);
    TMatrixD diff_vec_T(1,nbin);
    for (Int_t i=0;i!=nbin;i++){
      diff_vec(i,0) = pred_obs_vec(offset+i,0) - data_vec(offset+i,0);
      diff_vec_T(0,i) = pred_obs_vec(offset+i,0) - data_vec(offset+i,0);
      //std::cout << i << " " << pred_obs_vec(offset+i,0) << " " <<  data_vec(offset+i,0) << " " << diff_vec(i,0) << std::endl;
    }
    TMatrixD cov_temp_invert = cov_temp;
    cov_temp_invert.Invert();
    
    TMatrixD chi2 = diff_vec_T * cov_temp_invert * diff_vec;
    std::cout << "numuCC Gof: " <<  chi2(0,0) << std::endl;
    
    cov_temp.Draw("COLZ");
  }

  
  {
    // conditional covariance after numuCC measurements ...
    TMatrixD cov_xx(26+26,26+26);
    TMatrixD cov_yy(26+26+11*3, 26+26+11*3);
    TMatrixD cov_yx(26+26+11*3,26+26);
    TMatrixD cov_xy(26+26,26+26+11*3);
    for (Int_t i=0;i!=52;i++){
      for (Int_t j=0;j!=52;j++){
  	cov_xx(i,j) = cov_mat_tot(52+i,52+j); // numuCC
	cov_yy(i,j) = cov_mat_tot(i,j);  // nueCC

	cov_xy(i,j) = cov_mat_tot(52+i,j); // numu CC -- >nueCC
       	cov_yx(j,i) = cov_mat_tot(j,52+i); // nueCC --> numuCC
      }
      for (Int_t j=0;j!=11*3;j++){
	cov_yy(i,j+52) = cov_mat_tot(i,52+52+j);  // nueCC-> pi0
	cov_yy(j+52,i) = cov_mat_tot(j+52+52,i); // pi0 --> nueCC

	cov_xy(i,52+j) = cov_mat_tot(52+i,52+52+j); //numuCC --> pi0
	cov_yx(52+j,i) = cov_mat_tot(52+52+j,52+i); // pi0 --> numuCC
      }
    }

    for (Int_t i=0;i!=33;i++){
      for (Int_t j=0;j!=33;j++){
	cov_yy(52+i,52+j) = cov_mat_tot(104+i,104+j); // pi0
      }
    }
    TMatrixD cov_xx_inv = cov_xx;
    cov_xx_inv.Invert();
    TMatrixD cov_yy_inv = cov_yy;
    cov_yy_inv.Invert();

    TMatrixD vec_data_cont1(26+26,1);
    TMatrixD vec_data_cont2(26+26+11*3,1);

    TMatrixD vec_pred_cont1(26+26,1);
    TMatrixD vec_pred_cont2(26+26+11*3,1);

    for (Int_t i=0;i!=52;i++){
      vec_data_cont1(i,0) = data_vec(i+52,0);
      vec_pred_cont1(i,0) = pred_obs_vec(i+52,0);
      vec_data_cont2(i,0) = data_vec(i,0);
      vec_pred_cont2(i,0) = pred_obs_vec(i,0);
    }
    for (Int_t i=0;i!=33;i++){
      vec_data_cont2(52+i,0) = data_vec(104+i,0);
      vec_pred_cont2(52+i,0) = pred_obs_vec(104+i,0);
    }
    
    
    TMatrixD vec_after = vec_pred_cont2 + cov_yx * cov_xx_inv * (vec_data_cont1 - vec_pred_cont1);

    // for (Int_t i=0;i!=11;i++){
    //   std::cout << vec_pred_cont2(52+i,0) << " " << cov_yx(52+i,i) << " " << cov_xx_inv(i,i) << " " << vec_data_cont1(i,0) << " " << vec_pred_cont1(i,0) << std::endl;
    // }
    
    
    TMatrixD cov_after = cov_yy - cov_yx  * cov_xx_inv * cov_xy;

    // for (int i=0;i!=11;i++){
    //   std::cout << i << " " << cov_yy(52+i,52+i) << " " << cov_yx(52+i,i) << " " << cov_xx_inv(i,i) << " " << cov_xy(i,52+i) << " " << cov_after(52+i,52+i) << std::endl;
    // }
    
    TMatrixD diff_vec(52+33,1);
    TMatrixD diff_vec_T(1,52+33);

    // numuCC ...
    int nbin = 11*3;
    int offset = 52;

    for (Int_t i=0;i!=cov_after.GetNrows();i++){
      cov_after(i,i) += vec_after(i,0) - vec_pred_cont2(i,0);
      //      if (i>=52 && i<52+11)
	//std::cout << cov_after(i,i) << std::endl;
    }
    
    for (Int_t i=0;i!=nbin;i++){
      diff_vec(offset+i,0) = vec_data_cont2(offset+i,0) - vec_after(offset+i,0);
      diff_vec_T(0,offset+i) = vec_data_cont2(offset+i,0) - vec_after(offset+i,0);
      //      std::cout << i << " " << vec_data_cont2(offset+i,0)<< " " << vec_after(offset+i,0) << " " << vec_pred_cont2(offset+i,0) << " " << diff_vec(offset+i,0) << std::endl;
      
    }

    // for (Int_t i=0;i!=52;i++){
    //   for (Int_t j=0;j!=33;j++){
    // 	cov_after(i,j+52) = 0;
    // 	cov_after(j+52,i) = 0;
    //   }
    // }
    // for (Int_t i=0;i!=11;i++){
    //   for (Int_t j=0;j!=22;j++){
    // 	cov_after(52+i,52+11+j) = 0;
    // 	cov_after(52+11+j,52+i) = 0;
    //   }
    // }
    // for (Int_t i=0;i!=cov_after.GetNrows();i++){
    //   if (cov_after(i,i) ==0)  cov_after(i,i) = 1;
    // }
    
    
    TMatrixD cov_after_inv = cov_after;
    cov_after_inv.Invert();

    
    
    TMatrixD chi2 = diff_vec_T * cov_after_inv * diff_vec;
    std::cout << "Pi0 Gof: " << chi2(0,0) << std::endl;
    
  }

  
  // constrainting numu and pi0 ...
  {
     // conditional covariance after numuCC measurements ...
    TMatrixD cov_xx(26*4+11*3-52,26*4+11*3-52); 
    TMatrixD cov_yy(52, 52);
    TMatrixD cov_yx(52,26*4+11*3-52);
    TMatrixD cov_xy(26*4+11*3-52,52);

     for (Int_t i=0;i!=52;i++){
      for (Int_t j=0;j!=52;j++){
	cov_yy(i,j) = cov_mat_tot(i,j);  // nueCC
      }
    }

     for (Int_t i=0;i!=26*4+3*11-52;i++){
       for (Int_t j=0; j!= 26*4+3*11-52;j++){
	 cov_xx(i,j) = cov_mat_tot(i+52,j+52);
       }
     }
     
     for (Int_t i=0;i!=52;i++){
       for (Int_t j=0; j!= 26*4+3*11-52;j++){
	 cov_xy(j,i) = cov_mat_tot(i,j+52);
	 cov_yx(i,j) = cov_mat_tot(j+52,i);
       }
     }
     
    TMatrixD cov_xx_inv = cov_xx;
    cov_xx_inv.Invert();
    TMatrixD cov_yy_inv = cov_yy;
    cov_yy_inv.Invert();

    TMatrixD vec_data_cont1(26+26+11*3,1);
    TMatrixD vec_data_cont2(26+26,1);

    TMatrixD vec_pred_cont1(26+26+11*3,1);
    TMatrixD vec_pred_cont2(26+26,1);

    for (Int_t i=0;i!=52;i++){
      vec_data_cont2(i,0) = data_vec(i,0);
      vec_pred_cont2(i,0) = pred_obs_vec(i,0);
    }
    for (Int_t i=0;i!=33+52;i++){
      vec_data_cont1(i,0) = data_vec(52+i,0);
      vec_pred_cont1(i,0) = pred_obs_vec(52+i,0);
    }
    
    
    TMatrixD vec_after = vec_pred_cont2 + cov_yx * cov_xx_inv * (vec_data_cont1 - vec_pred_cont1);

    TMatrixD cov_after = cov_yy - cov_yx  * cov_xx_inv * cov_xy;

    TMatrixD diff_vec(52,1);
    TMatrixD diff_vec_T(1,52);
    for (Int_t i=0;i!=52;i++){
      diff_vec(i,0) = vec_data_cont2(i,0) - vec_after(i,0);
      diff_vec_T(0,i) = vec_data_cont2(i,0) - vec_after(i,0);

      // original ...
      cov_after(i,i) += vec_after(i,0) - vec_pred_cont2(i,0);
    }
    TMatrixD cov_after_inv = cov_after;
    cov_after_inv.Invert();
    
    TMatrixD chi2 = diff_vec_T * cov_after_inv * diff_vec;
    std::cout << "full nueCC Gof: " << chi2(0,0) << std::endl; // 
  }



  {
    TMatrixD cov_xx(26*4+3*11-8,26*4+3*11-8);
    TMatrixD cov_yy(8, 8);
    TMatrixD cov_yx(8, 26*4+3*11-8);
    TMatrixD cov_xy(26*4+3*11-8,8);

    for (Int_t i=0;i!=8;i++){
      for (Int_t j=0;j!=8;j++){
	cov_yy(i,j) = cov_mat_tot(i,j);  // nueCC
      }
    }

    for (Int_t i=0;i!=26*4+3*11-8;i++){
      for (Int_t j=0; j!= 26*4+3*11-8;j++){
	cov_xx(i,j) = cov_mat_tot(i+8,j+8);
      }
    }

    for (Int_t i=0;i!=8;i++){
      for (Int_t j=0; j!= 26*4+3*11-8;j++){
	cov_xy(j,i) = cov_mat_tot(i,j+8);
	cov_yx(i,j) = cov_mat_tot(j+8,i);
      }
    }

    TMatrixD cov_xx_inv = cov_xx;
    cov_xx_inv.Invert();
    TMatrixD cov_yy_inv = cov_yy;
    cov_yy_inv.Invert();

    TMatrixD vec_data_cont1(26*4+11*3-8,1);
    TMatrixD vec_data_cont2(8,1);

    TMatrixD vec_pred_cont1(26*4+11*3-8,1);
    TMatrixD vec_pred_cont2(8,1);

    for (Int_t i=0;i!=8;i++){
      vec_data_cont2(i,0) = data_vec(i,0);
      vec_pred_cont2(i,0) = pred_obs_vec(i,0);
    }
    for (Int_t i=0;i!=26*4+3*11-8;i++){
      vec_data_cont1(i,0) = data_vec(8+i,0);
      vec_pred_cont1(i,0) = pred_obs_vec(8+i,0);
    }
    
    
    TMatrixD vec_after = vec_pred_cont2 + cov_yx * cov_xx_inv * (vec_data_cont1 - vec_pred_cont1);
    TMatrixD cov_after = cov_yy - cov_yx  * cov_xx_inv * cov_xy;

    TMatrixD diff_vec(8,1);
    TMatrixD diff_vec_T(1,8);
    for (Int_t i=0;i!=8;i++){
      // vec_data_cont2(i,0) = 0;
      diff_vec(i,0) = vec_data_cont2(i,0) - vec_after(i,0);
      diff_vec_T(0,i) = vec_data_cont2(i,0) - vec_after(i,0);

      // original ...
      cov_after(i,i) += vec_after(i,0) - vec_pred_cont2(i,0);

      //std::cout << diff_vec(i,0) << " " << vec_data_cont2(i,0) << " " <<  vec_after(i,0) << std::endl;
    }
    TMatrixD cov_after_inv = cov_after;
    cov_after_inv.Invert();
    
    TMatrixD chi2 = diff_vec_T * cov_after_inv * diff_vec;
    std::cout << "low-energy nueCC Gof: " << chi2(0,0) << std::endl; // 

    // merge ... 
    TMatrixD transform(8,1);
    TMatrixD transform_T(1,8);
    for (Int_t i=0;i!=8;i++){
      // if (i>=6) continue;
      transform(i,0) = 1;
      transform_T(0,i) = 1;
    }
    TMatrixD cov_merge = transform_T * cov_after * transform;
    TMatrixD cov_merge_inv = cov_merge;
    cov_merge_inv.Invert();
    TMatrixD diff_merge = transform_T* diff_vec;
    std::cout << "Merge: " << diff_merge(0,0) * diff_merge(0,0) * cov_merge_inv(0,0) << std::endl;

    
    
  }

  double data_pot = 5.0549670e+19 ;
  
  {
    // 5.0549670e+19  Asimov sensitivity ...
    TMatrixD *mat_collapse_lee = new TMatrixD(26*4+11*3+26*2+26*4+11*3,26*4+33);
    *mat_collapse_lee = *mat_collapse_save;

    TMatrixD *mat_collapse_nolee = new TMatrixD(26*4+11*3+26*2+26*4+11*3,26*4+33);
    *mat_collapse_nolee = *mat_collapse_save;
    
    for (Int_t i=0;i!=52;i++){
      //     std::cout << (*mat_collapse)(26*4+11*3+i,i) << std::endl;
      (*mat_collapse_nolee)(26*4+11*3+i,i) = 0;
    }

    TMatrixD *mat_collapse_lee_T = new TMatrixD(mat_collapse->GetNcols(),mat_collapse->GetNrows());
    mat_collapse_lee_T->Transpose(*mat_collapse_lee);
    
    TMatrixD *mat_collapse_nolee_T = new TMatrixD(mat_collapse->GetNcols(),mat_collapse->GetNrows());
    mat_collapse_nolee_T->Transpose(*mat_collapse_nolee);
  
    //check input ...
    
    TMatrixD pred_obs_lee_vec = (*mat_collapse_lee_T) * (*pred_vec);
    TMatrixD pred_obs_nolee_vec = (*mat_collapse_nolee_T) * (*pred_vec);

    //pred_obs_lee_vec = data_vec;
    //    pred_obs_nolee_vec = data_vec;
    
    TMatrixD diff = pred_obs_lee_vec  - pred_obs_nolee_vec;

    // hack ...
    
    
    TMatrixD diff_T(1,diff.GetNrows());
    diff_T.Transpose(diff);
    
    
    if (lee == 0) // observe LEE
      {
	// no LEE
	TMatrixD cov_mat_tot = (*mat_collapse_nolee_T) * (*cov_mat_tot_before) * (*mat_collapse_nolee);
	cov_mat_tot += (*cov_mat_mcstat); // mc statistical uncertainties ...
	
	// add statistical uncertainties with CNP ...
	for (Int_t i=0;i!=cov_mat_tot.GetNcols();i++){
	  if ( pred_obs_lee_vec(i,0)==0 && pred_obs_nolee_vec(i,0) ==0){
	    for (Int_t j=0;j!=cov_mat_tot.GetNrows();j++){
	      cov_mat_tot(j,i) = 0;
	      cov_mat_tot(i,j) = 0;
	    }
	    cov_mat_tot(i,i) = 1e9; // very large uncertainties ...
	    
	  std::cout << "Bin: " << i << " not good" << std::endl;
	  }else{
	    
	    if (pred_obs_lee_vec(i,0) !=0)
	      cov_mat_tot(i,i) += 3./(1./pred_obs_lee_vec(i,0) + 2./pred_obs_nolee_vec(i,0));
	    else
	      cov_mat_tot(i,i) += pred_obs_nolee_vec(i,0)/2.;
	  }
	}
	
	TMatrixD cov_mat_tot_inv = cov_mat_tot;
	cov_mat_tot_inv.Invert();
	TMatrixD chi2 = diff_T * cov_mat_tot_inv * diff;
	
	std::cout << "Sensitivity to exclude SM at current POT: " << chi2(0,0) << std::endl;
      }
    
    if (lee == 1) // observe no LEE
      {
	// LEE
	TMatrixD cov_mat_tot = (*mat_collapse_lee_T) * (*cov_mat_tot_before) * (*mat_collapse_lee);
	cov_mat_tot += (*cov_mat_mcstat); // mc statistical uncertainties ... suppress by 50%
	
	// add statistical uncertainties with CNP ...
	for (Int_t i=0;i!=cov_mat_tot.GetNcols();i++){
	  if ( pred_obs_lee_vec(i,0)==0 && pred_obs_nolee_vec(i,0) ==0){
	    for (Int_t j=0;j!=cov_mat_tot.GetNrows();j++){
	      cov_mat_tot(j,i) = 0;
	      cov_mat_tot(i,j) = 0;
	    }
	    cov_mat_tot(i,i) = 1e9; // very large uncertainties ...
	    std::cout << "Bin: " << i << " not good" << std::endl;
	  }else{
	    if (pred_obs_nolee_vec(i,0) !=0)
	      cov_mat_tot(i,i) += 3./(2./pred_obs_lee_vec(i,0) + 1./pred_obs_nolee_vec(i,0));
	    else
	      cov_mat_tot(i,i) += pred_obs_lee_vec(i,0)/2.;
	  }
	}
	
	TMatrixD cov_mat_tot_inv = cov_mat_tot;
	cov_mat_tot_inv.Invert();
	TMatrixD chi2 = diff_T * cov_mat_tot_inv * diff;
	
	std::cout << "Sensitivity to exclude LEE at current POT: " << chi2(0,0) << std::endl;
	
      }
    
    if (lee == 0) // observe LEE
      {
	// no LEE
	TMatrixD cov_mat_tot = (*mat_collapse_nolee_T) * (*cov_mat_tot_before) * (*mat_collapse_nolee);
	cov_mat_tot += (*cov_mat_mcstat)*fac_mc_stat; // mc statistical uncertainties ... suppress by 50%
	
	cov_mat_tot *= pow(6.95e20/data_pot,2);
	
	// add statistical uncertainties with CNP ...
	for (Int_t i=0;i!=cov_mat_tot.GetNcols();i++){
	  if ( pred_obs_lee_vec(i,0)==0 && pred_obs_nolee_vec(i,0) ==0){
	    for (Int_t j=0;j!=cov_mat_tot.GetNrows();j++){
	      cov_mat_tot(j,i) = 0;
	      cov_mat_tot(i,j) = 0;
	    }
	    cov_mat_tot(i,i) = 1e9; // very large uncertainties ...
	    std::cout << "Bin: " << i << " not good" << std::endl;
	  }else{
	    if (pred_obs_lee_vec(i,0) !=0)
	      cov_mat_tot(i,i) += 3./(1./(pred_obs_lee_vec(i,0)*6.95e20/data_pot) + 2./(pred_obs_nolee_vec(i,0)*6.95e20/data_pot));
	    else
	      cov_mat_tot(i,i) += pred_obs_nolee_vec(i,0)*6.95e20/data_pot/2.;
	  }
	}
	
	diff *= 6.95e20/data_pot;
	diff_T *= 6.95e20/data_pot;
	
	TMatrixD cov_mat_tot_inv = cov_mat_tot;
	cov_mat_tot_inv.Invert();
	TMatrixD chi2 = diff_T * cov_mat_tot_inv * diff;
	
	std::cout << "Sensitivity to exclude SM at 6.95e20: " << chi2(0,0) << std::endl;
      }

    if (lee == 1) // observe no LEE
    {
      // LEE
      TMatrixD cov_mat_tot = (*mat_collapse_lee_T) * (*cov_mat_tot_before) * (*mat_collapse_lee);
      cov_mat_tot += (*cov_mat_mcstat)*fac_mc_stat; // mc statistical uncertainties ...

      cov_mat_tot *= pow(6.95e20/data_pot,2);
      
      // add statistical uncertainties with CNP ...
      for (Int_t i=0;i!=cov_mat_tot.GetNcols();i++){
	if ( pred_obs_lee_vec(i,0) ==0&& pred_obs_nolee_vec(i,0) ==0){
	  for (Int_t j=0;j!=cov_mat_tot.GetNrows();j++){
	    cov_mat_tot(j,i) = 0;
	    cov_mat_tot(i,j) = 0;
	  }
	  cov_mat_tot(i,i) = 1e9; // very large uncertainties ...
	  std::cout << "Bin: " << i << " not good" << std::endl;
	}else{
	  if (pred_obs_nolee_vec(i,0) !=0)
	    cov_mat_tot(i,i) += 3./(2./(pred_obs_lee_vec(i,0)*6.95e20/data_pot) + 1./(pred_obs_nolee_vec(i,0)*6.95e20/data_pot));
	  else
	    cov_mat_tot(i,i) += pred_obs_lee_vec(i,0)/2.*6.95e20/data_pot;
	}
      }

      TMatrixD cov_mat_tot_inv = cov_mat_tot;
      cov_mat_tot_inv.Invert();

      diff *= 6.95e20/data_pot;
      diff_T *= 6.95e20/data_pot;

      TMatrixD chi2 = diff_T * cov_mat_tot_inv * diff;

      std::cout << "Sensitivity to exclude LEE at 6.95e20: " << chi2(0,0) << std::endl;
    }
    


    
  }
  
  
}
