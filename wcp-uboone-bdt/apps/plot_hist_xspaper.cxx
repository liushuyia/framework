#include <iostream>
#include <sstream>
#include <string>
#include <map>
#include <utility>

#include "WCPLEEANA/master_cov_matrix.h"
#include "WCPLEEANA/bayes.h"
#include "WCPLEEANA/eval.h"
#include "WCPLEEANA/tagger.h"

#include "TROOT.h"
#include "TApplication.h"
#include "TCanvas.h"
#include "TImage.h"
#include "TGraphErrors.h"
#include "TGraphAsymmErrors.h" 
#include "TAxis.h"
#include "TLegend.h"
#include "TMath.h"
#include "TH1F.h"
#include "TFile.h"
#include "TTree.h"
#include "TStyle.h"
#include "THStack.h"
#include "TLegend.h"
#include "TString.h"
#include "TLine.h"
#include "TPDF.h"
#include "TF1.h"
#include "TMatrixD.h"

using namespace std;
using namespace LEEana;

float eff_error(float s, float sigma_s, float b, float sigma_b)
{
    if(s+b==0) return 0;
    return TMath::Sqrt( TMath::Power(sigma_s*b/((s+b)*(s+b)), 2) + TMath::Power(sigma_b*s/((s+b)*(s+b)),2) );
}


int main( int argc, char** argv )
{
  
  if (argc < 2){
    std::cout << "./plot_hist -r[#run, 0 for all] -e[1 for standard, 2 for Bayesian, 3 for total uncertainty] -l[LEE strength, check against total uncertainty config] -s[path-to-external covmatrix file] -c[check MC and DATA spectra against the ones from external file]" << std::endl;
  }
  
  int run = 1; // run 1 ...
  int flag_err = 1;// 1 for standard, 2 for Bayesian, 3 for breakdown
  std::map<int, double> sumtotalcov;
  std::map<int, std::pair<double, int>> GOF; 
  float lee_strength = 0; // no LEE strength ...
  int flag_display = 1;
  int flag_breakdown = 1;
  TString cov_inputfile = "";
  int flag_check = 0;
  int flag_truthlabel = 0;

  for (Int_t i=1;i!=argc;i++){
    switch(argv[i][1]){
    case 'r':{
      run = atoi(&argv[i][2]); // which run period
    }break;
    case 'e':{
      flag_err = atoi(&argv[i][2]); // error for plotting
    }break;
    case 'l':{
      lee_strength = atof(&argv[i][2]);
    }break;
    case 'd':{
      flag_display = atoi(&argv[i][2]);
    }break;
    case 'b':{
      flag_breakdown = atoi(&argv[i][2]);
    }break;
    case 's':{
      TString sss = argv[i];
      cov_inputfile = sss(2, sss.Length()-2);
    }break;
    case 'c':{
      flag_check = atoi(&argv[i][2]);
    }break;
    case 't':{
      flag_truthlabel = atoi(&argv[i][2]);
    }break;
    }
  }

  // ATTENTION!
  // mannually added pad3 and inset response matrix in pad1
  // configurations/tox_xs/xs_real_bin.txt
  Double_t xbins_true[11] = {200, 540, 705, 805, 920, 1050, 1200, 1375, 1570, 2050, 4000}; // MeV 
  float xbins_const[10] = {11909.4, 6539.56, 3809.78, 4065.94, 4074.76, 3933.4, 3507.89, 2557.06, 2272.93, 700.703};
  
  TFile* file1 = new TFile("./processed_checkout_rootfiles/checkout_prodgenie_bnb_nu_overlay_run1.root");
  TFile* file2 = new TFile("./processed_checkout_rootfiles/checkout_prodgenie_bnb_nu_overlay_run3.root");
  TTree* T_eval1 = (TTree*)file1->Get("wcpselection/T_eval");
  TTree* T_eval2 = (TTree*)file2->Get("wcpselection/T_eval");
  TTree* T_BDTvars1 = (TTree*)file1->Get("wcpselection/T_BDTvars");
  TTree* T_BDTvars2 = (TTree*)file2->Get("wcpselection/T_BDTvars");

  Bool_t match_isFC;
  Float_t weight_spline;
  Float_t weight_cv;
  Bool_t truth_isCC;
  Int_t truth_nuPdg;
  Bool_t truth_vtxInside;
  Float_t truth_nuEnergy;
  Float_t truth_energyInside;
  Float_t match_completeness_energy;


  T_eval1->SetBranchStatus("*",0);
  T_eval1->SetBranchStatus("match_isFC",1);
  T_eval1->SetBranchStatus("weight_spline",1);
  T_eval1->SetBranchStatus("weight_cv",1);
  T_eval1->SetBranchStatus("truth_isCC",1);
  T_eval1->SetBranchStatus("truth_nuPdg",1);
  T_eval1->SetBranchStatus("truth_vtxInside",1);
  T_eval1->SetBranchStatus("truth_nuEnergy",1);
  T_eval1->SetBranchStatus("truth_energyInside",1);
  T_eval1->SetBranchStatus("match_completeness_energy",1);

  T_eval1->SetBranchAddress("match_isFC", &match_isFC);
  T_eval1->SetBranchAddress("weight_spline",&weight_spline);
  T_eval1->SetBranchAddress("weight_cv", &weight_cv);
  T_eval1->SetBranchAddress("truth_isCC",&truth_isCC);
  T_eval1->SetBranchAddress("truth_nuPdg", &truth_nuPdg);
  T_eval1->SetBranchAddress("truth_vtxInside", &truth_vtxInside);
  T_eval1->SetBranchAddress("truth_nuEnergy",&truth_nuEnergy);
  T_eval1->SetBranchAddress("truth_energyInside",&truth_energyInside);
  T_eval1->SetBranchAddress("match_completeness_energy",&match_completeness_energy);
  
  T_eval2->SetBranchStatus("*",0);
  T_eval2->SetBranchStatus("match_isFC",1);
  T_eval2->SetBranchStatus("weight_spline",1);
  T_eval2->SetBranchStatus("weight_cv",1);
  T_eval2->SetBranchStatus("truth_isCC",1);
  T_eval2->SetBranchStatus("truth_nuPdg",1);
  T_eval2->SetBranchStatus("truth_vtxInside",1);
  T_eval2->SetBranchStatus("truth_nuEnergy",1);
  T_eval2->SetBranchStatus("truth_energyInside",1);
  T_eval2->SetBranchStatus("match_completeness_energy",1);

  T_eval2->SetBranchAddress("match_isFC", &match_isFC);
  T_eval2->SetBranchAddress("weight_spline",&weight_spline);
  T_eval2->SetBranchAddress("weight_cv", &weight_cv);
  T_eval2->SetBranchAddress("truth_isCC",&truth_isCC);
  T_eval2->SetBranchAddress("truth_nuPdg", &truth_nuPdg);
  T_eval2->SetBranchAddress("truth_vtxInside", &truth_vtxInside);
  T_eval2->SetBranchAddress("truth_nuEnergy",&truth_nuEnergy);
  T_eval2->SetBranchAddress("truth_energyInside",&truth_energyInside);
  T_eval2->SetBranchAddress("match_completeness_energy",&match_completeness_energy);

  Float_t numu_cc_flag;
  Float_t numu_score;

  T_BDTvars1->SetBranchStatus("*",0);
  T_BDTvars1->SetBranchStatus("numu_cc_flag",1);
  T_BDTvars1->SetBranchStatus("numu_score",1);

  T_BDTvars1->SetBranchAddress("numu_cc_flag",&numu_cc_flag);
  T_BDTvars1->SetBranchAddress("numu_score", &numu_score);

  T_BDTvars2->SetBranchStatus("*",0);
  T_BDTvars2->SetBranchStatus("numu_cc_flag",1);
  T_BDTvars2->SetBranchStatus("numu_score",1);

  T_BDTvars2->SetBranchAddress("numu_cc_flag",&numu_cc_flag);
  T_BDTvars2->SetBranchAddress("numu_score", &numu_score);


  TH1F* hpass_w = new TH1F("hpass_w", "", 10, xbins_true);
  TH1F* hpass_w2 = new TH1F("hpass_w2", "", 10, xbins_true);
  TH1F* hfail_w = new TH1F("hfail_w", "", 10, xbins_true);
  TH1F* hfail_w2 = new TH1F("hfail_w2", "", 10, xbins_true);
  TH1F* hFCpass_w = new TH1F("hFCpass_w", "", 10, xbins_true);
  TH1F* hFCpass_w2 = new TH1F("hFCpass_w2", "", 10, xbins_true);
  TH1F* hFCfail_w = new TH1F("hFCfail_w", "", 10, xbins_true);
  TH1F* hFCfail_w2 = new TH1F("hFCfail_w2", "", 10, xbins_true);

  for(int i=0; i<T_eval1->GetEntries(); i++)
  {
      T_eval1->GetEntry(i);
      T_BDTvars1->GetEntry(i);
      if(match_completeness_energy>0.1*truth_energyInside && truth_nuPdg==14 && truth_isCC==1 && truth_vtxInside==1 && truth_nuEnergy<=4000 && truth_nuEnergy>200)
      {
          if(numu_cc_flag>=0 && numu_score>0.9)
          {
            hpass_w->Fill(truth_nuEnergy, weight_cv*weight_spline);
            hpass_w2->Fill(truth_nuEnergy, pow(weight_cv*weight_spline, 2));
            if(match_isFC==1)
            {
                hFCpass_w->Fill(truth_nuEnergy, weight_cv*weight_spline);
                hFCpass_w2->Fill(truth_nuEnergy, pow(weight_cv*weight_spline, 2));
            }
            else
            {
                hFCfail_w->Fill(truth_nuEnergy, weight_cv*weight_spline);
                hFCfail_w2->Fill(truth_nuEnergy, pow(weight_cv*weight_spline, 2));
            }
          }
          else
          {
            hfail_w->Fill(truth_nuEnergy, weight_cv*weight_spline);
            hfail_w2->Fill(truth_nuEnergy, pow(weight_cv*weight_spline, 2));
            hFCfail_w->Fill(truth_nuEnergy, weight_cv*weight_spline);
            hFCfail_w2->Fill(truth_nuEnergy, pow(weight_cv*weight_spline, 2));
          }
      }
  }
  for(int i=0; i<T_eval2->GetEntries(); i++)
  {
      T_eval2->GetEntry(i);
      T_BDTvars2->GetEntry(i);
      if(match_completeness_energy>0.1*truth_energyInside && truth_nuPdg==14 && truth_isCC==1 && truth_vtxInside==1 && truth_nuEnergy<=4000 && truth_nuEnergy>200)
      {
          if(numu_cc_flag>=0 && numu_score>0.9)
          {
            hpass_w->Fill(truth_nuEnergy, weight_cv*weight_spline);
            hpass_w2->Fill(truth_nuEnergy, pow(weight_cv*weight_spline, 2));
            if(match_isFC==1)
            {
                hFCpass_w->Fill(truth_nuEnergy, weight_cv*weight_spline);
                hFCpass_w2->Fill(truth_nuEnergy, pow(weight_cv*weight_spline, 2));
            }
            else
            {
                hFCfail_w->Fill(truth_nuEnergy, weight_cv*weight_spline);
                hFCfail_w2->Fill(truth_nuEnergy, pow(weight_cv*weight_spline, 2));
            }
          }
          else
          {
            hfail_w->Fill(truth_nuEnergy, weight_cv*weight_spline);
            hfail_w2->Fill(truth_nuEnergy, pow(weight_cv*weight_spline, 2));
            hFCfail_w->Fill(truth_nuEnergy, weight_cv*weight_spline);
            hFCfail_w2->Fill(truth_nuEnergy, pow(weight_cv*weight_spline, 2));
          }
      }
  }

  TH1F* heff = new TH1F("heff", "", 10, xbins_true);
  heff->Add(hpass_w, 1);
  TH1F* htotal = (TH1F*)heff->Clone("htotal");
  htotal->Add(hfail_w, 1);
  heff->Divide(htotal);

  TH1F* hFCeff = new TH1F("hFCeff", "", 10, xbins_true);
  hFCeff->Add(hFCpass_w, 1);
  TH1F* hFCtotal = (TH1F*)hFCeff->Clone("hFCtotal");
  hFCtotal->Add(hFCfail_w, 1);
  hFCeff->Divide(hFCtotal);

  for(int i=1; i<=heff->GetNbinsX(); i++)
  {
      float s = hpass_w->GetBinContent(i);
      float sigma_s = TMath::Sqrt(hpass_w2->GetBinContent(i));
      float b = hfail_w->GetBinContent(i);
      float sigma_b = TMath::Sqrt(hfail_w2->GetBinContent(i));
      heff->SetBinError(i, eff_error(s, sigma_s, b, sigma_b));

      //float N = s + b;
      //heff->SetBinError(i, TMath::Sqrt(s/N*(1-s/N)/N));
      
  }

  for(int i=1; i<=hFCeff->GetNbinsX(); i++)
  {
      float s = hFCpass_w->GetBinContent(i);
      float sigma_s = TMath::Sqrt(hFCpass_w2->GetBinContent(i));
      float b = hFCfail_w->GetBinContent(i);
      float sigma_b = TMath::Sqrt(hFCfail_w2->GetBinContent(i));
      hFCeff->SetBinError(i, eff_error(s, sigma_s, b, sigma_b));
      
      //float N = s + b;
      //hFCeff->SetBinError(i, TMath::Sqrt(s/N*(1-s/N)/N));
  }

  TFile* file_xs = new TFile("merge_xs.root");
  TMatrixD* mat_response = (TMatrixD*)file_xs->Get("mat_R");
  TH2D* hist_response = new TH2D("hist_response", "", 10, xbins_true, 23, 200, 2500);
  for(int j=0; j<10; j++) // true, column
  {
      for(int i=2; i<25; i++)  // reco, FC, row
      {
          hist_response->SetBinContent(j+1, i-1, (*mat_response)(i, j)/xbins_const[j]);
      }
  }


  //// pad3 efficiency plot and pad1 inset response matrix 


  CovMatrix cov;

  // get data histograms ...
  // filetype, period, outfilename, external pot, fileno
  std::map<TString, std::tuple<int, int, TString, float, int, double, int> > map_inputfile_info = cov.get_map_inputfile_info();

  TFile *temp_file;
  TH1F *htemp;
  TTree *T;
  Double_t pot;
  std::map<TString, std::pair<TH1F*, double> > map_name_histogram;

  // data POT ...
  std::map<int, double> map_data_period_pot;
  //  std::vector<TH1F*> temp_histograms;
  
  // open all the histograms ...
  for (auto it = map_inputfile_info.begin(); it != map_inputfile_info.end(); it++){
    TString input_filename = it->first;
    int filetype = std::get<0>(it->second);
    int period = std::get<1>(it->second);
    TString out_filename = std::get<2>(it->second);
    int file_no = std::get<4>(it->second);
    temp_file = new TFile(out_filename);
    T = (TTree*)temp_file->Get("T");
    T->SetBranchAddress("pot",&pot);
    T->GetEntry(0);
    

    if (filetype==5 || filetype == 15){
      map_data_period_pot[period] = pot;
    }
    
    // std::cout << pot << std::endl;
    
    std::vector< std::tuple<TString,  int, float, float, TString, TString, TString, TString > > all_histo_infos;
    std::vector< std::tuple<TString,  int, float, float, TString, TString, TString, TString > > histo_infos = cov.get_histograms(input_filename,0);
    std::copy(histo_infos.begin(), histo_infos.end(), std::back_inserter(all_histo_infos));
    std::vector< std::tuple<TString,  int, float, float, TString, TString, TString, TString > > histo_infos_err2 = cov.get_histograms(input_filename,1);
    std::copy(histo_infos_err2.begin(), histo_infos_err2.end(), std::back_inserter(all_histo_infos));
    std::vector< std::tuple<TString,  int, float, float, TString, TString, TString, TString > > histo_infos_cros = cov.get_histograms(input_filename,2);
    std::copy(histo_infos_cros.begin(), histo_infos_cros.end(), std::back_inserter(all_histo_infos));
    
    for (size_t i=0;i!=all_histo_infos.size();i++){
      htemp = (TH1F*)temp_file->Get(std::get<0>(all_histo_infos.at(i)));
      //      std::cout << std::get<0>(all_histo_infos.at(i)) << " " << htemp->GetSum() << std::endl;
      //      temp_histograms.push_back(htemp);
      map_name_histogram[std::get<0>(all_histo_infos.at(i))] = std::make_pair(htemp, pot);
    }
  }

  
 
  
  // create histograms for data, create histograms for predictions
  // obsch --> histograms (data, prediction, prediction_err2
  std::map<int, std::vector<TH1F*> > map_obsch_histos;
  // obsch --> break down histograms (truth label: add_cut) prediction 
  std::map<int, std::vector<TH1F*> > map_obsch_subhistos;
  // Bayesian error needed ...
  // obsch --> bin with overflow bin --> vector of all channels (merge certain channels) --> mean and err2 
  std::map<int, std::vector< std::vector< std::tuple<double, double, double, int, double> > > > map_obsch_bayes;
  std::map<int, std::vector< std::vector< std::tuple<double, double, double, int, double> > > > map_obsch_infos;
  
  for (auto it = map_inputfile_info.begin(); it != map_inputfile_info.end(); it++){
    TString input_filename = it->first;
    int filetype = std::get<0>(it->second);
    int period = std::get<1>(it->second);
    TString out_filename = std::get<2>(it->second);
    int file_no = std::get<4>(it->second);
    if (filetype == 5 || filetype == 15){
      // name, nbin, lowlimit, highlimit, variable, channel cut, additional cut, weight
      std::vector< std::tuple<TString,  int, float, float, TString, TString, TString, TString > > histo_infos = cov.get_histograms(input_filename,0);
      for (auto it1 = histo_infos.begin(); it1 != histo_infos.end(); it1++){
	int obsch = cov.get_obsch_name(std::get<5>(*it1));	
	htemp = map_name_histogram[std::get<0>(*it1)].first;

	std::vector<TH1F*> vec_histos;
	
	TH1F *hdata = (TH1F*)htemp->Clone(Form("data_%d",obsch));
	hdata->Reset();
	TH1F *hpred = (TH1F*)htemp->Clone(Form("pred_%d",obsch));
	hpred->Reset();
	TH1F *hpred_err2 = (TH1F*)htemp->Clone(Form("pred_err2_%d",obsch));
	hpred_err2->Reset();

	vec_histos.push_back(hdata);
	vec_histos.push_back(hpred);
	vec_histos.push_back(hpred_err2);

	// get histograms ...
	map_obsch_histos[obsch] = vec_histos;
	//map_obsch_bayes[obsch].resize(htemp->GetNbinsX()+1);
	// for (Int_t i=0;i!=htemp->GetNbinsX()+1;i++){
	//   std::vector< std::tuple<double, double, double> > temp;
	//   map_obsch_bayes[obsch].push_back(temp);
	// }
	
	
	//std::cout << std::get<5>(*it1) << " " << obsch << " " << htemp->GetSum() << std::endl;
      }
      
      // break;
    }
  }
  
  // get data histograms ...
  cov.fill_data_histograms(run, map_obsch_histos, map_name_histogram);
  
  // get predictions and its uncertainties ...,
  cov.fill_pred_histograms(run, map_obsch_histos, map_obsch_bayes, map_obsch_infos, map_name_histogram, lee_strength, map_data_period_pot, flag_breakdown, map_obsch_subhistos);

    /* for (auto it = map_obsch_subhistos.begin(); it!= map_obsch_subhistos.end(); it++){ */
    /*         for(size_t i=0; i<it->second.size(); i++){ */
    /*             std::cout<<"DEBUG2: "<<it->first<<": "<<map_obsch_subhistos[it->first].at(i)->GetName()<<std::endl; */
    /*         } */
    /* } */ 

  // get Bayesian errrors ...

  if (flag_err==2){
    std::cout << lee_strength << " " << run << std::endl;
    
    for (auto it = map_obsch_histos.begin(); it!= map_obsch_histos.end(); it++){
      TH1F *h1 = it->second.at(1);
      TH1F *h2 = it->second.at(2);
      int obsch = it->first;

      std::vector<std::vector< std::tuple<double, double, double, int, double> > >  bayes_inputs = map_obsch_bayes[obsch];

      //std::cout << obsch << " " << bayes_inputs.size() << " " << bayes_inputs.at(0).size() << " " << h1->GetNbinsX() << std::endl;

      //if (obsch !=1) continue;
     double temp_sumtotalcov=0; 
      for (int i=0;i!=h1->GetNbinsX()+1;i++){
	Bayes bayes;
	//	if (i!=0) continue;
	//double temp = 0, temp1=0;
	double zero_weight = 0, zero_count = 0;
	double nonzero_meas = 0, nonzero_sigma2 = 0, nonzero_weight = 0;
	for (auto it1 = bayes_inputs.begin(); it1!=bayes_inputs.end(); it1++){
	  // bayes.add_meas_component(std::get<0>((*it1).at(i)), std::get<1>((*it1).at(i)), std::get<2>((*it1).at(i)));
	  // temp += std::get<0>((*it1).at(i));
	  // temp1 += std::get<1>((*it1).at(i));
	  //std::cout << i << " " << std::get<0>((*it1).at(i)) << " " << std::get<1>((*it1).at(i)) << " " << std::get<2>((*it1).at(i)) << " " << std::endl;

	  // approximation for zeros ...
	  if (std::get<0>((*it1).at(i)) == 0){
	    zero_weight += pow(std::get<2>((*it1).at(i)),2);
	    zero_count +=std::get<2>((*it1).at(i));
	  }else{
	    nonzero_meas += std::get<0>((*it1).at(i));
	    nonzero_sigma2 += std::get<1>((*it1).at(i));
	    nonzero_weight += std::get<2>((*it1).at(i));
	  }
	}
	if (zero_count != 0)	bayes.add_meas_component(0,0,zero_weight/zero_count);
	if (nonzero_meas != 0)  bayes.add_meas_component(nonzero_meas, nonzero_sigma2, nonzero_weight);
	
	bayes.do_convolution();
	
	double cov = bayes.get_covariance();
	//double cov1 = bayes.get_covariance_mc();
	std::cout << obsch << " " << i << " "	  << h1->GetBinContent(i+1) << " " << cov  << " "  << h2->GetBinContent(i+1) << std::endl;
	// std::cout << temp << " " << temp1 << " " << h1->GetBinContent(i+1) << " " << h2->GetBinContent(i+1) << std::endl;
        if(isnan(cov) || isinf(cov)) {
            std::cout<<"Abnormal Bayesian variance.\n";
            cov = bayes.get_covariance_mc();
            //cov = h1->SetBinError(i+1, h1->GetBinError(i));
        }
	h1->SetBinError(i+1,sqrt(cov));
    if(i!=h1->GetNbinsX()) temp_sumtotalcov += cov;
      }
      sumtotalcov[obsch] = temp_sumtotalcov;
      // obsch --> bin with overflow bin --> vector of all channels (merge certain channels) --> mean and err2 
      //std::map<int, std::vector< std::vector< std::tuple<double, double, double> > > > map_obsch_bayes;
    }
  }

  if (flag_err==3){
    std::cout <<"Total uncertainty from external covariance matrix: "<< cov_inputfile << std::endl;
    TFile* f_cov = new TFile(cov_inputfile, "READ");
    int flag_syst_flux_Xs = 0;
    int flag_syst_detector = 0;
    int flag_syst_additional = 0;
    int flag_syst_mc_stat = 0;
    double cov_lee_strength = 0;
    std::vector<double> *vc_val_GOF = new std::vector<double>;
    std::vector<int> *vc_val_GOF_NDF = new std::vector<int>;
    TTree* t_covconfig = (TTree*)f_cov->Get("tree");
    t_covconfig->SetBranchAddress("flag_syst_flux_Xs", &flag_syst_flux_Xs);
    t_covconfig->SetBranchAddress("flag_syst_detector", &flag_syst_detector);
    t_covconfig->SetBranchAddress("flag_syst_additional", &flag_syst_additional);
    t_covconfig->SetBranchAddress("flag_syst_mc_stat", &flag_syst_mc_stat);
    t_covconfig->SetBranchAddress("user_Lee_strength_for_output_covariance_matrix", &cov_lee_strength);
    t_covconfig->SetBranchAddress("vc_val_GOF", &vc_val_GOF);
    t_covconfig->SetBranchAddress("vc_val_GOF_NDF", &vc_val_GOF_NDF);
    t_covconfig->GetEntry(0);

    // fill GOF map
    for(size_t vv=0; vv<vc_val_GOF->size(); vv++)
    {
        GOF[vv] = std::make_pair(vc_val_GOF->at(vv), vc_val_GOF_NDF->at(vv));
    }

    // absolute cov matrix
    TMatrixD* matrix_absolute_cov = (TMatrixD*)f_cov->Get("matrix_absolute_cov_newworld");
    TMatrixD* matrix_absolute_detector_cov = (TMatrixD*)f_cov->Get("matrix_absolute_detector_cov_newworld");

    std::cout <<"User: "<< "chosen LEE strength: "<< lee_strength << " run option: " << run << std::endl;
    std::cout <<"Cov matrix config: \n"
                <<"\t syst_flux_Xs: "<< flag_syst_flux_Xs << std::endl  
                <<"\t syst_detector: "<< flag_syst_detector << std::endl  
                <<"\t syst_additional: "<< flag_syst_additional << std::endl  
                <<"\t syst_mc_stat: "<< flag_syst_mc_stat << std::endl  
                <<"\t LEE_strength: "<< cov_lee_strength << std::endl;
    if( abs(lee_strength-cov_lee_strength)>1e-6 ) {
        std::cout<<"ERROR: plot_hist -l option is inconsistent with external cov matrix LEE strength!\n";
        return 1;
    }
    
    // construct a map from (obsch, bin) to cov index
    std::map< std::pair<int, int>, int> obsch_bin_index;
    int index = 0;
    for(size_t i=0; i<map_obsch_histos.size(); i++)
    {
        // + overflow bin
        for(int j=1; j<=map_obsch_histos[i+1].at(1)->GetNbinsX()+1; j++)
        {
            index += 1;
            // cov index starts from 0
            obsch_bin_index.insert({{i+1, j}, index-1});
        }
    }

    TMatrixD* matrix_pred = (TMatrixD*)f_cov->Get("matrix_pred_newworld");
    TMatrixD* matrix_data = (TMatrixD*)f_cov->Get("matrix_data_newworld");
    for (auto it = map_obsch_histos.begin(); it!= map_obsch_histos.end(); it++){
      TH1F *h1 = it->second.at(1); // error --> total uncertainty
      TH1F *h2 = it->second.at(2); // bonus: error --> detector systematic uncertainty
      int obsch = it->first;
  
      TH1F* htemp_data = (TH1F*)h1->Clone("htemp_data");
      TH1F* htemp_pred = (TH1F*)h1->Clone("htemp_pred");
      htemp_data->Reset();
      htemp_pred->Reset();
      double temp_sumtotalcov=0;
      for (int i=0;i!=h1->GetNbinsX()+1;i++){
        int index = obsch_bin_index.find(std::make_pair(obsch, i+1))->second;
        double total_uncertainty = (*matrix_absolute_cov)(index, index); // only diagonal term
        //summation of total cov
        if(i!=h1->GetNbinsX()){ //no overflow bin in this calculation
        for(int j=0; j!=h1->GetNbinsX();j++){
            int jndex = obsch_bin_index.find(std::make_pair(obsch, j+1))->second;
            temp_sumtotalcov += (*matrix_absolute_cov)(index, jndex);
        }
        }
        //
        double detector_uncertainty = (*matrix_absolute_detector_cov)(index, index); // only diagonal term
        std::cout << obsch << " " << i << " "	  << h1->GetBinContent(i+1) << " " << total_uncertainty << " " <<sqrt(total_uncertainty)/h1->GetBinContent(i+1) << " "  << h2->GetBinContent(i+1) << " " << detector_uncertainty << std::endl;
	    h1->SetBinError(i+1,sqrt(total_uncertainty));
        h2->SetBinError(i+1,sqrt(detector_uncertainty));

        if(flag_check == 1){
            htemp_data->SetBinContent(i+1, (*matrix_data)(0, index));
            htemp_pred->SetBinContent(i+1, (*matrix_pred)(0, index));
        }
      }

      sumtotalcov[obsch] = temp_sumtotalcov;
      if(flag_check == 1){
        it->second.push_back(htemp_data);
        it->second.push_back(htemp_pred);
      }
    }
  }

  if (flag_display == 1){
    // plotting ...
    TApplication theApp("theApp",&argc,argv);
    theApp.SetReturnFromRun(true);
    //gStyle->SetOptStat(0);

    gROOT->ProcessLine(".x DrawOption.cc");
    gStyle->SetLegendBorderSize(0);
    gStyle->SetLegendFillColor(0);
    gStyle->SetLegendFont(132);
    gStyle->SetLegendTextSize(0.04);

    int nchannels = map_obsch_subhistos.size();
    TCanvas *canvas[nchannels];
    TGraphAsymmErrors *gr[nchannels];
    TGraphAsymmErrors *gratio_mc[nchannels];
    TGraphAsymmErrors *gratio_data[nchannels];
    THStack *hstack[nchannels]; 
    TLegend *legend[nchannels]; 
    TLegend *legend2[nchannels]; 
    TLegend *legend3[nchannels]; 
    for(auto it = map_obsch_subhistos.begin(); it!= map_obsch_subhistos.end(); it++){
        int obschannel = it->first;
        std::cout<<"Channel: "<<obschannel<<std::endl;
        canvas[obschannel-1] = new TCanvas(Form("canvas%d", obschannel), Form("channel%d", obschannel), 800, 1000);
        TPad *pad1 = new TPad("pad1", "", 0.01,0.5,0.99,0.99,0,0,0);
        TPad *pad2 = new TPad("pad2", "", 0.01,0.25,0.99,0.5,0,0,0);
        TPad *pad3 = new TPad("pad3", "", 0.01,0.01,0.99,0.25,0,0,0);
        pad1->SetBottomMargin(0.02);
        pad1->SetLeftMargin(0.12);
        pad1->SetRightMargin(0.05);
        pad2->SetTopMargin(0.04);
        pad2->SetLeftMargin(0.12);
        pad2->SetRightMargin(0.05);
        pad2->SetBottomMargin(0.23);
        pad3->SetTopMargin(0.02);
        pad3->SetLeftMargin(0.12);
        pad3->SetRightMargin(0.05);
        pad3->SetBottomMargin(0.21);
        pad1->Draw();
        pad2->Draw();
        pad3->Draw();
        TPad *pad1_2 = new TPad("pad1_2", "", 0.5,0.59,0.92,0.85,0,0,0);
        pad1_2->SetFrameFillStyle(4000);
        pad1_2->SetFillStyle(4000);
        pad1_2->SetTopMargin(0.1);
        pad1_2->SetBottomMargin(0.13);
        pad1_2->SetLeftMargin(0.18);
        pad1_2->SetRightMargin(0.12);
        pad1_2->Draw(); // inset 
        hstack[obschannel-1] = new THStack(Form("hs%d", obschannel),"");
        legend[obschannel-1] = new TLegend(0.2, 0.7, 0.85, 0.93);
        TH1F* hdata = (TH1F*)map_obsch_histos[obschannel].at(0)->Clone("hdata");
        TH1F* hXsecCosmic = (TH1F*)hdata->Clone("hXsecCosmic");
        TH1F* hXsecNumuCCinFV = (TH1F*)hdata->Clone("hXsecNumuCCinFV");
        TH1F* hXsecNC = (TH1F*)hdata->Clone("hXsecNC");
        TH1F* hXsecBkgCC = (TH1F*)hdata->Clone("hXsecBkgCC");
        TH1F* hext = (TH1F*)hdata->Clone("hext");
        hXsecCosmic->Reset();
        hXsecNumuCCinFV->Reset();
        hXsecNC->Reset();
        hXsecBkgCC->Reset();
        hext->Reset();
        //hack
        double scalePOT = 1.0; // overall POT scaling
        //scalePOT = 5.0/5.327;
        //end
        for(size_t i=0; i<it->second.size(); i++){
            TH1F* htemp = map_obsch_subhistos[obschannel].at(i);
            std::string histname = htemp->GetName();
            std::istringstream sss(histname);
            htemp->Scale(scalePOT);
            for(std::string line; std::getline(sss, line, '_');){
                if(line == "XsecCosmic") {
                    std::cout<<"XsecCosmic"<<" "<<histname<<std::endl;
                    hXsecCosmic->Add(htemp);
                    break;
                }
                if(line == "XsecNumuCCinFV") {
                    std::cout<<"XsecNumuCCinFV"<<" "<<histname<<std::endl;
                    hXsecNumuCCinFV->Add(htemp);
                    break;
                }
                if(line == "XsecNC") {
                    std::cout<<"XsecNC"<<" "<<histname<<std::endl;
                    hXsecNC->Add(htemp);
                    break;
                }
                if(line == "XsecBkgCC") {
                    std::cout<<"XsecBkgCC"<<" "<<histname<<std::endl;
                    hXsecBkgCC->Add(htemp);
                    break;
                }
                if(line == "ext") {
                    std::cout<<"ext"<<" "<<histname<<std::endl;
                    hext->Add(htemp);
                    break;
                }
            }
        }
        pad1->cd();

        gStyle->SetOptTitle(kFALSE);
        float datapot = 0;
        float datapot_numi = 0;
        if(run == 0){
            for(auto it=map_data_period_pot.begin(); it!=map_data_period_pot.end(); it++)
            {
                if(it->first<10) datapot += it->second;
                if(it->first>10) datapot_numi += it->second;
            }
        }else datapot = map_data_period_pot[run];

        gr[obschannel-1] = new TGraphAsymmErrors();
        legend[obschannel-1]->SetNColumns(3);
        // numi channels
        //legend[obschannel-1]->AddEntry((TObject*)0, Form("Scaled POT: %.3e", datapot*scalePOT), "");
        legend[obschannel-1]->AddEntry(gr[obschannel-1], Form("BNB data, %.1f", hdata->Integral()*scalePOT), "lp");
        //legend[obschannel-1]->AddEntry(gr[obschannel-1], Form("Scaled BNB data, %.1f", hdata->Integral()*scalePOT), "lp");

        TH1F* hmc = (TH1F*)map_obsch_histos[obschannel].at(1)->Clone("hmc");
        TH1F* hmc2 = (TH1F*)map_obsch_histos[obschannel].at(2)->Clone("hmc2");
        TH1F* hmcerror = (TH1F*)hmc->Clone("hmcerror");
        legend[obschannel-1]->AddEntry(hmcerror, "Pred. uncertainty", "lf");

        if(flag_truthlabel==0){
        // truth labels start
        hstack[obschannel-1]->Add(hext); 
        legend[obschannel-1]->AddEntry(hext, Form("Beam-off, %.1f", hext->Integral()), "F"); 
        hext->SetFillStyle(3004);
        hext->SetFillColorAlpha(kOrange+3, 0.5);
        hext->SetLineColor(kOrange+3);
        hext->SetLineWidth(1);
        
        hstack[obschannel-1]->Add(hXsecCosmic); 
        legend[obschannel-1]->AddEntry(hXsecCosmic, Form("Cosmic, %.1f", hXsecCosmic->Integral()), "F"); 
        hXsecCosmic->SetFillStyle(3224);
        hXsecCosmic->SetFillColorAlpha(kRed+2, 0.5);
        hXsecCosmic->SetLineColor(kRed+2);
        hXsecCosmic->SetLineWidth(1);

        hstack[obschannel-1]->Add(hXsecNC); 
        legend[obschannel-1]->AddEntry(hXsecNC, Form("NC,  %.1f", hXsecNC->Integral()), "F"); 
        hXsecNC->SetFillStyle(1001);
        hXsecNC->SetFillColorAlpha(kOrange+1, 0.5);
        hXsecNC->SetLineColor(38);
        hXsecNC->SetLineWidth(1);

        hstack[obschannel-1]->Add(hXsecBkgCC); 
        legend[obschannel-1]->AddEntry(hXsecBkgCC, Form("Other CC, %.1f", hXsecBkgCC->Integral()), "F"); 
        hXsecBkgCC->SetFillStyle(1001);
        hXsecBkgCC->SetFillColorAlpha(kGreen+1, 0.5);
        hXsecBkgCC->SetLineColor(30);
        hXsecBkgCC->SetLineWidth(1);

        hstack[obschannel-1]->Add(hXsecNumuCCinFV); 
        legend[obschannel-1]->AddEntry(hXsecNumuCCinFV, Form("#nu_{#mu} CC in FV, %.1f", hXsecNumuCCinFV->Integral()), "F"); 
        hXsecNumuCCinFV->SetFillStyle(1001);
        hXsecNumuCCinFV->SetFillColorAlpha(kAzure+6, 0.5);
        hXsecNumuCCinFV->SetLineColor(kAzure+6);
        hXsecNumuCCinFV->SetLineWidth(1);
        // truth labels end
        }

        hmc->Sumw2();
        hmc->Scale(scalePOT);
        hmc->Draw("hist");
        hmc->GetYaxis()->SetTitle("Event counts");
        hmc->GetYaxis()->SetTitleSize(0.06);
        hmc->GetYaxis()->SetTitleOffset(0.80);
        hmc->GetYaxis()->SetLabelSize(0.05);
        hmc->GetXaxis()->SetLabelColor(kWhite);
        //if(obschannel==9) hmc->GetXaxis()->SetRangeUser(0.5,1);
        float mcymax = hmc->GetBinContent(hmc->GetMaximumBin())*scalePOT;
        float dataymax = hdata->GetBinContent(hdata->GetMaximumBin())*scalePOT;
        if(dataymax>mcymax) mcymax = dataymax;
        hmc->SetMaximum(2.0*mcymax);
        hmc->GetYaxis()->SetRangeUser(-0.02*mcymax, 1.6*mcymax);
        hmc->SetLineColor(kBlack);
        hmc->SetLineWidth(5);


        hstack[obschannel-1]->Draw("hist same");
        hmcerror->Sumw2();
        hmcerror->Scale(scalePOT);
        hmcerror->Draw("same E2");
        hmcerror->SetFillColor(kGray+2);
        hmcerror->SetFillStyle(3002);
        hmcerror->SetLineWidth(0);
        hmcerror->SetLineColor(12);
        hmcerror->SetMarkerColor(0);
        hmcerror->SetMarkerSize(0);
        
        gratio_mc[obschannel-1] = new TGraphAsymmErrors();
        gratio_data[obschannel-1] = new TGraphAsymmErrors();
        float maxratio = 1.5;
        for(int i=0; i<hdata->GetNbinsX(); i++)
        {
            double x = hdata->GetBinCenter(i+1);
            double x_err = hdata->GetBinWidth(1)*0.5;
            double y = hdata->GetBinContent(i+1);
            double y_err = hdata->GetBinError(i+1);
            auto bayesError = cov.get_bayes_errors(y);
            double bayesError_low = bayesError.first;
            double bayesError_up = bayesError.second;
            double ymc = hmc->GetBinContent(i+1);
            double ymc_err = hmc->GetBinError(i+1);
            double ymc2_err = hmc2->GetBinError(i+1);
            gr[obschannel-1]->SetPoint(i,x,y*scalePOT);
            gratio_mc[obschannel-1]->SetPoint(i,x,1);
            if(ymc!=0){ 
                gratio_data[obschannel-1]->SetPoint(i,x,y/ymc); 
                if(maxratio<y/ymc) maxratio = y/ymc;
                gratio_mc[obschannel-1]->SetPointError(i, x_err, x_err, ymc_err/ymc, ymc_err/ymc);
            }
            else { 
                gratio_data[obschannel-1]->SetPoint(i, x, 10); // invalid value 
                gratio_mc[obschannel-1]->SetPointError(i, x_err, x_err, 0, 0);
            }
            if(flag_err==2 || flag_err==3){ //update data point errors
                gr[obschannel-1]->SetPointError(i, x_err, x_err, bayesError_low*scalePOT, bayesError_up*scalePOT);
                if(ymc!=0) gratio_data[obschannel-1]->SetPointError(i, x_err, x_err, bayesError_low/ymc, bayesError_up/ymc);
                else gratio_data[obschannel-1]->SetPointError(i, x_err, x_err, 0, 0);
            }
            if(flag_err==1){
                gr[obschannel-1]->SetPointError(i, x_err, x_err, y_err*scalePOT, y_err*scalePOT);
                if(ymc!=0) gratio_data[obschannel-1]->SetPointError(i, x_err, x_err, y_err/ymc, y_err/ymc);
                else gratio_data[obschannel-1]->SetPointError(i, x_err, x_err, 0, 0);
            }
        }
        gr[obschannel-1]->Draw("P same"); 
        gr[obschannel-1]->SetLineWidth(2);
        gr[obschannel-1]->SetMarkerStyle(20);
        gr[obschannel-1]->SetMarkerSize(1.5);
        gr[obschannel-1]->SetLineColor(kBlack);
        if(flag_check == 1){
            TH1F* hcheck_data = (TH1F*)map_obsch_histos[obschannel].at(3);
            TH1F* hcheck_pred = (TH1F*)map_obsch_histos[obschannel].at(4);
            hcheck_data->Draw("hist same");
            hcheck_data->SetLineColor(kBlack);
            hcheck_pred->Draw("hist same");
            hcheck_pred->SetLineColor(kRed);
        }


        //legend[obschannel-1]->SetFillStyle(0);
        double relerr_data = 1./TMath::Sqrt(hdata->Integral());
        double relerr_pred = TMath::Sqrt(sumtotalcov[obschannel])/hmc->Integral();
        double data_pred_ratio = hdata->Integral()/hmc->Integral();
        legend[obschannel-1]->SetHeader(Form("#SigmaDATA/#Sigma(MC+EXT)=%.2f#pm%.2f(data err)#pm%.2f(pred err)", data_pred_ratio, relerr_data*data_pred_ratio, relerr_pred*data_pred_ratio), "C");
        legend[obschannel-1]->AddEntry((TObject*)0, Form("%.3e POT", datapot*scalePOT), "");
        legend[obschannel-1]->SetTextSize(0.04);
        legend[obschannel-1]->SetFillStyle(0);
        legend[obschannel-1]->Draw();
        //pad1->Modified();

        pad1_2->cd();
        pad1_2->SetLogz();
        hist_response->GetZaxis()->SetRangeUser(0.001, 0.4);
        hist_response->Draw("colz");
        hist_response->GetYaxis()->SetNdivisions(505);
        hist_response->GetXaxis()->SetNdivisions(505);
        hist_response->GetXaxis()->SetTitle("True #nu energy [MeV]");
        hist_response->GetYaxis()->SetTitle("Reco #nu energy [MeV]");
        hist_response->GetYaxis()->SetLabelSize(0.06);
        hist_response->GetYaxis()->SetTitleSize(0.07);
        hist_response->GetYaxis()->SetTitleOffset(1.20);
        hist_response->GetXaxis()->SetLabelSize(0.06);
        hist_response->GetXaxis()->SetTitleSize(0.07);
        hist_response->GetXaxis()->SetTitleOffset(0.9);
        hist_response->GetZaxis()->SetLabelSize(0.05);

        pad2->cd();
        TH1F* hist = (TH1F*)hdata->Clone("hist");
        hist->Reset();
        hist->Scale(scalePOT);
        hist->GetYaxis()->SetRangeUser(0,int(1.5*maxratio)<2?int(1.5*maxratio):2);
        hist->GetYaxis()->SetNdivisions(505);
        hist->GetXaxis()->SetRangeUser(hmc->GetXaxis()->GetXmin(),hmc->GetXaxis()->GetXmax());
        hist->Draw("axis");
        hist->GetYaxis()->SetTitle("Data/Pred");
        hist->GetYaxis()->SetTitleSize(0.12);
        hist->GetYaxis()->SetLabelSize(0.10);
        hist->GetYaxis()->SetTitleOffset(0.40);
        hist->GetXaxis()->SetTitle("Reco neutrino energy [MeV]");
        hist->GetXaxis()->SetTitleSize(0.10);
        hist->GetXaxis()->SetLabelSize(0.10);
        hist->GetXaxis()->SetTickSize(0.06);
        hist->GetXaxis()->SetTitleOffset(0.93);

        gratio_mc[obschannel-1]->Draw("2 same");
        gratio_mc[obschannel-1]->SetFillColor(kRed-10);
        gratio_mc[obschannel-1]->GetYaxis()->SetRangeUser(0,int(1.5*maxratio)<2?int(1.5*maxratio):2);
        gratio_mc[obschannel-1]->GetXaxis()->SetRangeUser(hmc->GetXaxis()->GetXmin(),hmc->GetXaxis()->GetXmax());
        gratio_mc[obschannel-1]->GetXaxis()->SetTitleSize(0.1);
        gratio_mc[obschannel-1]->GetXaxis()->SetLabelSize(0.1);
        gratio_mc[obschannel-1]->GetYaxis()->SetTitleSize(0.1);
        gratio_mc[obschannel-1]->GetYaxis()->SetTitleOffset(0.35);
        gratio_mc[obschannel-1]->GetYaxis()->SetLabelSize(0.1);
        gratio_data[obschannel-1]->Draw("P same");
        gratio_data[obschannel-1]->GetYaxis()->SetRangeUser(0,int(1.5*maxratio)<2?int(1.5*maxratio):2);
        gratio_data[obschannel-1]->GetXaxis()->SetRangeUser(hmc->GetXaxis()->GetXmin(),hmc->GetXaxis()->GetXmax());
        gratio_data[obschannel-1]->SetLineWidth(2);
        gratio_data[obschannel-1]->SetMarkerStyle(20);
        gratio_data[obschannel-1]->SetMarkerSize(1.5);
        gratio_data[obschannel-1]->SetLineColor(kBlack);
       
        
        hist->Draw("axis same");

        TLine* line; 
        line = new TLine(hmc->GetXaxis()->GetXmin(),1,hmc->GetXaxis()->GetXmax(),1);
        //if(obschannel==5 || obschannel==6) line = new TLine(0,1,1200,1);
        line->Draw();
        line->SetLineWidth(2);
        line->SetLineStyle(kDashed);
        legend2[obschannel-1] = new TLegend(0.2, 0.7, 0.8, 0.95);
        legend2[obschannel-1]->SetNColumns(2);
        legend2[obschannel-1]->AddEntry((TObject*)0, Form("#chi^{2}/ndf=%.2f/%d", GOF[obschannel-1].first, GOF[obschannel-1].second), "");
        if(flag_err==1) legend2[obschannel-1]->AddEntry(gratio_mc[obschannel-1],"Pred stat. uncertainty", "F");
        if(flag_err==2) legend2[obschannel-1]->AddEntry(gratio_mc[obschannel-1],"Pred stat. uncertainty (Bayesian)", "F");
        if(flag_err==3) legend2[obschannel-1]->AddEntry(gratio_mc[obschannel-1],"Pred uncertainty", "F");
        //legend2[obschannel-1]->AddEntry(gratio_data[obschannel-1],"Data with stat. uncertainty", "lp");
        legend2[obschannel-1]->SetTextSize(0.10);
        legend2[obschannel-1]->SetFillStyle(0);
        legend2[obschannel-1]->Draw();
        pad2->Modified();

        pad3->cd();
        legend3[obschannel-1] = new TLegend(0.5, 0.35, 0.9, 0.75);
        legend3[obschannel-1]->SetNColumns(2);
        legend3[obschannel-1]->SetTextSize(0.10);
        heff->Draw("P");
        heff->SetLineWidth(2);
        heff->GetYaxis()->SetTitle("Efficiency");
        heff->GetYaxis()->SetRangeUser(0.0,1.0);
        heff->GetYaxis()->SetNdivisions(505);
        heff->GetYaxis()->SetTitleSize(0.13);
        heff->GetYaxis()->SetLabelSize(0.11);
        heff->GetYaxis()->SetTitleOffset(0.35);
        heff->GetXaxis()->SetTitle("True neutrino energy [MeV]");
        heff->GetXaxis()->SetTitleSize(0.11);
        heff->GetXaxis()->SetLabelSize(0.11);
        heff->GetXaxis()->SetTickSize(0.06);
        heff->GetXaxis()->SetTitleOffset(0.9);
        heff->SetLineColor(kRed);

        hFCeff->Draw("P same");
        hFCeff->SetLineWidth(2);
        hFCeff->SetLineColor(kBlue);
        
        legend3[obschannel-1]->AddEntry(heff, "FC+PC", "lp");
        legend3[obschannel-1]->AddEntry(hFCeff, "FC only", "lp");
        legend3[obschannel-1]->Draw();



        canvas[obschannel-1]->Print(Form("canvas%d.pdf", obschannel));
        
        if(obschannel==1) canvas[obschannel-1]->Print("selection.pdf(");
        else if(obschannel==nchannels) canvas[obschannel-1]->Print("selection.pdf)");
        else canvas[obschannel-1]->Print("selection.pdf");

    } 
    theApp.Run();
  }

}
