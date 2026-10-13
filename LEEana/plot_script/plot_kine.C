#include "/home/xqian/wire-cell/wcp-uboone-bdt/inc/WCPLEEANA/kine.h"

using namespace LEEana;

void plot_kine(){
  //TFile *file1 = new TFile("processed_checkout_rootfiles/checkout_data_bnb_run1_far_sideband_2021_02_03.root");
  //TFile *file1 = new TFile("processed_checkout_rootfiles/checkout_data_bnb_run1_5e19.root");
  //TFile *file1 = new TFile("processed_checkout_rootfiles/checkout_data_bnb_run1_5e19.root");
  TFile *file1 = new TFile("far_sideband_files/checkout_data_bnb_run1_far_sideband_2021_02_03.root");
  TTree *T_BDTvars = (TTree*)file1->Get("wcpselection/T_BDTvars");
  
  TTree *T_eval = (TTree*)file1->Get("wcpselection/T_eval");
  TTree *T_pot = (TTree*)file1->Get("wcpselection/T_pot");
  TTree *T_PFeval = (TTree*)file1->Get("wcpselection/T_PFeval");
  TTree *T_KINEvars = (TTree*)file1->Get("wcpselection/T_KINEvars");

  KineInfo kine;
  set_tree_address(T_KINEvars, kine);
  
  kine.kine_energy_particle = new std::vector<float>;
  kine.kine_energy_info = new std::vector<int>;
  kine.kine_particle_type = new std::vector<int>;
  kine.kine_energy_included = new std::vector<int>;

  double em_charge_scale = 0.95;
  
  for (int i=0;i!=100;i++){
    T_KINEvars->GetEntry(i);
    Double_t reco_Enu_corr = 0;
    if (kine.kine_reco_Enu > 0){
      for ( size_t j=0;j!= kine.kine_energy_particle->size();j++){
	
	if (kine.kine_energy_info->at(j) == 2 && kine.kine_particle_type->at(j) == 11){
	  reco_Enu_corr +=  kine.kine_energy_particle->at(j) * em_charge_scale;
	}else{
	  reco_Enu_corr +=  kine.kine_energy_particle->at(j);
	}
	
	//	std::cout << "p: " << kine.kine_energy_particle->at(j) << " " << kine.kine_energy_info->at(j) << " " << kine.kine_particle_type->at(j) << " " << kine.kine_energy_included->at(j) << std::endl; 
      }
      reco_Enu_corr += kine.kine_reco_add_energy;
      //      std::cout <<  kine.kine_energy_particle->size() << " " << kine.kine_energy_info->size() << " " << kine.kine_particle_type->size() << " " << kine.kine_energy_included->size() << std::endl;
      
      std::cout << i << " " << kine.kine_reco_Enu << " " << reco_Enu_corr << " " << kine.kine_reco_add_energy << std::endl;
    }
  }
}
