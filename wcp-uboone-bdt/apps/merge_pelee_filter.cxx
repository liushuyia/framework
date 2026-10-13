// cz: code modified from tutorials/tmva/TMVAClassification.C

#include <cstdlib>
#include <iostream>
#include <map>
#include <string>
#include <set>

#include "TChain.h"
#include "TFile.h"
#include "TTree.h"
#include "TString.h"
#include "TObjString.h"
#include "TSystem.h"
#include "TROOT.h"
#include "TMath.h"

#include "WCPLEEANA/tagger.h"

#include "TMVA/Factory.h"
#include "TMVA/DataLoader.h"
#include "TMVA/Tools.h"
#include "TMVA/TMVAGui.h"
#include "TMVA/Reader.h"


#include "WCPLEEANA/eval.h"

using namespace std;
using namespace LEEana;

#include "WCPLEEANA/bdt.h"
#include "WCPLEEANA/pot.h"
#include "WCPLEEANA/pfeval.h"
#include "WCPLEEANA/kine.h"
#include "WCPLEEANA/cuts.h"

int main( int argc, char** argv )
{
  if (argc < 5) {
    std::cout << "merge_pelee_filter #input_file #input_pelee_Np_file #input_pelee_0p_file #input_pelee_run_file #outfile" << std::endl;
    
    return -1;
  }

  TString input_file = argv[1];
  TString input_pelee_Np_file = argv[2];
  TString input_pelee_0p_file = argv[3];
  TString input_pelee_run_file = argv[4];
  TString outfile_name = argv[5];


  Int_t pl_run, pl_sub, pl_evt;
  Float_t pl_shr_energy_tot_cali, pl_shr_tkfit_dedx_Y;
  Float_t pl_reco_nu_vtx_sce_x, pl_reco_nu_vtx_sce_y, pl_reco_nu_vtx_sce_z;
  Double_t pl_reco_e;

  std::map<std::pair<int, int>, std::tuple<Float_t, Float_t, Float_t, Float_t, Float_t, Double_t> > map_pelee_np;
  std::map<std::pair<int, int>, std::tuple<Float_t, Float_t, Float_t, Float_t, Float_t, Double_t> > map_pelee_0p;
  std::set<std::pair<int, int> > set_pelee_rs;
  
  // read in PeLEE Np
  if (input_pelee_Np_file != "nan")
  {
    TFile *file_np = new TFile(input_pelee_Np_file);
    TTree *T_temp = (TTree*)file_np->Get("NeutrinoSelectionFilter");
  
    T_temp->SetBranchAddress("run",&pl_run);
    T_temp->SetBranchAddress("sub",&pl_sub);
    T_temp->SetBranchAddress("evt",&pl_evt);

    T_temp->SetBranchAddress("shr_energy_tot_cali",&pl_shr_energy_tot_cali);
    T_temp->SetBranchAddress("shr_tkfit_dedx_Y",&pl_shr_tkfit_dedx_Y);
    T_temp->SetBranchAddress("reco_nu_vtx_sce_x",&pl_reco_nu_vtx_sce_x);
    T_temp->SetBranchAddress("reco_nu_vtx_sce_y",&pl_reco_nu_vtx_sce_y);
    T_temp->SetBranchAddress("reco_nu_vtx_sce_z",&pl_reco_nu_vtx_sce_z);
    T_temp->SetBranchAddress("reco_e",&pl_reco_e);

    for (Int_t i=0;i!=T_temp->GetEntries();i++){
      T_temp->GetEntry(i);

      // std::cout << pl_shr_tkfit_dedx_Y << std::endl;
      
      map_pelee_np[std::make_pair(pl_run, pl_evt)] = std::make_tuple(pl_shr_energy_tot_cali, pl_shr_tkfit_dedx_Y,
								     pl_reco_nu_vtx_sce_x, pl_reco_nu_vtx_sce_y, pl_reco_nu_vtx_sce_z,
								     pl_reco_e);
    }
    
  }

  if (input_pelee_0p_file != "nan")
  // read in PeLEE 0p
  {
    TFile *file_0p = new TFile(input_pelee_0p_file);
    TTree *T_temp = (TTree*)file_0p->Get("NeutrinoSelectionFilter");

    T_temp->SetBranchAddress("run",&pl_run);
    T_temp->SetBranchAddress("sub",&pl_sub);
    T_temp->SetBranchAddress("evt",&pl_evt);

    T_temp->SetBranchAddress("shr_energy_tot_cali",&pl_shr_energy_tot_cali);
    T_temp->SetBranchAddress("shr_tkfit_dedx_Y",&pl_shr_tkfit_dedx_Y);
    T_temp->SetBranchAddress("reco_nu_vtx_sce_x",&pl_reco_nu_vtx_sce_x);
    T_temp->SetBranchAddress("reco_nu_vtx_sce_y",&pl_reco_nu_vtx_sce_y);
    T_temp->SetBranchAddress("reco_nu_vtx_sce_z",&pl_reco_nu_vtx_sce_z);
    T_temp->SetBranchAddress("reco_e",&pl_reco_e);

    for (Int_t i=0;i!=T_temp->GetEntries();i++){
      T_temp->GetEntry(i);
      map_pelee_0p[std::make_pair(pl_run, pl_evt)] = std::make_tuple(pl_shr_energy_tot_cali, pl_shr_tkfit_dedx_Y,
								     pl_reco_nu_vtx_sce_x, pl_reco_nu_vtx_sce_y, pl_reco_nu_vtx_sce_z,
								     pl_reco_e);
    }
  }

  // read in PeLEE run list
  if (input_pelee_run_file != "nan"){
    // to be added ...
    TFile *file_run = new TFile(input_pelee_run_file);
    TTree *T_temp = (TTree*)file_run->Get("SubRun");
    T_temp->SetBranchAddress("run",&pl_run);
    T_temp->SetBranchAddress("subRun",&pl_sub);
    for (Int_t i=0;i!=T_temp->GetEntries();i++){
      T_temp->GetEntry(i);
      set_pelee_rs.insert(std::make_pair(pl_run, pl_sub));
    }
  }


  std::cout << map_pelee_np.size() << " " << map_pelee_0p.size() << " " << set_pelee_rs.size() << std::endl;
  // return 0;
  
  TFile *file1 = new TFile(input_file);
  TTree *T_BDTvars = (TTree*)file1->Get("wcpselection/T_BDTvars");
  
  TTree *T_eval = (TTree*)file1->Get("wcpselection/T_eval");
  TTree *T_pot = (TTree*)file1->Get("wcpselection/T_pot");
  TTree *T_PFeval = (TTree*)file1->Get("wcpselection/T_PFeval");
  TTree *T_KINEvars = (TTree*)file1->Get("wcpselection/T_KINEvars");

  
  TFile *file2 = new TFile(outfile_name,"RECREATE");
  file2->mkdir("wcpselection");
  file2->cd("wcpselection");
  TTree *t4 = new TTree("T_BDTvars","T_BDTvars");
  TTree *t1 = new TTree("T_eval","T_eval");
  TTree *t2 = new TTree("T_pot","T_pot");
  TTree *t3 = new TTree("T_PFeval", "T_PFeval");
  TTree *t5 = new TTree("T_KINEvars", "T_KINEvars");

  bool flag_data = true;
  if (T_eval->GetBranch("weight_cv")) flag_data = false;
  
  EvalInfo eval;
  eval.file_type = new std::string();
  POTInfo pot;
  TaggerInfo tagger;
  PFevalInfo pfeval;
  KineInfo kine;

  kine.kine_energy_particle = new std::vector<float>;
  kine.kine_energy_info = new std::vector<int>;
  kine.kine_particle_type = new std::vector<int>;
  kine.kine_energy_included = new std::vector<int>;
    
  tagger.pio_2_v_dis2 = new std::vector<float>;
  tagger.pio_2_v_angle2 = new std::vector<float>;
  tagger.pio_2_v_acc_length = new std::vector<float>;
  tagger.pio_2_v_flag = new std::vector<float>;
  tagger.sig_1_v_angle = new std::vector<float>;
  tagger.sig_1_v_flag_single_shower = new std::vector<float>;
  tagger.sig_1_v_energy = new std::vector<float>;
  tagger.sig_1_v_energy_1 = new std::vector<float>;
  tagger.sig_1_v_flag = new std::vector<float>;
  tagger.sig_2_v_energy = new std::vector<float>;
  tagger.sig_2_v_shower_angle = new std::vector<float>;
  tagger.sig_2_v_flag_single_shower = new std::vector<float>;
  tagger.sig_2_v_medium_dQ_dx = new std::vector<float>;
  tagger.sig_2_v_start_dQ_dx = new std::vector<float>;
  tagger.sig_2_v_flag = new std::vector<float>;
  tagger.stw_2_v_medium_dQ_dx = new std::vector<float>;
  tagger.stw_2_v_energy = new std::vector<float>;
  tagger.stw_2_v_angle = new std::vector<float>;
  tagger.stw_2_v_dir_length = new std::vector<float>;
  tagger.stw_2_v_max_dQ_dx = new std::vector<float>;
  tagger.stw_2_v_flag = new std::vector<float>;
  tagger.stw_3_v_angle = new std::vector<float>;
  tagger.stw_3_v_dir_length = new std::vector<float>;
  tagger.stw_3_v_energy = new std::vector<float>;
  tagger.stw_3_v_medium_dQ_dx = new std::vector<float>;
  tagger.stw_3_v_flag = new std::vector<float>;
  tagger.stw_4_v_angle = new std::vector<float>;
  tagger.stw_4_v_dis = new std::vector<float>;
  tagger.stw_4_v_energy = new std::vector<float>;
  tagger.stw_4_v_flag = new std::vector<float>;
  tagger.br3_3_v_energy = new std::vector<float>;
  tagger.br3_3_v_angle = new std::vector<float>;
  tagger.br3_3_v_dir_length = new std::vector<float>;
  tagger.br3_3_v_length = new std::vector<float>;
  tagger.br3_3_v_flag = new std::vector<float>;
  tagger.br3_5_v_dir_length = new std::vector<float>;
  tagger.br3_5_v_total_length = new std::vector<float>;
  tagger.br3_5_v_flag_avoid_muon_check = new std::vector<float>;
  tagger.br3_5_v_n_seg = new std::vector<float>;
  tagger.br3_5_v_angle = new std::vector<float>;
  tagger.br3_5_v_sg_length = new std::vector<float>;
  tagger.br3_5_v_energy = new std::vector<float>;
  tagger.br3_5_v_n_main_segs = new std::vector<float>;
  tagger.br3_5_v_n_segs = new std::vector<float>;
  tagger.br3_5_v_shower_main_length = new std::vector<float>;
  tagger.br3_5_v_shower_total_length = new std::vector<float>;
  tagger.br3_5_v_flag = new std::vector<float>;
  tagger.br3_6_v_angle = new std::vector<float>;
  tagger.br3_6_v_angle1 = new std::vector<float>;
  tagger.br3_6_v_flag_shower_trajectory = new std::vector<float>;
  tagger.br3_6_v_direct_length = new std::vector<float>;
  tagger.br3_6_v_length = new std::vector<float>;
  tagger.br3_6_v_n_other_vtx_segs = new std::vector<float>;
  tagger.br3_6_v_energy = new std::vector<float>;
  tagger.br3_6_v_flag = new std::vector<float>;
  tagger.tro_1_v_particle_type = new std::vector<float>;
  tagger.tro_1_v_flag_dir_weak = new std::vector<float>;
  tagger.tro_1_v_min_dis = new std::vector<float>;
  tagger.tro_1_v_sg1_length = new std::vector<float>;
  tagger.tro_1_v_shower_main_length = new std::vector<float>;
  tagger.tro_1_v_max_n_vtx_segs = new std::vector<float>;
  tagger.tro_1_v_tmp_length = new std::vector<float>;
  tagger.tro_1_v_medium_dQ_dx = new std::vector<float>;
  tagger.tro_1_v_dQ_dx_cut = new std::vector<float>;
  tagger.tro_1_v_flag_shower_topology = new std::vector<float>;
  tagger.tro_1_v_flag = new std::vector<float>;
  tagger.tro_2_v_energy = new std::vector<float>;
  tagger.tro_2_v_stem_length = new std::vector<float>;
  tagger.tro_2_v_iso_angle = new std::vector<float>;
  tagger.tro_2_v_max_length = new std::vector<float>;
  tagger.tro_2_v_angle = new std::vector<float>;
  tagger.tro_2_v_flag = new std::vector<float>;
  tagger.tro_4_v_dir2_mag = new std::vector<float>;
  tagger.tro_4_v_angle = new std::vector<float>;
  tagger.tro_4_v_angle1 = new std::vector<float>;
  tagger.tro_4_v_angle2 = new std::vector<float>;
  tagger.tro_4_v_length = new std::vector<float>;
  tagger.tro_4_v_length1 = new std::vector<float>;
  tagger.tro_4_v_medium_dQ_dx = new std::vector<float>;
  tagger.tro_4_v_end_dQ_dx = new std::vector<float>;
  tagger.tro_4_v_energy = new std::vector<float>;
  tagger.tro_4_v_shower_main_length = new std::vector<float>;
  tagger.tro_4_v_flag_shower_trajectory = new std::vector<float>;
  tagger.tro_4_v_flag = new std::vector<float>;
  tagger.tro_5_v_max_angle = new std::vector<float>;
  tagger.tro_5_v_min_angle = new std::vector<float>;
  tagger.tro_5_v_max_length = new std::vector<float>;
  tagger.tro_5_v_iso_angle = new std::vector<float>;
  tagger.tro_5_v_n_vtx_segs = new std::vector<float>;
  tagger.tro_5_v_min_count = new std::vector<float>;
  tagger.tro_5_v_max_count = new std::vector<float>;
  tagger.tro_5_v_energy = new std::vector<float>;
  tagger.tro_5_v_flag = new std::vector<float>;
  tagger.lol_1_v_energy = new std::vector<float>;
  tagger.lol_1_v_vtx_n_segs = new std::vector<float>;
  tagger.lol_1_v_nseg = new std::vector<float>;
  tagger.lol_1_v_angle = new std::vector<float>;
  tagger.lol_1_v_flag = new std::vector<float>;
  tagger.lol_2_v_length = new std::vector<float>;
  tagger.lol_2_v_angle = new std::vector<float>;
  tagger.lol_2_v_type = new std::vector<float>;
  tagger.lol_2_v_vtx_n_segs = new std::vector<float>;
  tagger.lol_2_v_energy = new std::vector<float>;
  tagger.lol_2_v_shower_main_length = new std::vector<float>;
  tagger.lol_2_v_flag_dir_weak = new std::vector<float>;
  tagger.lol_2_v_flag = new std::vector<float>;
  tagger.cosmict_flag_10 = new std::vector<float>;
  tagger.cosmict_10_flag_inside = new std::vector<float>;
  tagger.cosmict_10_vtx_z = new std::vector<float>;
  tagger.cosmict_10_flag_shower = new std::vector<float>;
  tagger.cosmict_10_flag_dir_weak = new std::vector<float>;
  tagger.cosmict_10_angle_beam = new std::vector<float>;
  tagger.cosmict_10_length = new std::vector<float>;
  tagger.numu_cc_flag_1 = new std::vector<float>;
  tagger.numu_cc_1_particle_type = new std::vector<float>;
  tagger.numu_cc_1_length = new std::vector<float>;
  tagger.numu_cc_1_medium_dQ_dx = new std::vector<float>;
  tagger.numu_cc_1_dQ_dx_cut = new std::vector<float>;
  tagger.numu_cc_1_direct_length = new std::vector<float>;
  tagger.numu_cc_1_n_daughter_tracks = new std::vector<float>;
  tagger.numu_cc_1_n_daughter_all = new std::vector<float>;
  tagger.numu_cc_flag_2 = new std::vector<float>;
  tagger.numu_cc_2_length = new std::vector<float>;
  tagger.numu_cc_2_total_length = new std::vector<float>;
  tagger.numu_cc_2_n_daughter_tracks = new std::vector<float>;
  tagger.numu_cc_2_n_daughter_all = new std::vector<float>;
  tagger.pio_2_v_dis2 = new std::vector<float>;
  tagger.pio_2_v_angle2 = new std::vector<float>;
  tagger.pio_2_v_acc_length = new std::vector<float>;
  tagger.pio_2_v_flag = new std::vector<float>;
  tagger.sig_1_v_angle = new std::vector<float>;
  tagger.sig_1_v_flag_single_shower = new std::vector<float>;
  tagger.sig_1_v_energy = new std::vector<float>;
  tagger.sig_1_v_energy_1 = new std::vector<float>;
  tagger.sig_1_v_flag = new std::vector<float>;
  tagger.sig_2_v_energy = new std::vector<float>;
  tagger.sig_2_v_shower_angle = new std::vector<float>;
  tagger.sig_2_v_flag_single_shower = new std::vector<float>;
  tagger.sig_2_v_medium_dQ_dx = new std::vector<float>;
  tagger.sig_2_v_start_dQ_dx = new std::vector<float>;
  tagger.sig_2_v_flag = new std::vector<float>;
  tagger.stw_2_v_medium_dQ_dx = new std::vector<float>;
  tagger.stw_2_v_energy = new std::vector<float>;
  tagger.stw_2_v_angle = new std::vector<float>;
  tagger.stw_2_v_dir_length = new std::vector<float>;
  tagger.stw_2_v_max_dQ_dx = new std::vector<float>;
  tagger.stw_2_v_flag = new std::vector<float>;
  tagger.stw_3_v_angle = new std::vector<float>;
  tagger.stw_3_v_dir_length = new std::vector<float>;
  tagger.stw_3_v_energy = new std::vector<float>;
  tagger.stw_3_v_medium_dQ_dx = new std::vector<float>;
  tagger.stw_3_v_flag = new std::vector<float>;
  tagger.stw_4_v_angle = new std::vector<float>;
  tagger.stw_4_v_dis = new std::vector<float>;
  tagger.stw_4_v_energy = new std::vector<float>;
  tagger.stw_4_v_flag = new std::vector<float>;
  tagger.br3_3_v_energy = new std::vector<float>;
  tagger.br3_3_v_angle = new std::vector<float>;
  tagger.br3_3_v_dir_length = new std::vector<float>;
  tagger.br3_3_v_length = new std::vector<float>;
  tagger.br3_3_v_flag = new std::vector<float>;
  tagger.br3_5_v_dir_length = new std::vector<float>;
  tagger.br3_5_v_total_length = new std::vector<float>;
  tagger.br3_5_v_flag_avoid_muon_check = new std::vector<float>;
  tagger.br3_5_v_n_seg = new std::vector<float>;
  tagger.br3_5_v_angle = new std::vector<float>;
  tagger.br3_5_v_sg_length = new std::vector<float>;
  tagger.br3_5_v_energy = new std::vector<float>;
  tagger.br3_5_v_n_main_segs = new std::vector<float>;
  tagger.br3_5_v_n_segs = new std::vector<float>;
  tagger.br3_5_v_shower_main_length = new std::vector<float>;
  tagger.br3_5_v_shower_total_length = new std::vector<float>;
  tagger.br3_5_v_flag = new std::vector<float>;
  tagger.br3_6_v_angle = new std::vector<float>;
  tagger.br3_6_v_angle1 = new std::vector<float>;
  tagger.br3_6_v_flag_shower_trajectory = new std::vector<float>;
  tagger.br3_6_v_direct_length = new std::vector<float>;
  tagger.br3_6_v_length = new std::vector<float>;
  tagger.br3_6_v_n_other_vtx_segs = new std::vector<float>;
  tagger.br3_6_v_energy = new std::vector<float>;
  tagger.br3_6_v_flag = new std::vector<float>;
  tagger.tro_1_v_particle_type = new std::vector<float>;
  tagger.tro_1_v_flag_dir_weak = new std::vector<float>;
  tagger.tro_1_v_min_dis = new std::vector<float>;
  tagger.tro_1_v_sg1_length = new std::vector<float>;
  tagger.tro_1_v_shower_main_length = new std::vector<float>;
  tagger.tro_1_v_max_n_vtx_segs = new std::vector<float>;
  tagger.tro_1_v_tmp_length = new std::vector<float>;
  tagger.tro_1_v_medium_dQ_dx = new std::vector<float>;
  tagger.tro_1_v_dQ_dx_cut = new std::vector<float>;
  tagger.tro_1_v_flag_shower_topology = new std::vector<float>;
  tagger.tro_1_v_flag = new std::vector<float>;
  tagger.tro_2_v_energy = new std::vector<float>;
  tagger.tro_2_v_stem_length = new std::vector<float>;
  tagger.tro_2_v_iso_angle = new std::vector<float>;
  tagger.tro_2_v_max_length = new std::vector<float>;
  tagger.tro_2_v_angle = new std::vector<float>;
  tagger.tro_2_v_flag = new std::vector<float>;
  tagger.tro_4_v_dir2_mag = new std::vector<float>;
  tagger.tro_4_v_angle = new std::vector<float>;
  tagger.tro_4_v_angle1 = new std::vector<float>;
  tagger.tro_4_v_angle2 = new std::vector<float>;
  tagger.tro_4_v_length = new std::vector<float>;
  tagger.tro_4_v_length1 = new std::vector<float>;
  tagger.tro_4_v_medium_dQ_dx = new std::vector<float>;
  tagger.tro_4_v_end_dQ_dx = new std::vector<float>;
  tagger.tro_4_v_energy = new std::vector<float>;
  tagger.tro_4_v_shower_main_length = new std::vector<float>;
  tagger.tro_4_v_flag_shower_trajectory = new std::vector<float>;
  tagger.tro_4_v_flag = new std::vector<float>;
  tagger.tro_5_v_max_angle = new std::vector<float>;
  tagger.tro_5_v_min_angle = new std::vector<float>;
  tagger.tro_5_v_max_length = new std::vector<float>;
  tagger.tro_5_v_iso_angle = new std::vector<float>;
  tagger.tro_5_v_n_vtx_segs = new std::vector<float>;
  tagger.tro_5_v_min_count = new std::vector<float>;
  tagger.tro_5_v_max_count = new std::vector<float>;
  tagger.tro_5_v_energy = new std::vector<float>;
  tagger.tro_5_v_flag = new std::vector<float>;
  tagger.lol_1_v_energy = new std::vector<float>;
  tagger.lol_1_v_vtx_n_segs = new std::vector<float>;
  tagger.lol_1_v_nseg = new std::vector<float>;
  tagger.lol_1_v_angle = new std::vector<float>;
  tagger.lol_1_v_flag = new std::vector<float>;
  tagger.lol_2_v_length = new std::vector<float>;
  tagger.lol_2_v_angle = new std::vector<float>;
  tagger.lol_2_v_type = new std::vector<float>;
  tagger.lol_2_v_vtx_n_segs = new std::vector<float>;
  tagger.lol_2_v_energy = new std::vector<float>;
  tagger.lol_2_v_shower_main_length = new std::vector<float>;
  tagger.lol_2_v_flag_dir_weak = new std::vector<float>;
  tagger.lol_2_v_flag = new std::vector<float>;
  tagger.cosmict_flag_10 = new std::vector<float>;
  tagger.cosmict_10_flag_inside = new std::vector<float>;
  tagger.cosmict_10_vtx_z = new std::vector<float>;
  tagger.cosmict_10_flag_shower = new std::vector<float>;
  tagger.cosmict_10_flag_dir_weak = new std::vector<float>;
  tagger.cosmict_10_angle_beam = new std::vector<float>;
  tagger.cosmict_10_length = new std::vector<float>;
  tagger.numu_cc_flag_1 = new std::vector<float>;
  tagger.numu_cc_1_particle_type = new std::vector<float>;
  tagger.numu_cc_1_length = new std::vector<float>;
  tagger.numu_cc_1_medium_dQ_dx = new std::vector<float>;
  tagger.numu_cc_1_dQ_dx_cut = new std::vector<float>;
  tagger.numu_cc_1_direct_length = new std::vector<float>;
  tagger.numu_cc_1_n_daughter_tracks = new std::vector<float>;
  tagger.numu_cc_1_n_daughter_all = new std::vector<float>;
  tagger.numu_cc_flag_2 = new std::vector<float>;
  tagger.numu_cc_2_length = new std::vector<float>;
  tagger.numu_cc_2_total_length = new std::vector<float>;
  tagger.numu_cc_2_n_daughter_tracks = new std::vector<float>;
  tagger.numu_cc_2_n_daughter_all = new std::vector<float>;

  set_tree_address(T_BDTvars, tagger,2 );
  put_tree_address(t4, tagger,2);

  if (flag_data){
    set_tree_address(T_eval, eval,2);
    put_tree_address(t1, eval,2);

    set_tree_address(T_PFeval, pfeval,2);
    put_tree_address(t3, pfeval,2);
  }else{
    set_tree_address(T_eval, eval);
    put_tree_address(t1, eval);

    set_tree_address(T_PFeval, pfeval);
    put_tree_address(t3, pfeval);
  }

  set_tree_address(T_pot, pot);
  put_tree_address(t2, pot);

  set_tree_address(T_KINEvars, kine);
  put_tree_address(t5, kine);
    
  T_eval->SetBranchStatus("*",1);
  T_BDTvars->SetBranchStatus("*",1);

  bool flag_remove = true;
  Int_t pl_flag;
  // add the new information to eval tree
  t1->Branch("pl_shr_energy_tot_cali",&pl_shr_energy_tot_cali,"pl_shr_energy_tot_cali/F");
  t1->Branch("pl_shr_tkfit_dedx_Y",&pl_shr_tkfit_dedx_Y,"pl_shr_tkfit_dedx_Y/F");
  t1->Branch("pl_reco_nu_vtx_sce_x",&pl_reco_nu_vtx_sce_x,"pl_reco_nu_vtx_sce_x/F");
  t1->Branch("pl_reco_nu_vtx_sce_y",&pl_reco_nu_vtx_sce_y,"pl_reco_nu_vtx_sce_y/F");
  t1->Branch("pl_reco_nu_vtx_sce_z",&pl_reco_nu_vtx_sce_z,"pl_reco_nu_vtx_sce_z/F");
  t1->Branch("pl_reco_e",&pl_reco_e,"pl_reco_e/D");
  t1->Branch("pl_flag",&pl_flag,"pl_flag/I");
  
  for (int i=0;i!=T_BDTvars->GetEntries();i++){
    T_BDTvars->GetEntry(i);
    T_eval->GetEntry(i); tagger.match_isFC = eval.match_isFC;
    T_KINEvars->GetEntry(i); tagger.kine_reco_Enu = kine.kine_reco_Enu;
    T_PFeval->GetEntry(i);

    if (set_pelee_rs.size()!=0)
    // not in common list ...
      if (set_pelee_rs.find(std::make_pair(eval.run, eval.subrun) ) == set_pelee_rs.end() ) continue;

    // reset the files ...
    pl_flag = 0;
    pl_shr_energy_tot_cali = 0;
    pl_shr_tkfit_dedx_Y = 0;
    pl_reco_nu_vtx_sce_x = 0;
    pl_reco_nu_vtx_sce_y = 0;
    pl_reco_nu_vtx_sce_z = 0;
    pl_reco_e = 0;

    auto it1 = map_pelee_np.find(std::make_pair(eval.run, eval.event) );
    auto it2 = map_pelee_0p.find(std::make_pair(eval.run, eval.event) );
    if (it1 != map_pelee_np.end() ){
      pl_flag = 1; // Np
      pl_shr_energy_tot_cali = std::get<0>(it1->second);
      pl_shr_tkfit_dedx_Y = std::get<1>(it1->second);
      pl_reco_nu_vtx_sce_x = std::get<2>(it1->second);
      pl_reco_nu_vtx_sce_y = std::get<3>(it1->second);
      pl_reco_nu_vtx_sce_z = std::get<4>(it1->second);
      pl_reco_e = std::get<5>(it1->second);
      
    }else if (it2 != map_pelee_0p.end()){
      pl_flag = 2; // 0p

      pl_shr_energy_tot_cali = std::get<0>(it2->second);
      pl_shr_tkfit_dedx_Y = std::get<1>(it2->second);
      pl_reco_nu_vtx_sce_x = std::get<2>(it2->second);
      pl_reco_nu_vtx_sce_y = std::get<3>(it2->second);
      pl_reco_nu_vtx_sce_z = std::get<4>(it2->second);
      pl_reco_e = std::get<5>(it2->second);
    }
    
                
    t4->Fill();
    t1->Fill();
    t3->Fill();
    t5->Fill();
  }

  for (Int_t i=0;i!=T_pot->GetEntries();i++){
    T_pot->GetEntry(i);

    if (set_pelee_rs.size() != 0)
      // not in common list ...
      if (set_pelee_rs.find(std::make_pair(pot.runNo, pot.subRunNo) ) == set_pelee_rs.end() ) continue;
    
    t2->Fill();
  }


  
  file2->Write();
  file2->Close();


  return 0;

  
  
}
