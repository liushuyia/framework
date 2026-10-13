#include <cstdlib>
#include <iostream>
#include <fstream>
#include <map>
#include <string>

#include "TChain.h"
#include "TFile.h"
#include "TTree.h"
#include "TString.h"
#include "TObjString.h"
#include "TSystem.h"
#include "TROOT.h"
#include "TMath.h"

#include "WCPLEEANA/pot.h"

using namespace std;
using namespace LEEana;



int main(int argc, char** argv){

  if (argc < 2){
    std::cout << "pot_counting #extbnb_file " << std::endl;
    return -1;
  }

  int mode = 2;
  
  TString extbnb_file = argv[1];
 

  int run, subrun;
  double trigger_no;
  double total_bnb_trigger_no = 0;
  double total_extbnb_trigger_no = 0;

  std::map<int, std::tuple<float, float> > map_run_pte;
  
  
  std::map<std::pair<int, int>, double >  map_extbnb_infos;
  ifstream infile1("pot_extbnb.txt");
  while(!infile1.eof()){
    infile1 >> run >> subrun >> trigger_no;
    if (run >0)
      map_extbnb_infos[std::make_pair(run, subrun)] = trigger_no;
  }

  TFile *file2 = new TFile(extbnb_file);
  TTree *T_pot = (TTree*)file2->Get("wcpselection/T_pot");
  TTree *T_eval = (TTree*)file2->Get("wcpselection/T_eval");
  TTree *T_BDTvars = (TTree*)file2->Get("wcpselection/T_BDTvars");
  TTree *T_PFeval = (TTree*)file2->Get("wcpselection/T_PFeval");
  TTree *T_KINEvars = (TTree*)file2->Get("wcpselection/T_KINEvars");
  
  T_eval->SetBranchStatus("*",0);
  T_BDTvars->SetBranchStatus("*",0);
  T_PFeval->SetBranchStatus("*",0);
  T_KINEvars->SetBranchStatus("*",0);
  
  Bool_t match_isFC;
  Float_t match_energy;
  T_eval->SetBranchStatus("match_isFC",1); T_eval->SetBranchAddress("match_isFC",&match_isFC);
  T_eval->SetBranchStatus("run",1); T_eval->SetBranchAddress("run",&run);
  T_eval->SetBranchStatus("match_energy",1); T_eval->SetBranchAddress("match_energy", &match_energy);
  Float_t numu_score;
  Float_t nue_score;
  Float_t numu_cc_flag;
  Float_t nc_pio_score;
  Float_t nc_delta_score;
  T_BDTvars->SetBranchStatus("numu_score",1); T_BDTvars->SetBranchAddress("numu_score",&numu_score);
  T_BDTvars->SetBranchStatus("nue_score",1); T_BDTvars->SetBranchAddress("nue_score",&nue_score);
  T_BDTvars->SetBranchStatus("numu_cc_flag",1); T_BDTvars->SetBranchAddress("numu_cc_flag",&numu_cc_flag);
  if (T_BDTvars->GetBranch("nc_delta_score")){
    T_BDTvars->SetBranchStatus("nc_delta_score",1); T_BDTvars->SetBranchAddress("nc_delta_score",&nc_delta_score);
    T_BDTvars->SetBranchStatus("nc_pio_score",1); T_BDTvars->SetBranchAddress("nc_pio_score",&nc_pio_score);
  }
  Float_t kine_reco_Enu;
  T_KINEvars->SetBranchStatus("kine_reco_Enu",1); T_KINEvars->SetBranchAddress("kine_reco_Enu",&kine_reco_Enu);
  
  
  POTInfo pot;
  set_tree_address(T_pot, pot);
  float pass_ratio=1;
  if (T_pot->GetBranch("pass_ratio")) T_pot->SetBranchAddress("pass_ratio", &pass_ratio);
  
   
  for (Int_t i=0;i!=T_pot->GetEntries();i++){
    T_pot->GetEntry(i);
    auto it = map_extbnb_infos.find(std::make_pair(pot.runNo, pot.subRunNo));

    if (map_run_pte.find(pot.runNo) == map_run_pte.end())
      map_run_pte[pot.runNo] = std::make_tuple(0,0);

    auto it1 = map_run_pte.find(pot.runNo);
    
    if (it == map_extbnb_infos.end()  && pot.runNo >=4000  && pot.runNo <= 50000){
      std::cout << "Run: " << pot.runNo << " subRun: " << pot.subRunNo << "  not found!" << std::endl;
    }else{
      total_extbnb_trigger_no  += it->second * pass_ratio;
      std::get<0>(it1->second) += it->second * pass_ratio;
    }
  }

  int num_events = 0;
  for (Int_t i=0;i!= T_eval->GetEntries();i++){
    T_eval->GetEntry(i);
    T_BDTvars->GetEntry(i);
    T_KINEvars->GetEntry(i);

    

    
    // generic neutrino ...
    if (numu_cc_flag>=0
	//	  && match_energy < 300
	){
      auto it1 = map_run_pte.find(run);
      std::get<1>(it1->second) ++;
      //      num_events ++;
    } 
  }
  
  for (auto it = map_run_pte.begin(); it != map_run_pte.end(); it++){
    int run = it->first;
    float ntrigger = std::get<0>(it->second);
    float nobs = std::get<1>(it->second);
    if (run >=4000  && run <= 50000)
      std::cout << run << " " << ntrigger << " " << nobs  << std::endl;
    }
    //    std::cout << "EXTBNB: total trigger: " << total_extbnb_trigger_no << std::endl;
    //std::cout << "EXTBNB: total_pot: " << total_extbnb_trigger_no/total_bnb_trigger_no * total_bnb_pot << std::endl;
}
  
