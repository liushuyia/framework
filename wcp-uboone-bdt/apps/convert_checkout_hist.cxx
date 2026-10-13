#include <iostream>

#include <fstream>
#include <sstream>
#include <vector>
#include <cmath>
#include <iomanip>
#include <algorithm>
#include "TCanvas.h"
#include "TLegend.h"
#include "TLine.h"

#include "WCPLEEANA/master_cov_matrix.h"

#include "TROOT.h"
#include "TMath.h"
#include "TH1F.h"
#include "TFile.h"
#include "TTree.h"

#include "WCPLEEANA/cuts.h"
#include "WCPLEEANA/pot.h"
#include "WCPLEEANA/pfeval.h"
#include "WCPLEEANA/eval.h"

using namespace std;
using namespace LEEana;


// ============================================================
// NUE flux reweight switch
// true  = apply nue flux reweight
// false = do not apply nue flux reweight
// ============================================================
const bool USE_NUE_FLUX_REWEIGHT = true;

// ============================================================
// ROOTino bug correction switch
//
// true  = apply ROOTino correction
// false = do not apply ROOTino correction
//
// ROOTino correction is for:
//     SURPRISE + nu_overlay
//
// ============================================================
const bool USE_ROOTINO_REWEIGHT = true;

const TString ROOTINO_RUN4B_FILE = "/pnfs/uboone/persistent/users/uboonepro/surprise/retuple/BNB/checkout_MCC9.10_Run4b_v10_04_07_20_BNB_nu_overlay_retuple_retuple_hist.root";
const TString ROOTINO_RUN4C_FILE = "/pnfs/uboone/persistent/users/uboonepro/surprise/retuple/BNB/checkout_MCC9.10_Run4acd5_v10_04_07_20_BNB_nu_overlay_retuple_retuple_hist_4c.root";
const TString ROOTINO_RUN4D_FILE = "/pnfs/uboone/persistent/users/uboonepro/surprise/retuple/BNB/checkout_MCC9.10_Run4acd5_v10_04_07_20_BNB_nu_overlay_retuple_retuple_hist_4d.root";
const TString ROOTINO_RUN5_FILE = "/pnfs/uboone/persistent/users/uboonepro/surprise/retuple/BNB/checkout_MCC9.10_Run4acd5_v10_04_07_20_BNB_nu_overlay_retuple_retuple_hist_5.root";


bool load_nue_flux_target(const TString& csv_file, std::vector<double>& target, std::vector<double>& target_unweight, double& bin_width)
{
    std::ifstream fin(csv_file.Data());

    if (!fin.is_open()) 
    {
      std::cerr << "[NUE FLUX] ERROR: cannot open target file: " << csv_file << std::endl;
      return false;
    }

    std::string line;

    // --------------------------------------------------------
    // Read header
    // --------------------------------------------------------
    if (!std::getline(fin, line)) 
    {
      std::cerr << "[NUE FLUX] ERROR: empty target CSV." << std::endl;
      return false;
    }

    std::stringstream header_ss(line);

    std::vector<std::string> header;
    std::string item;

    while (std::getline(header_ss, item, ',')) 
    {
      header.push_back(item);
    }

    int idx_lo = -1;
    int idx_hi = -1;
    int idx_target_unweight = -1;
    int idx_target = -1;

    for (size_t i = 0; i < header.size(); ++i) 
    {
      if (header[i] == "bin_lo_gev")
          idx_lo = i;

      if (header[i] == "bin_hi_gev")
          idx_hi = i;

      if (header[i] == "target_no_xs_cv_weight_per_1e17pot")
          idx_target_unweight = i;

      if (header[i] == "target_xs_cv_weight_per_1e17pot")
          idx_target = i;
    }

    if (idx_lo < 0 || idx_hi < 0 || idx_target < 0 || idx_target_unweight < 0) 
    {
        std::cerr << "[NUE FLUX] ERROR: required columns " "not found in CSV." << std::endl;
        std::cerr << "Required:" << std::endl << "  bin_lo_gev" << std::endl << "  bin_hi_gev" << std::endl << "  target_xs_cv_weight_per_1e17pot" << "  target_no_xs_cv_weight_per_1e17pot" << std::endl;
        return false;
    }

    target.clear();
    target_unweight.clear();

    std::vector<double> bin_lo;
    std::vector<double> bin_hi;

    // --------------------------------------------------------
    // Read table
    // --------------------------------------------------------
    while (std::getline(fin, line))
    {
      if (line.empty())
          continue;

      std::stringstream ss(line);
      std::vector<std::string> values;

      while (std::getline(ss, item, ',')) 
      {
        values.push_back(item);
      }

      int max_idx = std::max(std::max(idx_target, idx_target_unweight), std::max(idx_lo, idx_hi));

      if ((int)values.size() <= max_idx)
          continue;

      double lo = std::stod(values[idx_lo]);
      double hi = std::stod(values[idx_hi]);
      double target_value_unweight = std::stod(values[idx_target_unweight]);
      double target_value = std::stod(values[idx_target]);

      bin_lo.push_back(lo);
      bin_hi.push_back(hi);
      target_unweight.push_back(target_value_unweight);
      target.push_back(target_value);
    }

    fin.close();

    if (target.empty()) 
    {
      std::cerr << "[NUE FLUX] ERROR: target table is empty." << std::endl;
      return false;
    }

    // --------------------------------------------------------
    // Check binning
    //
    //  - starts at 0
    //  - contiguous
    //  - uniform
    //
    // --------------------------------------------------------

    bin_width = bin_hi[0] - bin_lo[0];
    if (std::fabs(bin_lo[0]) > 1e-10) 
    {
      std::cerr << "[NUE FLUX] ERROR: first bin does " "not start at 0 GeV." << std::endl;
      return false;
    }

    for (size_t i = 0; i < bin_lo.size(); ++i) 
    {
      double width = bin_hi[i] - bin_lo[i];
      if (std::fabs(width - bin_width) > 1e-8) 
      {
        std::cerr << "[NUE FLUX] ERROR: bins are " "not uniform." << std::endl;
        return false;
      }
      if (i > 0) 
      {
        if (std::fabs(bin_lo[i] - bin_hi[i-1]) > 1e-8) 
        {
          std::cerr << "[NUE FLUX] ERROR: bins are " "not contiguous." << std::endl;
          return false;
        }
      }

      if (!std::isfinite(target[i]) ||  target[i] < 0.0) 
      {
        std::cerr << "[NUE FLUX] ERROR: invalid " "target value in bin " << i << std::endl;
        return false;
      }
    }

    return true;
}


// ============================================================
// Build the nue flux sampling weight for ONE ROOT file.
//
// This must be done using the ENTIRE file
// ============================================================

bool build_nue_flux_weights(TTree* T_eval, TTree* T_pot, EvalInfo& eval, POTInfo& pot, const std::vector<double>& target, const std::vector<double>& target_unweight, double bin_width, std::vector<double>& nue_flux_weight, double& file_POT, bool is_run123)
{
  const int n_bins = target.size();
  nue_flux_weight.assign(n_bins, 0.0);
  std::vector<double> weighted_counts(n_bins, 0.0);
  std::vector<Long64_t> raw_counts(n_bins, 0);

    // --------------------------------------------------------
    // Calculate POT.
    //
    // uses pot_tor875.
    // --------------------------------------------------------

    file_POT = 0.0;

    for (Long64_t i = 0; i < T_pot->GetEntries(); ++i) 
    {
      T_pot->GetEntry(i);
      file_POT += pot.pot_tor875;
    }

    if (file_POT <= 0.0) 
    {
      std::cerr << "[NUE FLUX] ERROR: " "file POT <= 0." << std::endl;
      return false;
    }

    // --------------------------------------------------------
    // First pass over ENTIRE T_eval.
    //
    Long64_t total_events = T_eval->GetEntries();

    Long64_t in_range_events = 0;
    Long64_t out_of_range_events = 0;
    Long64_t nonfinite_events = 0;

    Long64_t zero_weight_events = 0;

    for (Long64_t i = 0; i < T_eval->GetEntries(); ++i) 
    {
      T_eval->GetEntry(i);
      double enu = eval.truth_nuEnergy;
      if (!std::isfinite(enu)) 
      {
        nonfinite_events++;
        std::cerr << "[NUE FLUX] ERROR: "  "non-finite truth_nuEnergy " << "at entry " << i << std::endl;
        continue;
        // return false;
      }

      double enu_gev = enu / 1000.0;   /////////////////////// convert to GeV
      if (isnan(eval.weight_cv)||isinf(eval.weight_cv)|| eval.weight_cv<=0|| eval.weight_cv>1000)  eval.weight_cv=1;
      if (isnan(eval.weight_spline)||isinf(eval.weight_spline)||eval.weight_spline<=0||eval.weight_spline>1000) eval.weight_spline=1;
      float cv_spline = eval.weight_cv * eval.weight_spline;


      int bin = static_cast<int>(std::floor(enu_gev / bin_width));

    /////////////////////////////////////////
      //
      // in_range = (bin_idx >= 0) & (bin_idx < n_bins)
      // outside CSV range

      if (bin < 0 || bin >= n_bins) 
      {
        out_of_range_events++;
        continue;
      }

      // --------------------------------------------------------
      // Inside CSV range
      // --------------------------------------------------------
      in_range_events++;

      // Raw event count
      raw_counts[bin]++;

      if (isnan(eval.weight_cv)||isinf(eval.weight_cv)|| eval.weight_cv<=0|| eval.weight_cv>1000)  eval.weight_cv=1;
      if (isnan(eval.weight_spline)||isinf(eval.weight_spline)||eval.weight_spline<=0||eval.weight_spline>1000) eval.weight_spline=1;
      cv_spline = eval.weight_cv * eval.weight_spline;

      weighted_counts[bin] += cv_spline;
    }

    // --------------------------------------------------------
    // Convert to events / 1e17 POT
    //
    // file_rate = weighted_counts / (file_POT / 1e17)
    // --------------------------------------------------------

    std::vector<double> file_rate(n_bins, 0.0); 
    std::vector<double> file_rate_unweight(n_bins, 0.0); 
    if (is_run123)
    {
      for (int i = 0; i < n_bins; ++i) 
      {
        file_rate[i] = weighted_counts[i] / (file_POT / 1.0e17);
        if (weighted_counts[i] > 0.0) 
        {
          nue_flux_weight[i] = target[i] / file_rate[i];
        } else 
        {
          nue_flux_weight[i] = 0.0;
        }
      }
    }
    else
    {
      for (int i = 0; i < n_bins; ++i) 
      {
        file_rate_unweight[i] = raw_counts[i] / (file_POT / 1.0e17);
        if (raw_counts[i] > 0.0) 
        {
          nue_flux_weight[i] = target_unweight[i] / file_rate_unweight[i];
        } else 
        {
          nue_flux_weight[i] = 0.0;
        }
      }
    }

    int zero_event_bins = 0;
    int zero_weight_bins = 0;

    for (int i = 0; i < n_bins; ++i)
    {
      if (raw_counts[i] == 0)
          zero_event_bins++;

      if (nue_flux_weight[i] == 0.0)
          zero_weight_bins++;
    }

    

    int zero_target_bins = 0;
    for (int i = 0; i < n_bins; ++i)
    {
      if (target[i] == 0.0)
          zero_target_bins++;
    }
    // --------------------------------------------------------
    // Closure test
    //
    // post = sum(nue_flux_weight * cv_spline) / (POT / 1e17)
    //
    // Should equal target.
    // --------------------------------------------------------

    double post_integral = 0.0;
    double target_integral = 0.0;
    double post_integral_unweight = 0.0;
    double target_integral_unweight = 0.0;

    for (int i = 0; i < n_bins; ++i) 
    {
      post_integral += weighted_counts[i] * nue_flux_weight[i] / (file_POT / 1.0e17);
      target_integral += target[i];
      post_integral_unweight += raw_counts[i] * nue_flux_weight[i] / (file_POT / 1.0e17);
      target_integral_unweight += target_unweight[i];
    }

    
    return true;
}



// ============================================================
// Return nue flux sampling weight for one event.
// ============================================================

double get_nue_flux_weight(double enu_meV, const std::vector<double>& nue_flux_weight, double bin_width_GeV)
{
  if (!std::isfinite(enu_meV))
      return 0.0;
  // !!! enu_meV → GeV
  double enu_gev = enu_meV / 1000.0;
  int bin = static_cast<int>(std::floor(enu_gev / bin_width_GeV));
  if (bin < 0 || bin >= (int)nue_flux_weight.size()) 
  {
    return 0.0;
  }
  return nue_flux_weight[bin];
}

// ============================================================
// Count Delta(1600) events in eventweight_tree
//
// GTruth_ResNum == 9 -> Delta(1600)
//
// This is the raw event count in eventweight_tree.
// ============================================================
Long64_t count_ROOTino_Delta1600(TTree* T_eventweight)
{
    if (!T_eventweight)
    {
        std::cerr
            << "[ROOTINO] ERROR: null eventweight_tree."
            << std::endl;
        return 0;
    }

    if (!T_eventweight->GetBranch("GTruth_ResNum"))
    {
        std::cerr
            << "[ROOTINO] ERROR: GTruth_ResNum branch not found."
            << std::endl;
        return 0;
    }

    Int_t resnum = -1;

    T_eventweight->SetBranchStatus("*", 0);
    T_eventweight->SetBranchStatus("GTruth_ResNum", 1);
    T_eventweight->SetBranchAddress("GTruth_ResNum", &resnum);

    Long64_t n_delta1600 = 0;

    for (Long64_t i = 0; i < T_eventweight->GetEntries(); ++i)
    {
        T_eventweight->GetEntry(i);

        if (resnum == 9)
            ++n_delta1600;

    }

    return n_delta1600;
}

int main( int argc, char** argv )
{

  TString input_filename = argv[1];
  TString out_filename = argv[2];


  bool flag_data = true;

  bool flag_osc = false;
  for (Int_t i=1;i!=argc;i++){
    switch(argv[i][1]){
    case 'o':
      flag_osc = atoi(&argv[i][2]); // run oscillation
      break;
    }
  }
  
  
  TFile *file = new TFile(input_filename,"READ");
  
  CovMatrix cov;

  std::map<TString, std::tuple<int, int, TString, float, int, double, int> > map_inputfile_info = cov.get_map_inputfile_info();

  auto it = map_inputfile_info.find(input_filename);

  if (it == map_inputfile_info.end())
  {
    std::cerr << "ERROR: Cannot find input file in cv_input.txt: " << input_filename << std::endl;
    return 1;
  }

  int filetype = std::get<0>(it->second);
  int period   = std::get<1>(it->second);

  
  // ============================================================
  // NUE flux sampling reweighting
  // ============================================================

  // Only nue (filetype=1) is reweighted.
  // Run4A nue (period=4) keeps weight = 1.
  // bool is_intrinsic_nue = (filetype == 1);
  // bool is_run4a = (period == 4);

  // bool is_intrinsic_nue = input_filename.Contains("nue_overlay");
  // bool is_nu_overlay = input_filename.Contains("nu_overlay");
  // bool is_SURPRISE = input_filename.Contains("MCC9.10");
  // bool is_run4a = input_filename.Contains("run4a");
  // bool is_run1 = input_filename.Contains("run1");
  // bool is_run2 = input_filename.Contains("run2");
  // bool is_run3 = input_filename.Contains("run3");
  // bool is_run4b = input_filename.Contains("run4b");
  // bool is_run4c = input_filename.Contains("run4c");
  // bool is_run4d = input_filename.Contains("run4d");
  // bool is_run5 = input_filename.Contains("run5");
  // bool is_run123 = input_filename.Contains("run1") || input_filename.Contains("run2") || input_filename.Contains("run3");
  bool is_intrinsic_nue = (filetype == 1);
  bool is_nu_overlay = (filetype == 2);
  bool is_SURPRISE = input_filename.Contains("MCC9.10");
  bool is_run4a = (period == 4);
  bool is_run1 = (period == 1);
  bool is_run2 = (period == 2);
  bool is_run3 = (period == 3);
  bool is_run4b = (period == 42);
  bool is_run4c = (period == 43);
  bool is_run4d = (period == 44);
  bool is_run5 = (period == 5);
  bool is_run123 = (period == 1) || (period == 2) || (period == 3);
  bool use_nue_flux_reweight = (USE_NUE_FLUX_REWEIGHT && is_intrinsic_nue && !is_run4a);
  bool use_rootino_reweight = (USE_ROOTINO_REWEIGHT && is_SURPRISE && is_nu_overlay && is_run123);


  if (flag_osc) cov.add_osc_config();
  
  cov.print_rw(cov.get_rw_info());
  
  TTree *T_BDTvars = (TTree*)file->Get("wcpselection/T_BDTvars");
  TTree *T_eval = (TTree*)file->Get("wcpselection/T_eval");
  TTree *T_pot = (TTree*)file->Get("wcpselection/T_pot");
  TTree *T_PFeval = (TTree*)file->Get("wcpselection/T_PFeval");
  TTree *T_KINEvars = (TTree*)file->Get("wcpselection/T_KINEvars");

  // ==============================================
  // ROOTino eventweight tree
  // ==============================================
  TTree *T_eventweight = nullptr;

  if (use_rootino_reweight)
  {
    TDirectory *singlephotonana = dynamic_cast<TDirectory*>(file->Get("singlephotonana"));

    if (singlephotonana)
    {
        T_eventweight = dynamic_cast<TTree*>(singlephotonana->Get("eventweight_tree"));
    }

    if (!T_eventweight)
    {
        std::cerr
            << "[ROOTINO] ERROR: cannot find "
            << "singlephotonana/eventweight_tree"
            << " in SURPRISE file:"
            << input_filename
            << std::endl;

        return 1;
    }
  }
  
 
  if (T_eval->GetBranch("weight_cv")) flag_data = false;

  EvalInfo eval;
  POTInfo pot;
  TaggerInfo tagger;
  PFevalInfo pfeval;
  KineInfo kine;
  // ============================================================
  // ROOTino truth information
  // ============================================================
  Int_t GTruth_ResNum = -1;

  if (T_eventweight)
  {
    if (!T_eventweight->GetBranch("GTruth_ResNum"))
    {
        std::cerr
            << "[ROOTINO] ERROR: GTruth_ResNum branch "
            << "not found in eventweight_tree."
            << std::endl;

        return 1;
    }

    T_eventweight->SetBranchStatus("*", 0);
    T_eventweight->SetBranchStatus("GTruth_ResNum", 1);
    T_eventweight->SetBranchAddress("GTruth_ResNum", &GTruth_ResNum);
  }


#include "init.txt"

  set_tree_address(T_BDTvars, tagger,2 );
  if (flag_data){
    set_tree_address(T_eval, eval,2);
    set_tree_address(T_PFeval, pfeval,2);
  }else{
    set_tree_address(T_eval, eval);
    set_tree_address(T_PFeval, pfeval);
  }
  set_tree_address(T_pot, pot);
  set_tree_address(T_KINEvars, kine);


  std::vector<double> nue_flux_target;
  std::vector<double> nue_flux_target_unweight;
  std::vector<double> nue_flux_weight;

  double nue_flux_bin_width = 0.05; // 50 MeV
  double nue_flux_file_POT = 0.0;

  TString nue_flux_csv = "";
  std::vector<TString> possible_csv_paths;

  // possible_csv_paths.push_back("/home/lsy/workspace/LEEana/configurations/enu_reweighted_run4b4d_50MeV.csv");
  // possible_csv_paths.push_back("/nashome/s/shuyiliu/workspace/LEEana/configurations/enu_reweighted_run4b4d_50MeV.csv");
  possible_csv_paths.push_back("../LEEana/configurations/enu_reweighted_run4b4d_50MeV.csv");

  for (auto& path : possible_csv_paths) 
  {
    std::ifstream test(path.Data());
    if (test.good()) 
    {
      nue_flux_csv = path;
      test.close();
      break;
    }
    test.close();
  }

  if (use_nue_flux_reweight && nue_flux_csv == "") 
  {
    std::cerr << "[NUE FLUX] ERROR: " "cannot find " "enu_reweighted_run4b4d_50MeV.csv" << std::endl;
    return 1;
  }

  if (use_nue_flux_reweight) 
  {
    if (!load_nue_flux_target(nue_flux_csv, nue_flux_target, nue_flux_target_unweight, nue_flux_bin_width)) 
    {
      return 1;
    }
  }


  double total_pot = 0;
  for (Int_t i=0;i!=T_pot->GetEntries();i++){
    T_pot->GetEntry(i);
    total_pot += pot.pot_tor875;
  }
  double ext_pot = cov.get_ext_pot(input_filename);
  if (ext_pot != 0) total_pot = ext_pot;
  
  std::cout << "Total POT: " << total_pot << " external POT: " << ext_pot << std::endl;

  

  double ROOTino_current_POT = 0.0;
  Long64_t ROOTino_current_Delta1600 = 0;
  double ROOTino_POT45 = 0.0;

  Long64_t ROOTino_Delta4b = 0;
  Long64_t ROOTino_Delta4c = 0;
  Long64_t ROOTino_Delta4d = 0;
  Long64_t ROOTino_Delta5  = 0;

  Long64_t ROOTino_Delta45 = 0;

  double ROOTino_target_Delta1600 = 0.0;
  double ROOTino_weight = 1.0;

  if (use_rootino_reweight)
  {
    
    ROOTino_current_POT = 0.0;

    for (Long64_t i = 0; i < T_pot->GetEntries(); ++i)
    {
        T_pot->GetEntry(i);
        ROOTino_current_POT += pot.pot_tor875;
    }

    if (ROOTino_current_POT <= 0.0)
    {
        std::cerr
            << "[ROOTINO] ERROR: current input POT <= 0.\n"
            << "Input file = "
            << input_filename
            << std::endl;

        return 1;
    }

    ROOTino_current_Delta1600 = count_ROOTino_Delta1600(T_eventweight);
    if (T_eventweight)
    {
      T_eventweight->SetBranchStatus("*", 0);
      T_eventweight->SetBranchStatus("GTruth_ResNum", 1);
      T_eventweight->SetBranchAddress("GTruth_ResNum", &GTruth_ResNum);
    }


    TFile *rootino_file4b = TFile::Open(ROOTINO_RUN4B_FILE, "READ");
    TFile *rootino_file4c = TFile::Open(ROOTINO_RUN4C_FILE, "READ");
    TFile *rootino_file4d = TFile::Open(ROOTINO_RUN4D_FILE, "READ");
    TFile *rootino_file5 = TFile::Open(ROOTINO_RUN5_FILE, "READ");

    if (!rootino_file4b || rootino_file4b->IsZombie())
    {
        std::cerr
            << "[ROOTINO] ERROR opening Run4b:\n"
            << ROOTINO_RUN4B_FILE
            << std::endl;

        return 1;
    }

    if (!rootino_file4c || rootino_file4c->IsZombie())
    {
        std::cerr
            << "[ROOTINO] ERROR opening Run4c:\n"
            << ROOTINO_RUN4C_FILE
            << std::endl;

        return 1;
    }

    if (!rootino_file4d || rootino_file4d->IsZombie())
    {
        std::cerr
            << "[ROOTINO] ERROR opening Run4d:\n"
            << ROOTINO_RUN4D_FILE
            << std::endl;

        return 1;
    }

    if (!rootino_file5 || rootino_file5->IsZombie())
    {
        std::cerr
            << "[ROOTINO] ERROR opening Run5:\n"
            << ROOTINO_RUN5_FILE
            << std::endl;

        return 1;
    }

    TTree *R4b_pot = (TTree*)rootino_file4b->Get("wcpselection/T_pot");
    TTree *R4c_pot = (TTree*)rootino_file4c->Get("wcpselection/T_pot");
    TTree *R4d_pot = (TTree*)rootino_file4d->Get("wcpselection/T_pot");
    TTree *R5_pot = (TTree*)rootino_file5->Get("wcpselection/T_pot");


    if (!R4b_pot || !R4c_pot || !R4d_pot || !R5_pot)
    {
        std::cerr
            << "[ROOTINO] ERROR: cannot find "
            << "T_pot in Run4/5 reference files."
            << std::endl;

        return 1;
    }

    POTInfo R4b_pot_info;
    R4b_pot->SetBranchAddress("pot_tor875", &R4b_pot_info.pot_tor875);
    double POT4b = 0.0;
    for (Long64_t i = 0; i < R4b_pot->GetEntries(); ++i)
    {
        R4b_pot->GetEntry(i);
        POT4b += R4b_pot_info.pot_tor875;
    }

    POTInfo R4c_pot_info;
    R4c_pot->SetBranchAddress("pot_tor875", &R4c_pot_info.pot_tor875);
    double POT4c = 0.0;
    for (Long64_t i = 0; i < R4c_pot->GetEntries(); ++i)
    {
        R4c_pot->GetEntry(i);
        POT4c += R4c_pot_info.pot_tor875;
    }

    POTInfo R4d_pot_info;
    R4d_pot->SetBranchAddress("pot_tor875", &R4d_pot_info.pot_tor875);
    double POT4d = 0.0;
    for (Long64_t i = 0; i < R4d_pot->GetEntries(); ++i)
    {
        R4d_pot->GetEntry(i);
        POT4d += R4d_pot_info.pot_tor875;
    }

    POTInfo R5_pot_info;
    R5_pot->SetBranchAddress("pot_tor875", &R5_pot_info.pot_tor875);
    double POT5 = 0.0;
    for (Long64_t i = 0; i < R5_pot->GetEntries(); ++i)
    {
        R5_pot->GetEntry(i);
        POT5 += R5_pot_info.pot_tor875;
    }

    ROOTino_POT45 = POT4b + POT4c + POT4d + POT5;

    if (ROOTino_POT45 <= 0.0)
    {
        std::cerr
            << "[ROOTINO] ERROR: POT45 <= 0."
            << std::endl;

        return 1;
    }

    TDirectory *R4b_singlephotonana = dynamic_cast<TDirectory*>(rootino_file4b->Get("singlephotonana"));
    TDirectory *R4c_singlephotonana = dynamic_cast<TDirectory*>(rootino_file4c->Get("singlephotonana"));
    TDirectory *R4d_singlephotonana = dynamic_cast<TDirectory*>(rootino_file4d->Get("singlephotonana"));
    TDirectory *R5_singlephotonana = dynamic_cast<TDirectory*>(rootino_file5->Get("singlephotonana"));

    TTree *R4b_eventweight = nullptr;
    TTree *R4c_eventweight = nullptr;
    TTree *R4d_eventweight = nullptr;
    TTree *R5_eventweight  = nullptr;


    if (R4b_singlephotonana)
    {
        R4b_eventweight = dynamic_cast<TTree*>(R4b_singlephotonana->Get("eventweight_tree"));
    }

    if (R4c_singlephotonana)
    {
        R4c_eventweight = dynamic_cast<TTree*>(R4c_singlephotonana->Get("eventweight_tree"));
    }

    if (R4d_singlephotonana)
    {
        R4d_eventweight = dynamic_cast<TTree*>(R4d_singlephotonana->Get("eventweight_tree"));
    }

    if (R5_singlephotonana)
    {
        R5_eventweight = dynamic_cast<TTree*>(R5_singlephotonana->Get("eventweight_tree"));
    }


    if (!R4b_eventweight || !R4c_eventweight || !R4d_eventweight || !R5_eventweight)
    {
        std::cerr
            << "[ROOTINO] ERROR: cannot find "
            << "eventweight_tree in Run4/5."
            << std::endl;

        return 1;
    }

    ROOTino_Delta4b = count_ROOTino_Delta1600(R4b_eventweight);
    ROOTino_Delta4c = count_ROOTino_Delta1600(R4c_eventweight);
    ROOTino_Delta4d = count_ROOTino_Delta1600(R4d_eventweight);
    ROOTino_Delta5 = count_ROOTino_Delta1600(R5_eventweight);

    ROOTino_Delta45 = ROOTino_Delta4b + ROOTino_Delta4c + ROOTino_Delta4d + ROOTino_Delta5;

    ROOTino_target_Delta1600 = (ROOTino_current_POT / ROOTino_POT45) * static_cast<double>(ROOTino_Delta45);

    if (ROOTino_current_Delta1600 > 0)
    {
        ROOTino_weight = ROOTino_target_Delta1600 / static_cast<double>(ROOTino_current_Delta1600);
    }
    else
    {
        ROOTino_weight = 1.0;
        std::cerr
            << "[ROOTINO] WARNING: current "
            << "Delta(1600) count = 0. "
            << "ROOTino weight set to 1."
            << std::endl;
    }
  }


  // prepare histograms ...
  // declare histograms ...
  TH1F *htemp;
  std::map<TString, TH1F*> map_histoname_hist;
  std::vector< std::tuple<TString,  int, float, float, TString, TString, TString, TString > > all_histo_infos;
  
  std::vector< std::tuple<TString,  int, float, float, TString, TString, TString, TString > > histo_infos = cov.get_histograms(input_filename,0);
  std::copy(histo_infos.begin(), histo_infos.end(), std::back_inserter(all_histo_infos));
  //  std::cout << "CV:" << std::endl;
  for (auto it = histo_infos.begin(); it != histo_infos.end(); it++){
    TString histoname = std::get<0>(*it);
    Int_t nbin = std::get<1>(*it);
    float llimit = std::get<2>(*it);
    float hlimit = std::get<3>(*it);
    TString var_name = std::get<4>(*it);
    TString ch_name = std::get<5>(*it);
    TString add_cut = std::get<6>(*it);
    TString weight = std::get<7>(*it);
    
    //    std::cout << std::get<0>( *it)  << " " << std::get<1>(*it) << " " << std::get<4>(*it) << " " << std::get<5>(*it) << " " << std::get<6>(*it) << " " << std::get<7>(*it) << std::endl;
    htemp = new TH1F(histoname, histoname, nbin, llimit, hlimit);
    map_histoname_hist[histoname] = htemp;
  }
  //  std::cout << std::endl;

  std::vector< std::tuple<TString,  int, float, float, TString, TString, TString, TString > > histo_infos_err2 = cov.get_histograms(input_filename,1);
  std::copy(histo_infos_err2.begin(), histo_infos_err2.end(), std::back_inserter(all_histo_infos));
  //  std::cout << "Error2: " << std::endl;
  for (auto it = histo_infos_err2.begin(); it != histo_infos_err2.end(); it++){
    TString histoname = std::get<0>(*it);
    Int_t nbin = std::get<1>(*it);
    float llimit = std::get<2>(*it);
    float hlimit = std::get<3>(*it);
    TString var_name = std::get<4>(*it);
    TString ch_name = std::get<5>(*it);
    TString add_cut = std::get<6>(*it);
    TString weight = std::get<7>(*it);
    //std::cout << std::get<0>( *it) << " " << std::get<1>(*it) << " " << std::get<4>(*it) << " " << std::get<5>(*it) << " " << std::get<6>(*it) << " " << std::get<7>(*it) << std::endl;
    htemp = new TH1F(histoname, histoname, nbin, llimit, hlimit);
    map_histoname_hist[histoname] = htemp;
  }
  //  std::cout << std::endl;

  std::vector< std::tuple<TString,  int, float, float, TString, TString, TString, TString > > histo_infos_cros = cov.get_histograms(input_filename,2);
  std::copy(histo_infos_cros.begin(), histo_infos_cros.end(), std::back_inserter(all_histo_infos));
  //std::cout << "Cross: " << std::endl;
  for (auto it = histo_infos_cros.begin(); it != histo_infos_cros.end(); it++){
    TString histoname = std::get<0>(*it);
    Int_t nbin = std::get<1>(*it);
    float llimit = std::get<2>(*it);
    float hlimit = std::get<3>(*it);
    TString var_name = std::get<4>(*it);
    TString ch_name = std::get<5>(*it);
    TString add_cut = std::get<6>(*it);
    TString weight = std::get<7>(*it);
    //std::cout << std::get<0>( *it) << " " << std::get<1>(*it) << " " << std::get<4>(*it) << " " << std::get<5>(*it) << " " << std::get<6>(*it) << " " << std::get<7>(*it) << std::endl;
    htemp = new TH1F(histoname, histoname, nbin, llimit, hlimit);
    map_histoname_hist[histoname] = htemp;
  }
  //  std::cout << std::endl;
  
  // fill histogram ...
  T_BDTvars->SetBranchStatus("*",0);
  T_BDTvars->SetBranchStatus("numu_cc_flag",1);
  T_BDTvars->SetBranchStatus("numu_score",1);
  T_BDTvars->SetBranchStatus("nue_score",1);
  T_BDTvars->SetBranchStatus("cosmict_flag",1);
  T_BDTvars->SetBranchStatus("mip_vec_dQ_dx_0",1);
  T_BDTvars->SetBranchStatus("mip_vec_dQ_dx_1",1);
  T_BDTvars->SetBranchStatus("mip_vec_dQ_dx_2",1);
  T_BDTvars->SetBranchStatus("mip_vec_dQ_dx_3",1);
  T_BDTvars->SetBranchStatus("mip_vec_dQ_dx_4",1);
  T_BDTvars->SetBranchStatus("mip_vec_dQ_dx_5",1);
  T_BDTvars->SetBranchStatus("mip_vec_dQ_dx_6",1);
  T_BDTvars->SetBranchStatus("mip_vec_dQ_dx_7",1);
  T_BDTvars->SetBranchStatus("mip_vec_dQ_dx_8",1);
  T_BDTvars->SetBranchStatus("mip_vec_dQ_dx_9",1);
  T_BDTvars->SetBranchStatus("mip_vec_dQ_dx_10",1);
  T_BDTvars->SetBranchStatus("mip_vec_dQ_dx_11",1);
  T_BDTvars->SetBranchStatus("mip_vec_dQ_dx_12",1);
  T_BDTvars->SetBranchStatus("mip_vec_dQ_dx_13",1);
  T_BDTvars->SetBranchStatus("mip_vec_dQ_dx_14",1);
  T_BDTvars->SetBranchStatus("mip_vec_dQ_dx_15",1);
  T_BDTvars->SetBranchStatus("mip_vec_dQ_dx_16",1);
  T_BDTvars->SetBranchStatus("mip_vec_dQ_dx_17",1);
  T_BDTvars->SetBranchStatus("mip_vec_dQ_dx_18",1);
  T_BDTvars->SetBranchStatus("mip_vec_dQ_dx_19",1);
  T_BDTvars->SetBranchStatus("mip_energy",1);
  T_BDTvars->SetBranchStatus("mip_angle_beam",1);
  T_BDTvars->SetBranchStatus("spt_angle_vertical",1);
  T_BDTvars->SetBranchStatus("mip_quality_n_tracks",1);
  T_BDTvars->SetBranchStatus("mip_quality_n_showers",1);
  T_BDTvars->SetBranchStatus("gap_n_bad",1);
  T_BDTvars->SetBranchStatus("spt_angle_beam",1);
  T_BDTvars->SetBranchStatus("spt_angle_vertical",1);

  if (tagger.flag_nc_gamma_bdt){
    T_BDTvars->SetBranchStatus("nc_delta_score", 1);
    T_BDTvars->SetBranchStatus("nc_pio_score", 1);
  }
  
  T_eval->SetBranchStatus("*",0);
  T_eval->SetBranchStatus("match_energy",1);
  T_eval->SetBranchStatus("match_isFC",1);
  T_eval->SetBranchStatus("match_found",1);
  if (T_eval->GetBranch("match_found_asInt")) T_eval->SetBranchStatus("match_found_asInt",1); 
  T_eval->SetBranchStatus("stm_eventtype",1);
  T_eval->SetBranchStatus("stm_lowenergy",1);
  T_eval->SetBranchStatus("stm_LM",1);
  T_eval->SetBranchStatus("stm_TGM",1);
  T_eval->SetBranchStatus("stm_STM",1);
  T_eval->SetBranchStatus("stm_FullDead",1);
  T_eval->SetBranchStatus("stm_clusterlength",1);
  
  if (!flag_data){
    T_eval->SetBranchStatus("weight_spline",1);
    T_eval->SetBranchStatus("weight_cv",1);
    T_eval->SetBranchStatus("weight_lee",1);
    T_eval->SetBranchStatus("weight_change",1);
    // MC enable truth information ...
    T_eval->SetBranchStatus("truth_isCC",1);
    T_eval->SetBranchStatus("truth_nuPdg",1);
    T_eval->SetBranchStatus("truth_vtxInside",1);
    T_eval->SetBranchStatus("truth_nuEnergy",1);
    T_eval->SetBranchStatus("truth_energyInside",1);
    T_eval->SetBranchStatus("truth_vtxX",1);
    T_eval->SetBranchStatus("truth_vtxY",1);
    T_eval->SetBranchStatus("truth_vtxZ",1);
    T_eval->SetBranchStatus("match_completeness_energy",1);
  }

  
  T_KINEvars->SetBranchStatus("*",0);
  T_KINEvars->SetBranchStatus("kine_reco_Enu",1);
  T_KINEvars->SetBranchStatus("kine_energy_particle",1);
  T_KINEvars->SetBranchStatus("kine_particle_type",1);
  T_KINEvars->SetBranchStatus("kine_energy_info",1);
  T_KINEvars->SetBranchStatus("kine_energy_included",1);
  T_KINEvars->SetBranchStatus("kine_reco_add_energy",1);
  T_KINEvars->SetBranchStatus("kine_pio_mass",1);
  T_KINEvars->SetBranchStatus("kine_pio_flag",1);
  T_KINEvars->SetBranchStatus("kine_pio_vtx_dis",1);
  T_KINEvars->SetBranchStatus("kine_pio_energy_1",1);
  T_KINEvars->SetBranchStatus("kine_pio_theta_1",1);
  T_KINEvars->SetBranchStatus("kine_pio_phi_1",1);
  T_KINEvars->SetBranchStatus("kine_pio_dis_1",1);
  T_KINEvars->SetBranchStatus("kine_pio_energy_2",1);
  T_KINEvars->SetBranchStatus("kine_pio_theta_2",1);
  T_KINEvars->SetBranchStatus("kine_pio_phi_2",1);
  T_KINEvars->SetBranchStatus("kine_pio_dis_2",1);
  T_KINEvars->SetBranchStatus("kine_pio_angle",1);
  if (T_KINEvars->GetBranch("vlne_v4_numu_full_primaryE")) {
    T_KINEvars->SetBranchStatus("vlne_v4_numu_full_primaryE",1);
    T_KINEvars->SetBranchStatus("vlne_v4_numu_full_totalE",1);
    T_KINEvars->SetBranchStatus("vlne_v4_numu_partial_primaryE",1);
    T_KINEvars->SetBranchStatus("vlne_v4_numu_partial_totalE",1);
    // T_KINEvars->SetBranchStatus("vlne_nue_full_primaryE",1);
    // T_KINEvars->SetBranchStatus("vlne_nue_full_totalE",1);
    // T_KINEvars->SetBranchStatus("vlne_nue_partial_primaryE",1);
    // T_KINEvars->SetBranchStatus("vlne_nue_partial_totalE",1);
  }


  T_PFeval->SetBranchStatus("*",0);
  T_PFeval->SetBranchStatus("reco_nuvtxX",1);
  T_PFeval->SetBranchStatus("reco_nuvtxY",1);
  T_PFeval->SetBranchStatus("reco_nuvtxZ",1);
  T_PFeval->SetBranchStatus("reco_showervtxX",1);
  T_PFeval->SetBranchStatus("reco_showervtxY",1);
  T_PFeval->SetBranchStatus("reco_showervtxZ",1);
  T_PFeval->SetBranchStatus("reco_muonMomentum",1);
  T_PFeval->SetBranchStatus("reco_showerKE",1);
  if (!flag_data){
      T_PFeval->SetBranchStatus("nuvtx_diff",1);
      T_PFeval->SetBranchStatus("showervtx_diff",1);
      T_PFeval->SetBranchStatus("muonvtx_diff",1);
      T_PFeval->SetBranchStatus("truth_nuIntType",1);
      T_PFeval->SetBranchStatus("truth_muonMomentum",1);
      if(T_PFeval->GetBranch("truth_mother")){//prevents throwing an error for the non _PF files
        T_PFeval->SetBranchStatus("truth_Ntrack",1); 
        T_PFeval->SetBranchStatus("truth_pdg",1); 
        T_PFeval->SetBranchStatus("truth_mother",1); 
        T_PFeval->SetBranchStatus("truth_startMomentum",1); 
      }
      
  }
  if (pfeval.flag_NCDelta){
    
      if (!flag_data){
          T_PFeval->SetBranchStatus("truth_NCDelta",1);
          T_PFeval->SetBranchStatus("truth_NprimPio",1);
      }
  }
  
  if (pfeval.flag_recoprotonMomentum){
    T_PFeval->SetBranchStatus("reco_protonMomentum",1);
  }
  
  if (pfeval.flag_showerMomentum){
    T_PFeval->SetBranchStatus("reco_showerMomentum",1);
    T_PFeval->SetBranchStatus("reco_Nproton",1);
    if (!flag_data){
      T_PFeval->SetBranchStatus("truth_showerMomentum",1);
      T_PFeval->SetBranchStatus("truth_nuScatType",1);
      // oscillation formula ...
      T_PFeval->SetBranchStatus("truth_nu_momentum",1);
      T_PFeval->SetBranchStatus("neutrino_type",1);
      T_PFeval->SetBranchStatus("mcflux_ntype",1); // flux flavor before oscillation
      T_PFeval->SetBranchStatus("mcflux_dk2gen",1);
      T_PFeval->SetBranchStatus("mcflux_gen2vtx",1);
      T_PFeval->SetBranchStatus("mcflux_ndecay",1);
    }
  }

  std::cout << "Total entries: " << T_eval->GetEntries() << std::endl;

  if (use_nue_flux_reweight) 
  {
    if (is_run4a) 
    {
      std::cout << "nue flux sampling weight: run 4a had a different flux, skipping reweighting (all weights set to 1)." << std::endl;

    } else 
    {
      if (!build_nue_flux_weights(T_eval, T_pot, eval, pot, nue_flux_target, nue_flux_target_unweight, nue_flux_bin_width, nue_flux_weight, nue_flux_file_POT, is_run123)) 
      {
        return 1;
      }
    }
}


  static int debug_weight_count = 0;

  for (Int_t i=0;i!=T_eval->GetEntries();i++)
  {
    T_BDTvars->GetEntry(i);
    T_eval->GetEntry(i);
    T_KINEvars->GetEntry(i);
    T_PFeval->GetEntry(i);
    
    if (!is_preselection(eval)) continue;

    for (auto it = all_histo_infos.begin(); it != all_histo_infos.end(); it++)
    {
      TString histoname = std::get<0>(*it);
      Int_t nbin = std::get<1>(*it);
      float llimit = std::get<2>(*it);
      float hlimit = std::get<3>(*it);
      TString var_name = std::get<4>(*it);
      TString ch_name = std::get<5>(*it);
      TString add_cut = std::get<6>(*it);
      TString weight = std::get<7>(*it);
 
      htemp = map_histoname_hist[histoname];
      // get kinematics variable ...
      double val = get_kine_var(kine, eval, pfeval, tagger, flag_data, var_name);
      // get pass or not
      bool flag_pass = get_cut_pass(ch_name, add_cut, flag_data, eval, pfeval, tagger, kine);

      double osc_weight = 1.0;

     
      
      // std::cout << weight << std::endl;
      // get weight ...
      // double weight_val = get_weight(weight, eval, pfeval, kine, tagger, cov.get_rw_info(), flag_data);

    double original_weight = get_weight(weight, eval, pfeval, kine, tagger, cov.get_rw_info(), flag_data);
    double nue_weight = 1.0;

    if (use_nue_flux_reweight)
    {
      if (is_run4a)
      {
        nue_weight = 1.0;
      }
      else
      {
        nue_weight = get_nue_flux_weight(eval.truth_nuEnergy, nue_flux_weight, nue_flux_bin_width);
      }
    }

    // ============================================================
    // ROOTino reweight
    // ============================================================

    double ROOTino_reweight = 1.0;

    if (use_rootino_reweight)
    {
      T_eventweight->GetEntry(i);

      // GTruth_ResNum == 9 -> Delta(1600)
      bool is_Delta1600 = (GTruth_ResNum == 9);
      if (is_Delta1600)
      {
        ROOTino_reweight = ROOTino_weight;
      }
      else
      {
        ROOTino_reweight = 1.0;
      }
    }

    double final_weight = original_weight * nue_weight * ROOTino_reweight;
    double weight_val = final_weight;

  static Long64_t total_nue_reweighted_events = 0;
  static Long64_t zero_nue_weight_events = 0;
  static Long64_t out_of_range_nue_events = 0;

  if (use_nue_flux_reweight && !is_run4a)
  {
    total_nue_reweighted_events++;

    double enu_GeV = eval.truth_nuEnergy / 1000.0;

    int bin = static_cast<int>(std::floor(enu_GeV / nue_flux_bin_width));

    nue_weight = get_nue_flux_weight(eval.truth_nuEnergy, nue_flux_weight, nue_flux_bin_width);

    if (bin < 0 || bin >= (int)nue_flux_weight.size())
    {
        out_of_range_nue_events++;
    }

    if (nue_weight == 0.0)
    {
        zero_nue_weight_events++;
    }
  }

  
      if (flag_osc && cov.is_osc_channel(ch_name) && (!flag_data)){
	osc_weight = cov.get_osc_weight(eval, pfeval);
	weight_val *= osc_weight;
	if (weight == "cv_spline_cv_spline" || weight == "unity_unity" ||
	    weight == "spline_spline" )
	  weight_val *= osc_weight;
      }
      
      if (flag_pass)
	htemp->Fill(val,weight_val);
    }
  }
  

  // save histograms ...
  TFile *file1 = new TFile(out_filename,"RECREATE");
  file1->cd();
  TTree *T = new TTree("T","T");
  T->SetDirectory(file1);
  T->Branch("pot",&total_pot,"pot/D");
  T->Fill();
  
  for (auto it = map_histoname_hist.begin(); it!= map_histoname_hist.end(); it++){
    //std::cout<<"DEBUG: "<<it->first<<" "<<it->second->GetName()<<" "<<it->second->GetSum()<<"\n";
    it->second->SetDirectory(file1);
  }
  
  file1->Write();
  file1->Close();

  return 0;
}
