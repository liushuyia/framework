// modified from https://raw.githubusercontent.com/HaiwangYu/wcp-mc-eval/main/muon_energy.C
#include "TFile.h"
#include "TTree.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TCanvas.h"
#include "TLorentzVector.h"
#include "TEfficiency.h"

#include <vector>
#include <iostream>

using namespace std;


enum LIMITS {
    MAX_TRACKS = 30000,
};



void plot_nc_rez(){
    gInterpreter->GenerateDictionary("vector<vector<int> >", "vector");
    gInterpreter->GenerateDictionary("vector<string >", "vector");

    auto *tf = TFile::Open("checkout_prodgenie_bnb_nu_overlay_run123all_PF.root", "read");
    // auto *tf = TFile::Open("/data1/wgu/processed_checkout_rootfiles/checkout_prodgenie_bnb_nu_overlay_run123all_PF.root", "read");
    auto *dir = (TDirectoryFile *) tf->Get("wcpselection");

    auto *T_PFDump = (TTree *) dir->Get("T_PFeval");
    auto *T_BDTvars = (TTree *) dir->Get("T_BDTvars");
    auto *T_eval = (TTree *) dir->Get("T_eval");
    auto *T_KINEvars = (TTree *) dir->Get("T_KINEvars");
    // T_PFDump->AddFriend(T_BDTvars);
    // T_PFDump->AddFriend(T_eval);
    // T_PFDump->AddFriend(T_KINEvars);

    const int target_pdg = 13;  // 11, 13, 2212

    int truth_Ntrack;
    int truth_id[MAX_TRACKS];
    int truth_pdg[MAX_TRACKS];
    int truth_mother[MAX_TRACKS];
    float truth_startXYZT[MAX_TRACKS][4];
    float truth_endXYZT[MAX_TRACKS][4];
    float truth_startMomentum[MAX_TRACKS][4];
    float truth_endMomentum[MAX_TRACKS][4];
    std::vector<std::vector<int> > *truth_daughters = 0;
    std::vector<string > *truth_process = 0;

    int reco_Ntrack;
    int reco_id[MAX_TRACKS];
    int reco_pdg[MAX_TRACKS];
    int reco_process[MAX_TRACKS];
    int reco_mother[MAX_TRACKS];
    float reco_startXYZT[MAX_TRACKS][4];
    float reco_endXYZT[MAX_TRACKS][4];
    float reco_startMomentum[MAX_TRACKS][4];
    float reco_endMomentum[MAX_TRACKS][4];
    std::vector<std::vector<int> > *reco_daughters = 0;

    float truth_nuEnergy;
    float truth_energyInside;
    bool match_isFC;
    float match_completeness_energy;

    int truth_nuPdg;
    bool truth_isCC;
    float truth_nu_momentum[4];
    float truth_muonMomentum[4];
    float truth_corr_nuvtxX;
    float truth_corr_nuvtxY;
    float truth_corr_nuvtxZ;
    float truth_muonvtxX, truth_muonvtxY, truth_muonvtxZ;
    float truth_muonendX, truth_muonendY, truth_muonendZ;
    float reco_nuvtxX;
    float reco_nuvtxY;
    float reco_nuvtxZ;
    float reco_muonMomentum[4];

    float numu_cc_flag;
    float numu_score;
    float nue_score;
    float cosmict_flag;

    float stm_clusterlength;
    bool truth_isFC;
    bool truth_vtxInside;

    float kine_reco_Enu;
    float kine_reco_add_energy;
    std::vector<float> *kine_energy_particle = 0;
    std::vector<int> *kine_energy_info = 0;
    std::vector<int> *kine_particle_type = 0;
    int kine_pio_flag;
    float kine_pio_vtx_dis;
    float kine_pio_energy_1;
    float kine_pio_energy_2;
    float kine_pio_dis_1;
    float kine_pio_dis_2;
    float kine_pio_angle;
    float kine_pio_mass;

    T_PFDump->SetBranchAddress("truth_Ntrack", &truth_Ntrack);
    T_PFDump->SetBranchAddress("truth_id", &truth_id);
    T_PFDump->SetBranchAddress("truth_pdg", &truth_pdg);
    T_PFDump->SetBranchAddress("truth_mother", &truth_mother);
    T_PFDump->SetBranchAddress("truth_startXYZT", &truth_startXYZT);
    T_PFDump->SetBranchAddress("truth_endXYZT", &truth_endXYZT);
    T_PFDump->SetBranchAddress("truth_startMomentum", &truth_startMomentum);
    T_PFDump->SetBranchAddress("truth_endMomentum", &truth_endMomentum);
    T_PFDump->SetBranchAddress("truth_daughters", &truth_daughters);
    T_PFDump->SetBranchAddress("truth_process", &truth_process);

    T_PFDump->SetBranchAddress("reco_Ntrack", &reco_Ntrack);
    T_PFDump->SetBranchAddress("reco_id", &reco_id);
    T_PFDump->SetBranchAddress("reco_pdg", &reco_pdg);
    T_PFDump->SetBranchAddress("reco_mother", &reco_mother);
    T_PFDump->SetBranchAddress("reco_startXYZT", &reco_startXYZT);
    T_PFDump->SetBranchAddress("reco_endXYZT", &reco_endXYZT);
    T_PFDump->SetBranchAddress("reco_startMomentum", &reco_startMomentum);
    T_PFDump->SetBranchAddress("reco_endMomentum", &reco_endMomentum);
    T_PFDump->SetBranchAddress("reco_daughters", &reco_daughters);


    T_PFDump->SetBranchAddress("truth_nu_momentum", &truth_nu_momentum);
    T_PFDump->SetBranchAddress("truth_muonMomentum", &truth_muonMomentum);
    T_PFDump->SetBranchAddress("truth_corr_nuvtxX", &truth_corr_nuvtxX);
    T_PFDump->SetBranchAddress("truth_corr_nuvtxY", &truth_corr_nuvtxY);
    T_PFDump->SetBranchAddress("truth_corr_nuvtxZ", &truth_corr_nuvtxZ);
    T_PFDump->SetBranchAddress("truth_muonvtxX", &truth_muonvtxX);
    T_PFDump->SetBranchAddress("truth_muonvtxY", &truth_muonvtxY);
    T_PFDump->SetBranchAddress("truth_muonvtxZ", &truth_muonvtxZ);
    T_PFDump->SetBranchAddress("truth_muonendX", &truth_muonendX);
    T_PFDump->SetBranchAddress("truth_muonendY", &truth_muonendY);
    T_PFDump->SetBranchAddress("truth_muonendZ", &truth_muonendZ);
    T_PFDump->SetBranchAddress("reco_nuvtxX", &reco_nuvtxX);
    T_PFDump->SetBranchAddress("reco_nuvtxY", &reco_nuvtxY);
    T_PFDump->SetBranchAddress("reco_nuvtxZ", &reco_nuvtxZ); 
    T_PFDump->SetBranchAddress("reco_muonMomentum", &reco_muonMomentum);

    T_BDTvars->SetBranchAddress("numu_cc_flag", &numu_cc_flag);  // T_BDTvars
    T_BDTvars->SetBranchAddress("numu_score", &numu_score);
    T_BDTvars->SetBranchAddress("nue_score", &nue_score);
    T_BDTvars->SetBranchAddress("cosmict_flag", &cosmict_flag);

    T_eval->SetBranchAddress("truth_nuEnergy", &truth_nuEnergy);
    T_eval->SetBranchAddress("truth_energyInside", &truth_energyInside);
    T_eval->SetBranchAddress("match_isFC", &match_isFC);
    T_eval->SetBranchAddress("truth_nuPdg", &truth_nuPdg);
    T_eval->SetBranchAddress("truth_isCC", &truth_isCC);
    T_eval->SetBranchAddress("stm_clusterlength", &stm_clusterlength);  // T_eval
    T_eval->SetBranchAddress("truth_isFC", &truth_isFC);                // T_eval
    T_eval->SetBranchAddress("truth_vtxInside", &truth_vtxInside);
    T_eval->SetBranchAddress("match_completeness_energy", &match_completeness_energy);

    T_KINEvars->SetBranchAddress("kine_energy_particle", &kine_energy_particle);
    T_KINEvars->SetBranchAddress("kine_energy_info", &kine_energy_info);
    T_KINEvars->SetBranchAddress("kine_particle_type", &kine_particle_type);
    T_KINEvars->SetBranchAddress("kine_reco_Enu", &kine_reco_Enu); 
    T_KINEvars->SetBranchAddress("kine_reco_add_energy", &kine_reco_add_energy);
    T_KINEvars->SetBranchAddress("kine_pio_flag", &kine_pio_flag);
    T_KINEvars->SetBranchAddress("kine_pio_vtx_dis", &kine_pio_vtx_dis);
    T_KINEvars->SetBranchAddress("kine_pio_energy_1", &kine_pio_energy_1);
    T_KINEvars->SetBranchAddress("kine_pio_energy_2", &kine_pio_energy_2);
    T_KINEvars->SetBranchAddress("kine_pio_dis_1", &kine_pio_dis_1);
    T_KINEvars->SetBranchAddress("kine_pio_dis_2", &kine_pio_dis_2);
    T_KINEvars->SetBranchAddress("kine_pio_angle", &kine_pio_angle);
    T_KINEvars->SetBranchAddress("kine_pio_mass", &kine_pio_mass);

    float energyLim = 2.0;
    auto homega = new TH2F("homega","",100,0,energyLim,100,0,energyLim);
    homega->GetXaxis()->SetTitle("#omega = E_{#nu}^{init.} - E_{#nu}^{final}  (GeV)");
    homega->GetYaxis()->SetTitle("#omega^{reco} (GeV)");
    homega->GetXaxis()->SetTitleFont(22);
    homega->GetYaxis()->SetTitleFont(22);
    homega->GetYaxis()->SetTitleOffset(0.85);

    auto l1 = new TLine(0,0,energyLim,energyLim);
    l1->SetLineStyle(2);
    l1->SetLineColor(2);

    int numu_selection = 0;
    int nue_selection = 0;
    int ncpi0_selection = 0;
    int ncpi0_selection1 = 0;

    for (int ientry = 0; ientry < T_PFDump->GetEntries(); ++ientry) {

        T_PFDump->GetEntry(ientry);
        T_eval->GetEntry(ientry);
        T_BDTvars->GetEntry(ientry);
        T_KINEvars->GetEntry(ientry);

        // if (ientry % 1000 == 0) cout << "processing: " << ientry / 10000. * 100 << "%" << endl;
        if (ientry % 100000 == 0) cout << "processing: " << ientry  << " / " << T_PFDump->GetEntries() << endl;

        if (numu_cc_flag>=0 && numu_score>0.9) numu_selection ++;

        // if (ientry>100000) break;
       
        bool flag_NC = (!cosmict_flag) && numu_score<0;

        bool flag_pi0 = ( (kine_pio_flag==1 && kine_pio_vtx_dis < 9) || kine_pio_flag==2) && kine_pio_energy_1 > 40 && kine_pio_energy_2 > 25 && kine_pio_dis_1 < 110 && kine_pio_dis_2 < 120 && kine_pio_angle > 0 && kine_pio_angle < 174  && kine_pio_mass > 22 && kine_pio_mass < 300 ;

	bool flag_nueCC = numu_cc_flag>=0 && nue_score>7.0;
        if (flag_nueCC) nue_selection ++;

        if (flag_NC && flag_pi0 && (!flag_nueCC) ) { // NCpi0 selection

            ncpi0_selection ++;

            if (!truth_isCC && truth_vtxInside==1 && match_completeness_energy>0.1*truth_energyInside){ // truth NC

              // energy transfer for NC interaction (GeV)
              double omega = 1e-3*truth_nuEnergy - truth_endMomentum[0][3]; 
              homega->Fill(omega, kine_reco_Enu*1e-3);

              ncpi0_selection1 ++;
            }

	    // cout << "entry: " << ientry << " truth_nuPdg: " << truth_nuPdg << " truth_nuEnergy: " << 1e-3*truth_nuEnergy << " omega: " << omega << " kine_reco_Enu: " << kine_reco_Enu*1e-3 << " truth_isCC: " << truth_isCC << endl;

        }

    }

    homega->Draw("colz");
    l1->Draw("same");
    
    cout << "numu CC #: " << numu_selection << endl;
    cout << "nue CC #: " << nue_selection << endl;
    cout << "ncpi0 #: " << ncpi0_selection << endl;
    cout << "plot entries: " << homega->GetEntries() << endl;
    cout << "entries: " << ncpi0_selection1 << endl;
}
