if (use_rootino_reweight)
  {
    std::cout
        << "\n====================================================\n"
        << "[ROOTINO] Loading Run4b + Run4c + Run4d + Run5 "
        << "reference samples\n"
        << "====================================================\n";


    // --------------------------------------------------------
    // Open reference files
    // --------------------------------------------------------

    TFile *rootino_file1 = TFile::Open(ROOTINO_RUN1_FILE, "READ");
    TFile *rootino_file2 = TFile::Open(ROOTINO_RUN2_FILE, "READ");
    TFile *rootino_file3 = TFile::Open(ROOTINO_RUN3_FILE, "READ");
    TFile *rootino_file4b = TFile::Open(ROOTINO_RUN4B_FILE, "READ");
    TFile *rootino_file4c = TFile::Open(ROOTINO_RUN4C_FILE, "READ");
    TFile *rootino_file4d = TFile::Open(ROOTINO_RUN4D_FILE, "READ");
    TFile *rootino_file5 = TFile::Open(ROOTINO_RUN5_FILE, "READ");

    if (!rootino_file1 || rootino_file1->IsZombie())
    {
        std::cerr
            << "[ROOTINO] ERROR opening Run1:\n"
            << ROOTINO_RUN1_FILE
            << std::endl;
        return 1;
    }

    if (!rootino_file2 || rootino_file2->IsZombie())
    {
        std::cerr
            << "[ROOTINO] ERROR opening Run2:\n"
            << ROOTINO_RUN2_FILE
            << std::endl;
        return 1;
    }

    if (!rootino_file3 || rootino_file3->IsZombie())
    {
        std::cerr
            << "[ROOTINO] ERROR opening Run3:\n"
            << ROOTINO_RUN3_FILE
            << std::endl;
        return 1;
    }

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

    // --------------------------------------------------------
    // Get trees
    // --------------------------------------------------------
    TTree *R1_eval = (TTree*)rootino_file1->Get("wcpselection/T_eval");
    TTree *R1_pot = (TTree*)rootino_file1->Get("wcpselection/T_pot");
    TTree *R1_BDT = (TTree*)rootino_file1->Get("wcpselection/T_BDTvars");
    TTree *R1_PFeval = (TTree*)rootino_file1->Get("wcpselection/T_PFeval");
    TTree *R1_KINE = (TTree*)rootino_file1->Get("wcpselection/T_KINEvars");

    TTree *R2_eval = (TTree*)rootino_file2->Get("wcpselection/T_eval");
    TTree *R2_pot = (TTree*)rootino_file2->Get("wcpselection/T_pot");
    TTree *R2_BDT = (TTree*)rootino_file2->Get("wcpselection/T_BDTvars");
    TTree *R2_PFeval = (TTree*)rootino_file2->Get("wcpselection/T_PFeval");
    TTree *R2_KINE = (TTree*)rootino_file2->Get("wcpselection/T_KINEvars");

    TTree *R3_eval = (TTree*)rootino_file3->Get("wcpselection/T_eval");
    TTree *R3_pot = (TTree*)rootino_file3->Get("wcpselection/T_pot");
    TTree *R3_BDT = (TTree*)rootino_file3->Get("wcpselection/T_BDTvars");
    TTree *R3_PFeval = (TTree*)rootino_file3->Get("wcpselection/T_PFeval");
    TTree *R3_KINE = (TTree*)rootino_file3->Get("wcpselection/T_KINEvars");

    TTree *R4b_eval = (TTree*)rootino_file4b->Get("wcpselection/T_eval");
    TTree *R4b_pot = (TTree*)rootino_file4b->Get("wcpselection/T_pot");
    TTree *R4b_BDT = (TTree*)rootino_file4b->Get("wcpselection/T_BDTvars");
    TTree *R4b_PFeval = (TTree*)rootino_file4b->Get("wcpselection/T_PFeval");
    TTree *R4b_KINE = (TTree*)rootino_file4b->Get("wcpselection/T_KINEvars");

    TTree *R4c_eval = (TTree*)rootino_file4c->Get("wcpselection/T_eval");
    TTree *R4c_pot = (TTree*)rootino_file4c->Get("wcpselection/T_pot");
    TTree *R4c_BDT = (TTree*)rootino_file4c->Get("wcpselection/T_BDTvars");
    TTree *R4c_PFeval = (TTree*)rootino_file4c->Get("wcpselection/T_PFeval");
    TTree *R4c_KINE = (TTree*)rootino_file4c->Get("wcpselection/T_KINEvars");

    TTree *R4d_eval = (TTree*)rootino_file4d->Get("wcpselection/T_eval");
    TTree *R4d_pot = (TTree*)rootino_file4d->Get("wcpselection/T_pot");
    TTree *R4d_BDT = (TTree*)rootino_file4d->Get("wcpselection/T_BDTvars");
    TTree *R4d_PFeval = (TTree*)rootino_file4d->Get("wcpselection/T_PFeval");
    TTree *R4d_KINE =(TTree*)rootino_file4d->Get("wcpselection/T_KINEvars");

    TTree *R5_eval = (TTree*)rootino_file5->Get("wcpselection/T_eval");
    TTree *R5_pot = (TTree*)rootino_file5->Get("wcpselection/T_pot");
    TTree *R5_BDT = (TTree*)rootino_file5->Get("wcpselection/T_BDTvars");
    TTree *R5_PFeval = (TTree*)rootino_file5->Get("wcpselection/T_PFeval");
    TTree *R5_KINE = (TTree*)rootino_file5->Get("wcpselection/T_KINEvars");

    if (!R1_eval || !R1_pot || !R2_eval || !R2_pot || !R3_eval || !R3_pot || !R4b_eval || !R4b_pot || !R4c_eval || !R4c_pot || !R4d_eval || !R4d_pot || !R5_eval  || !R5_pot)
    {
        std::cerr
            << "[ROOTINO] ERROR: missing T_eval/T_pot "
            << "in one or more reference files."
            << std::endl;

        return 1;
    }

    // --------------------------------------------------------
    // Calculate POT of Run1
    // --------------------------------------------------------
    POTInfo R1_pot_info;
    R1_pot->SetBranchAddress("pot_tor875", &R1_pot_info.pot_tor875);
    double POT1 = 0.0;
    for (Long64_t j = 0; j < R1_pot->GetEntries(); ++j)
    {
        R1_pot->GetEntry(j);
        POT1 += R1_pot_info.pot_tor875;
    }

    // --------------------------------------------------------
    // Calculate POT of Run2
    // --------------------------------------------------------
    POTInfo R2_pot_info;
    R2_pot->SetBranchAddress("pot_tor875", &R2_pot_info.pot_tor875);
    double POT2 = 0.0;
    for (Long64_t j = 0; j < R2_pot->GetEntries(); ++j)
    {
        R2_pot->GetEntry(j);
        POT2 += R2_pot_info.pot_tor875;
    }

    // --------------------------------------------------------
    // Calculate POT of Run3
    // --------------------------------------------------------
    POTInfo R3_pot_info;
    R3_pot->SetBranchAddress("pot_tor875", &R3_pot_info.pot_tor875);
    double POT3 = 0.0;
    for (Long64_t j = 0; j < R3_pot->GetEntries(); ++j)
    {
        R3_pot->GetEntry(j);
        POT3 += R3_pot_info.pot_tor875;
    }

    // --------------------------------------------------------
    // Calculate POT of Run4b
    // --------------------------------------------------------
    POTInfo R4b_pot_info;
    R4b_pot->SetBranchAddress("pot_tor875", &R4b_pot_info.pot_tor875);
    double POT4b = 0.0;
    for (Long64_t j = 0; j < R4b_pot->GetEntries(); ++j)
    {
        R4b_pot->GetEntry(j);
        POT4b += R4b_pot_info.pot_tor875;
    }

    // --------------------------------------------------------
    // Run4c
    // --------------------------------------------------------

    POTInfo R4c_pot_info;
    R4c_pot->SetBranchAddress("pot_tor875", &R4c_pot_info.pot_tor875);
    double POT4c = 0.0;
    for (Long64_t j = 0; j < R4c_pot->GetEntries(); ++j)
    {
        R4c_pot->GetEntry(j);
        POT4c += R4c_pot_info.pot_tor875;
    }

    // --------------------------------------------------------
    // Run4d
    // --------------------------------------------------------

    POTInfo R4d_pot_info;
    R4d_pot->SetBranchAddress("pot_tor875", &R4d_pot_info.pot_tor875);
    double POT4d = 0.0;
    for (Long64_t j = 0; j < R4d_pot->GetEntries(); ++j)
    {
        R4d_pot->GetEntry(j);
        POT4d += R4d_pot_info.pot_tor875;
    }

    // --------------------------------------------------------
    // Run5
    // --------------------------------------------------------

    POTInfo R5_pot_info;
    R5_pot->SetBranchAddress("pot_tor875", &R5_pot_info.pot_tor875);
    double POT5 = 0.0;
    for (Long64_t j = 0; j < R5_pot->GetEntries(); ++j)
    {
        R5_pot->GetEntry(j);
        POT5 += R5_pot_info.pot_tor875;
    }

    // --------------------------------------------------------
    // Total Run4+Run5 POT
    // --------------------------------------------------------
    double POT45 =POT4b + POT4c + POT4d + POT5;

    if (POT45 <= 0.0)
    {
        std::cerr
            << "[ROOTINO] ERROR: POT45 <= 0."
            << std::endl;

        return 1;
    }

    double POTtarget = 0.0;
    if (is_run1)
    {
        POTtarget = POT1;
    }
    if (is_run2)
    {
        POTtarget = POT2;
    }
    if (is_run3)
    {
        POTtarget = POT3;
    }
    if (is_run123)
    {
        POTtarget = POT1 + POT2 + POT3;
    }
    if (ROOTINO_ALL_RUNS)
    {
        POTtarget = POT1 + POT2 + POT3 + POT45;
    }


    // --------------------------------------------------------
    // ROOTino replacement scale
    // --------------------------------------------------------

    double ROOTino_reference_scale = POTtarget / POT45;


    std::cout
        << "\n[ROOTINO] POT1 = "
        << POT1
        << "\n[ROOTINO] POT2 = "
        << POT2
        << "\n[ROOTINO] POT3 = "
        << POT3
        << "\n[ROOTINO] POT4b = "
        << POT4b
        << "\n[ROOTINO] POT4c = "
        << POT4c
        << "\n[ROOTINO] POT4d = "
        << POT4d
        << "\n[ROOTINO] POT5  = "
        << POT5
        << "\n[ROOTINO] POT45 = "
        << POT45
        << "\n[ROOTINO] Target POT = "
        << POTtarget
        << "\n[ROOTINO] Replacement scale = "
        << ROOTino_reference_scale
        << std::endl;
  }