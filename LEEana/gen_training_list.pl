#!/usr/bin/perl

system("./bin/gen_training_list ./old_files/run1_intrinsic_nue_POT1.2E23.root prodgenie_bnb_intrinsic_nue_overlay_run1");
system("./bin/gen_training_list ./old_files/run3_intrinsic_nue_POT8.3E22.root prodgenie_bnb_intrinsic_nue_overlay_run3");
system("./bin/gen_training_list ./old_files/run1_intrinsic_nue_LowE_POT6.1E23.root prodgenie_bnb_intrinsic_nue_overlay_LowE_run1");
system("./bin/gen_training_list ./old_files/run3_intrinsic_nue_LowE_POT6.0E23.root prodgenie_bnb_intrinsic_nue_overlay_LowE_run3");
system("./bin/gen_training_list ./old_files/run1_bnb_nu_POT1.2E21.root prodgenie_bnb_nu_overlay_run1");
system("./bin/gen_training_list ./old_files/run3_bnb_nu_POT1.2E21.root prodgenie_bnb_nu_overlay_run3");
system("./bin/gen_training_list ./old_files/run1_ext_bnb_C1_gt10_wcp_v00_14_00_POT1.2E20.root data_extbnb_run1");
system("./bin/gen_training_list ./old_files/run3_ext_bnb_F_G1_POT1.9E20.root data_extbnb_run3");
system("./bin/gen_training_list ./old_files/run1_bnb_nu_LowE_POT1.6E21.root prodgenie_bnb_nu_overlay_LowE_run1");
system("./bin/gen_training_list ./old_files/run3_bnb_nu_LowE_POT1.5E21.root prodgenie_bnb_nu_overlay_LowE_run3");

