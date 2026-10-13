#!/usr/bin/perl


#central histogram preparation ...
system("./convert_histo.pl");

#Flux sys, GEANT4 1-->16
system("./run_xf_sys.pl");

# MC stat ... 
system("./bin/merge_hist -r0 -l0 -e2 > ./mc_stat/xs_tot.log");

#central files ...
system("./bin/merge_hist -r0 -l1");
