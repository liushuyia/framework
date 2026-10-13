#!/usr/bin/perl

system("cp ../mc_stat/xs_tot.log ./mcstat/");
system("cp ../hist_rootfiles/DetVar/cov*.root ./DetVar/");
system("cp ../hist_rootfiles/XsFlux/cov_*.root ./XsFlux/");
system("cp ../merge_xs.root .");

