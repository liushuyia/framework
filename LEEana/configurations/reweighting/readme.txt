Usage:


Terse:
In Configure_Lee.h, set the reweight flags to true:
bool flag_syst_reweight        = 1;
bool flag_syst_reweight_cor    = 1;

Make sure you have a valid xf_cv_input.txt and rw_sys_input.txt and you have the required _PF files if needed for your reweighting. Rerun your ./convert_histo.pl and systematics, reweighting will be automatically applied. 
Then, to get the reweighting uncertainty, run:
bin/xf_cov_matrix -r18
bin/xf_cov_matrix -r19
If you are doing xs, instead of the running bin/xs_cov_matrix -r17, you run
bin/xs_cov_matrix -r0

Other than this, everything else proceeds as usual. Note that one need not apply the reweighting to the cv when calculating the rw uncertainty and vice versa; the code supports all options. Also not that all current reweighting functions require _PF files.


Detailed:
https://docs.google.com/document/d/19SUsQKqudMxFZ1j5ZcSxDXE5YRlT6d2suFc7GglWc98/edit?usp=sharing
