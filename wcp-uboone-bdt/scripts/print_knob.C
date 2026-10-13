void print_knob(){
  auto file = TFile::Open("/pnfs/uboone/persistent/users/jjo/WC_NuMI/files/data_mc_cv/checkout_prodgenie_numi_overlay_run1.root");
  auto t = (TTree*)file->Get("wcpweights/T_wgt");
  // t->Show(0);
  // t->Print();

  auto mcweight = new map<string,vector<float> >;
  t->SetBranchAddress("mcweight", &mcweight);

  int n=0;
  for(int i=0; i<t->GetEntries(); i++){
    t->GetEntry(i);
    if (mcweight->size()>0 && n==0) {
      cout << mcweight->size() << endl;
      n ++;

      for(auto const& e: *mcweight){
        cout << e.first << " " << e.second.size() << endl;
      }
    }
  }


  

}
