void test_xs(){
  TFile *file = new TFile("processed_checkout_rootfiles/checkout_prodgenie_bnb_nu_overlay_run1.root");
  TTree *t1 = (TTree*)file->Get("wcpselection/T_eval");
  TH1F *h1 = new TH1F("h1","h1",100,-10,10);
  TH1F *h2 = new TH1F("h2","h2",100,-10,10);

  TH1F *h3 = new TH1F("h3","h3",100,-10,10);
  TH1F *h4 = new TH1F("h4","h4",100,-10,10);

  t1->Project("h1","truth_isCC","weight_cv*(match_completeness_energy/truth_energyInside>=0.1 && truth_nuPdg==14 && truth_isCC==1 && truth_vtxInside==1 && truth_nuEnergy<=2500 && truth_nuEnergy>2000)");
  t1->Project("h2","truth_isCC","weight_cv*weight_cv*(match_completeness_energy/truth_energyInside>=0.1 && truth_nuPdg==14 && truth_isCC==1 && truth_vtxInside==1 && truth_nuEnergy<=2500 && truth_nuEnergy>2000)");
  std::cout << h1->GetSum()/498.661 << " +- " << sqrt(h2->GetSum())/498.661 << std::endl;

  t1->Project("h3","truth_isCC","weight_cv*(match_completeness_energy/truth_energyInside>=0.1 && truth_nuPdg==14 && truth_isCC==1 && truth_vtxInside==1 && truth_nuEnergy<=2000 && truth_nuEnergy>1600)");
  t1->Project("h4","truth_isCC","weight_cv*weight_cv*(match_completeness_energy/truth_energyInside>=0.1 && truth_nuPdg==14 && truth_isCC==1 && truth_vtxInside==1 && truth_nuEnergy<=2000 && truth_nuEnergy>1600)");
  std::cout << h3->GetSum()/1676.83 << " +- " << sqrt(h4->GetSum())/1676.83 << std::endl;
  // h1->Draw();
}
