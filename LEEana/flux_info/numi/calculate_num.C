void calculate_num(double init_E = 0., double final_E = 5, int nstep = 1000, double POT = 2.099e+20, double init_difval=0, double final_difval=0){


  // convert histogram to graph
  auto file1 = TFile::Open("output_uboone_fhc_run0_set1.root");
  TH1D* h1 = (TH1D*)file1->Get("numu/Detsmear/numu_CV_AV_TPC_5MeV_bin"); // nu / x POT / 5 MeV / m2
                                                                         // POT 5e5*492
                                                                         //
                                                                         // nu / POT / GeV / cm2

  TGraph* gflux_numu = new TGraph();
  for (int i=0; i<h1->GetNbinsX(); i++) {
    gflux_numu->SetPoint(i, h1->GetBinCenter(i+1), h1->GetBinContent(i+1)/(5e5*492.0)*200.0/1e4 );
  }


  // TFile *file = new TFile("numi_flux_graphs.root"); 
  // TGraph *gflux_numu = (TGraph*)file->Get("numu"); //  #nu / POT/ GeV / cm2
  // TGraph *gflux_numubar = (TGraph*)file->Get("numubar"); 

  double density = 1.3836; //g/cm^3
  double volume = 5.82515e7; //cm^3
  double m_mol = 39.95; //g/mol
  double NA = 6.022e23; // N/mol

  double target_N = density * volume * NA/m_mol;
  std::cout << "Number of Target Argon: " << target_N << std::endl;

  double pb = 1e-36; // pbarn
  std::cout << "POT                   : " << POT << std::endl;
  
  double flux_sum = 0;
  for (Int_t i=0;i!=nstep;i++){
    double x = init_E + (final_E-init_E)*(i+0.5)/nstep;
    flux_sum += gflux_numu->Eval(x)*(final_E-init_E)/nstep;
    // flux_sum=gFlux->Integral(init_E/5.0,final_E/5.0)/10000/5e5/492;
  }
  std::cout << "Flux                  : " << flux_sum << std::endl;

  std::cout << std::endl;
  double final_constant = flux_sum * POT * target_N * pb;

  if (final_difval!=init_difval and final_difval>0)
  std::cout << "Final constant        : " << final_constant * (final_difval - init_difval) << std::endl;
  else 
  std::cout << "Final constant        : " << final_constant << std::endl;
}
