void calc_bin_center(){
  auto file = TFile::Open("checkout_prodgenie_bnb_nu_overlay_run1.root");

  TTree* t1 = (TTree*)file->Get("wcpselection/T_eval");
  TTree* t2 = (TTree*)file->Get("wcpselection/T_PFeval");
  t1->SetBranchStatus("*",0);

  float weight_cv, weight_spline, truth_nuEnergy;
  int truth_nuPdg;
  bool truth_isCC, truth_vtxInside;
  t1->SetBranchStatus("weight_cv",1);
  t1->SetBranchStatus("weight_spline",1);
  t1->SetBranchStatus("truth_nuEnergy",1);
  t1->SetBranchStatus("truth_nuPdg",1);
  t1->SetBranchStatus("truth_isCC",1);
  t1->SetBranchStatus("truth_vtxInside",1);

  t1->SetBranchAddress("weight_cv", &weight_cv);
  t1->SetBranchAddress("weight_spline", &weight_spline);
  t1->SetBranchAddress("truth_nuEnergy", &truth_nuEnergy);
  t1->SetBranchAddress("truth_nuPdg", &truth_nuPdg);
  t1->SetBranchAddress("truth_isCC", &truth_isCC);
  t1->SetBranchAddress("truth_vtxInside", &truth_vtxInside);

  float truth_muonMomentum[4];
  t2->SetBranchStatus("*",0);
  t2->SetBranchStatus("truth_muonMomentum",1);
  t2->SetBranchAddress("truth_muonMomentum", &truth_muonMomentum[0]);

  TH1F* h1 = new TH1F("h1","h1",10,-0.5,9.5);
  TH1F* h2 = new TH1F("h2","h2",10,-0.5,9.5);
  for (int i=0; i<t1->GetEntries(); i++) {
    t1->GetEntry(i);
    t2->GetEntry(i);
//    cout << i << " " << truth_nuPdg << " " << truth_isCC << " " << truth_vtxInside << " " << truth_muonMomentum[3] << endl;

    if (truth_nuPdg==14 && truth_isCC==1 && truth_vtxInside==1) {
      int index = -1;

      // calculate MC weight for each bin
      if (truth_nuEnergy<=540 && truth_nuEnergy>200) index = 0;
      else if (truth_nuEnergy<=705 && truth_nuEnergy>540) index = 1;
      else if (truth_nuEnergy<=805 && truth_nuEnergy>705) index = 2;
      else if (truth_nuEnergy<=920 && truth_nuEnergy>805) index = 3;
      else if (truth_nuEnergy<=1050 && truth_nuEnergy>920) index = 4;
      else if (truth_nuEnergy<=1200 && truth_nuEnergy>1050) index = 5;
      else if (truth_nuEnergy<=1375 && truth_nuEnergy>1200) index = 6;
      else if (truth_nuEnergy<=1570 && truth_nuEnergy>1375) index = 7;
      else if (truth_nuEnergy<=2050 && truth_nuEnergy>1570) index = 8;
      else if (truth_nuEnergy<=4000 && truth_nuEnergy>2050) index = 9;

      if(index!=-1) {
        h1->Fill(index, weight_cv * weight_spline);
        // h2->Fill(index, weight_spline);
        h2->Fill(index);
      }

      // calculate bin center for Emuon
      // if (truth_muonMomentum[3]*1000<=226 && truth_muonMomentum[3]*1000>105.7) index = 0;
      // else if (truth_muonMomentum[3]*1000<=296 && truth_muonMomentum[3]*1000>226) index = 1;
      // else if (truth_muonMomentum[3]*1000<=386 && truth_muonMomentum[3]*1000>296) index = 2;
      // else if (truth_muonMomentum[3]*1000<=505 && truth_muonMomentum[3]*1000>386) index = 3;
      // else if (truth_muonMomentum[3]*1000<=577 && truth_muonMomentum[3]*1000>505) index = 4;
      // else if (truth_muonMomentum[3]*1000<=659 && truth_muonMomentum[3]*1000>577) index = 5;
      // else if (truth_muonMomentum[3]*1000<=753 && truth_muonMomentum[3]*1000>659) index = 6;
      // else if (truth_muonMomentum[3]*1000<=861 && truth_muonMomentum[3]*1000>753) index = 7;
      // else if (truth_muonMomentum[3]*1000<=984 && truth_muonMomentum[3]*1000>861) index = 8;
      // else if (truth_muonMomentum[3]*1000<=1285 && truth_muonMomentum[3]*1000>984) index = 9;
      // else if (truth_muonMomentum[3]*1000<=2506 && truth_muonMomentum[3]*1000>1285) index = 10;

      // if(index!=-1) {
      //   h1->Fill(index, truth_muonMomentum[3] * weight_cv * weight_spline);
      //   h2->Fill(index, weight_cv * weight_spline);
      // }

      //  calculate bin center for Ehad
      // float Ehadron = truth_nuEnergy - truth_muonMomentum[3]*1000;

      // if (Ehadron<=100 && Ehadron>30) index = 0;
      // else if (Ehadron<=150 && Ehadron>100) index = 1;
      // else if (Ehadron<=225 && Ehadron>150) index = 2;
      // else if (Ehadron<=275 && Ehadron>225) index = 3;
      // else if (Ehadron<=336 && Ehadron>275) index = 4;
      // else if (Ehadron<=411 && Ehadron>336) index = 5;
      // else if (Ehadron<=502 && Ehadron>411) index = 6;
      // else if (Ehadron<=614 && Ehadron>502) index = 7;
      // else if (Ehadron<=750 && Ehadron>614) index = 8;
      // else if (Ehadron<=1120 && Ehadron>750) index = 9;
      // else if (Ehadron<=2500 && Ehadron>1120) index = 10;

      // if(index!=-1) {
      //   h1->Fill(index, Ehadron * weight_cv * weight_spline);
      //   h2->Fill(index, weight_cv * weight_spline);
      // }


    } 

  }

  for (int i=0; i<10; i++) {
    cout << h1->GetBinContent(i+1) / h2->GetBinContent(i+1) << ", ";
  }
  cout << endl;

}
