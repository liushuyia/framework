/*
void plot_nueCC_genie_all(){

  // auto file = TFile::Open("/data1/wgu/LEEana-tmp/processed_checkout_rootfiles/prodgenie_bnb_intrinsic_nue_overlay_run3/UBGenieFluxSmallUni.root");
  auto file = TFile::Open("UBGenieFluxSmallUni_intrinsic_nue_overlay_run123.root");
  auto T_eval = (TTree*)file->Get("wcpselection/T_eval");
  auto T_BDTvars = (TTree*)file->Get("wcpselection/T_BDTvars");
  auto T_KINEvars = (TTree*)file->Get("wcpselection/T_KINEvars");
  auto T_weight = (TTree*)file->Get("wcpselection/T_weight");


  cout << T_eval->GetEntries() << endl;
  cout << T_BDTvars->GetEntries() << endl;
  cout << T_KINEvars->GetEntries() << endl;
  cout << T_weight->GetEntries() << endl;


  float weight_cv, weight_spline; // T_eval
  bool match_isFC; // T_eval
  float numu_cc_flag, nue_score; // T_BDTvars, numu_cc_flag>=0, nue_score>7.0
  float kine_reco_Enu; // T_KINEvars
  auto All_UBGenie = new std::vector<float>; // T_weight

  T_eval->SetBranchAddress("weight_cv", &weight_cv);
  T_eval->SetBranchAddress("weight_spline", &weight_spline);
  T_eval->SetBranchAddress("match_isFC", &match_isFC);
  T_BDTvars->SetBranchAddress("numu_cc_flag", &numu_cc_flag);
  T_BDTvars->SetBranchAddress("nue_score", &nue_score);
  T_KINEvars->SetBranchAddress("kine_reco_Enu", &kine_reco_Enu);
  T_weight->SetBranchAddress("All_UBGenie", &All_UBGenie);

  auto h1 = new TH1F("h1","#nu_{e}CC FC",25,0,25);
  auto h2 = new TH1F("h2","#nu_{e}CC PC",25,0,25);
  auto h3 = new TH1F("h3","#nu_{e}CC FC",25,0,25);
  auto h4 = new TH1F("h4","#nu_{e}CC PC",25,0,25);
  h1->SetLineColor(1);
  h2->SetLineColor(1);
  const int Nk = 600; // size of a knob
  // TH1F* hgenie1[Nk]; // ensemble for FC
  // TH1F* hgenie2[Nk]; // ensemble for PC
  vector<TH1F*> hgenie1(Nk);
  vector<TH1F*> hgenie2(Nk);
  for (int i=0; i<Nk; i++) {
    hgenie1[i] = new TH1F(Form("h1_%d",i),"#nu_{e}CC FC",25,0,25);
    hgenie2[i] = new TH1F(Form("h2_%d",i),"#nu_{e}CC PC",25,0,25);
    hgenie1[i]->SetLineColor(4);
    hgenie2[i]->SetLineColor(4);
  }
  auto hgenie_rel = new TH1F("hgenie_rel","All_UBGenie",100,0,2); // relative weight

  for (int i=0; i<T_eval->GetEntries(); i++)  {
  // for (int i=0; i<10000; i++)  {
    if (i%10000==0) cout << "entry: " << i << endl;
    T_eval->GetEntry(i);
    T_BDTvars->GetEntry(i);
    T_KINEvars->GetEntry(i);
    T_weight->GetEntry(i);

    // cout << weight_cv << " " << weight_spline << " " << match_isFC << " " << numu_cc_flag << " " << nue_score << " " << kine_reco_Enu << " " << endl; 
    // cout << "knob size: " << All_UBGenie->size() << endl;
    
    // select nueCC -> select FC/PC -> fill with weight
    if (numu_cc_flag>=0 and nue_score>7.0) {

      if (std::isnan(weight_spline) or std::isinf(weight_spline)) 
        weight_spline = 1.0;

      double event_weight = weight_cv * weight_spline;
      if (std::isnan(event_weight) or std::isinf(event_weight)) 
        event_weight = 1.0;

      float Enu = kine_reco_Enu * 1e-2; // x100MeV

      for(int j=0; j<Nk; j++) {

          double genie_rel_weight = All_UBGenie->at(j);
          if (std::isnan(genie_rel_weight) or std::isinf(genie_rel_weight)) 
            genie_rel_weight = 1.0;

          hgenie_rel->Fill(genie_rel_weight);
          if (match_isFC)
            hgenie1[j]->Fill(Enu, weight_spline * genie_rel_weight);
          else
            hgenie2[j]->Fill(Enu, weight_spline * genie_rel_weight);
      }

      // double genie_rel_weight = All_UBGenie->at(0);
      // if (std::isnan(genie_rel_weight) or std::isinf(genie_rel_weight)) 
      //   genie_rel_weight = 1.0;
      // if(match_isFC) h3->Fill(kine_reco_Enu, event_weight * genie_rel_weight);
      // else h4->Fill(kine_reco_Enu, event_weight * genie_rel_weight);

      if (match_isFC) {
        h1->Fill(Enu, event_weight);
      }
      else{
        h2->Fill(Enu, event_weight);
      }
    }
    
  } 

  // hgenie->Draw();

  cout << "FC integral: " << h1->Integral() << endl;
  cout << "PC integral: " << h2->Integral() << endl;

  auto ofile = new TFile("nueCC_All_UBGenie.root","recreate");
  h1->Write();
  h2->Write();
  // h3->Write();
  // h4->Write();
  for(int j=0; j<Nk; j++) {
    hgenie1[j]->Write();
    hgenie2[j]->Write();
  }
  hgenie_rel->Write();
  ofile->Close();

  
  auto c1 = new TCanvas("c1","",1200,600); 
  c1->Divide(2,1);
  h1->GetYaxis()->SetRangeUser(0,20000);
  h2->GetYaxis()->SetRangeUser(0,20000);
  
  c1->cd(1);
  h1->Draw("hist");
  for(int j=0; j<Nk; j++) {
    hgenie1[j]->Draw("same hist");
  }
  h1->Draw("same hist");
  
  c1->cd(2);
  h2->Draw("hist");
  for(int j=0; j<Nk; j++) {
    hgenie2[j]->Draw("same hist");
  }
  h2->Draw("same hist");

  c1->SaveAs("canvas.root");
}
*/


void plot_nueCC_genie_all(){
  auto c1 = new TCanvas("c1","",1200,600); 
  c1->Divide(2,1);

  auto file1 = TFile::Open("file_signalonly_nueCC_FC.root");
  auto h1_wo = (TH1D*)file1->Get("h1_pred_Y_noConstraint_clone");
  auto h1_wi = (TH1D*)file1->Get("h1_pred_Y_wiConstraint_clone");

  auto file2 = TFile::Open("file_signalonly_nueCC_PC.root");
  auto h2_wo = (TH1D*)file2->Get("h1_pred_Y_noConstraint_clone");
  auto h2_wi = (TH1D*)file2->Get("h1_pred_Y_wiConstraint_clone");
  
  auto file = TFile::Open("nueCC_All_UBGenie.root");
  auto hpredcv1 = (TH1F*)file->Get("h1"); // CV FC
  auto hpredcv2 = (TH1F*)file->Get("h2"); // CV PC
  auto hgenie = (TH1F*)file->Get("hgenie_rel"); // GENIE weights

  double norm1 = h1_wo->Integral(1,26) / hpredcv1->Integral(1,26);
  double norm2 = h2_wo->Integral(1,26) / hpredcv2->Integral(1,26);
  cout << "normalization FC: " << norm1 << endl;
  cout << "normalization PC: " << norm1 << endl;

  hpredcv1->Scale(norm1);
  hpredcv2->Scale(norm1); // use the same normalization

  c1->cd(1);
  hpredcv1->GetYaxis()->SetRangeUser(0,100);
  hpredcv1->Draw("hist");
  c1->cd(2);
  hpredcv2->GetYaxis()->SetRangeUser(0,100);
  hpredcv2->Draw("hist");

  for(int i=0; i<600; i++) {
    auto hpred1 = (TH1F*)file->Get(Form("h1_%d", i));
    auto hpred2 = (TH1F*)file->Get(Form("h2_%d", i));
    hpred1->Scale(norm1);
    hpred2->Scale(norm1); // use the same normalization

    c1->cd(1);
    hpred1->Draw("same hist");
    c1->cd(2);
    hpred2->Draw("same hist");
  }

  c1->cd(1);
  hpredcv1->Draw("same hist");
  h1_wi->SetLineStyle(2);
  h1_wi->SetLineColor(2);
  h1_wi->Draw("same hist");
  h1_wo->SetLineColor(2);
  h1_wo->Draw("same hist");

  auto lg = new TLegend(0.4,0.6,0.8,0.8);
  hgenie->SetFillColor(4);
  lg->AddEntry(hgenie,"GENIE All","f")->SetTextColor(4);
  lg->AddEntry(h1_wo,"Pred wo constraint","l")->SetTextColor(2);
  lg->AddEntry(h1_wi,"Pred wi constraint","l")->SetTextColor(2);
  lg->SetBorderSize(0);
  lg->Draw();

  c1->cd(2);
  hpredcv2->Draw("same hist");
  h2_wi->SetLineStyle(2);
  h2_wi->SetLineColor(2);
  h2_wi->Draw("same hist");
  h2_wo->SetLineColor(2);
  h2_wo->Draw("same hist");


}
