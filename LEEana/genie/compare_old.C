void compare(){
  TFile *file1 = new TFile("xsec_graphs_v2_tune1.root");
  TFile *file2 = new TFile("xsec_graphs_v3_10a_02_11a.root");

  TGraph *g1 = (TGraph*)file1->Get("nu_mu_Ar40/tot_cc");
  TGraph *g2 = (TGraph*)file2->Get("nu_mu_Ar40/tot_cc");
  TGraph *g3 = (TGraph*)file1->Get("nu_mu_Ar40/qel_cc_n");
  TGraph *g4 = (TGraph*)file2->Get("nu_mu_Ar40/qel_cc_n");

  TGraph *g5 = (TGraph*)file1->Get("nu_mu_Ar40/res_cc_n");
  TGraph *g6 = (TGraph*)file2->Get("nu_mu_Ar40/res_cc_n");

  TGraph *g7 = (TGraph*)file1->Get("nu_mu_Ar40/res_cc_p");
  TGraph *g8 = (TGraph*)file2->Get("nu_mu_Ar40/res_cc_p");

  TGraph *g10 = (TGraph*)file1->Get("nu_mu_Ar40/mec_cc");
  TGraph *g11 = (TGraph*)file2->Get("nu_mu_Ar40/mec_cc");

  TGraph *g12 = (TGraph*)file1->Get("nu_mu_Ar40/dis_cc_p");
  TGraph *g13 = (TGraph*)file2->Get("nu_mu_Ar40/dis_cc_p");
  TGraph *g14 = (TGraph*)file1->Get("nu_mu_Ar40/dis_cc_n");
  TGraph *g15 = (TGraph*)file2->Get("nu_mu_Ar40/dis_cc_n");
  
  g1->Draw("AL");
  g1->GetXaxis()->SetRangeUser(0,3);
  g1->GetYaxis()->SetRangeUser(0,120);
  g2->Draw("Lsame");
  g1->SetLineColor(1);
  g2->SetLineColor(1);
  g1->SetLineWidth(2);
  g2->SetLineWidth(2);
  g1->SetLineStyle(2);

  // QE
  g3->Draw("Lsame");
  g3->SetLineStyle(2);
  g4->Draw("Lsame");
  g3->SetLineColor(2);
  g4->SetLineColor(2);
  g3->SetLineWidth(2);
  g4->SetLineWidth(2);

  // RES
  g5->Draw("Lsame");
  g5->SetLineStyle(2);
  g6->Draw("Lsame");
  g5->SetLineColor(4);
  g6->SetLineColor(4);
  g5->SetLineWidth(2);
  g6->SetLineWidth(2);
  
  g7->Draw("Lsame");
  g7->SetLineStyle(2);
  g8->Draw("Lsame");
  g7->SetLineColor(3);
  g8->SetLineColor(3);
  g7->SetLineWidth(2);
  g8->SetLineWidth(2);
  
  
  // MEC
  g10->Draw("Lsame");
  g10->SetLineStyle(2);
  g11->Draw("Lsame");
  g10->SetLineColor(6);
  g11->SetLineColor(6);
  g10->SetLineWidth(2);
  g11->SetLineWidth(2);

 

  //DIS
  g12->Draw("Lsame");
  g12->SetLineStyle(2);
  g13->Draw("Lsame");
  g12->SetLineColor(8);
  g13->SetLineColor(8);
  g12->SetLineWidth(2);
  g13->SetLineWidth(2);

  g14->Draw("Lsame");
  g14->SetLineStyle(2);
  g15->Draw("Lsame");
  g14->SetLineColor(7);
  g15->SetLineColor(7);
  g14->SetLineWidth(2);
  g15->SetLineWidth(2);

  g1->GetXaxis()->SetTitle("E_{#nu} (GeV)");
  g1->GetYaxis()->SetTitle("10^{-38} cm^{2}/Ar");
  TLegend *le1 = new TLegend(0.6,0.6,0.89,0.89);
  le1->SetHeader("v3 10a_02_11a");
  le1->AddEntry(g2,"Total","l");
  le1->AddEntry(g4,"QE","l");
  le1->AddEntry(g6,"RES n","l");
  le1->AddEntry(g8,"RES p","l");
  le1->AddEntry(g11,"MEC","l");
  le1->AddEntry(g13,"DIS p","l");
  le1->AddEntry(g15,"DIS n","l");
  le1->Draw();
  
  
}
