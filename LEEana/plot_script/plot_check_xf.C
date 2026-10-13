void plot_check_xf(){
  TFile **file = new TFile*[17];
  for (Int_t i=0;i!=17;i++){
    file[i] = new TFile(Form("hist_rootfiles/XsFlux/cov_%d.root",i+1));
  }

  TH1F **h1 = new TH1F*[17];
  TH1F **h2 = new TH1F*[17];
  for (Int_t i=0;i!=17;i++){
    h1[i] = (TH1F*)file[i]->Get("pred_covch_4");
    h2[i] = (TH1F*)file[i]->Get("pred_covch_1");
  }
  
  TCanvas *c1 = new TCanvas("c1","c1",1200,800);
  c1->Divide(3,2);

   c1->cd(1);
  h1[0]->Draw("");  h1[0]->SetLineColor(1);
    
  for (Int_t i=1;i!=9;i++){
    h1[i]->SetLineColor(i+1);
    h1[i]->Draw("histsame");
  }

  c1->cd(2);
  h2[0]->Draw("");  h2[0]->SetLineColor(1);
    
  for (Int_t i=1;i!=9;i++){
    h2[i]->SetLineColor(i+1);
    h2[i]->Draw("histsame");
  }
  
  c1->cd(3);
  TGraph **g1 = new TGraph*[17];
  TMatrixD** frac_mat = new TMatrixD*[17];

  frac_mat[0] = (TMatrixD*)file[0]->Get(Form("frac_cov_xf_mat_%d",1));
  const int nbin = frac_mat[0]->GetNcols();
  double rel_err2[nbin];
  for (int i=0;i!=nbin;i++){
    rel_err2[i] = 0;
  }
  
  for(Int_t i=0;i!=17;i++){
    frac_mat[i] = (TMatrixD*)file[i]->Get(Form("frac_cov_xf_mat_%d",i+1));
    g1[i] = new TGraph();
    for (Int_t j=0;j!=frac_mat[i]->GetNcols();j++){
      g1[i]->SetPoint(j,j,sqrt((*frac_mat[i])(j,j)));
      rel_err2[j] += (*frac_mat[i])(j,j);
    }
  }
 

  g1[0]->Draw("A*");
  g1[0]->SetLineColor(1);
  g1[0]->SetMarkerColor(1);
  g1[0]->SetMarkerStyle(20);
   for (Int_t i=1;i!=17;i++){
     g1[i]->Draw("*same");
     g1[i]->SetLineColor(i+1);
     g1[i]->SetMarkerColor(i+1);
     g1[i]->SetMarkerStyle(20);
  }

   
   c1->cd(4);
   TLegend *le1 = new TLegend(0.1,0.1,0.89,0.89);
   le1->AddEntry(g1[0],"expskin_FluxUnisim","pl");
   le1->AddEntry(g1[1],"horncurrent_FluxUnisim","pl");
   le1->AddEntry(g1[2],"kminus_PrimaryHadronNormalization","pl");
   le1->AddEntry(g1[3],"kplus_PrimaryHadronFeynmanScaling","pl");
   le1->AddEntry(g1[4],"kzero_PrimaryHadronSanfordWang","pl");
   le1->AddEntry(g1[5],"nucleoninexsec_FluxUnisim","pl");
   le1->AddEntry(g1[6],"nucleonqexsec_FluxUnisim","pl");
   le1->AddEntry(g1[7],"nucleontotxsec_FluxUnisim","pl");
   le1->AddEntry(g1[8],"piminus_PrimaryHadronSWCentralSplineVariation","pl");
   le1->AddEntry(g1[9],"pioninexsec_FluxUnisim","pl");
   le1->AddEntry(g1[10],"pionqexsec_FluxUnisim","pl");
   le1->AddEntry(g1[11],"piontotxsec_FluxUnisim","pl");
   le1->AddEntry(g1[12],"piplus_PrimaryHadronSWCentralSplineVariation","pl");
   le1->AddEntry(g1[13],"reinteractions_piminus_Geant4","pl");
   le1->AddEntry(g1[14],"reinteractions_piplus_Geant4","pl");
   le1->AddEntry(g1[15],"reinteractions_proton_Geant4","pl");
   le1->AddEntry(g1[16],"All_UBGenie","pl");
   le1->Draw();

   c1->cd(5);
   TGraph **g2 = new TGraph*[9];
   for (Int_t i=0;i!=9;i++){
     g2[i] = new TGraph();
   }
   
   for (int i=0;i!=26;i++){     
     g2[0]->SetPoint(i, i+0.5, sqrt(rel_err2[i]));
     g2[1]->SetPoint(i, i+0.5, sqrt(rel_err2[i+26]));
     g2[2]->SetPoint(i, i+0.5, sqrt(rel_err2[i+26*2]));
     g2[3]->SetPoint(i, i+0.5, sqrt(rel_err2[i+26*3]));
     g2[7]->SetPoint(i, i+0.5, sqrt(rel_err2[i+26*4+21*3]));
     g2[8]->SetPoint(i, i+0.5, sqrt(rel_err2[i+26*4+21*3+26]));
   }

   for (Int_t i=0;i!=21;i++){
     g2[4]->SetPoint(i,i+0.5,sqrt(rel_err2[i+26*4]));
     g2[5]->SetPoint(i,i+0.5,sqrt(rel_err2[i+26*4+21]));
     g2[6]->SetPoint(i,i+0.5,sqrt(rel_err2[i+26*4+21*2]));
   }

   c1->cd(3);
   g2[2]->Draw("A*");
   g2[2]->SetMarkerStyle(21);
   g2[2]->SetMarkerColor(1);
   g2[3]->Draw("*same");
   g2[3]->SetMarkerStyle(21);
   g2[3]->SetMarkerColor(2);

   TLegend *le2 = new TLegend(0.6,0.6,0.89,0.89);
   le2->AddEntry(g2[2],"numuCC FC","p");
   le2->AddEntry(g2[3],"numuCC PC","p");
   le2->Draw();
   
   
   c1->cd(5);
   g2[0]->Draw("A*");
   g2[0]->SetMarkerStyle(20);
   g2[1]->Draw("*same");
   g2[1]->SetMarkerStyle(20);
   g2[1]->SetMarkerColor(2);
   g2[7]->Draw("*same");
   g2[7]->SetMarkerStyle(23);
   g2[7]->SetMarkerColor(4);
   g2[8]->Draw("*same");
   g2[8]->SetMarkerStyle(23);
   g2[8]->SetMarkerColor(6);

   TLegend *le3 = new TLegend(0.6,0.6,0.89,0.89);
   le3->AddEntry(g2[0],"int. nueCC FC","p");
   le3->AddEntry(g2[1],"int. nueCC PC","p");
   le3->AddEntry(g2[7],"LEE FC","p");
   le3->AddEntry(g2[8],"LEE PC","p");
   le3->Draw();
   
   c1->cd(6);
   g2[4]->Draw("A*");
   g2[4]->SetMarkerStyle(22);
   g2[5]->Draw("*same");
   g2[5]->SetMarkerStyle(22);
   g2[5]->SetMarkerColor(2);
   g2[6]->Draw("*same");
   g2[6]->SetMarkerStyle(22);
   g2[6]->SetMarkerColor(4);

   TLegend *le4 = new TLegend(0.6,0.6,0.89,0.89);
   le4->AddEntry(g2[4],"CC pio FC","p");
   le4->AddEntry(g2[5],"CC pio PC","p");
   le4->AddEntry(g2[6],"NC pio","p");
   le4->Draw();
   
}
