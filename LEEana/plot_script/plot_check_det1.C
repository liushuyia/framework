void plot_check_det1(){
  TFile **file = new TFile*[10];
  
  file[0] = new TFile("hist_rootfiles/DetVar/cov_LYDown.root"); //1
  file[5] = new TFile("hist_rootfiles/DetVar/cov_WMThetaXZ.root");  //6
  file[1] = new TFile("hist_rootfiles/DetVar/cov_LYRayleigh.root"); //2
  file[6] = new TFile("hist_rootfiles/DetVar/cov_WMThetaYZ.root");  //7
  file[2] = new TFile("hist_rootfiles/DetVar/cov_Recomb2.root"); //3
  file[7] = new TFile("hist_rootfiles/DetVar/cov_WMX.root");  //8
  file[3] = new TFile("hist_rootfiles/DetVar/cov_SCE.root");  //4
  file[8] = new TFile("hist_rootfiles/DetVar/cov_WMYZ.root");  //9
  file[4] = new TFile("hist_rootfiles/DetVar/cov_WMdEdx.root"); //5
  file[9] = new TFile("hist_rootfiles/DetVar/cov_LYatt.root"); //10

  TGraph **g1 = new TGraph*[10];
  
  for (Int_t i=0;i!=10;i++){
    TVectorD *v1 = (TVectorD*)file[i]->Get(Form("vec_mean_%d",i+1));
    TVectorD *v1d = (TVectorD*)file[i]->Get(Form("vec_mean_diff_%d",i+1));
    g1[i] = new TGraph();
    for (Int_t j=0;j!=26;j++){
      double y = ((*v1d)[52+j]+(*v1d)[52+26+j])/((*v1)[52+j]+(*v1)[52+26+j]);
      if (std::isnan(y) || std::isinf(y)) y=0;
      g1[i]->SetPoint(j,j+0.5,y);
    }
  }

  
  g1[3]->Draw("AL*");
  g1[3]->SetMarkerColor(1);
  g1[3]->SetMarkerStyle(20);
  g1[3]->GetYaxis()->SetRangeUser(-1,1);
  g1[3]->GetXaxis()->SetTitle("E_{#nu}^{reco} (x100 MeV)");
  g1[3]->GetYaxis()->SetTitle("Rel. Change");
  for (Int_t i=1;i!=10;i++){
    if (i==4) continue;
    if (i!=3) continue;
    g1[i]->Draw("*Lsame");
    if (i!=9){
      g1[i]->SetLineColor(i+1);
      g1[i]->SetMarkerColor(i+1);
    }else{
      g1[i]->SetLineColor(4+1);
      g1[i]->SetMarkerColor(4+1);
    }
    g1[i]->SetMarkerStyle(20);
  }

  // g1[0]->Draw("Lsame*");
  // g1[6]->Draw("Lsame*");

  TLegend *le1 = new TLegend(0.6,0.6,0.89,0.89);
   le1->AddEntry(g1[0],"LYDown","pl");
   le1->AddEntry(g1[1],"LYRayleigh","pl");
   le1->AddEntry(g1[2],"Recomb2","pl");
   le1->AddEntry(g1[3],"SCE","pl");
   le1->AddEntry(g1[4],"WMdEdx","pl");
   le1->AddEntry(g1[5],"WMThetaXZ","pl");
   le1->AddEntry(g1[6],"WMThetaYZ","pl");
   le1->AddEntry(g1[7],"WMX","pl");
   le1->AddEntry(g1[8],"WMYZ","pl");
   le1->AddEntry(g1[9],"LYatt","pl");
   le1->Draw();

   // g1[8]->Draw("A*");
   
}
