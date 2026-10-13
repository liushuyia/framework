void test_fake7(){
//========= Macro generated from object: /
//========= by ROOT version6.16/00
   
   Double_t _fx3001[25] = {
   50,
   150,
   250,
   350,
   450,
   550,
   650,
   750,
   850,
   950,
   1050,
   1150,
   1250,
   1350,
   1450,
   1550,
   1650,
   1750,
   1850,
   1950,
   2050,
   2150,
   2250,
   2350,
   2450};
   Double_t _fy3001[25] = {
   0,
   0,
   0.733103,
   1.35576,
   0.994291,
   0.809163,
   1.62622,
   1.36362,
   0.675926,
   1.30533,
   0.493868,
   0.477019,
   1.0353,
   1.569,
   0.458655,
   1.02011,
   0.653014,
   0.44936,
   0.931196,
   0.38886,
   0.695726,
   1.06825,
   1.11292,
   2.18271,
   1.2687};
   Double_t _felx3001[25] = {
   50,
   50,
   50,
   50,
   50,
   50,
   50,
   50,
   50,
   50,
   50,
   50,
   50,
   50,
   50,
   50,
   50,
   50,
   50,
   50,
   50,
   50,
   50,
   50,
   50};
   Double_t _fely3001[25] = {
   0,
   0,
   0.518382,
   0.479333,
   0.375807,
   0.286082,
   0.346711,
   0.330727,
   0.225309,
   0.337034,
   0.201621,
   0.21333,
   0.312156,
   0.405113,
   0.229327,
   0.360664,
   0.292037,
   0.259438,
   0.416443,
   0.274965,
   0.401678,
   0.534125,
   0.556462,
   0.891087,
   0.732485};
   Double_t _fehx3001[25] = {
   50,
   50,
   50,
   50,
   50,
   50,
   50,
   50,
   50,
   50,
   50,
   50,
   50,
   50,
   50,
   50,
   50,
   50,
   50,
   50,
   50,
   50,
   50,
   50,
   50};

   for (Int_t i=0;i!=25;i++){
     _fehx3001[i] = 0;
     _felx3001[i] = 0;
   }
   
   Double_t _fehy3001[25] = {
   0,
   0,
   0.518382,
   0.479333,
   0.375807,
   0.286082,
   0.346711,
   0.330727,
   0.225309,
   0.337034,
   0.201621,
   0.21333,
   0.312156,
   0.405113,
   0.229327,
   0.360664,
   0.292037,
   0.259438,
   0.416443,
   0.274965,
   0.401678,
   0.534125,
   0.556462,
   0.891087,
   0.732485};
   TGraphAsymmErrors *grae = new TGraphAsymmErrors(25,_fx3001,_fy3001,_felx3001,_fehx3001,_fely3001,_fehy3001);
   grae->SetName("");
   grae->SetTitle("");
   grae->SetFillStyle(1000);
   grae->SetLineWidth(2);
   grae->SetMarkerStyle(20);
   grae->SetMarkerSize(1.5);
   
   TH1F *Graph_Graph3001 = new TH1F("Graph_Graph3001","",100,0,2750);
   Graph_Graph3001->SetMinimum(0);
   Graph_Graph3001->SetMaximum(3.38117);
   Graph_Graph3001->SetDirectory(0);
   Graph_Graph3001->SetStats(0);
   Graph_Graph3001->SetLineWidth(2);
   Graph_Graph3001->SetMarkerStyle(20);
   Graph_Graph3001->GetXaxis()->SetNdivisions(509);
   Graph_Graph3001->GetXaxis()->SetLabelFont(132);
   Graph_Graph3001->GetXaxis()->SetLabelOffset(0.01);
   Graph_Graph3001->GetXaxis()->SetLabelSize(0.08);
   Graph_Graph3001->GetXaxis()->SetTitleSize(0.08);
   Graph_Graph3001->GetXaxis()->SetTitleOffset(0.95);
   Graph_Graph3001->GetXaxis()->SetTitleFont(132);
   Graph_Graph3001->GetYaxis()->SetNdivisions(509);
   Graph_Graph3001->GetYaxis()->SetLabelFont(132);
   Graph_Graph3001->GetYaxis()->SetLabelOffset(0.01);
   Graph_Graph3001->GetYaxis()->SetLabelSize(0.08);
   Graph_Graph3001->GetYaxis()->SetTitleSize(0.08);
   Graph_Graph3001->GetYaxis()->SetTitleOffset(0.95);
   Graph_Graph3001->GetYaxis()->SetTitleFont(132);
   Graph_Graph3001->GetZaxis()->SetLabelFont(132);
   Graph_Graph3001->GetZaxis()->SetLabelSize(0.08);
   Graph_Graph3001->GetZaxis()->SetTitleSize(0.08);
   Graph_Graph3001->GetZaxis()->SetTitleOffset(1.2);
   Graph_Graph3001->GetZaxis()->SetTitleFont(132);
   grae->SetHistogram(Graph_Graph3001);
   
   grae->Draw("A*");

   TF1 *f1 = new TF1("f1","1.0+[0]*sin([1]*x+[2])",0,2500);
   //TF1 *f1 = new TF1("f1","[0]",0,2500);
   f1->SetParameter(0,1);
   f1->SetParameter(1,0.02);
   f1->SetParameter(2,0);
   //   f1->SetParameter(2,3000);
   //   f1->Draw("Lsame");
   grae->Fit(f1,"","+",250,2450);
}
