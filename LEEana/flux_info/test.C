void test(){
  Double_t x = 250;
  Double_t y = 226;
  double_t z = 1031;

  
  TH1F *h1= new TH1F("h1","h1",100,0.95,1.05);
  for (Int_t i=0;i!=10000;i++){
    double ratio = (x+gRandom->Uniform(-1,1))/x
      * (y+gRandom->Uniform(-1,1))/y
      * (z+gRandom->Uniform(-1,1))/z;
    h1->Fill(ratio);
  }
  h1->Draw();
}
