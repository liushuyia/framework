#include "WCPLEEANA/GPKernel.h"
#include "WCPLEEANA/GPRegressor.h"
#include "TFile.h"
#include "TVectorD.h"
#include "TRandom3.h"
#include "TGraphErrors.h"

std::vector<double> convert_to_std_vector(TVectorD x){
  std::vector<double> xret;
  for(int i=0; i<x.GetNrows(); i++){
    xret.push_back(x[i]);
  }
  return xret;
}

int main(){

//Obsoleted, needs to be revised
/*
  std::vector<double> pars = {1., 1.};
  RationalQuadraticKernel kern = RationalQuadraticKernel(pars);
  
  // prepare data points
  double error = 0.25;
  std::vector<GPPoint> x;
  std::vector<double> xd, yd, eyd;
  TRandom3 rt; rt.SetSeed();
  for(int i=0; i<10; i++){
    double xt= i + 1.0;
    x.emplace_back(xt);
    xd.push_back(xt);
    double eyt= rt.Gaus(0,error);
    yd.push_back(5*sin(xt) + eyt);
    eyd.push_back(error);
  }
  TVectorD y(yd.size(), &yd[0]);
  // prepare new points
  std::vector<GPPoint> x_new;
  std::vector<double> xd_new, yd_new;
  for(int i=0; i<100; i++){
    double xt= i*0.1 + 1.0;
    x_new.emplace_back(xt);
    xd_new.push_back(xt);
    yd_new.push_back(5*sin(xt));
  }
  // prepare noise term (diagonal in cov)
  std::vector<double> nd(yd.size(), error*error);
  TVectorD noise(nd.size(), &nd[0]);
  TVectorD* noise_ptr = &noise;

  GPRegressor reg = GPRegressor(kern, true, noise_ptr);
  reg.Fit(x, y);
  reg.Predict(x_new);
  
  TVectorD yhat = reg.PosteriorMean();
  TVectorD ystd = reg.PosteriorStd();

  auto gh_meas  = new TGraphErrors(xd.size(), xd.data(), yd.data(), 0, eyd.data());
  auto yd_hat = convert_to_std_vector(yhat);
  auto yd_std = convert_to_std_vector(ystd);
  auto gh_pred_gp = new TGraphErrors(xd_new.size(), xd_new.data(), yd_hat.data(), 0, yd_std.data());
  gh_pred_gp->SetLineColor(4);
  gh_pred_gp->SetMarkerColor(4);
  auto gh_pred_orig = new TGraphErrors(xd_new.size(), xd_new.data(), yd_new.data());
  gh_pred_orig->SetLineColor(2);

  auto ofile = TFile::Open("output_gp.root","recreate");
  gh_meas->Write("gh_meas"); // measured points
  gh_pred_gp->Write("gh_pred_gp"); // gaussian process regression
  gh_pred_orig->Write("gh_pred_orig"); // original function
  ofile->Close();
*/
  return 0;
}
