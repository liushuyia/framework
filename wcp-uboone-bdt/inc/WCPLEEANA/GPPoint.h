#pragma once

#include <math.h>
#include "TMath.h"

class GPPoint
{
public:
  GPPoint(double x[5]) {
    SetX(x);
  }
  
  void SetX(double x[5]) { for (int i=0;i<5;i++) { this->x[i] = x[i]; } }
  
  double* X() { return x; }

  double Mag (GPPoint pt) {
    double* ptx = pt.X();
    double mag = 0;
    for (int i=0;i<5;i++) { mag += TMath::Power(ptx[i]-x[i], 2); }
    return mag;
  }
  
  bool operator==(GPPoint pt) const {
    double* ptx = pt.X();
    return (x[0] == ptx[0]) && (x[1] == ptx[1]) && (x[2] == ptx[2]) && (x[3] == ptx[3]) && (x[4] == ptx[4]);
  };

protected:
  double x[5];
};
