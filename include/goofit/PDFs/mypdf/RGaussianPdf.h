#pragma once

#include <goofit/PDFs/GooPdf.h>

namespace GooFit {

class RGaussianPdf : public GooPdf {
public:
  RGaussianPdf (std::string n, Observable _x, Variable m, Variable s); 
  __host__ fptype integrate (fptype lo, fptype hi) const; 
  __host__ virtual bool hasAnalyticIntegral () const {return true;} 



private:

};
} // namespace GooFit
