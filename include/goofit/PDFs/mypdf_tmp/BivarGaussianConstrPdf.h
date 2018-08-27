#pragma once

#include <goofit/PDFs/GooPdf.h>

namespace GooFit {


class BivarGaussianConstrPdf : public GooPdf {
public:
  BivarGaussianConstrPdf (std::string n, Variable* _xdummy, Observable* _x, Observable* _y, Variable* m1, Variable* s1, Variable* m2, Variable* s2, Variable* r); 
  __host__ fptype integrate (fptype lo, fptype hi) const; 
  __host__ virtual bool hasAnalyticIntegral () const {return true;} 



private:

};
} // namespace GooFit

