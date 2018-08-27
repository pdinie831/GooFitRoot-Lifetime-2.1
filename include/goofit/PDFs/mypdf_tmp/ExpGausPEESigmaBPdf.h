#pragma once

#include <goofit/PDFs/GooPdf.h>

namespace GooFit {

class ExpGausPEESigmaBPdf : public GooPdf {
public:
  ExpGausPEESigmaBPdf (std::string n, Observable* _x, Observable* _s, Variable* m,  Variable* t, Variable* lo,  Variable* hi); 
   __host__ fptype integrate (fptype lo, fptype hi) const; 
   __host__ virtual bool hasAnalyticIntegral () const {return true;} 

private:


};
} // namespace GooFit
