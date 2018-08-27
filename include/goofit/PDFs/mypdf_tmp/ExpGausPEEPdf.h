#pragma once

#include <goofit/PDFs/GooPdf.h>

namespace GooFit {

class ExpGausPEEPdf : public GooPdf {
public:
  ExpGausPEEPdf (std::string n, Observable* _x, Observable* _s, Variable* m,  Variable* t); 
   __host__ fptype integrate (fptype lo, fptype hi) const; 
   __host__ virtual bool hasAnalyticIntegral () const {return true;} 

private:


};
} // namespace GooFit

