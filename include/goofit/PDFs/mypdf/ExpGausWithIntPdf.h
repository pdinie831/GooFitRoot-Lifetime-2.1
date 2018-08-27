#pragma once

#include <goofit/PDFs/GooPdf.h>

namespace GooFit {

class ExpGausWithIntPdf : public GooPdf {
public:
  ExpGausWithIntPdf (std::string n, Observable* _x, Variable* m, Variable* s, Variable* t); 
   __host__ fptype integrate (fptype lo, fptype hi) const; 
   __host__ virtual bool hasAnalyticIntegral () const {return true;} 

private:

};
} // namespace GooFit
