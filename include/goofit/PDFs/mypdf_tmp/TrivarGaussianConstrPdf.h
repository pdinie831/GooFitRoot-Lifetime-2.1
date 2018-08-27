#pragma once

#include <goofit/PDFs/GooPdf.h>


namespace GooFit {

class TrivarGaussianConstrPdf : public GooPdf {
public:
  TrivarGaussianConstrPdf (std::string n, Variable* _xdummy, Observable* _x, Observable* _y, Observable* _z, 
                                                             Variable* m1, Variable* m2, Variable* m3, 
							     Variable* s11, Variable* s22, Variable* s33,
							     Variable* s12, Variable* s13, Variable* s23); 
  __host__ fptype integrate (fptype lo, fptype hi) const; 
  __host__ virtual bool hasAnalyticIntegral () const {return true;} 



private:

};
} // namespace GooFit

