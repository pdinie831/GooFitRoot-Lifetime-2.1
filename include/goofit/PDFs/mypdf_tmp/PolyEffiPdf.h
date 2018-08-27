#pragma once

#include <goofit/PDFs/GooPdf.h>

namespace GooFit {
class PolyEffiPdf : public GooPdf {
public:
  PolyEffiPdf (std::string n, Observable* _x, std::vector<Variable*> weights, Variable* x0 = 0, unsigned int lowestDegree = 0); 
  PolyEffiPdf (std::string n, std::vector<Observable*> obses, std::vector<Variable*> coeffs, std::vector<Variable*> offsets, unsigned int maxDegree); 
   __host__ fptype integrate (fptype lo, fptype hi) const; 
  //__host__ virtual bool hasAnalyticIntegral () const {return (1 == observables.size());} 
  __host__ fptype getCoefficient (int coef) const;
  __host__ virtual bool hasAnalyticIntegral () const {return true;}

private:
  Variable* center; 
};
} // namespace GooFit
