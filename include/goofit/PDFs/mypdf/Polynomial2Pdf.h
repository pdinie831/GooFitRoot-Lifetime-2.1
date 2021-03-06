#pragma once

#include <goofit/PDFs/GooPdf.h>

namespace GooFit {

class Polynomial2Pdf : public GooPdf {
  public:
    Polynomial2Pdf(std::string n, Observable _x, std::vector<Variable> weights, unsigned int lowestDegree = 0);
    Polynomial2Pdf(
        std::string n, Observable _x, std::vector<Variable> weights, Variable x0, unsigned int lowestDegree = 0);
    Polynomial2Pdf(std::string n,
                  std::vector<Observable> obses,
                  std::vector<Variable> coeffs,
                  std::vector<Variable> offsets,
                  unsigned int maxDegree);
    __host__ fptype integrate(fptype lo, fptype hi) const override;
//    __host__ virtual bool hasAnalyticIntegral () const {return (1 == observables.size());}
  __host__ virtual bool hasAnalyticIntegral () const {return false;} 
    __host__ fptype getCoefficient(int coef) const;

  private:
    std::unique_ptr<Variable> center;
};
} // namespace GooFit
