#pragma once

#include <goofit/PDFs/GooPdf.h>

namespace GooFit {

class BernsteinPdf : public GooPdf {
  public:
    BernsteinPdf(std::string n, Observable _x, std::vector<Variable> weights,std::vector<Variable> limits, unsigned int maxDegree = 0);
//
    BernsteinPdf(std::string n,
        	 std::vector<Observable> obses,
        	 std::vector<Variable> coeffs,
        	 std::vector<Variable> limits,
        	 unsigned int maxDegree1, 
        	 unsigned int maxDegree2, 
        	 unsigned int maxDegree3);
//		 
    BernsteinPdf(std::string n,
        	 std::vector<Observable> obses,
        	 std::vector<Variable> coeffs,
        	 std::vector<Variable> limits,
        	 unsigned int maxDegree1, 
        	 unsigned int maxDegree2, 
        	 unsigned int maxDegree3,
		 unsigned int dummy);
//
    BernsteinPdf(std::string n,
        	 std::vector<Observable> obses,
        	 std::vector<Variable>   coeffs,
        	 std::vector<Variable>   limits,
        	 std::vector<Variable>   binws,
        	 unsigned int maxDegree1,
        	 unsigned int maxDegree2,
        	 unsigned int maxDegree3);
    __host__ fptype integrate(fptype lo, fptype hi) const override;
//    __host__ virtual bool hasAnalyticIntegral () const {return (1 == observables.size());}
  __host__ virtual bool hasAnalyticIntegral () const {return true;} 
//    __host__ fptype getCoefficient(int coef) const;

  private:
    std::unique_ptr<Variable> center;
};
} // namespace GooFit
