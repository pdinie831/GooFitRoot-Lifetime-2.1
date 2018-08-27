#pragma once

#include <goofit/PDFs/GooPdf.h>

namespace GooFit {

class SimpleCheby2Pdf : public GooPdf {
public:
  SimpleCheby2Pdf (std::string n, Observable _x, Variable p0, Variable p1, Variable xmin, Variable xmax); 
  __host__ fptype integrate (fptype lo, fptype hi) const; 
//  __host__ virtual bool hasAnalyticIntegral () const {return true;} 



private:

};
} // namespace GooFit

