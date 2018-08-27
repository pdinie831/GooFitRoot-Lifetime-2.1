#pragma once

#include <goofit/PDFs/GooPdf.h>

namespace GooFit {

class ExpGausProdBPdf : public GooPdf {
public:
  ExpGausProdBPdf (std::string n, Observable _x, Observable _s, Variable m,  Variable t, Variable ss, Variable ms,  Variable ts,  Variable lo, Variable hi, Variable los, Variable his); 
   __host__ fptype integrate (fptype lo, fptype hi) const; 
   __host__ virtual bool hasAnalyticIntegral () const {return true;} 

private:


};
} // namespace GooFit
