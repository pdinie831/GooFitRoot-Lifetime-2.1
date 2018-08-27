#pragma once

#include <goofit/PDFs/GooPdf.h>

namespace GooFit {

class SigmoidGausPdf : public GooPdf {
public:
  SigmoidGausPdf (std::string n, Observable* _x, Variable* p0, Variable* p1, Variable* p2, Variable* p3, Variable* p4,Variable* mean,Variable* sigma) ; 

private:

};
} // namespace GooFit
