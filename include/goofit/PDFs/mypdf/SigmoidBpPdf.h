#pragma once

#include <goofit/PDFs/GooPdf.h>

namespace GooFit {

class SigmoidBpPdf : public GooPdf {
public:
  SigmoidBpPdf (std::string n, Observable* _x, Variable* p0, Variable* p1, Variable* p2, Variable* p3, Variable* p4, Variable* p5, Variable* p6, Variable* p7) ; 

private:

};
} // namespace GooFi
