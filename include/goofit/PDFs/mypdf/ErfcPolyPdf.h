#pragma once

#include <goofit/PDFs/GooPdf.h>

namespace GooFit {

class ErfcPolyPdf : public GooPdf {
public:
  ErfcPolyPdf (std::string n, Observable* _x,Variable* p0, Variable* p1, Variable* p2, Variable* p3, Variable* p4) ; 

private:

};
} // namespace GooFit
