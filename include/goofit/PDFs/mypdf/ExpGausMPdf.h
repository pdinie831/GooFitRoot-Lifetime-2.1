#pragma once

#include <goofit/PDFs/GooPdf.h>

namespace GooFit {

class ExpGausMPdf : public GooPdf {
public:
  ExpGausMPdf (std::string n, Observable* _x, Variable* s, Variable* m, Variable* t); 

private:

};
} // namespace GooFit
