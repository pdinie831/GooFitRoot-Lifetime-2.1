#pragma once

#include <goofit/PDFs/GooPdf.h>

namespace GooFit {

class ErfcMassPdf : public GooPdf {
public:
  ErfcMassPdf (std::string n, Observable _x,Variable p0, Variable p1, Variable p2, Variable p3) ; 

private:

};
} // namespace GooFit
