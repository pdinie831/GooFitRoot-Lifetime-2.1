#include <goofit/PDFs/mypdf/BinomialEffiPdf.h>

namespace GooFit {

__device__ fptype device_BinomialEffi(fptype *evt, fptype *p, unsigned int *indices) {
    fptype xRec     = evt[RO_CACHE(indices[2 + RO_CACHE(indices[0])])];
    fptype xGen     = evt[RO_CACHE(indices[2 + RO_CACHE(indices[0]) + 1])];

    fptype mean  = RO_CACHE(p[RO_CACHE(indices[1])]);
    fptype sigma = RO_CACHE(p[RO_CACHE(indices[2])]);

    fptype ret = exp(-0.5 * (x - mean) * (x - mean) / (sigma * sigma));

    return ret;
}

__device__ device_function_ptr ptr_to_BinomialEffi = device_BinomialEffi;

__host__ BinomialEffiPdf::BinomialEffiPdf(std::string n, Observable _r, Observable _g, Variable mean, Variable sigma)
    : GooPdf(n) {
    std::vector<unsigned int> pindices;
    pindices.push_back(registerObservable(_r));
    pindices.push_back(registerObservable(_g));
    pindices.push_back(registerParameter(mean));
    pindices.push_back(registerParameter(mean));
    pindices.push_back(registerParameter(sigma));
    GET_FUNCTION_ADDR(ptr_to_BinomialEffi);
    initialize(pindices);
}

__host__ fptype BinomialEffiPdf::integrate(fptype lo, fptype hi) const {
    return 1.0;
}

} // namespace GooFit
