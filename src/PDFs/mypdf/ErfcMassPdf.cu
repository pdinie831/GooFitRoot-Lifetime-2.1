#include <goofit/PDFs/mypdf/ErfcMassPdf.h>

namespace GooFit {

__device__ fptype device_ErfcMass (fptype* evt, fptype* p, unsigned int* indices) {
  fptype x     = evt[indices[2 + indices[0]]]; 
  fptype p0 = p[indices[1]];
  fptype p1 = p[indices[2]];
  fptype p2 = p[indices[3]];
  fptype p3 = p[indices[4]];

  fptype ret = (p2+p3*x)*erfc(p0*(x-p1));

  if (ret<=0) {
       printf("ErfcMass <=0!!!: x = %f , p0 = %f p1 = %f p2 = %f p3 = %f\n",x,p0,p1,p2,p3);
       return 0.;
  }else{
   return ret; 
  } 
}

__device__ device_function_ptr ptr_to_ErfcMass = device_ErfcMass; 

__host__ ErfcMassPdf::ErfcMassPdf (std::string n, Observable _x, Variable p0, Variable p1, Variable p2, Variable p3) 
  : GooPdf(n, _x)
{
  std::vector<unsigned int> pindices;
  pindices.push_back(registerParameter(p0));
  pindices.push_back(registerParameter(p1));
  pindices.push_back(registerParameter(p2));
  pindices.push_back(registerParameter(p3));
  GET_FUNCTION_ADDR(ptr_to_ErfcMass);
  initialize(pindices); 
}

} // namespace GooFit

