#include <goofit/PDFs/mypdf/ErfEffiBpPdf.h>

namespace GooFit {

__device__ fptype device_ErfEffiBp (fptype *evt, fptype *p, unsigned int *indices) {
  fptype x     = evt[indices[2 + indices[0]]]; 
  fptype p0 = p[indices[1]];
  fptype p1 = p[indices[2]];
  fptype p2 = p[indices[3]];
  fptype p3 = p[indices[4]];
  fptype p4 = p[indices[5]];
  fptype p5 = p[indices[6]];
  fptype p6 = p[indices[7]];
  fptype p7 = p[indices[8]];

  if (x<=0){
       printf("ErfEffiBp not defined for x<=0!!!: x = %f",x);
       return 0.;
  }
  fptype ret = (p2+p3*x+ p4*x*x + p5*x*x*x + p6*x*x*x*x + p7*x*x*x*x*x)*erf(p0*(x-p1));


  if (ret<=0) {
       printf("ErfEffiBp <=0!!!: ret = %f x = %f , p0 = %f p1 = %f p2 = %f p3 = %f p4 = %f p5 = %f p6 = %f p7 = %f\n",ret,x,p0,p1,p2,p3,p4,p5,p6,p7);
       return 0.;
  }
  return ret; 
}

__device__ device_function_ptr ptr_to_ErfEffiBp = device_ErfEffiBp; 

__host__ ErfEffiBpPdf::ErfEffiBpPdf (std::string n, Observable _x, Variable p0, Variable p1, Variable p2, Variable p3, Variable p4, Variable p5, Variable p6, Variable p7) 
  : GooPdf(n,_x )
{
  std::vector<unsigned int> pindices;
  pindices.push_back(registerParameter(p0));
  pindices.push_back(registerParameter(p1));
  pindices.push_back(registerParameter(p2));
  pindices.push_back(registerParameter(p3));
  pindices.push_back(registerParameter(p4));
  pindices.push_back(registerParameter(p5));
  pindices.push_back(registerParameter(p6));
  pindices.push_back(registerParameter(p7));
  GET_FUNCTION_ADDR(ptr_to_ErfEffiBp);
  initialize(pindices); 
}
} // namespace GooFit


