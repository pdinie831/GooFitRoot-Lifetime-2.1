#include "SimpleCheby2NPdf.h"

EXEC_TARGET fptype device_SimpleCheby2N (fptype* evt, fptype* p, unsigned int* indices) {
  fptype xx   = evt[indices[2 + indices[0]]];
  fptype p0   = p[indices[1]];
  fptype p1   = p[indices[2]];
  fptype xmin = p[indices[3]];
  fptype xmax = p[indices[4]];
  fptype x    = -1+2*(xx-xmin)/(xmax-xmin); 

  fptype ret = 1+p0*x+p1*(2*x*x-1);

  //if ((0 == THREADIDX) && (0 == BLOCKIDX)) cuPrintf("SimpleCheby2N Values %f %i %i %f %f %i\n", x, indices[1], indices[2], p0, p1, callnumber); 
  //cuPrintf("device_SimpleCheby2N %f %i %i %f %f %i %p %f\n", x, indices[1], indices[2], p0, p1, callnumber, indices, ret); 
  //if ((0 == THREADIDX) && (0 == BLOCKIDX))
  //printf("device_SimpleCheby2N %f %f %f %i %f\n", x, p0, p1, callnumber, ret);     
//  if ((0 == THREADIDX) && (0 == BLOCKIDX))
//   printf("device_SimpleCheby2N x=%f p0=%f p1=%f xmin=%f xmax=%f ret=%f\n", x, p0, p1,xmin,xmax, ret);

  return ret; 
}

MEM_DEVICE device_function_ptr ptr_to_SimpleCheby2N = device_SimpleCheby2N; 

__host__ SimpleCheby2NPdf::SimpleCheby2NPdf (std::string n, Observable* _x, Variable* p0, Variable* p1,Variable* xmin, Variable* xmax ) 
  : GooPdf(_x, n) 
{
  std::vector<unsigned int> pindices;
  pindices.push_back(registerParameter(p0));
  pindices.push_back(registerParameter(p1));
  pindices.push_back(registerParameter(xmin));
  pindices.push_back(registerParameter(xmax));
  GET_FUNCTION_ADDR(ptr_to_SimpleCheby2N);
  initialise(pindices); 
}

__host__ fptype SimpleCheby2NPdf::integrate (fptype lo, fptype hi) const {
  
  unsigned int* indices = host_indices+parameters; 
  fptype p0 = host_params[indices[1]];
  fptype p1 = host_params[indices[2]];
  
  fptype xmin =host_params[indices[3]]; 
  fptype xmax =host_params[indices[4]];
  
  fptype xMinL  = -1+2*(xmin-lo)/(xmax-xmin); 
  fptype xMaxL  =  1+2*(hi-xmax)/(xmax-xmin); 

//  fptype  xdelta1 = xMaxL - xMinL;
//  fptype  xdelta2 = xMaxL*xMaxL - xMinL*xMinL;
//  fptype  xdelta3 = xMaxL*xMaxL*xMaxL - xMinL*xMinL*xMinL;
  
  
  fptype xint = 0.5*(hi-lo)*((xMaxL + p0*xMaxL*xMaxL/2+p1*(2*xMaxL*xMaxL*xMaxL/3-xMaxL))- (xMinL +  p0*xMinL*xMinL/2+p1*(2*xMinL*xMinL*xMinL/3-xMinL)));

//   printf("integ_SimpleCheby2N p0=%f p1=%f xmin=%f xmax=%f lo=%f hi=%f int=%f\n",  p0, p1,xmin,xmax, lo,hi,xint);

  return xint; 
}

