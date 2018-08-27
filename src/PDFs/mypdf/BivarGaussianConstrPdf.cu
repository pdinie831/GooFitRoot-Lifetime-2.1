#include "BivarGaussianConstrPdf.h"

EXEC_TARGET fptype device_BivarGaussianConstr (fptype* evt, fptype* p, unsigned int* indices) {
//   fptype x = evt[indices[2 + indices[0]]]; 
//   fptype y = evt[indices[3 + indices[0]]]; 
  fptype x      = p[indices[1]];
  fptype y      = p[indices[2]];
  fptype mean1  = p[indices[3]];
  fptype sigma1 = p[indices[4]];
  fptype mean2  = p[indices[5]];
  fptype sigma2 = p[indices[6]];
  fptype rho    = p[indices[7]];
  
//  const fptype pi = 3.14159265359;
  fptype rhod = (1-rho*rho);
  if (rhod<=0) {
     printf("Error: rhod<0.!!! ==> x=%f  y=%f mean1=%f sigma1=%f mean2=%f sigma2=%f rho=%f\n",x,y,mean1,sigma1,mean2,sigma2,rho);
     return 0;
   } 
  fptype ret = (x-mean1)*(x-mean1)/(sigma1*sigma1)+(y-mean2)*(y-mean2)/(sigma2*sigma2)-2*rho*(x-mean1)*(y-mean2)/(sigma1*sigma2);
         ret = -0.5*ret/rhod;
//         ret = EXP(ret)/(2*pi*sqrt(rhod));

  //if ((0 == THREADIDX) && (0 == BLOCKIDX)) cuPrintf("BivarGaussianConstr Values %f %i %i %f %f %i\n", x, indices[1], indices[2], mean, sigma, callnumber); 
  //cuPrintf("device_BivarGaussianConstr %f %i %i %f %f %i %p %f\n", x, indices[1], indices[2], mean, sigma, callnumber, indices, ret); 
  //if ((0 == THREADIDX) && (0 == BLOCKIDX))
  //printf("device_BivarGaussianConstr %f %f %f %i %f\n", x, mean, sigma, callnumber, ret);     

//   if (ret<0) {
//      printf("Error: BivarGaussian<=0!!! ==> x=%f  y=%f mean1=%f sigma1=%f mean2=%f sigma2=%f rho=%f ret=%f\n",x,y,mean1,sigma1,mean2,sigma2,rho,ret);
//     ret=0;
//    } 
  return ret; 
}

MEM_DEVICE device_function_ptr ptr_to_BivarGaussianConstr = device_BivarGaussianConstr; 

__host__ BivarGaussianConstrPdf::BivarGaussianConstrPdf (std::string n, Variable* _xdummy, Observable* _x, Observable* _y, Variable* mean1, Variable* sigma1, Variable* mean2, Variable* sigma2, Variable* rho) 
  : GooPdf(_xdummy, n) 
{
//  registerObservable(_x);
//  registerObservable(_y);
  std::vector<unsigned int> pindices;
  pindices.push_back(registerParameter(_x));
  pindices.push_back(registerParameter(_y));
  pindices.push_back(registerParameter(mean1));
  pindices.push_back(registerParameter(sigma1));
  pindices.push_back(registerParameter(mean2));
  pindices.push_back(registerParameter(sigma2));
  pindices.push_back(registerParameter(rho));
  GET_FUNCTION_ADDR(ptr_to_BivarGaussianConstr);
  initialise(pindices); 
}

 __host__ fptype BivarGaussianConstrPdf::integrate (fptype lo, fptype hi) const {
//   //static const fptype root2 = sqrt(2.);
//   static const fptype rootPi = sqrt(atan2(0.0,-1.0));
//   static const fptype rootPiBy2 = rootPi / root2;
//   
//   unsigned int* indices = host_indices+parameters; 
//   fptype xscale = root2*host_params[indices[2]];
// 
//   /*
//   std::cout << "BivarGaussianConstr integral: " 
// 	    << xscale << " "
// 	    << host_params[indices[1]] << " "
// 	    << host_params[indices[2]] << " "
// 	    << ERF((hi-host_params[indices[1]])/xscale) << " "
// 	    << ERF((lo-host_params[indices[1]])/xscale) << " "
// 	    << rootPiBy2*host_params[indices[2]]*(ERF((hi-host_params[indices[1]])/xscale) -
// 						  ERF((lo-host_params[indices[1]])/xscale)) 
// 	    << std::endl; 
//   */
//   return rootPiBy2*host_params[indices[2]]*(ERF((hi-host_params[indices[1]])/xscale) - 
//   					    ERF((lo-host_params[indices[1]])/xscale));

  // Integral over all R. 
//   fptype sigma = host_params[indices[2]];
//   sigma *= root2*rootPi;
//   return sigma; 
return 1.;
}

