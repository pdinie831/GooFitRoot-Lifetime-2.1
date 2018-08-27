#include "TrivarGaussianConstrPdf.h"

EXEC_TARGET fptype device_TrivarGaussianConstr (fptype* evt, fptype* p, unsigned int* indices) {
//   fptype x = evt[indices[2 + indices[0]]]; 
//   fptype y = evt[indices[3 + indices[0]]]; 
  fptype x      = p[indices[1]];
  fptype y      = p[indices[2]];
  fptype z      = p[indices[3]];
  fptype mean1  = p[indices[4]];
  fptype mean2  = p[indices[5]];
  fptype mean3  = p[indices[6]];
  fptype sigma11 = p[indices[7]];
  fptype sigma22 = p[indices[8]];
  fptype sigma33 = p[indices[9]];
  fptype sigma12 = p[indices[10]];
  fptype sigma13 = p[indices[11]];
  fptype sigma23 = p[indices[12]];
  
//  const fptype pi = 3.14159265359;
  
  fptype dx      = x-mean1;
  fptype dy      = y-mean2;
  fptype dz      = x-mean3;
  
  fptype a11 = sigma33 * sigma22 - sigma23 * sigma23  ;
  fptype a12 = sigma13 * sigma23 - sigma33 * sigma12  ;
  fptype a13 = sigma12 * sigma23 - sigma13 * sigma22  ;

  fptype a22 = sigma33 * sigma11 - sigma13 * sigma13  ;
  fptype a23 = sigma12 * sigma13 - sigma11 * sigma23  ;
					       
  fptype a33 = sigma11 * sigma22 - sigma12 * sigma12  ;
  
  fptype Det  = (sigma11 * a11) + (sigma12 * a12) + (sigma13 * a13);
  
  fptype ret = (dx*dx*a11+dy*dy*a22+dz*dz*a33+2*dx*dy*a12+2*dx*dz*a13+2*y*z*a23)/Det;
  
  if (Det <=0) {
     printf("Error: Det<=0.!!! ==> x=%f  y=%f z=%f sigma11=%f sigma22=%f sigma33=%f sigma12=%f sigma13=%f sigma23=%f \n",x,y,z,sigma11,sigma22,sigma33,sigma12,sigma13,sigma23);
     return 0;
   } 

//  return -0.5*ret; 
  return EXP(-0.5*ret); 
}

MEM_DEVICE device_function_ptr ptr_to_TrivarGaussianConstr = device_TrivarGaussianConstr; 

__host__ TrivarGaussianConstrPdf::TrivarGaussianConstrPdf 
(std::string n, Variable* _xdummy, Observable* _x, Observable* _y, Observable* _z, 
                                   Variable* mean1, Variable* mean2, Variable* mean3, 
				   Variable* sigma11, Variable* sigma22, Variable* sigma33,
				   Variable* sigma12, Variable* sigma13, Variable* sigma23) 
  : GooPdf(_xdummy, n) 
{
//  registerObservable(_x);
//  registerObservable(_y);
  std::vector<unsigned int> pindices;
  pindices.push_back(registerParameter(_x));
  pindices.push_back(registerParameter(_y));
  pindices.push_back(registerParameter(_z));
  pindices.push_back(registerParameter(mean1));
  pindices.push_back(registerParameter(mean2));
  pindices.push_back(registerParameter(mean3));
  pindices.push_back(registerParameter(sigma11));
  pindices.push_back(registerParameter(sigma22));
  pindices.push_back(registerParameter(sigma33));
  pindices.push_back(registerParameter(sigma12));
  pindices.push_back(registerParameter(sigma13));
  pindices.push_back(registerParameter(sigma23));
  GET_FUNCTION_ADDR(ptr_to_TrivarGaussianConstr);
  initialise(pindices); 
}

 __host__ fptype TrivarGaussianConstrPdf::integrate (fptype lo, fptype hi) const {
//   //static const fptype root2 = sqrt(2.);
//   static const fptype rootPi = sqrt(atan2(0.0,-1.0));
//   static const fptype rootPiBy2 = rootPi / root2;
//   
//   unsigned int* indices = host_indices+parameters; 
//   fptype xscale = root2*host_params[indices[2]];
// 
//   /*
//   std::cout << "TrivarGaussianConstr integral: " 
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

