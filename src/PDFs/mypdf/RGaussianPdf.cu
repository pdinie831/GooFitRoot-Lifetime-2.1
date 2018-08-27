#include <goofit/PDFs/mypdf/RGaussianPdf.h>

namespace GooFit {

__device__ fptype device_RGaussian (fptype *evt, fptype *p, unsigned int *indices) {
     fptype x = evt[indices[2 + indices[0]]];
     fptype mean = p[indices[1]];
     fptype sigma = p[indices[2]];
//  fptype x     = evt[RO_CACHE(indices[2 + RO_CACHE(indices[0])])];
//  fptype mean  = RO_CACHE(p[RO_CACHE(indices[1])]);
//  fptype sigma = RO_CACHE(p[RO_CACHE(indices[2])]);

  fptype ret = exp(-0.5*(x-mean)*(x-mean)/(sigma*sigma));

  //if ((0 == THREADIDX) && (0 == BLOCKIDX)) cuPrintf("RGaussian Values %f %i %i %f %f %i\n", x, indices[1], indices[2], mean, sigma, callnumber); 
  //cuPrintf("device_RGaussian %f %i %i %f %f %i %p %f\n", x, indices[1], indices[2], mean, sigma, callnumber, indices, ret); 
  //if ((0 == THREADIDX) && (0 == BLOCKIDX))
  //printf("device_RGaussian %f %f %f %i %f\n", x, mean, sigma, callnumber, ret);     


  return ret; 
}

__device__ device_function_ptr ptr_to_RGaussian = device_RGaussian; 

__host__ RGaussianPdf::RGaussianPdf (std::string n, Observable _x, Variable mean, Variable sigma) 
  : GooPdf(n, _x) 
{
  std::vector<unsigned int> pindices;
  pindices.push_back(registerParameter(mean));
  pindices.push_back(registerParameter(sigma));
  GET_FUNCTION_ADDR(ptr_to_RGaussian);
  initialize(pindices); 
}

__host__ fptype RGaussianPdf::integrate (fptype lo, fptype hi) const {
  //static const fptype root2 = sqrt(2.);
  static const fptype rootPi = sqrt(atan2(0.0,-1.0));
  static const fptype rootPiBy2 = rootPi / root2;
  
  unsigned int* indices = host_indices+parameters; 
  fptype xscale = root2*host_params[indices[2]];

  /*
  std::cout << "RGaussian integral: " 
	    << xscale << " "
	    << host_params[indices[1]] << " "
	    << host_params[indices[2]] << " "
	    << ERF((hi-host_params[indices[1]])/xscale) << " "
	    << ERF((lo-host_params[indices[1]])/xscale) << " "
	    << rootPiBy2*host_params[indices[2]]*(ERF((hi-host_params[indices[1]])/xscale) -
						  ERF((lo-host_params[indices[1]])/xscale)) 
	    << std::endl; 
  */
  return rootPiBy2*host_params[indices[2]]*(erf((hi-host_params[indices[1]])/xscale) - 
  					    erf((lo-host_params[indices[1]])/xscale));

  // Integral over all R. 
//   fptype sigma = host_params[indices[2]];
//   sigma *= root2*rootPi;
//   return sigma; 
}
} // namespace GooFit

