#define MAX_THREADS_PER_BLOCK 512
#define MIN_BLOCKS_PER_MP     20
#include <goofit/PDFs/mypdf/FastBernsteinPdf.h>
#include <goofit/Variable.h>

//   __global__ void
//   __launch_bounds__(MAX_THREADS_PER_BLOCK, MIN_BLOCKS_PER_MP)

namespace GooFit {

/* a struct for storing double numbers */
// struct bernVal {
//   double bernFunc;
//   double bernIntg;
// };


__forceinline__ __device__ fptype device_coeffbinomial(fptype enne, fptype kappa){
 
        fptype factor=1.;
        for(fptype i = 1; i <=kappa; ++i) {
          factor *= (enne+1-i)/i; 
        }	 
 
        if (factor<=0 ){
	 printf("Error in FastBernsteinPdf coeffbinomial=> factor = %f enne=%f kappa=%f",factor,enne,kappa);
         return 0;
	} 
       return factor;
}
__forceinline__ __device__ fptype  device_bernsteinkn_func(fptype x, fptype enne, fptype kappa){
 
   return device_coeffbinomial(enne,kappa)*pow(x,kappa)*pow(1.0-x,enne-kappa);


}
/*  __device__ fptype  device_fastbernsteinkn_intg(fptype x, fptype enne, fptype kappa){
 
//  	if ((52 == THREADIDX) && (0 == BLOCKIDX)){
//       printf("==================================================\n");
//       printf("==================================================\n");
//      }
//       struct bernVal results;
//       results.bernFunc = 0;
//       results.bernIntg = 0;
      if (x<0 || x>1 ){
       printf(" Error in bernsteinkn_intg  x=%5.15f out of range [0,1]\n",x);
       return 0.;
      }
//      if (kappa>enne) return 0;
//      bernkn *= pow(x,kappa) ;
//      bernkn *= pow(1.0-x,enne-kappa) ;
      fptype integbernkn = 0;
      fptype ifactni = 0;
      fptype ifactik = 0;
      
       for(fptype i = kappa; i <=enne ; ++i) {
// n!/(i!(n-i)!)
//         ifactni=1;
//         for(float l = 1; l <=i; ++l) {
//           ifactni *= (maxDegree+1-l)/l; 
//         }	 

        ifactni =  device_coeffbinomial(enne,i);
// i!/(k!(i-k)!)
//         ifactik=1;
//         for(float l = 1; l <=k; ++l) {
//           ifactik *= (i+1-l)/l; 
//         }	
 
        ifactik =  device_coeffbinomial(i,kappa);
//
//        bernkn      += ifactni*ifactik*pow(menuno,i-kappa)*pow(x, i) ;
        integbernkn += ifactni*ifactik*pow(-1.0,i-kappa)/(i+1);
//	if ((52 == THREADIDX) && (0 == BLOCKIDX)){
//          printf("pow(x=%5.15f,i=%5.15f)=%5.15f\n",x,i,pow(x, i ));
//          printf("pow(-1=%5.15f,i-kappa)=%5.15f\n",menuno,pow(menuno,i-kappa));
//          printf("bernsteinkn=%5.15f integral=%5.15f ifactni=%5.15f ifactik=%5.15f \n",bernkn,integbernkn,ifactni,ifactik);
//          printf("bernsteinkn=%f integral=%f kappa=%f i=%f enne=%f ni=%f nk=%f\n",bernkn,integbernkn,kappa,i,enne,ifactni,ifactik);
//        }
       }

       if (integbernkn<=0 ){
//	if ((52 == THREADIDX) && (0 == BLOCKIDX)){
         printf(" Error in bernsteinkn_intg x=%5.15f integral = %5.15f THREADIDX=%d BLOCKIDX=%d\n", x,kappa,enne,integbernkn,THREADIDX,BLOCKIDX);
//	}
       }
//       results.bernFunc = bernkn;
//       results.bernIntg = integbernkn;
       return integbernkn;
}
 *///
//
//
//
__device__ fptype device_FastBernstein(fptype *evt, fptype *p, unsigned int *indices) {
    // Structure is nP lowestdegree c1 c2 c3 nO o1
     
     
//    struct bernVal bernknval;

    int numParams = (indices[0]) ;
    int maxDegree = (indices[1]);

    fptype x   = evt[(indices[2 + (indices[0])])];
    fptype ret = 0.0;
//    fptype integret = 0;
//    fptype bernkn = 0;
//    fptype integbernkn = 0;
    fptype intg_1 = 0.0;
    int    ipar = 2;
//    fptype ifactni=1;
//    fptype ifactik=1;
    fptype xmin=p[(indices[numParams-1])];
    fptype xmax=p[(indices[numParams])];
    x=(x-xmin)/(xmax-xmin);
//     printf("FastBernsteinPdf => limit xmin= %f xmax= %f\n",xmin,xmax);
//     return 0;


//     for(int i = 2; i < numParams; ++i) {
//         ret += (p[(indices[i])]) * pow(x, lowestDegree + i - 2);
//     }
    
      float k;
//      float i;
      for(k = 0; k <=maxDegree; ++k) {
       if (ipar>numParams-1){
        printf("Error in FastBernsteinPdf => ipar=%d > numParams=%d\n",ipar,numParams);
        return 0;
       }
       ret      += (p[(indices[ipar])]) * device_bernsteinkn_func(x,maxDegree,k);
       intg_1   += (p[(indices[ipar])]);
//       integret += (p[(indices[ipar])]) * device_bernsteinkn_intg(x,maxDegree,k);
//       printf("FastBernsteinPdf => %f integral = %f k=%d numparam=%d par=%f\n",ret,integret,k,numParams,(p[(indices[ipar])]));
       ipar++;
      }
      intg_1=(maxDegree+1)/intg_1/(xmax-xmin);
//       printf("FastBernsteinPdf => %f int = %f\n",ret,integret);
    if(ret<1.E-30)  return 1.E-30;
    return ret*intg_1;
//    return 0.;
}

/* __device__ fptype device_OffsetFastBernstein(fptype *evt, fptype *p, unsigned int *indices) {
    int numParams    = (indices[0]);
    int lowestDegree = (indices[1]);

    fptype x = evt[(indices[2 + numParams])];
    x -= (p[(indices[numParams])]);
    fptype ret = 0;

    for(int i = 2; i < numParams; ++i) {
        ret += (p[(indices[i])]) * pow(x, lowestDegree + i - 2);
    }

    return ret*ret;
}
*/
__device__ fptype device_MultiFastBernstein(fptype *evt, fptype *p, unsigned int *indices) {
//       if ((0 == THREADIDX) && (0 == BLOCKIDX)){
//        printf("==================================================\n");
//       }
// //     struct bernVal bernknvalx;
//     struct bernVal bernknvaly;
//     struct bernVal bernknvalz;
    int numObservables = (indices[(indices[0]) + 1]);
    int maxDegree1      = (indices[1]);
    int maxDegree2      = (indices[2]);
    int maxDegree3      = (indices[3]);
//      if ((0 == THREADIDX) && (0 == BLOCKIDX)){
//      printf("MultiFastBernstein 0=%d 0+1=%d 0+2=%d 0+3=%d\n",indices[(indices[0])],indices[(indices[0])+1],indices[(indices[0])+2],indices[(indices[0])+3]);
//      printf("MultiFastBernstein numObservables=%d maxDegree1=%d maxDegree2=%d maxDegree3=%d\n",numObservables,indices[1],indices[2],indices[3]);
//      }
    if (numObservables!=3) {
     printf("device_MultiFastBernstein error: Max Number of Observables is = 3!!! numObservables = %d\n",numObservables);
     return -100;
    }
 
    fptype x    = (evt[(indices[2 + (indices[0]) ])]); // x, y, z...
    fptype y    = (evt[(indices[2 + (indices[0]) + 1])]); // x, y, z...
    fptype z    = (evt[(indices[2 + (indices[0]) + 2])]); // x, y, z...
//      if ((0 == THREADIDX) && (0 == BLOCKIDX)){
//      printf("MultiFastBernstein x=%5.15f y=%5.15f z=%5.15f %d %d %d\n",x,y,z,numObservables,indices[1],indices[2],indices[3]);
//      }
    fptype xmin = (p[(indices[4 ])]);
    fptype xdif = (p[(indices[5 ])])-(p[(indices[4 ])]);
    x=(x-xmin)/xdif;
    fptype ymin = (p[(indices[6])]);
    fptype ydif = (p[(indices[7])])-(p[(indices[6])]);
    y=(y-ymin)/ydif;
    fptype zmin = (p[(indices[8])]);
    fptype zdif = (p[(indices[9])])-(p[(indices[8])]);
    z=(z-zmin)/zdif;
    
//        if ((0 == THREADIDX) && (0 == BLOCKIDX)){
// 	printf("MultiFastBernstein xmin=%5.15f xmax = %5.15f\n",xmin,xdif);
// 	printf("MultiFastBernstein ymin=%5.15f ymax = %5.15f\n",ymin,ydif); 
// 	printf("MultiFastBernstein zmin=%5.15f zmax = %5.15f\n",zmin,zdif);
// 	printf("MultiFastBernstein [0,1] x=%5.15f y=%5.15f z=%5.15f \n",x,y,z);
//        
//        }
       int ipar =4 + 2*numObservables;
//       int kk = 0;
//       int ii = 0;
//       int jj = 0;
       fptype func =0;
       fptype intg_1 =0;
       for(int i = 0; i <= maxDegree1 ; ++i) {
//       jj = 0;
         for(int j = 0; j <= maxDegree2 ; ++j) {
//	  std::cout<<"func = par["<<ipar<<"]*x^"<<kk<<"*y^"<<jj<<std::endl;
//          ii = 0;
          for(int k = 0; k <= maxDegree3 ; ++k) {
//	   std::cout<<"func = par["<<ipar<<"]*x^"<<ii<<"*y^"<<jj<<"*z^"<<kk<<std::endl;
    	   fptype bernknvalx =  device_coeffbinomial(maxDegree1,i)*pow(x,i)*pow(1.0-x,maxDegree1-i);
    	   fptype bernknvaly =  device_coeffbinomial(maxDegree2,j)*pow(y,j)*pow(1.0-y,maxDegree2-j);
    	   fptype bernknvalz =  device_coeffbinomial(maxDegree3,k)*pow(z,k)*pow(1.0-z,maxDegree3-k);
//     	   fptype bernknvalx =  device_bernsteinkn_func(x,maxDegree1,i);
//  	   fptype bernknvaly =  device_bernsteinkn_func(y,maxDegree2,j);
//  	   fptype bernknvalz =  device_bernsteinkn_func(z,maxDegree3,k);
//	   fptype bernknintx =  device_bernsteinkn_intg(x,maxDegree1,i);
//	   fptype bernkninty =  device_bernsteinkn_intg(y,maxDegree2,j);
//	   fptype bernknintz =  device_bernsteinkn_intg(z,maxDegree3,k);
//            func +=(p[(indices[ipar])])*bernknvalx*bernknvaly*bernknvalz;
//            intg +=(p[(indices[ipar])])*bernknintx*bernkninty*bernknintz;
           func   +=(p[(indices[ipar])])*bernknvalx*bernknvaly*bernknvalz;
           intg_1 +=(p[(indices[ipar])]);
// 	    if ((0 == THREADIDX) && (0 == BLOCKIDX)){
//  	     printf("MultiFastBernstein  par = %f       \n",(p[(indices[ipar])]));
// 	     printf("MultiFastBernstein  par = %f       B_(%d,%d,%d) = %f intg=%f\n",(p[(indices[ipar])]),ii,jj,kk,bernknvalx,bernknintx);
// 	    } 

//        if ((0 == THREADIDX) && (0 == BLOCKIDX)){
// 	printf("MultiFastBernstein MaxDegree=%d coefficient = %f   number = %d\n",maxDegree,(p[(indices[ipar])]),ipar-2-2*numObservables);
//        } 
	   
	   ipar++;
//           ii = (jj+kk+ii<maxDegree?++ii:0);
	  }
//          jj = (jj+kk+ii<maxDegree?++jj:0);
	  
	 
         }
//         kk= (jj+kk+ii<maxDegree?++kk:0);
       }
//       return  func*func;
//       return  func/(intg);
//      return  func/(intg)/xdif/ydif/zdif;
//      return  func/(intg)/xdif/ydif/zdif;
      if(func<1.E-30)  return 1.E-30;
      intg_1 = (1.0+maxDegree1)*(1.0+maxDegree2)*(1.0+maxDegree3)/intg_1;
      intg_1 = intg_1/(xdif*ydif*zdif);
      fptype ret = func*intg_1;
//      if(ret>1.0) return 0.99999;
      return  ret;
 }
// __device__ fptype device_MultiFastBernstein(fptype *evt, fptype *p, unsigned int *indices) {
//     // Structure is nP, maxDegree, offset1, offset2, ..., coeff1, coeff2, ..., nO, o1, o2, ...
// 
//     struct bernVal bernknval;
//     int numObservables = (indices[(indices[0]) + 1]);
//     int maxDegree      = (indices[1])+1;
//  
// //     printf("MultiFastBernstein CosL xmin = %f   xmax = %f\n",(p[(indices[4 + 0])]),(p[(indices[4 + 0]+1)]));
// //     printf("MultiFastBernstein CosK xmin = %f   xmax = %f\n",(p[(indices[4 + 1])]),(p[(indices[4 + 1]+1)]));
// //     printf("MultiFastBernstein Phi  xmin = %f   xmax = %f\n",(p[(indices[4 + 2])]),(p[(indices[4 + 2]+1)]));
// //     return 0;
//     
// //    printf("MultiFastBernstein  maxDegree = %d   numObservables = %d\n",maxDegree,numObservables);
// //    return 0;
//     
//     // Only appears in construction (maxDegree + 1) or (x > maxDegree), so
//     // may as well add the one and use >= instead.
// 
//     // Technique is to iterate over the full n-dimensional box, skipping matrix elements
//     // whose sum of indices is greater than maxDegree. Notice that this is increasingly
//     // inefficient as n grows, since a larger proportion of boxes will be skipped.
//     int numBoxes = 1;
//     
//     for(int i = 0; i < numObservables; ++i)
//         numBoxes *= maxDegree;
// 
// //     int coeffNumber = 2 + numObservables; // Index of first coefficient is 2 + nO, not 1 + nO, due to maxDegree. (nO
// //                                           // comes from offsets.)
//     int coeffNumber = 2 + 2*numObservables-1; // Index of first coefficient is 2 + nO, not 1 + nO, due to maxDegree. (nO
//                                           // comes from offsets.)
// //    fptype ret = (p[(indices[coeffNumber++])]); // Coefficient of constant term.
// //         printf("MultiFastBernstein  coefficient = %f   number = %d\n",(p[(indices[coeffNumber])]),(indices[coeffNumber]));
// //         printf("MultiFastBernstein  coefficient = %f   number = %d\n",(p[(indices[coeffNumber]+85)]),(indices[coeffNumber])+85);
// // 	return 0;
//     fptype ret = 0;
//     fptype intg= 0;
//     fptype coefficient = 0;
//     for(int i = 1; i <= numBoxes;
//         ++i) { // Notice skip of inmost 'box' in the pyramid, corresponding to all powers zero, already accounted for.
//         fptype currTerm  = 1;
//         int currIndex    = i;
//         int sumOfIndices = 0;
// 	fptype currIntg  = 1;
// 
//         // if ((gpuDebug & 1) && (THREADIDX == 50) && (BLOCKIDX == 3))
//         // if ((BLOCKIDX == internalDebug1) && (THREADIDX == internalDebug2))
//         // if ((1 > (int) floor(0.5 + evt[8])) && (gpuDebug & 1) && (paramIndices + debugParamIndex == indices))
//         // printf("[%i, %i] Start box %i %f %f:\n", BLOCKIDX, THREADIDX, i, ret, evt[8]);
//         for(int j = 0; j < numObservables; ++j) {
// //           if(sumOfIndices >= maxDegree) continue;
//             fptype x    = (evt[(indices[2 + (indices[0]) + j])]); // x, y, z...
//  	    fptype xmin = (p[(indices[4 + j])]);
//  	    fptype xmax = (p[(indices[4 + j]+1)]);
//             x=(x-xmin)/(xmax-xmin);
// //            fptype offset = (p[(indices[2 + j])]);                          // x0, y0, z0...
// //            x -= offset;
//             int currPower = currIndex % maxDegree;
//             currIndex /= maxDegree;
// //            currTerm *= pow(x, currPower);
//             bernknval = device_bernsteinkn(x,maxDegree,currPower);
//             currTerm *= bernknval.bernFunc;
// 	    currIntg *= bernknval.bernIntg;
// //            currTerm *= pow(x, currPower);
//             sumOfIndices += currPower;
//             // if ((gpuDebug & 1) && (THREADIDX == 50) && (BLOCKIDX == 3))
//             // if ((BLOCKIDX == internalDebug1) && (THREADIDX == internalDebug2))
//             // if ((1 > (int) floor(0.5 + evt[8])) && (gpuDebug & 1) && (paramIndices + debugParamIndex == indices))
//             // printf("  [%f -> %f^%i = %f] (%i %i) \n", evt[indices[2 + indices[0] + j]], x, currPower, pow(x,
//             // currPower), sumOfIndices, indices[2 + indices[0] + j]);
//         }
// 
//         // if ((gpuDebug & 1) && (THREADIDX == 50) && (BLOCKIDX == 3))
//         // if ((BLOCKIDX == internalDebug1) && (THREADIDX == internalDebug2))
//         // printf(") End box %i\n", i);
//         // All threads should hit this at the same time and with the same result. No branching.
//         if(sumOfIndices >= maxDegree) continue;
// 
//         coefficient = (p[(indices[coeffNumber++])]); // Coefficient from MINUIT
//         // if ((gpuDebug & 1) && (THREADIDX == 50) && (BLOCKIDX == 3))
//         // if ((BLOCKIDX == internalDebug1) && (THREADIDX == internalDebug2))
//         // if ((1 > (int) floor(0.5 + evt[8])) && (gpuDebug & 1) && (paramIndices + debugParamIndex == indices))
//         // printf("Box %i contributes %f * %f = %f -> %f\n", i, currTerm, p[indices[coeffNumber - 1]],
//         // coefficient*currTerm, (ret + coefficient*currTerm));
// 	 if ((0 == THREADIDX) && (0 == BLOCKIDX)){
// 	  printf("MultiFastBernstein MaxDegree=%d coefficient = %f   number = %d\n",maxDegree,(p[(indices[coeffNumber])]),coeffNumber-2-2*numObservables);
//       } 
//         currTerm *= coefficient;
//         currIntg *= coefficient;
//         ret += currTerm;
// 	intg+= currIntg;
//     }
//     // if ((1 > (int) floor(0.5 + evt[8])) && (gpuDebug & 1) && (paramIndices + debugParamIndex == indices))
//     // printf("Final FastBernstein: %f\n", ret);
// 
// // if (0 > ret/(intg)) ret = 0;
// // if (0 > ret) ret = -ret;
// // PADUL!!!
// // if (ret/(intg)>1) printf("Error in FastBernsteinPdf => %f",ret/(intg));
// // if (ret/(intg)<0) printf("Error in FastBernsteinPdf => %f",ret/(intg));
// //    return ret/(intg);
// //    return ret*ret/(intg*intg);
//     return ret/intg;
// //    return ret;
// }

__device__ device_function_ptr ptr_to_FastBernstein       = device_FastBernstein;
//__device__ device_function_ptr ptr_to_OffsetFastBernstein = device_OffsetFastBernstein;
__device__ device_function_ptr ptr_to_MultiFastBernstein  = device_MultiFastBernstein;

// Constructor for single-variate FastBernstein, with optional zero point.
// __host__ FastBernsteinPdf::FastBernsteinPdf(std::string n, Observable _x, std::vector<Variable> weights, unsigned int lowestDegree)
//     : GooPdf(n, _x) {
//     std::vector<unsigned int> pindices;
//     pindices.push_back(lowestDegree);
// 
//     for(auto &weight : weights) {
//         pindices.push_back(registerParameter(weight));
//     }
// 
//     GET_FUNCTION_ADDR(ptr_to_FastBernstein);
// 
//     initialize(pindices);
// }

//Constructor for single-variate FastBernstein, with optional zero point.
__host__ FastBernsteinPdf::FastBernsteinPdf(std::string n, Observable _x, std::vector<Variable> weights,std::vector<Variable> limits, unsigned int maxDegree)
    : GooPdf(n, _x) {
    std::vector<unsigned int> pindices;
    pindices.push_back(maxDegree);

    for(auto &weight : weights) {
        pindices.push_back(registerParameter(weight));
    }
    for(auto &limit : limits) {
        pindices.push_back(registerParameter(limit));
    }

     GET_FUNCTION_ADDR(ptr_to_FastBernstein);
//    GET_FUNCTION_ADDR(ptr_to_OffsetFastBernstein);

    initialize(pindices);
}
// 
 // Constructor for multivariate FastBernstein.
 __host__ FastBernsteinPdf::FastBernsteinPdf(std::string n,
				       std::vector<Observable> obses,
				       std::vector<Variable> coeffs,
				       std::vector<Variable> limits,
				       unsigned int maxDegree1,
				       unsigned int maxDegree2,
				       unsigned int maxDegree3 )
        : GooPdf(n) {
     unsigned int numParameters = 1;
 
     // For 1 observable, equal to n = maxDegree + 1.
     // For two, n*(n+1)/2, ie triangular number. This generalises:
     // 3: Pyramidal number n*(n+1)*(n+2)/(3*2)
     // 4: Hyperpyramidal number n*(n+1)*(n+2)*(n+3)/(4*3*2)
     // ...
     for(unsigned int i = 0; i < obses.size(); ++i) {
	 registerObservable(obses[i]);
//	 numParameters *= (maxDegree + 1 + i);
     }
//  
//      for(int i = observables.size(); i > 1; --i)
// 	 numParameters /= i;
 
//     int j=1;
//     numParameters = pow((maxDegree+1),coeffs.size());
     numParameters = (maxDegree1+1)*(maxDegree2+1)*(maxDegree3+1);
     while(numParameters > coeffs.size()) {
	 char varName[100];
	 sprintf(varName, "%s_extra_coeff_%i", getName().c_str(), static_cast<int>(coeffs.size()));
 
	 coeffs.emplace_back(varName, 10.,0.00001,0.,500.);
 
	 std::cout << "Warning: " << getName() << " created dummy variable " << varName
		   << "  to account for all terms.\n";
     }
 
     while(limits.size() < 2*obses.size()) {
	 char varName[100];
	 sprintf(varName, "%s_extra_limits_%i", getName().c_str(), static_cast<int>(limits.size()));
	 limits.emplace_back(varName, 0);
     }
 
     std::vector<unsigned int> pindices;
     pindices.push_back(maxDegree1);
     pindices.push_back(maxDegree2);
     pindices.push_back(maxDegree3);
 
     for(auto &limit : limits) {
	 pindices.push_back(registerParameter(limit));
     }
 
     for(auto &coeff : coeffs) {
	 pindices.push_back(registerParameter(coeff));
     }
 
     GET_FUNCTION_ADDR(ptr_to_MultiFastBernstein);
     initialize(pindices);
 }
//
 __host__ fptype FastBernsteinPdf::integrate(fptype lo, fptype hi) const {
       return 1.0;
 }

} // namespace GooFit
