#include <goofit/PDFs/mypdf/Polynomial2Pdf.h>
#include <goofit/Variable.h>

namespace GooFit {

__device__ fptype device_Polynomial2(fptype *evt, fptype *p, unsigned int *indices) {
    // Structure is nP lowestdegree c1 c2 c3 nO o1

    int numParams    = (indices[0]) + 1;
    int lowestDegree = (indices[1]);

    fptype x   = evt[(indices[2 + (indices[0])])];
    fptype ret = 0;

    for(int i = 2; i < numParams; ++i) {
        ret += (p[(indices[i])]) * pow(x, lowestDegree + i - 2);
    }

    return ret*ret;
}

__device__ fptype device_OffsetPolynomial2(fptype *evt, fptype *p, unsigned int *indices) {
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

__device__ fptype device_MultiPolynomial2(fptype *evt, fptype *p, unsigned int *indices) {
    int numObservables = (indices[(indices[0]) + 1]);
    int maxDegree      = (indices[1]);
    if (numObservables!=3) {
     printf("device_Polynomial2 error: Number of Observables should be = 3!!! numObservables = %d\n",numObservables);
    
     return -100;
    }
 
    fptype x    = (evt[(indices[2 + (indices[0]) ])]); // x, y, z...
    x -= (p[(indices[2])]); ;
    fptype y    = (evt[(indices[2 + (indices[0]) + 1])]); // x, y, z...
    y -= (p[(indices[3])]); ;
    fptype z    = (evt[(indices[2 + (indices[0]) + 2])]); // x, y, z...
    z -= (p[(indices[4])]); ;
    
       int ipar =2 + numObservables;
       int kk = 0;
       int ii = 0;
       int jj = 0;
       fptype func =0;
       for(int i = 0; i <= maxDegree ; ++i) {
       jj = 0;
         for(int j = i; j <= maxDegree ; ++j) {
//	  std::cout<<"func = par["<<ipar<<"]*x^"<<kk<<"*y^"<<jj<<std::endl;
          ii = 0;
          for(int k = j; k <= maxDegree ; ++k) {
//            if ((0 == THREADIDX) && (0 == BLOCKIDX)){
//  	    printf("polynomial2  par = %f	x^%d*y^%d*z^%d = %d\n",(p[(indices[ipar])]),ii,jj,kk);
//            } 
//	   std::cout<<"func = par["<<ipar<<"]*x^"<<ii<<"*y^"<<jj<<"*z^"<<kk<<std::endl;
           func +=(p[(indices[ipar])])*pow(x,ii)*pow(y,jj)*pow(z,kk);

//        if ((0 == THREADIDX) && (0 == BLOCKIDX)){
// 	printf("polynomial2 MaxDegree=%d coefficient = %f   number = %d\n",maxDegree,(p[(indices[ipar])]),ipar-2-2*numObservables);
//        } 
	   
	   ipar++;
           ii = (jj+kk+ii<maxDegree?++ii:0);
	  }
          jj = (jj+kk+ii<maxDegree?++jj:0);
         }
         kk= (jj+kk+ii<maxDegree?++kk:0);
       }
       return  func*func;
 
 
 }

/* __device__ fptype device_MultiPolynomial2(fptype *evt, fptype *p, unsigned int *indices) {
    // Structure is nP, maxDegree, offset1, offset2, ..., coeff1, coeff2, ..., nO, o1, o2, ...

    int numObservables = (indices[(indices[0]) + 1]);
    int maxDegree      = (indices[1]) + 1;
    // Only appears in construction (maxDegree + 1) or (x > maxDegree), so
    // may as well add the one and use >= instead.

    // Technique is to iterate over the full n-dimensional box, skipping matrix elements
    // whose sum of indices is greater than maxDegree. Notice that this is increasingly
    // inefficient as n grows, since a larger proportion of boxes will be skipped.
    int numBoxes = 1;
    
    for(int i = 0; i < numObservables; ++i)
        numBoxes *= maxDegree;

    int coeffNumber = 2 + numObservables; // Index of first coefficient is 2 + nO, not 1 + nO, due to maxDegree. (nO
                                          // comes from offsets.)
    fptype ret = (p[(indices[coeffNumber++])]); // Coefficient of constant term.
    fptype intg= ret*4*atan2(0.0, -1.0);

    for(int i = 1; i < numBoxes;
        ++i) { // Notice skip of inmost 'box' in the pyramid, corresponding to all powers zero, already accounted for.
        fptype currTerm  = 1;
        int currIndex    = i;
        int sumOfIndices = 0;
	fptype currIntg  = 1;
        fptype x1        = 0;
        fptype x2        = 0;

        // if ((gpuDebug & 1) && (THREADIDX == 50) && (BLOCKIDX == 3))
        // if ((BLOCKIDX == internalDebug1) && (THREADIDX == internalDebug2))
        // if ((1 > (int) floor(0.5 + evt[8])) && (gpuDebug & 1) && (paramIndices + debugParamIndex == indices))
        // printf("[%i, %i] Start box %i %f %f:\n", BLOCKIDX, THREADIDX, i, ret, evt[8]);
        for(int j = 0; j < numObservables; ++j) {
	    if(j<=1){
	     x1=-1;
	     x2= 1;
	    }
	    if(j>1){
	     x1=0;
	     x2=atan2(0.0, -1.0);
	    }
            fptype x      = (evt[(indices[2 + (indices[0]) + j])]); // x, y, z...
            fptype offset = (p[(indices[2 + j])]);                          // x0, y0, z0...
            x -= offset;
            int currPower = currIndex % maxDegree;
            currIndex /= maxDegree;
            currTerm *= pow(x, currPower);
            sumOfIndices += currPower;
	    currIntg *= (pow(x2, currPower+1)-pow(x1, currPower+1))/(currPower+1); 
            // if ((gpuDebug & 1) && (THREADIDX == 50) && (BLOCKIDX == 3))
            // if ((BLOCKIDX == internalDebug1) && (THREADIDX == internalDebug2))
            // if ((1 > (int) floor(0.5 + evt[8])) && (gpuDebug & 1) && (paramIndices + debugParamIndex == indices))
            // printf("  [%f -> %f^%i = %f] (%i %i) \n", evt[indices[2 + indices[0] + j]], x, currPower, pow(x,
            // currPower), sumOfIndices, indices[2 + indices[0] + j]);
        }

        // if ((gpuDebug & 1) && (THREADIDX == 50) && (BLOCKIDX == 3))
        // if ((BLOCKIDX == internalDebug1) && (THREADIDX == internalDebug2))
        // printf(") End box %i\n", i);
        // All threads should hit this at the same time and with the same result. No branching.
        if(sumOfIndices >= maxDegree)
            continue;

        fptype coefficient = (p[(indices[coeffNumber++])]); // Coefficient from MINUIT
//         if ((0 == THREADIDX) && (0 == BLOCKIDX)){
//          printf("Polynomial2 MaxDegree=%d coefficient = %f   number = %d\n",maxDegree,coefficient,coeffNumber-3-numObservables);
// 	}
        // if ((gpuDebug & 1) && (THREADIDX == 50) && (BLOCKIDX == 3))
        // if ((BLOCKIDX == internalDebug1) && (THREADIDX == internalDebug2))
        // if ((1 > (int) floor(0.5 + evt[8])) && (gpuDebug & 1) && (paramIndices + debugParamIndex == indices))
        // printf("Box %i contributes %f * %f = %f -> %f\n", i, currTerm, p[indices[coeffNumber - 1]],
        // coefficient*currTerm, (ret + coefficient*currTerm));
        currTerm *= coefficient;
        currIntg *= coefficient;
        ret += currTerm;
	intg+= currIntg;
    }

    // if ((1 > (int) floor(0.5 + evt[8])) && (gpuDebug & 1) && (paramIndices + debugParamIndex == indices))
    // printf("Final Polynomial2: %f\n", ret);

// if (0 > ret/(intg)) ret = 0;
// if (0 > ret) ret = -ret;
// PADUL!!!
// if (ret/(intg)>1) printf("Error in Polynomial2Pdf => %f",ret/(intg));
// if (ret/(intg)<0) printf("Error in Polynomial2Pdf => %f",ret/(intg));
//    return ret/(intg);
//    return ret*ret/(intg*intg);
    return ret*ret;
//    return ret;
}
 */
__device__ device_function_ptr ptr_to_Polynomial2       = device_Polynomial2;
__device__ device_function_ptr ptr_to_OffsetPolynomial2 = device_OffsetPolynomial2;
__device__ device_function_ptr ptr_to_MultiPolynomial2  = device_MultiPolynomial2;

// Constructor for single-variate Polynomial2, with optional zero point.
__host__
Polynomial2Pdf::Polynomial2Pdf(std::string n, Observable _x, std::vector<Variable> weights, unsigned int lowestDegree)
    : GooPdf(n, _x) {
    std::vector<unsigned int> pindices;
    pindices.push_back(lowestDegree);

    for(auto &weight : weights) {
        pindices.push_back(registerParameter(weight));
    }

    GET_FUNCTION_ADDR(ptr_to_Polynomial2);

    initialize(pindices);
}

// Constructor for single-variate Polynomial2, with optional zero point.
__host__ Polynomial2Pdf::Polynomial2Pdf(
    std::string n, Observable _x, std::vector<Variable> weights, Variable x0, unsigned int lowestDegree)
    : GooPdf(n, _x)
    , center(new Variable(x0)) {
    std::vector<unsigned int> pindices;
    pindices.push_back(lowestDegree);

    for(auto &weight : weights) {
        pindices.push_back(registerParameter(weight));
    }

    pindices.push_back(registerParameter(x0));
    GET_FUNCTION_ADDR(ptr_to_OffsetPolynomial2);

    initialize(pindices);
}

// Constructor for multivariate Polynomial2.
__host__ Polynomial2Pdf::Polynomial2Pdf(std::string n,
                                      std::vector<Observable> obses,
                                      std::vector<Variable> coeffs,
                                      std::vector<Variable> offsets,
                                      unsigned int maxDegree)
    : GooPdf(n) {
    unsigned int numParameters = 1;

    // For 1 observable, equal to n = maxDegree + 1.
    // For two, n*(n+1)/2, ie triangular number. This generalises:
    // 3: Pyramidal number n*(n+1)*(n+2)/(3*2)
    // 4: Hyperpyramidal number n*(n+1)*(n+2)*(n+3)/(4*3*2)
    // ...
    for(unsigned int i = 0; i < obses.size(); ++i) {
        registerObservable(obses[i]);
        numParameters *= (maxDegree + 1 + i);
    }

    for(int i = observables.size(); i > 1; --i)
        numParameters /= i;

    while(numParameters > coeffs.size()) {
        char varName[100];
        sprintf(varName, "%s_extra_coeff_%i", getName().c_str(), static_cast<int>(coeffs.size()));

        coeffs.emplace_back(varName, 0.1,-1000.,1000.);

        std::cout << "Warning: " << getName() << " created dummy variable " << varName
                  << "  to account for all terms.\n";
    }

    while(offsets.size() < obses.size()) {
        char varName[100];
        sprintf(varName, "%s_extra_offset_%i", getName().c_str(), static_cast<int>(offsets.size()));
        offsets.emplace_back(varName, 0);
    }

    std::vector<unsigned int> pindices;
    pindices.push_back(maxDegree);

    for(auto &offset : offsets) {
        pindices.push_back(registerParameter(offset));
    }

    for(auto &coeff : coeffs) {
        pindices.push_back(registerParameter(coeff));
    }

    GET_FUNCTION_ADDR(ptr_to_MultiPolynomial2);
    initialize(pindices);
}

__host__ fptype Polynomial2Pdf::integrate(fptype lo, fptype hi) const {
//     // This is *still* wrong. (13 Feb 2013.)
//      unsigned int *indices = host_indices + parameters;
//      fptype lowestDegree   = indices[1];
//  
//      if(center) {
// 	 hi -= host_params[indices[indices[0]]];
// 	 lo -= host_params[indices[indices[0]]];
//      }
//  
//      fptype ret = 0;
//  
//      for(int i = 2; i < indices[0] + (center ? 0 : 1); ++i) {
// 	 fptype powerPlusOne = lowestDegree + i - 2;
// 	 fptype curr	     = pow(hi, powerPlusOne);
// 	 curr -= pow(lo, powerPlusOne);
// 	 curr /= powerPlusOne;
// 	 ret += host_params[indices[i]] * curr;
//      }
//  
//      return ret;
       return 1;
}

__host__ fptype Polynomial2Pdf::getCoefficient(int coef) const {
    // NB! This function only works for single Polynomial2s.
    if(1 != observables.size()) {
        std::cout << "Warning: getCoefficient method of Polynomial2Pdf not implemented for multi-dimensional "
                     "Polynomial2s. Returning zero, which is very likely wrong.\n";
        return 0;
    }

    unsigned int *indices = host_indices + parameters;

    // True function is, say, ax^2 + bx + c.
    // We express this as (a'x^2 + b'x + c')*N.
    // So to get the true coefficient, multiply the internal
    // one by the normalisation. (In non-PDF cases the normalisation
    // equals one, which gives the same result.)

    // Structure is nP lowestdegree c1 c2 c3 nO o1
    if(coef < indices[1])
        return 0; // Less than least power.

    if(coef > indices[1] + (indices[0] - 1))
        return 0; // Greater than max power.

    fptype norm = normalize();
    norm        = (1.0 / norm);

    fptype param = host_params[indices[2 + coef - indices[1]]];
    return norm * param;
}
} // namespace GooFit
