# GooFitRoot-Lifetime
Private version of GooFit, based on GooFit Branch v2.1, Tag v2.1.1


GooFitRoot-Lifetime should be compiled with root. It has been tested with ROOT 6.04/18,
CUDA 8.0, with devtoolset-3 and python33 on SLC6.9 (GPU Tesla K40m) and Centos 7.5 (GPU GeForce GTX 750 Ti).

============ HOW TO INSTALL =====================================

```bash
git clone -b v2.1  git://github.com/GooFit/GooFit.git --recursive

cd GooFit

git fetch && git fetch --tags

git checkout v2.1.1

cd ..

wget https://github.com/pdinie831/GooFitRoot-Lifetime-2.1/archive/master.zip

unzip master.zip

cp -r GooFitRoot-Lifetime-2.1-master/* GooFit/

cd GooFit/
```

[Usually I install cmake in the following way:
 
```bash
mkdir cmake && wget -qO- "https://cmake.org/files/v3.9/cmake-3.9.4-Linux-x86_64.tar.gz" | tar --strip-components=1 -xz -C cmake
export PATH=`pwd`/cmake/bin:$PATH (or setenv PATH `pwd`/cmake/bin:$PATH, depending on the shell)
```
as suggested in: https://github.com/GooFit/GooFit/blob/virtuals/docs/SYSTEM_INSTALL.md ]

```bash
mkdir build

cd build; cmake ..  -DGOOFIT_PYTHON=OFF

make -j [insert the number of core/cpu available]
```

the new PDF sources are in src/PDFs/mypdf

the built analysis code is in build/example/test3DBp-2016-G21/testGoofit3DBp-2016

the code is in example/test3DBp-2016-G21/testGoofit3DBp-2016.cu



=================================================================


New PDFs have been added in this private version of Goofit.
The most important of these is  ExpGausProdBPdf: the product of two exponentially modified Gaussian 
distributions:

ExpGau(x,sigma;mean,lambda) X ExpGau(sigma;sigmas,means,lambdas).
  

This PDF has two Observable variables defined (x, sigma) and 9  parameters:
(mean,lambda,sigmas,means,lambdas,lo,hi,los,his), where
the parameters lo,hi are limits of the obsevable "x" (lo<x<hi), while los, his are limits of the  
observable "sigma" (los<sigma<his), and they should be fixed in the fit strategy. 


============ ExpGausProdBPdf explained =====================================

In B hadron Lifetime measurements ExpGau(x,sigma;mean,lambda) is the exponential function of the B 
meson decay convolved event-by-event with the resolution model (a gaussian), where:

variable x     == proper time.
variable sigma == proper time error keep as an estimation of the resolution.
parameter mean == mean of the gaussian, usually fixed to 0.
parameter lambda == 1/(Tau*c), where Tau is the Lifetime of the B meson and is a free parameter
of the fit.

It is possible to have different approches to the Lifetime measurements:

- 2D Fit, i.e. (joined) fit to mass and  proper time  of the reconstructed B meson.
- 3D Fit, i.e. (joined) fit to mass, proper time and time resolution of the reconstructed B meson.

ExpGausProdBPdf was written for the 3D Fit, where the ExpGau(x,sigma;mean,lambda) should be 
multiplied by a resolution model PDF. 
The distribution of the error of the proper time (== our estimed resolution) is, generally speaking, a 
sort of bell function with a tail, so it is possible to parametrize using (again...)  an exponentially 
modified Gaussian distribution. Other parameterizations could be Landau and Gamma distributions. 
The product of the two exponentially modified Gaussian distributions should be normalized carefully, 
since the resolution ("sigma") should be considered as a  "Conditional Observable" in the fit.
In ExpGausProdBPdf the convolution of the decay ExpGau(x,sigma;mean,lambda) is normalized only respect
to x (with lo<x<hi), the resolution model ExpGau(sigma;sigmas,means,lambdas) is normalized respect to 
sigma. This has been done avoiding the numerical normalization.

===========================================================================

Others PDFs:

	PolyEffiPdf
	ErfcPolyPdf
	SigmoidPdf
	SigmoidGausPdf
	
have been introduced to parametrize efficiency. 	
