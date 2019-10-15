/**************************** STOCR.H ************************ 2006-10-21 AF *
*
* This file defines additions to the C++ library of non-uniform random number
* generators for the R-language interface.
*
*
* class StocRBase
* ===============
* This class replaces the base classes for class StochasticLib3 when used for
* the R-language interface.
* Member functions:
*
* double Normal(double m, double s);
* Normal distribution with mean m and standard deviation s.
*
* int32 Hypergeometric (int32 n, int32 m, int32 N);
* Hypergeometric distribution. Taking n items out N, m of which are colored.
*
*
*
* source code:
* ============
* The code for EndOfProgram and FatalError is found in the file userintf.cpp.
* The code for the functions in StochasticLib1 is found in the file stoc1.cpp.
* The code for the functions in StochasticLib2 is found in the file stoc2.cpp.
* The code for the functions in StochasticLib3 is found in the file stoc3.cpp.
* The code for the functions in CWalleniusNCHypergeometric, 
* CMultiWalleniusNCHypergeometric and CMultiWalleniusNCHypergeometricMoments
* is found in the file wnchyppr.cpp.
* The code for the functions in CFishersNCHypergeometric and 
* CMultiFishersNCHypergeometric is found in the file fnchyppr.cpp
* LnFac is found in stoc1.cpp.
* Erf is found in wnchyppr.cpp.
*
*
* Examples:
* =========
*
* Documentation:
* ==============
* The file stocc.htm contains further instructions.
*
* The file distrib.pdf contains definitions of the standard statistic distributions:
* Bernoulli, Normal, Poisson, Binomial, Hypergeometric, Multinomial, MultiHypergeometric.
*
* The file sampmet.pdf contains theoretical descriptions of the methods used
* for sampling from these distributions.
*
* The file nchyp.pdf, available from www.agner.org/random/, contains
* definitions of the univariate and multivariate Wallenius and Fisher's 
* noncentral hypergeometric distributions and theoretical explanations of 
* the methods for calculating and sampling from these.
*
* © 2006 Agner Fog. GNU General Public License www.gnu.org/copyleft/gpl.html
*******************************************************************************/

#ifndef STOC_R_H
#define STOC_R_H

#include <R.h>
#include <Rinternals.h>

// Declaration specification for exported functions
#if defined(_WIN32) || defined(__WINDOWS__)
   #define REXPORTS extern "C" __declspec(dllexport)
#else
   #define REXPORTS extern "C"
#endif


/***********************************************************************
         Class StochasticLib1
***********************************************************************/

class StocRBase {
   // This class is used as base class for the random variate generating 
   // classes when used for the R-language interface
   // Encapsulates the random number generator in R.DLL.
public:
   StocRBase(int32 seed) {}                         // Constructor
   static void InitRan() {                          // Call this before first random number
      GetRNGstate();}                               // From R.DLL
   static void EndRan() {                           // Call this after last random number
      PutRNGstate();}                               // From R.DLL
   double Random() {                                // output random float number in the interval 0 <= x < 1
      return unif_rand();}                          // From R.DLL
   double Normal(double m, double s) {              // normal distribution
      return norm_rand()*s + m;}                    // From R.DLL
   int32 Hypergeometric(int32 n, int32 m, int32 N); // hypergeometric distribution (stocR.cpp)
protected:
   int32 HypInversionMod (int32 n, int32 M, int32 N);  // hypergeometric by inversion searching from mode
   int32 HypRatioOfUnifoms (int32 n, int32 M, int32 N);// hypergeometric by ratio of uniforms method
   static double fc_lnpk(int32 k, int32 N_Mn, int32 M, int32 n); // used by Hypergeometric
};

#endif
