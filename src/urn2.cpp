/*************************** urn2.cpp **********************************
* Author:        Agner Fog
* Date created:  2006
* Last modified: 2011-08-05
* Project:       BiasedUrn
* Source URL:    www.agner.org/random
*
* Description:
* R interface to multivariate noncentral hypergeometric distributions
*
* Copyright 2006-2011 by Agner Fog. 
* GNU General Public License http://www.gnu.org/licenses/gpl.html
*****************************************************************************/

#include <R.h>
#include <Rinternals.h>
#include "stocc.h"


/******************************************************************************
      dMFNCHypergeo
      Mass function for 
      Multivariate Fisher's NonCentral Hypergeometric distribution
******************************************************************************/
REXPORTS SEXP dMFNCHypergeo(
SEXP rx,         // Number of balls drawn of each color, vector or matrix
SEXP rm,         // Number of balls of each color in urn, vector
SEXP rn,         // Number of balls drawn from urn, scalar
SEXP rodds,      // Odds for each color, vector
SEXP rprecision  // Precision of calculation, scalar
) {

   // Check number of colors
   int colors = LENGTH(rm);
   if (colors < 1) error ("Number of colors too small");
   if (colors > MAXCOLORS) {
      error ("Number of colors (%i) exceeds maximum (%i).\n"
      "You may recompile the BiasedUrn package with a bigger value of MAXCOLORS in the file Makevars.",
      colors, MAXCOLORS);
   }   
   if (LENGTH(rn) != 1 || LENGTH(rprecision) != 1) error("Parameter n has wrong length");
   int nres;     // Number of results
   if (isMatrix(rx)) {
      nres = ncols(rx);
      if (nrows(rx) != colors) error("matrix x must have one row for each color and one column for each sample");
   }
   else {
      nres = 1;
      if (LENGTH(rx) != colors) error("Length of vectors x, m, and odds must be the same");
   }

   // Get parameter values
   int32 * px    =  INTEGER(rx);
   int32 * pm    =  INTEGER(rm);
   int     n     = *INTEGER(rn);
   double *podds =  REAL(rodds);
   double  prec  = *REAL(rprecision);
   int     N;                          // Total number of balls
   int     Nu;                         // Total number of balls with nonzero odds
   int     i, j;                       // Loop counter
   int     xsum;                       // Column sum of x = n

   // Check if odds = 1
   double OddsOne[MAXCOLORS];          // Used if odds = 1
   if (LENGTH(rodds) == 1 && *podds == 1.) {
      // Odds = scalar 1. Set to vector of all 1's
      for (i = 0; i < colors; i++) OddsOne[i] = 1.;
      podds = OddsOne;
   }
   else {
      if (LENGTH(rodds) != colors) error("Length of odds vector must match length of m vector");
   }

   // Get N = sum(m) and check validity of m and odds
   for (N = Nu = i = 0; i < colors; i++) {
      int32 m = pm[i];
      if (m < 0) error("m[%i] < 0", i+1);
      N += m;
      if (podds[i]) Nu += m;
      if ((unsigned int)N > 2000000000) error ("Integer overflow");
      if (!R_FINITE(podds[i]) || podds[i] < 0) error("Invalid value for odds[%i]", i+1);
   }

   // Check validity of scalar parameters
   if (n < 0)  error("Negative parameter n");
   if (n > N)  error ("n > sum(m): Taking more items than there are");
   if (n > Nu) error ("Not enough items with nonzero odds");
   if (!R_FINITE(prec) || prec < 0 || prec > 1) prec = 1E-7;

   // Allocate result vector
   SEXP result;  double * presult;
   PROTECT(result = allocVector(REALSXP, nres));
   presult = REAL(result);

   // Make object for calculating probabilities
   CMultiFishersNCHypergeometric mfnc(n, pm, podds, colors, prec);

   // Loop over x inputs
   for (i = 0; i < nres; i++) {
      // Calculate x sum and check each x
      for (xsum = j = 0; j < colors; j++) {
         xsum += px[j];
         /* Include this if you want error messages for x < 0 and x > m
         if (px[j] > pm[j]) {
            // Error
            if (nres == 1) error("x[%i] = %i is bigger than m[%i] = %i", j+1, px[j], j+1, pm[j]);
            else error("x[%i,%i] = %i is bigger than m[%i] = %i", j+1, i+1, px[j], j+1, pm[j]);
         }
         else if (px[j] < 0) {
            if (nres == 1) error("x[%i] = %i is negative", j+1, px[j]);
            else error("x[%i,%i] = %i is negative", j+1, i+1, px[j]);
         }
         */
      }
      // Check x sum
      if (xsum != n) {
         // Error
         if (nres == 1) error("sum(x) = %i must be equal to n = %i", xsum, n);
         else error("sum(x[,%i]) = %i must be equal to n = %i", i+1, xsum, n);
      }

      // Calculate probability
      presult[i] = mfnc.probability(px);         // Probability

      // Get next column
      px += colors;
   }

   // Return result
   UNPROTECT(1);
   return(result);
}


/******************************************************************************
      dMWNCHypergeo
      Mass function for 
      Multivariate Wallenius' NonCentral Hypergeometric distribution
******************************************************************************/
REXPORTS SEXP dMWNCHypergeo(
SEXP rx,         // Number of balls drawn of each color, vector or matrix
SEXP rm,         // Number of balls of each color in urn, vector
SEXP rn,         // Number of balls drawn from urn, scalar
SEXP rodds,      // Odds for each color, vector
SEXP rprecision  // Precision of calculation, scalar
) {

   // Check number of colors
   int colors = LENGTH(rm);
   if (colors < 1) error ("Number of colors too small");
   if (colors > MAXCOLORS) {
      error ("Number of colors (%i) exceeds maximum (%i).\n"
      "You may recompile the BiasedUrn package with a bigger value of MAXCOLORS in the file Makevars.",
      colors, MAXCOLORS);
   }   
   if (LENGTH(rn) != 1 || LENGTH(rprecision) != 1) error("Parameter n has wrong length");
   int nres;     // Number of results
   if (isMatrix(rx)) {
      nres = ncols(rx);
      if (nrows(rx) != colors) error("matrix x must have one row for each color and one column for each sample");
   }
   else {
      nres = 1;
      if (LENGTH(rx) != colors) error("Length of vectors x, m, and odds must be the same");
   }

   // Get parameter values
   int32 * px    =  INTEGER(rx);
   int32 * pm    =  INTEGER(rm);
   int     n     = *INTEGER(rn);
   double *podds =  REAL(rodds);
   double  prec  = *REAL(rprecision);
   int     N;                          // Total number of balls
   int     Nu;                         // Total number of balls with nonzero odds
   int     i, j;                       // Loop counter
   int     xsum;                       // Column sum of x = n

   // Check if odds = 1
   double OddsOne[MAXCOLORS];          // Used if odds = 1
   if (LENGTH(rodds) == 1 && *podds == 1.) {
      // Odds = scalar 1. Set to vector of all 1's
      for (i = 0; i < colors; i++) OddsOne[i] = 1.;
      podds = OddsOne;
   }
   else {
      if (LENGTH(rodds) != colors) error("Length of odds vector must match length of m vector");
   }

   // Get N = sum(m) and check validity of m and odds
   for (N = Nu = i = 0; i < colors; i++) {
      int32 m = pm[i];
      if (m < 0) error("m[%i] < 0", i+1);
      N += m;
      if (podds[i]) Nu += m;
      if ((unsigned int)N > 2000000000) error ("Integer overflow");
      if (!R_FINITE(podds[i]) || podds[i] < 0) error("Invalid value for odds[%i]", i+1);
   }

   // Check validity of scalar parameters
   if (n < 0)  error("Negative parameter n");
   if (n > N)  error ("n > sum(m): Taking more items than there are");
   if (n > Nu) error ("Not enough items with nonzero odds");
   if (!R_FINITE(prec) || prec < 0 || prec > 1) prec = 1E-7;

   // Allocate result vector
   SEXP result;  double * presult;
   PROTECT(result = allocVector(REALSXP, nres));
   presult = REAL(result);

   // Make object for calculating probabilities
   CMultiWalleniusNCHypergeometric mwnc(n, pm, podds, colors, prec);

   // Loop over x inputs
   for (i = 0; i < nres; i++) {
      // Calculate x sum and check each x
      for (xsum = j = 0; j < colors; j++) {
         xsum += px[j];
         /* Include this if you want error messages for x > m and x < 0
         if (px[j] > pm[j]) {
            // Error
            if (nres == 1) error("x[%i] = %i is bigger than m[%i] = %i", j+1, px[j], j+1, pm[j]);
            else error("x[%i,%i] = %i is bigger than m[%i] = %i", j+1, i+1, px[j], j+1, pm[j]);
         }
         else if (px[j] < 0) {
            if (nres == 1) error("x[%i] = %i is negative", j+1, px[j]);
            else error("x[%i,%i] = %i is negative", j+1, i+1, px[j]);
         }
         */
      }
      // Check x sum
      if (xsum != n) {
         // Error
         if (nres == 1) error("sum(x) = %i must be equal to n = %i", xsum, n);
         else error("sum(x[,%i]) = %i must be equal to n = %i", i+1, xsum, n);
      }

      // Calculate probability
      presult[i] = mwnc.probability(px);         // Probability

      // Get next column
      px += colors;
   }

   // Return result
   UNPROTECT(1);
   return(result);
}


/******************************************************************************
      rMFNCHypergeo
      Random variate generation function for
      Multivariate Fisher's NonCentral Hypergeometric distribution
******************************************************************************/
REXPORTS SEXP rMFNCHypergeo(
SEXP rnran,      // Number of random variates desired, scalar
SEXP rm,         // Number of balls of each color in urn, vector
SEXP rn,         // Number of balls drawn from urn, scalar
SEXP rodds,      // Odds for each color, vector
SEXP rprecision  // Precision of calculation, scalar
) {

   // Check number of colors
   int colors = LENGTH(rm);
   if (colors < 1) error ("Number of colors too small");
   if (colors > MAXCOLORS) {
      error ("Number of colors (%i) exceeds maximum (%i).\n"
      "You may recompile the BiasedUrn package with a bigger value of MAXCOLORS in the file Makevars.",
      colors, MAXCOLORS);
   }   
   if (LENGTH(rn) != 1) error("Parameter n has wrong length");
   if (LENGTH(rprecision) != 1) error("Parameter precision has wrong length");

   // Get parameter values
   int     nran = *INTEGER(rnran);  if (LENGTH(rnran) > 1) nran = LENGTH(rnran);
   int32 * pm    =  INTEGER(rm);
   int     n     = *INTEGER(rn);
   double *podds =  REAL(rodds);
   double  prec  = *REAL(rprecision);
   int     i;                          // Loop counter
   int     N;                          // Total number of balls
   int     Nu;                         // Total number of balls with nonzero odds

   // Check validity of scalar parameters
   if (n < 0)  error("Negative parameter n");
   if (nran <= 0)  error("Parameter nran must be positive");
   if (!R_FINITE(prec) || prec < 0 || prec > 1) prec = 1E-7;

   // Check if odds = 1
   double OddsOne[MAXCOLORS];          // Used if odds = 1
   if (LENGTH(rodds) == 1 && *podds == 1.) {
      // Odds = scalar 1. Set to vector of all 1's
      for (i = 0; i < colors; i++) OddsOne[i] = 1.;
      podds = OddsOne;
   }
   else {
      if (LENGTH(rodds) != colors) error("Length of odds vector must match length of m vector");
   }

   // Get N = sum(m) and check validity of m and odds
   for (N = Nu = i = 0; i < colors; i++) {
      int32 m = pm[i];
      if (m < 0) error("m[%i] < 0", i+1);
      N += m;
      if (podds[i]) Nu += m;
      if ((unsigned int)N > 2000000000) error ("Integer overflow");
      if (!R_FINITE(podds[i]) || podds[i] < 0) error("Invalid value for odds[%i]", i+1);
   }
   if (n > N)  error ("n > sum(m): Taking more items than there are");
   if (n > Nu) error ("Not enough items with nonzero odds");

   // Allocate result vector
   SEXP result;  int * presult;
   if (nran <= 1) { // One result. Make vector
      PROTECT(result = allocVector(INTSXP, colors));
   }
   else {           // Multiple results. Make matrix
      PROTECT(result = allocMatrix(INTSXP, colors, nran));
   }  
   
   presult = INTEGER(result);

   // Make object for generating variates
   StochasticLib3 sto(0);              // Seed is not used
   sto.SetAccuracy(prec);              // Set precision
   sto.InitRan();                      // Initialize RNG in R.dll

   // Generate variates one by one
   for (i = 0; i < nran; i++) {
      sto.MultiFishersNCHyp(presult, pm, podds, n, colors); // Generate variate
      presult += colors;               // Point to next column in matrix
   }

   sto.EndRan();                       // Return RNG state to R.dll

   // Return result
   UNPROTECT(1);
   return(result);
}


/******************************************************************************
      rMWNCHypergeo
      Random variate generation function for
      Multivariate Wallenius' NonCentral Hypergeometric distribution
******************************************************************************/
REXPORTS SEXP rMWNCHypergeo(
SEXP rnran,      // Number of random variates desired, scalar
SEXP rm,         // Number of balls of each color in urn, vector
SEXP rn,         // Number of balls drawn from urn, scalar
SEXP rodds,      // Odds for each color, vector
SEXP rprecision  // Precision of calculation, scalar
) {

   // Check number of colors
   int colors = LENGTH(rm);
   if (colors < 1) error ("Number of colors too small");
   if (colors > MAXCOLORS) {
      error ("Number of colors (%i) exceeds maximum (%i).\n"
      "You may recompile the BiasedUrn package with a bigger value of MAXCOLORS in the file Makevars.",
      colors, MAXCOLORS);
   }   
   if (LENGTH(rn) != 1) error("Parameter n has wrong length");
   if (LENGTH(rprecision) != 1) error("Parameter precision has wrong length");

   // Get parameter values
   int     nran  = *INTEGER(rnran);  if (LENGTH(rnran) > 1) nran = LENGTH(rnran);
   int32 * pm    =  INTEGER(rm);
   int     n     = *INTEGER(rn);
   double *podds =  REAL(rodds);
   double  prec  = *REAL(rprecision);
   int     i;                          // Loop counter
   int     N;                          // Total number of balls
   int     Nu;                         // Total number of balls with nonzero odds

   // Check validity of scalar parameters
   if (n < 0)  error("Negative parameter n");
   if (nran <= 0)  error("Parameter nran must be positive");
   if (!R_FINITE(prec) || prec < 0 || prec > 1) prec = 1E-7;

   // Check if odds = 1
   double OddsOne[MAXCOLORS];          // Used if odds = 1
   if (LENGTH(rodds) == 1 && *podds == 1.) {
      // Odds = scalar 1. Set to vector of all 1's
      for (i = 0; i < colors; i++) OddsOne[i] = 1.;
      podds = OddsOne;
   }
   else {
      if (LENGTH(rodds) != colors) error("Length of odds vector must match length of m vector");
   }

   // Get N = sum(m) and check validity of m and odds
   for (N = Nu = i = 0; i < colors; i++) {
      int32 m = pm[i];
      if (m < 0) error("m[%i] < 0", i+1);
      N += m;
      if (podds[i]) Nu += m;
      if ((unsigned int)N > 2000000000) error ("Integer overflow");
      if (!R_FINITE(podds[i]) || podds[i] < 0) error("Invalid value for odds[%i]", i+1);
   }
   if (n > N)  error ("n > sum(m): Taking more items than there are");
   if (n > Nu) error ("Not enough items with nonzero odds");

   // Allocate result vector
   SEXP result;  int * presult;
   if (nran <= 1) { // One result. Make vector
      PROTECT(result = allocVector(INTSXP, colors));
   }
   else {           // Multiple results. Make matrix
      PROTECT(result = allocMatrix(INTSXP, colors, nran));
   }  
   
   presult = INTEGER(result);

   // Make object for generating variates
   StochasticLib3 sto(0);              // Seed is not used
   sto.SetAccuracy(prec);              // Set precision
   sto.InitRan();                      // Initialize RNG in R.dll

   // Generate variates one by one
   for (i = 0; i < nran; i++) {
      sto.MultiWalleniusNCHyp(presult, pm, podds, n, colors); // Generate variate
      presult += colors;               // Point to next column in matrix
   }

   sto.EndRan();                       // Return RNG state to R.dll

   // Return result
   UNPROTECT(1);
   return(result);
}


/******************************************************************************
      momentsMFNCHypergeo
      Calculates the mean and variance of the
      Multivariate Fisher's NonCentral Hypergeometric distribution
******************************************************************************/
REXPORTS SEXP momentsMFNCHypergeo(
SEXP rm,         // Number of balls of each color in urn, vector
SEXP rn,         // Number of balls drawn from urn, scalar
SEXP rodds,      // Odds for each color, vector
SEXP rprecision  // Precision of calculation, scalar
) {

   // Check number of colors
   int colors = LENGTH(rm);
   if (colors < 1) error ("Number of colors too small");
   if (colors > MAXCOLORS) {
      error ("Number of colors (%i) exceeds maximum (%i).\n"
      "You may recompile the BiasedUrn package with a bigger value of MAXCOLORS in the file Makevars.",
      colors, MAXCOLORS);
   }   
   if (LENGTH(rn) != 1) error("Parameter n has wrong length");
   if (LENGTH(rprecision) != 1) error("Parameter precision has wrong length");

   // Get parameter values
   int32 * pm    =  INTEGER(rm);
   int     n     = *INTEGER(rn);
   double *podds =  REAL(rodds);
   double  prec  = *REAL(rprecision);
   int     i;                          // Loop counter
   int     N;                          // Total number of balls
   int     Nu;                         // Total number of balls with nonzero odds

   // Check validity of scalar parameters
   if (n < 0)  error("Negative parameter n");
   if (!R_FINITE(prec) || prec < 0) prec = 1;

   // Check if odds = 1
   double OddsOne[MAXCOLORS];          // Used if odds = 1
   if (LENGTH(rodds) == 1 && *podds == 1.) {
      // Odds = scalar 1. Set to vector of all 1's
      for (i = 0; i < colors; i++) OddsOne[i] = 1.;
      podds = OddsOne;
   }
   else {
      if (LENGTH(rodds) != colors) error("Length of odds vector must match length of m vector");
   }

   // Get N = sum(m) and check validity of m and odds
   for (N = Nu = i = 0; i < colors; i++) {
      int32 m = pm[i];
      if (m < 0) error("m[%i] < 0", i+1);
      N += m;
      if (podds[i]) Nu += m;
      if ((unsigned int)N > 2000000000) error ("Integer overflow");
      if (!R_FINITE(podds[i]) || podds[i] < 0) error("Invalid value for odds[%i]", i+1);
   }
   if (n > N)  error ("n > sum(m): Taking more items than there are");
   if (n > Nu) error ("Not enough items with nonzero odds");

   // Allocate result vector
   SEXP result;  double * presult;
   PROTECT(result = allocMatrix(REALSXP, colors, 2));
   
   presult = REAL(result);

   // Make object for calculating mean and variance
   CMultiFishersNCHypergeometric mfnc(n, pm, podds, colors, prec);

   if (prec >= 0.1) {
      // use approximate calculation methods
      mfnc.variance(presult + colors, presult);
   }
   else {
      // use exact calculation
      mfnc.moments(presult, presult + colors);
   }

   // Return result
   UNPROTECT(1);
   return(result);
}


/******************************************************************************
      momentsMWNCHypergeo
      Calculates the mean and variance of the
      Multivariate Wallenius' NonCentral Hypergeometric distribution
******************************************************************************/
REXPORTS SEXP momentsMWNCHypergeo(
SEXP rm,         // Number of balls of each color in urn, vector
SEXP rn,         // Number of balls drawn from urn, scalar
SEXP rodds,      // Odds for each color, vector
SEXP rprecision  // Precision of calculation, scalar
) {

   // Check number of colors
   int colors = LENGTH(rm);
   if (colors < 1) error ("Number of colors too small");
   if (colors > MAXCOLORS) {
      error ("Number of colors (%i) exceeds maximum (%i).\n"
      "You may recompile the BiasedUrn package with a bigger value of MAXCOLORS in the file Makevars.",
      colors, MAXCOLORS);
   }   
   if (LENGTH(rn) != 1) error("Parameter n has wrong length");
   if (LENGTH(rprecision) != 1) error("Parameter precision has wrong length");

   // Get parameter values
   int32 * pm    =  INTEGER(rm);
   int     n     = *INTEGER(rn);
   double *podds =  REAL(rodds);
   double  prec  = *REAL(rprecision);
   int     i;                          // Loop counter
   int     N;                          // Total number of balls
   int     Nu;                         // Total number of balls with nonzero odds

   // Check validity of scalar parameters
   if (n < 0)  error("Negative parameter n");
   if (!R_FINITE(prec) || prec < 0) prec = 1;

   // Check if odds = 1
   double OddsOne[MAXCOLORS];          // Used if odds = 1
   if (LENGTH(rodds) == 1 && *podds == 1.) {
      // Odds = scalar 1. Set to vector of all 1's
      for (i = 0; i < colors; i++) OddsOne[i] = 1.;
      podds = OddsOne;
   }
   else {
      if (LENGTH(rodds) != colors) error("Length of odds vector must match length of m vector");
   }

   // Get N = sum(m) and check validity of m and odds
   for (N = Nu = i = 0; i < colors; i++) {
      int32 m = pm[i];
      if (m < 0) error("m[%i] < 0", i+1);
      N += m;
      if (podds[i]) Nu += m;
      if ((unsigned int)N > 2000000000) error ("Integer overflow");
      if (!R_FINITE(podds[i]) || podds[i] < 0) error("Invalid value for odds[%i]", i+1);
   }
   if (n > N)  error ("n > sum(m): Taking more items than there are");
   if (n > Nu) error ("Not enough items with nonzero odds");

   // Allocate result vector
   SEXP result;  double * presult;
   PROTECT(result = allocMatrix(REALSXP, colors, 2));
   
   presult = REAL(result);

   // Make object for calculating mean and variance
   CMultiWalleniusNCHypergeometricMoments mwnc(n, pm, podds, colors, prec);

   if (prec >= 0.1) {
      // use approximate calculation methods
      mwnc.variance(presult + colors, presult);
   }
   else {
      // use exact calculation
      mwnc.moments(presult, presult + colors);
   }

   // Return result
   UNPROTECT(1);
   return(result);
}


/******************************************************************************
      oddsMFNCHypergeo
      Estimate odds ratio from mean for the
      Multivariate Fisher's NonCentral Hypergeometric distribution
******************************************************************************/
// Uses the multivariate extension of Cornfield's approximation.
// Precision is ignored
REXPORTS SEXP oddsMFNCHypergeo(
SEXP rmu,        // Number of balls drawn of each color, vector or matrix
SEXP rm,         // Number of balls of each color in urn, vector
SEXP rn,         // Number of balls drawn from urn, scalar
SEXP rprecision  // Precision of calculation, scalar
) {
   // Check number of colors
   int colors = LENGTH(rm);
   if (colors < 1) error ("Number of colors too small");
   if (colors > MAXCOLORS) {
      error ("Number of colors (%i) exceeds maximum (%i).\n"
      "You may recompile the BiasedUrn package with a bigger value of MAXCOLORS in the file Makevars.",
      colors, MAXCOLORS);
   }   
   int nres;     // Number of results
   if (isMatrix(rmu)) {
      nres = ncols(rmu);
      if (nrows(rmu) != colors) error("matrix mu must have one row for each color and one column for each sample");
   }
   else {
      nres = 1;
      if (LENGTH(rmu) != colors) error("Length of vectors mu and m must be the same");
   }

   // Get parameter values
   double *pmu   =  REAL(rmu);
   int32 * pm    =  INTEGER(rm);
   int     n     = *INTEGER(rn);
   double  prec  = *REAL(rprecision);
   int     N;                          // Total number of balls
   int     i, j;                       // Loop counter
   int     x1, x2;                     // x limits
   int     c0;                         // Reference color
   double  xd0, xd1, xd2;              // Used for searching for reference color
   double  mu;                         // Mean
   double  sum_mu = 0.;                // Sum of means
   int     err = 0;                    // Warning and error messages

   // Get N = sum(m) and check validity of m and odds
   for (N = i = 0; i < colors; i++) {
      int32 m = pm[i];
      if (m < 0) error("m[%i] < 0", i+1);
      N += m;
      if ((unsigned int)N > 2000000000) error ("Integer overflow");
      sum_mu += pmu[i];
   }
   if (n > 0 && fabs(sum_mu-n)/n > 0.1) {
      err |= 0x100;                    // sum of means should be equal to n
   }

   // Check validity of scalar parameters
   if (n < 0)  error("Negative parameter n");
   if (n > N)  error ("n > sum(m): Taking more items than there are");
   if (!R_FINITE(prec) || prec < 0 || prec > 1) prec = 0.1;
   if (prec < 0.05) warning ("Cannot obtain high precision");

   // Allocate result vector
   SEXP result;  double * presult;
   if (nres == 1) {
      PROTECT(result = allocVector(REALSXP, colors));
   }
   else {
      PROTECT(result = allocMatrix(REALSXP, colors, nres));
   }
   presult = REAL(result);

   // Loop over x inputs
   for (i = 0; i < nres; i++) {

      // Find the color with the highest variance to use as reference
      for (xd0 = 0., j = c0 = 0; j < colors; j++) {
         // Get limits for x[j]
         x1 = pm[j] + n - N;  if (x1 < 0) x1 = 0;
         x2 = n;  if (x2 > pm[j]) x2 = pm[j];
         // Find max distance of mu from limits
         xd1 = pmu[j] - x1;  xd2 = x2 - pmu[j];
         if (xd1 > xd2) xd1 = xd2;
         if (xd1 > xd0) {xd0 = xd1;  c0 = j;}
      }
      if (xd0 == 0.) {
         // All odds are indetermined
         err |= 0x10;
         for (j = 0; j < colors; j++) presult[j] = R_NaN;
      }
      else {
         // Use color c0 as reference
         presult[c0] = 1.;

         // Get odds for all colors except c0
         for (j = 0; j < colors; j++) {
            if (j != c0) {

               // Get limits for x[j]
               x1 = pm[j] + n - N;  if (x1 < 0) x1 = 0;
               x2 = n;  if (x2 > pm[j]) x2 = pm[j];
               mu = pmu[j];

               // Check limits
               if (x1 == x2) {
                  presult[j] = R_NaN;  err |= 1;         // Indetermined
                  continue;
               }
               if (mu <= double(x1)) {
                  if (mu == double(x1)) {
                     presult[j] = 0.;  err |= 2;         // Zero
                     continue;
                  }
                  presult[j] = R_NaN;  err |= 8;         // Out of range
                  continue;
               }
               if (mu >= double(x2)) {
                  if (mu == double(x2)) {
                     presult[j] = R_PosInf;  err |= 4;   // Infinite
                     continue;
                  }
                  presult[j] = R_NaN;  err |= 8;         // Out of range  
                  continue;
               }

               // Calculate odds relative to c0
               presult[j] = pmu[j] * (pm[c0] - pmu[c0]) / (pmu[c0] * (pm[j] - pmu[j]));
            }
         }
      }
      presult += colors;  pmu += colors;
   }

   // Check for errors
   if (err & 0x10) warning("All odds are indetermined");
   else if (err & 8) error("mu out of range");
   else if (err & 1) warning("odds is indetermined");
   else {
      if (err & 4) warning("odds is infinite");
      if (err & 2) warning("odds is zero with no precision");
   }
   if (err & 0x100) warning("Sum of means should be equal to n");

   // Return result
   UNPROTECT(1);
   return(result);
}


/******************************************************************************
      oddsMWNCHypergeo
      Estimate odds ratio from mean for the
      Multivariate Wallenius' NonCentral Hypergeometric distribution
******************************************************************************/
// Uses Manly's approximation. Precision is ignored
REXPORTS SEXP oddsMWNCHypergeo(
SEXP rmu,        // Number of balls drawn of each color, vector or matrix
SEXP rm,         // Number of balls of each color in urn, vector
SEXP rn,         // Number of balls drawn from urn, scalar
SEXP rprecision // Precision of calculation, scalar
) {
   // Check number of colors
   int colors = LENGTH(rm);
   if (colors < 1) error ("Number of colors too small");
   if (colors > MAXCOLORS) {
      error ("Number of colors (%i) exceeds maximum (%i).\n"
      "You may recompile the BiasedUrn package with a bigger value of MAXCOLORS in the file Makevars.",
      colors, MAXCOLORS);
   }   
   int nres;     // Number of results
   if (isMatrix(rmu)) {
      nres = ncols(rmu);
      if (nrows(rmu) != colors) error("matrix mu must have one row for each color and one column for each sample");
   }
   else {
      nres = 1;
      if (LENGTH(rmu) != colors) error("Length of vectors mu and m must be the same");
   }

   // Get parameter values
   double *pmu   =  REAL(rmu);
   int32 * pm    =  INTEGER(rm);
   int     n     = *INTEGER(rn);
   double  prec  = *REAL(rprecision);
   int     N;                          // Total number of balls
   int     i, j;                       // Loop counter
   int     x1, x2;                     // x limits
   int     c0;                         // Reference color
   double  xd0, xd1, xd2;              // Used for searching for reference color
   double  mu;                         // Mean
   double  sum_mu = 0.;                // Sum of means
   int     err = 0;                    // Warning and error messages

   // Get N = sum(m) and check validity of m and odds
   for (N = i = 0; i < colors; i++) {
      int32 m = pm[i];
      if (m < 0) error("m[%i] < 0", i+1);
      N += m;
      if ((unsigned int)N > 2000000000) error ("Integer overflow");
      sum_mu += pmu[i];
   }
   if (n > 0 && fabs(sum_mu-n)/n > 0.1) {
      err |= 0x100;                    // sum of means should be equal to n
   }

   // Check validity of scalar parameters
   if (n < 0)  error("Negative parameter n");
   if (n > N)  error ("n > sum(m): Taking more items than there are");
   if (!R_FINITE(prec) || prec < 0 || prec > 1) prec = 0.1;
   if (prec < 0.02) warning ("Cannot obtain high precision");

   // Allocate result vector
   SEXP result;  double * presult;
   if (nres == 1) {
      PROTECT(result = allocVector(REALSXP, colors));
   }
   else {
      PROTECT(result = allocMatrix(REALSXP, colors, nres));
   }
   presult = REAL(result);

   // Loop over x inputs
   for (i = 0; i < nres; i++) {

      // Find the color with the highest variance to use as reference
      for (xd0 = 0., j = c0 = 0; j < colors; j++) {
         // Get limits for x[j]
         x1 = pm[j] + n - N;  if (x1 < 0) x1 = 0;
         x2 = n;  if (x2 > pm[j]) x2 = pm[j];
         // Find max distance of mu from limits
         xd1 = pmu[j] - x1;  xd2 = x2 - pmu[j];
         if (xd1 > xd2) xd1 = xd2;
         if (xd1 > xd0) {xd0 = xd1;  c0 = j;}
      }
      if (xd0 == 0.) {
         // All odds are indetermined
         err |= 0x10;
         for (j = 0; j < colors; j++) presult[j] = R_NaN;
      }
      else {
         // Use color c0 as reference
         presult[c0] = 1.;

         // Get odds for all colors except c0
         for (j = 0; j < colors; j++) {
            if (j != c0) {

               // Get limits for x[j]
               x1 = pm[j] + n - N;  if (x1 < 0) x1 = 0;
               x2 = n;  if (x2 > pm[j]) x2 = pm[j];
               mu = pmu[j];

               // Check limits
               if (x1 == x2) {
                  presult[j] = R_NaN;  err |= 1;         // Indetermined
                  continue;
               }
               if (mu <= double(x1)) {
                  if (mu == double(x1)) {
                     presult[j] = 0.;  err |= 2;         // Zero
                     continue;
                  }
                  presult[j] = R_NaN;  err |= 8;         // Out of range
                  continue;
               }
               if (mu >= double(x2)) {
                  if (mu == double(x2)) {
                     presult[j] = R_PosInf;  err |= 4;   // Infinite
                     continue;
                  }
                  presult[j] = R_NaN;  err |= 8;         // Out of range  
                  continue;
               }

               // Calculate odds relative to c0
               presult[j] = log(1. - pmu[j] / pm[j]) / log(1. - pmu[c0] / pm[c0]);
            }
         }
      }
      presult += colors;  pmu += colors;
   }

   // Check for errors
   if (err & 0x10) warning("All odds are indetermined");
   else if (err & 8) error("mu out of range");
   else if (err & 1) warning("odds is indetermined");
   else {
      if (err & 4) warning("odds is infinite");
      if (err & 2) warning("odds is zero with no precision");
   }
   if (err & 0x100) warning("Sum of means should be equal to n");

   // Return result
   UNPROTECT(1);
   return(result);
}


/******************************************************************************
      numMFNCHypergeo
      Estimate number of balls of each color from experimental mean for
      Multivariate Fisher's NonCentral Hypergeometric distribution
******************************************************************************/
// Uses Cornfield's approximation. Precision is ignored.
// Calculation method: Solves the multivariate Cornfield's equation by
// Newton Raphson iteration with r as independent parameter.
REXPORTS SEXP numMFNCHypergeo(
SEXP rmu,        // Observed mean of x1
SEXP rn,         // Number of balls drawn from urn
SEXP rN,         // Number of balls in urn before sampling
SEXP rodds,      // Odds of getting a red ball among one red and one white
SEXP rprecision  // Precision of calculation
) {
   int nres;     // Number of results
   int colors;   // Number of colors

   // Check for vectors
   if (LENGTH(rn)         != 1 
   || LENGTH(rN)          != 1 
   || LENGTH(rprecision)  != 1 
   ) {
      error("Parameter has wrong length");
   }

   // Check mu matrix size
   if (isMatrix(rmu)) {
      nres = ncols(rmu);
      colors = nrows(rmu);
   }
   else {
      nres = 1;
      colors = LENGTH(rmu);
   }

   // Check number of colors
   if (colors < 1) error ("Number of colors too small");
   if (colors > MAXCOLORS) {
      error ("Number of colors (%i) exceeds maximum (%i).\n"
      "You may recompile the BiasedUrn package with a bigger value of MAXCOLORS in the file Makevars.",
      colors, MAXCOLORS);
   }

   // Get parameter values
   double *pmu   =  REAL(rmu);
   int     n     = *INTEGER(rn);
   int     N     = *INTEGER(rN);
   double *podds =  REAL(rodds);
   double  prec  = *REAL(rprecision);

   int     i, j;                       // Loop counter
   int     err, err1 = 0;              // Remember any error
   int     cu = 0;                     // Number of colors with nonzero odds
   double  smu;                        // Sum of means, reciprocal.
   double  mu[MAXCOLORS];              // Normalized means

   // Check if odds = 1
   double OddsOne[MAXCOLORS];          // Used if odds = 1
   if (LENGTH(rodds) == 1 && *podds == 1.) {
      // Odds = scalar 1. Set to vector of all 1's
      for (i = 0; i < colors; i++) OddsOne[i] = 1.;
      podds = OddsOne;
   }
   else {
      if (LENGTH(rodds) != colors) {
         // Size mismatch
         if (isMatrix(rmu)) {      
            error("matrix mu must have one row for each color and one column for each sample");
         }
         else {
            error("Length of vectors mu and odds must be the same");
         }
      }
   }

   // Check validity of parameters
   if (n < 0 || N < 0) error("Negative parameter");
   if ((unsigned int)N > 2000000000) error("Overflow");
   if (n > N) error ("n > N: Taking more items than there are");
   if (!R_FINITE(prec) || prec < 0 || prec > 1) prec = 0.1;
   if (prec < 0.05) warning ("Cannot obtain high precision");

   // Check validity of odds
   for (i = cu = 0; i < colors; i++) {
      if (!R_FINITE(podds[i]) || podds[i] < 0) error("Invalid value for odds[%i]", i+1);
      if (podds[i] > 0) cu++;
   }

   // Allocate result vector
   SEXP result;  double * presult;
   if (nres == 1) {
      PROTECT(result = allocVector(REALSXP, colors));
   }
   else {
      PROTECT(result = allocMatrix(REALSXP, colors, nres));
   }
   presult = REAL(result);

   // Loop for all mu inputs
   for (j = 0; j < nres; j++, presult += colors, pmu += colors) {
      err = 0;

      // Make results NAN in case of error exits below
      for (i = 0; i < colors; i++) presult[i] = R_NaN;  

      // Check limits
      if (n == 0) {
         err1 |= 1;         // Indetermined
         continue;
      }

      // Check sum of mu must equal n
      for (i = 0, smu = 0.; i < colors; i++) smu += pmu[i];
      if (smu <= 0.) {
         err1 |= 0x800;     // Sum of means must be positive
         break;
      }
      if (fabs(smu - n) > 0.02 * n) {
         err |= 0x100; // Warning: sum not approx. equal to n
      }
      smu = n / smu;
      for (i = 0; i < colors; i++) {
         mu[i] = pmu[i] * smu;    // Normalize mu
      }

      // More parameter checks
      if (n == N) {               // Results known exactly
         for (i = 0; i < colors; i++) {
            if (podds[i] == 0 && mu[i] != 0) {
               err1 |= 0x10;        // Out of range
            }
            else {
               presult[i] = mu[i];
            }
         }
         continue;
      }

      // Check odds
      if (cu < colors || colors < 2) {
         for (i = 0; i < colors; i++) {
            if (podds[i] == 0) {
               if (mu[i] != 0) err1 |= 0x10;   // Out of range
               else err1 |= 1;                 // Indetermined
            }
            else {
               if (cu == 1) presult[i] = N;    // Known exactly
            }
         }
         continue;
      }

      // check mu within bounds
      for (i = 0; i < colors; i++) {
         if (mu[i] <= 0.) {
            if (mu[i] == 0.) {
               presult[i] = 0;
               err |= 2;         // Zero
            }
            else {
               err |= 8;         // Out of range
            }
         }
         if (mu[i] >= double(n)) {
            if (mu[i] == double(n)) {
               presult[i] = N;
               err |= 4;
            }
            else {
               err |= 8;         // Out of range
            }
         }
      }

      if (err & 0x18) {
         // Results invalid
         err1 |= err;
         break;
      }

      // Calculate m[]
      double z;              // Newton Raphson function value
      double zd;             // Newton Raphson derivative of z
      double r, lastr;       // Independent parameter in Newton Raphson iteration
      int niter = 0;         // Number of iterations

      // Initial guess
      r = 1.;

      // Newton Raphson iteration
      do {
         lastr = r;
         // Calculate z and zd
         z = zd = 0.;
         for (i = 0; i < colors; i++) {
            z  += mu[i] * (1. + 1./(r*podds[i]));
            zd -= mu[i] / (podds[i]*r*r);
         }
         r -= (z - N) / zd;
         if (r <= 0.) {
            // r must be positive. Get r within range
            if (r < -lastr) {
               r = lastr * 0.125;
            }
            else {
               r = lastr * 0.5;
            }
         }
         if (++niter > 200) error ("Convergence problem");

      } while (fabs(r-lastr) > r * 1E-8);

      // Get results from r
      for (i = 0; i < colors; i++) {
         presult[i] = mu[i] * (r*podds[i] + 1.) / (r*podds[i]);
      }
      err1 |= err;
   }

   // Check for errors
   if (err1 & 0x808) error("Mean is out of range");
   else {
      if (err1 & 0x010) warning("Zero odds conflicts with nonzero mean");
      if (err1 & 1) warning("Number of items is indetermined");
      if (err1 & 0x100) warning("Sum of means is not equal to n");
   }

   // Return result
   UNPROTECT(1);
   return(result);
}


/******************************************************************************
      numMWNCHypergeo
      Estimate number of balls of each color from experimental mean for
      Multivariate Wallenius' NonCentral Hypergeometric distribution
******************************************************************************/
// Uses Manly's approximation. Precision is ignored.
// Calculation method: Solves Manly's equation by
// Newton Raphson iteration with theta as independent parameter.
REXPORTS SEXP numMWNCHypergeo(
SEXP rmu,        // Observed mean of x1
SEXP rn,         // Number of balls drawn from urn
SEXP rN,         // Number of balls in urn before sampling
SEXP rodds,      // Odds of getting a red ball among one red and one white
SEXP rprecision  // Precision of calculation
) {
   int nres;     // Number of results
   int colors;   // Number of colors

   // Check for vectors
   if (LENGTH(rn)         != 1 
   || LENGTH(rN)          != 1 
   || LENGTH(rprecision)  != 1 
   ) {
      error("Parameter has wrong length");
   }

   // Check mu matrix size
   if (isMatrix(rmu)) {
      nres = ncols(rmu);
      colors = nrows(rmu);
   }
   else {
      nres = 1;
      colors = LENGTH(rmu);
   }

   // Check number of colors
   if (colors < 1) error ("Number of colors too small");
   if (colors > MAXCOLORS) {
      error ("Number of colors (%i) exceeds maximum (%i).\n"
      "You may recompile the BiasedUrn package with a bigger value of MAXCOLORS in the file Makevars.",
      colors, MAXCOLORS);
   }

   // Get parameter values
   double *pmu   =  REAL(rmu);
   int     n     = *INTEGER(rn);
   int     N     = *INTEGER(rN);
   double *podds =  REAL(rodds);
   double  prec  = *REAL(rprecision);

   int     i, j;                       // Loop counter
   int     err, err1 = 0;              // Remember any error
   int     cu = 0;                     // Number of colors with nonzero odds
   double  smu;                        // Sum of means, reciprocal.
   double  mu[MAXCOLORS];              // Normalized means

   // Check if odds = 1
   double OddsOne[MAXCOLORS];          // Used if odds = 1
   if (LENGTH(rodds) == 1 && *podds == 1.) {
      // Odds = scalar 1. Set to vector of all 1's
      for (i = 0; i < colors; i++) OddsOne[i] = 1.;
      podds = OddsOne;
   }
   else {
      if (LENGTH(rodds) != colors) {
         // Size mismatch
         if (isMatrix(rmu)) {      
            error("matrix mu must have one row for each color and one column for each sample");
         }
         else {
            error("Length of vectors mu and odds must be the same");
         }
      }
   }

   // Check validity of parameters
   if (n < 0 || N < 0) error("Negative parameter");
   if ((unsigned int)N > 2000000000) error("Overflow");
   if (n > N) error ("n > N: Taking more items than there are");
   if (!R_FINITE(prec) || prec < 0 || prec > 1) prec = 0.1;
   if (prec < 0.02) warning ("Cannot obtain high precision");

   // Check validity of odds
   for (i = cu = 0; i < colors; i++) {
      if (!R_FINITE(podds[i]) || podds[i] < 0) error("Invalid value for odds[%i]", i+1);
      if (podds[i] > 0) cu++;
   }

   // Allocate result vector
   SEXP result;  double * presult;
   if (nres == 1) {
      PROTECT(result = allocVector(REALSXP, colors));
   }
   else {
      PROTECT(result = allocMatrix(REALSXP, colors, nres));
   }
   presult = REAL(result);

   // Loop for all mu inputs
   for (j = 0; j < nres; j++, presult += colors, pmu += colors) {
      err = 0;

      // Make results NAN in case of error exits below
      for (i = 0; i < colors; i++) presult[i] = R_NaN;  

      // Check limits
      if (n == 0) {
         err1 |= 1;         // Indetermined
         continue;
      }

      // Check sum of mu must equal n
      for (i = 0, smu = 0.; i < colors; i++) smu += pmu[i];
      if (smu <= 0.) {
         err1 |= 0x800;     // Sum of means must be positive
         break;
      }
      if (fabs(smu - n) > 0.02 * n) {
         err |= 0x100; // Warning: sum not approx. equal to n
      }
      smu = n / smu;
      for (i = 0; i < colors; i++) {
         mu[i] = pmu[i] * smu;    // Normalize mu
      }

      // More parameter checks
      if (n == N) {               // Results known exactly
         for (i = 0; i < colors; i++) {
            if (podds[i] == 0 && mu[i] != 0) {
               err1 |= 0x10;        // Out of range
            }
            else {
               presult[i] = mu[i];
            }
         }
         continue;
      }

      // Check odds
      if (cu < colors || colors < 2) {
         for (i = 0; i < colors; i++) {
            if (podds[i] == 0) {
               if (mu[i] != 0) err1 |= 0x10;   // Out of range
               else err1 |= 1;                 // Indetermined
            }
            else {
               if (cu == 1) presult[i] = N;    // Known exactly
            }
         }
         continue;
      }

      // check mu within bounds
      for (i = 0; i < colors; i++) {
         if (mu[i] <= 0.) {
            if (mu[i] == 0.) {
               presult[i] = 0;
               err |= 2;         // Zero
            }
            else {
               err |= 8;         // Out of range
            }
         }
         if (mu[i] >= double(n)) {
            if (mu[i] == double(n)) {
               presult[i] = N;
               err |= 4;
            }
            else {
               err |= 8;         // Out of range
            }
         }
      }

      if (err & 0x18) {
         // Results invalid
         err1 |= err;
         break;
      }

      // Calculate m[]
      double z;              // Newton Raphson function value
      double zd;             // Newton Raphson derivative of z
      double t, lastt;       // Independent parameter in Newton Raphson iteration
      double eot;            // exp(odds[i]*t)
      double eot1 = 1.;      // 1 - exp(odds[i]*t)
      int niter = 0;         // Number of iterations

      // Initial guess
      t = lastt = -1.;

      // Newton Raphson iteration
      do {
         // Calculate z and zd
         AGAIN:
         z = zd = 0.;
         for (i = 0; i < colors; i++) {
            eot = exp(podds[i]*t);
            eot1 = 1. - eot;
            if (eot1 <= 0. || eot <= 0.) {
               // Out of range
               lastt = t;
               t = 0.125 * lastt;
               goto AGAIN;
            }
            z  += mu[i] / eot1;
            zd += mu[i] * podds[i] * eot / (eot1*eot1);
         }
         lastt = t;
         t -= (z - N) / zd;
         if (t >= 0.) {
            // t must be negative. Get t within range
            if (t > -lastt) {
               t = lastt * 0.125;
            }
            else {
               t = lastt * 0.5;
            }
         }
         if (++niter > 200) error ("Convergence problem");

      } while (fabs(t-lastt) > -t * 1E-8);

      // Get results from t
      for (i = 0; i < colors; i++) {
         presult[i] = mu[i] / (1. - exp(podds[i]*t));
      }
      err1 |= err;
   }

   // Check for errors
   if (err1 & 0x808) error("Mean is out of range");
   else {
      if (err1 & 0x010) warning("Zero odds conflicts with nonzero mean");
      if (err1 & 1) warning("Number of items is indetermined");
      if (err1 & 0x100) warning("Sum of means is not equal to n");
   }

   // Return result
   UNPROTECT(1);
   return(result);
}
