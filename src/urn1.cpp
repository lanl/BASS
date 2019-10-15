/*************************** urn1.cpp **********************************
* Author:        Agner Fog
* Date created:  2006
* Last modified: 2011-08-05
* Project:       BiasedUrn
* Source URL:    www.agner.org/random
*
* Description:
* R interface to univariate noncentral hypergeometric distributions
*
* Copyright 2006-2011 by Agner Fog. 
* GNU General Public License http://www.gnu.org/licenses/gpl.html
*****************************************************************************/

#include <R.h>
#include <Rinternals.h>
#include "stocc.h"


/******************************************************************************
      dFNCHypergeo
      Mass function, Fisher's NonCentral Hypergeometric distribution
******************************************************************************/
REXPORTS SEXP dFNCHypergeo(
SEXP rx,         // Number of red balls drawn, scalar or vector
SEXP rm1,        // Number of red balls in urn
SEXP rm2,        // Number of white balls in urn
SEXP rn,         // Number of balls drawn from urn
SEXP rodds,      // Odds of getting a red ball among one red and one white
SEXP rprecision  // Precision of calculation
// ,SEXP rlog    // Will return log(p) if TRUE
) {
   // Check for vectors
   if (LENGTH(rx)         <  0 
    || LENGTH(rm1)        != 1 
    || LENGTH(rm2)        != 1 
    || LENGTH(rn)         != 1 
    || LENGTH(rodds)      != 1 
    || LENGTH(rprecision) != 1 
    // || LENGTH(rlog)       >  1
    ) {
       error("Parameter has wrong length");
   }
   // Get parameter values
   int     *px  =  INTEGER(rx);
   int     m1   = *INTEGER(rm1);
   int     m2   = *INTEGER(rm2);
   int     n    = *INTEGER(rn);
   double  odds = *REAL(rodds);
   double  prec = *REAL(rprecision);
   //int   ilog = *LOGICAL(rlog);
   int     nres = LENGTH(rx);          // Number of probability values to return
   int     N    = m1 + m2;             // Total number of balls
   double* buffer = 0;                 // Table of probabilities
   int     BufferLength;               // Length of table
   double  factor;                     // Scale factor
   int     x;                          // Temporary x
   int32   x1, x2;                     // Table limits
   int     xmin, xmax;                 // Absolute limits for x
   int     i;                          // Loop counter

   // Check validity of parameters
   if (!R_FINITE(odds) || odds < 0) error("Invalid value for odds");
   if (m1 < 0 || m2 < 0 || n < 0) error("Negative parameter");
   if ((unsigned int)N > 2000000000) error("Overflow");
   if (n > N) error ("n > m1 + m2: Taking more items than there are");
   if (n > m2 && odds == 0) error ("Not enough items with nonzero weight");
   if (!R_FINITE(prec) || prec < 0 || prec > 1) prec = 1E-7;

   // Allocate result vector
   SEXP result;  double * presult;
   PROTECT(result = allocVector(REALSXP, nres));
   presult = REAL(result);

   // Make object for calculating probabilities
   CFishersNCHypergeometric fnc(n, m1, N, odds, prec);

   // Check if it is advantageous to use MakeTable:
   if (nres > 1 &&
   (BufferLength = (int)fnc.MakeTable(buffer, 0, &x1, &x2),
   (uint32)nres > (uint32)BufferLength / 32)) {
      // Use MakeTable
      xmin = m1 + n - N;  if (xmin < 0) xmin = 0;  // Minimum x
      xmax = n;  if (xmax > m1) xmax = m1;         // Maximum x

      // Allocate buffer
      buffer = (double*)R_alloc(BufferLength, sizeof(double));

      // Make table of probabilities
      factor = 1. / fnc.MakeTable(buffer, BufferLength, &x1, &x2, prec*0.001);
      // Get probabilities from table
      for (i = 0; i < nres; i++) {
         x = px[i];
         if (x >= x1 && x <= x2) {
            // x within table
            presult[i] = buffer[x - x1] * factor;     // Get result from table
         }
         else if (x >= xmin && x <= xmax) {
            // Outside table. Result is very small but not 0
            presult[i] = fnc.probability(x);          // Calculate result
         }
         else {
            // Impossible value of x
            presult[i] = 0.;                          // Result is 0
         }
         // if (ilog) presult[i] = log(presult[i]);   // Log desired
      }
   }
   else {
      // Calculate probabilities one by one
      for (i = 0; i < nres; i++) {
         presult[i] = fnc.probability(px[i]);         // Probability
         //if (ilog) presult[i] = log(presult[i]);    // Log desired
      }
   }
   // Return result
   UNPROTECT(1);
   return(result);
}


/******************************************************************************
      dWNCHypergeo
      Mass function, Wallenius' NonCentral Hypergeometric distribution
******************************************************************************/
REXPORTS SEXP dWNCHypergeo(
SEXP rx,         // Number of red balls drawn, scalar or vector
SEXP rm1,        // Number of red balls in urn
SEXP rm2,        // Number of white balls in urn
SEXP rn,         // Number of balls drawn from urn
SEXP rodds,      // Odds of getting a red ball among one red and one white
SEXP rprecision  // Precision of calculation
// ,SEXP rlog    // Will return log(p) if TRUE
) {
   // Check for vectors
   if (LENGTH(rx)         <  0 
    || LENGTH(rm1)        != 1 
    || LENGTH(rm2)        != 1 
    || LENGTH(rn)         != 1 
    || LENGTH(rodds)      != 1 
    || LENGTH(rprecision) != 1 
    // || LENGTH(rlog)       >  1
    ) {
       error("Parameter has wrong length");
   }
   // Get parameter values
   int   * px   =  INTEGER(rx);
   int     m1   = *INTEGER(rm1);
   int     m2   = *INTEGER(rm2);
   int     n    = *INTEGER(rn);
   double  odds = *REAL(rodds);
   double  prec = *REAL(rprecision);
   //int   ilog = *LOGICAL(rlog);
   int     nres = LENGTH(rx);          // Number of probability values to return
   int     N    = m1 + m2;             // Total number of balls
   double* buffer = 0;                 // Table of probabilities
   int     BufferLength;               // Length of table
   int     x;                          // Temporary x
   int32   x1, x2;                     // Table limits
   int     xmin, xmax;                 // Absolute limits for x
   int     i;                          // Loop counter

   // Check validity of parameters
   if (!R_FINITE(odds) || odds < 0) error("Invalid value for odds");
   if (m1 < 0 || m2 < 0 || n < 0) error("Negative parameter");
   if ((unsigned int)N > 2000000000) error("Overflow");
   if (n > N) error ("n > m1 + m2: Taking more items than there are");
   if (n > m2 && odds == 0) error ("Not enough items with nonzero weight");
   if (!R_FINITE(prec) || prec < 0 || prec > 1) prec = 1E-7;

   // Allocate result vector
   SEXP result;  double * presult;
   PROTECT(result = allocVector(REALSXP, nres));
   presult = REAL(result);

   // Make object for calculating probabilities
   CWalleniusNCHypergeometric wnc(n, m1, N, odds, prec);

   // Check if it is advantageous to use MakeTable:
   if (nres > 1 &&
   (BufferLength = wnc.MakeTable(buffer, 0, &x1, &x2),
   x1)) {
      // Use MakeTable
      xmin = m1 + n - N;  if (xmin < 0) xmin = 0;  // Minimum x
      xmax = n;  if (xmax > m1) xmax = m1;         // Maximum x

      // Allocate buffer
      buffer = (double*)R_alloc(BufferLength, sizeof(double));
      // Make table of probabilities
      wnc.MakeTable(buffer, BufferLength, &x1, &x2, prec*0.001);
      // Get probabilities from table
      for (i = 0; i < nres; i++) {
         x = px[i];
         if (x >= x1 && x <= x2) {
            // x within table
            presult[i] = buffer[x - x1];              // Get result from table
         }
         else if (x >= xmin && x <= xmax) {
            // Outside table. Result is very small but not 0
            presult[i] = wnc.probability(x);          // Calculate result
         }
         else {
            // Impossible value of x
            presult[i] = 0.;                          // Result is 0
         }
         // if (ilog) presult[i] = log(presult[i]);   // Log desired
      }
   }
   else {
      // Calculate probabilities one by one
      for (i = 0; i < nres; i++) {
         presult[i] = wnc.probability(px[i]);
         //if (ilog) presult[i] = log(presult[i]);
      }
   }
   // Return result
   UNPROTECT(1);
   return(result);
}


/******************************************************************************
      pFNCHypergeo
      Cumulative distribution function for
      Fisher's NonCentral Hypergeometric distribution
******************************************************************************/
REXPORTS SEXP pFNCHypergeo(
SEXP rx,         // Number of red balls drawn, scalar or vector
SEXP rm1,        // Number of red balls in urn
SEXP rm2,        // Number of white balls in urn
SEXP rn,         // Number of balls drawn from urn
SEXP rodds,      // Odds of getting a red ball among one red and one white
SEXP rprecision, // Precision of calculation
SEXP rlower_tail // TRUE: P(X <= x), FALSE: P(X > x)
) {
   // Check for vectors
   if (LENGTH(rx)          <  0 
    || LENGTH(rm1)         != 1 
    || LENGTH(rm2)         != 1 
    || LENGTH(rn)          != 1 
    || LENGTH(rodds)       != 1 
    || LENGTH(rprecision)  != 1 
    || LENGTH(rlower_tail) != 1
    ) {
       error("Parameter has wrong length");
   }
   // Get parameter values
   int   * px   =  INTEGER(rx);
   int     m1   = *INTEGER(rm1);
   int     m2   = *INTEGER(rm2);
   int     n    = *INTEGER(rn);
   double  odds = *REAL(rodds);
   double  prec = *REAL(rprecision);
   int     lower_tail = *LOGICAL(rlower_tail);
   int     nres = LENGTH(rx);          // Number of probability values to return
   int     N    = m1 + m2;             // Total number of balls
   double* buffer = 0;                 // Table of probabilities
   int     BufferLength;               // Length of table
   double  factor;                     // Scale factor
   double  sum;                        // Used for summation
   double  p;                          // Probability
   int     x;                          // Temporary x
   int32   x1, x2;                     // Table limits
   int     xmin, xmax;                 // Absolute limits for x
   int     xmean;                      // Approximate mean of x
   int     i;                          // Loop counter

   // Check validity of parameters
   if (!R_FINITE(odds) || odds < 0) error("Invalid value for odds");
   if (m1 < 0 || m2 < 0 || n < 0) error("Negative parameter");
   if ((unsigned int)N > 2000000000) error("Overflow");
   if (n > N) error ("n > m1 + m2: Taking more items than there are");
   if (n > m2 && odds == 0) error ("Not enough items with nonzero weight");
   if (!R_FINITE(prec) || prec < 0 || prec > 1) prec = 1E-7;

   // min and max
   xmin = m1 + n - N;  if (xmin < 0) xmin = 0;  // Minimum x
   xmax = n;  if (xmax > m1) xmax = m1;         // Maximum x

   // Allocate result vector
   SEXP result;  double * presult;
   PROTECT(result = allocVector(REALSXP, nres));
   presult = REAL(result);

   // Make object for calculating probabilities
   CFishersNCHypergeometric fnc(n, m1, N, odds, prec);

   // Get necessary buffer length
   BufferLength = (int)fnc.MakeTable(buffer, 0, &x1, &x2, prec * 0.001);

   // Allocate buffer
   buffer = (double*)R_alloc(BufferLength, sizeof(double));

   // Make table of probabilities
   factor = 1. / fnc.MakeTable(buffer, BufferLength, &x1, &x2, prec * 0.001);

   // Get mean
   xmean = (int)(fnc.mean() + 0.5);           // Round mean

   // Check for consistency
   if (xmean < x1 || xmean > x2) {
      error("Inconsistency. mean = %i, lower limit = %i, upper limit = %i",
      xmean, x1, x2);
   }

   // Make left tail of table cumulative:
   for (x = x1, sum = 0; x <= xmean; x++) sum = buffer[x-x1] += sum;

   // Probabilities for x > xmean are calculated by summation from the
   // right in order to avoid loss of precision.
   // Make right tail of table cumulative from the right:
   for (x = x2, sum = 0; x > xmean; x--) sum = buffer[x-x1] += sum;   

   // Loop through x vector
   for (i = 0; i < nres; i++) {
      x = px[i];                       // Input x value
      if (x <= xmean) {
         // Left tail
         if (x < x1) {
            p = 0.;                    // Outside table
         }
         else {
            p = buffer[x-x1] * factor; // Probability from table
         }
         if (!lower_tail) p = 1. - p;  // Invert if right tail
         presult[i] = p;               // Store result
      }
      else {
         // Right tail
         if (x >= x2) {
            p = 0.;                    // Outside table
         }
         else {
            p = buffer[x-x1+1] * factor; // Probability from table
         }
         if (lower_tail) p = 1. - p;   // Invert if left tail
         presult[i] = p;               // Store result
      }
   }
   // Return result
   UNPROTECT(1);
   return(result);
}


/******************************************************************************
      pWNCHypergeo
      Cumulative distribution function for
      Wallenius' NonCentral Hypergeometric distribution
******************************************************************************/
REXPORTS SEXP pWNCHypergeo(
SEXP rx,         // Number of red balls drawn, scalar or vector
SEXP rm1,        // Number of red balls in urn
SEXP rm2,        // Number of white balls in urn
SEXP rn,         // Number of balls drawn from urn
SEXP rodds,      // Odds of getting a red ball among one red and one white
SEXP rprecision, // Precision of calculation
SEXP rlower_tail // TRUE: P(X <= x), FALSE: P(X > x)
) {
   // Check for vectors
   if (LENGTH(rx)          <  0 
    || LENGTH(rm1)         != 1 
    || LENGTH(rm2)         != 1 
    || LENGTH(rn)          != 1 
    || LENGTH(rodds)       != 1 
    || LENGTH(rprecision)  != 1 
    || LENGTH(rlower_tail) != 1
    ) {
       error("Parameter has wrong length");
   }
   // Get parameter values
   int   * px   =  INTEGER(rx);
   int     m1   = *INTEGER(rm1);
   int     m2   = *INTEGER(rm2);
   int     n    = *INTEGER(rn);
   double  odds = *REAL(rodds);
   double  prec = *REAL(rprecision);
   int     lower_tail = *LOGICAL(rlower_tail);
   int     nres = LENGTH(rx);          // Number of probability values to return
   int     N    = m1 + m2;             // Total number of balls
   double* buffer = 0;                 // Table of probabilities
   int     BufferLength;               // Length of table
   double  sum;                        // Used for summation
   double  p;                          // Probability
   int     x;                          // Temporary x
   int32   x1, x2;                     // Table limits
   int     xmin, xmax;                 // Absolute limits for x
   int     xmean;                      // Approximate mean of x
   int     i;                          // Loop counter

   // Check validity of parameters
   if (!R_FINITE(odds) || odds < 0) error("Invalid value for odds");
   if (m1 < 0 || m2 < 0 || n < 0) error("Negative parameter");
   if ((unsigned int)N > 2000000000) error("Overflow");
   if (n > N) error ("n > m1 + m2: Taking more items than there are");
   if (n > m2 && odds == 0) error ("Not enough items with nonzero weight");
   if (!R_FINITE(prec) || prec < 0 || prec > 1) prec = 1E-7;

   // min and max
   xmin = m1 + n - N;  if (xmin < 0) xmin = 0;  // Minimum x
   xmax = n;  if (xmax > m1) xmax = m1;         // Maximum x

   // Allocate result vector
   SEXP result;  double * presult;
   PROTECT(result = allocVector(REALSXP, nres));
   presult = REAL(result);

   // Make object for calculating probabilities
   CWalleniusNCHypergeometric wnc(n, m1, N, odds, prec);

   // Get necessary buffer length
   BufferLength = wnc.MakeTable(buffer, 0, &x1, &x2, prec * 0.001);

   // Allocate buffer
   buffer = (double*)R_alloc(BufferLength, sizeof(double));

   // Make table of probabilities
   wnc.MakeTable(buffer, BufferLength, &x1, &x2, prec * 0.001);

   // Get mean
   xmean = (int)(wnc.mean() + 0.5);           // Round mean

   // Check for consistency
   if (xmean < x1 || xmean > x2) {
      error("Inconsistency. mean = %i, lower limit = %i, upper limit = %i",
      xmean, x1, x2);
   }

   // Make left tail of table cumulative:
   for (x = x1, sum = 0; x <= xmean; x++) sum = buffer[x-x1] += sum;

   // Probabilities for x > xmean are calculated by summation from the
   // right in order to avoid loss of precision.
   // Make right tail of table cumulative from the right:
   for (x = x2, sum = 0; x > xmean; x--) sum = buffer[x-x1] += sum;   

   // Loop through x vector
   for (i = 0; i < nres; i++) {
      x = px[i];                       // Input x value
      if (x <= xmean) {
         // Left tail
         if (x < x1) {
            p = 0.;                    // Outside table
         }
         else {
            p = buffer[x-x1];          // Probability from table
         }
         if (!lower_tail) p = 1. - p;  // Invert if right tail
         presult[i] = p;               // Store result
      }
      else {
         // Right tail
         if (x >= x2) {
            p = 0.;                    // Outside table
         }
         else {
            p = buffer[x-x1+1];        // Probability from table
         }
         if (lower_tail) p = 1. - p;   // Invert if left tail
         presult[i] = p;               // Store result
      }
   }
   // Return result
   UNPROTECT(1);
   return(result);
}


/******************************************************************************
      qFNCHypergeo
      Quantile function for
      Fisher's NonCentral Hypergeometric distribution.
      Returns the lowest x for which P(X<=x) >= p when lower.tail = TRUE
      Returns the lowest x for which P(X >x) <= p when lower.tail = FALSE
******************************************************************************/
REXPORTS SEXP qFNCHypergeo(
SEXP rp,         // Cumulative probability
SEXP rm1,        // Number of red balls in urn
SEXP rm2,        // Number of white balls in urn
SEXP rn,         // Number of balls drawn from urn
SEXP rodds,      // Odds of getting a red ball among one red and one white
SEXP rprecision, // Precision of calculation
SEXP rlower_tail // TRUE: P(X <= x), FALSE: P(X > x)
) {
   // Check for vectors
   if (LENGTH(rp)          <  0 
    || LENGTH(rm1)         != 1 
    || LENGTH(rm2)         != 1 
    || LENGTH(rn)          != 1 
    || LENGTH(rodds)       != 1 
    || LENGTH(rprecision)  != 1 
    || LENGTH(rlower_tail) != 1
    ) {
       error("Parameter has wrong length");
   }
   // Get parameter values
   double* pp   = REAL(rp);
   int     m1   = *INTEGER(rm1);
   int     m2   = *INTEGER(rm2);
   int     n    = *INTEGER(rn);
   double  odds = *REAL(rodds);
   double  prec = *REAL(rprecision);
   int     lower_tail = *LOGICAL(rlower_tail);
   int     nres = LENGTH(rp);          // Number of probability values to return
   int     N    = m1 + m2;             // Total number of balls
   double* buffer = 0;                 // Table of probabilities
   int     BufferLength;               // Length of table
   double  factor;                     // Scale factor
   double  sum;                        // Used for summation
   double  p;                          // Probability
   int     x;                          // Temporary x
   int32   x1, x2;                     // Table limits
   int     i;                          // Loop counter
   unsigned int a, b, c;               // Used in binary search

   // Check validity of parameters
   if (!R_FINITE(odds) || odds < 0) error("Invalid value for odds");
   if (m1 < 0 || m2 < 0 || n < 0) error("Negative parameter");
   if ((unsigned int)N > 2000000000) error("Overflow");
   if (n > N) error ("n > m1 + m2: Taking more items than there are");
   if (n > m2 && odds == 0) error ("Not enough items with nonzero weight");
   if (!R_FINITE(prec) || prec < 0 || prec > 1) prec = 1E-7;

   // Allocate result vector
   SEXP result;  int * presult;
   PROTECT(result = allocVector(INTSXP, nres));
   presult = INTEGER(result);

   // Make object for calculating probabilities
   CFishersNCHypergeometric fnc(n, m1, N, odds, prec);

   // Get necessary buffer length
   BufferLength = (int)fnc.MakeTable(buffer, 0, &x1, &x2, prec * 0.001);

   // Allocate buffer
   buffer = (double*)R_alloc(BufferLength, sizeof(double));

   // Make table of probabilities
   factor = fnc.MakeTable(buffer, BufferLength, &x1, &x2, prec * 0.001);

   // Make table cumulative:
   for (x = x1, sum = 0; x <= x2; x++) sum = buffer[x-x1] += sum;

   // Loop through p vector
   for (i = 0; i < nres; i++) {
      p = pp[i];                       // Input p value
      if (!R_FINITE(p) || p < 0. || p > 1.) {
         presult[i] = NA_INTEGER;      // Invalid input. Return NA
      }
      else {
         if (!lower_tail) p = 1. - p;  // Invert if right tail
         p *= factor;                  // Table is scaled by factor

         // Binary search in table
         a = 0; b = x2 - x1 + 1;
         while (a < b) {
            c = (a + b) / 2;
            if (p <= buffer[c]) {
               b = c;
            }
            else {
               a = c + 1;
            }
         }
         x = x1 + a;  
         if (x > x2) x = x2;           // Prevent values > xmax that occur because of small imprecisions
         presult[i] = x;
      }
   }
   // Return result
   UNPROTECT(1);
   return(result);
}


/******************************************************************************
      qWNCHypergeo
      Quantile function for
      Wallenius' NonCentral Hypergeometric distribution.
      Returns the lowest x for which P(X<=x) >= p when lower.tail = TRUE
      Returns the lowest x for which P(X >x) <= p when lower.tail = FALSE
******************************************************************************/
REXPORTS SEXP qWNCHypergeo(
SEXP rp,         // Cumulative probability
SEXP rm1,        // Number of red balls in urn
SEXP rm2,        // Number of white balls in urn
SEXP rn,         // Number of balls drawn from urn
SEXP rodds,      // Odds of getting a red ball among one red and one white
SEXP rprecision, // Precision of calculation
SEXP rlower_tail // TRUE: P(X <= x), FALSE: P(X > x)
) {
   // Check for vectors
   if (LENGTH(rp)          <  0 
    || LENGTH(rm1)         != 1 
    || LENGTH(rm2)         != 1 
    || LENGTH(rn)          != 1 
    || LENGTH(rodds)       != 1 
    || LENGTH(rprecision)  != 1 
    || LENGTH(rlower_tail) != 1
    ) {
       error("Parameter has wrong length");
   }
   // Get parameter values
   double* pp   = REAL(rp);
   int     m1   = *INTEGER(rm1);
   int     m2   = *INTEGER(rm2);
   int     n    = *INTEGER(rn);
   double  odds = *REAL(rodds);
   double  prec = *REAL(rprecision);
   int     lower_tail = *LOGICAL(rlower_tail);
   int     nres = LENGTH(rp);          // Number of probability values to return
   int     N    = m1 + m2;             // Total number of balls
   double* buffer = 0;                 // Table of probabilities
   int     BufferLength;               // Length of table
   double  sum;                        // Used for summation
   double  p;                          // Probability
   int     x;                          // Temporary x
   int32   x1, x2;                     // Table limits
   int     i;                          // Loop counter
   unsigned int a, b, c;               // Used in binary search

   // Check validity of parameters
   if (!R_FINITE(odds) || odds < 0) error("Invalid value for odds");
   if (m1 < 0 || m2 < 0 || n < 0) error("Negative parameter");
   if ((unsigned int)N > 2000000000) error("Overflow");
   if (n > N) error ("n > m1 + m2: Taking more items than there are");
   if (n > m2 && odds == 0) error ("Not enough items with nonzero weight");
   if (!R_FINITE(prec) || prec < 0 || prec > 1) prec = 1E-7;

   // Allocate result vector
   SEXP result;  int * presult;
   PROTECT(result = allocVector(INTSXP, nres));
   presult = INTEGER(result);

   // Make object for calculating probabilities
   CWalleniusNCHypergeometric wnc(n, m1, N, odds, prec);

   // Get necessary buffer length
   BufferLength = wnc.MakeTable(buffer, 0, &x1, &x2, prec * 0.001);

   // Allocate buffer
   buffer = (double*)R_alloc(BufferLength, sizeof(double));

   // Make table of probabilities
   wnc.MakeTable(buffer, BufferLength, &x1, &x2, prec * 0.001);

   // Make table cumulative:
   for (x = x1, sum = 0; x <= x2; x++) sum = buffer[x-x1] += sum;

   // Loop through p vector
   for (i = 0; i < nres; i++) {
      p = pp[i];                       // Input p value
      if (!R_FINITE(p) || p < 0. || p > 1.) {
         presult[i] = NA_INTEGER;      // Invalid input. Return NA
      }
      else {
         if (!lower_tail) p = 1. - p;  // Invert if right tail

         // Binary search in table
         a = 0; b = x2 - x1 + 1;
         while (a < b) {
            c = (a + b) / 2;
            if (p <= buffer[c]) {
               b = c;
            }
            else {
               a = c + 1;
            }
         }
         x = x1 + a;  
         if (x > x2) x = x2;           // Prevent values > xmax that occur because of small imprecisions
         presult[i] = x;
      }
   }
   // Return result
   UNPROTECT(1);
   return(result);
}


/******************************************************************************
      rFNCHypergeo
      Random variate generation function for
      Fisher's NonCentral Hypergeometric distribution.
******************************************************************************/
REXPORTS SEXP rFNCHypergeo(
SEXP rnran,      // Number of random variates desired
SEXP rm1,        // Number of red balls in urn
SEXP rm2,        // Number of white balls in urn
SEXP rn,         // Number of balls drawn from urn
SEXP rodds,      // Odds of getting a red ball among one red and one white
SEXP rprecision  // Precision of calculation
) {
   // Check for vectors
   if (LENGTH(rnran)       != 1 
    || LENGTH(rm1)         != 1 
    || LENGTH(rm2)         != 1 
    || LENGTH(rn)          != 1 
    || LENGTH(rodds)       != 1 
    || LENGTH(rprecision)  != 1 
    ) {
       error("Parameter has wrong length");
   }
   // Get parameter values
   int     nran = *INTEGER(rnran);  if (LENGTH(rnran) > 1) nran = LENGTH(rnran);
   int     m1   = *INTEGER(rm1);
   int     m2   = *INTEGER(rm2);
   int     n    = *INTEGER(rn);
   double  odds = *REAL(rodds);
   double  prec = *REAL(rprecision);
   int     N    = m1 + m2;             // Total number of balls
   double* buffer = 0;                 // Table of probabilities
   int     BufferLength;               // Length of table
   double  sum;                        // Used for summation
   double  u;                          // Uniform random number
   int     x;                          // Temporary x
   int32   x1, x2;                     // Table limits
   unsigned int a, b, c;               // Used in binary search
   int     i;                          // Loop counter

   // Check validity of parameters
   if (!R_FINITE(odds) || odds < 0) error("Invalid value for odds");
   if (m1 < 0 || m2 < 0 || n < 0) error("Negative parameter");
   if (nran <= 0) error("Parameter nran must be positive");
   if ((unsigned int)N > 2000000000) error("Overflow");
   if (n > N) error ("n > m1 + m2: Taking more items than there are");
   if (n > m2 && odds == 0) error ("Not enough items with nonzero weight");
   if (!R_FINITE(prec) || prec < 0 || prec > 1) prec = 1E-7;

   // Allocate result vector
   SEXP result;  int * presult;
   PROTECT(result = allocVector(INTSXP, nran));
   presult = INTEGER(result);

   // Make object for generating variates
   StochasticLib3 sto(0);              // Seed is not used
   sto.SetAccuracy(prec);              // Set precision
   sto.InitRan();                      // Initialize RNG in R.dll

   if (nran > 4) {
      // Check necessary table length
      CFishersNCHypergeometric fnc(n, m1, N, odds, prec);
      BufferLength = (int)fnc.MakeTable(buffer, 0, &x1, &x2, prec * 0.001);

      if (BufferLength / 2 < nran) {
         // It is advantageous to make a table

         // Allocate buffer
         buffer = (double*)R_alloc(BufferLength, sizeof(double));

         // Make table of probabilities
         fnc.MakeTable(buffer, BufferLength, &x1, &x2, prec * 0.001);

         // Make table cumulative:
         for (x = x1, sum = 0; x <= x2; x++) sum = buffer[x-x1] += sum;

         // Loop for each variate
         for (i = 0; i < nran; i++) {

            // Make uniform random
            u = sto.Random() * sum;

            // Binary search in table
            a = 0;  b = x2 - x1 + 1;
            while (a < b) {
               c = (a + b) / 2;
               if (u < buffer[c]) {
                  b = c;
               }
               else {
                  a = c + 1;
               }
            }
            x = x1 + a;  
            if (x > x2) x = x2;   // Prevent values > xmax that occur because of small imprecisions
            presult[i] = x;
         }
         goto FINISHED_R;
      }
   }

   // Not using table.
   // Generate variates one by one
   for (i = 0; i < nran; i++) {
      presult[i] = sto.FishersNCHyp(n, m1, N, odds);
   }

   FINISHED_R:
   sto.EndRan();                       // Return RNG state to R.dll

   // Return result
   UNPROTECT(1);
   return(result);
}


/******************************************************************************
      rWNCHypergeo
      Random variate generation function for
      Wallenius' NonCentral Hypergeometric distribution.
******************************************************************************/
REXPORTS SEXP rWNCHypergeo(
SEXP rnran,      // Number of random variates desired
SEXP rm1,        // Number of red balls in urn
SEXP rm2,        // Number of white balls in urn
SEXP rn,         // Number of balls drawn from urn
SEXP rodds,      // Odds of getting a red ball among one red and one white
SEXP rprecision  // Precision of calculation
) {
   // Check for vectors
   if (LENGTH(rnran)       != 1
    || LENGTH(rm1)         != 1 
    || LENGTH(rm2)         != 1 
    || LENGTH(rn)          != 1 
    || LENGTH(rodds)       != 1 
    || LENGTH(rprecision)  != 1 
    ) {
       error("Parameter has wrong length");
   }
   // Get parameter values
   int     nran = *INTEGER(rnran);  if (LENGTH(rnran) > 1) nran = LENGTH(rnran);
   int     m1   = *INTEGER(rm1);
   int     m2   = *INTEGER(rm2);
   int     n    = *INTEGER(rn);
   double  odds = *REAL(rodds);
   double  prec = *REAL(rprecision);
   int     N    = m1 + m2;             // Total number of balls
   double* buffer = 0;                 // Table of probabilities
   int     BufferLength;               // Length of table
   double  sum;                        // Used for summation
   double  u;                          // Uniform random number
   int     x;                          // Temporary x
   int32   x1, x2;                     // Table limits
   unsigned int a, b, c;               // Used in binary search
   int     i;                          // Loop counter

   // Check validity of parameters
   if (!R_FINITE(odds) || odds < 0) error("Invalid value for odds");
   if (m1 < 0 || m2 < 0 || n < 0) error("Negative parameter");
   if (nran <= 0) error("Parameter nran must be positive");
   if ((unsigned int)N > 2000000000) error("Overflow");
   if (n > N) error ("n > m1 + m2: Taking more items than there are");
   if (n > m2 && odds == 0) error ("Not enough items with nonzero weight");
   if (!R_FINITE(prec) || prec < 0 || prec > 1) prec = 1E-7;

   // Allocate result vector
   SEXP result;  int * presult;
   PROTECT(result = allocVector(INTSXP, nran));
   presult = INTEGER(result);

   // Make object for generating variates
   StochasticLib3 sto(0);              // Seed is not used
   sto.SetAccuracy(prec);              // Set precision
   sto.InitRan();                      // Initialize RNG in R.dll

   if (nran > 4) {
      // Check necessary table length
      CWalleniusNCHypergeometric wnc(n, m1, N, odds, prec);
      BufferLength = (int)wnc.MakeTable(buffer, 0, &x1, &x2, prec * 0.001);

      if (BufferLength / 2 < nran) {
         // It is advantageous to make a table

         // Allocate buffer
         buffer = (double*)R_alloc(BufferLength, sizeof(double));

         // Make table of probabilities
         wnc.MakeTable(buffer, BufferLength, &x1, &x2, prec * 0.001);

         // Make table cumulative:
         for (x = x1, sum = 0; x <= x2; x++) sum = buffer[x-x1] += sum;

         // Loop for each variate
         for (i = 0; i < nran; i++) {

            // Make uniform random
            u = sto.Random() * sum;    // sum should be 1.0 but might be slightly less if tails are cut off in table

            // Binary search in table
            a = 0;  b = x2 - x1 + 1;
            while (a < b) {
               c = (a + b) / 2;
               if (u < buffer[c]) {
                  b = c;
               }
               else {
                  a = c + 1;
               }
            }
            x = x1 + a;  
            if (x > x2) x = x2;   // Prevent values > xmax that occur because of small imprecisions
            presult[i] = x;
         }
         goto FINISHED_R;
      }
   }

   // Not using table.
   // Generate variates one by one
   for (i = 0; i < nran; i++) {
      presult[i] = sto.WalleniusNCHyp(n, m1, N, odds);
   }

   FINISHED_R:
   sto.EndRan();                       // Return RNG state to R.dll

   // Return result
   UNPROTECT(1);
   return(result);
}


/******************************************************************************
      momentsFNCHypergeo
      Calculates the mean or variance of
      Fisher's NonCentral Hypergeometric distribution.
******************************************************************************/
// Uses simple approximations when precision >= 0.1. 
// Uses calculation by enumeration of all non-negligible x values when
// precision < 0.1.
// Note that several other approximations have been proposed in the literature.
// See e.g.:
// Levin, B. Biometrika, vol. 71, no. 3, 1984, pp. 630-632.
// Liao, J. Biometrics, vol. 48, no. 3, 1992, pp. 889-892.
// McCullagh, P. & Nelder, J.A.: Generalized Linear Models, 2'nd ed., 1989.

REXPORTS SEXP momentsFNCHypergeo(
SEXP rm1,        // Number of red balls in urn
SEXP rm2,        // Number of white balls in urn
SEXP rn,         // Number of balls drawn from urn
SEXP rodds,      // Odds of getting a red ball among one red and one white
SEXP rprecision, // Precision of calculation
SEXP rmoment     // 1 = mean, 2 = variance
) {
   // Check for vectors
   if (LENGTH(rm1)         != 1 
    || LENGTH(rm2)         != 1 
    || LENGTH(rn)          != 1 
    || LENGTH(rodds)       != 1 
    || LENGTH(rprecision)  != 1 
    ) {
       error("Parameter has wrong length");
   }
   // Get parameter values
   int     m1   = *INTEGER(rm1);
   int     m2   = *INTEGER(rm2);
   int     n    = *INTEGER(rn);
   double  odds = *REAL(rodds);
   double  prec = *REAL(rprecision);
   int     imoment = *INTEGER(rmoment);
   int     N    = m1 + m2;             // Total number of balls

   // Check validity of parameters
   if (!R_FINITE(odds) || odds < 0) error("Invalid value for odds");
   if (m1 < 0 || m2 < 0 || n < 0) error("Negative parameter");
   if ((unsigned int)N > 2000000000) error("Overflow");
   if (n > N) error ("n > m1 + m2: Taking more items than there are");
   if (n > m2 && odds == 0) error ("Not enough items with nonzero weight");
   if (imoment != 1 && imoment != 2) error ("Only moments 1 and 2 supported");
   if (!R_FINITE(prec) || prec < 0) prec = 1E-7;

   // Allocate result vector
   SEXP result;  double * presult;
   PROTECT(result = allocVector(REALSXP, 1));
   presult = REAL(result);

   // Make object for calculating mean and variance
   CFishersNCHypergeometric fnc(n, m1, N, odds, prec);

   // Check precision
   if (prec >= 0.1) {
      // Simple approximation allowed
      if (imoment == 1) {
         *presult = fnc.mean();
      }
      else {
         *presult = fnc.variance();
      }
   }
   else {
      // Exact calculation required
      // Values saved from last calculation:
      static int    old_m1   = 0;
      static int    old_m2   = 0;
      static int    old_n    = 0;
      static double old_odds = 0;
      static double old_prec = 0;
      static double old_mean = 0;
      static double old_var  = 0;

      if (m1 != old_m1 || m2 != old_m2 || n != old_n 
      || odds != old_odds || prec < old_prec) {
         // Parameters have changed. Cannot reuse results. 
         // Calculate mean and variance.

         // We are calculating both mean and variance in the same 
         // process. The values are stored for the next call in case
         // both mean and variance are requested
         fnc.moments(&old_mean, &old_var);

         // Store parameters for possible reuse in next call
         old_m1 = m1;  old_m2 = m2;  old_n = n;
         old_odds = odds;  old_prec = prec;
      }
      if (imoment == 1) {
         // Return mean
         *presult = old_mean;
      }
      else {
         // Return variance
         *presult = old_var;
      }
   }
   // Return result
   UNPROTECT(1);
   return(result);
}


/******************************************************************************
      momentsWNCHypergeo
      Calculates the mean or variance of
      Wallenius' NonCentral Hypergeometric distribution.
******************************************************************************/
// Uses simple approximations when precision >= 0.1. 
// Uses calculation by enumeration of all non-negligible x values when
// precision < 0.1.
REXPORTS SEXP momentsWNCHypergeo(
SEXP rm1,        // Number of red balls in urn
SEXP rm2,        // Number of white balls in urn
SEXP rn,         // Number of balls drawn from urn
SEXP rodds,      // Odds of getting a red ball among one red and one white
SEXP rprecision, // Precision of calculation
SEXP rmoment     // 1 = mean, 2 = variance
) {
   // Check for vectors
   if (LENGTH(rm1)         != 1 
    || LENGTH(rm2)         != 1 
    || LENGTH(rn)          != 1 
    || LENGTH(rodds)       != 1 
    || LENGTH(rprecision)  != 1 
    ) {
       error("Parameter has wrong length");
   }
   // Get parameter values
   int     m1   = *INTEGER(rm1);
   int     m2   = *INTEGER(rm2);
   int     n    = *INTEGER(rn);
   double  odds = *REAL(rodds);
   double  prec = *REAL(rprecision);
   int     imoment = *INTEGER(rmoment);
   int     N    = m1 + m2;             // Total number of balls

   // Check validity of parameters
   if (!R_FINITE(odds) || odds < 0) error("Invalid value for odds");
   if (m1 < 0 || m2 < 0 || n < 0) error("Negative parameter");
   if ((unsigned int)N > 2000000000) error("Overflow");
   if (n > N) error ("n > m1 + m2: Taking more items than there are");
   if (n > m2 && odds == 0) error ("Not enough items with nonzero weight");
   if (imoment != 1 && imoment != 2) error ("Only moments 1 and 2 supported");
   if (!R_FINITE(prec) || prec < 0) prec = 1E-7;

   // Allocate result vector
   SEXP result;  double * presult;
   PROTECT(result = allocVector(REALSXP, 1));
   presult = REAL(result);

   // Make object for calculating mean and variance
   CWalleniusNCHypergeometric wnc(n, m1, N, odds, prec);

   // Check precision
   if (prec >= 0.1) {
      // Simple approximation allowed
      if (imoment == 1) {
         *presult = wnc.mean();
      }
      else {
         *presult = wnc.variance();
      }
   }
   else {
      // Exact calculation required
      // Values saved from last calculation:
      static int    old_m1   = 0;
      static int    old_m2   = 0;
      static int    old_n    = 0;
      static double old_odds = 0;
      static double old_prec = 0;
      static double old_mean = 0;
      static double old_var  = 0;

      if (m1 != old_m1 || m2 != old_m2 || n != old_n 
      || odds != old_odds || prec < old_prec) {
         // Parameters have changed. Cannot reuse results. 
         // Calculate mean and variance.

         // We are calculating both mean and variance in the same 
         // process. The values are stored for the next call in case
         // both mean and variance are requested
         wnc.moments(&old_mean, &old_var);

         // Store parameters for possible reuse in next call
         old_m1 = m1;  old_m2 = m2;  old_n = n;
         old_odds = odds;  old_prec = prec;
      }
      if (imoment == 1) {
         // Return mean
         *presult = old_mean;
      }
      else {
         // Return variance
         *presult = old_var;
      }
   }
   // Return result
   UNPROTECT(1);
   return(result);
}


/******************************************************************************
      modeFNCHypergeo
      Calculates the mode of
      Fisher's NonCentral Hypergeometric distribution.
******************************************************************************/
REXPORTS SEXP modeFNCHypergeo(
SEXP rm1,        // Number of red balls in urn
SEXP rm2,        // Number of white balls in urn
SEXP rn,         // Number of balls drawn from urn
SEXP rodds       // Odds of getting a red ball among one red and one white
) {
   // Check for vectors
   if (LENGTH(rm1)         != 1 
    || LENGTH(rm2)         != 1 
    || LENGTH(rn)          != 1 
    || LENGTH(rodds)       != 1 
    ) {
       error("Parameter has wrong length");
   }
   // Get parameter values
   int     m1   = *INTEGER(rm1);
   int     m2   = *INTEGER(rm2);
   int     n    = *INTEGER(rn);
   double  odds = *REAL(rodds);
   int     N    = m1 + m2;             // Total number of balls

   // Check validity of parameters
   if (!R_FINITE(odds) || odds < 0) error("Invalid value for odds");
   if (m1 < 0 || m2 < 0 || n < 0) error("Negative parameter");
   if ((unsigned int)N > 2000000000) error("Overflow");
   if (n > N) error ("n > m1 + m2: Taking more items than there are");
   if (n > m2 && odds == 0) error ("Not enough items with nonzero weight");

   // Allocate result vector
   SEXP result;  int * presult;
   PROTECT(result = allocVector(INTSXP, 1));
   presult = INTEGER(result);

   // Calculate mode
   *presult = CFishersNCHypergeometric(n, m1, N, odds).mode();

   // Return result
   UNPROTECT(1);
   return(result);
}


/******************************************************************************
      modeWNCHypergeo
      Calculates the mode of
      Wallenius' NonCentral Hypergeometric distribution.
******************************************************************************/
REXPORTS SEXP modeWNCHypergeo(
SEXP rm1,        // Number of red balls in urn
SEXP rm2,        // Number of white balls in urn
SEXP rn,         // Number of balls drawn from urn
SEXP rodds,      // Odds of getting a red ball among one red and one white
SEXP rprecision  // Precision of calculation
) {
   // Check for vectors
   if (LENGTH(rm1)         != 1 
    || LENGTH(rm2)         != 1 
    || LENGTH(rn)          != 1 
    || LENGTH(rodds)       != 1 
    || LENGTH(rprecision)  != 1 
    ) {
       error("Parameter has wrong length");
   }
   // Get parameter values
   int     m1   = *INTEGER(rm1);
   int     m2   = *INTEGER(rm2);
   int     n    = *INTEGER(rn);
   double  odds = *REAL(rodds);
   double  prec = *REAL(rprecision);
   int     N    = m1 + m2;             // Total number of balls

   // Check validity of parameters
   if (!R_FINITE(odds) || odds < 0) error("Invalid value for odds");
   if (m1 < 0 || m2 < 0 || n < 0) error("Negative parameter");
   if ((unsigned int)N > 2000000000) error("Overflow");
   if (n > N) error ("n > m1 + m2: Taking more items than there are");
   if (n > m2 && odds == 0) error ("Not enough items with nonzero weight");
   if (!R_FINITE(prec) || prec < 0 || prec > 1) prec = 1E-7;

   // Allocate result vector
   SEXP result;  int * presult;
   PROTECT(result = allocVector(INTSXP, 1));
   presult = INTEGER(result);

   // Calculate mode
   *presult = CWalleniusNCHypergeometric(n, m1, N, odds, prec).mode();

   // Return result
   UNPROTECT(1);
   return(result);
}


/******************************************************************************
      oddsFNCHypergeo
      Estimate odds ratio from mean for
      Fisher's NonCentral Hypergeometric distribution.
******************************************************************************/
// Uses Cornfield's approximation. precision is ignored.
REXPORTS SEXP oddsFNCHypergeo(
SEXP rmu,        // Observed mean of x1
SEXP rm1,        // Number of red balls in urn
SEXP rm2,        // Number of white balls in urn
SEXP rn,         // Number of balls drawn from urn
SEXP rprecision  // Precision of calculation
) {
   // Check for vectors
   if (LENGTH(rmu)          < 1 
   || LENGTH(rm1)         != 1 
   || LENGTH(rm2)         != 1 
   || LENGTH(rn)          != 1 
   || LENGTH(rprecision)  != 1 
   ) {
      error("Parameter has wrong length");
   }
   // Get parameter values
   double *pmu  =  REAL(rmu);
   int     m1   = *INTEGER(rm1);
   int     m2   = *INTEGER(rm2);
   int     n    = *INTEGER(rn);
   double  prec = *REAL(rprecision);
   int     nres = LENGTH(rmu);
   int     N    = m1 + m2;             // Total number of balls
   int     i;                          // Loop counter
   int     err  = 0;                   // Remember any error

   // Check validity of parameters
   if (nres < 0) error("mu has wrong length");
   if (m1 < 0 || m2 < 0 || n < 0) error("Negative parameter");
   if ((unsigned int)N > 2000000000) error("Overflow");
   if (n > N) error ("n > m1 + m2: Taking more items than there are");
   if (!R_FINITE(prec) || prec < 0 || prec > 1) prec = 0.1;
   if (prec < 0.05) warning ("Cannot obtain high precision");

   // Allocate result vector
   SEXP result;  double * presult;
   PROTECT(result = allocVector(REALSXP, nres));
   presult = REAL(result);

   // Get xmin and xmax
   int xmin = m1 + n - N;  if (xmin < 0) xmin = 0;  // Minimum x
   int xmax = n;  if (xmax > m1) xmax = m1;         // Maximum x

   // Loop for all mu inputs
   for (i = 0; i < nres; i++) {
      double mu = pmu[i];

      // Check limits
      if (xmin == xmax) {
         presult[i] = R_NaN;  err |= 1;         // Indetermined
         continue;
      }
      if (mu <= double(xmin)) {
         if (mu == double(xmin)) {
            presult[i] = 0.;  err |= 2;         // Zero
            continue;
         }
         presult[i] = R_NaN;  err |= 8;         // Out of range
         continue;
      }
      if (mu >= double(xmax)) {
         if (mu == double(xmax)) {
            presult[i] = R_PosInf;  err |= 4;   // Infinite
            continue;
         }
         presult[i] = R_NaN;  err |= 8;         // Out of range  
         continue;
      }

      // Calculate odds ratio
      presult[i] = mu * (m2 - n + mu) / ((m1 - mu)*(n - mu));
   }
   // Check for errors
   if (err & 8) error("mu out of range");
   else if (err & 1) warning("odds is indetermined");
   else {
      if (err & 4) warning("odds is infinite");
      if (err & 2) warning("odds is zero with no precision");
   }

   // Return result
   UNPROTECT(1);
   return(result);
}


/******************************************************************************
      oddsWNCHypergeo
      Estimate odds ratio from mean for
      Wallenius' NonCentral Hypergeometric distribution.
******************************************************************************/
// Uses Manly's approximation. precision is ignored.
REXPORTS SEXP oddsWNCHypergeo(
SEXP rmu,        // Observed mean of x1
SEXP rm1,        // Number of red balls in urn
SEXP rm2,        // Number of white balls in urn
SEXP rn,         // Number of balls drawn from urn
SEXP rprecision  // Precision of calculation
) {
   // Check for vectors
   if (LENGTH(rmu)         < 1 
   || LENGTH(rm1)         != 1 
   || LENGTH(rm2)         != 1 
   || LENGTH(rn)          != 1 
   || LENGTH(rprecision)  != 1 
   ) {
      error("Parameter has wrong length");
   }
   // Get parameter values
   double *pmu  =  REAL(rmu);
   int     m1   = *INTEGER(rm1);
   int     m2   = *INTEGER(rm2);
   int     n    = *INTEGER(rn);
   double  prec = *REAL(rprecision);
   int     nres = LENGTH(rmu);
   int     N    = m1 + m2;             // Total number of balls
   int     i;                          // Loop counter
   int     err  = 0;                   // Remember any error

   // Check validity of parameters
   if (nres < 0) error("mu has wrong length");
   if (m1 < 0 || m2 < 0 || n < 0) error("Negative parameter");
   if ((unsigned int)N > 2000000000) error("Overflow");
   if (n > N) error ("n > m1 + m2: Taking more items than there are");
   if (!R_FINITE(prec) || prec < 0 || prec > 1) prec = 0.1;
   if (prec < 0.02) warning ("Cannot obtain high precision");

   // Allocate result vector
   SEXP result;  double * presult;
   PROTECT(result = allocVector(REALSXP, nres));
   presult = REAL(result);

   // Get xmin and xmax
   int xmin = m1 + n - N;  if (xmin < 0) xmin = 0;  // Minimum x
   int xmax = n;  if (xmax > m1) xmax = m1;         // Maximum x

   // Loop for all mu inputs
   for (i = 0; i < nres; i++) {
      double mu = pmu[i];

      // Check limits
      if (xmin == xmax) {
         presult[i] = R_NaN;  err |= 1;         // Indetermined
         continue;
      }
      if (mu <= double(xmin)) {
         if (mu == double(xmin)) {
            presult[i] = 0.;  err |= 2;         // Zero
            continue;
         }
         presult[i] = R_NaN;  err |= 8;         // Out of range
         continue;
      }
      if (mu >= double(xmax)) {
         if (mu == double(xmax)) {
            presult[i] = R_PosInf;  err |= 4;   // Infinite
            continue;
         }
         presult[i] = R_NaN;  err |= 8;         // Out of range  
         continue;
      }

      // Calculate odds ratio
      presult[i] = log(1. - mu / m1) / log(1. - (n-mu)/m2);
   }
   // Check for errors
   if (err & 8) error("mu out of range");
   else if (err & 1) warning("odds is indetermined");
   else {
      if (err & 4) warning("odds is infinite");
      if (err & 2) warning("odds is zero with no precision");
   }

   // Return result
   UNPROTECT(1);
   return(result);
}


/******************************************************************************
      numWNCHypergeo
      Estimate number of balls of each color from experimental mean for
      Wallenius' NonCentral Hypergeometric distribution.
******************************************************************************/
// Uses Manly's approximation. Precision is ignored.
/* Calculation method:
   Manly's approximate equation for the mean is transformed to:
   log(1-mu1/m1) = omega*(log(1-mu2/(N-m1))
   This equation is solved by Newton-Raphson iteration
*/
REXPORTS SEXP numWNCHypergeo(
SEXP rmu,        // Observed mean of x1
SEXP rn,         // Number of balls drawn from urn
SEXP rN,         // Number of balls in urn before sampling
SEXP rodds,      // Odds of getting a red ball among one red and one white
SEXP rprecision  // Precision of calculation
) {
   // Check for vectors
   if (LENGTH(rmu)         < 1 
   || LENGTH(rn)          != 1 
   || LENGTH(rN)          != 1 
   || LENGTH(rodds)       != 1 
   || LENGTH(rprecision)  != 1 
   ) {
      error("Parameter has wrong length");
   }
   // Get parameter values
   double *pmu  =  REAL(rmu);
   int     n    = *INTEGER(rn);
   int     N    = *INTEGER(rN);
   double  odds = *REAL(rodds);
   double  prec = *REAL(rprecision);
   int     nres = LENGTH(rmu);
   int     i;                          // Loop counter
   int     err  = 0;                   // Remember any error

   // Check validity of parameters
   if (nres < 0) error("mu has wrong length");
   if (n < 0 || N < 0) error("Negative parameter");
   if ((unsigned int)N > 2000000000) error("Overflow");
   if (n > N) error ("n > N: Taking more items than there are");
   if (!R_FINITE(odds) || odds < 0) error("Invalid value for odds");
   if (!R_FINITE(prec) || prec < 0 || prec > 1) prec = 0.1;
   if (prec < 0.02) warning ("Cannot obtain high precision");

   // Allocate result vector
   SEXP result;  double * presult;
   if (nres == 1) {
      PROTECT(result = allocVector(REALSXP, 2));
   }
   else {
      PROTECT(result = allocMatrix(REALSXP, 2, nres));
   }
   presult = REAL(result);

   // Loop for all mu inputs
   for (i = 0; i < nres; i++, presult += 2) {
      double mu = pmu[i];

      // Check limits
      if (n == 0) {
         presult[0] = presult[1] = R_NaN;  
         err |= 1;         // Indetermined
         continue;
      }
      if (odds == 0.) {
         presult[0] = presult[1] = R_NaN;
         if (mu == 0.) err |= 1;  // Indetermined
         else err |= 0x10;        // Out of range
         continue;
      }
      if (n == N) {        // Known exactly
         presult[0] = mu;
         presult[1] = N - mu;
         continue;
      }
      if (mu <= 0.) {
         if (mu == 0.) {
            presult[0] = 0;  presult[1] = N;
            err |= 2;         // Zero
            continue;
         }
         presult[0] = presult[1] = R_NaN;  
         err |= 8;         // Out of range
         continue;
      }
      if (mu >= double(n)) {
         if (mu == double(n)) {
            presult[0] = N;  presult[1] = 0;
            err |= 4;   // Infinite
            continue;
         }
         presult[0] = presult[1] = R_NaN;  
         err |= 8;         // Out of range  
         continue;
      }

      // Calculate m1
      double z, zd, m1, m2, lastm1, mu2 = n - mu;

      // Initial guess
      m1 = N * mu / n;
      m2 = N - m1;
      int niter = 0;

      // Newton Raphson iteration
      do {
         lastm1 = m1;
         z = log(1. - mu/m1) - odds*log(1. - mu2/m2);
         zd = mu/(m1*(m1-mu)) + odds*mu2/(m2*(m2-mu2));
         m1 -= z / zd;
         if (m1 <= mu) { // out of range
            m1 = (lastm1 + mu) * 0.5;
         }
         m2 = N - m1;
         if (m2 <= mu2) { // out of range
            m2 = (N - lastm1 + mu2) * 0.5;
            m1 = N - m2;
         }
         if (++niter > 200) error ("Convergence problem");

      } while (fabs(m1-lastm1) > N * 1E-10);

      presult[0] = m1;  presult[1] = N - m1;
   }

   // Check for errors
   if (err & 0x08) error("mu out of range");
   else {
      if (err & 0x10) warning("Zero odds conflicts with nonzero mean");
      if (err & 1) warning("odds is indetermined");
   }
   //else if (err & 6) warning("result is independent of odds");

   // Return result
   UNPROTECT(1);
   return(result);
}


/******************************************************************************
      numFNCHypergeo
      Estimate number of balls of each color from experimental mean for
      Fisher's NonCentral Hypergeometric distribution.
******************************************************************************/
// Uses Cornfield's approximation. Precision is ignored.
REXPORTS SEXP numFNCHypergeo(
SEXP rmu,        // Observed mean of x1
SEXP rn,         // Number of balls drawn from urn
SEXP rN,         // Number of balls in urn before sampling
SEXP rodds,      // Odds of getting a red ball among one red and one white
SEXP rprecision  // Precision of calculation
) {
   // Check for vectors
   if (LENGTH(rmu)         < 1 
   || LENGTH(rn)          != 1 
   || LENGTH(rN)          != 1 
   || LENGTH(rodds)       != 1 
   || LENGTH(rprecision)  != 1 
   ) {
      error("Parameter has wrong length");
   }
   // Get parameter values
   double *pmu  =  REAL(rmu);
   int     n    = *INTEGER(rn);
   int     N    = *INTEGER(rN);
   double  odds = *REAL(rodds);
   double  prec = *REAL(rprecision);
   int     nres = LENGTH(rmu);
   int     i;                          // Loop counter
   int     err  = 0;                   // Remember any error

   // Check validity of parameters
   if (nres < 0) error("mu has wrong length");
   if (n < 0 || N < 0) error("Negative parameter");
   if ((unsigned int)N > 2000000000) error("Overflow");
   if (n > N) error ("n > N: Taking more items than there are");
   if (!R_FINITE(odds) || odds < 0) error("Invalid value for odds");
   if (!R_FINITE(prec) || prec < 0 || prec > 1) prec = 0.1;
   if (prec < 0.02) warning ("Cannot obtain high precision");

   // Allocate result vector
   SEXP result;  double * presult;
   if (nres == 1) {
      PROTECT(result = allocVector(REALSXP, 2));
   }
   else {
      PROTECT(result = allocMatrix(REALSXP, 2, nres));
   }
   presult = REAL(result);

   // Loop for all mu inputs
   for (i = 0; i < nres; i++, presult += 2) {
      double mu = pmu[i];

      // Check limits
      if (n == 0) {
         presult[0] = presult[1] = R_NaN;  
         err |= 1;         // Indetermined
         continue;
      }
      if (odds == 0.) {
         presult[0] = presult[1] = R_NaN;
         if (mu == 0.) err |= 1;  // Indetermined
         else err |= 0x10;        // Out of range
         continue;
      }
      if (n == N) {        // Known exactly
         presult[0] = mu;
         presult[1] = N - mu;
         continue;
      }
      if (mu <= 0.) {
         if (mu == 0.) {
            presult[0] = 0;  presult[1] = N;
            err |= 2;         // Zero
            continue;
         }
         presult[0] = presult[1] = R_NaN;  
         err |= 8;         // Out of range
         continue;
      }
      if (mu >= double(n)) {
         if (mu == double(n)) {
            presult[0] = N;  presult[1] = 0;
            err |= 4;   // Infinite
            continue;
         }
         presult[0] = presult[1] = R_NaN;  
         err |= 8;         // Out of range  
         continue;
      }

      // Calculate m1
      double mu2 = n - mu, mu_o = mu / odds;;
      double m1 = (mu_o*(N-mu2) + mu*mu2) / (mu_o + mu2);
      presult[0] = m1;  presult[1] = N - m1;
   }

   // Check for errors
   if (err & 0x08) error("mu out of range");
   else {
      if (err & 0x10) warning("Zero odds conflicts with nonzero mean");
      if (err & 1) warning("odds is indetermined");
   }
   //else if (err & 6) warning("result is independent of odds");

   // Return result
   UNPROTECT(1);
   return(result);
}


/***********************************************************************
         DllMain
***********************************************************************/
// Define entry point DllMain if Windows and not Gnu compiler
#if defined (_WIN32) && ! defined (__GNUC__)
   extern "C"  __declspec(dllexport) 
   int __stdcall DllMain(int, int, void*) {
      return 1;
   }
#endif
