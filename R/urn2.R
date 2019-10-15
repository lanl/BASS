## Devin's hack to make this work when we have more than 32 variables:
# this is the only R code we need
# in src from BiasedUrn package, modify Makevars file to have a larger number than 32.  DO THIS BEFORE COMPILING

# Package BiasedUrn, file urn2.R
# R interface to multivariate noncentral hypergeometric distributions

# *****************************************************************************
#    dMWNCHypergeo
#    Mass function for
#    Multivariate Wallenius' NonCentral Hypergeometric distribution
# *****************************************************************************
#' @useDynLib BASS
dMWNCHypergeo <-
function(
   x,                   # Number of balls drawn of each color, vector or matrix
   m,                   # Number of balls of each color in urn, vector
   n,                   # Number of balls drawn from urn, scalar
   odds,                # Odds for each color, vector
   precision=1E-7) {    # Precision of calculation, scalar
   stopifnot(is.numeric(x), is.numeric(m), is.numeric(n), is.numeric(odds), is.numeric(precision));

   # Convert x to integer vector or matrix without loosing dimensions:
   if (is.matrix(x)) {
      xx <- matrix(as.integer(x), nrow=dim(x)[1], ncol=dim(x)[2]);
   }
   else {
      xx <- as.integer(x);
   }
   .Call("dMWNCHypergeo", xx, as.integer(m), as.integer(n),
   as.double(odds), as.double(precision), PACKAGE = "BASS");
}
