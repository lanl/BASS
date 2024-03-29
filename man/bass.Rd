% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/bass.R
\name{bass}
\alias{bass}
\title{Bayesian Adaptive Spline Surfaces (BASS)}
\usage{
bass(
  xx,
  y,
  maxInt = 3,
  maxInt.func = 3,
  maxInt.cat = 3,
  xx.func = NULL,
  degree = 1,
  maxBasis = 1000,
  npart = NULL,
  npart.func = NULL,
  nmcmc = 10000,
  nburn = 9000,
  thin = 1,
  g1 = 0,
  g2 = 0,
  s2.lower = 0,
  h1 = 10,
  h2 = 10,
  a.tau = 0.5,
  b.tau = NULL,
  w1 = 5,
  w2 = 5,
  beta.prior = "g",
  temp.ladder = NULL,
  start.temper = NULL,
  curr.list = NULL,
  save.yhat = TRUE,
  small = FALSE,
  verbose = TRUE,
  ret.str = F
)
}
\arguments{
\item{xx}{a data frame or matrix of predictors.  Categorical predictors should be included as factors.}

\item{y}{a response vector (scalar response) or matrix (functional response).  Note: If \code{sum(y^2)} is large (i.e. \code{1e10}), please center/rescale (and rescale \code{g1} and \code{g2} if necessary).}

\item{maxInt}{integer for maximum degree of interaction in spline basis functions.  Defaults to the number of predictors, which could result in overfitting.}

\item{maxInt.func}{(functional response only) integer for maximum degree of interaction in spline basis functions describing the functional response.}

\item{maxInt.cat}{(categorical input only) integer for maximum degree of interaction of categorical inputs.}

\item{xx.func}{a vector, matrix or data frame of functional variables.}

\item{degree}{degree of splines.  Stability should be examined for anything other than 1.}

\item{maxBasis}{maximum number of basis functions.  This should probably only be altered if you run out of memory.}

\item{npart}{minimum number of non-zero points in a basis function.  If the response is functional, this refers only to the portion of the basis function coming from the non-functional predictors. Defaults to 20 or 0.1 times the number of observations, whichever is smaller.}

\item{npart.func}{same as npart, but for functional portion of basis function.}

\item{nmcmc}{number of RJMCMC iterations.}

\item{nburn}{number of the \code{nmcmc} iterations to disregard.}

\item{thin}{keep every \code{thin} samples}

\item{g1}{shape for IG prior on \eqn{\sigma^2}.}

\item{g2}{scale for IG prior on \eqn{\sigma^2}.}

\item{s2.lower}{lower bound for s2. Turns IG prior for s2 into a truncated IG.}

\item{h1}{shape for gamma prior on \eqn{\lambda}.}

\item{h2}{rate for gamma prior on \eqn{\lambda}.  This is the primary way to control overfitting.  A large value of \code{h2} favors fewer basis functions.}

\item{a.tau}{shape for gamma prior on \eqn{\tau}.}

\item{b.tau}{rate for gamma prior on \eqn{\tau}. Defaults to one over the number of observations, which centers the prior for the basis function weights on the unit information prior.}

\item{w1}{nominal weight for degree of interaction, used in generating candidate basis functions.  Should be greater than 0.}

\item{w2}{nominal weight for variables, used in generating candidate basis functions.  Should be greater than 0.}

\item{beta.prior}{what type of prior to use for basis coefficients, "g" or "jeffreys"}

\item{temp.ladder}{temperature ladder used for parallel tempering.  The first value should be 1 and the values should increase.}

\item{start.temper}{when to start tempering (after how many MCMC iterations). Defaults to 1000 or half of burn-in, whichever is smaller.}

\item{curr.list}{list of starting models (one element for each temperature), could be output from a previous run under the same model setup.}

\item{save.yhat}{logical; should predictions of training data be saved?}

\item{small}{logical; if true, returns a smaller object by leaving out \code{curr.list} and other unnecessary objects.  Use in combination with \code{save.yhat} to get smaller memory footprint for very large models.}

\item{verbose}{logical; should progress be displayed?}

\item{ret.str}{logical; return data and prior structures}
}
\value{
An object of class 'bass'.  The other output will only be useful to the advanced user.  Rather, users may be interested in prediction and sensitivity analysis, which are obtained by passing the entire object to the predict.bass or sobol functions.
}
\description{
Fits a BASS model using RJMCMC.  Optionally uses parallel tempering to improve mixing.  Can be used with scalar or functional response.  Also can use categorical inputs.
}
\details{
Explores BASS model space by RJMCMC.  The BASS model has \deqn{y = f(x) + \epsilon,  ~~\epsilon \sim N(0,\sigma^2)} \deqn{f(x) = a_0 + \sum_{m=1}^M a_m B_m(x)} and \eqn{B_m(x)} is a BASS basis function (tensor product of spline basis functions). We use priors \deqn{a \sim N(0,\sigma^2/\tau (B'B)^{-1})} \deqn{M \sim Poisson(\lambda)} as well as the priors mentioned in the arguments above.
}
\examples{
\dontrun{
####################################################################################################
### univariate example
####################################################################################################
## simulate data (Friedman function)
f<-function(x){
  10*sin(pi*x[,1]*x[,2])+20*(x[,3]-.5)^2+10*x[,4]+5*x[,5]
}
sigma<-1 # noise sd
n<-500 # number of observations
x<-matrix(runif(n*10),n,10) #10 variables, only first 5 matter
y<-rnorm(n,f(x),sigma)

## fit BASS, no tempering
mod<-bass(x,y)
plot(mod)
## fit BASS, tempering
mod<-bass(x,y,temp.ladder=1.3^(0:8),start.temper=1000)
plot(mod)

## prediction
npred<-1000
xpred<-matrix(runif(npred*10),npred,10)
pred<-predict(mod,xpred,verbose=TRUE) # posterior predictive samples
true.y<-f(xpred)
plot(true.y,colMeans(pred),xlab='true values',ylab='posterior predictive means')
abline(a=0,b=1,col=2)

## sensitivity
sens<-sobol(mod)
plot(sens,cex.axis=.5)

####################################################################################################
### functional example
####################################################################################################
## simulate data (Friedman function with first variable as functional)
sigma<-1 # noise sd
n<-500 # number of observations
nfunc<-50 # size of functional variable grid
xfunc<-seq(0,1,length.out=nfunc) # functional grid
x<-matrix(runif(n*9),n,9) # 9 non-functional variables, only first 4 matter
X<-cbind(rep(xfunc,each=n),kronecker(rep(1,nfunc),x)) # to get y
y<-matrix(f(X),nrow=n)+rnorm(n*nfunc,0,sigma)

## fit BASS
mod<-bass(x,y,xx.func=xfunc)
plot(mod)

## prediction
npred<-100
xpred<-matrix(runif(npred*9),npred,9)
Xpred<-cbind(rep(xfunc,each=npred),kronecker(rep(1,nfunc),xpred))
ypred<-matrix(f(Xpred),nrow=npred)
pred<-predict(mod,xpred) # posterior predictive samples (each is a curve)
matplot(ypred,apply(pred,2:3,mean),type='l',xlab='observed',ylab='mean prediction')
abline(a=0,b=1,col=2)
matplot(t(ypred),type='l') # actual
matplot(t(apply(pred,2:3,mean)),type='l') # mean prediction

## sensitivity
sens<-sobol(mod,mcmc.use=1:10) # for speed, only use a few samples
plot(sens) # functional variable labelled "a"

sens.func<-sobol(mod,mcmc.use=1:10,func.var=1)
plot(sens.func)
}

## minimal example for CRAN testing
mod<-bass(1:2,1:2,nmcmc=2,nburn=1)
}
\seealso{
\link{predict.bass} for prediction and \link{sobol} for sensitivity analysis.
}
\keyword{analysis}
\keyword{data}
\keyword{functional}
\keyword{nonparametric}
\keyword{regression}
\keyword{splines}
