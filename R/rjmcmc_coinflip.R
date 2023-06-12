#' Make Weights for N2KD (coinflip) Proposal
#'
#' Finds probabilities for each coins based on count vector and tuning parameters
#'
#' @param eta vector of counts. The ith element of eta denotes the number of times predictor i is currently included in the model
#' @param p0 The expected interaction order
#' @param alpha Tuning parameter (non-negative). Larger values indicate more preference for previously selected variables (0 is uniform).
#' @param epsilon Tuning parameter (non-negative). A laplace smoother, larger values lead to more uniform sampling.
#' @param num_passes Number of passes for the deterministic weight generator.
#' @return A vector of coin flip weights
#' @details Generates coin flip weights
#' @export
#' @examples
#' p <- 20
#' p0 <- 8
#' eta <- rep(1, p)
#' w <- make_weights(eta, p0)
#' sum(w)
#' w
#'
#' eta[1:3] <- 1e5
#' w <- make_weights(eta, p0)
#' sum(w)
#' w
#'
#' eta <- abs(1e3*rt(p, 3))
#' w <- make_weights(eta, p0)
#' sum(w)
#' w
makeCoinWeights <- function(eta, p0, w3){
  epsilon    <- w3[1]
  alpha      <- w3[2]
  num_passes <- w3[3]
  p <- length(eta)
  v <- (eta^alpha + epsilon)/sum((eta^alpha + epsilon))*p0/p
  delta <- 0
  for(i in 1:num_passes){
    beta <- mean((log(p0 + delta) - log(p))/log(v))
    delta <- delta + p0 - sum(v^beta)
  }
  return(v^beta)
}


genCandBasisCoinflip<-function(minInt,maxExpectedInt,eta.vec,nint.proposal,p,xxt,degree,xx.unique.ind,vars.len,prior){
  # get number of variables in interaction
  p0 <- sample(minInt:maxExpectedInt, size=1, prob=nint.proposal)
  wts <- makeCoinWeights(eta.vec, p0, w3)
  chi.cand <- 0
  delayed_rejection_res <- 0
  while(sum(chi.cand) == 0){
    chi.cand <- rbinom(p, 1, wts)
    vars.cand <- which(chi.cand == 1)
    lprob_tmp <- 0
    for(jj in 1:maxExpectedInt){
      wts.j <- makeCoinWeights(eta.vec, jj, w3)
      term_j <- log(nint.proposal[jj]) + sum(log(wts.j[vars.cand])) + sum(log(1 - wts.j[-vars.cand]))
      lprob_tmp <- lprob_tmp + exp(term_j)
    }
    delayed_rejection_res <- delayed_rejection_res + log(lprob_tmp)
  }
  vars <- vars.cand
  n.int <- length(vars)
  if(n.int==0) #KR: This should never happen for coinflip
    return(list(basis=rep(1,ncol(xxt)),n.int=n.int,lbmcmp=0))
  #return(NULL)
  # get signs, vars, knots
  signs<-sample(c(-1,1),size=n.int,replace=T)
  if(n.int==1){
    #vars<-sample(1:p,size=1)
    #knotInd<-sample.int(vars.len[vars],size=1)
    knotInd<-sample(xx.unique.ind[[vars]],size=1)
  } else{
    #vars<-sort(sample(1:p,size=n.int,prob=z.vec,replace=F))
    #knotInd<-sapply(vars.len[vars],sample.int,size=1)
    knotInd<-sapply(xx.unique.ind[vars],sample,size=1)
  }
  #browser()
  # make basis, get reversibility term
  #knots<-sapply(1:n.int,function(nn) xxt.unique[[vars[nn]]][knotInd[nn]])
  knots<-xxt[cbind(vars,knotInd)]
  basis<-makeBasis(signs,vars,knots,xxt,degree)
  lpbmcmp<-logProbChangeModCoinflip(n.int,vars,p,vars.len,prior$nint.prior,prior$miC, delayed_rejection_res)

  return(list(basis=basis,n.int=n.int,signs=signs,vars=vars,knotInd=knotInd,knots=knots,lbmcmp=lpbmcmp))
}

## RJ reversibility term (and prior)
logProbChangeModCoinflip<-function(n.int,vars,p,vars.len,nint.prior,miC,delayedRejectionTerm){
  lprop <- (
    - sum(log(vars.len[vars]))  # probability of knots
    + delayedRejectionTerm      # probability of nint and vars for coinflipping procedure (with delayed rejection)
    - n.int*log(2)              # probability of signs
  )
  lprior <- (
    + sum(log(vars.len[vars]))   # probability of knots
    + n.int*log(2)               # probability of signs
    + log(nint.prior[n.int])     # probability of nint
    + lchoose(p, n.int)          # probability of vars
  )
  out <- lprop + lprior
  return(out)
}
