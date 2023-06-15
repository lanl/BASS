# Kellin Rumsey: 6/14/2023

#' @title Shapley-Taylor Effects with BASS
#'
#' @description Assess the sensitivity of a computer model using Shapley values
#' @param bassSob output from the \code{sobol} function.
#' @param max_order Generalize shapley effects will be computed for all interaction orders up to \code{maxOrder}.
#' @param relative_importance a function with one argument (j) giving the (unnormalized) relative importance for subsets of size j. Default is function(j) 1.
#' @param proportion should shapley effects be returned on variance or proportion scale?
#' @param mcmc.use an integer vector indexing which MCMC iterations to use for sensitivity analysis.
#' @param verbose logical; should progress be displayed?
#' @details Calculates the Shapley values based on the Sobol decomposition
#' @return A matrix containing estimated posteriors for the Shapley values.
#'
#' @keywords Shapley value
#' @seealso \link{bass} for model fitting and \link{sobol} for Sobol decomposition.
#' @export
#' @examples
#' f<-function(x){
#' 10*sin(2*pi*x[,1]*x[,2])+20*(x[,3]-.5)^2+10*x[,4]+5*x[,5]
#' }
#' sigma<-1 # noise sd
#' n<-50 # number of observations
#' x<-matrix(runif(n*10),n,10) #10 variables, only first 5 matter
#' y<-rnorm(n,f(x),sigma)
#'
#' mod <- bass(x, y)
#' sob <- sobol(mod)
#' shapley(mod)
#'
shapley_taylor <- function(bassSob, max_order=3, relative_importance=function(j) rep(1, length(j)), proportion=FALSE, mcmc.use=NULL, verbose=TRUE){
  if(is.null(mcmc.use)){
    mcmc.use <- 1:nrow(bassSob$S)
  }
  SS <- bassSob$S[mcmc.use,]
  TT <- bassSob$T[mcmc.use,]
  M <- length(mcmc.use)
  p <- ncol(TT)
  J <- p-1
  subset_sizes <- unlist(lapply(bassSob$names.ind, length))

  if(is.null(relative_importance)){
    relative_importance <- function(j) 1/p
  }
  if(proportion){
    vt <- 1
  }else{
    vt <- bassSob$var.tot[mcmc.use]
  }

  #Calculate the number of effects
  num_effects <- sum(choose(p, 1:max_order))

  if(verbose)
    cat('Shapley Start',BASS:::myTimestamp(),'Effects to compute:',num_effects,'\n')

  # Get containment indices
  CC <- matrix(0, nrow=p, ncol=sum(subset_sizes))
  cnt <- 1
  for(i in 1:length(bassSob$names.ind)){
    names_curr <- bassSob$names.ind[[i]]
    for(j in 1:subset_sizes[i]){
      subset <- string2intervec(names_curr[j])
      CC[subset,cnt] <- 1
      cnt <- cnt + 1
    }
  }

  # Start shapley calculations
  shapley <- list()
  cntr <- 1
  for(ord in 1:max_order){
    const <- 1/sum(relative_importance(0:(p-ord)))
    sets_ord <- combn(1:p, ord)
    shap_ord <- matrix(0, nrow=length(mcmc.use), ncol=ncol(sets_ord))
    for(k in 1:ncol(sets_ord)){
      I <- sets_ord[,k]
      shap_ord[,k] <- shapley_I(I, CC, SS, TT, vt, const, relative_importance)
      cntr <- cntr + 1
      if(verbose & cntr%%25==0)
        cat('Shapley Start',BASS:::myTimestamp(),'Effects computed:',cntr,'\n')
    }
    colnames(shap_ord) <- apply(sets_ord, 2, intervec2string)
    shapley[[ord]] <- shap_ord
  }
  return(shapley)
}

shapley_I <- function(I, CC, SS, TT, vt, const, relative_importance){
  n_I <- length(I)
  M <- nrow(SS)
  p <- ncol(TT)
  # CASE: j=0
  shap <- rep(0, M)
  res_I <- rep(1, ncol(CC))
  for(i in seq_along(I)){
    res_I <- res_I * CC[I[i],]
  }
  indx_term <- which(res_I == 1) # subsets which include I but none of j \in \mathcal J
  if(length(indx_term) > 0){
    term <- rep(0, M)
    for(ell in indx_term){
      term <- term + SS[,ell]
    }
    shap <- vt*const*relative_importance(0)*term
  }
  for(j in 1:(p-n_I)){
    sets_no_i <- combn((1:p)[-I], m=j)
    weight    <- const*relative_importance(j)/choose(p-n_I, j)

    for(k in 1:ncol(sets_no_i)){
      subset <- sets_no_i[,k]
      res <- res_I
      for(ell in subset){
        res <- res * (1 - CC[ell,])
      }
      indx_term <- which(res == 1) # subsets which include i but not j \in \mathcal J
      if(length(indx_term) > 0){
        term <- rep(0, M)
        for(ell in indx_term){
          term <- term + SS[,ell]
        }
        shap <- shap + vt*weight*term
      }
    }
  }
  return(shap)
}

