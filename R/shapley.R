# Kellin Rumsey: 6/14/2023


#' @title BASS Shapley Values
#'
#' @description Assess the sensitivity of a computer model using Shapley values
#' @param bassSob output from the \code{sobol} function.
#' @param relative_importance a function with one argument (j) giving the (unnormalized) relative importance for subsets of size j. Default is 1/ncol(X).
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
shapley <- function(bassSob, relative_importance=NULL, proportion=FALSE, mcmc.use=NULL, verbose=TRUE){
  S <- bassSob$S
  p <- ncol(bassSob$T)
  subset_sizes <- unlist(lapply(bassSob$names.ind, length))
  if(is.null(mcmc.use)){
    mcmc.use <- 1:nrow(S)
  }
  M <- length(mcmc.use)
  #J <- length(bassSob$names.ind) - 1
  #if(is.null(truncInt)){
  #  J <- min(J, truncInt)
  #}
  J <- p-1
  if(is.null(relative_importance)){
    relative_importance <- function(j) 1/p
    const <- 1
  }else{
    const <- 1/sum(relative_importance(0:(p-1)))
  }
  if(proportion){
    vt <- 1
  }else{
    vt <- bassSob$var.tot
  }
  if(verbose)
    cat('Shapley Start',BASS:::myTimestamp(),'Variables:',p,'\n')

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
  shapley <- matrix(0, nrow=length(mcmc.use), ncol=p)
  for(i in 1:p){
    if(verbose)
      cat('Shapley',BASS:::myTimestamp(),'Variable:',i,'\n')
    # empty set
    shapley[,i] <- bassSob$T[,i]*vt*const*relative_importance(0)
    for(j in 1:J){
      #indx_no_i <- !grepl(i, sob$names.ind[[j]])
      #sets_no_i <- bassSob$names.ind[[j]][indx_no_i] # These are sets \mathcal J
      sets_no_i <- combn((1:p)[-i], m=j)
      weight    <- const*relative_importance(j)/choose(p-1, j)
      for(k in 1:ncol(sets_no_i)){
        subset <- sets_no_i[,k]
        #subset_vec <- string2intervec(subset)
        res <- CC[i,]
        for(ell in subset){
          res <- res * (1 - CC[ell,])
        }
        indx_term <- which(res == 1) # subsets which include i but not j \in \mathcal J
        if(length(indx_term) > 0){
          term <- rep(0, M)
          for(ell in indx_term){
            term <- term + S[,ell]
          }
          shapley[,i] <- shapley[,i] + vt*weight*term
        }
      }
    }
  }
  return(shapley)
}

intervec2string <- function(xx){
  xx <- as.numeric(xx)
  xx <- sort(xx)
  str <- min(xx)
  for(i in seq_along(xx[-1]))
    str <- paste(str,'x',xx[i+1],sep='')
  return(str)
}

string2intervec <- function(str){
  brks <- c(0, stringr::str_locate_all(str, "x")[[1]][,2], nchar(str)+1)
  nn <- length(brks) - 1
  xx <- rep(NA, nn)
  for(i in 1:nn)
    xx[i] <- substr(str, brks[i]+1,brks[i+1]-1)
  return(as.numeric(xx))
}

