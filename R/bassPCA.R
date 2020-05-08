#######################################################
# Author: Devin Francom, Los Alamos National Laboratory
# Protected under GPL-3 license
# Los Alamos Computer Code release C19031
# github.com/lanl/BASS
#######################################################

########################################################################
## main BASS function
########################################################################


#' @title Bayesian Adaptive Spline Surfaces (BASS) with PCA decomposition of response
#'
#' @description Decomposes a multivariate or functional response onto a principal component basis and fits a BASS model to each basis coefficient.
#' @param xx a data frame or matrix of predictors with dimension n x p.  Categorical predictors should be included as factors.
#' @param y a response matrix (functional response) with dimension n x m.
#' @param dat optional (for more control) list with elements \code{xx} (same as above), \code{y} (same as above), \code{n.pc} (number of principal components used), \code{basis} (principal components with dimension m x \code{n.pc}), \code{newy} (reduced dimension \code{y} with dimension \code{n.pc} x n), \code{trunc.error} (optional truncation error with dimension n x m), \code{y.m} (vector mean removed before PCA with dimension m), \code{y.s} (vector sd scaled before PCA with dimension m).  If \code{dat} is specified, \code{xx}, \code{y} and \code{n.pc} do not need to be specified.
#' @param n.pc number of principal components to use
#' @param perc.var optionally specify percent of variance to explain instead of n.pc
#' @param n.cores integer number of cores (threads) to use
#' @param parType either "fork" or "socket".  Forking is typically faster, but not compatible with Windows. If \code{n.cores==1}, \code{parType} is ignored.
#' @param center logical whether to subtract the mean before getting the principal components, or else a numeric vector of dimension m for the center to be used
#' @param scale logical whether to divide by the standard deviation before getting the principal components, or else a numeric vector of dimension m for the scale to be used
#' @param ... arguements to be passed to \code{bass} function calls.
#' @details Gets the PCA decomposition of the response \code{y}, and fits a bass model to each PCA basis coefficient, \code{bass(dat$xx,dat$newy[i,],...)} for \code{i in 1 to n.pc}, possibly in parallel.
#' @return An object of class 'bassBasis' with two elements:
#'   \item{mod.list}{list (of length \code{n.pc}) of individual bass models}
#'   \item{dat}{same as dat above}
#' @keywords nonparametric regression splines functional data analysis
#' @seealso \link{predict.bassBasis} for prediction and \link{sobolBasis} for sensitivity analysis.
#' @export
#' @import utils
#' @example inst/examplesPCA.R
#'
bassPCA<-function(xx=NULL,y=NULL,dat=NULL,n.pc=NULL,perc.var=99,n.cores=1,parType="fork",center=T,scale=F,...){
  if(is.null(dat))
    dat<-bassPCAsetup(xx,y,n.pc,perc.var,center,scale)

  return(bassBasis(dat,n.cores,parType = parType,...))
}




bassPCAsetup<-function(xx,y,n.pc=NULL,perc.var=99,center=T,scale=F){
  if(perc.var>100 | perc.var<0)
    stop('perc.var must be between 0 and 100')
  n<-nrow(xx)
  y<-as.matrix(y)
  xx<-as.data.frame(xx)

  if(nrow(y)==1 | ncol(y)==1)
    stop('univariate y: use bass instead of bassPCA')

  if(nrow(y)!=nrow(xx))
    y<-t(y)
  if(nrow(y)!=nrow(xx))
    stop('x,y dimension mismatch')
  if(ncol(y)==nrow(y))
    warning("Caution: because y is square, please ensure that each row of x corresponds to a row of y (and not a column)")

  if(!is.null(n.pc)){
    if(n.pc>nrow(y))
      warning('n.pc too large, using all PCs intead')
  }

  if(class(center)=='logical' & length(center)==1){
    y.m<-colMeans(y)
    if(!center)
      y.m<-rep(0,ncol(y))
  } else if(class(center)=='numeric' & length(center)==ncol(y)){
    y.m<-center
  } else{
    stop("center parameter wrong dimension")
  }

  if(class(scale)=='logical' & length(scale)==1){
    y.s<-apply(y,2,sd)
    if(!scale)
      y.s<-rep(1,ncol(y))
  } else if(class(scale)=='numeric' & length(scale)==ncol(y)){
    y.s<-scale
  } else{
    stop("scale parameter wrong dimension")
  }

  yc<-t(scale(y,center=y.m,scale=y.s)) # maybe could get away with fewer transposes, but could require a lot of refactoring
  S<-svd(yc)


  if(is.null(n.pc)){
    ev<-S$d^2
    n.pc<-which(cumsum(ev/sum(ev))*100>perc.var)[1]
  }

  basis<-S$u[,1:n.pc,drop=F]%*%diag(S$d[1:n.pc],nrow=n.pc) # columns are basis functions
  newy<-t(S$v[,1:n.pc,drop=F])

  trunc.error<-basis%*%newy - yc

  ret<-list(xx=xx,y=y,n.pc=n.pc,basis=basis,newy=newy,trunc.error=trunc.error,y.m=y.m,y.s=y.s,ev=S$d^2)
  class(ret)<-'bassPCAsetup'
  return(ret)
}



#' @title Bayesian Adaptive Spline Surfaces (BASS) with basis decomposition of response
#'
#' @description Fits a BASS model to basis coefficients under the specified basis.
#' @param dat list that includes elements \code{xx}, \code{n.pc} (number of basis functions), \code{basis} (dimension m x \code{n.pc}), \code{newy} (dimension \code{n.pc} x n), \code{trunc.error} (optional truncation error with dimension n x m), \code{y.m} (vector mean removed before basis decomposition with dimension m), \code{y.s} (vector sd scaled before basis decomposition with dimension m).  See the documentation of \code{bassPCA} for more details.
#' @param n.cores integer number of cores (threads) to use
#' @param parType either "fork" or "socket".  Forking is typically faster, but not compatible with Windows. If \code{n.cores==1}, \code{parType} is ignored.
#' @param ... arguements to be passed to \code{bass} function calls.
#' @details Under a user defined basis decomposition, fits a bass model to each PCA basis coefficient independently, \code{bass(dat$xx,dat$newy[i,],...)} for \code{i in 1 to n.pc}, possibly in parallel.  The basis does not need to be orthogonal, but independent modeling of basis coefficients should be sensible.
#' @return An object of class 'bassBasis' with two elements:
#'   \item{mod.list}{list (of length \code{n.pc}) of individual bass models}
#'   \item{dat}{same as dat above}
#' @keywords nonparametric regression splines functional data analysis
#' @seealso \link{predict.bassBasis} for prediction and \link{sobolBasis} for sensitivity analysis.
#' @export
#' @import stats
#' @import utils
#' @example inst/examplesPCA.R
bassBasis<-function(dat,n.cores=1,parType='fork',...){


  if(n.cores>parallel::detectCores())
    warning(paste0("Specified n.cores = ",n.cores,'. Proceeding with n.cores = min(n.cores,dat$n.pc,detectCores()) = ',min(n.cores,dat$n.pc,parallel::detectCores())))
  n.cores<-min(n.cores,dat$n.pc,parallel::detectCores())

  if(n.cores==1){
    mod.list<-lapply(1:dat$n.pc,function(i) bass(dat$xx,dat$newy[i,],...))
  } else if(parType=='socket'){
    cl <- parallel::makeCluster(n.cores,setup_strategy = "sequential")
    mod.list<-parallel::parLapply(cl,1:dat$n.pc,function(i) bass(dat$xx,dat$newy[i,],...))
    parallel::stopCluster(cl)
  } else if(parType=='fork'){
    mod.list<-parallel::mclapply(1:dat$n.pc,function(i) bass(dat$xx,dat$newy[i,],...),mc.cores = n.cores,mc.preschedule = F)
  }

  ret<-list(mod.list=mod.list,dat=dat)
  class(ret)<-'bassBasis'
  return(ret)
}













################################################################################################################
## prediction



#' @title BASS Prediction
#'
#' @description Predict function for BASS.  Outputs the posterior predictive samples based on the specified MCMC iterations.
#' @param object a fitted model, output from the \code{bass} function.
#' @param newdata a matrix of new input values at which to predict.  The columns should correspond to the same variables used in the \code{bassBasis} or \code{bassPCA}  functions.
#' @param mcmc.use a vector indexing which MCMC iterations to use for prediction.
#' @param trunc.error logical, use basis truncation error when predicting?
#' @param nugget logical, use individual \code{bass} nugget variances when predicting?
#' @param n.cores number of cores, though 1 is often the fastest.
#' @param parType either "fork" or "socket".  Forking is typically faster, but not compatible with Windows. If \code{n.cores==1}, \code{parType} is ignored.
#' @param ... further arguments passed to or from other methods.
#' @details Prediction combined across \code{bass} models.
#' @return An array with first dimension corresponding to MCMC iteration, second dimension corresponding to the rows of \code{newdata}, and third dimension corresponding to the multivariate/functional response.
#' @seealso \link{bassPCA} and \link{bassBasis} for model fitting and \link{sobolBasis} for sensitivity analysis.
#' @export
#' @examples
#' # See examples in bass documentation.
#'
predict.bassBasis<-function(object,newdata,mcmc.use=NULL,trunc.error=FALSE,nugget=T,n.cores=1,parType="fork",...){

  if(is.null(mcmc.use)){ # if null, use all
    mcmc.use<-1:((object$mod.list[[1]]$nmcmc-object$mod.list[[1]]$nburn)/object$mod.list[[1]]$thin)
  }

  if(n.cores==1){
    # no parallel
    newy.pred<-array(unlist(lapply(1:object$dat$n.pc,function(i) predict1mod(object$mod.list[[i]],newdata,mcmc.use,nugget,...))),dim=c(length(mcmc.use),nrow(newdata),object$dat$n.pc))

    #browser()
    out<-array(unlist(lapply(1:length(mcmc.use),function(i) predict1mcmc(matrix(newy.pred[i,,],ncol=object$dat$n.pc,nrow=nrow(newdata)),object$dat))),dim=c(length(object$dat$y.m),nrow(newdata),length(mcmc.use)))
  } else if(parType=='socket'){

    # parLapply (socket)

    cl <- parallel::makeCluster(min(n.cores,object$dat$n.pc,parallel::detectCores()),setup_strategy = "sequential") # possibly a faster way to do this, but would need to keep cluster around
    parallel::clusterExport(cl,varlist=c("newdata"),envir=environment())

    newy.pred<-array(unlist(parallel::parLapply(cl,1:object$dat$n.pc,function(i) predict1mod(object$mod.list[[i]],newdata,mcmc.use,nugget,...))),dim=c(length(mcmc.use),nrow(newdata),object$dat$n.pc))


    out<-array(unlist(parallel::parLapply(cl,1:length(mcmc.use),function(i) predict1mcmc(matrix(newy.pred[i,,],ncol=object$dat$n.pc,nrow=nrow(newdata)),object$dat))),dim=c(length(object$dat$y.m),nrow(newdata),length(mcmc.use)))

    parallel::stopCluster(cl)
  } else if(parType=='fork'){
    # mclapply (fork - faster than socket, but not compatible with windows)
    newy.pred<-array(unlist(parallel::mclapply(1:object$dat$n.pc,function(i) predict1mod(object$mod.list[[i]],newdata,mcmc.use,nugget,...),mc.cores=n.cores)),dim=c(length(mcmc.use),nrow(newdata),object$dat$n.pc))

    out<-array(unlist(parallel::mclapply(1:length(mcmc.use),function(i) predict1mcmc(matrix(newy.pred[i,,],ncol=object$dat$n.pc,nrow=nrow(newdata)),object$dat),mc.cores=n.cores)),dim=c(length(object$dat$y.m),nrow(newdata),length(mcmc.use)))
  }


  out<-aperm(out,c(3,2,1))

  if(trunc.error)
    out<-out+array(truncErrSampN(length(mcmc.use)*nrow(newdata),object$dat$trunc.error),dim=c(length(mcmc.use),nrow(newdata),length(object$dat$y.m)))

  return(out) # should be nmcmc x npred x nfunc
}

predict1mcmc<-function(mat,dat){
  if(is.null(dim(mat)))
    mat<-t(mat)
  dat$basis%*%t(mat)*dat$y.s + dat$y.m
}

predict1mod<-function(mod,newdata,mcmc.use,nugget,...){
  pmat<-predict(mod,newdata,mcmc.use=mcmc.use,...)
  if(nugget)
    pmat<-pmat+rnorm(length(mcmc.use),0,sqrt(mod$s2[mcmc.use]))
  pmat
}

predict_fast.bassBasis<-function(object,newdata,n.cores=1,mcmc.use,trunc.error=FALSE,...){
  newy.pred<-array(unlist(parallel::mclapply(1:object$dat$n.pc,function(i) predict1mod_fast(object$mod.list[[i]],newdata,mcmc.use,...),mc.cores=min(n.cores,object$dat$n.pc))),dim=c(length(mcmc.use),nrow(newdata),object$dat$n.pc))

  out<-array(unlist(parallel::mclapply(1:length(mcmc.use),function(i) predict1mcmc(newy.pred[i,,],object$dat),mc.cores=min(n.cores,length(mcmc.use)))),dim=c(length(object$dat$y.m),nrow(newdata),length(mcmc.use)))

  out<-aperm(out,c(3,2,1))

  return(out) # should be nmcmc x npred x nfunc
}

predict1mod_fast<-function(mod,newdata,mcmc.use,...){
  pmat<-predict_fast(mod,newdata,mcmc.use=mcmc.use,...)
  #pmat<-pmat+rnorm(length(mcmc.use),0,sqrt(mod$s2[mcmc.use]))
  pmat
}


truncErrSampN<-function(n,te.mat){ # a function to quickly sample one of the truncation
  #   errors at each space and time, vectorized
  # te.mat is truncation error matrix, ncol is number of sims, nrow is length of EOFs
  #   (space time combinations)

  t(te.mat[,sample.int(nrow(te.mat),n,replace=T)])
}















################################################################################################################
## sobol


#' @title BASS Sensitivity Analysis
#'
#' @description Decomposes the variance of the BASS model into variance due to main effects, two way interactions, and so on, similar to the ANOVA decomposition for linear models.  Uses the Sobol' decomposition, which can be done analytically for MARS models.
#' @param mod output from the \code{bassBasis} or \code{bassPCA} function.
#' @param int.order an integer indicating the highest order of interactions to include in the Sobol decomposition.
#' @param prior a list with the same number of elements as there are inputs to mod.  Each element specifies the prior for the particular input.  Each prior is specified as a list with elements \code{dist} (one of \code{c("normal", "student", "uniform")}), \code{trunc} (a vector of dimension 2 indicating the lower and upper truncation bounds, taken to be the data bounds if omitted), and for "normal" or "student" priors, \code{mean} (scalar mean of the Normal/Student, or a vector of means for a mixture of Normals or Students), \code{sd} (scalar standard deviation of the Normal/Student, or a vector of standard deviations for a mixture of Normals or Students), \code{df} (scalar degrees of freedom of the Student, or a vector of degrees of freedom for a mixture of Students), and \code{weights} (a vector of weights that sum to one for the mixture components, or the scalar 1).  If unspecified, a uniform is assumed with the same bounds as are represented in the input to mod.
#' @param mcmc.use an integer indicating which MCMC iteration to use for sensitivity analysis. Defaults to the last iteration.
#' @param nind number of Sobol indices to keep (will keep the largest nind).
#' @param n.cores number of cores to use (nearly linear speedup for adding cores).
#' @param parType either "fork" or "socket".  Forking is typically faster, but not compatible with Windows. If \code{n.cores==1}, \code{parType} is ignored.
#' @param plot logical; whether to plot results.
#' @param verbose logical; print progress.
#' @details Performs analytical Sobol' decomposition for each MCMC iteration in mcmc.use (each corresponds to a MARS model), yeilding a posterior distribution of sensitivity indices.  Can obtain Sobol' indices as a function of one functional variable.
#' @return If non-functional (\code{func.var = NULL}), a list with two elements:
#'  \item{S}{a data frame of sensitivity indices with number of rows matching the length of \code{mcmc.use}.  The columns are named with a particular main effect or interaction.  The values are the proportion of variance in the model that is due to each main effect or interaction.}
#'  \item{T}{a data frame of total sensitivity indices with number of rows matching the length of \code{mcmc.use}.  The columns are named with a particular variable.}
#'  Otherwise, a list with four elements:
#'  \item{S}{an array with first dimension corresponding to MCMC samples (same length as \code{mcmc.use}), second dimension corresponding to different main effects and interactions (labeled in \code{names.ind}), and third dimension corresponding to the grid used for the functional variable.  The elements of the array are sensitivity indices.}
#'  \item{S.var}{same as \code{S}, but scaled in terms of total variance rather than percent of variance.}
#'  \item{names.ind}{a vector of names of the main effects and interactions used.}
#'
#' @keywords Sobol decomposition
#' @seealso \link{bassPCA} and \link{bassBasis} for model fitting and \link{predict.bassBasis} for prediction.
#' @export
#' @examples
#' # See examples in bass documentation.

sobolBasis<-function(mod,int.order,prior=NULL,mcmc.use=NULL,nind=NULL,n.cores=1,parType='fork',plot=F,verbose=T){

  if(is.null(mcmc.use))
    mcmc.use<-length(mod$mod.list[[1]]$s2)


  bassMod<-mod$mod.list[[1]] # for structuring everything, assuming that model structures are the same for different PCs
  pdescat<-sum(bassMod$pdes)+sum(bassMod$pcat) # sums make NULLs 0s

  if(is.null(prior))
    prior<-list()

  if(length(prior)<pdescat){
    for(i in (length(prior)+1):pdescat)
      prior[[i]]<-list(dist=NA)
  }

  #browser()
  for(i in 1:pdescat){
    if(is.null(prior[[i]]))
      prior[[i]]<-list(dist=NA)

    if(is.na(prior[[i]]$dist)){
      prior[[i]]<-list()
      prior[[i]]$dist<-'uniform'
      #prior[[i]]$trunc<-bassMod$range.des[,i] - not right index when there are categorical vars
    }
  }


  if(bassMod$func){
    if(is.null(prior.func)){
      prior.func<-list()
      for(i in 1:bassMod$pfunc){
        prior.func[[i]]<-list()
        prior.func[[i]]$dist<-'uniform'
        #prior.func[[i]]$trunc<-bassMod$range.func[,i]
      }
    }
    for(i in 1:length(prior.func))
      class(prior.func[[i]])<-prior.func[[i]]$dist
  }

  for(i in 1:length(prior))
    class(prior[[i]])<-prior[[i]]$dist # class will be used for integral functions, should be uniform, normal, or student

  if(bassMod$cat){
    which.cat<-which(bassMod$cx=='factor')
    prior.cat<-list()
    for(i in 1:length(which.cat)){
      prior.cat[i]<-prior[which.cat[i]]
    }
    prior[which.cat]<-NULL
  } else{
    prior.cat<-NULL
  }

  #browser()
  if(bassMod$des){
    for(i in 1:length(prior)){
      if(is.null(prior[[i]]$trunc)){
        prior[[i]]$trunc<-c(0,1)
      } else{
        prior[[i]]$trunc<-scale.range(prior[[i]]$trunc,bassMod$range.des[,i])
      }

      if(prior[[i]]$dist %in% c('normal','student')){
        prior[[i]]$mean<-scale.range(prior[[i]]$mean,bassMod$range.des[,i])
        prior[[i]]$sd<-prior[[i]]$sd/(bassMod$range.des[2,i]-bassMod$range.des[1,i])
        if(prior[[i]]$dist == 'normal'){
          prior[[i]]$z<-pnorm((prior[[i]]$trunc[2]-prior[[i]]$mean)/prior[[i]]$sd) - pnorm((prior[[i]]$trunc[1]-prior[[i]]$mean)/prior[[i]]$sd)
        } else{
          prior[[i]]$z<-pt((prior[[i]]$trunc[2]-prior[[i]]$mean)/prior[[i]]$sd,prior[[i]]$df) - pt((prior[[i]]$trunc[1]-prior[[i]]$mean)/prior[[i]]$sd,prior[[i]]$df)
        }
        cc<-sum(prior[[i]]$weights*prior[[i]]$z)
        prior[[i]]$weights<-prior[[i]]$weights/cc#prior[[i]]$z # change weights with truncation # divide by cc instead to keep the same prior shape
        # does the truncation change the distribution shape in the non-truncated regions??
        #browser()
      }
    }
  }





  #browser()





  tl<-list(prior=prior)

















  pc.mod<-mod$mod.list
  pcs<-mod$dat$basis

  if(verbose)
    cat('Start',timestamp(quiet = T),'\n')
  p<-pc.mod[[1]]$p
  if(int.order>p){
    int.order<-p
    warning("int.order > number of inputs, changing to int.order = number of inputs")
  }

  u.list<-lapply(1:int.order,function(i) combn(1:p,i))
  ncombs.vec<-unlist(lapply(u.list,ncol))
  ncombs<-sum(ncombs.vec)
  nxfunc<-nrow(pcs)
  #sob<-matrix(nrow=nxfunc,ncol=ncombs)
  sob<-ints<-list()

  n.pc<-ncol(pcs)
  w0<-unlist(lapply(1:n.pc,function(pc) get.f0(prior,pc.mod,pc,mcmc.use)))

  #browser()

  f0r2<-(pcs%*%w0)^2

  max.nbasis<-max(unlist(lapply(pc.mod,function(x) x$nbasis[mcmc.use])))
  C1Basis.array<-array(dim=c(n.pc,p,max.nbasis))
  for(i in 1:n.pc){
    nb<-pc.mod[[i]]$nbasis[mcmc.use]
    mcmc.mod.usei<-pc.mod[[i]]$model.lookup[mcmc.use]
    for(j in 1:p){
      for(k in 1:nb){
        C1Basis.array[i,j,k]<-C1Basis(prior,pc.mod,j,k,i,mcmc.mod.usei)
      }
    }
    #print(i)
  }

  # browser()
  #
  # C2Basis.array<-array(dim=c(n.pc,n.pc,p,max.nbasis,max.nbasis))
  # for(i1 in 1:n.pc){
  #   nb1<-pc.mod[[i1]]$nbasis[mcmc.use]
  #   mcmc.mod.usei1<-pc.mod[[i1]]$model.lookup[mcmc.use]
  #   for(i2 in 1:n.pc){
  #     nb2<-pc.mod[[i2]]$nbasis[mcmc.use]
  #     mcmc.mod.usei2<-pc.mod[[i2]]$model.lookup[mcmc.use]
  #     for(j in 1:p){
  #       for(k1 in 1:nb1){
  #         for(k2 in 1:nb2){
  #           C2Basis.array[i1,i2,j,k1,k2]<-C2Basis(pc.mod,j,k1,k2,i1,i2,mcmc.mod.usei1,mcmc.mod.usei2) #C2Basis(pc.mod,l,mi,mj,i,j,mcmc.mod.usei,mcmc.mod.usej)
  #         }
  #       }
  #     }
  #   }
  #   print(i1)
  # }


  #browser()
  u.list1<-list()
  for(i in 1:int.order)
    u.list1<-c(u.list1,split(u.list[[i]], col(u.list[[i]])))
  #require(parallel)
  #browser()
  if(verbose)
    cat('Integrating',timestamp(quiet = T),'\n')

  u.list.temp<-c(list(1:p),u.list1)





  if(n.cores==1){
    # no parallel
    ints1.temp<-lapply(u.list.temp,function(x) func.hat(prior,x,pc.mod,pcs,mcmc.use,f0r2,C1Basis.array))
  } else if(parType=='socket'){

    # parLapply (socket)

    cl <- parallel::makeCluster(min(n.cores,parallel::detectCores()),setup_strategy = "sequential") # possibly a faster way to do this, but would need to keep cluster around
    parallel::clusterExport(cl,varlist=c("prior","x","pc.mod","pcs","mcmc.use","f0r2","C1Basis.array"),envir=environment())

    ints1.temp<-parallel::parLapply(cl,u.list.temp,function(x) func.hat(prior,x,pc.mod,pcs,mcmc.use,f0r2,C1Basis.array))

    parallel::stopCluster(cl)
  } else if(parType=='fork'){
    # mclapply (fork - faster than socket, but not compatible with windows)
    ints1.temp<-parallel::mclapply(u.list.temp,function(x) func.hat(prior,x,pc.mod,pcs,mcmc.use,f0r2,C1Basis.array),mc.cores=n.cores,mc.preschedule = F)
  }




  V.tot<-ints1.temp[[1]]
  ints1<-ints1.temp[-1]

  #ints1<-mclapply(u.list1,function(x) func.hat(prior,x,pc.mod,pcs,mcmc.use,f0r2,C1Basis.array),mc.cores=n.cores,mc.preschedule = preschedule)
  ints<-list()
  ints[[1]]<-do.call(cbind,ints1[1:ncol(u.list[[1]])])
  if(int.order>1){
    for(i in 2:int.order)
      ints[[i]]<-do.call(cbind,ints1[sum(ncombs.vec[1:(i-1)])+1:ncol(u.list[[i]])])
  }

  # for(i in 1:length(u.list))
  #   ints[[i]]<-apply(u.list[[i]],2,function(x) func.hat(x,pc.mod,pcs,mcmc.use,f0r2)) # the heavy lifting


  sob[[1]]<-ints[[1]]
  # matplot(t(apply(sob[[1]],1,cumsum)),type='l')
  # matplot(t(apply(sens.func$S.var[1,1:5,],2,cumsum)),type='l',add=T)

  #V.tot<-func.hat(prior,1:p,pc.mod,pcs,mcmc.use,f0r2) # need to add this to the above



  # plot(V.tot)
  # points(apply(sens.func$S.var[1,,],2,sum),col=2)

  if(verbose)
    cat('Shuffling',timestamp(quiet = T),'\n')
  if(length(u.list)>1){
    for(i in 2:length(u.list)){
      sob[[i]]<-matrix(nrow=nxfunc,ncol=ncol(ints[[i]]))
      for(j in 1:ncol(u.list[[i]])){
        cc<-rep(0,nxfunc)
        for(k in 1:(i-1)){
          ind<-which(apply(u.list[[k]],2,function(x) all(x%in%u.list[[i]][,j])))
          cc<-cc+(-1)^(i-k)*rowSums(ints[[k]][,ind])
        }
        sob[[i]][,j]<-ints[[i]][,j]+cc
      }
    }
  }


  # sens.func.use<-lapply(strsplit(sens.func$names.ind,'x'),as.numeric)
  # sl<-sapply(sens.func.use,length)
  # ind.list<-list()
  # sob.small<-list()
  # for(i in 1:length(u.list)){
  #   ind.list[[i]]<-NA
  #   k<-0
  #   for(j in which(sl==i)){
  #     k<-k+1
  #     ind.list[[i]][k]<-which(apply(u.list[[i]],2,function(x) all(x==sens.func.use[[j]])))
  #   }
  #   sob.small[[i]]<-sob[[i]][,ind.list[[i]]]
  # }
  #
  # sob.small<-do.call(cbind,sob.small)
  # matplot(t(apply(sob.small,1,cumsum)),type='l')
  # matplot(t(apply(sens.func$S.var[1,,],2,cumsum)),type='l',add=T)

  #browser()

  if(is.null(nind))
    nind<-ncombs


  sob.comb.var<-do.call(cbind,sob)

  vv<-colMeans(sob.comb.var)
  ord<-order(vv,decreasing = T)
  cutoff<-vv[ord[nind]]
  if(nind>length(ord))
    cutoff<-min(vv)
  use<-sort(which(vv>=cutoff))


  V.other<-V.tot-rowSums(sob.comb.var[,use])
  use<-c(use,ncombs+1)

  sob.comb.var<-t(cbind(sob.comb.var,V.other))
  sob.comb<-t(t(sob.comb.var)/c(V.tot))

  sob.comb.var<-sob.comb.var[use,,drop=F]
  sob.comb<-sob.comb[use,,drop=F]

  dim(sob.comb)<-c(1,length(use),nxfunc)
  dim(sob.comb.var)<-c(1,length(use),nxfunc)


  names.ind<-c(unlist(lapply(u.list,function(x) apply(x,2,paste,collapse='x',sep=''))),'other')
  names.ind<-names.ind[use]

  if(verbose)
    cat('Finish',timestamp(quiet = T),'\n')

  #browser()

  ret<-list(S=sob.comb,S.var=sob.comb.var,Var.tot=V.tot,names.ind=names.ind,xx=seq(0,1,length.out = nxfunc),func=T)
  class(ret)<-'bassSob'

  if(plot)
    plot(ret)

  return(ret)

}

################################################################################
## Functions
################################################################################
func.hat<-function(prior,u,pc.mod,pcs,mcmc.use,f0r2,C1Basis.array){ # could speed this up
  #browser()
  res<-rep(0,nrow(pcs))
  n.pc<-length(pc.mod)
  for(i in 1:n.pc){
    res<-res+pcs[,i]^2*Ccross(prior,pc.mod,i,i,u,mcmc.use,C1Basis.array)

    if(i<n.pc){
      for(j in (i+1):n.pc){
        res<-res+2*pcs[,i]*pcs[,j]*Ccross(prior,pc.mod,i,j,u,mcmc.use,C1Basis.array)
        #print(c(i,j))
      }
    }
  }
  return(res-f0r2)
}

Ccross<-function(prior,pc.mod,i,j,u,mcmc.use=1,C1Basis.array){ # inner product of main effects from different eof models
  p<-pc.mod[[1]]$p
  mcmc.mod.usei<-pc.mod[[i]]$model.lookup[mcmc.use]
  mcmc.mod.usej<-pc.mod[[j]]$model.lookup[mcmc.use]

  Mi<-pc.mod[[i]]$nbasis[mcmc.use]
  Mj<-pc.mod[[j]]$nbasis[mcmc.use]
  mat<-matrix(nrow=Mi,ncol=Mj)

  #CC<-C2Basis.temp<-CCu<-matrix(1,nrow=Mi,ncol=Mj)

  a0i<-pc.mod[[i]]$beta[mcmc.use,1]
  a0j<-pc.mod[[j]]$beta[mcmc.use,1]
  f0i<-get.f0(prior,pc.mod,i,mcmc.use)
  f0j<-get.f0(prior,pc.mod,j,mcmc.use)

  out<- a0i*a0j + a0i*(f0j-a0j) + a0j*(f0i-a0i)
  #browser()

  if(Mi>0 & Mj>0){
    ai<-pc.mod[[i]]$beta[mcmc.use,1+1:Mi]
    aj<-pc.mod[[j]]$beta[mcmc.use,1+1:Mj]

    for(mi in 1:Mi){
      for(mj in 1:Mj){
        temp1<-ai[mi]*aj[mj]
        temp2<-temp3<-1
        for(l in (1:p)[-u]){
          #temp2<-temp2*C1Basis(pc.mod,l,mi,i,mcmc.mod.usei)*C1Basis(pc.mod,l,mj,j,mcmc.mod.usej) # make a C1Basis lookup table instead (this is the bottleneck)
          temp2<-temp2*C1Basis.array[i,l,mi]*C1Basis.array[j,l,mj]
          #browser()
        }
        #CC[mi,mj]<-temp2
        for(l in u){
          temp3<-temp3*C2Basis(prior,pc.mod,l,mi,mj,i,j,mcmc.mod.usei,mcmc.mod.usej) # would be nice to use a lookup table here too, but its too big
        }
        #C2Basis.temp[mi,mj]<-temp3
        #CCu[mi,mj]<-temp4
        out<-out+temp1*temp2*temp3#(temp3-1) not -1 since we subtract f0^2 later
        #print(out)
        #mat[mi,mj]<-temp
        #print(c(temp1*temp2*temp3))
      }
    }
  }
  #out<-out+ai%*%(CC*C2Basis.temp/CCu)%*%aj
  if(length(out)==0)
    browser()
  return(out)
}


C1Basis<-function(prior,pc.mod,l,m,pc,mcmc.mod.use){ # l = variable, m = basis function, pc = eof index
  if(l<=pc.mod[[pc]]$pdes){
    int.use.l<-which(pc.mod[[pc]]$vars.des[mcmc.mod.use,m,]==l)
    if(length(int.use.l)==0)
      return(1)
    s<-pc.mod[[pc]]$signs[mcmc.mod.use,m,int.use.l]
    t.ind<-pc.mod[[pc]]$knotInd.des[mcmc.mod.use,m,int.use.l]
    t<-pc.mod[[pc]]$xx.des[t.ind,l]
    q<-pc.mod[[pc]]$degree
    #return((1/(q+1)*((s+1)/2-s*t))*s^2)

    if(s==0)
      return(0)
    cc<-const(signs=s,knots=t,degree=q)
    if(s==1){
      a<-max(prior[[l]]$trunc[1],t)
      b<-prior[[l]]$trunc[2]
      if(b<t)
        return(0)
      out<-intabq1(prior[[l]],a,b,t,q)/cc
      #return(intabq1(tl$prior[[k]],a,b,t,q)/cc)
    } else{
      a<-prior[[l]]$trunc[1]
      b<-min(prior[[l]]$trunc[2],t)
      if(t<a)
        return(0)
      out<-intabq1(prior[[l]],a,b,t,q)*(-1)^q/cc
      #return(intabq1(tl$prior[[k]],a,b,t,q)*(-1)^q/cc)
    }
    if(out< -1e-15)
      browser()
    return(out)


  } else{
    l.cat<-l-pc.mod[[pc]]$pdes # assumes that des vars come before cat vars, which I think we do internally.
    int.use.l<-which(pc.mod[[pc]]$vars.cat[mcmc.mod.use,m,]==l.cat)
    if(length(int.use.l)==0)
      return(1)
    lD1<-pc.mod[[pc]]$sub.size[mcmc.mod.use,m,int.use.l]
    nlevels<-pc.mod[[pc]]$nlevels[l.cat]
    return(lD1/nlevels)
  }
}






C2Basis<-function(prior,pc.mod,l,m1,m2,pc1,pc2,mcmc.mod.use1,mcmc.mod.use2){
  if(l<=pc.mod[[pc1]]$pdes){ # could do pc1 or pc2, they have the same vars
    int.use.l1<-which(pc.mod[[pc1]]$vars.des[mcmc.mod.use1,m1,]==l)
    int.use.l2<-which(pc.mod[[pc2]]$vars.des[mcmc.mod.use2,m2,]==l)
    if(length(int.use.l1)==0 & length(int.use.l2)==0)
      return(1)
    if(length(int.use.l1)==0)
      return(C1Basis(prior,pc.mod,l,m2,pc2,mcmc.mod.use2))
    if(length(int.use.l2)==0)
      return(C1Basis(prior,pc.mod,l,m1,pc1,mcmc.mod.use1))

    #if(pc1==pc2 & m1==m2)
    #  return(C1Basis(prior,pc.mod,l,m1,pc1,mcmc.mod.use1)^2) ## is this right??

    q<-pc.mod[[pc1]]$degree
    s1<-pc.mod[[pc1]]$signs[mcmc.mod.use1,m1,int.use.l1]
    s2<-pc.mod[[pc2]]$signs[mcmc.mod.use2,m2,int.use.l2]
    t.ind1<-pc.mod[[pc1]]$knotInd.des[mcmc.mod.use1,m1,int.use.l1]
    t.ind2<-pc.mod[[pc2]]$knotInd.des[mcmc.mod.use2,m2,int.use.l2]
    t1<-pc.mod[[pc1]]$xx.des[t.ind1,l]
    t2<-pc.mod[[pc2]]$xx.des[t.ind2,l]



    if(t2<t1){
      temp<-t1
      t1<-t2
      t2<-temp
      temp<-s1
      s1<-s2
      s2<-temp
    }
    #browser()
    return(C22Basis(prior[[l]],t1,t2,s1,s2,q,m1,m2,pc1,pc2))
  } else{
    l.cat<-l-pc.mod[[pc1]]$pdes

    int.use.l1<-which(pc.mod[[pc1]]$vars.cat[mcmc.mod.use1,m1,]==l.cat)
    int.use.l2<-which(pc.mod[[pc2]]$vars.cat[mcmc.mod.use2,m2,]==l.cat)

    if(length(int.use.l1)==0 & length(int.use.l2)==0)
      return(1)
    if(length(int.use.l1)==0)
      return(C1Basis(prior,pc.mod,l,m2,pc2,mcmc.mod.use2))
    if(length(int.use.l2)==0)
      return(C1Basis(prior,pc.mod,l,m1,pc1,mcmc.mod.use1))

    #browser()
    sub1<-pc.mod[[pc1]]$sub.list[[mcmc.mod.use1]][[m1]][[int.use.l1]]
    sub2<-pc.mod[[pc2]]$sub.list[[mcmc.mod.use2]][[m2]][[int.use.l2]]
    if(is.na(sub1[1]) & is.na(sub2[1]))
      browser()
    nlevels<-pc.mod[[pc1]]$nlevels[l.cat]
    return(length(intersect(sub1,sub2))/nlevels)
  }
}


C22Basis<-function(prior,t1,t2,s1,s2,q,m1,m2,pc1,pc2){ # t1<t2
  cc<-const(signs=c(s1,s2),knots=c(t1,t2),degree=q)
  if((s1*s2)==0){
    return(0)
  }
  # if(m1==m2 & pc1==pc2){ #t1=t2, s1=s2 - NOT TRUE, since these could be different eof models
  #   return(1/(2*q+1)*((s1+1)/2-s1*t1)^(2*q+1)/cc)
  #   intabq1(prior[[l]],a,b,t,q)/cc
  #   if(s1==1){
  #
  #   } else{
  #
  #   }
  # } else{
    if(s1==1){
      if(s2==1){
        return(intabq2(prior,t2,1,t1,t2,q)/cc)
      } else{
        return(intabq2(prior,t1,t2,t1,t2,q)*(-1)^q/cc)
      }
    } else{
      if(s2==1){
        return(0)
      } else{
        return(intabq2(prior,0,t1,t1,t2,q)/cc)
      }
    }
  #}
}







get.f0<-function(prior,pc.mod,pc,mcmc.use){ # mcmc.mod.use is mcmc index not model index
  mcmc.mod.use<-pc.mod[[pc]]$model.lookup[mcmc.use]
  out<-pc.mod[[pc]]$beta[mcmc.use,1] # intercept
  if(pc.mod[[pc]]$nbasis[mcmc.use] > 0){
    for(m in 1:pc.mod[[pc]]$nbasis[mcmc.use]){
      out1<-pc.mod[[pc]]$beta[mcmc.use,1+m]
      for(l in 1:pc.mod[[pc]]$p){
        out1<-out1*C1Basis(prior,pc.mod,l,m,pc,mcmc.mod.use)
      }
      out<-out+out1
    }
  }
  return(out)
}



##################################################################################################################################################################
##################################################################################################################################################################
## modularized calibration
rmnorm<-function(mu, S){
  mu+c(rnorm(length(mu))%*%chol(S))
}
calibrate.bassBasis<-function(mod,y,a,b,nmcmc,verbose=T){ # assumes inputs to mod are standardized to (0,1), equal variance for all y values (should change to sim covariance)
  p<-ncol(mod$dat$xx)
  ny<-length(y)
  ns<-mod$mod.list[[1]]$nmcmc-mod$mod.list[[1]]$nburn # number of emu mcmc samples

  theta<-matrix(nrow=nmcmc,ncol=p)
  s2<-rep(NA,nmcmc)

#browser()
  theta[1,]<-.5
  pred.curr<-predict(mod,theta[1,,drop=F],mcmc.use=sample(ns,size=1),trunc.error=F)
  s2[1]<-1/rgamma(1,ny/2+a,b+sum((y-pred.curr)^2))

  eps<-1e-10
  cc<-2.4^2/p
  S<-diag(p)*eps
  count<-0
  for(i in 2:nmcmc){
    s2[i]<-1/rgamma(1,ny/2+a,b+sum((y-pred.curr)^2))

    theta[i,]<-theta[i-1,]
    if(i>300){
      mi<-1#max(1,i-300)
      S<-cov(theta[mi:(i-1),])*cc+diag(eps*cc,p)
    }
    theta.cand<-rmnorm(theta[i-1,],S)
    if(any(theta.cand<0 | theta.cand>1))
      alpha<- -9999
    else{
      pred.cand<-predict(mod,t(theta.cand),mcmc.use=sample(ns,size=1),trunc.error=F)
      alpha<- -.5/s2[i]*(sum((y-pred.cand)^2)-sum((y-pred.curr)^2))
    }
    if(log(runif(1))<alpha){
      theta[i,]<-theta.cand
      count<-count+1
    }

    pred.curr<-predict(mod,theta[i,,drop=F],mcmc.use=sample(ns,size=1),trunc.error=F)

    if(verbose & i%%100==0){
      pr<-c('MCMC iteration',i,myTimestamp(),'count:',count)
      cat(pr,'\n')
    }
  }

  return(list(theta=theta,s2=s2,count=count))
}



calibrateIndep.bassBasis<-function(mod,y,a,b,nmcmc,verbose=T){ # assumes inputs to mod are standardized to (0,1), equal variance for all y values (should change to sim covariance)
  p<-ncol(mod$dat$xx)
  ny<-length(y)
  ns<-mod$mod.list[[1]]$nmcmc-mod$mod.list[[1]]$nburn # number of emu mcmc samples

  theta<-matrix(nrow=nmcmc,ncol=p)
  s2<-rep(NA,nmcmc)

  #browser()
  theta[1,]<-.5
  pred.curr<-predict(mod,theta[1,,drop=F],mcmc.use=sample(ns,size=1),trunc.error=F)
  s2[1]<-1/rgamma(1,ny/2+a,b+sum((y-pred.curr)^2))

  count<-rep(0,p)
  for(i in 2:nmcmc){
    s2[i]<-1/rgamma(1,ny/2+a,b+sum((y-pred.curr)^2))

    theta[i,]<-theta[i-1,]

    for(j in 1:p){
      theta.cand<-theta[i,]
      theta.cand[j]<-runif(1)
      pred.cand<-predict(mod,t(theta.cand),mcmc.use=sample(ns,size=1),trunc.error=F)
      alpha<- -.5/s2[i]*(sum((y-pred.cand)^2)-sum((y-pred.curr)^2))
      if(log(runif(1))<alpha){
        theta[i,]<-theta.cand
        count[j]<-count[j]+1
      }
    }

    pred.curr<-predict(mod,theta[i,,drop=F],mcmc.use=sample(ns,size=1),trunc.error=F)

    if(verbose & i%%100==0){
      pr<-c('MCMC iteration',i,myTimestamp(),'count:',count)
      cat(pr,'\n')
    }
  }

  return(list(theta=theta,s2=s2,count=count))
}



##################################################################################################################################################################
##################################################################################################################################################################



plot.prior<-function(prior,plot=TRUE,n=1000,...){
  xx<-seq(prior$trunc[1],prior$trunc[2],length.out=n)
  if(prior$dist=='uniform'){
    out<-dunif(xx,prior$trunc[1],prior$trunc[2])
    z<-1
  }
  if(prior$dist=='normal'){
    out<-0
    z<-0
    for(i in 1:length(prior$weights)){
      zi<-pnorm(prior$trunc[2],prior$mean[i],prior$sd[i]) - pnorm(prior$trunc[1],prior$mean[i],prior$sd[i])
      z<-z+zi*prior$weights[i]
      out<-out+prior$weights[i]*dnorm(xx,prior$mean[i],prior$sd[i])
    }
  }
  if(prior$dist=='student'){
    out<-0
    z<-0
    for(i in 1:length(prior$weights)){
      zi<-pt((prior$trunc[2]-prior$mean[i])/prior$sd[i],prior$df[i]) - pt((prior$trunc[1]-prior$mean[i])/prior$sd[i],prior$df[i])
      z<-z+zi*prior$weights[i]
      out<-out+prior$weights[i]*(dt((xx-prior$mean[i])/prior$sd[i],prior$df[i])/prior$sd[i])
    }
  }
  if(plot)
    plot(xx,out/z,...)
  return(cbind(xx,out/z))
}






sample.prior<-function(prior,n){
  p<-length(prior)
  out<-matrix(nrow=n,ncol=p)
  for(i in 1:p){
    if(prior[[i]]$dist=='uniform'){
      out[,i]<-runif(n,prior[[i]]$trunc[1],prior[[i]]$trunc[2])
    } else{
    ncomp<-length(prior[[i]]$weights)
    comp<-sample(1:ncomp,size=n,prob=prior[[i]]$weights,replace=T)
    if(prior[[i]]$dist=='normal')
      #out[,i]<-rnorm(n,prior[[i]]$mean[comp],prior[[i]]$sd[comp])
      out[,i]<-suppressWarnings(truncdist::rtrunc(n,spec='norm',a=(prior[[i]]$trunc[1]-prior[[i]]$mean[comp])/prior[[i]]$sd[comp],b=(prior[[i]]$trunc[2]-prior[[i]]$mean[comp])/prior[[i]]$sd[comp])*prior[[i]]$sd[comp]+prior[[i]]$mean[comp])
    if(prior[[i]]$dist=='student')
      out[,i]<-truncdist::rtrunc(n,spec='t',df=prior[[i]]$df[comp],a=(prior[[i]]$trunc[1]-prior[[i]]$mean[comp])/prior[[i]]$sd[comp],b=(prior[[i]]$trunc[2]-prior[[i]]$mean[comp])/prior[[i]]$sd[comp])*prior[[i]]$sd[comp]+prior[[i]]$mean[comp]
    }
  }
  out
}


