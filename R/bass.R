#######################################################
# Author: Devin Francom, Los Alamos National Laboratory
# Protected under GPL-3 license
# Los Alamos Computer Code release C19031
# github.com/lanl/BASS
#######################################################

########################################################################
## main BASS function
########################################################################

#' @title Bayesian Adaptive Spline Surfaces (BASS)
#'
#' @description Fits a BASS model using RJMCMC.  Optionally uses parallel tempering to improve mixing.  Can be used with scalar or functional response.  Also can use categorical inputs.
#' @param xx a data frame or matrix of predictors.  Categorical predictors should be included as factors.
#' @param y a response vector (scalar response) or matrix (functional response).  Note: If \code{sum(y^2)} is large (i.e. \code{1e10}), please center/rescale (and rescale \code{g1} and \code{g2} if necessary).
#' @param maxInt integer for maximum degree of interaction in spline basis functions.  Defaults to the number of predictors, which could result in overfitting.
#' @param maxInt.func (functional response only) integer for maximum degree of interaction in spline basis functions describing the functional response.
#' @param maxInt.cat (categorical input only) integer for maximum degree of interaction of categorical inputs.
#' @param xx.func a vector, matrix or data frame of functional variables.
#' @param degree degree of splines.  Stability should be examined for anything other than 1.
#' @param maxBasis maximum number of basis functions.  This should probably only be altered if you run out of memory.
#' @param npart minimum number of non-zero points in a basis function.  If the response is functional, this refers only to the portion of the basis function coming from the non-functional predictors. Defaults to 20 or 0.1 times the number of observations, whichever is smaller.
#' @param npart.func same as npart, but for functional portion of basis function.
#' @param nmcmc number of RJMCMC iterations.
#' @param nburn number of the \code{nmcmc} iterations to disregard.
#' @param thin keep every \code{thin} samples
#' @param g1 shape for IG prior on \eqn{\sigma^2}.
#' @param g2 scale for IG prior on \eqn{\sigma^2}.
#' @param s2.lower lower bound for s2. Turns IG prior for s2 into a truncated IG.
#' @param h1 shape for gamma prior on \eqn{\lambda}.
#' @param h2 rate for gamma prior on \eqn{\lambda}.  This is the primary way to control overfitting.  A large value of \code{h2} favors fewer basis functions.
#' @param a.tau shape for gamma prior on \eqn{\tau}.
#' @param b.tau rate for gamma prior on \eqn{\tau}. Defaults to one over the number of observations, which centers the prior for the basis function weights on the unit information prior.
#' @param w1 nominal weight for degree of interaction, used in generating candidate basis functions.  Should be greater than 0.
#' @param w2 nominal weight for variables, used in generating candidate basis functions.  Should be greater than 0.
#' @param beta.prior what type of prior to use for basis coefficients, "g" or "jeffreys"
#' @param temp.ladder temperature ladder used for parallel tempering.  The first value should be 1 and the values should increase.
#' @param start.temper when to start tempering (after how many MCMC iterations). Defaults to 1000 or half of burn-in, whichever is smaller.
#' @param curr.list list of starting models (one element for each temperature), could be output from a previous run under the same model setup.
#' @param save.yhat logical; should predictions of training data be saved?
#' @param small logical; if true, returns a smaller object by leaving out \code{curr.list} and other unnecessary objects.  Use in combination with \code{save.yhat} to get smaller memory footprint for very large models.
#' @param verbose logical; should progress be displayed?
#' @param ret.str logical; return data and prior structures
#' @details Explores BASS model space by RJMCMC.  The BASS model has \deqn{y = f(x) + \epsilon,  ~~\epsilon \sim N(0,\sigma^2)} \deqn{f(x) = a_0 + \sum_{m=1}^M a_m B_m(x)} and \eqn{B_m(x)} is a BASS basis function (tensor product of spline basis functions). We use priors \deqn{a \sim N(0,\sigma^2/\tau (B'B)^{-1})} \deqn{M \sim Poisson(\lambda)} as well as the priors mentioned in the arguments above.
#' @return An object of class 'bass'.  The other output will only be useful to the advanced user.  Rather, users may be interested in prediction and sensitivity analysis, which are obtained by passing the entire object to the predict.bass or sobol functions.
#' @keywords nonparametric regression splines functional data analysis
#' @seealso \link{predict.bass} for prediction and \link{sobol} for sensitivity analysis.
#' @export
#' @import stats
#' @import utils
#' @example inst/examples.R
#'
bass<-function(xx,y,maxInt=3,maxInt.func=3,maxInt.cat=3,xx.func=NULL,degree=1,maxBasis=1000,npart=NULL,npart.func=NULL,nmcmc=10000,nburn=9000,thin=1,g1=0,g2=0,s2.lower=0,h1=10,h2=10,a.tau=.5,b.tau=NULL,w1=5,w2=5,beta.prior='g',temp.ladder=NULL,start.temper=NULL,curr.list=NULL,save.yhat=TRUE,small=FALSE,verbose=TRUE,ret.str=F){

  cl<-match.call()
  ########################################################################
  ## setup

  ## check inputs

  if(!posInt(maxInt))
    stop('invalid maxInt')
  if(!posInt(maxInt.func))
    stop('invalid maxInt.func')
  if(!posInt(maxInt.cat))
    stop('invalid maxInt.cat')
  #if(!posInt(degree))
  #  stop('invalid degree')
  if(!is.null(npart)){
    if(!posInt(npart))
      stop('invalid npart')
  }
  if(!is.null(npart.func)){
    if(!posInt(npart.func))
      stop('invalid npart.func')
  }
  if(!is.null(start.temper)){
    if(!posInt(start.temper))
      stop('invalid start.temper')
  }
  if(!posInt(nmcmc))
    stop('invalid nmcmc')
  #if(!posInt(nburn))
  #  stop('invalid nburn')
  if(!posInt(thin))
    stop('invalid thin')
  if(nburn>=nmcmc)
    stop('nmcmc must be greater than nburn')
  if(thin>(nmcmc-nburn))
    stop('combination of thin, nmcmc and nburn results in no samples')
  if(any(c(g1,g2)<0))
    stop('g1 and g2 must be greater than or equal to 0')
  if(any(c(h1,h2,a.tau,b.tau,w1,w2)<=0))
    stop('h1,h2,a.tau,b.tau,w1,w2 must be greater than 0')
  if(s2.lower<0)
    stop('s2.lower must be >= 0')

  ## process data
  if(any(is.na(xx)) | any(is.na(y)))
    stop('Current version does not allow missing data')

  y<-as.matrix(y)
  xx<-as.data.frame(xx)
  dx<-dim(xx)
  dxf<-dim(xx.func)
  dy<-dim(y)
  if(any(dy==1))
    y<-c(y)
  dy<-dim(y)

  if(is.null(dy)){
    func<-F
    pfunc<-0
    if(!is.null(xx.func))
      warning('xx.func ignored because there is no functional variable')
    if(length(y)!=dx[1])
      stop('dimension mismatch between xx and y')
  } else {
    func<-T
    if(is.null(xx.func))
      stop('missing xx.func')
    xx.func<-as.matrix(xx.func)
    dxf<-dim(xx.func)
    if(dy[1]!=dx[1]){
      y<-t(y)
      dy<-dim(y)
    }
    if(dy[1]!=dx[1])
      stop('dimension mismatch between xx and y')
    if(dy[2]!=dxf[1])
      xx.func<-t(xx.func)
    dxf<-dim(xx.func)
    if(dy[2]!=dxf[1])
      stop('dimension mismatch between xx.func and y')
    pfunc<-dxf[2]
    range.func<-apply(xx.func,2,range)
    xx.func<-apply(xx.func,2,scale_range)
  }

  if(func){
    if(dx[1]==dxf[1]) # this is dangerous because we would automatically correct it if we could tell it was wrong, but can't tell if it is wrong here
      warning('Possible dimension problem: make sure rows of y correspond to functional data')
  }

  des<-T
  cx<-sapply(xx,class)
  cx.factor<- cx == 'factor'
  if(any(cx.factor)){
    cat<-T
    if(all(cx.factor))
      des<-F
    xx.des<-as.matrix(xx[,!cx.factor,drop=F])
    xx.cat<-xx[,cx.factor,drop=F]
  } else{
    cat<-F
    xx.des<-as.matrix(xx)
    xx.cat<-NULL
  }
  if(des){
    range.des<-apply(xx.des,2,range)
    xx.des<-apply(xx.des,2,scale_range)
  }
  des.vars<-which(!cx.factor)
  cat.vars<-which(cx.factor)
  pdes<-length(des.vars)
  pcat<-length(cat.vars)

  type<-''
  if(des)
    type<-paste(type,'des',sep='_')
  if(cat)
    type<-paste(type,'cat',sep='_')
  if(func)
    type<-paste(type,'func',sep='_')




  # so cases are des, cat, des_cat, des_func, cat_func, des_cat_func

  ## handle tempering arguements

  if(is.null(temp.ladder)){
    temp.ladder<-1
  }
  if(max(temp.ladder)>(dx[1]/2)){
    temp.ladder<-temp.ladder[temp.ladder<(dx[1]/2)]
    if(length(temp.ladder)==0)
      stop('invalid temp.ladder (temperatures too high)')
  }
  if(min(temp.ladder)!=1)
    warning('min(temp.ladder) should equal 1')
  ntemps<-length(temp.ladder)
  if(ntemps==1){
    start.temper<-nmcmc
  }
  if(is.null(start.temper))
    start.temper<-min(1000,ceiling(nburn*.5))
  temp.val<-matrix(nrow=nmcmc,ncol=ntemps)

  if(any(temp.ladder<=0))
    stop('temp.ladder must be greater than 0 (should be greater than 1)')
  if(any(temp.ladder<1))
    warning('temp.ladder should be greater than 1')

  ## make a data object

  data<-list()
  data$y<-y
  if(des){
    data$xxt.des<-t(xx.des)
    data$vars.len.des<-NA
    data$xxt.des.unique<-list()
    data$unique.ind.des<-list()
    for(i in 1:pdes){
      data$xxt.des.unique[[i]]<-unique(data$xxt.des[i,])
      data$unique.ind.des[[i]]<-which(!duplicated(data$xxt.des[i,])) # gets the first instance of each unique value
      data$vars.len.des[i]<-length(data$xxt.des.unique[[i]])
    }
  }
  if(func){
    data$xxt.func<-t(xx.func)
    data$vars.len.func<-NA
    data$xxt.func.unique<-list()
    data$unique.ind.func<-list()
    for(i in 1:pfunc){
      data$xxt.func.unique[[i]]<-unique(data$xxt.func[i,])
      data$unique.ind.func[[i]]<-which(!duplicated(data$xxt.func[i,])) # gets the first instance of each unique value
      data$vars.len.func[i]<-length(data$xxt.func.unique[[i]])
    }
  }
  if(cat){
    data$levels<-lapply(xx.cat,levels)
    data$nlevels<-sapply(data$levels,length)
    data$xx.cat<-xx.cat
  }
  data$pdes<-pdes
  data$pfunc<-pfunc
  data$pcat<-pcat
  data$p<-dx[2]
  data$ndes<-dx[1]
  data$nfunc<-dxf[1]
  data$des<-des
  data$func<-func
  data$cat<-cat

  data$n<-prod(data$ndes,data$nfunc)
  data$ssy<-sum(data$y^2)
  data$death.prob.next<-1/3
  data$birth.prob<-1/3
  data$birth.prob.last<-1/3
  data$death.prob<-1/3
  data$itemp.ladder<-1/temp.ladder


  ## make a prior object

  npart.des<-npart
  if(is.null(npart.des)){
    npart.des<-min(20,.1*data$ndes)
  }
  if(is.null(npart.func) & func){
    npart.func<-min(20,.1*data$nfunc)
  }

  maxBasis<-min(maxBasis,data$n) # can't have more basis functions than data points
  maxInt.des<-min(maxInt,pdes) # can't have more interactions than variables
  maxInt.cat<-min(maxInt.cat,pcat)
  maxInt.func<-min(maxInt.func,pfunc)

  # prior object
  prior<-list()
  prior$maxInt.des<-maxInt.des
  prior$maxInt.cat<-maxInt.cat
  prior$maxInt.func<-maxInt.func
  prior$q<-degree
  prior$npart.des<-npart.des
  prior$npart.func<-npart.func
  prior$h1<-h1
  prior$h2<-h2
  prior$g1<-g1
  prior$g2<-g2
  prior$s2.lower<-s2.lower
  prior$a.beta.prec<-a.tau
  if(is.null(b.tau)){
    prior$b.beta.prec<-2/data$n
  } else{
  prior$b.beta.prec<-b.tau
  }
  prior$maxBasis<-maxBasis
  prior$minInt<-0
  if(des+cat+func==1) # if there is only one part, can't have minInt of 0
    prior$minInt<-1
  prior$miC<-abs(prior$minInt-1)
  prior$beta.gprior.ind<-as.numeric(beta.prior=='g')
  prior$beta.jprior.ind<-as.numeric(beta.prior=='jeffreys')


  ## make an object to store current MCMC state (one for each temperature)

  if(is.null(curr.list)){
    curr.list<-list()
    for(i in 1:ntemps){
      curr.list[[i]]<-list()

      if(des){
        curr.list[[i]]$I.star.des<-rep(w1,prior$maxInt.des+prior$miC)
        curr.list[[i]]$I.vec.des<-curr.list[[i]]$I.star.des/sum(curr.list[[i]]$I.star.des)
        curr.list[[i]]$z.star.des<-rep(w2,data$pdes)
        curr.list[[i]]$z.vec.des<-curr.list[[i]]$z.star.des/sum(curr.list[[i]]$z.star.des)
        curr.list[[i]]$des.basis<-matrix(rep(1,data$ndes))
      }
      if(cat){
        curr.list[[i]]$I.star.cat<-rep(w1,prior$maxInt.cat+prior$miC)
        curr.list[[i]]$I.vec.cat<-curr.list[[i]]$I.star.cat/sum(curr.list[[i]]$I.star.cat)
        curr.list[[i]]$z.star.cat<-rep(w2,data$pcat)
        curr.list[[i]]$z.vec.cat<-curr.list[[i]]$z.star.cat/sum(curr.list[[i]]$z.star.cat)
        curr.list[[i]]$cat.basis<-matrix(rep(1,data$ndes))
      }
      if(func){
        curr.list[[i]]$I.star.func<-rep(w1,prior$maxInt.func+prior$miC)
        curr.list[[i]]$I.vec.func<-curr.list[[i]]$I.star.func/sum(curr.list[[i]]$I.star.func)
        curr.list[[i]]$z.star.func<-rep(w2,data$pfunc)
        curr.list[[i]]$z.vec.func<-curr.list[[i]]$z.star.func/sum(curr.list[[i]]$z.star.func)
        curr.list[[i]]$func.basis<-matrix(rep(1,data$nfunc))
      }

      if(des & cat)
        curr.list[[i]]$dc.basis<-curr.list[[i]]$des.basis*curr.list[[i]]$cat.basis

      curr.list[[i]]$s2<-1
      curr.list[[i]]$lam<-1
      curr.list[[i]]$beta.prec<-1*prior$beta.gprior.ind
      curr.list[[i]]$nbasis<-0
      curr.list[[i]]$nc<-1

      curr.list[[i]]$knots.des<-matrix(numeric(0),ncol=maxInt.des)
      curr.list[[i]]$knotInd.des<-matrix(integer(0),ncol=maxInt.des)
      curr.list[[i]]$signs.des<-matrix(integer(0),ncol=maxInt.des)
      curr.list[[i]]$vars.des<-matrix(integer(0),ncol=maxInt.des)
      curr.list[[i]]$n.int.des<-0

      curr.list[[i]]$sub.list<-list()
      curr.list[[i]]$sub.size<-matrix(integer(0),ncol=maxInt.cat)
      curr.list[[i]]$vars.cat<-matrix(integer(0),ncol=maxInt.cat)
      curr.list[[i]]$n.int.cat<-0

      curr.list[[i]]$knots.func<-matrix(numeric(0),ncol=maxInt.func)
      curr.list[[i]]$knotInd.func<-matrix(integer(0),ncol=maxInt.func)
      curr.list[[i]]$signs.func<-matrix(integer(0),ncol=maxInt.func)
      curr.list[[i]]$vars.func<-matrix(integer(0),ncol=maxInt.func)
      curr.list[[i]]$n.int.func<-0

      curr.list[[i]]$Xty<-rep(NA,maxBasis+2)
      curr.list[[i]]$Xty[1]<-sum(data$y)
      curr.list[[i]]$XtX<-matrix(NA,nrow=maxBasis+2,ncol=maxBasis+2)
      curr.list[[i]]$XtX[1,1]<-data$n
      curr.list[[i]]$R<-chol(curr.list[[i]]$XtX[1,1])
      curr.list[[i]]$R.inv.t<-t(solve(curr.list[[i]]$R))
      curr.list[[i]]$bhat<-mean(data$y)
      curr.list[[i]]$qf<-crossprod(curr.list[[i]]$R%*%curr.list[[i]]$bhat)
      curr.list[[i]]$count<-rep(0,3)
      curr.list[[i]]$cmod<-F
      curr.list[[i]]$step<-NA
      curr.list[[i]]$temp.ind<-i
      curr.list[[i]]$type<-type
    }
  }

  # define functions according to type.  Doing eval parse every time the functions are used is slow.
  funcs<-list()
  funcs$birth<-eval(parse(text=paste('birth',type,sep='')))
  funcs$death<-eval(parse(text=paste('death',type,sep='')))
  funcs$change<-eval(parse(text=paste('change',type,sep='')))
  funcs$getYhat<-eval(parse(text=paste('getYhat',type,sep='')))


  ## prepare storage objects for mcmc draws

  nmod.max<-(nmcmc-nburn)/thin # max number of models (models don't necessarily change every iteration)
  if(des){
    signs.des<-knotInd.des<-vars.des<-array(dim=c(nmod.max,maxBasis,maxInt.des)) # truncate when returning at end of function
    n.int.des<-matrix(nrow=nmod.max,ncol=maxBasis) # degree of interaction
  }
  if(cat){
    sub.list<-list() # this is big...
    sub.size<-vars.cat<-array(dim=c(nmod.max,maxBasis,maxInt.cat))
    n.int.cat<-matrix(nrow=nmod.max,ncol=maxBasis)
  }
  if(func){
    signs.func<-knotInd.func<-vars.func<-array(dim=c(nmod.max,maxBasis,maxInt.func))
  # arrays use less space, esp integer arrays
    n.int.func<-matrix(nrow=nmod.max,ncol=maxBasis)
  }

  beta<-matrix(nrow=nmod.max,ncol=maxBasis+1) # +1 for intercept, nmcmc-nburn instead of nmod because beta updates every iteration
  nbasis<-s2<-lam<-beta.prec<-NA
  cmod<-F # indicator for whether we have changed models since last storing
  model.lookup<-NA # lookup table between models and mcmc iterations

  log.post.cold<-rep(NA,nmcmc) # log posterior for cold chain (the one we care about)

  if(save.yhat){
    yhat.sum<-0 # if we don't want to store all yhat draws, can still get running average
    yhat<-array(dim=c(nmod.max,data$ndes,data$nfunc))
  }

  # temperature index
  cold.chain<-1 # to start, the cold chain is curr.list[[1]]
  temp.ind<-1:ntemps # we will change this vector as we swap temperatures
  count.swap<-count.swap1000<-count.swap.prop<-rep(0,ntemps-1) # number of swaps between each set of neighbors
  swap<-NA # to keep track of swaps
  #require(parallel) # for tempering


  ########################################################################
  ## MCMC


  if(verbose)
    cat('MCMC Start',myTimestamp(),'nbasis:',curr.list[[cold.chain]]$nbasis,'\n')
  n.models<-keep.sample<-0 # indexes for storage
  for(i in 2:nmcmc){

    ## update model for each temperature - can be parallel

    curr.list<-lapply(curr.list,updateMCMC,prior=prior,data=data,funcs=funcs)
    #curr.list<-parLapply(cluster,curr.list,updateMCMC)
    #curr.list<-parallel::mclapply(curr.list,updateMCMC,prior=prior,data=data,funcs=funcs,mc.preschedule=T,mc.cores=1)
    #curr.list<-parLapplyLB(cl,curr.list,updateMCMC,prior=prior,data=data,funcs=funcs)
    # TODO: DO SOMETHING LIKE THIS BUT KEEP EVERYTHING SEPARATE ON THE CLUSTER, all we need is lpost, cmod - MPI

    ## parallel tempering swap

    # if(i>start.temper){# & (i%%20==0)){ #only start after a certain point, and only try every 20
    #   # sample temp.ind.swap from 1:(ntemps-1), then swap with temp.ind.swap+1
    #   temp.ind.swap1<-sample(1:(ntemps-1),size=1) # corresponds to temperature temp.ladder[temp.ind.swap1]
    #   temp.ind.swap2<-temp.ind.swap1+1 # always use the neighboring chain on the right
    #   chain.ind1<-which(temp.ind==temp.ind.swap1) # which chain has temperature temp.ladder[temp.ind.swap1]
    #   chain.ind2<-which(temp.ind==temp.ind.swap2)
    #   alpha.swap<-(data$itemp.ladder[temp.ind.swap1]-data$itemp.ladder[temp.ind.swap2])*(curr.list[[chain.ind2]]$lpost-curr.list[[chain.ind1]]$lpost)
    #   if(is.nan(alpha.swap) | is.na(alpha.swap)){
    #     alpha.swap<- -9999
    #     warning('large values of temp.ladder too large')
    #   }
    #   count.swap.prop[temp.ind.swap1]<-count.swap.prop[temp.ind.swap1]+1
    #   #browser()
    #   if(log(runif(1)) < alpha.swap){
    #     # swap temperatures
    #     temp.ind[chain.ind1]<-temp.ind.swap2
    #     temp.ind[chain.ind2]<-temp.ind.swap1
    #     curr.list[[chain.ind1]]$temp.ind<-temp.ind.swap2
    #     curr.list[[chain.ind2]]$temp.ind<-temp.ind.swap1
    #
    #     count.swap[temp.ind.swap1]<-count.swap[temp.ind.swap1]+1
    #     count.swap1000[temp.ind.swap1]<-count.swap1000[temp.ind.swap1]+1
    #     swap[i]<-temp.ind.swap1
    #     if(temp.ind.swap1==1){
    #       cmod<-T # we changed models
    #       cold.chain<-chain.ind2 #which(temp.ind==1)
    #     }
    #   }
    # }

    if(i>start.temper){# & (i%%20==0)){ #only start after a certain point, and only try every 20
      # sample temp.ind.swap from 1:(ntemps-1), then swap with temp.ind.swap+1
      for(dummy in 1:ntemps){
        ts<-sort(sample(1:ntemps,size=2))
        temp.ind.swap1<-ts[1]#sample(1:(ntemps-1),size=1) # corresponds to temperature temp.ladder[temp.ind.swap1]
        temp.ind.swap2<-ts[2]#temp.ind.swap1+1 # always use the neighboring chain on the right
        chain.ind1<-which(temp.ind==temp.ind.swap1) # which chain has temperature temp.ladder[temp.ind.swap1]
        chain.ind2<-which(temp.ind==temp.ind.swap2)
        alpha.swap<-(data$itemp.ladder[temp.ind.swap1]-data$itemp.ladder[temp.ind.swap2])*(curr.list[[chain.ind2]]$lpost-curr.list[[chain.ind1]]$lpost)
        if(is.nan(alpha.swap) | is.na(alpha.swap)){
          alpha.swap<- -9999
          warning('large values of temp.ladder too large')
        }
        count.swap.prop[temp.ind.swap1]<-count.swap.prop[temp.ind.swap1]+1
        #browser()
        if(log(runif(1)) < alpha.swap){
          # swap temperatures
          temp.ind[chain.ind1]<-temp.ind.swap2
          temp.ind[chain.ind2]<-temp.ind.swap1
          curr.list[[chain.ind1]]$temp.ind<-temp.ind.swap2
          curr.list[[chain.ind2]]$temp.ind<-temp.ind.swap1

          count.swap[temp.ind.swap1]<-count.swap[temp.ind.swap1]+1
          count.swap1000[temp.ind.swap1]<-count.swap1000[temp.ind.swap1]+1
          swap[i]<-temp.ind.swap1
          if(temp.ind.swap1==1){
            cmod<-T # we changed models
            cold.chain<-chain.ind2 #which(temp.ind==1)
          }
        }
      }
    }

    log.post.cold[i]<-curr.list[[cold.chain]]$lpost
    temp.val[i,]<-temp.ind



    ## write current model if past burnin and model is unique
    if((i>nburn) & (((i-nburn)%%thin)==0)){
      # these things are updated every time
      keep.sample<-keep.sample+1 # indexes samples
      nb<-curr.list[[cold.chain]]$nbasis
      nbasis[keep.sample]<-nb
      beta[keep.sample,1:(nb+1)]<-curr.list[[cold.chain]]$beta
      s2[keep.sample]<-curr.list[[cold.chain]]$s2
      lam[keep.sample]<-curr.list[[cold.chain]]$lam
      beta.prec[keep.sample]<-curr.list[[cold.chain]]$beta.prec
      if(save.yhat){
        yhat.current<-funcs$getYhat(curr.list[[cold.chain]],nb)
        if(func){
          yhat[keep.sample,,]<-yhat.current
        } else{
          yhat[keep.sample,]<-yhat.current
        }
        yhat.sum<-yhat.sum+yhat.current

      }
      # save cold chain basis parms if they are different from previous (cmod=T)
      if(cmod || curr.list[[cold.chain]]$cmod){ # can I actually get curr.list[[cold.chain]]$cmod easily from the core it is on?
        n.models<-n.models+1 # indexes models
        if(nb>0){
          if(des){
            vars.des[n.models,1:nb,]<-as.integer(curr.list[[cold.chain]]$vars.des)
            signs.des[n.models,1:nb,]<-as.integer(curr.list[[cold.chain]]$signs.des)
            knotInd.des[n.models,1:nb,]<-as.integer(curr.list[[cold.chain]]$knotInd.des)
            n.int.des[n.models,1:nb]<-as.integer(curr.list[[cold.chain]]$n.int.des)
          }
          if(cat){
            vars.cat[n.models,1:nb,]<-as.integer(curr.list[[cold.chain]]$vars.cat)
            sub.size[n.models,1:nb,]<-as.integer(curr.list[[cold.chain]]$sub.size)
            sub.list[[n.models]]<-curr.list[[cold.chain]]$sub.list
            n.int.cat[n.models,1:nb]<-as.integer(curr.list[[cold.chain]]$n.int.cat)
          }
          if(func){
            vars.func[n.models,1:nb,]<-as.integer(curr.list[[cold.chain]]$vars.func)
            signs.func[n.models,1:nb,]<-as.integer(curr.list[[cold.chain]]$signs.func)
            knotInd.func[n.models,1:nb,]<-as.integer(curr.list[[cold.chain]]$knotInd.func)
            n.int.func[n.models,1:nb]<-as.integer(curr.list[[cold.chain]]$n.int.func)
          }
        }
        cmod<-F # reset change model indicator after writing current model
        curr.list[[cold.chain]]$cmod<-F
      }
      model.lookup[keep.sample]<-n.models # update lookup table
    }

    #if(calibrate){
    #  theta[i,]<-sampleTheta(curr.list[[cold.chain]])
    #  delta[i,]<-sampleDelta(curr.list[[cold.chain]])
    #}

    if(verbose & i%%1000==0){
      pr<-c('MCMC iteration',i,myTimestamp(),'nbasis:',curr.list[[cold.chain]]$nbasis)
      if(i>start.temper)
        pr<-c(pr,'tempering acc',round(count.swap1000/1000*(ntemps-1),3)) # swap acceptance rate
        #pr<-c(pr,'tempering acc',round(count.swap/count.swap.prop,3)) # swap acceptance rate
      cat(pr,'\n')
      count.swap1000<-rep(0,ntemps-1)
    }

  }

  ########################################################################
  ## return

  out.yhat<-list()
  if(save.yhat){
    out.yhat<-list(yhat.mean=yhat.sum/nmod.max,yhat=yhat)
  }
  out.str<-list()
  if(ret.str)
    out.str<-list(data=data,prior=prior,funcs=funcs)

  out<-list(
       call=cl,
       beta=beta,
       s2=s2,
       lam=lam,
       nbasis=nbasis,
       degree=degree,
       nmcmc=nmcmc,
       nburn=nburn,
       thin=thin,
       p=data$p,
       beta.prec=beta.prec,
       y=y,
       log.post.cold=log.post.cold,
       swap=swap,
       count.swap=count.swap,
       count.swap.prop=count.swap.prop,
       temp.val=temp.val,
       temp.ladder=temp.ladder,
       n.models=n.models,
       model.lookup=model.lookup,
       des=des,func=func,cat=cat,type=type,cx=cx
  )
  if(!small){
    out$curr.list<-curr.list # for restarting
  }
  mb<-max(nbasis)


  out.des<-list()
  if(des){
    out.des<-list(
      knotInd.des=knotInd.des[1:n.models,1:mb,,drop=F],
      signs.des=signs.des[1:n.models,1:mb,,drop=F],
      vars.des=vars.des[1:n.models,1:mb,,drop=F],
      n.int.des=n.int.des[1:n.models,1:mb,drop=F],
      maxInt.des=maxInt.des,
      pdes=pdes,
      xx.des=xx.des,range.des=range.des,
      unique.ind.des=data$unique.ind.des
    )
    if(!small){
      out.des$des.basis<-curr.list[[cold.chain]]$des.basis
    }
  }

  out.cat<-list()
  if(cat){
    out.cat<-list(
      vars.cat=vars.cat[1:n.models,1:mb,,drop=F],
      sub.size=sub.size[1:n.models,1:mb,,drop=F],
      sub.list=sub.list,
      n.int.cat=n.int.cat[1:n.models,1:mb,drop=F],
      maxInt.cat=maxInt.cat,
      pcat=pcat,
      xx.cat=xx.cat,
      nlevels=data$nlevels
    )
    if(!small){
      out.cat$cat.basis<-curr.list[[cold.chain]]$cat.basis
    }
  }

  out.func<-list()
  if(func){
    out.func<-list(
      knotInd.func=knotInd.func[1:n.models,1:mb,,drop=F],
      signs.func=signs.func[1:n.models,1:mb,,drop=F],
      vars.func=vars.func[1:n.models,1:mb,,drop=F],
      n.int.func=n.int.func[1:n.models,1:mb,drop=F],
      maxInt.func=maxInt.func,
      pfunc=pfunc,
      xx.func=xx.func,range.func=range.func,
      unique.ind.func=data$unique.ind.func
    )
    if(!small){
      out.func$func.basis=curr.list[[cold.chain]]$func.basis
    }
  }

  #stopCluster(cluster)
  ret<-c(out.yhat,out,out.des,out.cat,out.func,out.str)
  class(ret)<-'bass'
  return(ret)
}
