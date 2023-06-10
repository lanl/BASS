#' @title Calibrate a bassPCA or bassBasis Model to Data
#'
#' @description Robust modular calibration of a bassPCA or bassBasis emulator using adaptive Metropolis, tempering, and decorrelation steps in an effort to be free of any user-required tuning.
#' @param y vector of calibration data
#' @param mod a emulator of class bassBasis, whose predictions should match y (i.e., predictions from mod should be the same length as y)
#' @param type one of c(1,2). 1 indicates a model that uses independent truncation error variance, no measurement error correlation, and discrepancy on a basis while type 2 indicates a model that uses a full truncation error covariance matrix, a full measurement error correlation matrix, a fixed full discrepancy covariance matrix, and a fixed discrepancy mean. 1 is for situations where computational efficiency is important (because y is dense), while 2 is only for cases where y is a short vector.
#' @param sd.est vector of prior estimates of measurement error standard deviation
#' @param s2.df vector of degrees of freedom for measurement error sd prior estimates
#' @param s2.ind index vector, same length as y, indicating which sd.est goes with which y
#' @param meas.error.cor a fixed correlation matrix for the measurement errors
#' @param bounds a 2xp matrix of bounds for each input parameter, where p is the number of input parameters.
#' @param discrep.mean discrepancy mean (fixed), only used if type=2
#' @param discrep.mat discrepancy covariance (fixed, for type 2) or basis (if not square, for type 1)
#' @param nmcmc number of MCMC iterations.
#' @param temperature.ladder an increasing vector, all greater than 1, for tempering. Geometric spacing is recommended, so that you have (1+delta)^(0:ntemps), where delta is small (typically between 0.05 and 0.2) and ntemps is the number of elements in the vector.
#' @param decor.step.every integer number of MCMC iterations between decorrelation steps.
#' @param verbose logical, whether to print progress.
#'
#' @details Fits a modular Bayesian calibration model, with \deqn{y = Kw(\theta) + Dv + \epsilon,  ~~\epsilon \sim N(0,\sigma^2 R)} \deqn{f(x) = a_0 + \sum_{m=1}^M a_m B_m(x)} and \eqn{B_m(x)} is a BASS basis function (tensor product of spline basis functions). We use priors \deqn{a \sim N(0,\sigma^2/\tau (B'B)^{-1})} \deqn{M \sim Poisson(\lambda)} as well as the priors mentioned in the arguments above.
#' @return An object
#' @seealso \link{predict.bassBasis} for prediction and \link{sobolBasis} for sensitivity analysis.
#' @export
#' @import utils
#' @example inst/examplesPCA.R
#'
calibrate.bassBasis<-function(y,
                          mod,
                          type,
                          sd.est,
                          s2.df,
                          s2.ind,
                          meas.error.cor,
                          bounds,
                          discrep.mean,
                          discrep.mat,
                          nmcmc=10000,
                          temperature.ladder=1.05^(0:30),
                          decor.step.every=100,
                          verbose=T
){

  tl<-temperature.ladder
  decor<-decor.step.every
  func<-function(mod,xx,ii)
    predict(mod,xx,mcmc.use=ii,nugget=F,trunc.error=F)
  vars<-do.call(cbind,lapply(mod$mod.list,function(a) a$s2))
  basis<-mod$dat$basis

  if(type==1)
    trunc.error.cov<-diag(diag(cov(t(mod$dat$trunc.error))))
  if(type==2)
    trunc.error.cov<-cov(t(mod$dat$trunc.error))

  rmnorm<-function(mu, S){
    mu+c(rnorm(length(mu))%*%chol(S))
  }
  if(any(s2.df==0)){
    ldig.kern<-function(x,a,b)
      -log(x+1)
  } else{
    ldig.kern<-function(x,a,b)
      -(a+1)*log(x)-b/x
  }
  unscale.range<-function(x,r){
    x*(r[2]-r[1])+r[1]
  }
  rigammaTemper<-function(n,shape,scale,itemper){
    1/rgamma(n,itemper*(shape+1)-1,rate=itemper*scale)
  }
  myTimestamp<-function(){
    x<-Sys.time()
    paste('#--',format(x,"%b %d %X"),'--#')
  }

  cor2cov<-function(R,S) # https://stats.stackexchange.com/questions/62850/obtaining-covariance-matrix-from-correlation-matrix
    outer(S,S) * R

  p<-ncol(bounds)
  a<-s2.df/2
  b<-a*sd.est^2
  ns2<-length(unique(s2.ind))
  n.s2.ind<-rep(0,ns2)
  for(j in 1:ns2)
    n.s2.ind[j]<-sum(s2.ind==j)

  s2.first.ind<-NA
  for(j in 1:ns2)
    s2.first.ind[j]<-which(s2.ind==j)[1]

  ny<-length(y)
  ntemps<-length(tl)

  class<-'mult'
  nd<-1
  if(ncol(discrep.mat)<ny){
    class<-'func'
    nd<-ncol(discrep.mat)
  }

  theta<-array(dim=c(nmcmc,ntemps,p))
  s2<-array(dim=c(nmcmc,ntemps,ns2))
  discrep.vars<-array(dim=c(nmcmc,ntemps,nd))
  #temp.ind<-matrix(nrow=nmcmc,ncol=ntemps)

  itl<-1/tl # inverse temperature ladder

  tran<-function(th){
    #th2<-pnorm(th)
    for(i in 1:ncol(th)){
      th[,i]<-unscale.range(th[,i],bounds[,i])
    }
    th
  }

  theta[1,,]<-runif(prod(dim(theta[1,,,drop=F])))
  s2[1,,]<-1
  discrep.vars[1,,]<-0

  my.solve<-function(x){
    u <- chol(x)
    chol2inv(u)
  }

  swm<-function(Ainv,U,Cinv,V){ # sherman woodbury morrison (A+UCV)^-1
    Ainv - Ainv %*% U %*% my.solve(Cinv + V %*% Ainv %*% U) %*% V %*% Ainv
  }
  swm.ldet<-function(Ainv,U,Cinv,V,Aldet,Cldet){ # sherman woodbury morrison |A+UCV|
    determinant(Cinv + V %*% Ainv %*% U, logarithm=T)$mod + Aldet + Cldet
  }

  curr<-list()
  for(t in 1:ntemps)
    curr[[t]]<-list()
  dat<-list(y=y,basis=basis,s2.ind=s2.ind)



  if(class=='mult'){



    dat$trunc.error.cov<-trunc.error.cov
    dat$meas.error.cor<-meas.error.cor
    dat$discrep.cov<-discrep.mat
    dat$discrep.mean<-discrep.mean

    for(t in 1:ntemps){
      curr[[t]]$s2<-s2[1,t,]
    }

    lik.cov.inv<-function(dat,curr){#trunc.error.cov,Sigma,discrep.mat,discrep.vars,basis,emu.vars){
      Sigma<-cor2cov(dat$meas.error.cor,sqrt(curr$s2[dat$s2.ind]))
      mat<-chol(dat$trunc.error.cov+Sigma+dat$discrep.cov+dat$basis%*%diag(curr$emu.vars,mod$dat$n.pc)%*%t(dat$basis))
      inv<-chol2inv(mat)
      ldet<-2*sum(log(diag(mat)))
      return(list(inv=inv, ldet=ldet))
    }


    llik<-function(dat,curr){#y,pred,discrep,cov.inv,ldet){ # use ldet=0 if it doesn't matter
      vec <- dat$y-curr$pred-dat$discrep.mean
      -.5*(curr$cov$ldet + t(vec)%*%curr$cov$inv%*%(vec))
    }

  }

  if(class=='func'){



    dat$trunc.error.var<-diag(trunc.error.cov) # assumed diagonal
    dat$D<-discrep.mat
    dat$discrep<-dat$D%*%discrep.vars[1,t,]
    dat$discrep.tau<-1
    dat$nd<-ncol(dat$D)

    for(t in 1:ntemps){
      curr[[t]]$s2<-s2[1,t,]
      #curr[[t]]$vars.sigma<-curr[[t]]$s2[s2.ind]
      curr[[t]]$discrep<-dat$D%*%discrep.vars[1,t,]
    }

    lik.cov.inv<-function(dat,curr){#trunc.error.cov,Sigma,discrep.mat,discrep.vars,basis,emu.vars){ # Sigma, trunc.error.cov are diagonal
      vec<-dat$trunc.error.var+curr$s2[dat$s2.ind]
      Ainv<-diag(1/vec)
      Aldet<-sum(log(vec))
      inv<-swm(Ainv,basis,diag(1/curr$emu.vars),t(basis))
      ldet<-swm.ldet(Ainv,basis,diag(1/curr$emu.vars),t(basis),Aldet,sum(log(curr$emu.vars)))
      return(list(inv=inv, ldet=ldet))
    }

    llik<-function(dat,curr){#y,pred,discrep,cov.inv,ldet){ # use ldet=0 if it doesn't matter
      vec <- dat$y-curr$pred-curr$discrep
      -.5*(curr$cov$ldet + t(vec)%*%curr$cov$inv%*%(vec))
    }

  }


  bigMat.curr<-Sigma.curr<-list()
  nmcmc.emu<-nrow(vars)

  ii<-NA
  ii[1]<-sample(nmcmc.emu,1)
  pred.curr<-matrix(func(mod,tran(matrix(theta[1,,],ncol=p)),ii[1]),nrow=ntemps)

  for(t in 1:ntemps){
    curr[[t]]$pred<-pred.curr[t,]
    curr[[t]]$emu.vars<-vars[ii[1],]
  }


  eps<-1e-13
  tau<-rep(-4,ntemps) # scaling
  tau.ls2<-rep(0,ntemps)
  cc<-2.4^2/p
  S<-mu<-cov<-S.ls2<-mu.ls2<-cov.ls2<-list()
  for(t in 1:ntemps){
    S[[t]]<-diag(p)*1e-6
    S.ls2[[t]]<-diag(ns2)*1e-6
  }
  count<-matrix(0,nrow=ntemps,ncol=ntemps)
  count.decor<-matrix(0,nrow=p,ncol=ntemps)
  count100<-count.s2<-rep(0,ntemps)


  for(i in 2:nmcmc){

    theta[i,,]<-theta[i-1,,] # current set at previous (update below)
    s2[i,,]<-s2[i-1,,]

    ########################################################
    ## update s2


    ii[i]<-sample(nmcmc.emu,1)
    pred.curr<-matrix(func(mod,tran(matrix(theta[i-1,,],ncol=p)),ii[i]),nrow=ntemps)
    for(t in 1:ntemps){
      curr[[t]]$pred<-pred.curr[t,]
      curr[[t]]$emu.vars<-vars[ii[i],]
      curr[[t]]$cov<-lik.cov.inv(dat,curr[[t]])
      curr[[t]]$llik<-llik(dat,curr[[t]])
    }


    if(i==300){ # start adapting
      for(t in 1:ntemps){
        if(ns2>1){
          mu.ls2[[t]]<-colMeans(log(s2[1:(i-1),t,]))
          cov.ls2[[t]]<-cov(log(s2[1:(i-1),t,]))
        } else{
          mu.ls2[[t]]<-matrix(mean(log(s2[1:(i-1),t,])))
          cov.ls2[[t]]<-matrix(var(log(s2[1:(i-1),t,])))
        }

        S.ls2[[t]]<- (cov.ls2[[t]]*cc+diag(eps*cc,ns2))*exp(tau.ls2[t])
      }
    }
    if(i>300){ # adaptation updates
      for(t in 1:ntemps){
        #browser()
        mu.ls2[[t]]<-mu.ls2[[t]]+(log(s2[(i-1),t,])-mu.ls2[[t]])/(i-1)
        cov.ls2[[t]]<-(i-2)/(i-1)*cov.ls2[[t]] + (i-2)/(i-1)^2*tcrossprod(log(s2[(i-1),t,])-mu.ls2[[t]])
        S.ls2[[t]]<-(cov.ls2[[t]]*cc+diag(eps*cc,ns2))*exp(tau.ls2[t])
      }
    }

    #ls2.cand<-matrix(nrow=ntemps,ncol=ns2)
    for(t in 1:ntemps){
      ls2.cand<-rmnorm(log(curr[[t]]$s2),S.ls2[[t]]) # generate candidate

      cand<-curr[[t]]
      cand$s2<-exp(ls2.cand)
      cand$cov<-lik.cov.inv(dat,cand)
      cand$llik<-llik(dat,cand)

      alpha<- itl[t]*(
        + cand$llik + sum(ldig.kern(cand$s2,a,b)) + sum(log(cand$s2))
        - curr[[t]]$llik - sum(ldig.kern(curr[[t]]$s2,a,b)) - sum(log(curr[[t]]$s2))
      ) # log lik + log prior + log jacobian

      #if(i>5000)
      #  browser()
      if(log(runif(1))<alpha){
        curr[[t]]<-cand
        s2[i,t,]<-cand$s2
        count.s2[t]<-count.s2[t]+1
      }

    }

    #if(i>5000)
    #  browser()

    ########################################################
    ## update discrep.vars with gibbs step
    #???~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    if(class=='func'){
      for(t in 1:ntemps){
        discrep.S<-my.solve(diag(dat$nd)/dat$discrep.tau + t(dat$D) %*% curr[[t]]$cov$inv %*% dat$D)
        discrep.m<-t(dat$D) %*% curr[[t]]$cov$inv %*% (dat$y-curr[[t]]$pred)
        discrep.vars[i,t,]<-c(rmnorm(discrep.S%*%discrep.m, discrep.S/itl[t]))
        curr[[t]]$discrep.vars<-discrep.vars[i,t,]
        curr[[t]]$discrep<-dat$D%*%discrep.vars[i,t,]
      }
    }


    ########################################################
    ## adaptive block update for theta within each temperature - covariance of previous samples scaled (by exp(tau)) based on acceptance rate for last 100 samples, since tempering makes large gaps

    if(i==300){ # start adapting
      for(t in 1:ntemps){
        mu[[t]]<-colMeans(theta[1:(i-1),t,])
        cov[[t]]<-cov(theta[1:(i-1),t,])
        S[[t]]<-(cov[[t]]*cc+diag(eps*cc,p))*exp(tau[t])
      }
    }
    if(i>300){ # adaptation updates
      for(t in 1:ntemps){
        #browser()
        mu[[t]]<-mu[[t]]+(theta[(i-1),t,]-mu[[t]])/(i-1)
        cov[[t]]<-(i-2)/(i-1)*cov[[t]] + (i-2)/(i-1)^2*tcrossprod(theta[(i-1),t,]-mu[[t]])
        S[[t]]<-(cov[[t]]*cc+diag(eps*cc,p))*exp(tau[t])
      }
    }


    theta.cand<-matrix(nrow=ntemps,ncol=p)
    for(t in 1:ntemps)
      theta.cand[t,]<-rmnorm(theta[i-1,t,],S[[t]]) # generate candidate for each temperature


    pred.cand<-matrix(func(mod,tran(matrix(theta.cand,ncol=p)),ii[i]),nrow=ntemps)

    # Dinv.cand is the same as Dinv.curr

    for(t in 1:ntemps){ # loop over temperatures, do a block MCMC update

      cand<-curr[[t]]
      cand$pred<-pred.cand[t,]
      cand$theta<-theta.cand[t,]


      if(any(cand$theta<0 | cand$theta>1))
        alpha<- -9999
      else{

        cand$llik<-llik(dat,cand)

        alpha<- itl[t]*(
          + cand$llik
          - curr[[t]]$llik
        ) # log lik + uniform prior

      }

      if(log(runif(1))<alpha){
        curr[[t]]<-cand
        theta[i,t,]<-theta.cand[t,]
        count[t,t]<-count[t,t]+1
        pred.curr[t,]<-pred.cand[t,]
        count100[t]<-count100[t]+1
      }
    }

    ########################################################
    ## acceptance-rate based adaptation of covariance scale
    if(i%%100==0){
      delta<-min(.1,1/sqrt(i)*5)
      for(t in 1:ntemps){
        if(count100[t]<23){
          tau[t]<-tau[t]-delta
        } else if(count100[t]>23){
          tau[t]<-tau[t]+delta
        }
      } # could be vectorized, but probably not expensive
      count100<-count100*0
    }


    ########################################################
    ## decorrelation step for theta, especially to decorrelate the thetas that dont change anything

    if(i%%decor==0){ # every so often, do a decorrelation step (use single-site independence sampler)
      for(k in 1:p){
        theta.cand<-theta[i,,,drop=F] # most up-to-date value
        theta.cand[1,,k]<-runif(ntemps) # independence sampler candidate (vectorize over temperatures)
        pred.cand<-matrix(func(mod,tran(matrix(theta.cand[1,,],ncol=p)),ii[i]),nrow=ntemps)

        for(t in 1:ntemps){
          cand<-curr[[t]]
          cand$pred<-pred.cand[t,]
          cand$theta<-theta.cand[1,t,]
          cand$llik<-llik(dat,cand)

          alpha<- itl[t]*(
            + cand$llik
            - curr[[t]]$llik
          ) # log lik + uniform prior

          if(log(runif(1))<alpha){
            curr[[t]]<-cand
            theta[i,t,k]<-theta.cand[1,t,k]
            count.decor[k,t]<-count.decor[k,t]+1
            pred.curr[t,]<-pred.cand[t,]
          }
        }
      }
    }





    ########################################################
    ## tempering swaps


    if(i>1000 & ntemps>1){ # tempering swap
      for(dummy in 1:ntemps){ # repeat tempering step a bunch of times
        sw<-sort(sample(ntemps,size=2))

        alpha<-(itl[sw[2]]-itl[sw[1]])*(
          curr[[sw[1]]]$llik + sum(ldig.kern(curr[[sw[1]]]$s2,a,b)) # plus discrep terms
          - curr[[sw[2]]]$llik - sum(ldig.kern(curr[[sw[2]]]$s2,a,b))
        )

        if(log(runif(1))<alpha){
          temp<-theta[i,sw[1],]
          theta[i,sw[1],]<-theta[i,sw[2],]
          theta[i,sw[2],]<-temp
          temp<-s2[i,sw[1],]
          s2[i,sw[1],]<-s2[i,sw[2],]
          s2[i,sw[2],]<-temp
          temp<-discrep.vars[i,sw[1],]
          discrep.vars[i,sw[1],]<-discrep.vars[i,sw[2],]
          discrep.vars[i,sw[2],]<-temp
          count[sw[1],sw[2]]<-count[sw[1],sw[2]]+1
          temp<-pred.curr[sw[1],] # not sampling posterior predictive each time, for speed
          pred.curr[sw[1],]<-pred.curr[sw[2],]
          pred.curr[sw[2],]<-temp
          temp<-curr[[sw[1]]]
          curr[[sw[1]]]<-curr[[sw[2]]]
          curr[[sw[2]]]<-temp
        }
      }
    }










    if(verbose & i%%100==0){
      pr<-c('MCMC iteration',i,myTimestamp(),'count:',diag(count))
      cat(pr,'\n')
    }
  }

  ##th2<-pnorm(theta)
  #for(ii in 1:p){
  #  theta[,,ii]<-unscale.range(theta[,,ii],bounds[,ii])
  #}
  return(list(theta=theta,s2=s2,count=count,count.decor=count.decor,tau=tau,ii=ii,curr=curr,dat=dat,discrep.vars=discrep.vars))
}




















#
#
#
# ##################################################################################################################################################################
# ##################################################################################################################################################################
# ## modularized calibration
# rmnorm<-function(mu, S){
#   mu+c(rnorm(length(mu))%*%chol(S))
# }
# calibrate<-function(mod,y,sd.est,s2.df,bounds,nmcmc,tl=1,verbose=T,decor=100,pred.ncores=1,pred.parType='fork',trunc.error=F){ # assumes inputs to mod are standardized to (0,1), equal variance for all y values (should change to sim covariance)
#   p<-ncol(bounds)
#   a<-s2.df/2
#   b<-a*sd.est^2
#
#   # could allow s2.ind vector, like in python code
#   # or could allow s2.mult for a functional setup
#
#   ny<-length(y)
#   ntemps<-length(tl)
#
#   theta<-array(dim=c(nmcmc,ntemps,p))
#   s2<-matrix(nrow=nmcmc,ncol=ntemps)
#   #temp.ind<-matrix(nrow=nmcmc,ncol=ntemps)
#
#   itl<-1/tl # inverse temperature ladder
#
#   tran<-function(th){
#     #th2<-pnorm(th)
#     for(i in 1:ncol(th)){
#       th[,i]<-unscale.range(th[,i],bounds[,i])
#     }
#     th
#   }
#
#   #browser()
#   theta[1,,]<-runif(prod(dim(theta[1,,,drop=F])))
#   # pred.curr<-matrix(
#   #   predict(mod,
#   #           tran(matrix(theta[1,,],ncol=p)),
#   #           mcmc.use=sample(ns,size=1),
#   #           trunc.error=trunc.error,
#   #           nugget = T,
#   #           n.cores=pred.ncores,
#   #           parType = pred.parType),
#   #   nrow=ntemps)
#
#   pred.curr<-matrix(mod(tran(matrix(theta[1,,],ncol=p))),nrow=ntemps)
#   s2[1,]<-rigammaTemper(ntemps, ny/2+a, b+colSums((t(pred.curr)-y)^2)/2, itl)
#
#
#
#   eps<-1e-13
#   tau<-rep(-4,ntemps) # scaling
#   cc<-2.4^2/p
#   S<-mu<-cov<-S.rev<-mu.rev<-cov.rev<-list()
#   for(t in 1:ntemps)
#     S[[t]]<-diag(p)*1e-6
#   count<-matrix(0,nrow=ntemps,ncol=ntemps)
#   count.decor<-matrix(0,nrow=p,ncol=ntemps)
#   count100<-rep(0,ntemps)
#
#
#   for(i in 2:nmcmc){
#
#     theta[i,,]<-theta[i-1,,] # current set at previous (update below)
#
#     ########################################################
#     ## adaptive block update for theta within each temperature - covariance of previous samples scaled (by exp(tau)) based on acceptance rate for last 100 samples, since tempering makes large gaps
#
#     if(i==300){ # start adapting
#       for(t in 1:ntemps){
#         mu[[t]]<-colMeans(theta[1:(i-1),t,])
#         cov[[t]]<-cov(theta[1:(i-1),t,])
#         S[[t]]<-cov(theta[1:(i-1),t,])*cc+diag(eps*cc,p)
#       }
#     }
#     if(i>300){ # adaptation updates
#       for(t in 1:ntemps){
#         mu[[t]]<-mu[[t]]+(theta[(i-1),t,]-mu[[t]])/(i-1)
#         cov[[t]]<-(i-2)/(i-1)*cov[[t]] + (i-2)/(i-1)^2*tcrossprod(theta[(i-1),t,]-mu[[t]])
#         S[[t]]<-(cov[[t]]*cc+diag(eps*cc,p))*exp(tau[t])
#       }
#     }
#
#     # if(i>300 & i<1000){ # start adapting
#     #   for(t in 1:ntemps){
#     #     mu[[t]]<-colMeans(theta[1:(i-1),t,])
#     #     cov[[t]]<-cov(theta[1:(i-1),t,])
#     #     S[[t]]<-cov[[t]]*cc+diag(eps*cc,p)
#     #   }
#     # }
#     # if(i>1000){ # radius-based adaptation
#     #   for(t in 1:ntemps){
#     #     #browser()
#     #     dist<-sqrt(colSums((t(theta[1:(i-2),t,])-theta[i-1,t,])^2))
#     #     use<-which(dist<=quantile(dist,.15))
#     #     mu[[t]]<-colMeans(theta[use,t,])
#     #     cov[[t]]<-cov(theta[use,t,])
#     #     S[[t]]<-cov[[t]]*cc+diag(eps*cc,p)
#     #   }
#     # }
#
#     theta.cand<-matrix(nrow=ntemps,ncol=p)
#     for(t in 1:ntemps)
#       theta.cand[t,]<-rmnorm(theta[i-1,t,],S[[t]]) # generate candidate for each temperature
#
#     # if(i>1000){ # for reversibility when radius-based adapting
#     #   for(t in 1:ntemps){
#     #     #browser()
#     #     dist<-sqrt(colSums((t(theta[1:(i-2),t,])-theta.cand[t,])^2))
#     #     use<-which(dist<=quantile(dist,.15))
#     #     mu.rev[[t]]<-colMeans(theta[use,t,])
#     #     cov.rev[[t]]<-cov(theta[use,t,])
#     #     S.rev[[t]]<-cov.rev[[t]]*cc+diag(eps*cc,p)
#     #   }
#     # }
#
#     # pred.cand<-matrix(
#     #   predict(mod,
#     #           tran(matrix(theta.cand,ncol=p)),
#     #           mcmc.use=sample(ns,size=1),
#     #           trunc.error=trunc.error,
#     #           nugget=T,
#     #           n.cores=pred.ncores,
#     #           parType = pred.parType),
#     #   nrow=ntemps) # get BASS prediction at each candidate (major speedup by vectorizing across temperatures)
#
#     pred.cand<-matrix(mod(tran(matrix(theta.cand,ncol=p))),nrow=ntemps)
#
#     for(t in 1:ntemps){ # loop over temperatures, do a block MCMC update
#       if(any(theta.cand[t,]<0 | theta.cand[t,]>1))
#         alpha<- -9999
#       else{
#         alpha<- (-.5/s2[i-1,t]*itl[t]*(sum((y-pred.cand[t,])^2)-sum((y-pred.curr[t,])^2)) # posterior
#                 #-mnormt::dmnorm(theta.cand[t,],theta[i-1,t,],S[[t]],log = T) + mnormt::dmnorm(theta[i-1,t,],theta.cand[t,],S.rev[[t]],log = T)  # proposal when using radius-based adapting
#         )
#       }
#       if(log(runif(1))<alpha){
#         theta[i,t,]<-theta.cand[t,]
#         count[t,t]<-count[t,t]+1
#         pred.curr[t,]<-pred.cand[t,]
#         count100[t]<-count100[t]+1
#       }
#     }
#
#     ########################################################
#     ## acceptance-rate based adaptation of covariance scale
#     if(i%%100==0){
#       delta<-min(.1,1/sqrt(i)*5)
#       for(t in 1:ntemps){
#         if(count100[t]<23){
#           tau[t]<-tau[t]-delta
#         } else if(count100[t]>23){
#           tau[t]<-tau[t]+delta
#         }
#       } # could be vectorized, but probably not expensive
#       count100<-count100*0
#     }
#
#
#     ########################################################
#     ## decorrelation step for theta, especially to decorrelate the thetas that dont change anything
#
#     if(i%%decor==0){ # every so often, do a decorrelation step (use single-site independence sampler)
#       for(k in 1:p){
#         theta.cand<-theta[i,,,drop=F] # most up-to-date value
#         theta.cand[1,,k]<-runif(ntemps) # independence sampler candidate (vectorize over temperatures)
#         # pred.cand<-matrix(
#         #   predict(mod,
#         #           tran(matrix(theta.cand[1,,],ncol=p)),
#         #           mcmc.use=sample(ns,size=1),
#         #           trunc.error=trunc.error,
#         #           nugget=T,
#         #           n.cores=pred.ncores,
#         #           parType = pred.parType),
#         #   nrow=ntemps)
#         pred.cand<-matrix(mod(tran(matrix(theta.cand[1,,],ncol=p))),nrow=ntemps)
#         for(t in 1:ntemps){
#           alpha<- -.5/s2[i-1,t]*itl[t]*(sum((y-pred.cand[t,])^2)-sum((y-pred.curr[t,])^2)) # could do with colsums, but this is cheap
#           if(log(runif(1))<alpha){
#             theta[i,t,k]<-theta.cand[1,t,k]
#             count.decor[k,t]<-count.decor[k,t]+1
#             pred.curr[t,]<-pred.cand[t,]
#           }
#         }
#       }
#     }
#
#     ########################################################
#     ## update s2 with gibbs step
#
#     s2[i,]<-rigammaTemper(ntemps,ny/2+a, b+colSums((t(pred.curr)-y)^2)/2, itl) # update error variance
#
#     #if(i==10000){
#     #  browser()
#     #}
#
#     #pred.curr<-predict(mod,theta[i,,],mcmc.use=sample(ns,size=1),trunc.error=F,nugget=T)[1,,]
#
#     ########################################################
#     ## tempering swaps
#
#     if(i>1000 & ntemps>1){ # tempering swap
#       for(dummy in 1:ntemps){ # repeat tempering step a bunch of times
#         sw<-sort(sample(ntemps,size=2))
#
#         alpha<-(itl[sw[2]]-itl[sw[1]])*(
#           -ny/2*log(s2[i,sw[1]]) - .5/s2[i,sw[1]]*sum((y-pred.curr[sw[1],])^2) -(a+1)*log(s2[i,sw[1]])-b/s2[i,sw[1]]
#           +ny/2*log(s2[i,sw[2]]) + .5/s2[i,sw[2]]*sum((y-pred.curr[sw[2],])^2) +(a+1)*log(s2[i,sw[2]])+b/s2[i,sw[2]]
#         )
#
#         if(log(runif(1))<alpha){
#           temp<-theta[i,sw[1],]
#           theta[i,sw[1],]<-theta[i,sw[2],]
#           theta[i,sw[2],]<-temp
#           temp<-s2[i,sw[1]]
#           s2[i,sw[1]]<-s2[i,sw[2]]
#           s2[i,sw[2]]<-temp
#           count[sw[1],sw[2]]<-count[sw[1],sw[2]]+1
#           temp<-pred.curr[sw[1],] # not sampling posterior predictive each time, for speed
#           pred.curr[sw[1],]<-pred.curr[sw[2],]
#           pred.curr[sw[2],]<-temp
#         }
#       }
#     }
#
#     ########################################################
#     ## take a new posterior predictive sample from the emulator (helps not get stuck in a mode from a particularly good sample)
#     # pred.curr<-matrix(
#     #   predict(mod,
#     #           tran(matrix(theta[i,,],ncol=p)),
#     #           mcmc.use=sample(ns,size=1),
#     #           trunc.error=trunc.error,
#     #           nugget=T,
#     #           n.cores=pred.ncores,
#     #           parType = pred.parType),
#     #   nrow=ntemps)
#
#     #if(!all(pred.curr == matrix(mod(tran(theta[i,,])),nrow=ntemps)))
#     #  browser()
#
#
#     if(verbose & i%%100==0){
#       pr<-c('MCMC iteration',i,myTimestamp(),'count:',diag(count))
#       cat(pr,'\n')
#     }
#   }
#
#   #th2<-pnorm(theta)
#   for(ii in 1:p){
#     theta[,,ii]<-unscale.range(theta[,,ii],bounds[,ii])
#   }
#   return(list(theta=theta,s2=s2,count=count,count.decor=count.decor,tau=tau))
# }
#
#
#
#
#
# # calibrateIndep<-function(mod,y,a,b,nmcmc,verbose=T){ # assumes inputs to mod are standardized to (0,1), equal variance for all y values (should change to sim covariance)
# #   p<-ncol(mod$dat$xx)
# #   ny<-length(y)
# #   ns<-mod$mod.list[[1]]$nmcmc-mod$mod.list[[1]]$nburn # number of emu mcmc samples
# #
# #   theta<-matrix(nrow=nmcmc,ncol=p)
# #   s2<-rep(NA,nmcmc)
# #
# #   #browser()
# #   theta[1,]<-.5
# #   pred.curr<-predict(mod,theta[1,,drop=F],mcmc.use=sample(ns,size=1),trunc.error=F)
# #   s2[1]<-1/rgamma(1,ny/2+a,b+sum((y-pred.curr)^2))
# #
# #   count<-rep(0,p)
# #   for(i in 2:nmcmc){
# #     s2[i]<-1/rgamma(1,ny/2+a,b+sum((y-pred.curr)^2))
# #
# #     theta[i,]<-theta[i-1,]
# #
# #     for(j in 1:p){
# #       theta.cand<-theta[i,]
# #       theta.cand[j]<-runif(1)
# #       pred.cand<-predict(mod,t(theta.cand),mcmc.use=sample(ns,size=1),trunc.error=F)
# #       alpha<- -.5/s2[i]*(sum((y-pred.cand)^2)-sum((y-pred.curr)^2))
# #       if(log(runif(1))<alpha){
# #         theta[i,]<-theta.cand
# #         count[j]<-count[j]+1
# #       }
# #     }
# #
# #     pred.curr<-predict(mod,theta[i,,drop=F],mcmc.use=sample(ns,size=1),trunc.error=F)
# #
# #     if(verbose & i%%100==0){
# #       pr<-c('MCMC iteration',i,myTimestamp(),'count:',count)
# #       cat(pr,'\n')
# #     }
# #   }
# #
# #   return(list(theta=theta,s2=s2,count=count))
# # }
#
# ldig.kernal<-function(x,a,b)
#   (-a-1)*log(x) - b/x
#
# calibrate.probit<-function(mod,y,s2.est,s2.df,bounds,nmcmc,tl=1,verbose=T,decor=100,pred.ncores=1,pred.parType='fork',trunc.error=F){ # assumes inputs to mod are standardized to (0,1), equal variance for all y values (should change to sim covariance)
#   # if(class(mod)=='bass'){
#   #   p<-mod$p
#   #   ns<-length(mod$s2) # number of emu mcmc samples
#   # } else if(class(mod)=='bassBasis'){
#   #   p<-ncol(mod$dat$xx)
#   #   ns<-length(mod$mod.list[[1]]$s2) # number of emu mcmc samples
#   # }
#
#   p<-ncol(bounds)
#   a<-s2.df/2
#   b<-a*sd.est^2
#
#   lpost<-NA
#   # could allow s2.ind vector, like in python code
#   # or could allow s2.mult for a functional setup
#
#   ny<-length(y)
#   ntemps<-length(tl)
#
#   theta<-array(dim=c(nmcmc,ntemps,p))
#   s2<-matrix(nrow=nmcmc,ncol=ntemps)
#   #temp.ind<-matrix(nrow=nmcmc,ncol=ntemps)
#
#   itl<-1/tl # inverse temperature ladder
#
#   tran<-function(th){
#     th<-pnorm(th)
#     for(i in 1:ncol(th)){
#       th[,i]<-unscale.range(th[,i],bounds[,i])
#     }
#     th
#   }
#
#   #browser()
#   theta[1,,]<-rnorm(prod(dim(theta[1,,,drop=F])))
#   # pred.curr<-matrix(
#   #   predict(mod,
#   #           tran(matrix(theta[1,,],ncol=p)),
#   #           mcmc.use=sample(ns,size=1),
#   #           trunc.error=trunc.error,
#   #           nugget = T,
#   #           n.cores=pred.ncores,
#   #           parType = pred.parType),
#   #   nrow=ntemps)
#   pred.curr<-matrix(mod(tran(matrix(theta[1,,],ncol=p))),nrow=ntemps)
#   s2[1,]<-1/rgammaTemper(ntemps, ny/2+a, b+colSums((t(pred.curr)-y)^2)/2, itl)
#
#
#
#   eps<-1e-13
#   tau<-rep(-0,ntemps) # scaling
#   cc<-2.4^2/p
#   S<-mu<-cov<-S.rev<-mu.rev<-cov.rev<-list()
#   for(t in 1:ntemps)
#     S[[t]]<-diag(p)*1e-6
#   count<-matrix(0,nrow=ntemps,ncol=ntemps)
#   count.decor<-matrix(0,nrow=p,ncol=ntemps)
#   count100<-rep(0,ntemps)
#
#
#   for(i in 2:nmcmc){
#
#     theta[i,,]<-theta[i-1,,] # current set at previous (update below)
#
#     ########################################################
#     ## adaptive block update for theta within each temperature - covariance of previous samples scaled (by exp(tau)) based on acceptance rate for last 100 samples, since tempering makes large gaps
#
#     if(i==300){ # start adapting
#       for(t in 1:ntemps){
#         mu[[t]]<-colMeans(theta[1:(i-1),t,])
#         cov[[t]]<-cov(theta[1:(i-1),t,])
#         S[[t]]<-cov(theta[1:(i-1),t,])*cc+diag(eps*cc,p)
#       }
#     }
#     if(i>300){ # adaptation updates
#       for(t in 1:ntemps){
#         mu[[t]]<-mu[[t]]+(theta[(i-1),t,]-mu[[t]])/(i-1)
#         cov[[t]]<-(i-2)/(i-1)*cov[[t]] + (i-2)/(i-1)^2*tcrossprod(theta[(i-1),t,]-mu[[t]])
#         S[[t]]<-(cov[[t]]*cc+diag(eps*cc,p))*exp(tau[t])
#       }
#     }
#
#     theta.cand<-matrix(nrow=ntemps,ncol=p)
#     for(t in 1:ntemps)
#       theta.cand[t,]<-rmnorm(theta[i-1,t,],S[[t]]) # generate candidate for each temperature
#
#     # pred.cand<-matrix(
#     #   predict(mod,
#     #           tran(matrix(theta.cand,ncol=p)),
#     #           mcmc.use=sample(ns,size=1),
#     #           trunc.error=trunc.error,
#     #           nugget=T,
#     #           n.cores=pred.ncores,
#     #           parType = pred.parType),
#     #   nrow=ntemps) # get BASS prediction at each candidate (major speedup by vectorizing across temperatures)
#     pred.cand<-matrix(mod(tran(matrix(theta.cand,ncol=p))),nrow=ntemps)
#
#     for(t in 1:ntemps){ # loop over temperatures, do a block MCMC update
#         #alpha<- (-.5/s2[i-1,t]*itl[t]*(sum((y-pred.cand[t,])^2)-sum((y-pred.curr[t,])^2)) + itl[t]*(sum(dnorm(theta.cand[t,],log=T)) - sum(dnorm(theta[i-1,t,],log=T))) )
#         alpha<-itl[t]*(
#           sum(dnorm(y,pred.cand[t,],sqrt(s2[i-1,t]),log=T))
#           -sum(dnorm(y,pred.curr[t,],sqrt(s2[i-1,t]),log=T))
#           +sum(dnorm(theta.cand[t,],log=T))
#           -sum(dnorm(theta[i-1,t,],log=T))
#         )
#       if(log(runif(1))<alpha){
#         theta[i,t,]<-theta.cand[t,]
#         count[t,t]<-count[t,t]+1
#         pred.curr[t,]<-pred.cand[t,]
#         count100[t]<-count100[t]+1
#       }
#     }
#
#     ########################################################
#     ## acceptance-rate based adaptation of covariance scale
#     if(i%%100==0){
#       delta<-min(.1,1/sqrt(i)*5)
#       for(t in 1:ntemps){
#         if(count100[t]<23){
#           tau[t]<-tau[t]-delta
#         } else if(count100[t]>23){
#           tau[t]<-tau[t]+delta
#         }
#       } # could be vectorized, but probably not expensive
#       count100<-count100*0
#     }
#
#
#     ########################################################
#     ## decorrelation step for theta, especially to decorrelate the thetas that dont change anything
#
#     if(i%%decor==0 & i>1000){ # every so often, do a decorrelation step (use single-site independence sampler)
#       for(k in 1:p){
#         theta.cand<-theta[i,,,drop=F] # most up-to-date value
#         theta.cand[1,,k]<-rnorm(ntemps) # independence sampler candidate (vectorize over temperatures)
#         #pred.cand<-matrix(
#           # predict(mod,
#           #         tran(matrix(theta.cand[1,,],ncol=p)),
#           #         mcmc.use=sample(ns,size=1),
#           #         trunc.error=trunc.error,
#           #         nugget=T,
#           #         n.cores=pred.ncores,
#           #         parType = pred.parType),
#           # nrow=ntemps)
#           pred.cand<-matrix(mod(tran(matrix(theta.cand[1,,],ncol=p))),nrow=ntemps)
#
#         for(t in 1:ntemps){
#           #alpha<- -.5/s2[i-1,t]*itl[t]*(sum((y-pred.cand[t,])^2)-sum((y-pred.curr[t,])^2)) + itl[t]*(sum(dnorm(theta.cand[1,t,],log=T)) - sum(dnorm(theta[i,t,],log=T))) # could do with colsums, but this is cheap
#           alpha<-itl[t]*(
#             sum(dnorm(y,pred.cand[t,],sqrt(s2[i-1,t]),log=T))
#             -sum(dnorm(y,pred.curr[t,],sqrt(s2[i-1,t]),log=T))
#             +dnorm(theta.cand[1,t,k],log=T)
#             -dnorm(theta[i,t,k],log=T)
#             )- dnorm(theta.cand[1,t,k],log=T)+ dnorm(theta[i,t,k],log=T)
#
#           if(log(runif(1))<alpha){
#             #browser()
#             theta[i,t,k]<-theta.cand[1,t,k]
#             count.decor[k,t]<-count.decor[k,t]+1
#             pred.curr[t,]<-pred.cand[t,]
#           }
#         }
#       }
#     }
#
#     ########################################################
#     ## update s2 with gibbs step
#
#     #browser()
#     s2[i,]<-1/rgammaTemper(ntemps,ny/2+a, b+colSums((t(pred.curr)-y)^2)/2, itl) # update error variance
#
#     ########################################################
#     ## tempering swaps
#
#     if(i>1000 & ntemps>1){ # tempering swap
#       #browser()
#       for(dummy in 1:ntemps){ # repeat tempering step a bunch of times
#         sw<-sort(sample(ntemps,size=2))
#
#         # alpha<-(itl[sw[2]]-itl[sw[1]])*(
#         #   -ny/2*log(s2[i,sw[1]]) - .5/s2[i,sw[1]]*sum((y-pred.curr[sw[1],])^2) -(a+1)*log(s2[i,sw[1]])-b/s2[i,sw[1]]
#         #   +ny/2*log(s2[i,sw[2]]) + .5/s2[i,sw[2]]*sum((y-pred.curr[sw[2],])^2) +(a+1)*log(s2[i,sw[2]])+b/s2[i,sw[2]]
#         #   + sum(dnorm(theta[i,sw[1],],log=T)) - sum(dnorm(theta[i,sw[2],],log=T))
#         # )
#
#         alpha<-(itl[sw[2]]-itl[sw[1]])*(
#           sum(dnorm(y,pred.curr[sw[1],],sqrt(s2[i,sw[1]]),log=T))
#           -sum(dnorm(y,pred.curr[sw[2],],sqrt(s2[i,sw[2]]),log=T))
#           +sum(dnorm(theta[i,sw[1],],log=T))
#           -sum(dnorm(theta[i,sw[2],],log=T))
#           +ldig.kernal(s2[i,sw[1]],a,b)
#           -ldig.kernal(s2[i,sw[2]],a,b)
#         )
#
#         if(log(runif(1))<alpha){
#           if(sw[1]==0 & any(abs(theta[i,sw[2],])>4.5)){
#             print('bad')
#           }
#           #temp<-theta[i,sw[1],]
#           #theta[i,sw[1],]<-theta[i,sw[2],]
#           #theta[i,sw[2],]<-temp
#           theta[i,sw,]<-theta[i,sw[2:1],]
#           #temp<-s2[i,sw[1]]
#           #s2[i,sw[1]]<-s2[i,sw[2]]
#           #s2[i,sw[2]]<-temp
#           s2[i,sw]<-s2[i,sw[2:1]]
#           count[sw[1],sw[2]]<-count[sw[1],sw[2]]+1
#           #temp<-pred.curr[sw[1],] # not sampling posterior predictive each time, for speed
#           #pred.curr[sw[1],]<-pred.curr[sw[2],]
#           #pred.curr[sw[2],]<-temp
#           pred.curr[sw,]<-pred.curr[sw[2:1],]
#         }
#       }
#     }
#
#     ########################################################
#     ## take a new posterior predictive sample from the emulator (helps not get stuck in a mode from a particularly good sample)
#     # pred.curr<-matrix(
#     #   predict(mod,
#     #           tran(matrix(theta[i,,],ncol=p)),
#     #           mcmc.use=sample(ns,size=1),
#     #           trunc.error=trunc.error,
#     #           nugget=T,
#     #           n.cores=pred.ncores,
#     #           parType = pred.parType),
#     #   nrow=ntemps)
#
#     lpost[i]<-(sum(dnorm(y,pred.curr[1,],sqrt(s2[i,1]),log=T))
#     +sum(dnorm(theta[i,1,],log=T))
#     +ldig.kernal(s2[i,1],a,b)
#     )
#
#
#     if(verbose & i%%1000==0){
#       pr<-c('MCMC iteration',i,myTimestamp(),'count:',diag(count))
#       cat(pr,'\n')
#     }
#   }
#
#   #th2<-pnorm(theta)
#   #for(ii in 1:p){
#   #  theta[,,ii]<-unscale.range(pnorm(theta[,,ii]),bounds[,ii])
#   #}
#   return(list(theta=theta,s2=s2,count=count,count.decor=count.decor,tau=tau,lpost=lpost))
# }
#
#
#
#
#
#
#
# #calibrate.full<-function(mod,y,s2.est,s2.df,bounds,nmcmc,tl=1,verbose=T,decor=100,pred.ncores=1,pred.parType='fork',trunc.error=F){ # assumes inputs to mod are standardized to (0,1), equal variance for all y values (should change to sim covariance)
#
# #}
#
#
# calibrate.naive.cut<-function(mod,y,s2.est,s2.df,bounds,nmcmc,tl=1,verbose=T,decor=100,pred.ncores=1,pred.parType='fork',trunc.error=F,stoch=F){ # assumes inputs to mod are standardized to (0,1), equal variance for all y values (should change to sim covariance)
#
#
#   p<-ncol(bounds)
#   a<-s2.df/2
#   b<-a*sd.est^2
#
#   # could allow s2.ind vector, like in python code
#   # or could allow s2.mult for a functional setup
#
#   ny<-length(y)
#   ntemps<-length(tl)
#
#   theta<-array(dim=c(nmcmc,ntemps,p))
#   s2<-matrix(nrow=nmcmc,ncol=ntemps)
#   #temp.ind<-matrix(nrow=nmcmc,ncol=ntemps)
#
#   itl<-1/tl # inverse temperature ladder
#
#   tran<-function(th){
#     #th2<-pnorm(th)
#     for(i in 1:ncol(th)){
#       th[,i]<-unscale.range(th[,i],bounds[,i])
#     }
#     th
#   }
#
#   #browser()
#   theta[1,,]<-runif(prod(dim(theta[1,,,drop=F])))
#   # pred.curr<-matrix(
#   #   predict(mod,
#   #           tran(matrix(theta[1,,],ncol=p)),
#   #           mcmc.use=sample(ns,size=1),
#   #           trunc.error=trunc.error,
#   #           nugget = T,
#   #           n.cores=pred.ncores,
#   #           parType = pred.parType),
#   #   nrow=ntemps)
#
#   pred.curr<-matrix(mod(tran(matrix(theta[1,,],ncol=p))),nrow=ntemps)
#   s2[1,]<-rigammaTemper(ntemps, ny/2+a, b+colSums((t(pred.curr)-y)^2)/2, itl)
#
#
#
#   eps<-1e-13
#   tau<-rep(-4,ntemps) # scaling
#   cc<-2.4^2/p
#   S<-mu<-cov<-S.rev<-mu.rev<-cov.rev<-list()
#   for(t in 1:ntemps)
#     S[[t]]<-diag(p)*1e-6
#   count<-matrix(0,nrow=ntemps,ncol=ntemps)
#   count.decor<-matrix(0,nrow=p,ncol=ntemps)
#   count100<-rep(0,ntemps)
#
#
#   for(i in 2:nmcmc){
#
#     theta[i,,]<-theta[i-1,,] # current set at previous (update below)
#
#     ########################################################
#     ## adaptive block update for theta within each temperature - covariance of previous samples scaled (by exp(tau)) based on acceptance rate for last 100 samples, since tempering makes large gaps
#
#     if(i==300){ # start adapting
#       for(t in 1:ntemps){
#         mu[[t]]<-colMeans(theta[1:(i-1),t,])
#         cov[[t]]<-cov(theta[1:(i-1),t,])
#         S[[t]]<-cov(theta[1:(i-1),t,])*cc+diag(eps*cc,p)
#       }
#     }
#     if(i>300){ # adaptation updates
#       for(t in 1:ntemps){
#         mu[[t]]<-mu[[t]]+(theta[(i-1),t,]-mu[[t]])/(i-1)
#         cov[[t]]<-(i-2)/(i-1)*cov[[t]] + (i-2)/(i-1)^2*tcrossprod(theta[(i-1),t,]-mu[[t]])
#         S[[t]]<-(cov[[t]]*cc+diag(eps*cc,p))*exp(tau[t])
#       }
#     }
#
#
#
#     theta.cand<-matrix(nrow=ntemps,ncol=p)
#     for(t in 1:ntemps)
#       theta.cand[t,]<-rmnorm(theta[i-1,t,],S[[t]]) # generate candidate for each temperature
#
#
#     # pred.cand<-matrix(
#     #   predict(mod,
#     #           tran(matrix(theta.cand,ncol=p)),
#     #           mcmc.use=sample(ns,size=1),
#     #           trunc.error=trunc.error,
#     #           nugget=T,
#     #           n.cores=pred.ncores,
#     #           parType = pred.parType),
#     #   nrow=ntemps) # get BASS prediction at each candidate (major speedup by vectorizing across temperatures)
#
#     pred.cand<-matrix(mod(tran(matrix(theta.cand,ncol=p))),nrow=ntemps)
#
#     for(t in 1:ntemps){ # loop over temperatures, do a block MCMC update
#       if(any(theta.cand[t,]<0 | theta.cand[t,]>1))
#         alpha<- -9999
#       else{
#         alpha<- (-.5/s2[i-1,t]*itl[t]*(sum((y-pred.cand[t,])^2)-sum((y-pred.curr[t,])^2)) # posterior
#                  #-mnormt::dmnorm(theta.cand[t,],theta[i-1,t,],S[[t]],log = T) + mnormt::dmnorm(theta[i-1,t,],theta.cand[t,],S.rev[[t]],log = T)  # proposal when using radius-based adapting
#         )
#       }
#       if(log(runif(1))<alpha){
#         theta[i,t,]<-theta.cand[t,]
#         count[t,t]<-count[t,t]+1
#         pred.curr[t,]<-pred.cand[t,]
#         count100[t]<-count100[t]+1
#       }
#     }
#
#     ########################################################
#     ## acceptance-rate based adaptation of covariance scale
#     if(i%%100==0){
#       delta<-min(.1,1/sqrt(i)*5)
#       for(t in 1:ntemps){
#         if(count100[t]<23){
#           tau[t]<-tau[t]-delta
#         } else if(count100[t]>23){
#           tau[t]<-tau[t]+delta
#         }
#       } # could be vectorized, but probably not expensive
#       count100<-count100*0
#     }
#
#
#     ########################################################
#     ## decorrelation step for theta, especially to decorrelate the thetas that dont change anything
#
#     if(i%%decor==0){ # every so often, do a decorrelation step (use single-site independence sampler)
#       for(k in 1:p){
#         theta.cand<-theta[i,,,drop=F] # most up-to-date value
#         theta.cand[1,,k]<-runif(ntemps) # independence sampler candidate (vectorize over temperatures)
#         # pred.cand<-matrix(
#         #   predict(mod,
#         #           tran(matrix(theta.cand[1,,],ncol=p)),
#         #           mcmc.use=sample(ns,size=1),
#         #           trunc.error=trunc.error,
#         #           nugget=T,
#         #           n.cores=pred.ncores,
#         #           parType = pred.parType),
#         #   nrow=ntemps)
#         pred.cand<-matrix(mod(tran(matrix(theta.cand[1,,],ncol=p))),nrow=ntemps)
#         for(t in 1:ntemps){
#           alpha<- -.5/s2[i-1,t]*itl[t]*(sum((y-pred.cand[t,])^2)-sum((y-pred.curr[t,])^2)) # could do with colsums, but this is cheap
#           if(log(runif(1))<alpha){
#             theta[i,t,k]<-theta.cand[1,t,k]
#             count.decor[k,t]<-count.decor[k,t]+1
#             pred.curr[t,]<-pred.cand[t,]
#           }
#         }
#       }
#     }
#
#     ########################################################
#     ## update s2 with gibbs step
#
#     s2[i,]<-rigammaTemper(ntemps,ny/2+a, b+colSums((t(pred.curr)-y)^2)/2, itl) # update error variance
#
#     #pred.curr<-predict(mod,theta[i,,],mcmc.use=sample(ns,size=1),trunc.error=F,nugget=T)[1,,]
#
#     ########################################################
#     ## tempering swaps
#
#     if(i>1000 & ntemps>1){ # tempering swap
#       for(dummy in 1:ntemps){ # repeat tempering step a bunch of times
#         sw<-sort(sample(ntemps,size=2))
#
#         alpha<-(itl[sw[2]]-itl[sw[1]])*(
#           -ny/2*log(s2[i,sw[1]]) - .5/s2[i,sw[1]]*sum((y-pred.curr[sw[1],])^2) -(a+1)*log(s2[i,sw[1]])-b/s2[i,sw[1]]
#           +ny/2*log(s2[i,sw[2]]) + .5/s2[i,sw[2]]*sum((y-pred.curr[sw[2],])^2) +(a+1)*log(s2[i,sw[2]])+b/s2[i,sw[2]]
#         )
#
#         if(log(runif(1))<alpha){
#           temp<-theta[i,sw[1],]
#           theta[i,sw[1],]<-theta[i,sw[2],]
#           theta[i,sw[2],]<-temp
#           temp<-s2[i,sw[1]]
#           s2[i,sw[1]]<-s2[i,sw[2]]
#           s2[i,sw[2]]<-temp
#           count[sw[1],sw[2]]<-count[sw[1],sw[2]]+1
#           temp<-pred.curr[sw[1],] # not sampling posterior predictive each time, for speed
#           pred.curr[sw[1],]<-pred.curr[sw[2],]
#           pred.curr[sw[2],]<-temp
#         }
#       }
#     }
#
#     ########################################################
#     ## take a new posterior predictive sample from the emulator (helps not get stuck in a mode from a particularly good sample)
#     # pred.curr<-matrix(
#     #   predict(mod,
#     #           tran(matrix(theta[i,,],ncol=p)),
#     #           mcmc.use=sample(ns,size=1),
#     #           trunc.error=trunc.error,
#     #           nugget=T,
#     #           n.cores=pred.ncores,
#     #           parType = pred.parType),
#     #   nrow=ntemps)
#
#     #if(!all(pred.curr == matrix(mod(tran(theta[i,,])),nrow=ntemps)))
#     #  browser()
#
#
#     if(verbose & i%%100==0){
#       pr<-c('MCMC iteration',i,myTimestamp(),'count:',diag(count))
#       cat(pr,'\n')
#     }
#   }
#
#   #th2<-pnorm(theta)
#   for(ii in 1:p){
#     theta[,,ii]<-unscale.range(theta[,,ii],bounds[,ii])
#   }
#   return(list(theta=theta,s2=s2,count=count,count.decor=count.decor,tau=tau))
#
#
# }
