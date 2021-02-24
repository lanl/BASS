##################################################################################################################################################################
##################################################################################################################################################################
## modularized calibration
rmnorm<-function(mu, S){
  mu+c(rnorm(length(mu))%*%chol(S))
}
calibrate<-function(mod,y,a,b,nmcmc,tl=1,verbose=T,decor=100,pred.ncores=1,pred.parType='fork'){ # assumes inputs to mod are standardized to (0,1), equal variance for all y values (should change to sim covariance)
  if(class(mod)=='bass'){
    p<-mod$p
    ns<-length(mod$s2) # number of emu mcmc samples
  } else if(class(mod)=='bassBasis'){
    p<-ncol(mod$dat$xx)
    ns<-length(mod$mod.list[[1]]$s2) # number of emu mcmc samples
  }

  ny<-length(y)
  ntemps<-length(tl)

  theta<-array(dim=c(nmcmc,ntemps,p))
  s2<-matrix(nrow=nmcmc,ncol=ntemps)
  #temp.ind<-matrix(nrow=nmcmc,ncol=ntemps)

  itl<-1/tl # inverse temperature ladder

  theta[1,,]<-runif(prod(dim(theta[1,,,drop=F])))
  pred.curr<-matrix(
    predict(mod,
            matrix(theta[1,,],ncol=p),
            mcmc.use=sample(ns,size=1),
            trunc.error=F,
            nugget = T,
            n.cores=pred.ncores,
            parType = pred.parType),
    nrow=ntemps)
  s2[1,]<-1/rgammaTemper(ntemps, ny/2+a, b+colSums((t(pred.curr)-y)^2)/2, itl)



  eps<-1e-13
  tau<-rep(-4,ntemps) # scaling
  cc<-2.4^2/p
  S<-mu<-cov<-S.rev<-mu.rev<-cov.rev<-list()
  for(t in 1:ntemps)
    S[[t]]<-diag(p)*1e-6
  count<-matrix(0,nrow=ntemps,ncol=ntemps)
  count.decor<-matrix(0,nrow=p,ncol=ntemps)
  count100<-rep(0,ntemps)


  for(i in 2:nmcmc){

    theta[i,,]<-theta[i-1,,] # current set at previous (update below)

    ########################################################
    ## adaptive block update for theta within each temperature - covariance of previous samples scaled (by exp(tau)) based on acceptance rate for last 100 samples, since tempering makes large gaps

    if(i==300){ # start adapting
      for(t in 1:ntemps){
        mu[[t]]<-colMeans(theta[1:(i-1),t,])
        cov[[t]]<-cov(theta[1:(i-1),t,])
        S[[t]]<-cov(theta[1:(i-1),t,])*cc+diag(eps*cc,p)
      }
    }
    if(i>300){ # adaptation updates
      for(t in 1:ntemps){
        mu[[t]]<-mu[[t]]+(theta[(i-1),t,]-mu[[t]])/(i-1)
        cov[[t]]<-(i-2)/(i-1)*cov[[t]] + (i-2)/(i-1)^2*tcrossprod(theta[(i-1),t,]-mu[[t]])
        S[[t]]<-(cov[[t]]*cc+diag(eps*cc,p))*exp(tau[t])
      }
    }

    # if(i>300 & i<1000){ # start adapting
    #   for(t in 1:ntemps){
    #     mu[[t]]<-colMeans(theta[1:(i-1),t,])
    #     cov[[t]]<-cov(theta[1:(i-1),t,])
    #     S[[t]]<-cov[[t]]*cc+diag(eps*cc,p)
    #   }
    # }
    # if(i>1000){ # radius-based adaptation
    #   for(t in 1:ntemps){
    #     #browser()
    #     dist<-sqrt(colSums((t(theta[1:(i-2),t,])-theta[i-1,t,])^2))
    #     use<-which(dist<=quantile(dist,.15))
    #     mu[[t]]<-colMeans(theta[use,t,])
    #     cov[[t]]<-cov(theta[use,t,])
    #     S[[t]]<-cov[[t]]*cc+diag(eps*cc,p)
    #   }
    # }

    theta.cand<-matrix(nrow=ntemps,ncol=p)
    for(t in 1:ntemps)
      theta.cand[t,]<-rmnorm(theta[i-1,t,],S[[t]]) # generate candidate for each temperature

    # if(i>1000){ # for reversibility when radius-based adapting
    #   for(t in 1:ntemps){
    #     #browser()
    #     dist<-sqrt(colSums((t(theta[1:(i-2),t,])-theta.cand[t,])^2))
    #     use<-which(dist<=quantile(dist,.15))
    #     mu.rev[[t]]<-colMeans(theta[use,t,])
    #     cov.rev[[t]]<-cov(theta[use,t,])
    #     S.rev[[t]]<-cov.rev[[t]]*cc+diag(eps*cc,p)
    #   }
    # }

    pred.cand<-matrix(
      predict(mod,
              matrix(theta.cand,ncol=p),
              mcmc.use=sample(ns,size=1),
              trunc.error=F,
              nugget=T,
              n.cores=pred.ncores,
              parType = pred.parType),
      nrow=ntemps) # get BASS prediction at each candidate (major speedup by vectorizing across temperatures)

    for(t in 1:ntemps){ # loop over temperatures, do a block MCMC update
      if(any(theta.cand[t,]<0 | theta.cand[t,]>1))
        alpha<- -9999
      else{
        alpha<- (-.5/s2[i-1,t]*itl[t]*(sum((y-pred.cand[t,])^2)-sum((y-pred.curr[t,])^2)) # posterior
                #-mnormt::dmnorm(theta.cand[t,],theta[i-1,t,],S[[t]],log = T) + mnormt::dmnorm(theta[i-1,t,],theta.cand[t,],S.rev[[t]],log = T)  # proposal when using radius-based adapting
        )
      }
      if(log(runif(1))<alpha){
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
        pred.cand<-matrix(
          predict(mod,
                  matrix(theta.cand[1,,],ncol=p),
                  mcmc.use=sample(ns,size=1),
                  trunc.error=F,
                  nugget=T,
                  n.cores=pred.ncores,
                  parType = pred.parType),
          nrow=ntemps)
        for(t in 1:ntemps){
          alpha<- -.5/s2[i-1,t]*itl[t]*(sum((y-pred.cand[t,])^2)-sum((y-pred.curr[t,])^2)) # could do with colsums, but this is cheap
          if(log(runif(1))<alpha){
            theta[i,t,k]<-theta.cand[1,t,k]
            count.decor[k,t]<-count.decor[k,t]+1
            pred.curr[t,]<-pred.cand[t,]
          }
        }
      }
    }

    ########################################################
    ## update s2 with gibbs step

    s2[i,]<-1/rgammaTemper(ntemps,ny/2+a, b+colSums((t(pred.curr)-y)^2)/2, itl) # update error variance

    #pred.curr<-predict(mod,theta[i,,],mcmc.use=sample(ns,size=1),trunc.error=F,nugget=T)[1,,]

    ########################################################
    ## tempering swaps

    if(i>1000 & ntemps>1){ # tempering swap
      for(dummy in 1:ntemps){ # repeat tempering step a bunch of times
        sw<-sort(sample(ntemps,size=2))

        alpha<-(itl[sw[2]]-itl[sw[1]])*(
          -ny/2*log(s2[i,sw[1]]) - .5/s2[i,sw[1]]*sum((y-pred.curr[sw[1],])^2) -(a+1)*log(s2[i,sw[1]])-b/s2[i,sw[1]]
          +ny/2*log(s2[i,sw[2]]) + .5/s2[i,sw[2]]*sum((y-pred.curr[sw[2],])^2) +(a+1)*log(s2[i,sw[2]])+b/s2[i,sw[2]]
        )

        if(log(runif(1))<alpha){
          temp<-theta[i,sw[1],]
          theta[i,sw[1],]<-theta[i,sw[2],]
          theta[i,sw[2],]<-temp
          temp<-s2[i,sw[1]]
          s2[i,sw[1]]<-s2[i,sw[2]]
          s2[i,sw[2]]<-temp
          count[sw[1],sw[2]]<-count[sw[1],sw[2]]+1
          temp<-pred.curr[sw[1],] # not sampling posterior predictive each time, for speed
          pred.curr[sw[1],]<-pred.curr[sw[2],]
          pred.curr[sw[2],]<-temp
        }
      }
    }

    ########################################################
    ## take a new posterior predictive sample from the emulator (helps not get stuck in a mode from a particularly good sample)
    pred.curr<-matrix(
      predict(mod,
              matrix(theta[i,,],ncol=p),
              mcmc.use=sample(ns,size=1),
              trunc.error=F,
              nugget=T,
              n.cores=pred.ncores,
              parType = pred.parType),
      nrow=ntemps)



    if(verbose & i%%100==0){
      pr<-c('MCMC iteration',i,myTimestamp(),'count:',diag(count))
      cat(pr,'\n')
    }
  }

  return(list(theta=theta,s2=s2,count=count,count.decor=count.decor,tau=tau))
}





# calibrateIndep<-function(mod,y,a,b,nmcmc,verbose=T){ # assumes inputs to mod are standardized to (0,1), equal variance for all y values (should change to sim covariance)
#   p<-ncol(mod$dat$xx)
#   ny<-length(y)
#   ns<-mod$mod.list[[1]]$nmcmc-mod$mod.list[[1]]$nburn # number of emu mcmc samples
#
#   theta<-matrix(nrow=nmcmc,ncol=p)
#   s2<-rep(NA,nmcmc)
#
#   #browser()
#   theta[1,]<-.5
#   pred.curr<-predict(mod,theta[1,,drop=F],mcmc.use=sample(ns,size=1),trunc.error=F)
#   s2[1]<-1/rgamma(1,ny/2+a,b+sum((y-pred.curr)^2))
#
#   count<-rep(0,p)
#   for(i in 2:nmcmc){
#     s2[i]<-1/rgamma(1,ny/2+a,b+sum((y-pred.curr)^2))
#
#     theta[i,]<-theta[i-1,]
#
#     for(j in 1:p){
#       theta.cand<-theta[i,]
#       theta.cand[j]<-runif(1)
#       pred.cand<-predict(mod,t(theta.cand),mcmc.use=sample(ns,size=1),trunc.error=F)
#       alpha<- -.5/s2[i]*(sum((y-pred.cand)^2)-sum((y-pred.curr)^2))
#       if(log(runif(1))<alpha){
#         theta[i,]<-theta.cand
#         count[j]<-count[j]+1
#       }
#     }
#
#     pred.curr<-predict(mod,theta[i,,drop=F],mcmc.use=sample(ns,size=1),trunc.error=F)
#
#     if(verbose & i%%100==0){
#       pr<-c('MCMC iteration',i,myTimestamp(),'count:',count)
#       cat(pr,'\n')
#     }
#   }
#
#   return(list(theta=theta,s2=s2,count=count))
# }

