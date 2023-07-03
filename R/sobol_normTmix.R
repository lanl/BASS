#######################################################
# Author: Devin Francom, Los Alamos National Laboratory
# Protected under GPL-3 license
# Los Alamos Computer Code release C19031
# github.com/lanl/BASS
# Full copyright in the README.md in the repository
#######################################################

############################################################
## get Sobol decomposition
############################################################

#' @title BASS Sensitivity Analysis
#'
#' @description Decomposes the variance of the BASS model into variance due to main effects, two way interactions, and so on, similar to the ANOVA decomposition for linear models.  Uses the Sobol' decomposition, which can be done analytically for MARS models.
#' @param bassMod a fitted model output from the \code{bass} function.
#' @param prior a list of priors; uniform, truncated mixture of Normals or Ts for continuous; vector of category weights for categorical.  Default is uniform over range of data.
#' @param prior.func prior for functional variable.  In almost all cases, keep this as the uniform default.
#' @param mcmc.use an integer vector indexing which MCMC iterations to use for sensitivity analysis.
#' @param func.var an integer indicating which functional variable to make sensitivity indices a function of.  Disregard if \code{bassMod} is non-functional or if scalar sensitivity indices are desired.
#' @param xx.func.var grid for functional variable specified by \code{func.var}.  Disregard if \code{func.var} is not specified.  If \code{func.var} is specified and \code{xx.func.var} not specified, the grid used to fit \code{bass} will be used.
#' @param verbose logical; should progress be displayed?
#' @param getEffects logical; should Sobols ANOVA decomposition be computed?
#' @details Performs analytical Sobol' decomposition for each MCMC iteration in mcmc.use (each corresponds to a MARS model), yeilding a posterior distribution of sensitivity indices.  Can obtain Sobol' indices as a function of one functional variable.
#' @return If non-functional (\code{func.var = NULL}), a list with two elements:
#'  \item{S}{a data frame of sensitivity indices with number of rows matching the length of \code{mcmc.use}.  The columns are named with a particular main effect or interaction.  The values are the proportion of variance in the model that is due to each main effect or interaction.}
#'  \item{T}{a data frame of total sensitivity indices with number of rows matching the length of \code{mcmc.use}.  The columns are named with a particular variable.}
#'  Otherwise, a list with four elements:
#'  \item{S}{an array with first dimension corresponding to MCMC samples (same length as \code{mcmc.use}), second dimension corresponding to different main effects and interactions (labeled in \code{names.ind}), and third dimension corresponding to the grid used for the functional variable.  The elements of the array are sensitivity indices.}
#'  \item{S.var}{same as \code{S}, but scaled in terms of total variance rather than percent of variance.}
#'  \item{names.ind}{a vector of names of the main effects and interactions used.}
#'  \item{xx}{the grid used for the functional variable.}
#'
#' @keywords Sobol decomposition
#' @seealso \link{bass} for model fitting and \link{predict.bass} for prediction.
#' @import hypergeo
#' @export
#' @examples
#' # See examples in bass documentation.
#'
sobol<-function(bassMod,prior=NULL,prior.func=NULL,mcmc.use=NULL,func.var=NULL,xx.func.var=NULL,verbose=TRUE,getEffects=FALSE){
  if(!inherits(bassMod,'bass'))
    stop('First input needs to be a bass object')
  if(bassMod$p==1 & !bassMod$func)
    stop('Sobol only used for multiple input models')
  if(bassMod$p==1 & sum(func.var)>0)
    stop('Cannot decompose the variance in terms of only one variable')
  mcmc.use.poss<-1:((bassMod$nmcmc-bassMod$nburn)/bassMod$thin)
  if(any(!(mcmc.use%in%mcmc.use.poss))){
    mcmc.use<-mcmc.use.poss
    warning('disregarding mcmc.use because of bad values')
  }
  if(any(is.null(mcmc.use))){
    mcmc.use<-mcmc.use.poss
  }

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
        prior[[i]]$trunc<-bassMod$range.des[,i]#c(0,1)
      }

      prior[[i]]$trunc<-scale_range(prior[[i]]$trunc,bassMod$range.des[,i])
      if(prior[[i]]$trunc[1]<0 | prior[[i]]$trunc[2]>1)
        warning('truncation range larger than training range...it is unwise to ask an emulator to extrapolate.')
      #browser()

      if(prior[[i]]$dist %in% c('normal','student')){
        prior[[i]]$mean<-scale_range(prior[[i]]$mean,bassMod$range.des[,i])
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

  if(bassMod$func){
    for(i in 1:length(prior.func)){
      if(is.null(prior.func[[i]]$trunc)){
        prior.func[[i]]$trunc<-bassMod$range.func[,i]#c(0,1)
      }

      prior.func[[i]]$trunc<-scale_range(prior.func[[i]]$trunc,bassMod$range.func[,i])
      if(prior.func[[i]]$trunc[1]<0 | prior.func[[i]]$trunc[2]>1)
        stop('truncation range of functional variable larger than training range...it is unwise to ask an emulator to extrapolate.')
      #browser()

      if(prior.func[[i]]$dist %in% c('normal','student')){
        prior.func[[i]]$mean<-scale_range(prior.func[[i]]$mean,bassMod$range.func[,i])
        prior.func[[i]]$sd<-prior.func[[i]]$sd/(bassMod$range.func[2,i]-bassMod$range.func[1,i])
        if(prior.func[[i]]$dist == 'normal'){
          prior.func[[i]]$z<-pnorm((prior.func[[i]]$trunc[2]-prior.func[[i]]$mean)/prior.func[[i]]$sd) - pnorm((prior.func[[i]]$trunc[1]-prior.func[[i]]$mean)/prior.func[[i]]$sd)
        } else{
          prior.func[[i]]$z<-pt((prior.func[[i]]$trunc[2]-prior.func[[i]]$mean)/prior.func[[i]]$sd,prior.func[[i]]$df) - pt((prior.func[[i]]$trunc[1]-prior.func[[i]]$mean)/prior.func[[i]]$sd,prior.func[[i]]$df)
        }
        cc<-sum(prior.func[[i]]$weights*prior.func[[i]]$z)
        prior.func[[i]]$weights<-prior.func[[i]]$weights/cc#prior[[i]]$z # change weights with truncation # divide by cc instead to keep the same prior shape
        # does the truncation change the distribution shape in the non-truncated regions??
        #browser()
      }
    }
  }

  # if(bassMod$func){
  #   for(i in 1:length(prior.func)){
  #     if(is.null(prior.func[[i]]$trunc)){
  #       prior.func[[i]]$trunc<-c(0,1)
  #     } else{
  #       prior.func[[i]]$trunc<-scale_range(prior.func[[i]]$trunc,bassMod$range.func[,i])
  #     }
  #
  #     if(prior.func[[i]]$dist %in% c('normal','student')){
  #       prior.func[[i]]$mean<-scale_range(prior.func[[i]]$mean,bassMod$range.func[,i])
  #       prior.func[[i]]$sd<-prior.func[[i]]$sd/(bassMod$range.func[2,i]-bassMod$range.func[1,i])
  #       if(prior.func[[i]]$dist == 'normal'){
  #         prior.func[[i]]$z<-pnorm((prior.func[[i]]$trunc[2]-prior.func[[i]]$mean)/prior.func[[i]]$sd) - pnorm((prior.func[[i]]$trunc[1]-prior.func[[i]]$mean)/prior.func[[i]]$sd)
  #       } else{
  #         prior.func[[i]]$z<-pt((prior.func[[i]]$trunc[2]-prior.func[[i]]$mean)/prior.func[[i]]$sd,prior.func[[i]]$df) - pt((prior.func[[i]]$trunc[1]-prior.func[[i]]$mean)/prior.func[[i]]$sd,prior.func[[i]]$df)
  #       }
  #       prior.func[[i]]$weights<-prior.func[[i]]$weights/prior.func[[i]]$z # change weights with truncation
  #     }
  #   }
  # }
  #
  #
  # prior.func[[func.var]]<-NULL
  # prior<-c(prior,prior.func)
  # check to see if this worked how we expected, pass prior and prior.cat around


  prior<-c(prior,prior.func)

  #  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  if(getEffects){
    if(bassMod$cat | bassMod$func){
      getEffects<-F
      warning('getEffects not yet implemented for functional response or categorical inputs.')
    }
  }

  if(is.null(func.var)){
    func<-F
  } else{
    func<-T
    if(!bassMod$func){
      func<-F
      warning('disregarding func.var because bassMod parameter is not functional')
    }
  }

  if(!is.null(xx.func.var) & !func)
    warning('disregarding xx.func.var because bassMod parameter is not functional')

  if(func){
    if(!(func.var%in%(1:ncol(bassMod$xx.func))))
      stop('func.var in wrong range of values')
    if(is.null(xx.func.var)){
      xx.func.var<-bassMod$xx.func[,func.var,drop=F]
    } else{
      rr<-range(xx.func.var)
      if(rr[1]<bassMod$range.func[1,func.var] | rr[2]>bassMod$range.func[2,func.var])
        warning(paste('range of func.var in bass function (',bassMod$range.func[1,func.var],',',bassMod$range.func[2,func.var],') is smaller than range of xx.func.var (',rr[1],',',rr[2],'), indicating some extrapolation',sep=''))
      xx.func.var<-scale_range(xx.func.var,bassMod$range.func[,func.var])
    }
    return(sobol_des_func(bassMod=bassMod,mcmc.use=mcmc.use,verbose=verbose,func.var=func.var,xx.func.var=xx.func.var,prior=prior,prior.cat=prior.cat))
  } else{
    return(sobol_des(bassMod=bassMod,mcmc.use=mcmc.use,verbose=verbose,prior=prior,prior.cat=prior.cat,getEffects=getEffects)) # applies to both des & func as long as functional sobol indices are not desired
  }
}






## get sobol indices - no functional

sobol_des<-function(bassMod,mcmc.use,verbose,prior,prior.cat,getEffects){
  models<-bassMod$model.lookup[mcmc.use] # only do the heavy lifting once for each model
  uniq.models<-unique(models)
  nmodels<-length(uniq.models)
  maxInt.tot<-sum(bassMod$maxInt.des)+sum(bassMod$maxInt.cat)+sum(bassMod$maxInt.func)
  maxBasis<-max(bassMod$nbasis)
  q<-bassMod$degree
  pdes<-sum(bassMod$pdes)
  pcat<-sum(bassMod$pcat)
  pfunc<-sum(bassMod$pfunc)
  p<-pdes+pcat+pfunc



  #prior<-c(prior,prior.func)

  ################################################
  # get combs including functional variables
  #browser()
  allCombs<-getCombs(bassMod,uniq.models,nmodels,maxBasis,maxInt.tot)
  combs<-allCombs$combs # which main effects and interactions included
  names.ind<-allCombs$names.ind # which main effects and interactions included
  num.ind<-allCombs$num.ind # number of terms in each interaction level
  cs.num.ind<-allCombs$cs.num.ind # cumsum of num.ind
  ################################################
  sob<-array(0,dim=c(length(mcmc.use),sum(num.ind)))
  var.tot.store<-f0.store<-rep(0,length(mcmc.use))

  ngrid<-100
  #xx=seq(0,1,length.out=ngrid) # make a different size xx for each variable...to debug
  xxt<-lapply(1:p,function(i) t(seq(0,1,length.out=ngrid+i)))
  effects<-list()
  if(getEffects){
    for(iint in 1:min(length(combs),2)){
      effects[[iint]]<-array(dim=c(length(mcmc.use),ncol(combs[[iint]]),rep(ngrid,iint)))
    } # do above somewhere
  }

  if(verbose)
    cat('Sobol Start',myTimestamp(),'Models:',length(unique(models)),'\n')

  i<-1
  mod.count<-0
  for(mod in uniq.models){ #do this in parallel?
    mod.count<-mod.count+1 # number of models
    mcmc.use.mod<-mcmc.use[models==mod] # which part of mcmc.use does this correspond to?
    mod.ind<-i:(i+length(mcmc.use.mod)-1) # index for which mcmc samples are this model?
    M<-bassMod$nbasis[mcmc.use.mod][1] # number of basis functions in this model
    if(M>0){
      # for each model, m, tl stores everything necessary to get sobol indices
      tl<-get_tl(bassMod,mcmc.use.mod,M,mod,p,q,cs.num.ind,combs) # tl initially stores coefs, knots, signs, variables, which combinations included, etc.
      tl$prior<-prior
      tl$prior.cat<-prior.cat
      if(tl$cat)
        tl<-add_tlCat(tl,bassMod,mcmc.use.mod,mod) # if there are categorical variables, add the associated info
      tl<-add_tl(tl,p) # calculates and adds the C1 and C2 integrals to tl
      lens<-apply(tl$Kind,1,function(x) length(na.omit(x))) # number of interactions in each basis function
#browser()
      var.tot<-Vu(1:p,tl) # total variance, one for each mcmc iteration (since coefs change each mcmc iteration rather than each model)


      # vt<-myCC<-myC2<-matrix(nrow=M,ncol=M)
      # for(i in 1:M){
      #   for(j in 1:M){
      #     prod1<-prod2<-1
      #     for(k in 1:p){
      #       tt1<-C1(k,i,tl)
      #       tt2<-C1(k,j,tl)
      #       tt3<-C2(k,i,j,tl)
      #       # if(!any(na.omit(tl$Kind[i,]==k)))
      #       #   tt1=1
      #       # if(!any(na.omit(tl$Kind[j,]==k)))
      #       #   tt2=1
      #       # if(!any(na.omit(intersect(tl$Kind[i,], tl$Kind[j,])==k)))
      #       #   tt3=1
      #       prod1<-prod1*tt1*tt2
      #       prod2<-prod2*tt3#~~~~~~~~~~~~~~~~~~~~~~~~~ finish this debugging!!!!
      #     }
      #     myCC[i,j]<-prod1
      #     myC2[i,j]<-prod2
      #     vt[i,j]=prod2-prod1#prod1*(prod2/prod1-1)
      #   }
      # }
      #
      # u=1:p
      # CCu<-apply(tl$C1.all.prod3[,,u,drop=F],1:2,prod) # product over u's TODO: utilize symmetry
      # C2.temp<-apply(tl$C2.all2[,,u,drop=F],1:2,prod)
      # mat<-tl$CC*(C2.temp/CCu-1)
      # mat[tl$CC==0]<-0
      #
      # C2.temp[CCu==0]
      # CCu[C2.temp==0]
      #
      # CCu[1,4]
      # myCC[1,4]
      # tl$CC[1,4]
      # C2.temp[1,4]
      # myC2[1,4]
      #
      # VuMat(1:p,tl)[4,1]
      # vt[4,1]
      # range(VuMat(1:p,tl)-vt)
      # (tl$a%*%vt%*%t(tl$a))
      # (tl$a%*%mat%*%t(tl$a))
      # var.tot
      # browser()

      vars.used<-unique(unlist(na.omit(c(tl$Kind)))) # which variables are used in this model?
      vars.used<-sort(vars.used)

      tl$integrals<-matrix(0,nrow=length(mcmc.use.mod),ncol=max(cs.num.ind)) # where we store all the integrals (not normalized by subtraction) - matches dim of sob

      ## get main effect variances

      tl$integrals[,which(combs[[1]]%in%vars.used)]<-apply(t(vars.used),2,Vu,tl=tl) # one for each effect for each mcmc iteration, only done for the used vars, but inserted into the correct column of tl$integrals
      sob[mod.ind,1:cs.num.ind[1]]<-tl$integrals[,1:cs.num.ind[1]] # main effects are done

      ## get interactions

      # consider changing lens to n.int
      if(max(lens)>1){ # if there are any interactions
        for(l in 2:max(lens)){ # must go in order for it to work (tl$integrals is made sequentially)
          int.l.ind<-(cs.num.ind[l-1]+1):cs.num.ind[l] # index for which columns of sob will have this level of interaction, and they are names.ind[[l]]
          basis.int.l.use<-matrix(nrow=M,ncol=num.ind[l]) # basis.int.l.use[i,j]=1 if interaction combs[[l]][,j] occurs in ith basis function, shows us which interactions (of order l) are used in which basis functions
          for(m in 1:M){
            #basis.int.l.use[m,]<-unlist(lapply(combs.list[[l]],function(el){prod(el%in%tl$Kind[m,])})) #all the variables in question must be in the basis
            basis.int.l.use[m,]<-apply(combs[[l]],2,function(el){prod(el%in%tl$Kind[m,])})
          }
          use<-colSums(basis.int.l.use)>0 # indicator for which interactions (of order l) are used at all in this model
          tl$integrals[,int.l.ind[use]]<-apply(combs[[l]][,use,drop=F],2,Vu,tl=tl) # perform the necessary integration, but only for interactions that are actually used
          sob.l<-matrix(0,nrow=length(mcmc.use.mod),ncol=num.ind[l])
          for(l2 in (1:num.ind[l])[use]){ # only go through the interactions that are actually in Kind (but still allow for 2-way when actual is 3-way, etc.)
            sob.l[,l2]<-VuInt(combs[[l]][,l2],tl) # do the normalizing (subtract out lower order terms to make this orthogonal to them)
          }
          sob[mod.ind,int.l.ind]<-sob.l # variance explained
        }
      }
      sob[mod.ind,]<-sob[mod.ind,]/var.tot # scale so that percent variance explained
      i<-i+length(mcmc.use.mod) # for index
      var.tot.store[mod.ind]<-var.tot # store the total variance for this model
      #browser()
      C1.temp<-(1/(tl$q+1)*((tl$s+1)/2-tl$s*tl$t))*tl$s^2
      C1.temp[C1.temp==0]<-1

      f0.store[mcmc.use.mod]<-bassMod$beta[mcmc.use.mod,1]+bassMod$beta[mcmc.use.mod,2:(tl$M+1)]%*%matrix(apply(C1.temp,1,prod))
      #browser()


      main.effects<-T
      # if(length(combs[[1]])==0)
      #   main.effects<-F

      two.ints<-F
      if(length(combs)>1)
        two.ints<-T

      if(getEffects){
##################################################################
      # main effects - trying to get fi - f0, can leave off the intercepts (a0 in Chen 2004)

        f0<-bassMod$beta[mcmc.use.mod,2:(tl$M+1)]%*%matrix(apply(C1.temp,1,prod)) # without a0
        if(main.effects){
        for(ef in 1:ncol(combs[[1]])){ # go through each effect
          effects[[1]][mod.ind,ef,]<- -f0
          for(m in 1:tl$M){ # go through each basis function
            pp.use<-combs[[1]][,ef] # which variable is this effect
            if(tl$s[m,pp.use]!=0){ # if the basis function uses the variable
              effects[[1]][mod.ind,ef,]<-effects[[1]][mod.ind,ef,]+tcrossprod(bassMod$beta[mcmc.use.mod,m+1],makeBasis(tl$s[m,pp.use],1,tl$t[m,pp.use],xxt,1)*prod(C1.temp[m,-pp.use]))
            } else{ # if the basis function does not use the variable (integrate over all)
              effects[[1]][mod.ind,ef,]<-effects[[1]][mod.ind,ef,]+bassMod$beta[mcmc.use.mod,m+1]*prod(C1.temp[m,])
            }
          }
          #matplot(t(effects[[1]][1,,]),type='l')
          #matplot(cbind(10*sin(pi*xx)^2/(pi*xx) - a1, 10*sin(pi*xx)^2/(pi*xx) - a1, 20*(xx-.5)^2 - 5/3, 10*xx - 5, 5*xx-5/2),type='l')
          #matplot(t(effects[[1]][,1,]),type='l')
        }

        }

        if(two.ints){
        ## 2-way interactions
        for(ef in 1:ncol(combs[[2]])){
          for(kk in 1:length(mcmc.use.mod))
            effects[[2]][mod.ind[kk],ef,,]<- -f0[kk]

          pp.use<-combs[[2]][,ef]
          pp.use1<-which(combs[[1]]%in%combs[[2]][1,ef])
          pp.use2<-which(combs[[1]]%in%combs[[2]][2,ef])

            effects[[2]][mod.ind,ef,,]<- sweep(effects[[2]][mod.ind,ef,,,drop=F], c(1,2,3), effects[[1]][mod.ind,pp.use1,,drop=F])

            effects[[2]][mod.ind,ef,,]<- sweep(effects[[2]][mod.ind,ef,,,drop=F], c(1,2,4), effects[[1]][mod.ind,pp.use2,,drop=F])

            #image.plot(effects[[2]][mcmc.use.mod[1],ef,,])


          for(m in 1:tl$M){

              b1<-makeBasis(tl$s[m,pp.use[1]],1,tl$t[m,pp.use[1]],xxt,1)
              b2<-makeBasis(tl$s[m,pp.use[2]],1,tl$t[m,pp.use[2]],xxt,1)
              if(all(b1==0))
                b1=b1+1
              if(all(b2==0))
                b2=b2+1

              effects[[2]][mod.ind,ef,,]<-effects[[2]][mod.ind,ef,,] + drop(bassMod$beta[mcmc.use.mod,m+1]%o%(tcrossprod(b1,b2)*prod(C1.temp[m,-pp.use])))

            }
          }
        #browser()
        #image.plot(effects[[2]][1,1,,],zlim=c(-10,10))
        #image.plot(matrix(10*sin(2*pi*xx2[,1]*xx2[,2]),ncol=100))

        #xx2<-expand.grid(t(xxt),t(xxt))
        #a1<--5/pi*sum( (-4*pi^2)^(1:50)/((2*(1:50))*factorial(2*(1:50))) )
        #image.plot(matrix(10*sin(2*pi*xx2[,1]*xx2[,2]) - 10*sin(pi*xx2[,1])^2/(pi*xx2[,1]) - 10*sin(pi*xx2[,2])^2/(pi*xx2[,2]) + a1,nrow=100),zlim=c(-10,10))

        }
        }
          ##################################################################


      if(verbose & mod.count%%10==0)
        cat('Sobol',myTimestamp(),'Model:',mod.count,'\n')


    }
  }

  sob<-as.data.frame(sob)
  names(sob)<-unlist(names.ind) # give labels


  if(verbose)
    cat('Total Sensitivity',myTimestamp(),'\n')

  tot<-getTot(combs,sob,names.ind) # get total indices

  ## reorder socolumns of sob & tot matrices to match original data order

  sob.reorder<-NA
  sob.reorder[1:length(names.ind[[1]])]<-mixOrd(allCombs$dispNames[[1]])
  if(length(names.ind)>1){
    for(l in 2:length(names.ind)){
      sob.reorder[(cs.num.ind[l-1]+1):(cs.num.ind[l])]<-cs.num.ind[l-1]+mixOrd(allCombs$dispNames[[l]])
    }
  }
  sob<-sob[,sob.reorder,drop=F]
  names(sob)<-unlist(allCombs$dispNames)[sob.reorder]
  tot<-tot[,sob.reorder[1:length(names.ind[[1]])],drop=F]
  names(tot)<-allCombs$dispNames[[1]][sob.reorder[1:length(names.ind[[1]])]]

  if(any(sob<0))
    browser()

  ret<-list(S=sob,T=tot,func=F,var.tot=var.tot.store,f0=f0.store,ints=tl$integrals,prior=prior,effects=effects,names.ind=names.ind)
  class(ret)<-'bassSob'

  return(ret)
}






## get sobol indices - functional

sobol_des_func<-function(bassMod,mcmc.use,verbose,func.var,xx.func.var,prior,prior.cat){
  models<-bassMod$model.lookup[mcmc.use] # only do the heavy lifting once for each model
  uniq.models<-unique(models)
  nmodels<-length(uniq.models)
  maxInt.tot<-sum(bassMod$maxInt.des)+sum(bassMod$maxInt.cat)+sum(bassMod$maxInt.func)
  maxBasis<-max(bassMod$nbasis)
  q<-bassMod$degree
  pdes<-sum(bassMod$pdes)
  pcat<-sum(bassMod$pcat)
  pfunc<-sum(bassMod$pfunc)
  p<-pdes+pcat+pfunc

  ################################################
  # get combs including functional variables
  allCombs<-getCombs(bassMod,uniq.models,nmodels,maxBasis,maxInt.tot,func.var)
  combs<-allCombs$combs
  names.ind<-allCombs$names.ind
  num.ind<-allCombs$num.ind
  cs.num.ind<-allCombs$cs.num.ind
  ################################################
  sob<-sob2<-array(0,dim=c(length(mcmc.use),sum(num.ind),length(xx.func.var)))

  #prior.func[[func.var]]<-NULL
  #prior<-c(prior,prior.func)

  if(verbose)
    cat('Sobol Start',myTimestamp(),'Models:',length(unique(models)),'\n')

  i<-1
  mod.count<-0
  for(mod in uniq.models){ #do this in parallel?
    mod.count<-mod.count+1
    mcmc.use.mod<-mcmc.use[models==mod] # which part of mcmc.use does this correspond to?
    mod.ind<-i:(i+length(mcmc.use.mod)-1)
    M<-bassMod$nbasis[mcmc.use.mod][1] # number of basis functions in this model
    if(M>0){
      tl<-get_tl(bassMod,mcmc.use.mod,M,mod,p,q,cs.num.ind,combs,func.var,xx.func.var)
      tl$prior<-prior
      tl$prior.cat<-prior.cat
      if(tl$cat)
        tl<-add_tlCat(tl,bassMod,mcmc.use.mod,mod)
      tl<-add_tl(tl,p)
      lens<-apply(tl$Kind,1,function(x) length(na.omit(x)))

      tl$tfunc.basis<-makeBasisMatrixVar(mod,M,vars=bassMod$vars.func,signs=bassMod$signs.func,knots.ind=bassMod$knotInd.func,q=bassMod$degree,xxt=t(tl$xx),n.int=bassMod$n.int.func,xx.train=bassMod$xx.func,var=func.var)[-1,,drop=F]

      var.tot<-Vu_des_func(1:p,tl) # total variance
      vars.used<-unique(unlist(na.omit(c(tl$Kind)))) # which variables are used?
      vars.used<-sort(vars.used)

      tl$integrals<-array(0,dim=c(length(mcmc.use.mod),max(cs.num.ind),length(xx.func.var))) # where we store all the integrals (not normalized by subtraction)
      jj=0
      for(pp in vars.used){
        #jj=jj+1
        tl$integrals[,which(combs[[1]]==pp),]<-Vu_des_func(pp,tl)
      }



      sob[mod.ind,1:cs.num.ind[1],]<-tl$integrals[,1:cs.num.ind[1],]

      if(max(lens)>1){ # if there are any interactions
        for(l in 2:max(lens)){ # must go in order for it to work (tl$integrals is made sequentially)
          int.l.ind<-(cs.num.ind[l-1]+1):cs.num.ind[l]
          #basis.int.l.use<-matrix(nrow=M,ncol=length(combs.list[[l]]))
          basis.int.l.use<-matrix(nrow=M,ncol=num.ind[l])

          for(m in 1:M){
            #basis.int.l.use[m,]<-unlist(lapply(combs.list[[l]],function(el){prod(el%in%tl$Kind[m,])})) # all the variables in question must be in the basis
            basis.int.l.use[m,]<-apply(combs[[l]],2,function(el){prod(el%in%tl$Kind[m,])})
          }
          use<-colSums(basis.int.l.use)>0
          for(pp in which(use)){
            tl$integrals[,int.l.ind[pp],]<-Vu_des_func(combs[[l]][,pp,drop=F],tl) # perform the necessary integration
          }

          sob.l<-array(0,dim=c(length(mcmc.use.mod),num.ind[l],length(xx.func.var)))
          for(l2 in (1:num.ind[l])[use]){ # only go through the interactions that are actually in Kind (but still allow for 2-way when actual is 3-way, etc.)
            sob.l[,l2,]<-VuInt_des_func(combs[[l]][,l2],tl) # do the normalizing
          }
          sob[mod.ind,int.l.ind,]<-sob.l
        }
      }

      kk<-0
      for(ii in mod.ind){
        kk=kk+1
        sob2[ii,,]<-t(t(sob[ii,,])/var.tot[kk,]) # sobol indices
      }
      i<-i+length(mcmc.use.mod) # for index

      if(verbose & mod.count%%10==0)
        cat('Sobol',myTimestamp(),'Model:',mod.count,'\n')
    }
  }




  # reorder for display
  sob.reorder<-NA
  sob.reorder[1:length(names.ind[[1]])]<-mixOrd(allCombs$dispNames[[1]])
  if(length(names.ind)>1){
    for(l in 2:length(names.ind)){
      sob.reorder[(cs.num.ind[l-1]+1):(cs.num.ind[l])]<-cs.num.ind[l-1]+mixOrd(allCombs$dispNames[[l]])
    }
  }
  sob<-sob[,sob.reorder,,drop=F]
  sob2<-sob2[,sob.reorder,,drop=F]


  #if(any(sob2<0))
  #  browser()

  ret<-list(S=sob2,S.var=sob,names.ind=unlist(allCombs$dispNames)[sob.reorder],xx=unscale.range(tl$xx,bassMod$range.func),func=T,prior=prior)
  class(ret)<-'bassSob'
  return(ret)
}

# get total sensitivity indices
getTot<-function(combs,sob,names.ind){
  vars.use<-unique(unlist(combs))
  puse<-length(vars.use)
  ncombs<-sapply(combs,ncol)
  tot<-sob[,1:ncombs[1]] # start with main effects, then add to them
  ncombs[1]<-0
  if(length(combs)>1){ # if there are interactions
    int.use<-2:length(combs) # this should work because the lower order combs have to be there
    for(pp in 1:puse){ # go through variables
      for(l in int.use){ # go through interactions
        tot[,pp]<-tot[,pp]+rowSums(
          sob[,
              puse + # after main effects
              ncombs[l-1] + # after interactions of order l-1
              which(apply(t(combs[[l]]),1,function(r){vars.use[pp]%in%r})) # which combs use this variable
            ,drop=F]
        )
      }
    }
  }
  return(tot)
}



########################################################################
## processing functions
########################################################################

## make for only one variable
makeBasisMatrixVar<-function(i,nbasis,vars,signs,knots.ind,q,xxt,n.int,xx.train,var){
  n<-ncol(xxt)
  tbasis.mat<-matrix(nrow=nbasis+1,ncol=n)
  tbasis.mat[1,]<-1
  if(nbasis>0){
    for(m in 1:nbasis){
      if(all(na.omit(vars[i,m,])!=var)){
        tbasis.mat[m+1,]<-1
      } else{
        use<-which(vars[i,m,]==var)
        knots<-xx.train[cbind(knots.ind[i,m,use],vars[i,m,use])]
        tbasis.mat[m+1,]<-makeBasis(signs[i,m,use],1,knots,xxt,q)
      }
    }
  }
  return(tbasis.mat)
}

## get all the variable combinations used in the models, storing in proper structures
getCombs<-function(bassMod,uniq.models,nmodels,maxBasis,maxInt.tot,func.var=NULL){
  #browser()
  des.labs<-which(bassMod$cx%in%c('numeric','integer'))
  cat.labs<-which(bassMod$cx=='factor')
  func.labs<-letters[0:sum(bassMod$pfunc)]
  labs<-c(des.labs,func.labs,cat.labs) # this is the order things end up in
  vf<-bassMod$vars.func[uniq.models,,]
  sub<-0
  if(!is.null(func.var)){
    vf[vf==func.var]<-NA
    sub=-1
  }
  n.un<-array(c(as.integer(bassMod$vars.des[uniq.models,,]),as.integer(vf+sum(bassMod$pdes)),as.integer(bassMod$vars.cat[uniq.models,,]+sum(bassMod$pdes)+sum(bassMod$pfunc))),dim=c(nmodels,maxBasis,maxInt.tot)) # all the variables/interactions used (with non-overlapping variable indices)
  #n.un<-apply(n.un,1:2,sort) # this has list properties, n.un[1,] is a list (for the first model) with elements for basis function interactions - something wrong here when there are no NAs
  x<-apply(n.un,3,function(x) x)
  n.un<-lapply(1:nrow(x),function(i) sort(x[i,]))
  n.un<-unique(c(n.un)) # unique combinations
  n.un[sapply(n.un,length)==0]<-NULL # get rid of list elements with nothing in them
  int.lower<-list() # will hold lower order combinations
  for(ii in 1:length(n.un)){
    pp<-length(n.un[[ii]])
    if(pp==1){
      int.lower<-c(int.lower,n.un[[ii]])
    } else{
      int.lower<-c(int.lower,do.call(c,sapply(1:pp,function(x) combn(n.un[[ii]],x,simplify=F)))) # get all the combinations, add them to the list
    }
  }
  #browser()
  int.lower<-lapply(int.lower,as.integer)
  n.un<-unique(c(n.un,unique(int.lower))) # this now has all the lower order combinations
  ord.intSize<-order(sapply(n.un,length))
  n.un<-n.un[ord.intSize] # order by interaction size
  int.begin.ind<-NA # a[i] is the index in n.un where ith order interactions begin
  for(ii in 1:maxInt.tot)
    int.begin.ind[ii]<-which(sapply(n.un,length)==ii)[1]
  int.begin.ind[maxInt.tot+1]<-length(n.un)+1 # need the top level to help with indexing later

  combs<-names.ind<-dispNames<-dispCombs<-mat<-list()
  ints.used<-(1:maxInt.tot)[!is.na(int.begin.ind[-(maxInt.tot+1)])] # which interactions we actually use in the models
  int.begin.ind<-na.omit(int.begin.ind) # will help index things
  k<-0
  for(ii in ints.used){ # go through used interactions
    k<-k+1
    if(!is.na(int.begin.ind[ii])){ # probably don't need this if anymore
      mat<-do.call(rbind,n.un[int.begin.ind[k]:(int.begin.ind[k+1]-1)]) # get all the ineractions of order ii
      mat<-mat[do.call(order, as.data.frame(mat)),,drop=F] # sort them (order works the way we want for ordering vectors)
      combs[[ii]]<-t(mat) # the combinations of order ii
      names.ind[[ii]]<-apply(combs[[ii]],2,paste,collapse='x') # labels for output
      dispCombs[[ii]]<-matrix(labs[combs[[ii]]],nrow = nrow(combs[[ii]])) # takes into account categorical and functional labels
      dd<-dim(dispCombs[[ii]])
      dispCombs[[ii]]<-apply(dispCombs[[ii]],2,mixSort) # sorts numbers and letters the way we want
      dim(dispCombs[[ii]])<-dd
      dispNames[[ii]]<-apply(dispCombs[[ii]],2,paste,collapse='x') # same as names.ind, but better indexing for display (with functional, categorical variables)
      #combs.list[[ii]]<-split(combs[[ii]],rep(1:ncol(combs[[ii]]),each=nrow(combs[[ii]]))) # list version of combs[[ii]]
    }
  }

  num.ind<-sapply(combs,ncol) # num.ind[i] is number of interactions of order i
  cs.num.ind<-cumsum(num.ind) # used for indexing
  return(list(combs=combs,names.ind=names.ind,num.ind=num.ind,cs.num.ind=cs.num.ind,dispCombs=dispCombs,dispNames=dispNames))
}

## process model information into a temporary list
get_tl<-function(bassMod,mcmc.use.mod,M,mod,p,q,cs.num.ind,combs,func.var=NULL,xx.func.var=NULL){
  a<-bassMod$beta[mcmc.use.mod,2:(M+1),drop=F] # basis coefficients excluding intercept
  vf<-bassMod$vars.func[mod,1:M,]
  pfunc<-sum(bassMod$pfunc)
  if(!is.null(func.var)){
    vf[vf==func.var]<-NA # if there is a functional variable specified, we make it NA so that we don't integrate over it
    pfunc<-pfunc-1 # disregarding the func.var so we don't integrate over it
  }

  if(bassMod$des){
    Kind<-cbind(bassMod$vars.des[mod,1:M,],vf+bassMod$pdes) # Kind[i,] is the variables used in the ith basis function, including relevant functional ones
  } else if(bassMod$func){
    Kind<-matrix(vf)
  } else{
    Kind<-NA
  }

  if(M==1){
    Kind<-t(Kind)
  }
  p.df<-sum(bassMod$pdes)+sum(bassMod$pfunc)
  t<-s<-matrix(0,nrow=M,ncol=sum(bassMod$pdes)+sum(bassMod$pfunc))
  if(p.df>0){
    for(k in 1:M){ # these matrices mimic the output of earth
      if(bassMod$des){
        n.int.des<-bassMod$n.int.des[mod,k]
        knotInd.des<-bassMod$knotInd.des[mod,k,1:n.int.des]
        vars.des<-bassMod$vars.des[mod,k,1:n.int.des]
        t[k,vars.des]<-bassMod$xx.des[cbind(knotInd.des,vars.des)]
        s[k,vars.des]<-bassMod$signs.des[mod,k,1:n.int.des]
      }
      if(pfunc>0){
        if(bassMod$n.int.func[mod,k]>0){
          n.int.func<-bassMod$n.int.func[mod,k]
          knotInd.func<-bassMod$knotInd.func[mod,k,1:n.int.func]
          vars.func<-bassMod$vars.func[mod,k,1:n.int.func]
          t[k,vars.func+sum(bassMod$pdes)]<-bassMod$xx.func[cbind(knotInd.func,vars.func)]
          s[k,vars.func+sum(bassMod$pdes)]<-bassMod$signs.func[mod,k,1:n.int.func]
        }
      }
    }
  }
  if(!is.null(func.var)){
    s[,bassMod$pdes+func.var]<-t[,bassMod$pdes+func.var]<-0
  }
  tl<-list(s=s,t=t,q=q,a=a,M=M,Kind=Kind,cs.num.ind=cs.num.ind,combs=combs,xx=xx.func.var,pfunc=sum(bassMod$pfunc),cat=bassMod$cat,pdes=sum(bassMod$pdes)) #temporary list
  return(tl)
}




## process model information into a temporary list - categorical part
add_tlCat<-function(tl,bassMod,mcmc.use.mod,mod){
  tl$pcat<-bassMod$pcat
  tl$sub.cnt<-matrix(0,nrow=tl$M,ncol=tl$pcat)
  tl$sub<-list()
  for(mm in 1:tl$M){
    vars<-na.omit(bassMod$vars.cat[mod,mm,])
    tl$sub.cnt[mm,vars]<-bassMod$sub.size[mod,mm,bassMod$vars.cat[mod,mm,]%in%vars]/bassMod$nlevels[vars]
    tl$sub[[mm]]<-list()
    for(k in 1:tl$pcat){
      tl$sub[[mm]][[k]]<-NA
      if(k %in% vars)
        tl$sub[[mm]][[k]]<-bassMod$sub.list[[mod]][[mm]][[which(vars==k)]]
    }
  }
  p.df<-sum(bassMod$pdes)+sum(tl$pfunc)
  if(p.df>0)
    tl$Kind<-cbind(tl$Kind,bassMod$vars.cat[mod,1:tl$M,]+sum(bassMod$pdes)+sum(tl$pfunc))
  if(p.df==0)
    tl$Kind<-bassMod$vars.cat[mod,1:tl$M,]
  if(is.null(nrow(tl$Kind)))
    tl$Kind<-matrix(tl$Kind)
  tl$nlevels<-bassMod$nlevels
  return(tl)
}

## process model information into a temporary list - evaluate integrals (from Chen 2004, but vectorized as much as possible)
add_tl<-function(tl,p){
  #browser()
  p.df<-sum(tl$pdes)+sum(tl$pfunc)
  p.use<-p
  if(p.df==0){
    C1.all.cat<-tl$sub.cnt
    C1.all<-C1.all.cat
  } else{
    #C1.all<-C1(tl)#(1/(tl$q+1)*((tl$s+1)/2-tl$s*tl$t))*tl$s^2 # so I don't need C function anymore
    C1.all<-C1All(tl) # C1.all[i,j] is C1 (from Chen, 2004) for variable j in basis function i
    if(tl$cat){
      C1.all.cat<-tl$sub.cnt
      C1.all.cat[C1.all.cat==0]<-1
      C1.all<-cbind(C1.all,C1.all.cat)
    }
  }

  #C1.all.temp<-replace(C1.all,which(C1.all==0,arr.ind=T),1) # for products, make 0's into 1's, see Eq 50 of Chen 2004 to understand why (if we didn't the products wouldn't work)
  #browser()
  C1.all.temp<-C1.all
  C1.all.prod<-apply(C1.all.temp,1,prod) # product over basis functions (from 1 to M)
  tl$CC<-tcrossprod(C1.all.prod) # Eq 35, this is the product from 1:p (M in their notation) for all combinations of basis functions
  C2.all<-C1.all.prod2<-both.ind<-array(0,dim=c(tl$M,tl$M,p.use))

  for(ii in 1:tl$M){
    for(jj in ii:tl$M){
      bb<-intersect(na.omit(tl$Kind[ii,]),na.omit(tl$Kind[jj,])) # variables that basis functions ii and jj have in common



      ## test
      C1.all.prod2[ii,jj,]<-C1.all.prod2[jj,ii,]<-C1.all.temp[ii,]*C1.all.temp[jj,]
      ##



      #if(length(bb)>0){
        #C1.all.prod2[ii,jj,]<-C1.all.prod2[jj,ii,]<-C1.all[ii,]*C1.all[jj,] # pairwise products of C1.all, without final product like tl$CC.  Nothing is ever multiplied by 0 because of if statement above
        #if(ii==11 & jj==11)
        #  browser()
        #both.ind[ii,jj,bb]<-both.ind[jj,ii,bb]<-1
        #bb.cat<-bb[bb>p.df]
        #bb.des<-bb[bb<=p.df]

      #browser()

        temp<-rep(0,p.use)
        #if(length(bb.des)>0){
        if(p.df>0)
          temp[1:p.df]<-apply(t(1:p.df),2,C2,m=ii,n=jj,tl=tl)
        #}
        #if(length(bb.cat)>0){
        if(tl$cat){
          temp[(p.df+1):p.use]<-apply(t(1:tl$pcat),2,C2Cat,m=ii,n=jj,tl=tl)
          #if(temp[(p.df+1):p.use]==0)
            #browser()
        }
        C2.all[ii,jj,]<-C2.all[jj,ii,]<-temp
      #}
    }
  }
  tl$C1.all.prod3<-C1.all.prod2
  tl$C1.all.prod3[C1.all.prod2==0]<-1 # again, for products...otherwise you have to deal with divide by 0 by putting mat[tl$CC==0]<-0 inside Vu function
  tl$C2.all2<-C2.all
  #tl$C2.all2[!as.logical(both.ind)]<-1
  #tl$C2.all2[tl$C2.all2==0]<-1

  #browser()
  #if(min(tl$C2.all2)<1e-13)
  #  browser()
  return(tl)
}

## sorting function with mixed numerical and character
mixSort<-function(x){
  ind<-is.na(suppressWarnings(as.numeric(x)))
  num<-which(!ind)
  char<-which(ind)
  return(c(sort(as.numeric(x[num])),x[char]))
}

## ordering function with mixed numerical and character
mixOrd<-function(x){
  ind<-is.na(suppressWarnings(as.numeric(x)))
  ord<-1:length(x)
  num<-which(!ind)
  char<-which(ind)
  return(c(ord[num][order(as.numeric(x[num]))],ord[char]))
}


########################################################################
## functions for Sobol decomposition - these all use scaling from const function
########################################################################


# C1<-function(tl){ # can have a different one for each variable
#   (1/(tl$q+1)*((tl$s+1)/2-tl$s*tl$t))*tl$s^2
# }
# C1All<-function(tl)
#   C1(tl)

C1All<-function(tl){
  M<-tl$M
  puse<-ncol(tl$s)
  out<-matrix(nrow=M,ncol=puse)
  for(m in 1:M){
    for(p in 1:puse){
      out[m,p]<-C1(p,m,tl)
    }
  }
  out
}


################################################################################
# integral from a to b of (x-t)^q * prior(x) when q positive integer
intabq1 <- function (prior,a,b,t,q) {
  UseMethod("intabq1", prior)
}

intabq1.uniform<-function(prior,a,b,t,q){
  1/(q+1)*((b-t)^(q+1)-(a-t)^(q+1)) * 1/(prior$trunc[2]-prior$trunc[1])
  #int<-integrate(function(x) (x-t)*dunif(x,prior$trunc[1],prior$trunc[2]),lower=a,upper=b)$value
}
intabq1.normal<-function(prior,a,b,t,q){
  if(q!=1)
    stop('degree other than 1 not supported for normal priors')
  out<-0
  for(k in 1:length(prior$weights)){
    zk<-pnorm(b,prior$mean[k],prior$sd[k]) - pnorm(a,prior$mean[k],prior$sd[k])
    #tnorm.mean.zk<-prior$mean[k]*zk - prior$sd[k]*(dnorm(b,prior$mean[k],prior$sd[k]) - dnorm(a,prior$mean[k],prior$sd[k]))

    ast<-(a-prior$mean[k])/prior$sd[k]
    bst<-(b-prior$mean[k])/prior$sd[k]
    dnb<-dnorm(bst)
    dna<-dnorm(ast)
    tnorm.mean.zk<-prior$mean[k]*zk - prior$sd[k]*(dnb - dna)
    out<-out+prior$weights[k]*(tnorm.mean.zk - t*zk)
    # in parens should match integrate(function(x){(x-t)*dnorm(x,prior$mean[k],prior$sd[k])},lower=a,upper=b)
  }
  out
}
intabq1.student<-function(prior,a,b,t,q){
  if(q!=1)
    stop('degree other than 1 not supported for student priors')
  out<-0
  for(k in 1:length(prior$weights)){
    # zk<-pnorm(b,prior$mean[k],prior$sd[k]) - pnorm(a,prior$mean[k],prior$sd[k])
    # ast<-(a-prior$mean[k])/prior$sd[k]
    # bst<-(b-prior$mean[k])/prior$sd[k]
    # dnb<-dnorm(bst)
    # dna<-dnorm(ast)
    # tnorm.mean.zk<-prior$mean[k]*zk - prior$sd[k]*(dnb - dna)
    # out<-out+prior$weights[k]*(tnorm.mean.zk - t*zk)
    # in parens should match integrate(function(x){(x-t)*dnorm(x,prior$mean[k],prior$sd[k])},lower=a,upper=b)
    #int<-prior$mean[k]+prior$sd[k]*integrate(function(x) (x-t)*dt(x,prior$df[k]),lower=(a-prior$mean[k])/prior$sd[k],upper=(b-prior$mean[k])/prior$sd[k])$value



    zk<-pt((b-prior$mean[k])/prior$sd[k],prior$df[k]) - pt((a-prior$mean[k])/prior$sd[k],prior$df[k])
    # ast<-(a-prior$mean[k])/prior$sd[k]
    # bst<-(b-prior$mean[k])/prior$sd[k]
    # dnb<-dt(bst)/prior$sd[k]
    # dna<-dt(ast)/prior$sd[k]
    # tnorm.mean.zk<-prior$mean[k]*zk - prior$sd[k]*(dnb - dna)
    # out<-out+prior$weights[k]*(tnorm.mean.zk - t*zk)

    #int<-integrate(function(x) (x-t)*dt((x-prior$mean[k])/prior$sd[k],prior$df[k])/prior$sd[k],lower=a,upper=b)$value
    #browser()
    int<-intx1Student(b,prior$mean[k],prior$sd[k],prior$df[k],t) - intx1Student(a,prior$mean[k],prior$sd[k],prior$df[k],t)
    #integrate(function(x) (x-t)*dt.scaled(x,prior$df[k],prior$mean[k],prior$sd[k]),lower=a,upper=b)
    out<-out+prior$weights[k]*int ## TODO: handle truncation, mixture...~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~v
  }
  #browser()
  out
}

.S3method("intabq1", "uniform")
.S3method("intabq1", "normal")
.S3method("intabq1", "student")

intx1Student<-function(x,m,s,v,t){
  #browser()
  # a<-1/2
  # b<-(1 + v)/2
  # c<-3/2
  # xx<- -(m - x)^2/(s^2 *v)
  # f21<-gsl::hyperg_2F1(a,c-b,c,1-1/(1-xx))/(1-xx)^a # instead of hyperg_2F1(a,b,c,x)
  # if(is.nan(f21))
  #   f21<-gsl::hyperg_2F1(a,b,c,xx)
  temp<-(s^2*v)/(m^2 + s^2*v - 2*m*x + x^2)
  -((v/(v + (m - x)^2/s^2))^(v/2) *
      sqrt(temp) *
      sqrt(1/temp) *
      (s^2*v* (sqrt(1/temp) -
                 (1/temp)^(v/2)) +
         (t-m)*(-1 + v)*(-m + x) *
         (1/temp)^(v/2) *
         robust2f1(1/2,(1 + v)/2,3/2,-(m - x)^2/(s^2 *v)) )) /
    (s *(-1 + v)* sqrt(v) *beta(v/2, 1/2))
}


C1<-function(k,m,tl,tq=F){ #k is variable, m is basis function # deals with sign & truncation
  t<-tl$t[m,k]
  s<-tl$s[m,k]
  #browser()
  if(s==0)
    return(1)
  q<-tl$q
  cc<-const(signs=s,knots=t,degree=q)
  #int<-integrate(function(x) pos(s*(x-t))*dunif(x,tl$prior[[k]]$trunc[1],tl$prior[[k]]$trunc[2]),lower=0,upper=1)$value
  #browser()
  #return(int/cc)
  if(tq){
    q<-2*q
    cc<-cc^2
  }
  #browser()
  if(s==1){
    a<-max(tl$prior[[k]]$trunc[1],t)
    b<-tl$prior[[k]]$trunc[2]
    if(b<t)
      return(0)
    out<-intabq1(tl$prior[[k]],a,b,t,q)/cc
    #return(intabq1(tl$prior[[k]],a,b,t,q)/cc)
  } else{
    a<-tl$prior[[k]]$trunc[1]
    b<-min(tl$prior[[k]]$trunc[2],t)
    if(t<a)
      return(0)
    out<-intabq1(tl$prior[[k]],a,b,t,q)*(-1)^q/cc
    #return(intabq1(tl$prior[[k]],a,b,t,q)*(-1)^q/cc)
  }
  if(out<0)
    browser()
  return(out)
}


## refer to francom 2016 paper
pCoef<-function(i,q){
  factorial(q)^2*(-1)^i/(factorial(q-i)*factorial(q+1+i))
}

################################################################################
## integral from a to b of [(x-t1)(x-t2)]^q * prior(x) when q positive integer
intabq2 <- function (prior,a,b,t1,t2,q) {
  UseMethod("intabq2", prior)
}
intabq2.uniform<-function(prior,a,b,t1,t2,q){
  (sum(pCoef(0:q,q)*(b-t1)^(q-0:q)*(b-t2)^(q+1+0:q)) - sum(pCoef(0:q,q)*(a-t1)^(q-0:q)*(a-t2)^(q+1+0:q))) * 1/(prior$trunc[2]-prior$trunc[1])
  #integrate(function(x) (x-t1)*(x-t2)*dunif(x,prior$trunc[1],prior$trunc[2]),lower=a,upper=b)$value
}
intabq2.normal<-function(prior,a,b,t1,t2,q){
  if(q>1)
    stop('spline degree >1 not supported yet')
  out<-0
  for(k in 1:length(prior$weights)){
    zk<-pnorm(b,prior$mean[k],prior$sd[k]) - pnorm(a,prior$mean[k],prior$sd[k])
    if(zk<.Machine$double.eps)
      next
    ast<-(a-prior$mean[k])/prior$sd[k]
    bst<-(b-prior$mean[k])/prior$sd[k]
    dnb<-dnorm(bst)
    dna<-dnorm(ast)
    tnorm.mean.zk<-prior$mean[k]*zk - prior$sd[k]*(dnb - dna)
    tnorm.var.zk<-zk*prior$sd[k]^2*(1 + (ast*dna-bst*dnb)/zk - ((dna-dnb)/zk)^2) + tnorm.mean.zk^2/zk # variance + expectation^2
    out<-out+prior$weights[k]*(tnorm.var.zk - (t1+t2)*tnorm.mean.zk + t1*t2*zk)
    # in parens should be integrate(function(x){(x-t1)*(x-t2)*dnorm(x,prior$mean[k],prior$sd[k])},lower=a,upper=b)
    if(out<0 & abs(out)<1e-12)
      out<-0
  }
  #browser()
  out
}
intabq2.student<-function(prior,a,b,t1,t2,q){
  out<-0
  for(k in 1:length(prior$weights)){
    #int<-prior$mean[k]+prior$sd[k]*integrate(function(x) (x-t1)*(x-t2)*dt(x,prior$df[k]),lower=(a-prior$mean[k])/prior$sd[k],upper=(b-prior$mean[k])/prior$sd[k])$value

   # browser()
    #int<-integrate(function(x) (x-t1)*(x-t2)*dt((x-prior$mean[k])/prior$sd[k],prior$df[k])/prior$sd[k],lower=a,upper=b)$value
    int<-intx2Student(b,prior$mean[k],prior$sd[k],prior$df[k],t1,t2) - intx2Student(a,prior$mean[k],prior$sd[k],prior$df[k],t1,t2)
    out<-out+prior$weights[k]*int
  }
  out
}

.S3method("intabq2", "uniform")
.S3method("intabq2", "normal")
.S3method("intabq2", "student")

# robust2f1<-function(a,b,c,x){
#   if(abs(x)<1)
#     return(gsl::hyperg_2F1(a,b,c,x))
#   return(gsl::hyperg_2F1(a,c-b,c,1-1/(1-x))/(1-x)^a)
# }

robust2f1<-function(a,b,c,x){
  if(abs(x)<1)
    return(hypergeo::f15.3.8(a,b,c,x))
  return(hypergeo::f15.3.8(a,c-b,c,1-1/(1-x))/(1-x)^a)
}
intx2Student<-function(x,m,s,v,t1,t2){
  #x=b;m=prior$mean[k];s=prior$sd[k];v=prior$df[k]
  temp<-(s^2*v)/(m^2 + s^2*v - 2*m*x + x^2)
  ((v/(v + (m - x)^2/s^2))^(v/2) *
      sqrt(temp) *
      sqrt(1/temp) *
      (-3*(-t1-t2+2*m)*s^2*v* (sqrt(1/temp) -
                 (1/temp)^(v/2)) +
         3*(-t1+m)*(-t2+m)*(-1 + v)*(-m + x) *
         (1/temp)^(v/2) *
         robust2f1(1/2,(1 + v)/2,3/2,-(m - x)^2/(s^2 *v)) +
         (-1+v)*(-m+x)^3*(1/temp)^(v/2) *
         robust2f1(3/2,(1 + v)/2,5/2,-(m - x)^2/(s^2 *v)) )) /
    (3*s *(-1 + v)* sqrt(v) *beta(v/2, 1/2))
}


## integral of two pieces of tensor that have same variable - deals with sign, truncation
C2<-function(k,m,n,tl){ # k is variable, n & m are basis indices
  #browser()
  q<-tl$q
  t1<-tl$t[n,k]
  s1<-tl$s[n,k]
  t2<-tl$t[m,k]
  s2<-tl$s[m,k]
  cc<-const(signs=c(s1,s2),knots=c(t1,t2),degree=q)
  # if((s1*s2)==0){
  #   return(0)
  # }

  ## test
  if(s1==0 & s2==0){
    #browser()
    return(1)
  }
  if(s1==0 & s2!=0){
    return(C1(k,m,tl))
  }
  if(s1!=0 & s2==0){
    return(C1(k,n,tl))
  }
  ##

  if(t2<t1){
    t1<-tl$t[m,k]
    s1<-tl$s[m,k]
    t2<-tl$t[n,k]
    s2<-tl$s[n,k]
  }

  #cc1<-const(signs=s1,knots=t1,degree=q)
  #cc2<-const(signs=s2,knots=t2,degree=q)
  #int<-integrate(function(x) pos(s1*(x-t1))/cc1*pos(s2*(x-t2))/cc2*dunif(x,tl$prior[[k]]$trunc[1],tl$prior[[k]]$trunc[2]),lower=0,upper=1,stop.on.error = F)$value
  #browser()
  #return(int)
  #if(m==n){ #t1=t2, s1=s2 # didn't need this part
    #return(1/(2*q+1)*((s1+1)/2-s1*t1)^(2*q+1)/cc)
   # browser()
  #  return(C1(k,m,tl,tq=T)) # this is the same as if you let it run below, maybe faster?

  #} #else{
    if(s1==1){
      if(s2==1){
        a<-max(t2,tl$prior[[k]]$trunc[1])
        b<-tl$prior[[k]]$trunc[2]
        if(a>=b)
          return(0)
        out<-intabq2(tl$prior[[k]],a,b,t1,t2,q)/cc
        #return(intabq2(tl$prior[[k]],a,b,t1,t2,q)/cc)
      } else{
        a<-max(t1,tl$prior[[k]]$trunc[1])
        b<-min(t2,tl$prior[[k]]$trunc[2])
        if(a>=b)
          return(0)
        out<-intabq2(tl$prior[[k]],a,b,t1,t2,q)*(-1)^q/cc
        #return(intabq2(tl$prior[[k]],a,b,t1,t2,q)*(-1)^q/cc)
      }
    } else{
      if(s2==1){
        return(0)
      } else{
        a<-tl$prior[[k]]$trunc[1]
        b<-min(t1,tl$prior[[k]]$trunc[2])
        if(a>=b)
          return(0)
        out<-intabq2(tl$prior[[k]],a,b,t1,t2,q)/cc
        #return(intabq2(tl$prior[[k]],a,b,t1,t2,q)/cc)
      }
    }
  #}

  if(abs(out)<.Machine$double.eps)
    out<-0
  if(out<0)
    browser()
  return(out)
}

## same as C2, but categorical
C2Cat<-function(k,m,n,tl){ # k is variable (categorical), m & n are basis functions

  if(tl$sub.cnt[n,k]==0 & tl$sub.cnt[m,k]==0){
    #browser()
    return(1)
  }

  if(tl$sub.cnt[n,k]==0){
    #browser()
    return(tl$sub.cnt[m,k])
  }

  if(tl$sub.cnt[m,k]==0)
    return(tl$sub.cnt[n,k])

  # return(length(na.omit(intersect(tl$sub[[m]][[k]],tl$sub[[n]][[k]])))/tl$nlevels[k])
  out<-length((intersect(tl$sub[[m]][[k]],tl$sub[[n]][[k]])))/tl$nlevels[k]
  #print(out)
  #if(out==0)
  #  browser()
  return(out)
}

## matrix used in sobol main effect variances - where most of the time is spent
VuMat<-function(u,tl){
  #browser()
  CCu<-apply(tl$C1.all.prod3[,,u,drop=F],1:2,prod) # product over u's TODO: utilize symmetry
  C2.temp<-apply(tl$C2.all2[,,u,drop=F],1:2,prod)
  mat<-tl$CC*(C2.temp/CCu-1) # Eq 35 of Chen 2004
  #mat[tl$CC==0]<-0
  return(mat)
}

## sobol main effect variances
Vu<-function(u,tl){
  mat<-VuMat(u,tl)
  out<-apply(tl$a,1,function(x) t(x)%*%mat%*%x) # Eq 35 of Chen 2004.  Should mat have non-negative eigenvalues?
  if(any(out<0))
    browser()
  return(out)
}

## functional sobol main effect variances
Vu_des_func<-function(u,tl){
  mat<-VuMat(u,tl)
  nx<-length(tl$xx)
  nmodels<-length(tl$a[,1])
  out<-matrix(nrow=nmodels,ncol=nx)
  for(i in 1:nmodels){
    for(j in 1:nx){
      tt<-tl$a[i,]*tl$tfunc.basis[,j]
      out[i,j]<-t(tt)%*%mat%*%tt # Eq 35 of Chen 2004, not integrating over one of the variables
    }
  }
  if(any(out<0))
    browser()
  return(out)
}

## sobol interaction variances
VuInt<-function(u,tl){
  add<-0
  len<-length(u)
  for(l in 1:len){
    ind<-((sum(tl$cs.num.ind[l-1])+1):tl$cs.num.ind[l])[apply(tl$combs[[l]],2,function(x) all(x%in%u))] # this gets index for which combs are subsets of u, sum(cs.num.ind[l-1]) makes it 0 when it should be
    add<-add+(-1)^(len-l)*rowSums(tl$integrals[,ind,drop=F])
  }
  add[abs(add)<1.5e-13]<-0
  if(any(add<0))
    browser()
  return(add)
}

## sobol interaction variances - functional
VuInt_des_func<-function(u,tl){
  add<-0
  len<-length(u)
  #browser()
  for(l in 1:len){
    ind<-((sum(tl$cs.num.ind[l-1])+1):tl$cs.num.ind[l])[apply(tl$combs[[l]],2,function(x) all(x%in%u))] # this gets index for which combs are subsets of u, sum() makes it 0 when it should be
    add<-add+(-1)^(len-l)*apply(tl$integrals[,ind,,drop=F],c(1,3),sum)
  }
  add[abs(add)<1e-13]<-0
  if(any(add<0))
    browser()
  return(add)
}

