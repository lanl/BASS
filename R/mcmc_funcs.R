########################################################################
## functions used in MCMC
########################################################################

## RJ reversibility term (and prior)
logProbChangeMod<-function(n.int,vars,I.vec,z.vec,p,vars.len,maxInt,miC){ 
  if(n.int==1){ #for acceptance ratio
    out<-log(I.vec[n.int+miC])-log(2*p*vars.len[vars]) + #proposal
      log(2*p*vars.len[vars])+log(maxInt) # prior
  } else{
    # perms<-permutations(n.int,n.int,vars)
    # sum.perm<-sum(apply(perms,1,function(row){1/prod(1-cumsum(z.vec[row][-n.int]))}))
    # lprob.vars.noReplace<-sum(log(z.vec[vars]))+log(sum.perm)
    #require(BiasedUrn) # this is much faster than above (esp for large maxInt)
    x<-rep(0,p)
    x[vars]<-1
    #lprob.vars.noReplace<-log(BiasedUrn::dMWNCHypergeo(x,rep(1,p),n.int,z.vec)) - do this in combination with imports: BiasedUrn in DESCRIPTION file, but that has a limit to MAXCOLORS
    lprob.vars.noReplace<-log(dMWNCHypergeo(x,rep(1,p),n.int,z.vec))
    out<-log(I.vec[n.int+miC])+lprob.vars.noReplace-n.int*log(2)-sum(log(vars.len[vars])) + # proposal
      +n.int*log(2)+sum(log(vars.len[vars]))+lchoose(p,n.int)+log(maxInt) # prior
  }
  return(out)
}

## RJ reversibility term (and prior) for categorical
logProbChangeModCat<-function(n.int,vars,I.vec,z.vec,p,nlevels,sub.size,maxInt,miC){
  if(n.int==1){ #for acceptance ratio
    out<-log(I.vec[n.int+miC])-log(p*(nlevels[vars]-1))-lchoose(nlevels[vars],sub.size[1:n.int]) + # proposal
      log(p*(nlevels[vars]-1))+lchoose(nlevels[vars],sub.size[1:n.int])+log(maxInt) # prior
  } else{
    x<-rep(0,p)
    x[vars]<-1
    lprob.vars.noReplace<-log(dMWNCHypergeo(x,rep(1,p),n.int,z.vec))
    out<-log(I.vec[n.int+miC])+lprob.vars.noReplace-n.int*sum(log(nlevels[vars]-1))-sum(lchoose(nlevels[vars],sub.size[1:n.int])) + # proposal
      n.int*sum(log(nlevels[vars]-1))+sum(lchoose(nlevels[vars],sub.size[1:n.int]))+lchoose(p,n.int)+log(maxInt) # prior
  }
  if(length(out)>1)
    browser()
  if(is.na(out))
    browser()
  return(out)
}


## log posterior
lp<-function(curr,prior,data){ 
  tt<-(
    - (curr$s2.rate+prior$g2)/curr$s2
    -(data$n/2+1+(curr$nbasis+1)/2 +prior$g1)*log(curr$s2) # changed -g1 to +g1
    + sum(log(abs(diag(curr$R)))) # .5*determinant of XtX
    + (prior$a.beta.prec+(curr$nbasis+1)/2-1)*log(curr$beta.prec) - prior$b.beta.prec*curr$beta.prec
    - (curr$nbasis+1)/2*log(2*pi)
    + (prior$h1+curr$nbasis-1)*log(curr$lam) - curr$lam*(prior$h2+1) # curr$nbasis-1 because poisson prior is excluding intercept (for curr$nbasis instead of curr$nbasis+1)
    #-lfactorial(curr$nbasis) # added, but maybe cancels with prior
  )
  if(curr$nbasis==0){
    return(tt)
  }
  #priors for basis parameters
  if(F){#(data$des){ # should these be involved in tempering??
    tt<-tt+(
      - sum(curr$n.int.des)*log(2) # signs for each basis function
      - sum(lchoose(data$pdes,curr$n.int.des)) # variables for each basis function
      - sum(log(data$vars.len.des[na.omit(c(curr$vars.des))])) # knots for each basis function
      - curr$nbasis*log(prior$maxInt.des) # degree of interaction for each basis function
    )
  }
  if(F){#(data$cat){
    tt<-tt+(
      - sum(sapply(1:curr$nbasis,function(i) curr$n.int.cat[i]*sum(log(data$nlevels[na.omit(curr$vars.cat[i,])]-1))))
      - sum(sapply(1:curr$nbasis,function(i) sum(lchoose(data$nlevels[na.omit(c(curr$vars.cat[i,]))],curr$sub.size[i,1:curr$n.int.cat[i]]))))
      - sum(lchoose(data$pcat,curr$n.int.cat))
      - curr$nbasis*log(prior$maxInt.cat)
    )
  }
  if(F){#(data$func){
    tt<-tt+(
      - sum(curr$n.int.func)*log(2)
      - sum(lchoose(data$pfunc,curr$n.int.func))
      - sum(log(data$vars.len.func[na.omit(c(curr$vars.func))]))
      - curr$nbasis*log(prior$maxInt.func)
    )
  }
  
  return(tt)
}


## get quadratic form that shows up in RJ acceptance probability
getQf<-function(XtX,Xty){ 
  R<-tryCatch(chol(XtX), error=function(e) matrix(F))
  if(R[1,1]){
    dr<-diag(R)
    if(length(dr)>1){
      if(max(dr[-1])/min(dr)>1e3) # TODO: this is a hack, otherwise we get huge variance inflation in beta
        return(NULL)
    }
    bhat<-backsolve(R,forwardsolve(R,Xty,transpose=T,upper.tri=T))
    qf<-crossprod(bhat,Xty)# same as sum((R%*%bhat)^2)
    return(list(R=R,bhat=bhat,qf=qf))
  } else{
    return(NULL)
  }
}

