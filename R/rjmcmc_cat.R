#######################################################
# Author: Devin Francom, Los Alamos National Laboratory
# Protected under GPL-3 license
# Los Alamos Computer Code release C19031
# github.com/lanl/BASS
#######################################################

########################################################################
## perform RJMCMC step (birth, death, or change)
########################################################################
birth_cat<-function(curr,prior,data){
  cand.cat<-genCandBasisCat(minInt=prior$minInt,maxInt=prior$maxInt.cat,I.vec=curr$I.vec.cat,z.vec=curr$z.vec.cat,p=data$pcat,xx=data$xx.cat,nlevels=data$nlevels,levels=data$levels,prior)

  if(sum(cand.cat$basis!=0)<prior$npart.des){
    return(curr)
  }

  ata<-crossprod(cand.cat$basis)
  Xta<-crossprod(curr$cat.basis,cand.cat$basis)
  aty<-crossprod(cand.cat$basis,data$y)

  curr$Xty[curr$nc+1]<-aty
  curr$XtX[1:curr$nc,curr$nc+1]<-Xta
  curr$XtX[curr$nc+1,curr$nc+1]<-ata

  qf.cand.list<-getQf(curr$XtX[1:(curr$nc+1),1:(curr$nc+1)],curr$Xty[1:(curr$nc+1)])

  fullRank<-!is.null(qf.cand.list$qf)
  if(!fullRank){
    return(curr)
  }

  ## calculate log acceptance probability
  alpha<- data$itemp.ladder[curr$temp.ind]*(
    .5/curr$s2*(qf.cand.list$qf-curr$qf)/(1+curr$beta.prec)
    + log(curr$lam) - log(curr$nc) + log(data$death.prob.next/data$birth.prob)
    - cand.cat$lbmcmp + .5*log(curr$beta.prec+prior$beta.jprior.ind) - .5*log(1+curr$beta.prec)
    + prior$beta.jprior.ind*(
      .5*log(2*pi)#*curr$s2)
      + .5*sum(log(diag(qf.cand.list$R)))
      -.5*sum(log(diag(curr$R)))
    )
    )

  ## assign new values
  if(log(runif(1)) < alpha){
    curr<-addBasis(curr,cand.cat,qf.cand.list,prior)
    curr<-addBasisCat(curr,cand.cat,qf.cand.list,prior)
  }
  return(curr)
}


death_cat<-function(curr,prior,data){
  basis<-sample(1:curr$nbasis,size=1)
  ind<-(1:curr$nc)[-(basis+1)]

  qf.cand.list<-getQf(curr$XtX[ind,ind],curr$Xty[ind])
  fullRank<-!is.null(qf.cand.list$qf)
  if(!fullRank){
    return(curr) # TODO: not sure why I need this, I shouldn't need it in theory
  }

  I.star.cat<-curr$I.star.cat
  I.star.cat[curr$n.int.cat[basis]]<-I.star.cat[curr$n.int.cat[basis]]-1
  I.vec.cat<-I.star.cat/sum(I.star.cat)
  z.star.cat<-curr$z.star.cat
  z.star.cat[curr$vars.cat[basis,1:curr$n.int.cat[basis]]]<-z.star.cat[curr$vars.cat[basis,1:curr$n.int.cat[basis]]]-1
  z.vec.cat<-z.star.cat/sum(z.star.cat)

  lpbmcmp<-logProbChangeModCat(curr$n.int.cat[basis],curr$vars.cat[basis,1:curr$n.int.cat[basis]],I.vec.cat,z.vec.cat,data$pcat,data$nlevels,curr$sub.size[basis,],prior$maxInt.cat,prior$miC)

  # calculate log acceptance probability
  alpha<- data$itemp.ladder[curr$temp.ind]*(
    .5/curr$s2*(qf.cand.list$qf-curr$qf)/(1+curr$beta.prec)
    - log(curr$lam) + log(data$birth.prob.last/data$death.prob)
    + log(curr$nbasis) + lpbmcmp
    - .5*log(curr$beta.prec+prior$beta.jprior.ind) + .5*log(1+curr$beta.prec)
    + prior$beta.jprior.ind*(
      -.5*log(2*pi)#*curr$s2)
      +.5*sum(log(diag(qf.cand.list$R)))
      -.5*sum(log(diag(curr$R)))
    )
    )

  if(log(runif(1)) < alpha){
    curr<-deleteBasis(curr,basis,ind,qf.cand.list,I.star.cat,I.vec.cat,z.star.cat,z.vec.cat)
    curr<-deleteBasisCat(curr,basis,ind,qf.cand.list,I.star.cat,I.vec.cat,z.star.cat,z.vec.cat)
  }
  return(curr)
}


change_cat<-function(curr,prior,data){
  basis<-sample(1:curr$nbasis,size=1)
  int.change<-sample(1:(curr$n.int.cat[basis]),size=1)
  use<-1:curr$n.int.cat[basis]
  cand.cat<-genBasisChangeCat(curr,basis,int.change,data$xx.cat,data$nlevels,data$levels,curr$sub.size[basis,use],curr$sub.list[[basis]],vars=curr$vars.cat[basis,use])

  if(sum(cand.cat$basis!=0)<prior$npart.des){
    return(curr)
  }

  XtX.cand<-curr$XtX[1:curr$nc,1:curr$nc]
  XtX.cand[basis+1,]<-XtX.cand[,basis+1]<-crossprod(curr$cat.basis,cand.cat$basis)
  XtX.cand[basis+1,basis+1]<-crossprod(cand.cat$basis)
  Xty.cand<-curr$Xty[1:curr$nc]
  Xty.cand[basis+1]<-crossprod(cand.cat$basis,data$y)

  qf.cand.list<-getQf(XtX.cand,Xty.cand)

  fullRank<-!is.null(qf.cand.list$qf)
  if(!fullRank){
    return(curr)
  }

  alpha<-data$itemp.ladder[curr$temp.ind]*(
    .5/curr$s2*(qf.cand.list$qf-curr$qf)/(1+curr$beta.prec)
    + prior$beta.jprior.ind*(
      .5*sum(log(diag(qf.cand.list$R)))-.5*sum(log(diag(curr$R)))
    )
  )

  if(log(runif(1))<alpha){
    curr<-changeBasis(curr,cand.cat,basis,qf.cand.list,XtX.cand,Xty.cand)
    curr<-changeBasisCat(curr,cand.cat,basis,qf.cand.list,XtX.cand,Xty.cand)
  }
  return(curr)
}







