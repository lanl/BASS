#######################################################
# Author: Devin Francom, Los Alamos National Laboratory
# Protected under GPL-3 license
# Los Alamos Computer Code release C19031
# github.com/lanl/BASS
#######################################################

########################################################################
## perform RJMCMC step (birth, death, or change)
########################################################################
birth_des_cat<-function(curr,prior,data){

  cand.des<-genCandBasis(minInt=prior$minInt,maxInt=prior$maxInt.des,I.vec=curr$I.vec.des,z.vec=curr$z.vec.des,p=data$pdes,xxt=data$xxt.des,q=prior$q,xx.unique.ind=data$unique.ind.des,vars.len=data$vars.len.des,prior)

  cand.cat<-genCandBasisCat(minInt=prior$minInt,maxInt=prior$maxInt.cat,I.vec=curr$I.vec.cat,z.vec=curr$z.vec.cat,p=data$pcat,xx=data$xx.cat,nlevels=data$nlevels,levels=data$levels,prior)

  if(cand.des$n.int + cand.cat$n.int == 0) # intercept
    return(curr)

  dc<-cand.des$basis*cand.cat$basis

  if(sum(dc!=0)<prior$npart.des){
    return(curr)
  }

  ata<-crossprod(dc)
  Xta<-crossprod(curr$dc.basis,dc)
  aty<-crossprod(dc,data$y)

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
    - cand.des$lbmcmp - cand.cat$lbmcmp
    + .5*log(curr$beta.prec+prior$beta.jprior.ind) - .5*log(1+curr$beta.prec)
    + prior$beta.jprior.ind*(
      .5*log(2*pi)#*curr$s2)
      + .5*sum(log(diag(qf.cand.list$R)))
      -.5*sum(log(diag(curr$R)))
    )
    )


  ## assign new values
  if(log(runif(1)) < alpha){
    curr<-addBasis(curr,cand.des,qf.cand.list,prior)
    curr<-addBasisDes(curr,cand.des,qf.cand.list,prior)
    curr<-addBasisCat(curr,cand.cat,qf.cand.list,prior)
    curr<-addBasisDC(curr,dc)
    # if type has cat and des, want to update curr$dc.basis also
  }
  return(curr)
}


death_des_cat<-function(curr,prior,data){
  basis<-sample(1:curr$nbasis,size=1)
  ind<-(1:curr$nc)[-(basis+1)]

  qf.cand.list<-getQf(curr$XtX[ind,ind],curr$Xty[ind])
  fullRank<-!is.null(qf.cand.list$qf)
  if(!fullRank){
    return(curr) # TODO: not sure why I need this, I shouldn't need it in theory
  }
  I.star.des<-curr$I.star.des
  I.star.des[curr$n.int.des[basis]+1]<-I.star.des[curr$n.int.des[basis]+1]-1
  I.vec.des<-I.star.des/sum(I.star.des)
  z.star.des<-curr$z.star.des
  z.star.des[curr$vars.des[basis,1:curr$n.int.des[basis]]]<-z.star.des[curr$vars.des[basis,1:curr$n.int.des[basis]]]-1
  z.vec.des<-z.star.des/sum(z.star.des)

  I.star.cat<-curr$I.star.cat
  I.star.cat[curr$n.int.cat[basis]+1]<-I.star.cat[curr$n.int.cat[basis]+1]-1
  I.vec.cat<-I.star.cat/sum(I.star.cat)
  z.star.cat<-curr$z.star.cat
  z.star.cat[curr$vars.cat[basis,1:curr$n.int.cat[basis]]]<-z.star.cat[curr$vars.cat[basis,1:curr$n.int.cat[basis]]]-1
  z.vec.cat<-z.star.cat/sum(z.star.cat)


  lpbmcmp<-0
  if(curr$n.int.des[basis]>0){
    lpbmcmp<-lpbmcmp+logProbChangeMod(curr$n.int.des[basis],curr$vars.des[basis,1:curr$n.int.des[basis]],I.vec.des,z.vec.des,data$pdes,data$vars.len.des,prior$maxInt.des,prior$miC)
  }

  if(curr$n.int.cat[basis]>0){#n.int,vars,I.vec,z.vec,p,nlevels,sub.size,maxInt
    lpbmcmp<-lpbmcmp+logProbChangeModCat(curr$n.int.cat[basis],curr$vars.cat[basis,1:curr$n.int.cat[basis]],I.vec.cat,z.vec.cat,data$pcat,data$nlevels,curr$sub.size[basis,],prior$maxInt.cat,prior$miC)
  }

  #if(is.na(lpbmcmp))
    #browser()

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
    curr<-deleteBasis(curr,basis,ind,qf.cand.list,I.star.des,I.vec.des,z.star.des,z.vec.des)
    curr<-deleteBasisDes(curr,basis,ind,qf.cand.list,I.star.des,I.vec.des,z.star.des,z.vec.des)
    curr<-deleteBasisCat(curr,basis,ind,qf.cand.list,I.star.cat,I.vec.cat,z.star.cat,z.vec.cat)
    curr<-deleteBasisDC(curr,basis)
  }
  return(curr)
}


change_des_cat<-function(curr,prior,data){
  basis<-sample(1:curr$nbasis,size=1)


  type.change<-sample(c('des','cat'),size=1,prob=c(curr$n.int.des[basis],curr$n.int.cat[basis]))

  if(type.change=='des'){
    int.change<-sample(1:(curr$n.int.des[basis]),size=1)
    use<-1:curr$n.int.des[basis]
    cand.des<-genBasisChange(curr,basis,int.change,data$xxt.des,prior$q,knots=curr$knots.des[basis,use],knotInd=curr$knotInd.des[basis,use],signs=curr$signs.des[basis,use],vars=curr$vars.des[basis,use],xx.unique.ind=data$unique.ind.des)
    dc<-cand.des$basis*curr$cat.basis[,basis+1]
  } else{
    int.change<-sample(1:(curr$n.int.cat[basis]),size=1)
    use<-1:curr$n.int.cat[basis]
    #curr,basis,int.change,xx,nlevels,levels,sub.size,sub,vars
    cand.cat<-genBasisChangeCat(curr,basis,int.change,data$xx.cat,data$nlevels,data$levels,curr$sub.size[basis,use],curr$sub.list[[basis]],vars=curr$vars.cat[basis,use])
    dc<-cand.cat$basis*curr$des.basis[,basis+1]
  }

  if(sum(dc!=0)<prior$npart.des){
    return(curr)
  }

  XtX.cand<-curr$XtX[1:curr$nc,1:curr$nc]
  XtX.cand[basis+1,]<-XtX.cand[,basis+1]<-crossprod(curr$dc.basis,dc)
  XtX.cand[basis+1,basis+1]<-crossprod(dc)
  Xty.cand<-curr$Xty[1:curr$nc]
  Xty.cand[basis+1]<-crossprod(dc,data$y)

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
    curr<-changeBasis(curr,cand.des,basis,qf.cand.list,XtX.cand,Xty.cand)
    curr<-changeBasisDC(curr,dc,basis)
    if(type.change=='des')
      curr<-changeBasisDes(curr,cand.des,basis,qf.cand.list,XtX.cand,Xty.cand)
    if(type.change=='cat')
      curr<-changeBasisCat(curr,cand.cat,basis,qf.cand.list,XtX.cand,Xty.cand)
  }
  return(curr)
}







