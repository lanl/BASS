#######################################################
# Author: Devin Francom, Los Alamos National Laboratory
# Protected under GPL-3 license
# Los Alamos Computer Code release C19031
# github.com/lanl/BASS
#######################################################

########################################################################
## perform RJMCMC step (birth, death, or change)
########################################################################
birth_des_func<-function(curr,prior,data){

  cand.des<-genCandBasis(minInt=prior$minInt,maxInt=prior$maxInt.des,I.vec=curr$I.vec.des,z.vec=curr$z.vec.des,p=data$pdes,xxt=data$xxt.des,q=prior$q,xx.unique.ind=data$unique.ind.des,vars.len=data$vars.len.des,prior)
  if(sum(cand.des$basis!=0)<prior$npart.des)
      return(curr)

  cand.func<-genCandBasis(minInt=prior$minInt,maxInt=prior$maxInt.func,I.vec=curr$I.vec.func,z.vec=curr$z.vec.func,p=data$pfunc,xxt=data$xxt.func,q=prior$q,xx.unique.ind=data$unique.ind.func,vars.len=data$vars.len.func,prior)
  if(sum(cand.func$basis!=0)<prior$npart.func)
    return(curr)

  if(cand.des$n.int + cand.func$n.int == 0) # intercept
    return(curr)

  ata<-crossprod(cand.des$basis)*crossprod(cand.func$basis)
  Xta<-crossprod(curr$des.basis,cand.des$basis)*crossprod(curr$func.basis,cand.func$basis)
  aty<-tcrossprod(crossprod(cand.des$basis,data$y),cand.func$basis)

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
    - cand.des$lbmcmp - cand.func$lbmcmp
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
    curr<-addBasisFunc(curr,cand.func,qf.cand.list,prior)
    # if type has cat and des, want to update curr$dc.basis also
  }
  return(curr)
}


death_des_func<-function(curr,prior,data){
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
  if(curr$n.int.des[basis]>0)
    z.star.des[curr$vars.des[basis,1:curr$n.int.des[basis]]]<-z.star.des[curr$vars.des[basis,1:curr$n.int.des[basis]]]-1
  z.vec.des<-z.star.des/sum(z.star.des)

  I.star.func<-curr$I.star.func
  I.star.func[curr$n.int.func[basis]+1]<-I.star.func[curr$n.int.func[basis]+1]-1
  I.vec.func<-I.star.func/sum(I.star.func)
  z.star.func<-curr$z.star.func
  if(curr$n.int.func[basis]>0)
    z.star.func[curr$vars.func[basis,1:curr$n.int.func[basis]]]<-z.star.func[curr$vars.func[basis,1:curr$n.int.func[basis]]]-1
  z.vec.func<-z.star.func/sum(z.star.func)

  lpbmcmp<-0
  if(curr$n.int.des[basis]>0){
    lpbmcmp<-lpbmcmp+logProbChangeMod(curr$n.int.des[basis],curr$vars.des[basis,1:curr$n.int.des[basis]],I.vec.des,z.vec.des,data$pdes,data$vars.len.des,prior$maxInt.des,prior$miC)
  }

  if(curr$n.int.func[basis]>0){
    lpbmcmp<-lpbmcmp+logProbChangeMod(curr$n.int.func[basis],curr$vars.func[basis,1:curr$n.int.func[basis]],I.vec.func,z.vec.func,data$pfunc,data$vars.len.func,prior$maxInt.func,prior$miC)
  }

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
    curr<-deleteBasisFunc(curr,basis,ind,qf.cand.list,I.star.func,I.vec.func,z.star.func,z.vec.func)
  }
  return(curr)
}


change_des_func<-function(curr,prior,data){
  basis<-sample(1:curr$nbasis,size=1)
  type.change<-sample(c('des','func'),size=1,prob=c(curr$n.int.des[basis],curr$n.int.func[basis]))

  if(type.change=='des'){
    int.change<-sample(1:(curr$n.int.des[basis]),size=1)
    use<-1:curr$n.int.des[basis]
    cand.des<-genBasisChange(curr,basis,int.change,data$xxt.des,prior$q,knots=curr$knots.des[basis,use],knotInd=curr$knotInd.des[basis,use],signs=curr$signs.des[basis,use],vars=curr$vars.des[basis,use],xx.unique.ind=data$unique.ind.des)
    cand.func<-list(basis=curr$func.basis[,basis+1])
  } else{
    int.change<-sample(1:(curr$n.int.func[basis]),size=1)
    use<-1:curr$n.int.func[basis]
    cand.func<-genBasisChange(curr,basis,int.change,data$xxt.func,prior$q,knots=curr$knots.func[basis,use],knotInd=curr$knotInd.func[basis,use],signs=curr$signs.func[basis,use],vars=curr$vars.func[basis,use],xx.unique.ind=data$unique.ind.func)
    cand.des<-list(basis=curr$des.basis[,basis+1])
  }

  if(sum(cand.des$basis!=0)<prior$npart.des){
    return(curr)
  }
  if(sum(cand.func$basis!=0)<prior$npart.func){
    return(curr)
  }

  XtX.cand<-curr$XtX[1:curr$nc,1:curr$nc]
  XtX.cand[basis+1,]<-XtX.cand[,basis+1]<-crossprod(curr$des.basis,cand.des$basis)*crossprod(curr$func.basis,cand.func$basis)
  XtX.cand[basis+1,basis+1]<-crossprod(cand.des$basis)*crossprod(cand.func$basis)
  Xty.cand<-curr$Xty[1:curr$nc]
  Xty.cand[basis+1]<-tcrossprod(crossprod(cand.des$basis,data$y),cand.func$basis)

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
    if(type.change=='des')
      curr<-changeBasisDes(curr,cand.des,basis,qf.cand.list,XtX.cand,Xty.cand)
    if(type.change=='func')
      curr<-changeBasisFunc(curr,cand.func,basis,qf.cand.list,XtX.cand,Xty.cand)
  }
  return(curr)
}







