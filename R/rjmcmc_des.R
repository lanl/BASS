#######################################################
# Author: Devin Francom, Los Alamos National Laboratory
# Protected under GPL-3 license
# Los Alamos Computer Code release C19031
# github.com/lanl/BASS
#######################################################

########################################################################
## perform RJMCMC step (birth, death, or change)
########################################################################
birth_des<-function(curr,prior,data){
  cand.des<-genCandBasis(minInt=prior$minInt,maxInt=prior$maxInt.des,I.vec=curr$I.vec.des,z.vec=curr$z.vec.des,p=data$pdes,xxt=data$xxt.des,q=prior$q,xx.unique.ind=data$unique.ind.des,vars.len=data$vars.len.des,prior)

  if(sum(cand.des$basis!=0)<prior$npart.des){
    return(curr)
  }

  ata<-crossprod(cand.des$basis)
  Xta<-crossprod(curr$des.basis,cand.des$basis)
  aty<-crossprod(cand.des$basis,data$y)

  curr$Xty[curr$nc+1]<-aty
  curr$XtX[1:curr$nc,curr$nc+1]<-Xta
  curr$XtX[curr$nc+1,curr$nc+1]<-ata

  qf.cand.list<-getQf(curr$XtX[1:(curr$nc+1),1:(curr$nc+1)],curr$Xty[1:(curr$nc+1)])

  fullRank<-!is.null(qf.cand.list$qf)
  if(!fullRank){
    return(curr)
  }


  ## calculate log acceptance probability
  alpha<- data$itemp.ladder[curr$temp.ind]*(.5/curr$s2*(qf.cand.list$qf-curr$qf)/(1+curr$beta.prec) + log(curr$lam) - log(curr$nc) + log(data$death.prob.next/data$birth.prob) - cand.des$lbmcmp + .5*log(curr$beta.prec) - .5*log(1+curr$beta.prec))
  #cat(- cand.des$lbmcmp,' ')

  ## assign new values
  if(log(runif(1)) < alpha){
    curr<-addBasis(curr,cand.des,qf.cand.list,prior)
    curr<-addBasisDes(curr,cand.des,qf.cand.list,prior)
    # if type has cat and des, want to update curr$dc.basis also
  }
  return(curr)
}


death_des<-function(curr,prior,data){
  basis<-sample(1:curr$nbasis,size=1)
  ind<-(1:curr$nc)[-(basis+1)]

  qf.cand.list<-getQf(curr$XtX[ind,ind],curr$Xty[ind])
  fullRank<-!is.null(qf.cand.list$qf)
  if(!fullRank){
    return(curr) # TODO: not sure why I need this, I shouldn't need it in theory
  }
  I.star.des<-curr$I.star.des
  I.star.des[curr$n.int.des[basis]]<-I.star.des[curr$n.int.des[basis]]-1
  I.vec.des<-I.star.des/sum(I.star.des)
  z.star.des<-curr$z.star.des
  z.star.des[curr$vars.des[basis,1:curr$n.int.des[basis]]]<-z.star.des[curr$vars.des[basis,1:curr$n.int.des[basis]]]-1
  z.vec.des<-z.star.des/sum(z.star.des)

  lpbmcmp<-logProbChangeMod(curr$n.int.des[basis],curr$vars.des[basis,1:curr$n.int.des[basis]],I.vec.des,z.vec.des,data$pdes,data$vars.len.des,prior$maxInt.des,prior$miC)

  # calculate log acceptance probability
  alpha<- data$itemp.ladder[curr$temp.ind]*(.5/curr$s2*(qf.cand.list$qf-curr$qf)/(1+curr$beta.prec) - log(curr$lam) + log(data$birth.prob.last/data$death.prob) + log(curr$nbasis) + lpbmcmp - .5*log(curr$beta.prec) + .5*log(1+curr$beta.prec))

  if(log(runif(1)) < alpha){
    curr<-deleteBasis(curr,basis,ind,qf.cand.list,I.star.des,I.vec.des,z.star.des,z.vec.des)
    curr<-deleteBasisDes(curr,basis,ind,qf.cand.list,I.star.des,I.vec.des,z.star.des,z.vec.des)
  }
  return(curr)
}


change_des<-function(curr,prior,data){
  basis<-sample(1:curr$nbasis,size=1)
  int.change<-sample(1:(curr$n.int.des[basis]),size=1)
  use<-1:curr$n.int.des[basis]
  cand.des<-genBasisChange(curr,basis,int.change,data$xxt.des,prior$q,knots=curr$knots.des[basis,use],knotInd=curr$knotInd.des[basis,use],signs=curr$signs.des[basis,use],vars=curr$vars.des[basis,use],xx.unique.ind=data$unique.ind.des)



  if(sum(cand.des$basis!=0)<prior$npart.des){
    return(curr)
  }

  XtX.cand<-curr$XtX[1:curr$nc,1:curr$nc]
  XtX.cand[basis+1,]<-XtX.cand[,basis+1]<-crossprod(curr$des.basis,cand.des$basis)
  XtX.cand[basis+1,basis+1]<-crossprod(cand.des$basis)
  Xty.cand<-curr$Xty[1:curr$nc]
  Xty.cand[basis+1]<-crossprod(cand.des$basis,data$y)

  qf.cand.list<-getQf(XtX.cand,Xty.cand)

  fullRank<-!is.null(qf.cand.list$qf)
  if(!fullRank){
    return(curr)
  }

  alpha<-data$itemp.ladder[curr$temp.ind]*.5/curr$s2*(qf.cand.list$qf-curr$qf)/(1+curr$beta.prec)

  if(log(runif(1))<alpha){
    curr<-changeBasis(curr,cand.des,basis,qf.cand.list,XtX.cand,Xty.cand)
    curr<-changeBasisDes(curr,cand.des,basis,qf.cand.list,XtX.cand,Xty.cand)
  }
  return(curr)
}







