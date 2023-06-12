#######################################################
# Author: Devin Francom, Los Alamos National Laboratory
# Protected under GPL-3 license
# Los Alamos Computer Code release C19031
# github.com/lanl/BASS
#######################################################

########################################################################
## perform RJMCMC step (birth, death, or change)
########################################################################
birth_coinflip_des<-function(curr,prior,data){
  cand.des<-genCandBasisCoinflip(minInt=prior$minInt,maxExpectedInt=prior$maxExpectedInt.des,eta.vec=curr$eta.star.des,
                         nint.proposal=curr$nint.proposal,p=data$pdes,xxt=data$xxt.des,degree=prior$q,
                         xx.unique.ind=data$unique.ind.des,vars.len=data$vars.len.des,prior=prior,w3=curr$w3)

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
  alpha<- data$itemp.ladder[curr$temp.ind]*(
    .5/curr$s2*(qf.cand.list$qf-curr$qf)/(1+curr$beta.prec)
    + log(curr$lam) - log(curr$nc)
    + log(data$death.prob.next/data$birth.prob) - cand.des$lbmcmp
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
    # if type has cat and des, want to update curr$dc.basis also
  }
  return(curr)
}


death_coinflip_des<-function(curr,prior,data){
  basis<-sample(1:curr$nbasis,size=1)
  ind<-(1:curr$nc)[-(basis+1)]

  qf.cand.list<-getQf(curr$XtX[ind,ind],curr$Xty[ind])
  fullRank<-!is.null(qf.cand.list$qf)
  if(!fullRank){
    return(curr) # TODO: not sure why I need this, I shouldn't need it in theory
  }
  #KR:
  #I.star.des<-curr$I.star.des
  #I.star.des[curr$n.int.des[basis]]<-I.star.des[curr$n.int.des[basis]]-1
  #I.vec.des<-I.star.des/sum(I.star.des)
  #z.star.des<-curr$z.star.des
  #z.star.des[curr$vars.des[basis,1:curr$n.int.des[basis]]]<-z.star.des[curr$vars.des[basis,1:curr$n.int.des[basis]]]-1
  #z.vec.des<-z.star.des/sum(z.star.des)

  eta.star.des<-curr$eta.star.des
  eta.star.des[curr$n.int.des[basis]]<-eta.star.des[curr$n.int.des[basis]]-1
  vars.cand <- curr$vars.des[basis,1:curr$n.int.des[basis]]
  lprob_tmp <- 0
  for(jj in 1:prior$maxExpectedInt){
    wts.j <- makeCoinWeights(eta.star.des, jj, curr$w3)
    term_j <- log(curr$nint.proposal[jj]) + sum(log(wts.j[vars.cand])) + sum(log(1 - wts.j[-vars.cand]))
    lprob_tmp <- lprob_tmp + exp(term_j)
  }
  delayedRejectionTerm <- log(lprob_tmp)

  lpbmcmp<-logProbChangeModCoinflip(curr$n.int.des[basis],vars.cand,data$pdes,data$vars.len.des,prior$nint.prior,prior$miC,delayedRejectionTerm)

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
  }
  return(curr)
}


change_coinflip_des<-function(curr,prior,data){
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

  alpha<-data$itemp.ladder[curr$temp.ind]*(
    .5/curr$s2*(qf.cand.list$qf-curr$qf)/(1+curr$beta.prec)
    + prior$beta.jprior.ind*(
      .5*sum(log(diag(qf.cand.list$R)))-.5*sum(log(diag(curr$R)))
    )
  )

  if(log(runif(1))<alpha){
    curr<-changeBasis(curr,cand.des,basis,qf.cand.list,XtX.cand,Xty.cand)
    curr<-changeBasisDes(curr,cand.des,basis,qf.cand.list,XtX.cand,Xty.cand)
  }
  return(curr)
}







