#######################################################
# Author: Devin Francom, Los Alamos National Laboratory
# Protected under GPL-3 license
# Los Alamos Computer Code release C19031
# github.com/lanl/BASS
#######################################################

########################################################################
## generate a cadidate basis (birth, change)
########################################################################
genCandBasis<-function(minInt,maxInt,I.vec,z.vec,p,xxt,q,xx.unique.ind,vars.len,prior){
  # get number of variables in interaction
  n.int<-sample(minInt:maxInt,size=1,prob=I.vec)
  if(n.int==0)
    return(list(basis=rep(1,ncol(xxt)),n.int=n.int,lbmcmp=0))
    #return(NULL)
  # get signs, vars, knots
  signs<-sample(c(-1,1),size=n.int,replace=T)
  if(n.int==1){
    vars<-sample(1:p,size=1)
    #knotInd<-sample.int(vars.len[vars],size=1)
    knotInd<-sample(xx.unique.ind[[vars]],size=1)
  } else{
    vars<-sort(sample(1:p,size=n.int,prob=z.vec,replace=F))
    #knotInd<-sapply(vars.len[vars],sample.int,size=1)
    knotInd<-sapply(xx.unique.ind[vars],sample,size=1)
  }
  # make basis, get reversibility term
  #knots<-sapply(1:n.int,function(nn) xxt.unique[[vars[nn]]][knotInd[nn]])
  knots<-xxt[cbind(vars,knotInd)]
  basis<-makeBasis(signs,vars,knots,xxt,q)
  lpbmcmp<-logProbChangeMod(n.int,vars,I.vec,z.vec,p,vars.len,maxInt,prior$miC)

  return(list(basis=basis,n.int=n.int,signs=signs,vars=vars,knotInd=knotInd,knots=knots,lbmcmp=lpbmcmp))
}

genCandBasisCat<-function(minInt,maxInt,I.vec,z.vec,p,xx,nlevels,levels,prior){
  # get number of variables in interaction
  n.int<-sample(minInt:maxInt,size=1,prob=I.vec)
  if(n.int==0)
    return(list(basis=rep(1,nrow(xx)),n.int=n.int,lbmcmp=0,sub=list(NA)))
  # get vars, subsets
  if(n.int==1){
    vars<-sample(1:p,size=1)
  } else{
    vars<-sort(sample(1:p,size=n.int,prob=z.vec,replace=F))
  }
  sub.size<-NA # for each of vars, number of categories included in subset
  sub<-list() # actual subsets
  for(ii in 1:n.int){
    sub.size[ii]<-sample(1:(nlevels[vars[ii]]-1),size=1) # sample the size of the subset
    sub[[ii]]<-sample(levels[[vars[ii]]],size=sub.size[ii]) # sample the subset
  }
  # make basis, get reversibility term
  basis<-makeBasisCat(vars,sub,xx)
  lpbmcmp<-logProbChangeModCat(n.int,vars,I.vec,z.vec,p,nlevels,sub.size,maxInt,prior$miC)
  if(is.na(lpbmcmp))
    browser()
  return(list(basis=basis,n.int=n.int,vars=vars,sub.size=sub.size,sub=sub,lbmcmp=lpbmcmp))
}



genBasisChange<-function(curr,basis,int.change,xxt,q,knots,knotInd,signs,vars,xx.unique.ind){
  signs[int.change]<-sample(c(-1,1),size=1)
  #knotInd[int.change]<-sample.int(vars.len[vars[int.change]],size=1)
  knotInd[int.change]<-sample(xx.unique.ind[[vars[int.change]]],size=1)
  #knots[int.change]<-xxt.unique[[vars[int.change]]][knotInd[int.change]]
  knots[int.change]<-xxt[vars[int.change],knotInd[int.change]]

  basis<-makeBasis(signs,vars,knots,xxt,q)
  return(list(knots=knots,knotInd=knotInd,signs=signs,basis=basis))
}

genBasisChangeCat<-function(curr,basis,int.change,xx,nlevels,levels,sub.size,sub,vars){
  sub.size[int.change]<-sample(1:(nlevels[vars[int.change]]-1),size=1)
  sub[[int.change]]<-sample(levels[[vars[int.change]]],size=sub.size[int.change])

  basis<-makeBasisCat(vars,sub,xx)
  return(list(sub.size=sub.size,sub=sub,basis=basis))
}

########################################################################
## write to curr
########################################################################
addBasis<-function(curr,cand,qf.cand.list,prior){
  # basis characteristics
  curr$nbasis<-curr$nbasis+1
  curr$nc<-curr$nbasis+1

  # updates to quantities used elsewhere (XtX & Xty are already updated)
  curr$qf<-qf.cand.list$qf
  curr$bhat<-qf.cand.list$bhat
  curr$R<-qf.cand.list$R
  curr$R.inv.t<-backsolve(curr$R,diag(curr$nc))

  # diagnostics
  curr$cmod<-T
  curr$step<-1
  curr$count[1]<-curr$count[1]+1

  return(curr)
}
addBasisDes<-function(curr,cand,qf.cand.list,prior){
  curr$n.int.des[curr$nbasis]<-cand$n.int
  fill<-rep(NA,prior$maxInt.des-cand$n.int)
  curr$knots.des<-rbind(curr$knots.des,c(cand$knots,fill))
  curr$knotInd.des<-rbind(curr$knotInd.des,c(cand$knotInd,fill))
  curr$signs.des<-rbind(curr$signs.des,c(cand$signs,fill))
  curr$vars.des<-rbind(curr$vars.des,c(cand$vars,fill))

  curr$I.star.des[cand$n.int+prior$miC]<-curr$I.star.des[cand$n.int+prior$miC]+1
  curr$I.vec.des<-curr$I.star.des/sum(curr$I.star.des)
  curr$z.star.des[cand$vars]<-curr$z.star.des[cand$vars]+1
  curr$z.vec.des<-curr$z.star.des/sum(curr$z.star.des)

  # basis functions
  curr$des.basis<-cbind(curr$des.basis,cand$basis)
  return(curr)
}
addBasisFunc<-function(curr,cand,qf.cand.list,prior){
  curr$n.int.func[curr$nbasis]<-cand$n.int
  fill<-rep(NA,prior$maxInt.func-cand$n.int)
  curr$knots.func<-rbind(curr$knots.func,c(cand$knots,fill))
  curr$knotInd.func<-rbind(curr$knotInd.func,c(cand$knotInd,fill))
  curr$signs.func<-rbind(curr$signs.func,c(cand$signs,fill))
  curr$vars.func<-rbind(curr$vars.func,c(cand$vars,fill))

  curr$I.star.func[cand$n.int+prior$miC]<-curr$I.star.func[cand$n.int+prior$miC]+1
  curr$I.vec.func<-curr$I.star.func/sum(curr$I.star.func)
  curr$z.star.func[cand$vars]<-curr$z.star.func[cand$vars]+1
  curr$z.vec.func<-curr$z.star.func/sum(curr$z.star.func)

  # basis functions
  curr$func.basis<-cbind(curr$func.basis,cand$basis)
  return(curr)
}
addBasisCat<-function(curr,cand,qf.cand.list,prior){
  curr$n.int.cat[curr$nbasis]<-cand$n.int
  fill<-rep(NA,prior$maxInt.cat-cand$n.int)
  curr$sub.size<-rbind(curr$sub.size,c(cand$sub.size,fill))
  #browser()
  curr$sub.list[[curr$nbasis]]<-cand$sub #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ maybe should use a design matrix style for factors instead of this list...wouldn't need to use %in%
  curr$vars.cat<-rbind(curr$vars.cat,c(cand$vars,fill))

  curr$I.star.cat[cand$n.int+prior$miC]<-curr$I.star.cat[cand$n.int+prior$miC]+1
  curr$I.vec.cat<-curr$I.star.cat/sum(curr$I.star.cat)
  curr$z.star.cat[cand$vars]<-curr$z.star.cat[cand$vars]+1
  curr$z.vec.cat<-curr$z.star.cat/sum(curr$z.star.cat)

  # basis functions
  curr$cat.basis<-cbind(curr$cat.basis,cand$basis)
  return(curr)
}

addBasisDC<-function(curr,dc){
  curr$dc.basis<-cbind(curr$dc.basis,dc)
  return(curr)
}

deleteBasis<-function(curr,basis,ind,qf.cand.list,I.star,I.vec,z.star,z.vec){
  # basis characteristics
  curr$nbasis<-curr$nbasis-1
  curr$nc<-curr$nbasis+1

  # updates to quantities used elsewhere
  curr$Xty[1:curr$nc]<-curr$Xty[ind]
  curr$XtX[1:curr$nc,1:curr$nc]<-curr$XtX[ind,ind]
  curr$qf<-qf.cand.list$qf
  curr$bhat<-qf.cand.list$bhat
  curr$R<-qf.cand.list$R
  curr$R.inv.t<-backsolve(curr$R,diag(curr$nc))

  # diagnostics
  curr$cmod<-T
  curr$step<-2
  curr$count[2]<-curr$count[2]+1

  return(curr)
}

deleteBasisDes<-function(curr,basis,ind,qf.cand.list,I.star,I.vec,z.star,z.vec){
  curr$n.int.des<-curr$n.int.des[-basis]
  curr$knots.des<-curr$knots.des[-basis,,drop=F]
  curr$knotInd.des<-curr$knotInd.des[-basis,,drop=F]
  curr$signs.des<-curr$signs.des[-basis,,drop=F]
  curr$vars.des<-curr$vars.des[-basis,,drop=F]
  curr$I.star.des<-I.star
  curr$I.vec.des<-I.vec
  curr$z.star.des<-z.star
  curr$z.vec.des<-z.vec
  # basis functions
  curr$des.basis<-curr$des.basis[,-(basis+1),drop=F]
  return(curr)
}
deleteBasisFunc<-function(curr,basis,ind,qf.cand.list,I.star,I.vec,z.star,z.vec){
  curr$n.int.func<-curr$n.int.func[-basis]
  curr$knots.func<-curr$knots.func[-basis,,drop=F]
  curr$knotInd.func<-curr$knotInd.func[-basis,,drop=F]
  curr$signs.func<-curr$signs.func[-basis,,drop=F]
  curr$vars.func<-curr$vars.func[-basis,,drop=F]
  curr$I.star.func<-I.star
  curr$I.vec.func<-I.vec
  curr$z.star.func<-z.star
  curr$z.vec.func<-z.vec
  # basis functions
  curr$func.basis<-curr$func.basis[,-(basis+1),drop=F]
  return(curr)
}
deleteBasisCat<-function(curr,basis,ind,qf.cand.list,I.star,I.vec,z.star,z.vec){
  curr$n.int.cat<-curr$n.int.cat[-basis]
  curr$sub.size<-curr$sub.size[-basis,,drop=F]
  curr$sub.list[[basis]]<-NULL
  curr$vars.cat<-curr$vars.cat[-basis,,drop=F]
  curr$I.star.cat<-I.star
  curr$I.vec.cat<-I.vec
  curr$z.star.cat<-z.star
  curr$z.vec.cat<-z.vec
  # basis functions
  curr$cat.basis<-curr$cat.basis[,-(basis+1),drop=F]
  return(curr)
}
deleteBasisDC<-function(curr,basis.ind){
  curr$dc.basis<-curr$dc.basis[,-(basis.ind+1),drop=F]
  return(curr)
}

changeBasis<-function(curr,cand,basis,qf.cand.list,XtX.cand,Xty.cand){
  # updates to quantities used elsewhere
  curr$Xty[basis+1]<-Xty.cand[basis+1]
  curr$XtX[1:curr$nc,1:curr$nc]<-XtX.cand
  curr$qf<-qf.cand.list$qf
  curr$bhat<-qf.cand.list$bhat
  curr$R<-qf.cand.list$R
  curr$R.inv.t<-backsolve(curr$R,diag(curr$nc))

  # diagnostics
  curr$cmod<-T
  curr$step<-3
  curr$count[3]<-curr$count[3]+1

  return(curr)
}
changeBasisDes<-function(curr,cand,basis,qf.cand.list,XtX.cand,Xty.cand){
  # basis characteristics
  curr$knots.des[basis,1:curr$n.int.des[basis]]<-cand$knots
  curr$knotInd.des[basis,1:curr$n.int.des[basis]]<-cand$knotInd
  curr$signs.des[basis,1:curr$n.int.des[basis]]<-cand$signs
  # basis functions
  curr$des.basis[,basis+1]<-cand$basis
  return(curr)
}
changeBasisFunc<-function(curr,cand,basis,qf.cand.list,XtX.cand,Xty.cand){
  # basis characteristics
  curr$knots.func[basis,1:curr$n.int.func[basis]]<-cand$knots
  curr$knotInd.func[basis,1:curr$n.int.func[basis]]<-cand$knotInd
  curr$signs.func[basis,1:curr$n.int.func[basis]]<-cand$signs
  # basis functions
  curr$func.basis[,basis+1]<-cand$basis
  return(curr)
}
changeBasisCat<-function(curr,cand,basis,qf.cand.list,XtX.cand,Xty.cand){
  # basis characteristics
  curr$sub.size[basis,1:curr$n.int.cat[basis]]<-cand$sub.size
  curr$sub.list[[basis]]<-cand$sub
  # basis functions
  curr$cat.basis[,basis+1]<-cand$basis
  return(curr)
}
changeBasisDC<-function(curr,cand,basis){
  curr$dc.basis[,basis+1]<-cand
  return(curr)
}


# check monotonicity of functional output - for now only with des (no cat) and degree=1
checkMono<-function(tdes.basis.ext,func.basis,bhat){

  if(length(bhat)>1)
    out<-func.basis%*%diag(c(bhat))%*%tdes.basis.ext
  else
    out<-func.basis%*%bhat%*%tdes.basis.ext

  #browser()
  min(diff(out)) >= 0
}
# checkMono<-function(cand.des,cand.func,curr,data){
#
#   basis<-makeBasis(cand.des$signs,cand.des$vars,cand.des$knots,data$ext,1)
#
#   tdes.basis.ext<-rbind(makeBasisMatrixCurr(curr,data),makeBasis(cand.des$signs,cand.des$vars,cand.des$knots,data$ext,1))
#   out<-cbind(curr$func.basis,cand.func$basis)%*%diag(curr$bhat)%*%tdes.basis.ext
#   min(apply(out,2,function(x) min(diff(x)))) < 0
# }

checkMono2<-function(curr,data,k){
  # make des.basis.ext (n.ext x nbasis) using data$ext
  # make func.deriv (nfunc x nbasis) matrix
  # func.deriv%*%diag(beta)%*%t(des.basis.ext)

  # shortcut: use regular func.basis, get differences
  tdes.basis.ext<-makeBasisMatrixCurr(curr,data)

  if(length(curr$beta)>1){
    out<-curr$func.basis%*%diag(c(curr$beta))%*%tdes.basis.ext
  } else{
    out<-curr$func.basis%*%curr$beta%*%tdes.basis.ext
  }

  if(k>10)
    browser()

  min(diff(out)) >= 0


}
makeBasisMatrixCurr<-function(curr,data){#i,nbasis,vars,signs,knots.ind,q,xxt,n.int,xx.train){
  xxt.ext<-data$ext
  n<-ncol(xxt.ext)
  tbasis.mat<-matrix(nrow=curr$nbasis+1,ncol=n)
  tbasis.mat[1,]<-1
  if(curr$nbasis>0){
    for(m in 1:curr$nbasis){
      if(curr$n.int.des[m]==0){
        tbasis.mat[m+1,]<-1
      } else{
        use<-1:curr$n.int.des[m]
        #knots<-xx.train[cbind(knots.ind[i,m,use],vars[i,m,use])] # get knots from knots.ind
        tbasis.mat[m+1,]<-makeBasis(curr$signs.des[m,use],curr$vars.des[m,use],curr$knots.des[m,use],xxt.ext,1)
      }
    }
  }
  return(tbasis.mat)
}
