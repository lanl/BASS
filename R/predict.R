#######################################################
# Author: Devin Francom, Los Alamos National Laboratory
# Protected under GPL-3 license
# Los Alamos Computer Code release C19031
# github.com/lanl/BASS
#######################################################

###############################################################
## predict methods
###############################################################
scale_range_mat<-function(x,r){
  #sweep(sweep(x,2,r[1,]),2,r[2,]-r[1,],FUN='/')
  t((t(x)-r[1,])/c(diff(r)))
  #(x - matrix(r[1,], dim(x)[1], dim(x)[2], byrow = TRUE))/matrix(diff(r), dim(x)[1], dim(x)[2], byrow = TRUE)
}

#' @title BASS Prediction
#'
#' @description Predict function for BASS.  Outputs the posterior predictive samples based on the specified MCMC iterations.
#' @param object a fitted model, output from the \code{bass} function.
#' @param newdata a matrix of new input values at which to predict.  The columns should correspond to the same variables used in the \code{bass} function.
#' @param newdata.func a matrix of new values of the functional variable.  If none, the same values will be used as in the training data.
#' @param mcmc.use a vector indexing which MCMC iterations to use for prediction.
#' @param verbose logical; should progress be displayed?
#' @param nugget logical; should predictions include error? If FALSE, predictions will be for mean.
#' @param ... further arguments passed to or from other methods.
#' @details Efficiently predicts when two MCMC iterations have the same basis functions (but different weights).
#' @return If model output is a scalar, this returns a matrix with the same number of rows as \code{newdata} and columns corresponding to the the MCMC iterations \code{mcmc.use}.  These are samples from the posterior predictive distribution.  If model output is functional, this returns an array with first dimension corresponding to MCMC iteration, second dimension corresponding to the rows of \code{newdata}, and third dimension corresponding to the rows of \code{newdata.func}.
#' @seealso \link{bass} for model fitting and \link{sobol} for sensitivity analysis.
#' @export
#' @examples
#' # See examples in bass documentation.
#'
predict.bass<-function(object,newdata,newdata.func=NULL,mcmc.use=NULL,verbose=FALSE,nugget=FALSE,...){
  if(is.null(mcmc.use)){ # if null, use all
    mcmc.use<-1:((object$nmcmc-object$nburn)/object$thin)
  }
  if(object$func){
    if(is.null(newdata.func))
      newdata.func<-object$xx.func
    else{
      dxf<-dim(newdata.func)
      if(is.null(dxf))
        newdata.func<-matrix(newdata.func)
      for(i in 1:ncol(newdata.func)){
        newdata.func[,i]<-scale_range(newdata.func[,i],object$range.func[,i])
      }
    }
  } else{
    newdata.func<-t(1) # placeholder
  }
  tnewdata.func<-t(newdata.func)

  dx<-dim(newdata)
  if(is.null(dx)){
    newdata<-data.frame(newdata)
    dx<-dim(newdata)
  }
  pd<-sum(object$pdes)+sum(object$pcat)
  if(dx[2]!=pd){
    newdata<-t(newdata)
    dx<-dim(newdata)
    if(dx[2]!=pd)
      stop('number of variables in newdata does not match number used in object')
  }
  newdata<-as.data.frame(newdata)
  cx<-sapply(newdata,class)
  cx.factor<- cx == 'factor'
  object.cx.factor<- object$cx == 'factor'
  #if(!all(cx==object$cx))
  #  stop('number/order of columns of newdata does not match number/order of inputs used to train object')
  if(!all(cx.factor == object.cx.factor))
    stop('number/order of columns of newdata does not match number/order of inputs used to train object')


  newdata.des<-newdata[,!cx.factor,drop=F]
  newdata.cat<-newdata[,cx.factor,drop=F]

  if(ncol(newdata.des)>0){
    # for(i in 1:ncol(newdata.des)){
    #   newdata.des[,i]<-scale_range(newdata.des[,i],object$range.des[,i])
    # }
    # browser()
    newdata.des<-scale_range_mat(newdata.des,object$range.des)
  }
  tnewdata.des<-t(newdata.des)
  out<-array(dim=c(length(mcmc.use),nrow(newdata),nrow(newdata.func)))
  k<-0
  models<-object$model.lookup[mcmc.use]

  if(verbose)
    cat('Predict Start',myTimestamp(),'Models:',length(unique(models)),'\n')

  #func<-eval(parse(text=paste('mult',object$type,sep='')))
  func<-get(paste('mult',object$type,sep=''))

  mod.ind<-0
  for(j in unique(models)){ # loop though models, could be parallel?
    mod.ind<-mod.ind+1
    mcmc.use.j<-mcmc.use[models==j]
    ind<-k+(1:length(mcmc.use.j)) # index for storage
    k<-k+length(ind) # used for start of index
    out[ind,,]<-func(model=j,mcmc.use.mod=mcmc.use.j,object=object,tnewdata.des=tnewdata.des,newdata.cat=newdata.cat,tnewdata.func=tnewdata.func)

    if(verbose & mod.ind%%100==0)
      cat('Predict',myTimestamp(),'Model:',mod.ind,'\n')
  }

  if(nugget)
    return(drop(out)+rnorm(n=prod(dim(out)),sd=sqrt(object$s2[mcmc.use])))

  return(drop(out))
}



predict_fast<-function(object,newdata,newdata.func=NULL,mcmc.use=NULL,verbose=FALSE,...){
  if(is.null(mcmc.use)){ # if null, use all
    mcmc.use<-1:((object$nmcmc-object$nburn)/object$thin)
  }
  if(object$func){
    if(is.null(newdata.func))
      newdata.func<-object$xx.func
    else{
      dxf<-dim(newdata.func)
      if(is.null(dxf))
        newdata.func<-matrix(newdata.func)
      for(i in 1:ncol(newdata.func)){
        newdata.func[,i]<-scale_range(newdata.func[,i],object$range.func[,i])
      }
    }
  } else{
    newdata.func<-t(1) # placeholder
  }
  tnewdata.func<-t(newdata.func)

  dx<-dim(newdata)
  if(is.null(dx)){
    newdata<-data.frame(newdata)
    dx<-dim(newdata)
  }
  pd<-sum(object$pdes)+sum(object$pcat)
  if(dx[2]!=pd){
    newdata<-t(newdata)
    dx<-dim(newdata)
    if(dx[2]!=pd)
      stop('number of variables in newdata does not match number used in object')
  }
  newdata<-as.data.frame(newdata)
  cx<-sapply(newdata,class)
  cx.factor<- cx == 'factor'
  if(!all(cx==object$cx))
    stop('number/order of columns of newdata does not match number/order of inputs used to train object')

  newdata.des<-newdata[,!cx.factor,drop=F]
  newdata.cat<-newdata[,cx.factor,drop=F]

  #if(ncol(newdata.des)>0){
  #  for(i in 1:ncol(newdata.des)){
  #    newdata.des[,i]<-scale_range(newdata.des[,i],object$range.des[,i])
  #  }
  #}
  tnewdata.des<-t(newdata.des)
  out<-array(dim=c(length(mcmc.use),nrow(newdata),nrow(newdata.func)))
  k<-0
  models<-object$model.lookup[mcmc.use]

  if(verbose)
    cat('Predict Start',myTimestamp(),'Models:',length(unique(models)),'\n')

  func<-mult_des#eval(parse(text=paste('mult',object$type,sep='')))

  mod.ind<-0
  for(j in unique(models)){ # loop though models, could be parallel?
    mod.ind<-mod.ind+1
    mcmc.use.j<-mcmc.use[models==j]
    ind<-k+(1:length(mcmc.use.j)) # index for storage
    k<-k+length(ind) # used for start of index
    out[ind,,]<-func(model=j,mcmc.use.mod=mcmc.use.j,object=object,tnewdata.des=tnewdata.des,newdata.cat=newdata.cat,tnewdata.func=tnewdata.func)

    if(verbose & mod.ind%%100==0)
      cat('Predict',myTimestamp(),'Model:',mod.ind,'\n')
  }
  return(drop(out))
}



## make basis functions for model i - continuous portion
makeBasisMatrix<-function(i,nbasis,vars,signs,knots.ind,q,xxt,n.int,xx.train){
  n<-ncol(xxt)
  tbasis.mat<-matrix(nrow=nbasis+1,ncol=n)
  tbasis.mat[1,]<-1
  if(nbasis>0){
    for(m in 1:nbasis){
      if(n.int[i,m]==0){
        tbasis.mat[m+1,]<-1
      } else{
        use<-1:n.int[i,m]
        knots<-xx.train[cbind(knots.ind[i,m,use],vars[i,m,use])] # get knots from knots.ind
        tbasis.mat[m+1,]<-makeBasis(signs[i,m,use],vars[i,m,use],knots,xxt,q)
      }
    }
  }
  return(tbasis.mat)
}

# trying to speed up prediction for large number of basis functions: vectorizing like this doesn't help
# makeBasis_vec<-Vectorize(makeBasis,c('signs','vars','knots'))
#
# makeBasisMatrix<-function(i,nbasis,vars,signs,knots.ind,q,xxt,n.int,xx.train){
#   n<-ncol(xxt)
#   tbasis.mat<-matrix(nrow=nbasis+1,ncol=n)
#   tbasis.mat[1,]<-1
#   signs.list<-knots.list<-vars.list<-list()
#   if(nbasis>0){
#     for(m in 1:nbasis){
#       if(n.int[i,m]==0){
#         tbasis.mat[m+1,]<-1
#       } else{
#         use<-1:n.int[i,m]
#         knots.list[[m]]<-xx.train[cbind(knots.ind[i,m,use],vars[i,m,use])] # get knots from knots.ind
#         signs.list[[m]]<-signs[i,m,use]
#         vars.list[[m]]<-vars[i,m,use]
#         #tbasis.mat[m+1,]<-makeBasis(signs[i,m,use],vars[i,m,use],knots,xxt,q)
#       }
#     }
#     #browser()
#     tbasis.mat<-t(cbind(1,makeBasis_vec(signs.list,vars.list,knots.list,xxt,q)))
#   }
#   return(tbasis.mat)
# }

## make basis functions for model i - categorical portion
makeBasisMatrixCat<-function(i,nbasis,vars,xx,n.int,sub){
  n<-nrow(xx)
  tbasis.mat<-matrix(nrow=nbasis+1,ncol=n)
  tbasis.mat[1,]<-1
  for(m in 1:nbasis){
    if(n.int[i,m]==0){
      tbasis.mat[m+1,]<-1
    } else{
      use<-1:n.int[i,m]
      tbasis.mat[m+1,]<-makeBasisCat(vars[i,m,use],sub[[i]][[m]],xx)
    }
  }
  return(tbasis.mat)
}

## do multiplication to get yhat under the different scenarios
mult_des<-function(model,mcmc.use.mod,object,tnewdata.des,newdata.cat,tnewdata.func){
  M<-object$nbasis[mcmc.use.mod[1]]
  #browser()
  tmat.des<-makeBasisMatrix(model,M,object$vars,object$signs,object$knotInd,object$degree,tnewdata.des,object$n.int,object$xx.des)
  out<-object$beta[mcmc.use.mod,1:(M+1),drop=F]%*%tmat.des
  return(out)
}

mult_cat<-function(model,mcmc.use.mod,object,tnewdata.des,newdata.cat,tnewdata.func){
  M<-object$nbasis[mcmc.use.mod[1]]
  tmat.cat<-makeBasisMatrixCat(model,M,object$vars.cat,newdata.cat,object$n.int.cat,object$sub.list)
  out<-object$beta[mcmc.use.mod,1:(M+1),drop=F]%*%tmat.cat
  return(out)
}

mult_des_cat<-function(model,mcmc.use.mod,object,tnewdata.des,newdata.cat,tnewdata.func){
  M<-object$nbasis[mcmc.use.mod[1]]
  tmat.des<-makeBasisMatrix(model,M,object$vars.des,object$signs.des,object$knotInd.des,object$degree,tnewdata.des,object$n.int.des,object$xx.des)
  tmat.cat<-makeBasisMatrixCat(model,M,object$vars.cat,newdata.cat,object$n.int.cat,object$sub.list)
  out<-object$beta[mcmc.use.mod,1:(M+1),drop=F]%*%(tmat.des*tmat.cat)
  return(out)
}

mult_des_func<-function(model,mcmc.use.mod,object,tnewdata.des,newdata.cat,tnewdata.func){
  M<-object$nbasis[mcmc.use.mod[1]]
  tmat.des<-makeBasisMatrix(model,M,object$vars.des,object$signs.des,object$knotInd.des,object$degree,tnewdata.des,object$n.int.des,object$xx.des)
  tmat.func<-makeBasisMatrix(model,M,object$vars.func,object$signs.func,object$knotInd.func,object$degree,tnewdata.func,object$n.int.func,object$xx.func)
  out<-array(dim=c(length(mcmc.use.mod),ncol(tnewdata.des),ncol(tnewdata.func)))
  for(i in 1:length(mcmc.use.mod)){
    out[i,,]<-crossprod(diag(c(object$beta[mcmc.use.mod[i],1:(M+1)]),M+1)%*%tmat.des,tmat.func)
  }
  return(out)
}

mult_cat_func<-function(model,mcmc.use.mod,object,tnewdata.des,newdata.cat,tnewdata.func){
  M<-object$nbasis[mcmc.use.mod[1]]
  tmat.cat<-makeBasisMatrixCat(model,M,object$vars.cat,newdata.cat,object$n.int.cat,object$sub.list)
  tmat.func<-makeBasisMatrix(model,M,object$vars.func,object$signs.func,object$knotInd.func,object$degree,tnewdata.func,object$n.int.func,object$xx.func)
  out<-array(dim=c(length(mcmc.use.mod),nrow(newdata.cat),ncol(tnewdata.func)))
  for(i in 1:length(mcmc.use.mod)){
    out[i,,]<-crossprod(diag(c(object$beta[mcmc.use.mod[i],1:(M+1)]),M+1)%*%tmat.cat,tmat.func)
  }
  return(out)
}

mult_des_cat_func<-function(model,mcmc.use.mod,object,tnewdata.des,newdata.cat,tnewdata.func){
  M<-object$nbasis[mcmc.use.mod[1]]
  tmat.des<-makeBasisMatrix(model,M,object$vars.des,object$signs.des,object$knotInd.des,object$degree,tnewdata.des,object$n.int.des,object$xx.des)
  tmat.cat<-makeBasisMatrixCat(model,M,object$vars.cat,newdata.cat,object$n.int.cat,object$sub.list)
  tmat.func<-makeBasisMatrix(model,M,object$vars.func,object$signs.func,object$knotInd.func,object$degree,tnewdata.func,object$n.int.func,object$xx.func)
  out<-array(dim=c(length(mcmc.use.mod),ncol(tnewdata.des),ncol(tnewdata.func)))
  for(i in 1:length(mcmc.use.mod)){
    out[i,,]<-crossprod(diag(c(object$beta[mcmc.use.mod[i],1:(M+1)]),M+1)%*%(tmat.des*tmat.cat),tmat.func)
  }
  return(out)
}
