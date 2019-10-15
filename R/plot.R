###############################################################
## plot objects
###############################################################

#' @title BASS Plot Diagnostics
#'
#' @description Generate diagnostic plots for BASS model fit.
#' @param x a \code{bass} object.
#' @param quants quantiles for intervals, if desired.  NULL if not desired.
#' @param ... graphical parameters.
#' @details The first two plots are trace plots for diagnosing convergence.  The third plot is posterior predicted vs observed, with intervals for predictions.  The fourth plot is a histogram of the residuals (of the posterior mean model), with a red curve showing the assumed Normal density (using posterior mean variance).  If \code{bass} was run with \code{save.yhat = FALSE}, the third and fourth plots are omitted.
#' @export 
#' @import graphics
#' @seealso \link{bass}, \link{predict.bass}, \link{sobol}
#' @examples
#' # See examples in bass documentation.
#'
plot.bass<-function(x,quants=c(.025,.975),...){
  if(class(x)!='bass')
    stop('x must be an object of class bass')
  pred<-T
  if(is.null(x$yhat.mean))
    pred<-F
  op<-par(no.readonly=T)
  if(pred)
    par(mfrow=c(2,2))
  else
    par(mfrow=c(1,2))
  plot(x$nbasis,type='l',ylab='number of basis functions',xlab='MCMC iteration (post-burn)')
  plot(x$s2,type='l',ylab='error variance',xlab='MCMC iteration (post-burn)')
  if(pred){
    margin<-2
    if(x$func)
      margin<-2:3
    s<-sqrt(x$s2)
    if(!is.null(quants)){
      qq1<-apply(x$yhat+qnorm(quants[2])*s,margin,quantile,probs=quants[2])
      qq2<-apply(x$yhat+qnorm(quants[1])*s,margin,quantile,probs=quants[1])
      ylim=range(c(qq1,qq2))
      ylab='interval'
    } else{
      ylim=range(x$yhat.mean)
      ylab='mean'
    }
    plot(x$y,x$yhat.mean,ylim=ylim,ylab=paste('posterior predictive',ylab),xlab='observed',main='Training Fit',type='n',...)
    if(!is.null(quants))
      segments(x$y,qq1,x$y,qq2,col='lightgrey')
    points(x$y,x$yhat.mean)
    abline(a=0,b=1,col=2)
  
    hist(x$y-x$yhat.mean,freq=F,main='Posterior mean residuals',xlab='residuals')
    curve(dnorm(x,sd=mean(s)),col=2,add=T)
  }
  par(op)
}



#' @title BASS Plot Diagnostics
#'
#' @description Generate diagnostic plots for BASS model fit.
#' @param x a \code{bass} object.
#' @param quants quantiles for intervals, if desired.  NULL if not desired.
#' @param ... graphical parameters.
#' @details The first two plots are trace plots for diagnosing convergence.  The third plot is posterior predicted vs observed, with intervals for predictions.  The fourth plot is a histogram of the residuals (of the posterior mean model), with a red curve showing the assumed Normal density (using posterior mean variance).  If \code{bass} was run with \code{save.yhat = FALSE}, the third and fourth plots are omitted.
#' @export 
#' @import graphics
#' @seealso \link{bass}, \link{predict.bass}, \link{sobol}
#' @examples
#' # See examples in bass documentation.
#'
plot.bassOB<-function(x,quants=c(.025,.975),...){
  if(class(x)!='bassOB')
    stop('x must be an object of class bassOB')
  pred<-T
  if(is.null(x$yhat.mean))
    pred<-F
  op<-par(no.readonly=T)
  if(pred)
    par(mfrow=c(2,2))
  else
    par(mfrow=c(1,2))
  plot(x$nbasis,type='l',ylab='number of basis functions',xlab='MCMC iteration (post-burn)')
  plot(x$s2,type='l',ylab='error variance',xlab='MCMC iteration (post-burn)')
  if(pred){
    margin<-2
    if(x$func)
      margin<-2:3
    s<-sqrt(x$s2)
    if(!is.null(quants)){
      qq1<-apply(x$yhat+qnorm(quants[2])*s,margin,quantile,probs=quants[2])
      qq2<-apply(x$yhat+qnorm(quants[1])*s,margin,quantile,probs=quants[1])
      ylim=range(c(qq1,qq2))
      ylab='interval'
    } else{
      ylim=range(x$yhat.mean)
      ylab='mean'
    }
    plot(x$y,x$yhat.mean,ylim=ylim,ylab=paste('posterior predictive',ylab),xlab='observed',main='Training Fit',type='n',...)
    if(!is.null(quants))
      segments(x$y,qq1,x$y,qq2,col='lightgrey')
    points(x$y,x$yhat.mean)
    abline(a=0,b=1,col=2)
    
    hist(x$y-x$yhat.mean,freq=F,main='Posterior mean residuals',xlab='residuals')
    curve(dnorm(x,sd=mean(s)),col=2,add=T)
  }
  par(op)
}


#' @title Plot BASS sensitivity indices
#'
#' @description Generate plots for sensitivity analysis of BASS.
#' @param x a \code{bassSob} object, returned from \code{sobol}.
#' @param ... graphical parameters.
#' @details If \code{func.var} in the call to \code{sobol} was \code{NULL}, this returns boxplots of sensitivity indices and total sensitivity indices.  If there were functional variables, they are labeled with letters alphabetically.  Thus, if I fit a model with 4 categorical/continuous inputs and 2 functional inputs, the functional inputs are labeled a and b.  If \code{func.var} was not \code{NULL}, then posterior mean functional sensitivity indices are plotted, along with the functional partitioned variance.  Variables and interactions that are excluded did not explain any variance.
#' @export
#' @seealso \link{bass}, \link{predict.bass}, \link{sobol}
#' @examples
#' # See examples in bass documentation.
#'
plot.bassSob<-function(x,var.names=NULL,...){
  op<-par(no.readonly=T)
  par(mfrow=c(1,2),xpd=T)
  if(x$func){
    ord<-order(x$xx)
    x.mean<-apply(x$S,2:3,mean)
    matplot(x$xx[ord],t(apply(x.mean,2,cumsum))[ord,],type='l',xlab='x',ylab='proportion variance',ylim=c(0,1),main='Sensitivity',...)
    lab.x<-apply(x.mean,1,which.max)
    cs<-rbind(0,apply(x.mean,2,cumsum))
    cs.diff<-apply(x.mean,2,function(x) diff(cumsum(c(0,x))))
    text(x=x$xx[lab.x],y=cs[cbind(1:length(lab.x),lab.x)] + (cs.diff/2)[cbind(1:length(lab.x),lab.x)],x$names.ind,...)
    x.mean.var<-apply(x$S.var,2:3,mean)
    matplot(x$xx[ord],t(apply(x.mean.var,2,cumsum))[ord,],type='l',xlab='x',ylab='variance',main='Variance Decomposition',...)
  } else{
    boxplot(x$S,las=2,ylab='proportion variance',main='Sensitivity',range=0,...)
    boxplot(x$T,main='Total Sensitivity',range=0,...)
  }
  par(op)
}

