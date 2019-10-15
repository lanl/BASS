#' @title Print BASS Details
#'
#' @description Print some of the details of a BASS model. 
#' @param x a \code{bass} object, returned from \code{bass}.
#' @param ... further arguments passed to or from other methods.
#' @export
#'
print.bass<-function(x,...){
  cat("\nCall:\n",deparse(x$call),'\n')
  
  
  # numeric inputs
  p.cat<-sum(x$cx=='factor')
  p.num<-x$p-p.cat
  p.func<-ncol(x$xx.func)
  reps<-nrow(x$xx.des)
  func.grid<-nrow(x$xx.func)
  
  if(p.cat==0)
    cat("\n              Number of variables: ",x$p,sep='')
  if(p.cat>0)
    cat("\n              Number of variables: ",x$p,' (',p.cat,' categorical)',sep='')
  cat("\n                      Sample size: ",reps,sep='')
  if(!is.null(p.func)){
    cat("\n   Number of functional variables:",p.func)
    cat("\n             Functional grid size:",func.grid)
  }
  
  cat('\n\n')
  
}


#' @title Summarize BASS Details
#'
#' @description Summarize some of the details of a BASS model. 
#' @param object a \code{bass} object, returned from \code{bass}.
#' @param ... further arguments passed to or from other methods.
#' @export
#'
summary.bass<-function(object,...){
  cat("\nCall:\n",deparse(object$call),'\n')
  
  # numeric inputs
  p.cat<-sum(object$cx=='factor')
  p.num<-object$p-p.cat
  p.func<-ncol(object$xx.func)
  reps<-nrow(object$xx.des)
  func.grid<-nrow(object$xx.func)
  
  if(p.cat==0)
    cat("\n              Number of variables: ",object$p,sep='')
  if(p.cat>0)
    cat("\n              Number of variables: ",object$p,' (',p.cat,' categorical)',sep='')
  cat("\n                      Sample size: ",reps,sep='')
  if(!is.null(p.func)){
    cat("\n   Number of functional variables:",p.func)
    cat("\n             Functional grid size:",func.grid)
  }
  
  cat("\n\nNumber of basis functions (range):",range(object$nbasis))
  cat("\n    Posterior mean error variance:",mean(object$s2),'\n\n')
  
}


