#' @title Write/read BASS object to JSON
#'
#' @description Write/read an object of class \code{bass}, \code{bassBasis}, or \code{bassSob} to a JSON file.
#' @param x an object of class \code{bass}, \code{bassBasis}, or \code{bassSob}, returned from \code{bass}, \code{bassPCA}, \code{bassBasis}, \code{sobol}, or \code{sobolBasis}.
#' @param file json file to write or read.
#' @param ... further arguments passed to or from other methods.
#' @name writeBassJSON
#' @export
NULL
#' @rdname writeBassJSON
writeBassJSON<-function(x,file,...){
  cl<-class(x)

  if(cl=='bass'){
    class(x)<-'list'
    class(x$call)<-'list'
    for(i in 1:length(x$call))
      class(x$call[[i]])<-'character'
    jj<-jsonlite::toJSON(x,digits=NA,...)
    jsonlite::write_json(jj,file)
  }

  if(cl=='bassSob'){
    class(x)<-'list'
    for(i in 1:length(x$prior))
      class(x$prior[[i]])<-'list'
    jj<-jsonlite::toJSON(x,digits=NA,...)
    jsonlite::write_json(jj,file)
  }

  if(cl=='bassBasis'){
    class(x)<-'list'
    class(x$dat)<-'list'
    names(x$mod.list)<-paste0('n',1:length(x$mod.list)) # may need to do something like this for tempering, as well
    for(j in 1:length(x$mod.list)){
      class(x$mod.list[[j]])<-'list'
      class(x$mod.list[[j]]$call)<-'list'
      for(i in 1:length(x$mod.list[[j]]$call))
        class(x$mod.list[[j]]$call[[i]])<-'character'
    }
    jj<-jsonlite::toJSON(x,digits=NA,...)
    jsonlite::write_json(jj,file)
  }

  #return(NULL)
}
#' @rdname writeBassJSON
#' @export
readBassJSON<-function(file){
  jj<-jsonlite::fromJSON(file)#,simplifyDataFrame=F,simplifyVector=T,simplifyMatrix=F)
  x<-jsonlite::fromJSON(jj)#,simplifyDataFrame=F,simplifyVector=T,simplifyMatrix=F)

  if(!is.null(x$beta)){
    class(x)<-'bass'
    for(i in 1:length(x$call))
      x$call[[i]]<-as.name(x$call[[i]])
    x$call<-as.call(x$call)
  }

  if(!is.null(x$S)){
    class(x)<-'bassSob'
    for(i in 1:length(x$prior))
      class(x$prior[[i]])<-x$prior[[i]]$dist
  }

  if(!is.null(x$mod.list)){
    class(x)<-'bassBasis'
    for(j in 1:length(x$mod.list)){
      class(x$mod.list[[j]])<-'bass'
      #for(i in 1:length(x$mod.list[[j]]$call))
      #  x$mod.list[[j]]$call[[i]]<-as.name(x$mod.list[[j]]$call[[i]])
      #x$mod.list[[j]]$call<-as.call(x$mod.list[[j]]$call)
    }
  }

  return(x)
}

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

#' @title Print BASS Details
#'
#' @description Print some of the details of a BASS model.
#' @param x a \code{bassBasis} object, returned from \code{bassPCA} or \code{bassBasis}.
#' @param ... further arguments passed to or from other methods.
#' @export
#'
print.bassBasis<-function(x,...){
  # numeric inputs
  p.cat<-sum(x$mod.list[[1]]$cx=='factor')
  p.num<-x$mod.list[[1]]$p-p.cat
  reps<-nrow(x$mod.list[[1]]$xx.des)

  if(p.cat==0)
    cat("\n              Number of variables: ",x$mod.list[[1]]$p,sep='')
  if(p.cat>0)
    cat("\n              Number of variables: ",x$mod.list[[1]]$p,' (',p.cat,' categorical)',sep='')
  cat("\n                      Sample size: ",reps,sep='')
  cat("\n   Number of output basis functions: ",length(x$mod.list))

  cat('\n\n')

}


#' @title Summarize BASS Details
#'
#' @description Summarize some of the details of a BASS model.
#' @param object a \code{bassBasis} object, returned from \code{bassPCA} or \code{bassBasis}.
#' @param ... further arguments passed to or from other methods.
#' @export
#'
summary.bassBasis<-function(object,...){
  x<-object
  p.cat<-sum(x$mod.list[[1]]$cx=='factor')
  p.num<-x$mod.list[[1]]$p-p.cat
  reps<-nrow(x$mod.list[[1]]$xx.des)

  if(p.cat==0)
    cat("\n              Number of variables: ",x$mod.list[[1]]$p,sep='')
  if(p.cat>0)
    cat("\n              Number of variables: ",x$mod.list[[1]]$p,' (',p.cat,' categorical)',sep='')
  cat("\n                      Sample size: ",reps,sep='')
  cat("\n   Number of output basis functions: ",length(x$mod.list))

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


