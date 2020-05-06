#######################################################
# Author: Devin Francom, Los Alamos National Laboratory
# Protected under GPL-3 license
# Los Alamos Computer Code release C19031
# github.com/lanl/BASS
#######################################################

########################################################################
## miscellaneous functions
########################################################################

## sample a tempered gamma
rgammaTemper<-function(n,shape,rate,itemper){
  rgamma(n,itemper*(shape-1)+1,itemper*rate)
}
## sample a tempered IG
rigammaTemper<-function(n,shape,scale,itemper){
  1/rgamma(n,itemper*(shape+1)-1,rate=itemper*scale)
}

## sample a truncated tempered IG
rtigammaTemper<-function(n,shape,scale,itemper,lower){
  1/rtgamma(n,1/lower,itemper*(shape+1)-1,rate=itemper*scale)
}

## sample from an upper-truncated gamma
rtgamma<-function(n,upper,shape,rate){
  out<-rep(upper,n)
  if(pgamma(upper,shape=shape,rate=rate)>0) # if cdf at upper bound is positive, sample, otherwise use upper bound
    out<-truncdist::rtrunc(n,'gamma',b=upper,shape=shape,rate=rate)
  return(out)
}

## scale a vector to be between 0 and 1
scale.range<-function(x,r=NULL){ # x is a vector
  if(is.null(r))
    r<-range(x)
  if((r[2]-r[1])==0)
    return(x-r[1])
  return((x-r[1])/(r[2]-r[1]))
}
## rescale a vector between 0 and 1 to range r
unscale.range<-function(x,r){
  x*(r[2]-r[1])+r[1]
}

## get yhat under the different scenarios
getYhat_des<-function(curr,nb){
  curr$des.basis%*%curr$beta
}
getYhat_cat<-function(curr,nb){
  curr$cat.basis%*%curr$beta
}
getYhat_des_cat<-function(curr,nb){
  curr$dc.basis%*%curr$beta
}
getYhat_des_func<-function(curr,nb){
  tcrossprod(curr$des.basis%*%diag(c(curr$beta),nb+1),curr$func.basis)
}
getYhat_cat_func<-function(curr,nb){
  tcrossprod(curr$cat.basis%*%diag(c(curr$beta),nb+1),curr$func.basis)
}
getYhat_des_cat_func<-function(curr,nb){
  tcrossprod(curr$dc.basis%*%diag(c(curr$beta),nb+1),curr$func.basis)
}

getYhat_des2<-function(des.basis,beta){
  des.basis%*%beta
}
getYhat_des_func2<-function(des.basis,func.basis,beta){
  tcrossprod(des.basis%*%diag(beta),func.basis)
}

## for checking inputs
posInt<-function(x){
  x==as.integer(x) & x>0
}

## replacement for timestamp(), since that seems to give Rstudio trouble on Windows
myTimestamp<-function(){
  x<-Sys.time()
  paste('#--',format(x,"%b %d %X"),'--#')
}
