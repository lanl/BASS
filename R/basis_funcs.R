#######################################################
# Author: Devin Francom, Los Alamos National Laboratory
# Protected under GPL-3 license
# Los Alamos Computer Code release C19031
# github.com/lanl/BASS
#######################################################

########################################################################
## make basis functions
########################################################################

## makes negative values 0
pos<-function(vec){
  #replace(vec,vec<0,0)
  (abs(vec)+vec)/2
}

## largest value of basis function, assuming x's in [0,1], used for scaling
const<-function(signs,knots,degree){
  cc<-prod(((signs+1)/2 - signs*knots))^degree
  if(cc==0)
    return(1)
  return(cc)
} # since a product, can find for functional & categorical pieces separately, take product

## make basis function (from continuous variables)
makeBasis<-function(signs,vars,knots,datat,degree){
  cc<-const(signs,knots,degree)
  temp1<-pos(signs*(datat[vars,,drop=F]-knots))^degree # this only works for t(data)...
  if(length(vars)==1){
    return(c(temp1)/cc)
  } else{
    temp2<-1
    for(pp in 1:length(vars)){ # faster than apply
      temp2<-temp2*temp1[pp,]
    }
    return(temp2/cc)
  }
}

## make basis function (from categorical variables)
makeBasisCat<-function(vars,sub,data){
  temp<-1
  for(ii in 1:length(vars)){
    temp<-temp*as.numeric(data[,vars[ii]] %in% sub[[ii]])
  }
  return(temp)
}
