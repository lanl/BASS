\dontrun{
####################################################################################################
### univariate example
####################################################################################################
## simulate data (Friedman function)
f<-function(x){
  10*sin(pi*x[,1]*x[,2])+20*(x[,3]-.5)^2+10*x[,4]+5*x[,5]
}
sigma<-1 # noise sd
n<-500 # number of observations
x<-matrix(runif(n*10),n,10) #10 variables, only first 5 matter
y<-rnorm(n,f(x),sigma)

## fit BASS, no tempering
mod<-bass(x,y)
plot(mod)
## fit BASS, tempering
mod<-bass(x,y,temp.ladder=1.3^(0:8),start.temper=1000)
plot(mod)

## prediction
npred<-1000
xpred<-matrix(runif(npred*10),npred,10)
pred<-predict(mod,xpred,verbose=TRUE) # posterior predictive samples
true.y<-f(xpred)
plot(true.y,colMeans(pred),xlab='true values',ylab='posterior predictive means')
abline(a=0,b=1,col=2)

## sensitivity
sens<-sobol(mod)
plot(sens,cex.axis=.5)

####################################################################################################
### functional example
####################################################################################################
## simulate data (Friedman function with first variable as functional)
sigma<-1 # noise sd
n<-500 # number of observations
nfunc<-50 # size of functional variable grid
xfunc<-seq(0,1,length.out=nfunc) # functional grid
x<-matrix(runif(n*9),n,9) # 9 non-functional variables, only first 4 matter
X<-cbind(rep(xfunc,each=n),kronecker(rep(1,nfunc),x)) # to get y
y<-matrix(f(X),nrow=n)+rnorm(n*nfunc,0,sigma)

## fit BASS
mod<-bass(x,y,xx.func=xfunc)
plot(mod)

## prediction
npred<-100
xpred<-matrix(runif(npred*9),npred,9)
Xpred<-cbind(rep(xfunc,each=npred),kronecker(rep(1,nfunc),xpred))
ypred<-matrix(f(Xpred),nrow=npred)
pred<-predict(mod,xpred) # posterior predictive samples (each is a curve)
matplot(ypred,apply(pred,2:3,mean),type='l',xlab='observed',ylab='mean prediction')
abline(a=0,b=1,col=2)
matplot(t(ypred),type='l') # actual
matplot(t(apply(pred,2:3,mean)),type='l') # mean prediction

## sensitivity
sens<-sobol(mod,mcmc.use=1:10) # for speed, only use a few samples
plot(sens) # functional variable labelled "a"

sens.func<-sobol(mod,mcmc.use=1:10,func.var=1)
plot(sens.func)
}

## minimal example for CRAN testing
mod<-bass(1:2,1:2,nmcmc=2,nburn=1)