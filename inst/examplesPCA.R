\dontrun{
  ## simulate data (Friedman function)
f<-function(x){
  10*sin(pi*x[,1]*x[,2])+20*(x[,3]-.5)^2+10*x[,4]+5*x[,5]
}
## simulate data (Friedman function with first variable as functional)
sigma<-.1 # noise sd
n<-500 # number of observations
nfunc<-50 # size of functional variable grid
xfunc<-seq(0,1,length.out=nfunc) # functional grid
x<-matrix(runif(n*9),n,9) # 9 non-functional variables, only first 4 matter
X<-cbind(rep(xfunc,each=n),kronecker(rep(1,nfunc),x)) # to get y
y<-matrix(f(X),nrow=n)+rnorm(n*nfunc,0,sigma)

## fit BASS
library(parallel)
mod<-bassPCA(x,y,n.pc=5,n.cores=min(5,parallel::detectCores()))
plot(mod$mod.list[[1]])
plot(mod$mod.list[[2]])
plot(mod$mod.list[[3]])
plot(mod$mod.list[[4]])
plot(mod$mod.list[[5]])

hist(mod$dat$trunc.error)

## prediction
npred<-100
xpred<-matrix(runif(npred*9),npred,9)
Xpred<-cbind(rep(xfunc,each=npred),kronecker(rep(1,nfunc),xpred))
ypred<-matrix(f(Xpred),nrow=npred)
pred<-predict(mod,xpred,mcmc.use=1:1000) # posterior predictive samples (each is a curve)
matplot(ypred,apply(pred,2:3,mean),type='l',xlab='observed',ylab='mean prediction')
abline(a=0,b=1,col=2)
matplot(t(ypred),type='l') # actual
matplot(t(apply(pred,2:3,mean)),type='l') # mean prediction

## sensitivity
sens<-sobolBasis(mod,int.order = 2,ncores = max(parallel::detectCores()-2,1),
                 mcmc.use=1000) # for speed, only use a few samples
plot(sens)


## calibration
x.true<-runif(9,0,1) # what we are trying to learn
yobs<-f(cbind(xfunc,kronecker(rep(1,nfunc),t(x.true)))) +
  rnorm(nfunc,0,.1) # calibration data (with measurement error)
plot(yobs)

cal<-calibrate.bassBasis(y=yobs,mod=mod,
                         discrep.mean=rep(0,nfunc),
                         discrep.mat=diag(nfunc)[,1:2]*.0000001,
                         sd.est=.1,
                         s2.df=50,
                         s2.ind=rep(1,nfunc),
                         meas.error.cor=diag(nfunc),
                         bounds=rbind(rep(0,9),rep(1,9)),
                         nmcmc=10000,
                         temperature.ladder=1.05^(0:30),type=1)

nburn<-5000
uu<-seq(nburn,10000,5)

pairs(rbind(cal$theta[uu,1,],x.true),col=c(rep(1,length(uu)),2),ylim=c(0,1),xlim=c(0,1))

pred<-apply(predict(mod,cal$theta[uu,1,],nugget = T,trunc.error = T,
      mcmc.use = cal$ii[uu]),3,function(x) diag(x)+rnorm(length(uu),0,sqrt(cal$s2[uu,1,1])))
qq<-apply(pred,2,quantile,probs=c(.025,.975))
matplot(t(qq),col='lightgrey',type='l')
lines(yobs,lwd=3)


}
## minimal example for CRAN testing
mod<-bassPCA(1:10,matrix(1:20,10),n.pc=2,nmcmc=2,nburn=1)
