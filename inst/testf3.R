###############################################################
## f3: modified Friedman (1991) function (Francom et al., 2018), continuous variables with functional variable
###############################################################
library(BASS)

f<-function(x){
  10*sin(2*pi*x[,1]*x[,2])+20*(x[,3]-.5)^2+10*x[,4]+5*x[,5]
}
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


# true values
Si<-function(x,j=1:50){sum((-1)^(j-1)*(x^(2*j-1))/((2*j-1)*factorial(2*j-1)))}
Si2pi<-Si(2*pi);Si4pi<-Si(4*pi);a1<--5/pi*sum( (-4*pi^2)^(1:50)/((2*(1:50))*factorial(2*(1:50))) )
a2<-5/3;a3<-5;a4<-5/2;a5<-a2+a3+a4;a6<-a1
f0<-a1+a2+a3+a4
D0<-50-25*Si4pi/(2*pi)+5+125/3+2*(a1*a2+a1*a3+a1*a4+a2*a3+a2*a4+a3*a4) - f0^2
D1<-50/pi*(2*Si2pi-Si4pi)-a1^2
D2<-D1
D3<-5-a2^2
D4<-100/3-a3^2
D5<-25/3-a4^2
D12<-50-25*Si4pi/(2*pi)-100/pi*(2*Si2pi-Si4pi) + a1^2
c(D1,D2,D3,D4,D5,D12)/D0

xx<-xfunc
f0.func<-10*sin(pi*xx)^2/(pi*xx)+a2+a3+a4
D0.func<-50-25*sin(4*pi*xx)/(2*pi*xx)+5+125/3+2*(a2+a3+a4)*10*sin(pi*xx)^2/(pi*xx)+2*a2*a3+2*a2*a4+2*a3*a4-f0.func^2
D2.func<-50-25*sin(4*pi*xx)/(2*pi*xx)-100*sin(pi*xx)^4/(pi^2*xx^2)
D3.func<-rep(D3,nfunc)
D4.func<-rep(D4,nfunc)
D5.func<-rep(D5,nfunc)

Dmat<-cbind(D2.func/D0.func,D3.func/D0.func,D4.func/D0.func,D5.func/D0.func)
S.func<-t(apply(Dmat,1,cumsum))
matplot(S.func,type='l')


# save small model
bm<-bass(x,y,xx.func = xfunc,save.yhat = F,nmcmc=10000,nburn=9998,small=T)
ss<-BASS::sobol(bm)
ss.func<-BASS::sobol(bm,func.var = 1)
pred<-predict(bm,xpred)
saveRDS(bm,'tests/f3_mod.rda')
saveRDS(ss,'tests/f3_sob.rda')
saveRDS(ss.func,'tests/f3_sobFunc.rda')
saveRDS(xpred,'tests/f3_testX.rda')
saveRDS(pred,'tests/f3_testPred.rda')


