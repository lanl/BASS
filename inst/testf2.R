###############################################################
## f2: modified Friedman (1991) function (Francom et al., 2018), continuous variables
###############################################################
library(BASS)

f<-function(x){
  10*sin(2*pi*x[,1]*x[,2])+20*(x[,3]-.5)^2+10*x[,4]+5*x[,5]
}
sigma<-1 # noise sd
n<-500 # number of observations
x<-matrix(runif(n*10),n,10) #10 variables, only first 5 matter
y<-rnorm(n,f(x),sigma)

## fit BASS, no tempering
mod<-bass(x,y)
plot(mod)

## prediction
npred<-100
xpred<-matrix(runif(npred*10),npred,10)
pred<-predict(mod,xpred,verbose=TRUE) # posterior predictive samples
true.y<-f(xpred)
plot(true.y,colMeans(pred),xlab='true values',ylab='posterior predictive means')
abline(a=0,b=1,col=2)

## sensitivity
sens<-sobol(mod)
plot(sens,cex.axis=.5)

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


# save small model
bm<-bass(x,y,save.yhat = F,nmcmc=10000,nburn=9998,small=T)
ss<-BASS::sobol(bm)
pred<-predict(bm,xpred)
saveRDS(bm,'tests/f2_mod.rda')
saveRDS(ss,'tests/f2_sob.rda')
saveRDS(xpred,'tests/f2_testX.rda')
saveRDS(pred,'tests/f2_testPred.rda')


