###############################################################
## f1: McKay 1997 function ("Nonparametric variance-based methods of assessing uncertainty importance"), legendre polynomials - continuous and categorical variables
###############################################################
library(BASS)

f<-function(x,t){
  #t<-t*4+1
  (t==1)*x + (t==2)*.5*(3*x^2-1) + (t==3)*.5*(5*x^3-3*x) + (t==4)*1/8*(35*x^4-30*x^2+3) + (t==5)*1/8*(63*x^5-70*x^3+15*x)
}


n<-500
x<-cbind(runif(n,-1,1),sample(1:5,size=n,replace=T),sample(1:5,size=n,replace=T),runif(n))
xf<-as.data.frame(x)
xf[,2]<-as.factor(xf[,2])
xf[,3]<-as.factor(xf[,3])
y<-f(x[,1],x[,2])
plot(x[,1],y)
bm<-bass(xf,y)
plot(bm)
var(predict(bm,xf,mcmc.use=1))

ss<-BASS::sobol(bm)
mean(ss$var.tot) # answer=3034/17325
plot(ss)

mean(ss$S$'1') # answer=.2
mean(ss$S$'1x2') # answer=.8

ntest<-100
xtest<-cbind(runif(ntest,-1,1),sample(1:5,size=ntest,replace=T),sample(1:5,size=ntest,replace=T),runif(ntest))
xftest<-as.data.frame(xtest)
xftest[,2]<-as.factor(xftest[,2])
xftest[,3]<-as.factor(xftest[,3])
plot(colMeans(predict(bm,xftest)),f(xtest[,1],xtest[,2])); abline(a=0,b=1,col=2)

# save small model
bm<-bass(xf,y,save.yhat = F,nmcmc=10000,nburn=9998,small=T)
ss<-BASS::sobol(bm)
pred<-predict(bm,xftest)
saveRDS(bm,'tests/f1_mod.rda')
saveRDS(ss,'tests/f1_sob.rda')
saveRDS(xftest,'tests/f1_testX.rda')
saveRDS(pred,'tests/f1_testPred.rda')


