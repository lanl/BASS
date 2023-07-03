###############################################################
## f5: Ishigami function, continuous variables (https://www.sfu.ca/~ssurjano/ishigami.html)
###############################################################
library(BASS)



f<-function(x,a=7,b=.1)
  sin(x[,1])*(1+b*x[,3]^4)+a*sin(x[,2])^2


b=.1
a=7


V1=(5+b*pi^4)^2/50
V2=a^2/8
V3=0
V13=8*b^2*pi^8/225
Vt=V1+V2+V3+V13
p1<-c(V1,V2,V3,V13)/Vt

V1=(1+3*b)^2*(exp(2)-1)/exp(2)/2
V2=a^2*(exp(4)-1)^2/exp(8)/8
V3=0
V13=48*b^2*(exp(2)-1)/exp(2)
Vt=V1+V2+V3+V13
p2<-c(V1,V2,V3,V13)/Vt

V1=(pi^2-8)*(5+b*pi^4)^2/(50*pi^2)
V2=a^2/8
V3=64*b^2*pi^6/225
V13=8*b^2*pi^6*(pi^2-8)/225
Vt=V1+V2+V3+V13
p3<-c(V1,V2,V3,V13)/Vt

n<-1000
p<-3
library(lhs)
xst<-maximinLHS(n,p)

low<- -4
high<- 4
rr<-c(low,high)

#f<-function(x)
#  x[,1] * x[,2] * x[,3]

x<-do.call(cbind,lapply(1:ncol(xst),function(i) BASS:::unscale.range(xst[,i],rr)))
#x<-rbind(x,as.matrix(expand.grid(rr,rr,rr)))
y<-f(x)

#library(rgl)
#plot3d(cbind(x[,c(1,2)],y))

#pairs(cbind(x,y))

library(BASS)
mod<-bass(x,y)
plot(mod)


prior1<-prior2<-prior3<-list()
for(i in 1:p){
  prior1[[i]]<-list(dist='uniform',trunc=c(-pi,pi))
  prior2[[i]]<-list(dist='normal',mean=0,sd=1,weights=1,trunc=mod$range.des[,i])
  prior3[[i]]<-list(dist='uniform',trunc=c(0,pi))
}

#BASS:::plot_prior(prior1[[1]])
#BASS:::plot_prior(prior2[[1]])
#curve(dnorm(x,0,1)*2.93,col=2,add=T)

BASS:::plot_prior(prior2[[1]])


ss1<-BASS::sobol(mod,prior1)
ss2<-BASS::sobol(mod,prior2)
ss3<-BASS::sobol(mod,prior3)

boxplot(ss1$S,range=0)
points(p1,col=2)
boxplot(ss2$S,range=0)
points(p2,col=2)
boxplot(ss3$S,range=0)
points(p3,col=2)






ntest<-100
xtest<-matrix(runif(ntest*p),ncol=p)
plot(colMeans(predict(mod,xtest)),f(xtest)); abline(a=0,b=1,col=2)

# save small model
bm<-bass(x,y,save.yhat = F,nmcmc=10000,nburn=9998,small=T)
ss1<-BASS::sobol(bm,prior1)
ss2<-BASS::sobol(bm,prior2)
ss3<-BASS::sobol(bm,prior3)
pred<-predict(bm,xtest)
saveRDS(bm,'tests/f5_mod.rda')
saveRDS(ss1,'tests/f5_sob1.rda')
saveRDS(ss2,'tests/f5_sob2.rda')
saveRDS(ss3,'tests/f5_sob3.rda')
saveRDS(xtest,'tests/f5_testX.rda')
saveRDS(pred,'tests/f5_testPred.rda')
saveRDS(prior1,'tests/f5_sobPrior1.rda')
saveRDS(prior2,'tests/f5_sobPrior2.rda')
saveRDS(prior3,'tests/f5_sobPrior3.rda')

