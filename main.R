######################################################################
# Project: Estimation and Inference for Generalized Geoadditive Models
# Estimation for negative binomial in the HorseShoe domain
# Author: Shan Yu, Li Wang and GuanNan Wang
# Date: 08/18/2019
######################################################################
rm(list = ls())
# library
library(mgcv)
library(splines)
library(bindata)
library(parallel)
library(BPST)
# install.packages('~/Dropbox/GGAM/Code/GgAM_0906', repos = NULL, type="source")
library(GgAM)
# source('../Data_Generator_poi_0603.R')

# true beta at gridded points
X0=matrix(rep(seq(-0.5,0.5,0.01),each=2),ncol=2,byrow=TRUE)
m1=F2(X0[,1])
m2=F3(X0[,2])
# beta0.3=beta.3(X0[,3])
fsb=list(fs.boundary())[[1]]
xx=seq(-1,4,0.05)
yy=seq(-1,1,0.05)
S0=expand.grid(xx,yy)
v=S0[,1]
w=S0[,2]
m=fs.test(v,w,b=1)
names(fsb)=c("v","w")
ind=inSide(fsb,x=v,y=w)
S0=S0[ind,]

m0=m[ind]
m0=m0-mean(m0)
m0=m0/2
# triangulation information
library(Triangulation)
Tri.hs=TriMesh(horseshoe,n=4)
V=Tri.hs$V
Tr=Tri.hs$Tr

n=1000L # 2000 # sample size
p.x=3L # dimension for x
s.x=3L # dimension for x with signals
mu0=0 # intercept: 1 for "binomial, [-1.8,-.9] for "poisson"
penalty="grALasso" # "grALasso" "grSCAD"
family='gaussian' # poisson # binomial # negbin(c(2,8)) #
sig=0.5
Simulation<-function(iter){
  iter=1
  set.seed(2*iter)
  print(iter)
  data=Data_Generator(F1,F2,F3,n,p.x,s.x,mu0,sig,scale=1,family)
  y.sam=data$y
  x.sam=data$x
  s.sam=data$s
  dat_poi=as.data.frame(cbind(y.sam,x.sam,s.sam))
  colnames(dat_poi)<-c('y','x1','x2','x3','s1','s2')
  formula=y~x1+u(x2,N=3)+u(x3,N=3)+b(s1,s2,V=V,Tr=Tr,d=2,r=1)
  fit=plbpsm(formula=formula,data=as.data.frame(dat_poi),family = family,
             group=c(0,0))
  fit$se_beta
  plot(fit)
  # family = negbin(fit$est_theta)
  # fit=ggam(Y=Data$Y,X=Data$X,S=Data$S,Tr=Tr,V=V,N=4,backfitting=TRUE,
  #          lambda=lambda,off=0,family=nb_bps(),nb=TRUE)
  # pred=ggam.pred(fit,X0,S0,V,Tr,testlevel=0.05)
  # pred=predict(fit,as.data.frame(pop),type = 'terms')
  # m1_hat=pred[,"u(x1)"]
  # m2_hat=pred[,"u(x2)"]
  # m_hat=pred[,"b(s1,s2)"]
  B0=basis(as.matrix(V),as.matrix(Tr),2,1,as.matrix(S0))
  newB=B0$B
  newind=B0$Ind.inside
  m_hat=as.matrix(newB)%*%fit$coefficients_bivariate
  beta_est1=fit$coefficients[1]
  ind=seq(10,90,4)
  bands=summary(fit,X0=X0,Ind=ind)
  b1=bands$bands[[1]]
  b2=bands$bands[[2]]
  m1_hat=bands$bands[[1]]$m1_sbk
  m2_hat=bands$bands[[2]]$m1_sbk
  
  aa=c((1-beta_est1)^2,
       sum((m1-m1_hat)^2)*0.01,
       sum((m2-m2_hat)^2)*0.01,
       sum((m0[newind]-m_hat)^2)*0.05^2)
  
  # coverage
  
  cover=c(prod((m1[ind]>=b1$lower& m1[ind]<= b1$upper))==1,
          prod((m2[ind]>=b2$lower& m2[ind]<= b2$upper))==1)
  all=c(aa,cover,fit$coefficients[1],fit$se_beta)
  # results
  names(all)=c('beta1.MSE','alpha1.MISE',"alpha2.MISE",'b.MISE','alpha1.cover','alpha2.cover',"est.beta1",
               'std.err.beta1')
  print(all)
  # SCBs plots
  # plot(b1$x0,b1$est,type='l',ylim=c(min(b1$lower),max(b1$upper)),xlab='',main='u(x2)',ylab='')
  # polygon(c(b1$x0, rev(b1$x0)), c(b1$upper, rev(b1$lower)),
  #         col = "grey", border = NA)
  # lines(b1$x0,b1$est,type='l')
  # lines(b1$x0,F2(b1$x0),lty=2)
  # plot(b2$x0,b2$est,type='l',ylim=c(min(b2$lower),max(b2$upper)),xlab='',main='u(x3)',ylab='')
  # polygon(c(b2$x0, rev(b2$x0)), c(b2$upper, rev(b2$lower)),
  #         col = "grey", border = NA)
  # lines(b2$x0,b2$est,type='l')
  # lines(b2$x0,F3(b2$x0),lty=2)
  
  return(all)
}
# Simulation(2)
family='gaussian'
result1000_gaussian=mclapply(1:100,Simulation,mc.cores=1)
result1000_gaussian2=mclapply(1:500,Simulation,mc.cores=1)
result2=do.call(rbind,result1000_gaussian2)
colMeans(result2) #
result=do.call(rbind,result1000_gaussian)
colMeans(result2) #
sd(result2[,7])
# 0.05788071
family='poisson'
result1000_poi=mclapply(1:500,Simulation,mc.cores=1)
n=2000
result2000_gaussian2=mclapply(1:500,Simulation,mc.cores=1)

# record standard error & others
result=do.call(rbind,result1000_gaussian2)
colMeans(result)
res_1000_poi=matrix(colMeans(result),nr=1)

saveRDS(result,'scb_gaus_1000_0912.rds')
p1=matrix(c(sd(result[,7]),
            mean(result[,8]),
            median(result[,8]),
            IQR(result[,8])/1.349),nc=1)
p1
result=do.call(rbind,result2000_gaussian2)
colMeans(result)

res_2000_poi=matrix(colMeans(result),nr=1)

saveRDS(result,'scb_gaus_2000_0912.rds')
p2=matrix(c(sd(result[,7]),
            mean(result[,8]),
            median(result[,8]),
            IQR(result[,8])/1.349),nc=1)
result_gaus=cbind(p1,p2)
saveRDS(result_gaus,'se_gaus_0912.rds')

# cv.error #
n=1000 # 2000
require(caret)
pred_all=c()
cverror_all=c()
for (i in 1:500){
  set.seed(i)
  print(i)
  data=Data_Generator(F1,F2,F3,n,p.x,s.x,mu0,sig,scale=1,family)
  y.sam=data$y
  x.sam=data$x
  s.sam=data$s
  dat_poi_GGAMS=as.data.frame(cbind(y.sam,x.sam,s.sam))
  colnames(dat_poi_GGAMS)<-c('y','x1','x2','x3','s1','s2')
  formula=y~x1+u(x2)+u(x3)+b(s1,s2,V=V,Tr=Tr,d=2,r=1)
  flds <- createFolds(1:n, k = 10, list = TRUE, returnTrain = FALSE)
  pred=rep(NA,n)
  for(kk in 1:length(flds)){
    dat.est=dat_poi_GGAMS[-flds[[kk]],]
    res_poi2=plbpsm(formula=formula,data=as.data.frame(dat.est),
                    family = family,	group=c(0,0))
    # a=predict(res_poi2,as.data.frame(dat_poi_GGAMS[c(((f-1)*100+1):(f*10)),]))
    pred[flds[[kk]]]=predict(res_poi2, as.data.frame(dat_poi_GGAMS[flds[[kk]],]))
  }
  cverror=sum((pred-dat_poi_GGAMS$y)^2,na.rm = TRUE)/sum(!is.na(pred))
  cverror_all=c(cverror_all,cverror)
  pred_all=cbind(pred_all,pred)
  print(cverror)
}
mean(cverror_all)
