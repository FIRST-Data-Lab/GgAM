library(mgcv)
library(BPST)
library(Triangulation)
library(ggplot2)


########## sydney  ############
# Triangulation #
library(HRW)
data('SydneyRealEstateBdry')
TriSydney=TriMesh(as.matrix(SydneyRealEstateBdry),n=15) #10,15,20
#B0=basis(as.matrix(TT$V),Tr=as.matrix(TT$Tr),d=2,r=1,Z=cbind(new.dataset$XCOORD,new.dataset$YCOORD))
setwd("C:\\Users\\JueWang\\Desktop\\Summer 2019\\GgAM Applications\\Application\\Sydney")
chernih <- read.csv("chernih.csv")
head(chernih)
new.dataset=chernih[,c('LogSalePrice','XCOORD','YCOORD','LotSize',
                       'Income','PM10','MainRoad','Coastline')]
new.dataset$LotSize=log(new.dataset$LotSize+0.01)
new.dataset$Coastline=log(new.dataset$Coastline+0.01)
# change
new.dataset$MainRoad=log(new.dataset$MainRoad+0.01)

### SMILE GGAM ###
# edf transform for SBL estimates
ed_trans=function(x){
  Fn=ecdf(x)
  Fn(x)
}
sydneydat=as.data.frame(new.dataset)
X_ed=apply(sydneydat[,c('LotSize','Income','MainRoad','Coastline')],2,ed_trans)

# sydneydat2=data.frame(sydneydat$LogSalePrice,X_ed,sydneydat$XCOORD, sydneydat$YCOORD)
# names(sydneydat2)<-c("LogSalePrice",'logLotSize','Income','logMainRoad','logCoastline','XCOORD',"YCOORD")
names(sydneydat)<-c("LogSalePrice",'XCOORD',"YCOORD",'logLotSize','Income', "PM10",'logMainRoad','logCoastline')
#SydneyHousing=sydneydat2
b_sydney=basis(V = TriSydney$V, Tr = TriSydney$Tr, d = 2, r = 1,Z=cbind(sydneydat$XCOORD,sydneydat$YCOORD))
saveRDS(b_sydney,'b_sydney.rds')

X_ed3=apply(sydneydat[b_sydney$Ind.inside,c('LotSize','Income','MainRoad','Coastline')],2,ed_trans)
sydneydat3=data.frame(sydneydat$LogSalePrice[b_sydney$Ind.inside],X_ed3,sydneydat$XCOORD[b_sydney$Ind.inside],
                      sydneydat$YCOORD[b_sydney$Ind.inside])
names(sydneydat3)<-c("LogSalePrice",'logLotSize','Income','logMainRoad','logCoastline','XCOORD',"YCOORD")

formula <- LogSalePrice ~ u(logLotSize) + u(Income) + u(logMainRoad) + u(logCoastline) + b(XCOORD,YCOORD, V = TriSydney$V, Tr = TriSydney$Tr, d = 2, r = 1)
formula2 <- LogSalePrice ~ u(logLotSize) + u(Income) + u(logMainRoad) + u(logCoastline) +
  b(XCOORD,YCOORD,V=TriSydney$V,Tr=TriSydney$Tr,d=2,r=1,B=b_sydney$B,Q2=b_sydney$Q2,K=b_sydney$K,ind=1:length(b_sydney$Ind.inside))
SmileSydney <- plbpsm(formula = formula2, data = as.data.frame(sydneydat[b_sydney$Ind.inside,]),ecdf=TRUE)
SmileSydney2 <- plbpsm(formula = formula, data = as.data.frame(sydneydat3))
sum((SmileSydney2$fitted.values.sbl-SmileSydney2$y)^2,
    na.rm = TRUE)/sum(!is.na(SmileSydney2$fitted.values.sbl))
SmileSydney$linear.names; SmileSydney$nonlinear.names
SmileSydney$coefficients
sum((SmileSydney$fitted.values.sbl-SmileSydney$y)^2,
    na.rm = TRUE)/sum(!is.na(SmileSydney$fitted.values.sbl))
### 0.08020877
resgam2=gam(formula=formula_gam_sydney2,data=as.data.frame(sydneydat3))
sum((resgam2$fitted.values-resgam2$y)^2,
    na.rm = TRUE)/sum(!is.na(resgam2$fitted.values))
### 0.1043919
resgam1=gam(formula=formula_gam_sydney1,data=as.data.frame(sydneydat3))
sum((resgam1$fitted.values-resgam1$y)^2,
    na.rm = TRUE)/sum(!is.na(resgam1$fitted.values))
### 0.1040206

# cv error #
set.seed(4321)
cv_sydey_5=cv.plbpsm(formula = formula, data = as.data.frame(sydneydat3),
                     kfold = 5, family = 'gaussian')
set.seed(4321)
cv_sydey_10=cv.plbpsm(formula = formula, data = as.data.frame(sydneydat3[b_sydney$Ind.inside,]),
                      kfold = 10, family = 'gaussian')
cv_sydey_5$cverror
### 0.1031939
cv_sydey_10$cverror
### 0.1029653
saveRDS(cv_sydey_5,'cv5_sydeny_0827.rds')
saveRDS(cv_sydey_10,'cv10_sydney_0827.rds')
saveRDS(SmileSydney2,'sydney_0827_v2.rds')
saveRDS(SmileSydney,'sydney_0827_v1.rds')
# gam cv #
cv5_sydeny_0827 <- readRDS("C:/Users/JueWang
                           /Desktop/Summer 2019/GgAM Applications/application_results_0827/cv5_sydeny_0827.rds")
flds=cv5_sydeny_0827$foldind
formula_gam_sydney1=LogSalePrice~s(logLotSize)+s(Income)+s(logMainRoad)+s(logCoastline)+
  s(XCOORD,YCOORD)
pred=rep(NA,dim(sydneydat3)[1])
for(kk in 1:length(flds)){
  dat.est=sydneydat3[-flds[[kk]],]
  resgam=gam(formula=formula_gam_sydney1,data=as.data.frame(dat.est))
  # a=predict(res_gauss2,as.data.frame(dat_gauss_GGAMS[c(((f-1)*100+1):(f*10)),]))
  pred[flds[[kk]]]=predict(resgam,as.data.frame(sydneydat3[flds[[kk]],]),type='response')
}
cverror_gam_sydney1=sum((pred-sydneydat3[,as.character(formula_gam_sydney1[[2]])] )^2,na.rm = TRUE)/sum(!is.na(pred))
### 0.1090294 ###
formula_gam_sydney2=LogSalePrice~s(logLotSize)+s(Income)+logMainRoad+s(logCoastline)+
  s(XCOORD,YCOORD)
pred=rep(NA,dim(sydneydat3)[1])
for(kk in 1:length(flds)){
  dat.est=sydneydat3[-flds[[kk]],]
  resgam=gam(formula=formula_gam_sydney2,data=as.data.frame(dat.est))
  # a=predict(res_gauss2,as.data.frame(dat_gauss_GGAMS[c(((f-1)*100+1):(f*10)),]))
  pred[flds[[kk]]]=predict(resgam,as.data.frame(sydneydat3[flds[[kk]],]),type='response')
}
cverror_gam_sydney2=sum((pred-sydneydat3[,as.character(formula_gam_sydney2[[2]])] )^2,na.rm = TRUE)/sum(!is.na(pred))
### 0.1047512

# summary & plots #
summary(SmileSydney)
SummarySydney <- summary(SydneySmile)
Summary$bands[[1]]$est
Summary$bands[[1]]$upper
Summary$bands[[1]]$lower
plot(SmileSydney, ab.line=TRUE,add.rug=TRUE)

res.plot <- plot(SmileSydney, n1 = 450,n2 = 450)
res.plot$basisplot$p.beta+geom_polygon(data=SydneyRealEstateBdry,aes(x=longitude,y=latitude),inherit.aes=F, colour='black', fill=NA, lwd=1)
saveRDS(res.plot,'plot_sydney_0827.rds')
