predict_backfitting<- function(object,newdata,type,ColNames,...) {

  term.smooth=sapply(object$basis_info, class)[1,]
  Terms <- delete.response(object$pterms)
  Xpp <- model.matrix(Terms,newdata)

  smooth.bivariate=which(term.smooth=="bivariate.smooth")
  if(length(smooth.bivariate)>0){
    basisinfo=object$basis_info[[smooth.bivariate]]
    loc=newdata[,basisinfo$term]
    B0=basis(as.matrix(basisinfo$V),as.matrix(basisinfo$Tr),basisinfo$d,basisinfo$r,as.matrix(loc))
    newB=B0$B
    newind=B0$Ind.inside
    ## prediction of the nonlinear part
    newbeta=newB%*%object$coefficients_bivariate
    mm.test=newbeta
    Xp=Xpp[newind,,drop=FALSE]
    yy=object$y[basisinfo$ind]
  } else {
    newbeta=0
    mm.test=0
    Xp=Xpp
    yy=object$y
    newind=1:dim(newdata)[1]
  }
print(dim(Xp))
print(length(Xp))

  if (object$MI){
    if (length(object$ind.l)>0){
      Xp=cbind(Xp,newdata[newind, sapply(object$basis_info_MI[object$ind.l],function(X){X$term})])
    }
  }

  if (!is.null(dim(Xp))) {
    if (length(smooth.bivariate)>0 && object$intercept){
      Xp2=Xp[,-1,drop=FALSE]
    } else {
      Xp2=Xp
    }
    print(dim(Xp2)[2])
    if (dim(Xp2)[2]>0){
      mhat.test=Xp2*object$coefficients[1:(dim(Xp2)[2])]
    } else {
      mhat.test=NULL
    }
    first.para=alpha0=dim(as.matrix(Xp2))[2]
  } else {
    if (length(smooth.bivariate)>0 && object$intercept){
      Xp2=NULL
      mhat.test=NULL
      first.para=0
      } else {
      Xp2=as.matrix(Xp)
      mhat.test=Xp2*object$coefficients[1:(dim(Xp2)[2])]
      first.para=alpha0=dim(as.matrix(Xp2))[2]
    }
  }

  # Xp2=as.matrix(Xp[,-1])
  # mhat.test=Xp2%*%object$coefficients[1:(dim(Xp2)[2])]
  if (object$MI){
    smooth.univariate=object$ind.nl
  } else {
    smooth.univariate=which(term.smooth=="univariate.smooth")
  }
  #first.para=alpha0=dim(as.matrix(Xp2))[2]
  pred<-NULL
  if (!is.null(mhat.test)){
    temp<-matrix(NA,ncol=dim(mhat.test)[2],nrow=dim(newdata)[1])

  }
  if (!is.null(dim(mhat.test)[2])){
    temp[newind,]<-as.matrix(mhat.test)
  } else if (length(mhat.test)>0){
    temp[newind]<-mhat.test} else #if (is.null(mhat.test))
      {
    temp<-NULL}
  pred<-cbind(pred,temp)

  if (length(smooth.univariate)>0){
  for(k in 1:length(smooth.univariate)){
    x0=newdata[newind,object$basis_info[[smooth.univariate[k]]]$term]

    #### back-fitting estimation ####
    alpha0=first.para+k
    initial=runif(length(x0))-0.5
    m_sbk=SBK_locp(yy,object$X2,initial,alpha0,x0,mhat=as.matrix(object$mhat),object$h_opt_all[k],family=object$family)
    mhat.test=cbind(mhat.test,m_sbk)
    if (type=='terms'){
      temp<-rep(NA,dim(newdata)[1])
      temp[newind]<-m_sbk
      pred<-cbind(pred,temp)
    }
  }
  }
  print(dim(mhat.test))
  if (type=='terms'){
    temp2<-rep(NA,dim(newdata)[1])
    temp2[newind]<-mm.test
    pred<-cbind(pred,temp2)
    colnames(pred)<-ColNames

  }

  mhat.test=cbind(mhat.test,mm.test)
  if (type=='response'){
     eta_test=rowSums(as.matrix(mhat.test))
  y_pred=object$family$linkinv(eta_test)
  pred=rep(NA,dim(newdata)[1])
  pred[newind]<-y_pred
  }
  return(list(pred=pred,newind=newind))
}
