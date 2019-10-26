SBK_locp=function(y,X,initial,alpha0,xfixed,mhat,halpha,family,theta=0){
  mhatsbk=t(initial)
  n=ifelse(is.null(dim(X)),length(X),nrow(X))
  m=length(xfixed)
  XX=if (is.null(dim(X))){X} else {X[,alpha0]}

  xalpha=data.matrix(XX)

  # calculate distance matrix
  Sample=matrix(rep(xalpha,each=m),ncol=m,byrow=TRUE)
  fix=matrix(rep(xfixed,each=n),nrow=n)
  ualpha=Sample-fix
  # kernel: biweight
  kalpha=15/16*((1-(ualpha/halpha)^2+abs(1-(ualpha/halpha)^2))/2)^2/halpha
  mhat2=mhat[,-alpha0]
  thetaalpha=rowSums(as.matrix(mhat2))
  sapply(1:m,loc_fit,y=y,ualpha,kalpha,thetaalpha=thetaalpha,theta=theta,family=family)
}
