#' @importFrom mgcv gam
#' @importFrom stats approxfun density predict
key_functionsSBLL=function(y,X,alpha,mhat,sigma_2,xfixed,family,basisinfo,theta=0){

  #### f_k density ####
  #df=approxfun(density(X[,alpha],kernel='biweight'))
  df=if(is.null(dim(X))){approxfun(density(X,kernel='biweight'))} else {
    approxfun(density(X[,alpha],kernel='biweight'))}

  df_fixed=df(xfixed)
  linkinv=family$linkinv
  mu.eta=family$mu.eta
  if(theta==0) {variance=family$variance} else {
    variance=function(mu,theta) mu+mu^2/theta
  }
  if (family$link=='identity'){
    gprime=function(mu) 1
  }  else if (family$link=='logit'){
    gprime=function(mu) 1/(mu-mu^2)
  } else if (family$link=='log'){
    gprime=function(mu) 1/mu
  } else if (family$link=='inverse'){
    gprime=function(mu) -1/mu^2

  }


  #### E \rho_2 ####
  #muhat=linkinv(rowSums(mhat))
  mhat=as.matrix(mhat)
  muhat=if(dim(mhat)[2]>2){linkinv(rowSums(mhat))} else {linkinv(mhat)}

  if(theta==0) erho2=1/(sigma_2*gprime(muhat)^2*variance(muhat)) else
    erho2=1/(sigma_2*gprime(muhat)^2*variance(muhat,theta))
  # xalpha=data.matrix(X[,alpha])
  # n=nrow(X)
  # m=length(xfixed)
  # # calculate distance matrix
  # Sample=matrix(rep(xalpha,each=m),ncol=m,byrow=TRUE)
  # fix=matrix(rep(xfixed,each=n),nrow=n)
  # ualpha=Sample-fix
  # # kernel: biweight
  # halpha=0.2
  # kalpha=15/16*((1-(ualpha/halpha)^2+abs(1-(ualpha/halpha)^2))/2)^2/halpha
  # esigma=t(erho2)%*%kalpha/colSums(kalpha)
  # esigma=as.vector(esigma)
  XX=if (is.null(dim(X))){X} else {X[,alpha]}

  dd=data.frame(erho2,XX)
  names(dd)=c('V1','V2')
  fit=gam(V1~s(V2),data=dd)
  dd0=data.frame(xfixed)
  names(dd0)='V2'
  esigma=predict(fit,newdata=dd0)
  #### \beta \prime \prime ####
  dd1=data.frame(y,XX)
  names(dd1)=c('V1','V2')
  #fit=gam(V1~s(V2),data=dd1,family=quasipoisson(),offset=rowSums(mhat[,-alpha]))

  # eps=1e-7
  # X0=predict(fit,dd0,type='lpmatrix')
  # newDFeps_p=dd0+eps
  # X1=predict(fit,newDFeps_p,type='lpmatrix')
  # newDFeps_m=dd0-eps
  # X_1=predict(fit,newDFeps_m,type='lpmatrix')
  # # design matrix for second derivative
  # Xpp=(X1+X_1-2*X0)/eps^2
  # # second derivative
  # beta_dev2=as.vector(Xpp%*%coef(fit))
  # # # plot(xfixed,fd_d2)
  # # #

  ### CHANGE! prior knots###
  knots=Basis_generator(XX,basisinfo$N,basisinfo$q,basisinfo$KnotsLocation,knots=basisinfo
                        $prior.knots)$knots
  if (!is.null(dim(X))){
    ZZ=truncated_spline(X[,alpha],knots)
  } else {
    ZZ=truncated_spline(matrix(X,ncol=1),knots)
  }
  Z=ZZ$XX
  thetaalpha=if(dim(mhat)[2]>2){rowSums(mhat[,-alpha])} else {mhat[,-alpha]}

  if (family$family[1]=='gaussian'){ # problems!
    W=t(Z)%*%Z
    W2=t(Z)%*%(y-thetaalpha)
    b_new=ginv(W)%*%W2
    # need to be changed
  } else{
    step=0
    b_new=rep(0,dim(Z)[2])
    b=rep(1,dim(Z)[2])
    while(sum((b_new-b)^2)>0.00001 & step <= 10){
      step=step+1
      b=b_new
      eta=Z%*%b_new+thetaalpha
      mu=linkinv(eta)
      mevg=mu.eta(eta)
      if(theta==0) {var=variance(mu)} else {
        variance=function(mu,theta) mu+mu^2/theta
        var=variance(mu,theta)
      }
      Y_iter=(eta-thetaalpha)+(y-mu)/mu.eta(eta)
      W_iter=as.vector((mevg^2)/var)
      temp1=W_iter*Z
      temp2=W_iter*Y_iter
      temp3=crossprod(Z,temp1)
      b_new=solve(temp3,crossprod(Z,temp2))
      if(step==0) print('not convergent')
    }
  }

  XX2=truncated_spline(xfixed,knots)$XX2
  beta_dev2=as.vector(XX2%*%b_new)
  # #
  # #### variance function ####
  # # v2=(y-variance(muhat))^2
  # # #esigma=XX%*%solve(t(XX)%*%XX,t(XX)%*%sigma)
  # # v2_hat=t(v2)%*%kalpha/colSums(kalpha)
  # # v2_hat=as.vector(v2_hat)
  list(df_fixed=df_fixed,esigma=esigma,beta_dev2=beta_dev2)
}
