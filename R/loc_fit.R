loc_fit=function(iter,y,ualpha,kalpha,thetaalpha,family,theta=0){
  u=ualpha[,iter]
  K=kalpha[,iter]
  Z=cbind(1,u)
  b=c(1,1)
  b_new=c(0,0)
  step=0

  linkinv=family$linkinv
  mu.eta=family$mu.eta
  if(theta==0) {variance=family$variance} else {
    variance=function(mu,theta) mu+mu^2/theta
  }
  while(sum((b_new-b)^2)>0.00001 & step <= 10){
    step=step+1
    b=b_new
    eta=Z%*%b+thetaalpha
    mu=linkinv(eta)
    mevg=mu.eta(eta)
    if(theta==0) var=variance(mu) else
      var=variance(mu,theta)
    Y_iter=(eta-thetaalpha)+(y-mu)/mu.eta(eta)
    W_iter=as.vector((mevg^2)/var*K)
    temp1=W_iter*Z
    temp2=W_iter*Y_iter
    temp3=crossprod(Z,temp1)
    b_new=solve(temp3,crossprod(Z,temp2))
    if(step==0) print('not convergent')
  }
  b_new[1]
}
