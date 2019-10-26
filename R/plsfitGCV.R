#' Penalized Least Square Fit under GCV
#'
#' This is an internal function of package \code{ggam}.

#'
#' @param B The bernstein basis matrix.
#' @param Q2 The \code{Q2} matrix from QR decomposition of the transpose of the constraint matrix.
#' @param P The penalty matrix.
#' @param lambda The smoothing penalty parameter.
#' @param Y Response variable.
#' @param fx indicates whether the term is a fixed d.f. regression
#' spline (\code{TRUE}) or a penalized regression spline (\code{FALSE}).
#' @param Z The parametric model matrix. set to '\code{NULL}' if it is not provided.
#' @param ... other arguments.

#' @details
#' The method is a computationally efficient means of applying \code{GCV} to the problem of smoothing parameter selection:
#' \deqn{\min _ { \boldsymbol { \beta } , \boldsymbol { \gamma } } \frac { 1 } { 2 } \left\{ \| \mathbf { Y } - \mathbf { Z } \boldsymbol { \beta } - \mathbf { B } \boldsymbol { \gamma } \| ^ { 2 } + \lambda \boldsymbol { \gamma } ^ { \top } \mathbf { P } \gamma \right\}}{\min _ { \boldsymbol { \beta } , \boldsymbol { \gamma } } \frac { 1 } { 2 } \left\{ \| \mathbf { Y } - \mathbf { Z } \boldsymbol { \beta } - \mathbf { B } \boldsymbol { \gamma } \| ^ { 2 } + \lambda \boldsymbol { \gamma } ^ { \top } \mathbf { P } \gamma \right\}}
#' subject to constraints \eqn{\mathbf { H } \gamma = \mathbf { 0 }}{\mathbf { H } \gamma = \mathbf { 0 }}.
#' \code{Z} is a parametrix design matrix, \eqn{\beta}{\beta} a parameter vector, \eqn{Y}{Y} a data vector,
#' \eqn{\gamma}{\gamma} is the berstein coefficients, \eqn{B}{B} is the Bernsterin basis matrix,
#' \eqn{H}{H} is contraint matrix.
#' @return
#' A list of fit information.
#'
#' @examples
#' library(Matrix)
#' library(BPST)
#' data("eg1pop_dat")
#' eg1_V1=eg1pop_dat[['V1']]
#' eg1_T1=eg1pop_dat[['T1']]
#' eg1pop_rho03=eg1pop_dat[['rho03']]
#' sam=eg1pop_rho03[sample(1:dim(eg1pop_rho03)[1],100),]
#' B0=basis(eg1_V1,eg1_T1, d=2, r=1, sam[,3:4])
#' B=B0$Bi
#' ind=B0$Ind.inside
#' Q2=B0$Q2
#' K=B0$K
#' P=t(Q2)%*%K%*%Q2
#' Z=sam[ind,c(5:12)]
#' Y=sam[ind,'Y']
#' lambda=10^(seq(-2,5,by=1))
#' plsfitGCV(as.matrix(B),Q2,P,lambda,Y,fx=FALSE,Z=Z)
#'
#'### without parametric part
#' plsfitGCV(as.matrix(B),Q2,P,lambda,Y,fx=FALSE)
#' @export

######################################################################
plsfitGCV=function(B,Q2,P,lambda,Y,fx,Z=NULL,...){
  
  if (!is.null(B)){
    BQ2 <- as.matrix(B%*%Q2)
    J=ncol(Q2)
  } else {BQ2 <- NULL}
  n=length(Y)

  if (!is.null(BQ2)){
    #if(!is.null(Z)){
    #BQ2=B%*%Q2
    W=cbind(as.matrix(BQ2),Z)
    W=as.matrix(W)
    WW=t(W)%*%W
    rhs=t(W)%*%Y
    if (!is.null(Z)){
      np=ncol(Z)
      VV=WW+diag(c(rep(10^-10,J),rep(0,np)))
      Ainv=chol(VV)
      A=solve(t(Ainv))
      #P=as.matrix(t(Q2)%*%K%*%Q2)
      D=matrix(0,(J+np),(J+np))
      # print(dim(D))
      # print(dim(P))
      D[1:J,1:J]=as.matrix(P)
    } else {
      np=0
      W=B%*%Q2
      W=as.matrix(W)
      WW=t(W)%*%W
      rhs=t(W)%*%Y
      VV=WW+diag(rep(10^-10,J))
      Ainv=chol(VV)
      A=solve(t(Ainv))
      #D=as.matrix(t(Q2)%*%K%*%Q2)
      D=P
    }
    ADA=A%*%D%*%t(A)
    eigs=eigen(ADA)
    C=eigs$values
  } else {
    W=Z
    W=as.matrix(W)
    WW=t(W)%*%W
    rhs=t(W)%*%Y
    if (!is.null(Z)){
      np=ncol(Z)
      VV=WW
      # Ainv=chol(VV)
      # A=solve(t(Ainv))
      # P=as.matrix(t(Q2)%*%K%*%Q2)
      # D=matrix(0,(J+np),(J+np))
      # D[1:J,1:J]=P
    } else {
      np=0
      W=NULL
      rhs=NULL
    }
  }
   if (!is.null(lambda)){
    lambda=as.matrix(lambda)
    nl=nrow(lambda)
    alpha_all=c()
    theta_all=c()
    gamma_all=c()
    beta_all=c()
    res_all=c()
    sse_all=c()
    df_all=c()
    gcv_all=c()
    cv_all=c()
    bic_all=c()
    dfs_all=c()
    Yhat_all=c()
    for(il in 1:nl){
      Lam=lambda[il]
      Dlam=Lam*D
      lhs=WW+Dlam
      theta=solve(lhs)%*%rhs
      #
      if(np!=0) {
        alpha=theta[-(1:J)]
        theta=theta[1:J]
        gamma=Q2%*%theta
        beta=B%*%gamma
        Yhat=beta+as.matrix(Z)%*%as.matrix(alpha)
      } else {
        alpha=NULL
        gamma=Q2%*%theta
        beta=B%*%gamma
        Yhat=beta
      }
      alpha_all=cbind(alpha_all,alpha)
      theta_all=cbind(theta_all,theta)
      gamma_all=cbind(gamma_all,gamma)
      beta_all=cbind(beta_all,beta)
      Yhat_all=cbind(Yhat_all,Yhat)
      res=Y-Yhat
      res_all=cbind(res_all,res)
      sse=sum((Y-Yhat)^2)
      sse_all=c(sse_all,sse)
      df=sum(1/(1+C*Lam))
      dfs=1/(1+C*Lam)
      dfs_all=cbind(dfs_all,dfs)
      df_all=c(df_all,df)
      gcv=n*sse/(n-df)^2
      tem=1/(1+C*Lam)
      gcv_all=c(gcv_all,gcv)
      bic=log(sse/n)+df*log(n)/n
      bic_all=c(bic_all,bic)
    }
    j=which.min(gcv_all)
    #j=which.min(bic_all)
    lambdac=lambda[j,]
    alpha=alpha_all[,j]
    theta=theta_all[,j]
    gamma=gamma_all[,j]
    beta=beta_all[,j]
    dfs=dfs_all[,j]
    sse=sse_all[j]
    Yhat=Yhat_all[,j]
    gcv=gcv_all[j]
    bic=bic_all[j]
    df=df_all[j]
   } else {
     lhs=WW
     theta=solve(lhs)%*%rhs
     #
     if(np!=0 && !is.null(BQ2)) {
       alpha=theta[-(1:J)]
       theta=theta[1:J]
       gamma=Q2%*%theta
       beta=B%*%gamma
       Yhat=beta+as.matrix(Z)%*%as.matrix(alpha)
     } else if (np==0 && !is.null(BQ2)){
       alpha=NULL
       gamma=Q2%*%theta
       beta=B%*%gamma
       Yhat=beta
     } else if (np!=0 && is.null(BQ2)){
       alpha=theta
       theta=NULL
       gamma=NULL
       beta=rep(0,length(Y))
       Yhat=as.matrix(Z)%*%as.matrix(alpha)
     } else{
       # everything is NULL
     }
     lambdac=NULL
     res=Y-Yhat
     sse=sum((Y-Yhat)^2)
     df=sum(1/(1+0))
     dfs=1/(1+0)
     gcv=n*sse/(n-df)^2
     bic=log(sse/n)+df*log(n)/n
  }
  # se_beta=rep(0,np)
  # if(se=="TRUE"){
  #   sb=beta_se(V,Tr,B,Q2,K,lambdac,lam2,X,Z,Y,method)
  #   se_beta=sb$se_beta
  #   Ve=sb$Ve
  # }
  list(alpha_hat=alpha,beta_hat=beta,gamma_hat=gamma,theta_hat=theta,sse=sse,gcv=gcv,bic=bic,lam0=lambdac,Yhat=Yhat,res=res,edf=df, edfs=dfs)
}

