#' Penalized Least Square Fit under CV
#'
#' This is an internal function of package \code{ggam}.

#' @importFrom MASS ginv
#' @importFrom Matrix chol
#' @param B The bernstein basis matrix.
#' @param Q2 The \code{Q2} matrix from QR decomposition of the transpose of the constraint matrix.
#' @param K The energy matrix to construct penalty matrix.
#' @param lambda The smoothing penalty parameter.
#' @param Y Response variable.
#' @param fx indicates whether the term is a fixed d.f. regression.
#'   spline (\code{TRUE}) or a penalized regression spline (\code{FALSE}).
#' @param fold number of folders to do cross validation.
#' @param Z The parametric model matrix. set to '\code{NULL}' if it is not provided.

#' @details
#' The method is a computationally efficient means of applying cross validation to the problem of smoothing parameter selection:
#' \deqn{\min _ { \boldsymbol { \beta } , \boldsymbol { \gamma } } \frac { 1 } { 2 } \left\{ \| \mathbf { Y } - \mathbf { Z } \boldsymbol { \beta } - \mathbf { B } \boldsymbol { \gamma } \| ^ { 2 } + \lambda \boldsymbol { \gamma } ^ { \top } \mathbf { P } \gamma \right\}}{\min _ { \boldsymbol { \beta } , \boldsymbol { \gamma } } \frac { 1 } { 2 } \left\{ \| \mathbf { Y } - \mathbf { Z } \boldsymbol { \beta } - \mathbf { B } \boldsymbol { \gamma } \| ^ { 2 } + \lambda \boldsymbol { \gamma } ^ { \top } \mathbf { P } \gamma \right\}}
#' subject to constraints \eqn{\mathbf { H } \gamma = \mathbf { 0 }}{\mathbf { H } \gamma = \mathbf { 0 }}.
#' \code{Z} is a parametrix design matrix, \eqn{\beta}{\beta} a parameter vector, \eqn{Y}{Y} a data vector,
#' \eqn{\gamma}{\gamma} is the berstein coefficients, \eqn{B}{B} is the Bernsterin basis matrix,
#' \eqn{H}{H} is contraint matrix.
#' @return
#' A list of fit information.
#'
#' @examples
#'library(GgAM)
#'library(Matrix)
#'library(BPST)
#'data("eg2pop_dat")
#'eg2_V20=eg2pop_dat[['V20']]
#'eg2_T20=eg2pop_dat[['T20']]
#'eg2pop=eg2pop_dat[['pop']]
#'d=2
#'r=1
#'sam=eg2pop[sample(1:dim(eg2pop)[1],100),]
#'B0=basis(eg2_V20,eg2_T20, d, r, sam[,3:4])
#'B=B0$Bi
#'ind=B0$Ind.inside
#' Q2=B0$Q2
#' K=B0$K
#' Z=sam[ind,c(5:12)]
#' Y=sam[ind,'Y']
#' lambda=10^(seq(-2,5,by=1))
#' plsfitCV(as.matrix(B),Q2,K,lambda,Y,fx=FALSE,Z=Z)
#'### without parametric part
#' plsfitCV(as.matrix(B),Q2,K,lambda,Y,fx=FALSE)
#' @export
######################################################################
plsfitCV=function(B,Q2,K,lambda,Y,fx,fold=5,Z=NULL){
  n=length(Y)
  J=ncol(Q2)
  if(!is.null(Z)){
    np=ncol(Z)
    BQ2=B%*%Q2
    W=cbind(BQ2,Z)
    WW=t(W)%*%W
    rhsorig=t(W)%*%Y
    VV=WW+diag(c(rep(10^-10,J),rep(0,np)))
    Ainv=chol(VV)
    A=solve(t(Ainv))

    P=as.matrix(t(Q2)%*%K%*%Q2)
    D=matrix(0,(J+np),(J+np))
    D[1:J,1:J]=P
    ADA=A%*%D%*%t(A)
    eigs=eigen(ADA)
    C=eigs$values
    #alpha_all=matrix(rep(0,(d1+d2)*nl),ncol=nl)
  } else{
    np=0
    W=B%*%Q2
    W=as.matrix(W)
    WW=t(W)%*%W
    rhsorig=t(W)%*%Y
    VV=WW+diag(rep(10^-10,J))
    Ainv=chol(VV)
    A=solve(t(Ainv))

    D=as.matrix(t(Q2)%*%K%*%Q2)
    ADA=A%*%D%*%t(A)
    eigs=eigen(ADA)
    C=eigs$values
    #alpha_all=matrix(rep(0,d2*nl),ncol=nl)
  }
  if (fx==TRUE){
    Lam=NULL
    lhs=WW#+diag(c(rep(10^-10,J),rep(0,np)))
    theta=solve(lhs)%*%rhsorig
    if(np!=0) {
      alpha=theta[-(1:J)]
      theta=theta[1:J]
      gamma=Q2%*%theta
      beta=B%*%gamma
      Yhat=beta+Z%*%alpha
    } else {
      alpha=NULL
      gamma=Q2%*%theta
      beta=B%*%gamma
      Yhat=beta
    }
    res=Y-Yhat
    sse=sum((Y-Yhat)^2)
    df=sum(1/(1+C*0))
    dfs=1/(1+C*0)
    bic=log(sse/n)+df*log(n)/n
    cv=NULL
  }

  if (fx==FALSE){
    lambda=as.matrix(lambda)
    nl=nrow(lambda)
    cv_all=c()
    for(il in 1:nl){
      indp=sample(n,n,replace=FALSE)
      nk=fold
      sk=round(n/nk)
      Y.sspe=0
      for(ik in 1:nk){
        if(ik==nk){
          indk=sort(indp[((ik-1)*sk+1):n])
        }else{
          indk=sort(indp[((ik-1)*sk+1):(ik*sk)])
        }
        indnk=sort(setdiff(indp,indk))
        Ynk=Y[indnk]
        Bnk=B[indnk,]
        Znk=Z[indnk,]
        BQ2=Bnk%*%Q2
        if (!is.null(Z)){
          Znk=Z[indnk,]
          W=cbind(BQ2,Znk)
          WWnk=t(W)%*%W
          rhs=t(W)%*%Ynk
        } else {
          W=Bnk%*%Q2
          WWnk=t(W)%*%W
          rhs=t(W)%*%Ynk
        }
        Lam=lambda[il]
        Dlam=Lam*D
        lhs=WWnk+Dlam
        theta=solve(lhs)%*%rhs
        if (np!=0){
          alphak=theta[-(1:J)]
          thetak=theta[1:J]
          gammak=Q2%*%thetak
          Yk=Y[indk]
          Bk=B[indk,]
          Zk=Z[indk,]
          mk.hat=Bk%*%gammak
          Yk.hat=Zk%*%alphak+mk.hat
        } else {
          thetak=theta
          gammak=Q2%*%thetak
          Yk=Y[indk]
          Bk=B[indk,]
          mk.hat=Bk%*%gammak
          Yk.hat=mk.hat
        }
        Y.sspe=Y.sspe+sum((Yk.hat-Yk)^2,na.rm=TRUE)
      }
      cv_all=c(cv_all,Y.sspe/n)
    }
    j=which.min(cv_all)
    cv=cv_all[j]

    Lam=lambda[j,]
    Dlam=Lam*D
    lhs=WW+Dlam
    theta=solve(lhs)%*%rhsorig
    if (np!=0){
      alpha=theta[-(1:J)]
      theta=theta[1:J]
      gamma=Q2%*%theta
      beta=B%*%gamma
      Yhat=beta+Z%*%alpha
    } else {
      alpha=NULL
      gamma=Q2%*%theta
      beta=B%*%gamma
      Yhat=beta
    }
    res=Y-Yhat
    sse=sum((Y-Yhat)^2)
    df=sum(1/(1+C*Lam))
    dfs=1/(1+C*Lam)
    gcv=n*sse/(n-df)^2
    bic=log(sse/n)+df*log(n)/n
  }
  list(alpha_hat=alpha,beta_hat=beta,gamma_hat=gamma,theta_hat=theta,sse=sse,cv=cv,bic=bic,lam0=Lam,Yhat=Yhat,res=res,edf=df, edfs=dfs)
}

