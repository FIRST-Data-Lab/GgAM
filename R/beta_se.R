#' Standard Error Estimation for Linear Coefficients
#'
#' This is an internal function of package \code{ggam}.
#' @importFrom psych tr
#' @importFrom MASS ginv
#' @param B The bernstein basis matrix.
#' @param Q2 The \code{Q} matrix from QR decomposition of the constraint matrix.
#' @param K The energy matrix.
#' @param UB The univariate basis function matrix constructed.
#' @param lam1 The smoothing penalty parameter.
#' @param lam2 The variable selection penalty parameter.
#' @param Z The parametric matrix.
#' @param Y Response variable.
#' @param ind.c The indexed for parametric coefficient for which standard error want to be calculated.
#' @param VS variable section is conducted or not.
#' @param ... other arguments.
#' @details A sandwich formula is developed to find the standard error for \eqn{\beta}{\beta}.
#' The detailed algorithm is in the paper Wang et al. (2018).
#'
#' @return
#' \item{se_beta}{The standard error of linear coefficients.}
#' \item{Ve}{The estimated covariance matrix.}
#'
#' @examples
#' library(MASS)
#' library(grpreg)
#' library(Matrix)
#' library(BPST)
#' data("eg1pop_dat")
#' eg1_V1=eg1pop_dat[['V1']]
#' eg1_T1=eg1pop_dat[['T1']]
#' eg1pop_rho03=eg1pop_dat[['rho03']]
#' sam=eg1pop_rho03[sample(1:dim(eg1pop_rho03)[1],100),]
#' B0=basis(eg1_V1,eg1_T1, d=2, r=1, sam[,3:4])
#' B=B0$B
#' ind=B0$Ind.inside
#' Q2=B0$Q2
#' K=B0$K
#' Z=sam[ind,c(5:12)]
#' Y=sam[ind,'Y']
#' lam1=0.1
#' lam2=0
#' beta_se(UB=NULL,as.matrix(B),Q2,K,lam1,lam2,Z,Y,ind.c=1:6)
#' @export

beta_se=function(UB,B,Q2,K,lam1,lam2,Z,Y,ind.c,VS=FALSE,...){
  if (VS){
    Xp=Z
  } else {
    Xp=Z[,ind.c,drop=FALSE]}
  if (!is.null(B) && !is.null(Xp) && is.null(UB)){
    lam1=ifelse(length(lam1)>0,lam1,0)
    n=length(Y)
    B=as.matrix(B)
    temp1=t(Q2)%*%(as.matrix(t(B)%*%B)+as.matrix(lam1*K))%*%Q2
    temp1=as.matrix(temp1)
    temp2=ginv(temp1)
    H_B=B%*%Q2%*%temp2%*%t(Q2)%*%t(B)
    temp1=t(Xp)%*%(diag(dim(H_B)[1])-H_B)%*%Xp
    temp2=solve(temp1,t(Xp))
    Xp_hat=as.matrix(H_B)%*%as.matrix(Xp)
    Hmtx=t(Xp-Xp_hat)%*%as.matrix(Xp)
    #Sigma_lambda=scadP(alpha,lam2,method)$Sigma_lambda
    rhs=Hmtx#+n*Sigma_lambda
    rhs=as.matrix(rhs)
    temp3=as.matrix(Xp-Xp_hat)%*%ginv(rhs)%*%t(Xp-Xp_hat)
    Slambda=temp3+H_B
    Y_hat=Slambda%*%Y
    temp4=sum((Y-Y_hat)^2)
    Slambda=as.matrix(Slambda)
    temp5=tr(Slambda)
    temp6=n-temp5
    #sigma2_hat=temp4/temp6
    se_beta=sqrt(diag(tcrossprod(temp2)))
    #Ve=sigma2_hat*ginv(rhs)%*%
    #t(Xp-Xp_hat)%*%as.matrix(Xp-Xp_hat)%*%ginv(rhs)
    Ve=ginv(rhs)%*%
      t(Xp-Xp_hat)%*%as.matrix(Xp-Xp_hat)%*%ginv(rhs)
    # print(se_beta)
    # print(Ve)
    print('aa')
    
  } else if (!is.null(B) && !is.null(Xp) && !is.null(UB)) {
    BQ2=as.matrix(B%*%Q2)
    BB <- cbind(UB,BQ2)
    V22 <- crossprod(BB,BB)
    P <- t(Q2)%*%K%*%Q2
    v22_inv <-solve(V22+adiag(matrix(0,ncol=ncol(UB),nrow=ncol(UB)),
                              as.matrix((lam1*P))))
    temp3 <- Xp
    #print(dim(Xp))
    Xphat <- as.matrix(BB)%*%v22_inv%*%crossprod(as.matrix(BB),as.matrix(temp3))
    temp4 <- (Xp-Xphat)
    
    Ve <- solve(crossprod(as.matrix(Xp-Xphat),as.matrix(temp4))) #version 1
    se_beta <- sqrt(diag(Ve))
    
    #print(se_beta)
    #se_beta <- Ve <- NULL
  } else if (!is.null(Xp) && !is.null(UB) && is.null(B)){
    # V22 <- crossprod(UB,UB)
    # v22_inv <- solve(V22)
    # temp3 <- Xp
    # Xphat <- as.matrix(UB)%*%v22_inv%*%crossprod(as.matrix(UB),as.matrix(temp3))
    # temp4 <- (Xp-Xphat)
    # Ve <- solve(crossprod(as.matrix(Xp-Xphat),as.matrix(temp4))) #version 1
    # se_beta <- sqrt(diag(Ve))
    se_beta <- Ve <- NULL
  } else if (!is.null(Xp) && is.null(UB) && is.null(B)){
    Ve <- t(Xp)%*%Xp
    se_beta <- sqrt(diag(Ve))
  } else {
    se_beta <- Ve <- NULL
  }
  return(list(se_beta=se_beta,Ve=Ve))
}