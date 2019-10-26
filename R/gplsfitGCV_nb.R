#' Generalized Penalized Least Square Fit under GCV for negative binomial family
#'
#' This is an internal function of package \code{ggam}.

#' @param B The bernstein basis matrix.
#' @param Q2 The \code{Q2} matrix from QR decomposition of the transpose of the constraint matrix.
#' @param P The penalty matrix.
#' @param UB The univariate basis function matrix constructed.
#' @param lambda The smoothing penalty parameter.
#' @param family The family object, specifying the distribution and link to use.
#' @param offset Can be used to supply a model offset for use in fitting. Note that this offset
#' will always be completely ignored when predicting.
#' @param Y Response variable.
#' @param fx indicates whether the term is a fixed d.f. regression
#' spline (\code{TRUE}) or a penalized regression spline (\code{FALSE}).
#' @param X The parametric model matrix. set to '\code{NULL}' if it is not provided.
#' @param control A list of fit control parameters to replace defaults returned by \code{\link{plbpsm.control}}.
#' Any control parameters not supplied stay at their default values.
#' @param r.theta All the theta values given.
#' @param ind_c The vector of index to indicate the parametric part.
#' @param fixedSteps How many steps to take: useful when only using this routine to get rough starting
#' values for other methods.
#' @param ... other arguments passed onto \code{\link{gplsfitGCV}}.

#' @details
#' In this function, the estimator of \eqn{\theta}{\theta} is chosen to ensure that
#' the Pearson estimate of the scale parameter is as close as possible to 1. The other parts follow from the
#' routine of \code{\link{gplsfitGCV}}.
#' @return
#' A list of fit information.
#' @export
gplsfitGCV_nb=function(Y,B,Q2,P,UB,lambda,family,offset,r.theta=c(2,8),fx,
                    control,X=NULL,ind_c=1:ncol(X),fixedSteps=(control$maxstep+1),...){
  if (length(r.theta)==1){
    result=gplsfitGCV(Y,B,Q2,P,UB,lambda,family,offset,r.theta,fx=fx,control=control,X=X,ind_c=ind_c,...)
    middle=r.theta
  } else {
    lower=r.theta[1]
    upper=r.theta[2]
    result1=gplsfitGCV(Y,B,Q2,P,UB,lambda,family,offset,upper,fx=fx,control=control,X=X,ind_c=ind_c,...)
    y_hat=result1$Yhat
    variance=function(mu,theta) mu+mu^2/theta
    tt=(Y-y_hat)^2/variance(y_hat,upper)
    U_sigma=1/(length(Y)-result1$df)*sum(tt)

    result2=gplsfitGCV(Y,B,Q2,P,UB,lambda,family,offset,lower,fx=fx,control=control,X=X,ind_c=ind_c,...)
    y_hat=result2$Yhat
    tt=(Y-y_hat)^2/variance(y_hat,lower)
    L_sigma=1/(length(Y)-result2$df)*sum(tt)


    M_sigma=2
    step=1

    if ( (L_sigma-1)>0 & (U_sigma-1) >0 ) {
      M_sigma=1
      result=result2
      middle=lower
    }

    if ( (L_sigma-1)<0 & (U_sigma-1) < 0 ) {
      M_sigma=1
      result=result1
      middle=upper
    }
    while(abs(M_sigma-1)>0.01 & step <= 10){
      middle=(upper+lower)/2
      result=gplsfitGCV(Y,B,Q2,P,UB,lambda,family,offset,middle,fx=fx,control=control,X=X,ind_c=ind_c,...)
      y_hat=result$Yhat
      tt=(Y-y_hat)^2/variance(y_hat,middle)
      M_sigma=1/(length(Y)-result$df)*sum(tt)
      if((M_sigma-1)*(U_sigma-1)<0){
        lower=middle
        L_sigma=M_sigma
      } else {
        upper=middle
        U_sigma=M_sigma
      }
      print(c(middle,M_sigma))
      step=step+1
    }
  }

  list(boundary=result$boundary,alpha_hat=result$alpha_hat,theta_hat=result$theta_hat,
       gamma_hat=result$gamma_hat,beta_hat=result$beta_hat,
       lam0=result$lam0,gcv=result$gcv,df=result$df,est_theta=middle,middle=middle,
       eta_hat=result$eta_hat,Yhat=result$Yhat,sse=result$sse,res=result$res,w=result$w,se_beta=result$se_beta, z=result$z,Ve=result$Ve)
}
