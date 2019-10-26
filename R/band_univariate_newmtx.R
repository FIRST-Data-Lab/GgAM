#' @importFrom  stats runif
band_univariate_newmtx <- function (object,mhat,X,alpha,alpha0,smooth.univariate,h_opt,x0=NULL,
                               Ind=NULL,testlevel=0.05,...) {

  #select optimal bandwidth
  if (is.null(Ind)){
    Ind=1:length(x0)
  }
 if (is.null(x0)){
   x0=seq(-0.5,0.5,0.01)

   #x0=seq(0.01,0.99,0.01)
 }
  Band=NULL
  sigma_2=object$sigma_2

  #### back-fitting estimation ####
  y=object$y
  if (!is.null(object$coefficients_bivariate)){
    ind=object$basis_info[[length(object$basis_info)]]$ind
    y=y[ind]
  }
  n=dim(mhat)[1]

  if (is.null(h_opt)){
    XX=if (is.null(dim(X))){X} else {X[,alpha0]}
    #### need to know the original function???
    if (substr(object$family$family[1], 1, 17) == "Negative Binomial") {
      fun1=key_functionsSBLL(y,X,alpha0,mhat,sigma_2,XX,family=object$family
                             ,basisinfo=object$basis_info[[smooth.univariate[alpha]]],theta=object$est_theta)
    } else { #could delete
      fun1=key_functionsSBLL(y,X,alpha0,mhat,sigma_2,XX,family=object$family
                             ,basisinfo=object$basis_info[[smooth.univariate[alpha]]])
    }
    h_opt=
      length(y)^(-0.2)*(5*7*sum(fun1$esigma/fun1$df_fixed)/sum(fun1$beta_dev2^2))^(0.2)
    if(h_opt>0.8) h_opt=0.8
  }

  initial=runif(length(x0))-0.5
  if (substr(object$family$family[1], 1, 17) == "Negative Binomial") {
    theta <- object$est_theta
  } else {
    theta=0
  }
  m1_sbk=SBK_locp(y,X,initial,alpha0,x0,as.matrix(mhat),h_opt,family=object$family,theta)

  #confidence band
  aa=key_functionsSBLL(y,X,alpha0,mhat,sigma_2,x0[Ind],family=object$family
                       ,basisinfo=object$basis_info[[smooth.univariate[alpha]]])
  a_h=sqrt(-2*log(h_opt))
  q_testlevel=a_h+(log(3^0.5)-log(2*pi)-log(-log(1-testlevel)/2))/a_h
  upper=q_testlevel/sqrt(n*h_opt)*sqrt(1/aa$esigma/aa$df_fixed*5/7)+m1_sbk[Ind]
  lower=-q_testlevel/sqrt(n*h_opt)*sqrt(1/aa$esigma/aa$df_fixed*5/7)+m1_sbk[Ind]
  #Band=c(Band,prod(tt[Ind]>=lower & tt[Ind]<=upper)==1)
  # }

  ret<-list(upper=upper,lower=lower,x0=x0[Ind],est=m1_sbk[Ind],m1_sbk=m1_sbk,label=
              object$basis_info[[smooth.univariate[alpha]]]$label,
            term=object$basis_info[[smooth.univariate[alpha]]]$term,h_opt=h_opt)
  class(ret)<-"band.univariate"
  ret
} ## end
