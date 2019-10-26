#' @importFrom  stats runif
band_univariate_orig <- function (object,mhat,X,alpha,alpha0,smooth.univariate,h_opt,
                                  Ind=NULL,testlevel=0.05,...) {

  sigma_2=object$sigma_2

  #### back-fitting estimation ####
  y=object$y
  if (!is.null(object$coefficients_bivariate)){
    ind=object$basis_info[[length(object$basis_info)]]$ind
    #X=X[ind,]
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
  initial=runif(n)-0.5

  if (is.null(Ind) && length(initial)>10000){
    Ind=seq(1,length(initial),by=30)
  } else {Ind=seq(1,length(initial))}
  if (substr(object$family$family[1], 1, 17) == "Negative Binomial") {
    theta <- object$est_theta
  } else {
    theta=0
  }

  m1_sbk=SBK_locp(y,X,initial[Ind],alpha0,xfixed=X[Ind,alpha0],as.matrix(mhat),h_opt,family=object$family,theta)

  #confidence band
  aa=key_functionsSBLL(y,X,alpha0,mhat,sigma_2,X[Ind,alpha0],family=object$family
                       ,basisinfo=object$basis_info[[smooth.univariate[alpha]]])
  a_h=sqrt(-2*log(h_opt))
  q_testlevel=a_h+(log(3^0.5)-log(2*pi)-log(-log(1-testlevel)/2))/a_h
  upper=q_testlevel/sqrt(n*h_opt)*sqrt(1/aa$esigma/aa$df_fixed*5/7)+m1_sbk
  lower=-q_testlevel/sqrt(n*h_opt)*sqrt(1/aa$esigma/aa$df_fixed*5/7)+m1_sbk
  #Band=c(Band,prod(tt[Ind]>=lower & tt[Ind]<=upper)==1)
  # }
  nameshere=object$basis_info[[smooth.univariate[alpha]]]$term
  B0=object$basis_info[[length(object$basis_info)]]
  dat3=object$model[B0$ind,]
  dat4=dat3[Ind,]
  Index=order(dat4[,nameshere])

  ret<-list(Ind=Ind,Index=Index,upper=upper[Index],lower=lower[Index],x0=dat4[Index,nameshere],est=m1_sbk[Index],m1_sbk=m1_sbk,label=
              object$basis_info[[smooth.univariate[alpha]]]$label,dat=dat4,
            term=object$basis_info[[smooth.univariate[alpha]]]$term,h_opt=h_opt)
  class(ret)<-"band.univariate"
  ret
} ## end
