#' Variable selection in Partial Linear Bivariate penalized Spline fit in generalized spatial model
#'
#' This is an internal function of package \code{GgAM}. Bivariate penalized least squares problem
#' is solved with penalty parameter chosen by GCV or CV. Variable selection by using
#' adaptive LASSO or group SCAD is applied in parametric coefficients.
#'
#' @importFrom grpreg grpreg
#' @importFrom stats gaussian model.frame glm coef runif

#' @param G An object of the type returned by \code{plbpsm} when \code{fit=FALSE}.
#' @param criterion The criterion to choose the penalty parameter lambda. \code{"GCV"} to use
#' generalized cross validation method and \code{"CV"} for cross validation
#' @param family The family object, specifying the distribution and link to use.
#' @param method 'ALASSO' or 'SCAD' to penalize the coefficients for parametric part.
#' @param ind_c The given index of covariates that are selected.
#' @param VS '\code{TRUE}' for using ALASSO/SCAD to select linear variables.
#' @param control A list of fit control parameters to replace defaults returned by \code{\link{plbpsm.control}}.
#' Any control parameters not supplied stay at their default values.
#' @param MI whether model identification is conducted or not.
#' @param ... other arguments.
#' @details
#' This is an internal function of package \code{GgAM}. We propose Iteratively Reweighted Least square based
#' algorithm to get the poilot estimation and then use it to get a a spline-backfitted local polynomial estimation.
#' The smoothing penalty parameter could
#' be chosen by \code{GCV} or \code{CV} using the routines: \code{\link{gplsfitGCV}}.
#' @return A list of fit information.
#' @export
ggrplsfit <- function(G, criterion, method, family, ind_c, VS, control = plbpsm.control(), MI,...)
{
  # Deal with NB
  if (substr(family$family[1], 1, 17) == "Negative Binomial") {
    theta <- family$getTheta()
    if (length(theta) == 1) {
      G$sig2 <- 1
    } else {
      if (length(theta) > 2)
      find.theta <- TRUE
    }
    nb.link <- family$link
    gambpsname <- 'gplsfitGCV_nb'
  } else {
    gambpsname <- 'gplsfitGCV'
    theta <- 0
  }

  Terms <- delete.response(G$pterms)
  Xpp <- model.matrix(Terms,G$mf,contrasts=G$contrasts)
  if (length(G$basis_info)>0){
    tem <- "bivariate.smooth"%in%sapply(G$basis_info, class)
  }else {
    Gb <-NULL
    tem <- FALSE}
    if (tem){
      Gb <- G$basis_info[[length(G$basis_info)]]
      V <- Gb$V
      Tr <- Gb$Tr
      J <- ncol(Gb$Q2)
      ### standardize before doing ALASSO/SCAD
      B <- Gb$B
      Q2 <- Gb$Q2
      K <- Gb$K
      P <- t(Q2)%*%K%*%Q2
      lambda <- G$lambda
      Xpp <- Xpp[,-1,drop=FALSE]
      Xp <- Xpp[Gb$ind,,drop=FALSE]
    }else {
      BQ2 <- NULL
      B <- NULL
      Q2 <- NULL
      P <- NULL
      Gb <-NULL
      Xp <- Xpp
      lambda <- NULL}
  nvars <- NCOL(G$X)
  y <- G$y # original data

  X <- G$X # original design matrix
  if (MI){
    X <- cbind(G$Xp,X)
  }
  if (nvars == 0)
    stop("Model seems to contain no terms")
  if (NCOL(y) > 1)
    stop("y must be univariate unless binomial")

  if ((G$intercept & G$m_b!=0)){
    X <- X[,-1]
    if(is.null(dim(X))){
      X <- as.matrix(X)
    }
  } else {X <- X}
  if (length(X)==0){VS=FALSE}


  # univariate basis info
  UB <- c()
  if (MI){
    if (length(G$ind.nl)>0){
      UB <- G$UB
    } else {UB <- NULL}
  } else {
    if (length(G$basis_info)>0){
      for (j in 1:length(G$basis_info)){
        if (class(G$basis_info[[j]])[1]=='univariate.smooth'){
          UB <- cbind(UB,G$basis_info[[j]]$B)
        }}
    } else {UB <- NULL}
  }

  offset <- G$offset
  weights <- G$w
  n.score <- sum(weights != 0)
  ### Define some functions to use in 'family'
  variance <- family$variance
  dev.resids <- family$dev.resids
  aic <- family$aic
  linkinv <- family$linkinv
  linkfun <- family$linkfun
  mu.eta <- family$mu.eta

  if (NCOL(y) > 1)
    stop("y must be univariate unless binomial")

  if (!is.function(variance) || !is.function(linkinv))
    stop("illegal `family' argument")

  gplsfit.control<-list(trace=control$trace,maxstep=control$maxstep,epsilon=control$epsilon)

  ### main algorithm for grplsfit starts here

  Y <- y
  n <- nobs <- length(Y)

  good <- weights > 0
  if (all(!good)) {
    warning(gettextf("No observations informative."))
  }
  if (!is.null(dim(X))){
    X <- X[good,,drop = FALSE]
    np <- ncol(X)
    #X <- apply(X,2,scale)
  } else {
    X <- X[good,drop = FALSE]
    np <- ncol(as.matrix(X))
    #X <- apply(as.matrix(X),2,scale)
  }
  Xp <- cbind(Xp,G$linear_cov)
  #print(head(Xp))
  colnames(Xp)<-c(colnames(Xpp),G$linear.names)
  #print(dim(Xp))

  # start of GCV
  if (criterion=='GCV'){
    # No variable selection
    if(VS==FALSE){
      ### no selection of variables so no ind.c
      if (length(Xp)>0){
        ind.c <- 1:(dim(Xp)[2])
      } else {
          ind.c <- 0}
        lam2 <- 0
        mfit <- do.call(gambpsname,list(y,B,Q2,P,UB,lambda,family,offset,theta,fx=G$fx,control=gplsfit.control,X=X,ind_c=ind.c,...))
        print(ind.c)
        print(mfit$se_beta)
        lam1 <- mfit$lam0
        tt <- (y-mfit$Yhat)^2/variance(mfit$Yhat)
        sigma_2 <- 1/(length(y)-mfit$df)*sum(tt)
        alpha_hat <- mfit$alpha_hat
        beta_hat <- mfit$beta_hat
        gamma_hat <- mfit$gamma_hat
        theta_hat <- mfit$theta_hat
        sse <- mfit$sse
        res <- mfit$res
        gcv <- mfit$gcv
        cv <- NULL
        Yhat <- mfit$Yhat
        edf <- mfit$df
        if (length(Xp[,ind.c])>0 && length(mfit$alpha_hat)>0){
            se_beta <- mfit$se_beta[1:length(ind.c)]
            Ve <- mfit$Ve[1:length(ind.c),1:length(ind.c)]
          } else {
            se_beta <- NULL
            Ve <- NULL
          }
    }

    ### variable selection
    if(VS==TRUE){
     # add more in future
    }
  }# end of GCV


  # get the final estimation here
  eta <- mfit$eta_hat

  # sbl estimation starts here
 if (G$backfitting){
   term.smooth <- sapply(G$basis_info, class)[1,]
   if (MI){
     #term.smooth <- sapply(object$basis_info2, class)[1,]
     smooth.univariate <- G$ind.nl
   } else {
     smooth.univariate <- which(term.smooth=="univariate.smooth")
   }
   smooth.bivariate <- which(term.smooth=="bivariate.smooth")
   X3 <- X2 <- Xp
   first.para <- dim(as.matrix(X2))[2]
   mhat <- as.matrix(X2)*alpha_hat[1:(dim(as.matrix(X2))[2])]

   ### univariate part
   if (length(smooth.univariate)>0){
     for (k in smooth.univariate){
       bsinf <- G$basis_info[[k]]
       last.para <- first.para+bsinf$n.para
       tem <- bsinf$B%*%alpha_hat[(first.para+1):last.para]
       temm <- G$mf[,bsinf$term]
       first.para <- last.para
       if (!is.null(Gb)){
         mhat <- cbind(mhat,tem)
         # G$Xorig
         X2 <- cbind(X2,temm[Gb$ind])} else {
           mhat <- cbind(mhat,tem)
           # G$Xorig
           X2 <- cbind(X2,temm)
         }
     }
   }
   ### bivariate part
   if(length(smooth.bivariate)>0){
     bsinf <- G$basis_info[[smooth.bivariate]]
     ind <- bsinf$ind
     temm <- bsinf$B%*%gamma_hat
     mhat <- cbind(mhat,temm)
   }

    first.para <- alpha0 <- dim(as.matrix(X3))[2]
    mhat.sbl <- c()
    h_opt_all <- c()
    if (length(smooth.univariate)>0){
      for (k in 1:length(smooth.univariate)){
        alpha0 <- first.para+k
        XX <- if (is.null(dim(X2))){X2} else {X2[,alpha0]}
        fun1 <- key_functionsSBLL(y,X2,alpha0,mhat,sigma_2,XX,family=family
                               ,basisinfo=G$basis_info[[smooth.univariate[k]]])
        h_opt <-
          length(y)^(-0.2)*(5*7*sum(fun1$esigma/fun1$df_fixed)/sum(fun1$beta_dev2^2))^(0.2)
        x0 <- XX
        initial <- runif(length(x0))-0.5
        if (length(initial)>10000){
          Ind <- seq(1,length(initial),by=30)}  else {Ind <- seq(1,length(initial),by=1)}
        m_sbk <- SBK_locp(y,X2,initial[Ind],alpha0,x0[Ind],mhat,h_opt,family=family)
        mhat.sbl <- cbind(mhat.sbl,m_sbk)
        h_opt_all <- c(h_opt_all,h_opt)
      }
    }
    if (!is.null(dim(Xp))) {
      Xp2 <- Xp
      #print(dim(Xp2))
      para <- if (!is.null(dim(Xp2)) && dim(Xp2)[2]>0) {as.matrix(Xp2)*alpha_hat[1:(dim(as.matrix(Xp2))[2])]} else {NULL}
    } else {
      para <- NULL}
    eta_sbl <- rowSums(as.matrix(cbind(para[Ind],as.matrix(mhat.sbl),beta_hat[Ind])))
    # print(para)
    # print(mhat.sbl)
    # print(as.matrix(cbind(para[Ind],mhat.sbl,beta_hat[Ind])))
    y_sbl <- linkinv(eta_sbl)
    # eta_sbl <- rep(0,length(beta_hat))
    # y_sbl <- linkinv(eta_sbl)
  } else {
      y_sbl <- rep(NA,length(y))
      Ind <- 1:length(y)
      mhat.sbl <- NULL
      h_opt_all <- NULL
      mhat <- NULL
      X3 <- X2 <- NULL
  }
  if (any(!is.finite(alpha_hat))) {
    warning(gettextf("Non-finite coefficients"))}
  eta <- as.matrix(eta)
  mu <- linkinv(eta <- eta + offset)
  eta <- linkfun(mu)
  dev <- sum(dev.resids(Y, mu, weights))
  if (mfit$boundary)
    warning("Algorithm stopped at boundary value")
  # check with specific family
  eps <- 10 * .Machine$double.eps
  if (family$family[1] == "binomial") {
    if (any(mu > 1 - eps) || any(mu < eps))
      warning("fitted probabilities numerically 0 or 1 occurred")
  }
  if (family$family[1] == "poisson") {
    if (any(mu < eps))
      warning("fitted rates numerically 0 occurred")
  }

  ### rewrite family with theta
  if (substr(family$family[1], 1, 17) == "Negative Binomial"){
      family<-do.call("negbin",list(theta=mfit$est_theta,link=nb.link))
  }

  residuals <- rep(NA, nobs)
  ### not sure
  residuals[good] <- mfit$z - (eta-offset)[good]
  ### For those points not in the triangles, return "NA"
  prior_weights_all <- G$y <- residuals_all <- weights_all <- eta_all <- fitted.values_all <- fitted.values_all_sbl <- rep(NA,dim(G$mf)[1])
  if (!is.null(Gb$ind)){
    G$y[Gb$ind]<-y
    fitted.values_all[Gb$ind]<-mu
    fitted.values_all_sbl3 <- fitted.values_all_sbl[Gb$ind]
    fitted.values_all_sbl3[Ind]<-y_sbl
    fitted.values_all_sbl[Gb$ind]<-fitted.values_all_sbl3
    residuals_all[Gb$ind] <- residuals
    #residuals_all <- NULL
    eta_all[Gb$ind] <- eta
    weights_all[Gb$ind] <- if (is.null(mfit$w)){rep(1,length(y))} else {mfit$w}
    prior_weights_all[Gb$ind] <- weights

  } else {
    G$y<-y
    fitted.values_all<-mu
    fitted.values_all_sbl3 <- fitted.values_all_sbl
    fitted.values_all_sbl3[Ind]<-y_sbl
    fitted.values_all_sbl<-fitted.values_all_sbl3
    residuals_all <- residuals
    eta_all <- eta
    weights_all <- if (is.null(mfit$w)){rep(1,length(y))} else {mfit$w}
    prior_weights_all<- weights
  }
  dev_sbl <- sum(dev.resids(Y, y_sbl, weights))
  wtdmu <- if (G$intercept)
    sum(weights * Y)/sum(weights) else linkinv(offset)
  nulldev <- sum(dev.resids(Y, wtdmu, weights))
  n.ok <- nobs - sum(weights == 0)
  nulldf <- n.ok - as.integer(G$intercept)
  rank2<-ifelse(is.null(Gb),0,dim(Gb$Q2)[1])
  # aic.model <- aic(Y, n, mu, weights, dev) + 2 * sum(G$edf)

  object <- list(coefficients = as.vector(alpha_hat), coefficients_bivariate = as.vector(gamma_hat),
              residuals = residuals_all,fitted.values.sbl=fitted.values_all_sbl,Xp=X3,dev_sbl=dev_sbl,#fitted.values_all_sbl2=fitted.values_all_sbl2,
              fitted.values = fitted.values_all, family = family, linear.predictors = eta_all,
              deviance = dev, null.deviance = nulldev, edf=edf, prior.weights=prior_weights_all, weights = weights_all,#edfs=edfs,
              df.null = nulldf,converged = mfit$conv,backfitting=G$backfitting,
              y = G$y, est_theta=mfit$est_theta,middle=mfit$middle,#scale=G$scale,
              r = G$r, nsdf = G$nsdf, Ve = Ve, iter=mfit$iter,
              gcv_opt= gcv, h_opt_all=h_opt_all,cv_opt=cv, # aic = aic.model,
              boundary = mfit$boundary,mhat.sbl=mhat.sbl,mhat=mhat,X2=X2,
              theta_hat=as.vector(theta_hat),ind_c=ind.c,
              se_beta=se_beta,lam1=lam1,lam2=lam2, Ve=Ve,intercept=G$intercept,sigma_2=sigma_2,VS=VS)
}
