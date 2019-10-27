#' Variable selection with Bivariate penalized Spline (GCV/CV) fit
#'
#' This is an internal function of package \code{ggam}. Bivariate penalized least squares problem
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
#' @param MI Whether model identification is conducted or not.
#' @param ... other arguments.
#' @details
#' This is an internal function of package \code{ggam}. We propose a coordinate descent based
#' algorithm to perform the variable selection efficiently.
#' The smoothing penalty parameter could
#' be chosen by \code{GCV} or \code{CV} using the routines: \code{\link{plsfitGCV}} and
#' \code{\link{plsfitCV}}. In this function, the user can also choose whether to do variable selection
#' or not.
#' @return A list of fit information.
#' @export
######################################################################
######################################################################
######################################################################
grplsfit <- function(G, criterion, method, family, ind_c, VS, MI,...)
{
  Terms <- delete.response(G$pterms)
  Xpp <- model.matrix(Terms,G$mf,contrasts=G$contrasts)
  if (G$m_b>0){
    for (i in 1:length(G$basis_info)){
      if (class(G$basis_info[[i]])[1]%in% 'bivariate.smooth'){
        Gb <- G$basis_info[[i]]
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
        Xp <- Xpp[Gb$ind,]
      }
    }
  } else {
    BQ2 <- NULL
    B <- NULL
    Q2 <- NULL
    P <- NULL
    Gb <-NULL
    Xp <- Xpp
    lambda <- NULL}
  y <- G$y
  X <- G$X
  if (MI){
    X <- cbind(G$Xp,X)
  }
  nvars <- NCOL(G$X)

  if (nvars == 0)
    stop("Model seems to contain no terms")
  if (NCOL(y) > 1)
    stop("y must be univariate unless binomial")

  if (G$intercept & !is.null(Gb)){
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
      for (j in G$ind.nl){
        if (class(G$basis_info[[j]])[1]=='univariate.smooth'){
          UB <- cbind(UB,G$basis_info[[j]]$B)
        }
      }
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
  if (!is.function(variance) || !is.function(linkinv))
    stop("illegal `family' argument")

  ### not sure want to include or not, deal with starting values, with like 'newton' method?

  # valideta <- family$valideta
  # if (is.null(valideta))
  #   valideta <- function(eta) TRUE
  # validmu <- family$validmu
  # if (is.null(validmu))
  #   validmu <- function(mu) TRUE


  ### initialize or not?
  # mu=initialize(y)
  # eta=linkfun(mu)
  mu <- y
  eta <- linkfun(y)
  # if (!(validmu(mu) && valideta(eta)))
  #   stop("Can't find valid starting values: please specify some")

  good <- weights > 0
  varmu <- variance(mu)[good]
  if (any(is.na(varmu)))
    stop("NAs in V(mu)")
  if (any(varmu == 0))
    stop("0s in V(mu)")
  mu.eta.val <- mu.eta(eta)
  if (any(is.na(mu.eta.val[good])))
    stop("NAs in d(mu)/d(eta)")

  good <- (weights > 0) & (mu.eta.val != 0)
  if (all(!good)) {
    warning(gettextf("No observations informative."))
  }
  mevg <- mu.eta.val[good]
  mug <- mu[good]
  yg <- y[good]
  weg <- weights[good]
  var.mug <- variance(mug)

  y <- z <- (eta - offset)[good] + (yg - mug)/mevg
  w <- sqrt((weg * mevg^2)/var.mug)
  if (!is.null(dim(X))){
    X <- X[good,,drop = FALSE]
    np <- ncol(X)
    X <- if(VS){apply(X,2,scale)} else {X}
  } else {
    X <- X[good,drop = FALSE]
    np <- ncol(as.matrix(X))
    X <- if(VS) {apply(as.matrix(X),2,scale)} else {X}
  }

  if (sum(!is.finite(G$y)) + sum(!is.finite(G$w)) > 0)
    stop("iterative weights or data non-finite in grplsfit")

  ### main algorithm for grplsfit starts here

  Y <- y
  n <- nobs <- length(Y)

  Xp <- cbind(Xp,G$linear_cov)
  colnames(Xp)<-c(colnames(Xpp),G$linear.names)
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
      mfit <- plsfitGCV(B,Q2,P,lambda,Y,fx=G$fx,Z=X)
      lam1 <- mfit$lam0
      ###
      alpha_hat <- mfit$alpha_hat
      beta_hat <- mfit$beta_hat
      gamma_hat <- mfit$gamma_hat
      theta_hat <- mfit$theta_hat
      sse <- mfit$sse
      res <- mfit$res
      gcv <- mfit$gcv
      tt <- (y-mfit$Yhat)^2/variance(mfit$Yhat)
      #print(variance(mfit$Yhat))
      sigma_2 <- 1/(length(y)-mfit$edf)*sum(tt)
      cv <- NULL
      bic <- mfit$bic
      Yhat <- mfit$Yhat
      edf <- mfit$edf
      edfs <- mfit$edfs
      if (length(Xp)>0){
        if (length(mfit$alpha_hat)>0){
          sb <- beta_se(UB,B,Q2,K,lam1,lam2,X,Y,ind.c,...)
          se_beta <- sqrt(sigma_2)*sb$se_beta
          Ve <- sqrt(sigma_2)*sigma_2*sb$Ve
        } else {
          se_beta <- NULL
          Ve <- NULL
        }
      } else {
        se_beta <- NULL
        Ve <- NULL
      }
    }

    ### variable selection
    if(VS==TRUE){
      if (!is.null(Gb)){
        lfit <- plsfitGCV(B,Q2,P,lambda,Y,fx=G$fx)
        lam1 <- lfit$lam0
        mfit <- lfit
        BQ2 <- B%*%Q2
        if (G$fx==FALSE){
          BQ2 <- as.matrix(BQ2)
          B <- as.matrix(B)
          HB <- BQ2%*%solve(t(Q2)%*%(t(B)%*%B+lam1*K)%*%Q2)%*%t(BQ2)
        } else {
          HB <- BQ2%*%solve(t(Q2)%*%(t(B)%*%B)%*%Q2)%*%t(BQ2)
        }
        Zstar <- as.matrix(X-HB%*%X)
        Ystar <- Y-lfit$beta_hat
      } else {
        lam1 <- NULL
        Zstar <- as.matrix(X)
        Ystar <- Y
      }
      ## check whether ind_c is given or not
      if (is.null(ind_c)){
        if(method=="SCAD"){
          mfit <- grpreg(Zstar,Ystar,group=1:np,penalty="grSCAD",family=G$family[[1]])
        }
        if(method=="ALASSO"){
          mfit <- grpreg(Zstar,Ystar,group=1:np,penalty="grLasso",family=G$family[[1]])
          temp <- as.matrix(mfit$beta[-1,])
          df.x <- apply(temp!=0,2,sum)
          BIC <- n*log(mfit$loss/n)+df.x*log(n)
          BIC[is.infinite(BIC)] <- 10^6
          alpha_hat <- mfit$beta[,which.min(BIC)]
          alpha.c <- alpha_hat[1]
          alpha <- alpha_hat[-1]
          temp <- 1/abs(alpha^3)
          weight <- ifelse(temp==Inf,10^3,temp)
          mfit <- grpreg(Zstar,Ystar,group=1:np,penalty="grLasso",family=G$family[[1]],
                      group.multiplier=weight)
        }
        temp2 <- as.matrix(mfit$beta[-1,])
        df.x <- apply(temp2!=0,2,sum)
        BIC <- n*log(mfit$loss/n)+df.x*log(n)
        BIC[is.infinite(BIC)] <- 10^6
        lam2 <- mfit$lambda[which.min(BIC)]
        alpha <- mfit$beta[,which.min(BIC)]
        alpha.c <- alpha[1]
        alpha <- alpha[-1]
        ind.c <- (1:np)[alpha!=0]
      } else {
        ##
        lam2 <- NULL
        ind.c <- ind_c
        alpha.c <- 0
      }
      Zc <- as.matrix(X[,ind.c])

      ##
      if (is.null(Gb)){
        # need to be changed according to how to do variable selection when UB is available.
        #lam1 <- 0
        alpha_hat <- alpha
        Yhat <- X%*%alpha+alpha.c
        res <- Y-Yhat
        sse <- sum((Y-Yhat)^2)
        tt <- (y-Yhat)^2/variance(Yhat)
        sigma_2 <- mfit$sigma_2
        gcv <- NULL
        cv <- NULL
        gamma_hat <- NULL
        theta_hat <- NULL
        bic <- BIC
        beta_hat <- 0
        ## need to change: se in grpreg
        ## if (se=='TRUE)
        se_beta <- rep(0,length(ind.c))
        Ve <- matrix(0,ncol=length(ind.c),nrow=length(ind.c))
        edf <- NULL
        #edfs <- NULL
      } else {
        # need to be changed according to how to do variable selection when UB is available.
        mfit <- plsfitGCV(B,Q2,P,G$lambda,Y,fx=G$fx,Z=Zc,...)
        lam1 <- mfit$lam0
        alpha_hat <- rep(0,np)
        alpha_hat[ind.c] <- mfit$alpha_hat
        beta_hat <- mfit$beta_hat
        gamma_hat <- mfit$gamma_hat
        theta_hat <- mfit$theta_hat
        sse <- mfit$sse
        res <- mfit$res
        gcv <- mfit$gcv
        cv <- NULL
        bic <- mfit$bic
        Yhat <- mfit$Yhat
        tt <- (y-mfit$Yhat)^2/variance(mfit$Yhat)
        sigma_2 <- 1/(length(y)-mfit$edf)*sum(tt)
        #df for s(u,v)
        edf <- mfit$edf
        #edfs <- mfit$edfs
        #print(length(mfit$alpha_hat))
        if (length(mfit$alpha_hat)>0){
          sb <- beta_se(UB,B,Q2,K,lam1,lam2,Zc,Y,VS=VS)
          se_beta <- sb$se_beta
          Ve <- sb$Ve
        } else {
          se_beta <- NULL
          Ve <- NULL
         }
      }
    }
  }# end of GCV

  # start of CV, need to test!
  if (criterion=='CV'){
    # No variable selection
    if(VS==FALSE){
      ### no selection of variables so no ind.c
      ind.c <- 1:length(which(!G$term.names %in% '(Intercept)'))
      lam2 <- 0
      if (!is.null(Gb)){
        if (length(X)>0){mfit <- plsfitCV(B,Q2,K,lambda,Y,Z=X)
        } else {mfit <- plsfitCV(B,Q2,K,lambda,Y,fx=G$fx)}
        lam1 <- mfit$lam0
        ###
        alpha_hat <- mfit$alpha_hat
        #alpha.c <- 0
        beta_hat <- mfit$beta_hat
        gamma_hat <- mfit$gamma_hat
        theta_hat <- mfit$theta_hat
        sse <- mfit$sse
        res <- mfit$res
        cv <- mfit$cv
        gcv <- NULL
        bic <- mfit$bic
        Yhat <- mfit$Yhat
        edf <- mfit$edf
        #edfs <- mfit$edfs
        if (length(X)>0){
          if (length(mfit$alpha_hat)>0){
            sb <- beta_se(B,Q2,K,lam1,lam2,X,Y)
            se_beta <- sb$se_beta
            Ve <- sb$Ve
          } else {
            se_beta <- NULL
            Ve <- NULL
          }
        } else {
          se_beta <- NULL
          Ve <- NULL
        }
      }
      else {
        # no bivariate term, no penalty on that
        lam1  <-  NULL
        # here X is without intercept
        # if (G$intercept){
          #mfit <- glm(Y~X,family=G$family[[1]])
          mfit <- plsfitGCV(B,Q2,P,lambda,Y,G$fx,X,...)
          alpha.c <- mfit$alpha_hat[1]
          alpha_hat <- mfit$alpha_hat[-1]
          # } else {
          #   mfit <- glm(Y~X-1,family=G$family[[1]])
          #   alpha.c <- NULL
          #   alpha_hat <- coef(mfit)
          # }
        Yhat <- mfit$fitted.values
        res <- Y-Yhat
        sse <- sum(res^2)
        cv <- NULL
        gcv <- NULL
        bic <- NULL
        theta_hat <- NULL
        gamma_hat <- NULL
        edf <- NULL
        beta_hat <- 0
        #need to change (se in glm)
        # if (length(X)>0){
        #     se_beta <- summary(mfit)$coefficients[, 2]
        #     Ve <- summary(mfit)$cov.unscaled
        # }else {
          se_beta <- NULL
          Ve <- NULL
          #}
      }
    }

    ### variable selection
    if(VS==TRUE){
      if (!is.null(Gb)){
        lfit <- plsfitCV(B,Q2,K,lambda,Y,fx=G$fx)
        lam1 <- lfit$lam0
        mfit <- lfit
        BQ2 <- B%*%Q2
        if (G$fx==FALSE){
          HB <- BQ2%*%solve(t(Q2)%*%(t(B)%*%B+lam1*K)%*%Q2)%*%t(BQ2)
        } else {
          HB <- BQ2%*%solve(t(Q2)%*%(t(B)%*%B)%*%Q2)%*%t(BQ2)
        }
        Zstar <- as.matrix(X-HB%*%X)
        Ystar <- Y-lfit$beta_hat
      } else {
        lam1 <- NULL
        Zstar <- as.matrix(X)
        Ystar <- Y
      }
      ## check whether ind_c is given or not
      if (is.null(ind_c)){
        if(method=="SCAD"){
          mfit <- grpreg(Zstar,Ystar,group=1:np,penalty="grSCAD",family=G$family[[1]])
        }
        if(method=="ALASSO"){
          mfit <- grpreg(Zstar,Ystar,group=1:np,penalty="grLasso",family=G$family[[1]])
          df.x <- apply(mfit$beta[-1,]!=0,2,sum)
          BIC <- n*log(mfit$loss/n)+df.x*log(n)
          BIC[is.infinite(BIC)] <- 10^6
          alpha_hat <- mfit$beta[,which.min(BIC)]
          alpha.c <- alpha_hat[1]
          alpha <- alpha_hat[-1]
          temp <- 1/abs(alpha^3)
          weight <- ifelse(temp==Inf,10^3,temp)
          mfit <- grpreg(Zstar,Ystar,group=1:np,penalty="grLasso",family=G$family[[1]],
                      group.multiplier=weight)
        }
        df.x <- apply(mfit$beta[-1,]!=0,2,sum)
        BIC <- n*log(mfit$loss/n)+df.x*log(n)
        BIC[is.infinite(BIC)] <- 10^6
        lam2 <- mfit$lambda[which.min(BIC)]
        alpha <- mfit$beta[,which.min(BIC)]
        alpha.c <- alpha[1]
        alpha <- alpha[-1]
        ind.c <- (1:np)[alpha!=0]
      } else {
        ##
        lam2 <- NULL
        ind.c <- ind_c
        alpha.c <- 0
      }
      Zc <- as.matrix(X[,ind.c])
      ##
      if (is.null(Gb)){
        #lam1 <- 0
        alpha_hat <- alpha
        Yhat <- X%*%alpha+alpha.c
        res <- Y-Yhat
        tt <- (y-Yhat)^2/variance(Yhat)
        sigma_2 <- 1/(length(y)-mfit$edf)*sum(tt)
        sse <- sum((Y-Yhat)^2)
        cv <- NULL
        gcv <- NULL
        gamma_hat <- NULL
        theta_hat <- NULL
        bic <- BIC
        beta_hat <- 0
        ## need to change: se in grpreg
        ## if (se=='TRUE)
        se_beta <- rep(0,length(ind.c))
        Ve <- matrix(0,ncol=length(ind.c),nrow=length(ind.c))
        edf <- NULL
        #edfs <- NULL
      } else {
        mfit <- plsfitCV(B,Q2,K,G$lambda,Y,G$fx,Z=Zc)
        lam1 <- mfit$lam0
        alpha_hat <- rep(0,np)
        alpha_hat[ind.c] <- mfit$alpha_hat
        beta_hat <- mfit$beta_hat
        gamma_hat <- mfit$gamma_hat
        theta_hat <- mfit$theta_hat
        sse <- mfit$sse
        res <- mfit$res
        cv <- mfit$cv
        gcv <- NULL
        bic <- mfit$bic
        Yhat <- mfit$Yhat
        tt <- (y-mfit$Yhat)^2/variance(mfit$Yhat)
        sigma_2 <- 1/(length(y)-mfit$edf)*sum(tt)
        #df for s(u,v)
        edf <- mfit$edf
        #edfs <- mfit$edfs
        if (length(mfit$alpha_hat)>0){
          sb <- beta_se(B,Q2,K,lam1,lam2,Zc,Y,ind.c,...)
          se_beta <- sb$se_beta
          Ve <- sb$Ve
        } else {
          se_beta <- NULL
          Ve <- NULL
        }
      }
    }
  } # end of CV

  # get the final estimation here
  eta <- mfit$Yhat
  if (any(!is.finite(alpha_hat))) {
    warning(gettextf("Non-finite coefficients"))}

  mu <- linkinv(eta <- eta + offset)
  eta <- linkfun(mu)
  dev <- sum(dev.resids(Y, mu, weights))

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
  residuals <- rep(NA, nobs)
  residuals[good] <- z - (eta - offset)[good]

  ### For those points not in the triangles, return "NA"
  G$y <- residuals_all <- weights_all <- eta_all <- fitted.values_all <- Yhat_all <- fitted.values_all_sbl <- rep(NA,dim(G$mf)[1])
  if (!is.null(Gb$ind)){
    G$y[Gb$ind]<-y
    fitted.values_all[Gb$ind]<-mu
    fitted.values_all_sbl3 <- fitted.values_all_sbl[Gb$ind]
    fitted.values_all_sbl3[Ind]<-y_sbl
    fitted.values_all_sbl[Gb$ind]<-fitted.values_all_sbl3
    #fitted.values_all_sbl[Gb$ind]<-y_sbl
    residuals_all[Gb$ind] <- residuals
    eta_all[Gb$ind] <- eta
    weights_all[Gb$ind] <- if (is.null(weights)){rep(1,length(y))} else {weights}
  } else {
    G$y<-y
    fitted.values_all<-mu
    fitted.values_all_sbl3 <- fitted.values_all_sbl
    fitted.values_all_sbl3[Ind]<-y_sbl
    fitted.values_all_sbl<-fitted.values_all_sbl3
    residuals_all <- residuals
    eta_all <- eta
    weights_all <- if (is.null(weights)){rep(1,length(y))} else {weights}
  }

  wtdmu <- if (G$intercept)
    sum(weights * Y)/sum(weights) else linkinv(offset)
  nulldev <- sum(dev.resids(Y, wtdmu, weights))
  n.ok <- nobs - sum(weights == 0)
  nulldf <- n.ok - as.integer(G$intercept)
  rank2<-ifelse(is.null(Gb),0,dim(Gb$Q2)[1])
  # aic.model <- aic(Y, n, mu, weights, dev) + 2 * sum(G$edf)

  object <- list(coefficients = as.vector(alpha_hat), coefficients_bivariate = as.vector(gamma_hat),
              residuals = residuals_all,dev_sbl=dev_sbl,
              fitted.values = fitted.values_all, family = family, linear.predictors = eta_all,
              deviance = dev, null.deviance = nulldev, edf=edf, prior.weights = weights_all,#edfs=edfs,
              df.null = nulldf,backfitting=G$backfitting,Xp=X3,
              y = G$y, sigma_2=sigma_2,#sig2 = G$sig2,
              r = G$r, nsdf = G$nsdf, Ve = Ve, nsdf=G$nsdf,
              gcv_opt= gcv, cv_opt=cv, # aic = aic.model,
              mhat.sbl=mhat.sbl,mhat=mhat,X2=X2,h_opt_all=h_opt_all,
              fitted.values.sbl=fitted.values_all_sbl,
              theta_hat=as.vector(theta_hat),ind_c=ind.c,
              se_beta=se_beta,lam1=lam1,lam2=lam2, Ve=Ve,intercept=G$intercept,VS=VS)
}
