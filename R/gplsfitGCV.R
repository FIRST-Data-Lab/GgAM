#' Generalized Penalized Least Square Fit under GCV
#'
#' This is an internal function of package \code{ggam}.
#'
#' @importFrom magic adiag
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
#' @param control A list of fit control parameters to replace defaults returned by \code{\link{plbpsm.control}}.
#' Any control parameters not supplied stay at their default values.
#' spline (\code{TRUE}) or a penalized regression spline (\code{FALSE}).
#' @param X The parametric model matrix. set to '\code{NULL}' if it is not provided.
#' @param ind_c The vector of index to indicate the parametric part.
#' @param start Initial values for model coefficients
#' @param etastart Initial values for linear predictor.
#' @param mustart Initial values for the expected response.
#' @param theta The given theta values in negative binomial family.
#' @param fixedSteps How many steps to take: useful when only using this routine to get rough starting
#' values for other methods.
#' @param ... other arguments.

#' @details
#' See section 4 'Implementation' in Shan et al. (2018).
#' @return
#' A list of fit information.
#' @examples
#' library(BPST)
#' data("eg_poi")
#' eg1_V1 <- eg_poi[['V1']]
#' eg1_T1 <- eg_poi[['T1']]
#' sam <- eg_poi[['sam_poi']]
#' d <- 2
#' r <- 1
#' B0 <- basis(eg1_V1,eg1_T1, d, r, sam[,c('loc1','loc2')])
#' B <- B0$Bi
#' ind <- B0$Ind.inside
#' Q2 <- B0$Q2
#' K <- B0$K
#' Z <- sam[ind,c(5:12)]
#' Y <- sam[ind,'y']
#' lambda_start <- 0.01
#' lambda_end <- 32
#' nlambda <- 10
#' lambda <- exp(seq(log(lambda_start),log(lambda_end),length.out=nlambda))
#' X <- sam[,1:15]
#' P <- t(Q2)%*%K%*%Q2
#' gplsfitGCV(Y,as.matrix(B),Q2,P,UB=NULL,lambda=lambda,family=poisson(),offset=0,fx=FALSE,
#' control = plbpsm.control(),X=as.matrix(X))
#' @export
gplsfitGCV=function(Y, B, Q2, P, UB = NULL, lambda, family, offset, theta = 0, fx, control,
                         start = NULL, etastart = NULL, mustart = NULL, X = NULL, ind_c = 1:ncol(X),
                         fixedSteps = (control$maxstep + 1),...){
  conv <- FALSE
  y <- Y
  Xp<- X[,ind_c]
  if (!is.null(B)){
    BQ2 <- as.matrix(B%*%Q2)
  } else {BQ2 <- NULL}
  ### Define some functions to use in different families
  variance <- family$variance
  dev.resids <- family$dev.resids
  aic <- family$aic
  linkinv <- family$linkinv
  linkfun <- family$linkfun
  mu.eta <- family$mu.eta
  if (!is.function(variance) || !is.function(linkinv))
    stop("illegal `family' argument")
  valideta <- family$valideta
  if (is.null(valideta))
    valideta <- function(eta) TRUE
  validmu <- family$validmu
  if (is.null(validmu))
    validmu <- function(mu) TRUE

  ### set up
  nobs<-length(Y)
  weights <- rep(1,nobs)
  nl <- length(lambda)
  nobs <- length(Y)
  if (!is.null(BQ2)){
    if(!is.null(X)){
      d1 <- dim(X)[2]
      d2 <- dim(BQ2)[2]
      zero1 <- matrix(rep(0,d1*d2),ncol = d2)
      zero2 <- matrix(rep(0,d1*(d1+d2)),ncol=d1)
      D <- rbind(zero1,P)
      D <- cbind(zero2,D)
      Z <- cbind(X,BQ2)
    } else{
      d1 <- 0
      d2 <- dim(BQ2)[2]
      D <- P
      Z <- BQ2
      #alpha_all <- matrix(rep(0,d2*nl),ncol=nl)
    }
  } else {
    if(!is.null(X)){
      d1 <- dim(X)[2]
      d2 <- 0
      # zero1 <- matrix(rep(0,d1*d2),ncol = d2)
      # zero2 <- matrix(rep(0,d1*(d1+d2)),ncol=d1)
      # D <- rbind(zero1,P)
      # D <- cbind(zero2,D)
      Z <- cbind(X,BQ2)
    } else{
      d1 <- 0
      d2 <- 0
      Z <- NULL
    }
  }

  nvars <- ncol(Z)

  ### From Wood
  ### mustart given or not
  if (is.null(mustart))   # new from version 1.5.0
  { eval(family$initialize)} else { mukeep <- mustart
  eval(family$initialize)
  mustart <- mukeep
  }

  ### etastart given or not
  coefold <- NULL                 # 1.5.0
  etastart2 <- if (!is.null(etastart)) {
    etastart} else if (!is.null(start)){
      if (length(as.matrix(start)) != nvars){
        stop(gettextf("Length of start should equal %d and correspond to initial coefs.",
                      nvars))} else {
                        coefold <- start
                        offset + as.vector(if (NCOL(X) == 1){
                          X * start
                        } else {X %*% start})
                      }
    } else {
      family$linkfun(mustart)}

  ### recalculate mustart
  mustart2 <- linkinv(etastart2)
  if (!(validmu(mustart2) && valideta(etastart2)))
    stop("Can't find valid starting values: please specify some")
  devold <- sum(dev.resids(Y, mustart2, weights))
  boundary <- FALSE

  ## fx = TRUE, no penalization on bivariate spline

  ## fx = FALSE
  #if (fx == FALSE){
    ###
    if (!is.null(lambda)){
    W_iter_all <- c()
    se_beta_2_all <- c()
    se_beta_1_all <- c()
    alpha_all <- c()
    var_beta_1_all <- list()
    var_beta_2_all <- list()
    gcv_all <- c()
      for(il in 1:nl){
        mu <- mustart2
        eta <- etastart2
        for (iter in 1:fixedSteps) {
          # while(delta1>epsilon & sum(is.infinite(delta2))==0 & step<=maxstep ){
          #   step <- step+1
          good <- weights > 0
          if (theta!=0){
            variance <- function(mu,theta) mu+mu^2/theta
            varmu <- variance(mu,theta)[good]
          } else {
            varmu <- variance(mu)[good]
          }
          if (any(is.na(varmu)))
            stop("NAs in V(mu)")
          if (any(varmu == 0))
            stop("0s in V(mu)")
          mu.eta.val <- mu.eta(eta)
          if (any(is.na(mu.eta.val[good])))
            stop("NAs in d(mu)/d(eta)")
          good <- (weights > 0) & (mu.eta.val != 0) # note good modified here => must re-calc each iter
          if (all(!good)) {
            conv <- FALSE
            warning(gettextf("No observations informative at iteration %d",
                             iter))
            break
          }
          mevg <- mu.eta.val[good]
          mug <- mu[good]
          yg <- Y[good]
          weg <- weights[good]
          #var.mug <- variance(mug)

          ### need to be changed according the Wood.
          if(theta!=0) {
            #variance=function(mu,theta) mu+mu^2/theta
            var <- variance(mu,theta)
          }else{
            var <- variance(mu)
          }
          #mevg <- mu.eta(eta)
          Y_iter <- z <- (eta - offset)[good] + (yg - mug)/mevg
          W_iter <- weg*(mevg^2)/var
          W_iter <- as.vector(W_iter*weg)
          temp1 <- as.matrix(W_iter*Z)
          temp2 <- W_iter*Y_iter
          Z <- as.matrix(Z)
          temp3 <- crossprod(Z,temp1)+lambda[il]*D
          alpha_old <- solve(temp3,crossprod(Z,temp2))
          if (any(!is.finite(alpha_old))) {
            conv <- FALSE
            warning(gettextf("Non-finite coefficients at iteration %d",iter))
            break
          }
          #print(alpha_old[1:3])
          start <- alpha_old
          eta <- Z%*%alpha_old+offset
          eta <- as.matrix(eta)
          mu <- linkinv(eta)
          dev <- sum(dev.resids(Y, mu, weights))
          #print(dev)
          if (!is.finite(dev)) {
            if (is.null(coefold))
              stop("no valid set of coefficients has been found:please supply starting values",
                   call. = FALSE)
            warning("Step size truncated due to divergence",call.=FALSE)
            ii <- 1
            while (!is.finite(dev)) {
              if (ii > control$maxstep)
                stop("inner loop 1; can't correct step size")
              ii <- ii + 1
              start <- (start + coefold)/2
              eta<-drop(X %*% start)
              mu <- linkinv(eta <- eta + offset)
              #eta <- linkfun(mu)
              dev <- sum(dev.resids(Y, mu, weights))
            }
            boundary <- TRUE
          }
          if (!(valideta(eta) && validmu(mu))) {
            warning("Step size truncated: out of bounds.",call.=FALSE)
            ii <- 1
            while (!(valideta(eta) && validmu(mu))) {
              if (ii > control$maxstep)
                stop("inner loop 2; can't correct step size")
              ii <- ii + 1
              start <- (start + coefold)/2
              eta<-drop(X %*% start)
              mu <- linkinv(eta <- eta + offset)
              eta<-linkfun(mu)
            }
            boundary <- TRUE
            dev <- sum(dev.resids(Y, mu, weights))
            if (control$trace)
              cat("Step halved: new deviance =", dev, "\n")
          }## end of check
          if (abs(dev - devold)/(0.1 + abs(dev)) < control$epsilon || #olm ||
              iter >= fixedSteps) {
            conv <- TRUE
            coef <- start #1.5.0
            break
          }else {
            devold <- dev
            coefold <- coef<-start
          }
        }
        alpha_all <- cbind(alpha_all,alpha_old)
        W_iter_all <- c(W_iter_all,W_iter)
        #calculate gcv
        DD <- solve(temp3)
        #Slambda <- Z%*%tcrossprod(DD,Z)
        Slambda <- Z%*%tcrossprod(DD,temp1)
        df <- sum(diag(Slambda))
        yhat <- Z%*%alpha_old

        # calculate gcv
        gcv <- nobs*sum(W_iter*(Y_iter-yhat)^2)/(nobs-df)^2;
        gcv_all <- c(gcv_all,gcv)
        #calculate standard deviation
        if(d1!=0) {
          beta_hat <- alpha_old[1:d1]
          theta_hat <- alpha_old[(d1+1):(d1+d2)]
          gamma_hat <- Q2%*%theta_hat
          etahat <- B%*%as.vector(gamma_hat)+as.matrix(X)%*%beta_hat
        } else {
          theta_hat <- alpha_old
          gamma_hat <- Q2%*%theta_hat
          beta_hat <- NULL
          etahat <- B%*%as.vector(gamma_hat)
        }
        etahat <- as.matrix(etahat)
        YHAT <- linkinv(etahat)
        if(theta!=0){
          tt <- (y-YHAT)^2/variance(YHAT,theta)
        } else {tt <- (y-YHAT)^2/variance(YHAT)}
        sigma2 <- 1/(length(y)-df)*sum(tt)

        BQ2 <- as.matrix(BQ2)
        BB <- cbind(UB,BQ2)
        W_iter <- W_iter/sigma2
        # print(length(Xp))
        # print(head(B))
        # print(111)
        if (!is.null(B) && !is.null(Xp)){
          V22 <- crossprod(BB,W_iter*BB)
          v22_inv <- if (!is.null(UB)) {
            solve(V22+adiag(matrix(0,ncol=ncol(UB),nrow=ncol(UB)),
                                                       as.matrix((lambda[il]*P))))
            } else {solve(V22+as.matrix(lambda[il]*P))}
          #v22_inv <- solve(V22+as.matrix(lambda[il]*P))
          if (length(Xp)>0){
            temp3 <- W_iter*Xp
            Xphat <- as.matrix(BB)%*%v22_inv%*%crossprod(as.matrix(BB),as.matrix(temp3))
            temp4 <- W_iter*(Xp-Xphat)
            var_beta_1 <- var_beta_1_all[[il]] <- solve(crossprod(as.matrix(Xp-Xphat),as.matrix(temp4))) #version 1
            #var_beta_2 <- var_beta_2_all[[il]] <- solve(crossprod(Xp,temp3)) #version 2
            se_beta_1_all <- cbind(se_beta_1_all,sqrt(diag(var_beta_1)))
            #se_beta_2_all <- cbind(se_beta_2_all,sqrt(diag(var_beta_2)))
            # print(sqrt(diag(var_beta_1)))
            # print(222)
          }
        } else if
        (!is.null(Xp) && is.null(UB)){
            var_beta_1 <- t(Xp) %*% W_iter%*%Xp
            se_beta_1 <- sqrt(diag(var_beta_1))
            }else{
          var_beta_1 <-se_beta_1 <- NULL
        }
      }# end of 1:nl
      #plot(log(lambda),gcv_all,type='l')
      j <- which.min(gcv_all)
      W_iter <- W_iter_all[j]
      gcv <- gcv_all[j]
      lambdac <- lambda[j]
      alpha_hat <- alpha_all[,j]
    } #end of with BPST

    else {
      mu <- mustart2
      eta <- etastart2
      for (iter in 1:fixedSteps) {
        good <- weights > 0
        if (theta!=0){
          variance <- function(mu,theta) mu+mu^2/theta
          varmu <- variance(mu,theta)[good]
        } else {
          varmu <- variance(mu)[good]
        }
        if (any(is.na(varmu)))
          stop("NAs in V(mu)")
        if (any(varmu == 0))
          stop("0s in V(mu)")
        mu.eta.val <- mu.eta(eta)
        if (any(is.na(mu.eta.val[good])))
          stop("NAs in d(mu)/d(eta)")
        good <- (weights > 0) & (mu.eta.val != 0) # note good modified here => must re-calc each iter
        if (all(!good)) {
          conv <- FALSE
          warning(gettextf("No observations informative at iteration %d",
                           iter))
          break
        }
        mevg <- mu.eta.val[good]
        mug <- mu[good]
        yg <- Y[good]
        weg <- weights[good]
        #var.mug <- variance(mug)

        ### need to be changed according the Wood.
        if(theta!=0) {
          #variance=function(mu,theta) mu+mu^2/theta
          var <- variance(mu,theta)
        }else{
          var <- variance(mu)
        }
        #mevg <- mu.eta(eta)
        Y_iter <- z <- (eta - offset)[good] + (yg - mug)/mevg
        W_iter <- weg*(mevg^2)/var
        W_iter <- as.vector(W_iter*weg)
        temp1 <- as.matrix(W_iter*Z)
        temp2 <- W_iter*Y_iter
        Z <- as.matrix(Z)
        temp3 <- crossprod(Z,temp1)
        alpha_old <- solve(temp3,crossprod(Z,temp2))
        if (any(!is.finite(alpha_old))) {
          conv <- FALSE
          warning(gettextf("Non-finite coefficients at iteration %d",iter))
          break
        }
        #print(alpha_old[1:3])
        start <- alpha_old
        eta <- Z%*%alpha_old+offset
        eta <- as.matrix(eta)
        mu <- linkinv(eta)
        dev <- sum(dev.resids(Y, mu, weights))
        #print(dev)
        if (!is.finite(dev)) {
          if (is.null(coefold))
            stop("no valid set of coefficients has been found:please supply starting values",
                 call. = FALSE)
          warning("Step size truncated due to divergence",call.=FALSE)
          ii <- 1
          while (!is.finite(dev)) {
            if (ii > control$maxstep)
              stop("inner loop 1; can't correct step size")
            ii <- ii + 1
            start <- (start + coefold)/2
            eta<-drop(X %*% start)
            mu <- linkinv(eta <- eta + offset)
            #eta <- linkfun(mu)
            dev <- sum(dev.resids(Y, mu, weights))
          }
          boundary <- TRUE
        }
        if (!(valideta(eta) && validmu(mu))) {
          warning("Step size truncated: out of bounds.",call.=FALSE)
          ii <- 1
          while (!(valideta(eta) && validmu(mu))) {
            if (ii > control$maxstep)
              stop("inner loop 2; can't correct step size")
            ii <- ii + 1
            start <- (start + coefold)/2
            eta<-drop(X %*% start)
            mu <- linkinv(eta <- eta + offset)
            eta<-linkfun(mu)
          }
          boundary <- TRUE
          dev <- sum(dev.resids(Y, mu, weights))
          if (control$trace)
            cat("Step halved: new deviance =", dev, "\n")
        }## end of check
        if (abs(dev - devold)/(0.1 + abs(dev)) < control$epsilon || #olm ||
            iter >= fixedSteps) {
          conv <- TRUE
          coef <- start #1.5.0
          break
        }else {
          devold <- dev
          coefold <- coef<-start
        }
      }
      #calculate gcv
      DD <- solve(temp3)
      #Slambda <- Z%*%tcrossprod(DD,Z)
      Slambda <- Z%*%tcrossprod(DD,temp1)
      df <- sum(diag(Slambda))
      yhat <- Z%*%alpha_old
      alpha_hat <- alpha_old
      lambdac <- NULL

      # calculate gcv
      gcv <- nobs*sum(W_iter*(Y_iter-yhat)^2)/(nobs-df)^2;

      #calculate standard deviation
      if (!is.null(UB) && !is.null(Xp)){
        V22 <- crossprod(UB,W_iter*UB)
        v22_inv <- solve(V22)
        if (length(Xp)>0){
          temp3 <- W_iter*Xp
          Xphat <- as.matrix(UB)%*%v22_inv%*%crossprod(as.matrix(UB),as.matrix(temp3))
          temp4 <- W_iter*(Xp-Xphat)
          var_beta_1 <- solve(crossprod(as.matrix(Xp-Xphat),as.matrix(temp4))) #version 1
          se_beta_1 <- sqrt(diag(var_beta_1))
          #print(se_beta_1)
        }
      } else if (!is.null(Xp)){
        var_beta_1 <- t(Xp) %*% W_iter%*%Xp
        se_beta_1 <- sqrt(diag(var_beta_1))
      } else {
        var_beta_1 <-se_beta_1 <- NULL}
      }
  #}#end of fx == FALSE

  if(d1!=0 && d2!=0) {
    beta_hat <- alpha_hat[1:d1]
    theta_hat <- alpha_hat[(d1+1):(d1+d2)]
    gamma_hat <- Q2%*%theta_hat
    etahat <- B%*%as.vector(gamma_hat)+as.matrix(X)%*%beta_hat
    if (length(Xp)>0){
      #print(var_beta_1_all)
          var_beta_1 <- var_beta_1_all[[j]]
    se_beta_1 <- se_beta_1_all[,j]
    } else {
      se_beta_1 <- var_beta_1 <- NULL
    }
    beta_hat2 <- B%*%as.vector(gamma_hat)
  } else if (d1!=0 && d2==0) {
    beta_hat <- alpha_hat
    theta_hat <- NULL
    gamma_hat <- NULL
    etahat <- as.matrix(X)%*%beta_hat
    # if (length(Xp)>0){
    #   se_beta_1 <- var_beta_1 <- NULL
    # }
    beta_hat2 <- rep(0,length(y))
    } else if (d1==0 && d2!=0){
    theta_hat <- alpha_hat
    gamma_hat <- Q2%*%theta_hat
    beta_hat <- NULL
    beta_hat2 <- B%*%as.vector(gamma_hat)
    etahat <- B%*%as.vector(gamma_hat)
  }
  etahat <- as.matrix(etahat)
  Yhat <- linkinv(etahat)

  ### need to be changed
  res <- Y_iter-yhat
  #sse <- sum((Y-Yhat)^2)
  if (!conv)
  { warning("Algorithm did not converge")
  }

  list(boundary = boundary,est_theta = theta,alpha_hat = beta_hat,theta_hat = theta_hat,gamma_hat = gamma_hat,beta_hat = beta_hat2,
       lam0 = lambdac,gcv = gcv,df = df,eta_hat = etahat,Yhat = Yhat,yiter = Y_iter,res = res,w = W_iter,se_beta = se_beta_1, z = z,Ve = var_beta_1,conv = conv,hat = Slambda)
}
