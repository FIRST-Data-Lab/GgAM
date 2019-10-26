#' Summary for a PLBPSM fit
#'
#' Takes a fitted \code{plbpsm} object produced by \code{plbpsm()} and produces various useful summaries from it.
#' @importFrom stats pnorm pt pchisq pf var
#' @param object a fitted \code{plbpsm} object as produced by \code{plbpsm()}.
#' @param x a \code{summary.plbpsm} object produced by \code{summary.plbpsm()}.
#' @param dispersion a value for the dispersion parameter: not normally used.
#'@param digits controls number of digits printed in output.
#'@param signif.stars Should significance stars be printed alongside output.
#'@param h_opt the bandwidth given for Spline-backfitting local estimator, default is \code{NULL}.
#'@param X0 the new predict matrix for obtaining simultaneous confidence band.
#' @param ... other arguments.
#' @method summary plbpsm
#' @return \code{summary.plbpsm} produces a list of summary information for a fitted \code{plbpsm} object.
#'\item{p.coeff}{is an array of estimates of the strictly parametric model coefficients.}
#'\item{p.t}{is an array of the \code{p.coeff}'s divided by their standard errors.}
#'  \item{p.pv}{is an array of p-values for the null hypothesis that the corresponding parameter is zero.
#'  Calculated with reference to the t distribution with the estimated residual
#'  degrees of freedom for the model fit if the dispersion parameter has been
#'  estimated, and the standard normal if not.}
#'  \item{m}{The number of smooth terms in the model.}
#'  \item{se}{array of standard error estimates for all parameter estimates.}
#'  \item{r.sq}{The adjusted r-squared for the model. Defined as the proportion of variance explained, where original variance and
#'  residual variance are both estimated using unbiased estimators. This quantity can be negative if your model is worse than a one
#'  parameter constant model, and can be higher for the smaller of two nested models! The proportion null deviance
#'  explained is probably more appropriate for non-normal errors. Note that \code{r.sq} does not include any offset in the one parameter model.}
#' \item{dev.expl}{The proportion of the null deviance explained by the model. The null deviance is computed taking account of any offset, so
#' \code{dev.expl} can be substantially lower than \code{r.sq} when an offset is present.}
#'  \item{edf}{array of estimated degrees of freedom for the model terms.}
#' \item{residual.df}{estimated residual degrees of freedom.}
#'  \item{n}{number of data.}
#'  \item{np}{number of model coefficients (regression coefficients, not smoothing parameters or other parameters of likelihood).}
#' \item{criterion}{The criterion to choose the penalty parameter lambda. \code{"GCV"} to use
#' generalized cross validation method and \code{"CV"} for cross validation}
#' \item{family}{The family object, specifying the distribution and link to use.}
#' \item{method}{'ALASSO' or 'SCAD' to penalize the coefficients for parametric part.}
#'  \item{formula}{the original PLBPSM formula.}
#'  \item{dispersion}{the scale parameter.}
#'  \item{pTerms.df}{the degrees of freedom associated with each parametric term
#'  (excluding the constant).}
#'  \item{pTerms.chi.sq}{a Wald statistic for testing the null hypothesis that the
#'  each parametric term is zero.}
#' \item{pTerms.pv}{p-values associated with the tests that each term is
#'  zero. For penalized fits these are approximate. The reference distribution
#'  is an appropriate chi-squared when the
#' scale parameter is known, and is based on an F when it is not.}
#'  \item{cov.scaled}{The estimated covariance matrix of the parameters.}
#'  \item{p.table}{significance table for parameters}
#' \item{p.Terms}{significance table for parametric model terms}
#' \item{gcv_opt}{The optimized gcv score.}
#' \item{cv_opt}{The optimized cv score.}
#' \item{bands}{A list of confidence bands for univaratie functions estimates.}
#' \item{mhat}{The estimated values for each linear or nonlinear term.}
#'@examples
#'  library(MASS)
#'  library(grpreg)
#'  library(BPST)
#'  data("eg1pop_dat")
#'  eg1_V2=eg1pop_dat[['V2']]
#'  eg1_T2=eg1pop_dat[['T2']]
#' eg1pop_rho03=eg1pop_dat[['rho03']]
#' sam=eg1pop_rho03[sample(1:dim(eg1pop_rho03)[1],100),]
#' lambda=10^(seq(-2,5,by=1))
#' data=sam
#' formula=Y~z1+z2+z3+z4+z5+z6+z7+z8+b(x1,x2,V=eg1_V2,Tr=eg1_T2,d=2,r=1,lambda=lambda)
#' # example 1
#' res=plbpsm(formula=formula,data=as.data.frame(data),VS=TRUE)
#' # example 12: ALASSO
#' res12=plbpsm(formula=formula,data=as.data.frame(data),VS=TRUE)
#' res10=plbpsm(formula=formula,data=as.data.frame(data),drop.intercept=TRUE)
#' # compare results under different settings
#' summary(res)
#' summary(res10)
#' summary(res12)
#'
#' ### GGAM-SMILE ###
#' data(eg1pop_poi2)
#' n=100
#' Npop=nrow(eg1pop_poi2)
#' ind.pop=(1:Npop)
#' sam.ind=sort(sample(ind.pop,n,replace=FALSE))
#' sam=eg1pop_poi2[sam.ind,]
#' data=sam
#' formula=Y~z1+u(z2)+u(z3)+b(x1,x2,V=eg1_V2,Tr=eg1_T2,d=2,r=1)
#' res_eg1_poi_add=plbpsm(formula=formula,data=as.data.frame(data),family='poisson')
#' summary(res_eg1_poi_add)
#' res_ggams=summary(res_eg1_poi_add)
#'
#' # The following is the SBL estimator for u(z2)
#' res_ggams$bands[[1]]$est
#'@export
summary.plbpsm <- function (object, h_opt=NULL,X0=NULL,dispersion = NULL, ...) {
  ## summary method for plbpsm object - provides approximate p values for parametric part

  ######### calculate bands #######################
  if (length(object$basis_info)>0){
      term.smooth=sapply(object$basis_info, class)[1,]
      if (object$MI){
        #term.smooth=sapply(object$basis_info2, class)[1,]
        smooth.univariate=object$ind.nl
      } else {
        smooth.univariate=which(term.smooth=="univariate.smooth")
      }
      if (object$backfitting){
        if (length(smooth.univariate)>0){
          bands=list()
          for (alpha in 1:length(smooth.univariate)){
            ### change???
            alpha0=alpha+dim(as.matrix(object$Xp))[2]
            if (is.null(X0)){
              bands[[alpha]]=band_univariate_orig(object,object$mhat,object$X2,alpha=alpha,alpha0=alpha0,
                                                  smooth.univariate=smooth.univariate, h_opt=object$h_opt_all[alpha],...)
            } else{
              bands[[alpha]]=band_univariate_newmtx(object,object$mhat,object$X2,alpha=alpha,alpha0=alpha0,
                                                    smooth.univariate=smooth.univariate, h_opt=object$h_opt_all[alpha],x0=X0[,alpha],...)
            }
          }
        } else {
          bands=NULL
        }
      }else {
        bands=NULL
      }
  } else {
    bands=NULL
    mhat=NULL
    }

  #################
  pinv<-function(V,M,rank.tol=1e-6) {
    ## a local pseudoinverse function
    D <- eigen(V,symmetric=TRUE)
    M1<-length(D$values[D$values>rank.tol*D$values[1]])
    if (M>M1) M<-M1 # avoid problems with zero eigen-values

    if (M+1<=length(D$values)) D$values[(M+1):length(D$values)]<-1
    D$values<- 1/D$values
    if (M+1<=length(D$values)) D$values[(M+1):length(D$values)]<-0
    res <- D$vectors%*%(D$values*t(D$vectors))  ##D$u%*%diag(D$d)%*%D$v
    attr(res,"rank") <- M
    res
  } ## end of pinv

  # if (is.null(object$R)) {
  #   warning("p-values for any terms that can be penalized to zero will be unreliable: refit model to fix this.")
  #   useR <- FALSE
  # } else useR <- TRUE

  # if (p.type < -1) useR <- FALSE
  #
  # if (p.type!=0) warning("p.type!=0 is deprecated, and liable to be removed in future")

  p.table <- pTerms.table <- NULL

  covmat <- object$Ve
  # name <- names(object$edf) # change in plbpsm
  # dimnames(covmat) <- list(name, name)
  # covmat.unscaled <- covmat/object$sig2
  #est.disp <- object$scale.estimated
  # if (!is.null(dispersion)) {
  #   covmat <- dispersion * covmat.unscaled
  #   object$Ve <- object$Ve*dispersion/object$sig2 ## freq
  #   object$Vp <- object$Vp*dispersion/object$sig2 ## Bayes
  #   est.disp <- FALSE
  # } else dispersion <- object$sig2
  #

  ## Now the individual parameteric coefficient p-values...
  p.ind=which(!grepl("\\.",names(object$coefficients)))
  new.p.ind=p.ind[object$ind_c]
  if (length(covmat)>0){
    se <- diag(as.matrix(covmat))^0.5
  } else {se=1}
  residual.df<-sum(!is.na(object$fitted.values))-sum(object$edf)
  if (length(new.p.ind) > 0) { # individual parameters
    ### not consider here
    if (length(object$nsdf)>1) { ## several linear predictors
      pstart <- attr(object$nsdf,"pstart")
      ind <- rep(0,0)
      for (i in 1:length(object$nsdf)) if (object$nsdf[i]>0) ind <-
        c(ind,pstart[i]:(pstart[i]+object$nsdf[i]-1))
    } else {
      pstart <- 1;ind <- new.p.ind
      } ## only one lp
    #if (length(object$coefficients)>0){
      #if (names(object$coefficients)[1]=='(Intercept)'){ind_c=object$ind_c+1} else {ind_c=object$ind_c}
  #} else {ind_c=object$ind_c}
      p.coeff <- object$coefficients[new.p.ind]
      p.se <- se[seq(1,length(new.p.ind))]
      p.t<-p.coeff/p.se

    ##
    est.disp=FALSE
    if (!est.disp & length(p.ind)>0) {
      p.pv <- 2*pnorm(abs(p.t),lower.tail=FALSE)
      p.table <- cbind(p.coeff, p.se, p.t, p.pv)
      dimnames(p.table) <- list(names(p.coeff), c("Estimate", "Std. Error", "z value", "Pr(>|z|)"))
    } else {
      p.pv <- 2*pt(abs(p.t),df=residual.df,lower.tail=FALSE)
      p.table <- cbind(p.coeff, p.se, p.t, p.pv)
      dimnames(p.table) <- list(names(p.coeff), c("Estimate", "Std. Error", "t value", "Pr(>|t|)"))
    }
  } else {
    p.coeff <- p.t <- p.pv <- array(0,0)
    est.disp=FALSE}

  ## Next the p-values for parametric terms, so that factors are treated whole...

  pterms <- if (is.list(object$pterms)) object$pterms else list(object$pterms)
  if (!is.list(object$assign)) object$assign <- list(object$assign)
  #npt <- length(unlist(lapply(pterms,attr,"term.labels")))
  # if (names(object$coefficients)[1]=='intercept'){
  #   ind_c=object$ind_c+1
  # }
  npt <- length(object$ind_c)
  if (npt>0)  pTerms.df <- pTerms.chi.sq <- pTerms.pv <- array(0,npt)
  term.labels <- rep("",0)
  k <- 0 ## total term counter
  for (j in 1:length(pterms)) {
    ##term.labels <- attr(object$pterms,"term.labels")
    #tlj <- attr(pterms[[j]],"term.labels")
    tlj <- colnames(object$Ve)
    nt <- length(tlj)
    if (j>1 && nt>0) tlj <- paste(tlj,j-1,sep=".")
    term.labels <- c(term.labels,tlj)
    if (nt>0) { # individual parametric terms
      nt=np=length(object$ind_c)
      #np <- length(object$assign[[j]])
      #ind <- pstart[j] - 1 + 1:np
      ind=1:length(object$ind_c)
      Vb <- covmat[ind,ind,drop=FALSE]
      bp <- array(object$coefficients[new.p.ind],np)
      # pTerms.pv <- if (j==1) array(0,nt) else c(pTerms.pv,array(0,nt))
      # #attr(pTerms.pv,"names") <- term.labels
      # attr(pTerms.pv,"names") <- colnames(covmat)
      # pTerms.df <- pTerms.chi.sq <- pTerms.pv
      for (i in 1:length(new.p.ind)) {
        k <- k + 1
        #intercept
        # if(names(object$coefficients)[1]=='(Intercept)'){
        #   ind <- object$assign[[j]]==object$ind_c[i]-1
        # } else {ind <- object$assign[[j]]==object$ind_c[i]}
        b <- bp[i];V <- Vb[i,i]
        ## pseudo-inverse needed in case of truncation of parametric space
        if (length(b)==1) {
          V <- 1/V
          pTerms.df[k] <- nb <- 1
          pTerms.chi.sq[k] <- V*b*b
        } else {
          V <- pinv(V,length(b),rank.tol=.Machine$double.eps^.5)
          pTerms.df[k] <- nb <- attr(V,"rank")
          pTerms.chi.sq[k] <- t(b)%*%V%*%b
        }
        if (!est.disp)
          pTerms.pv[k] <- pchisq(pTerms.chi.sq[k],df=nb,lower.tail=FALSE)
        else
          pTerms.pv[k] <- pf(pTerms.chi.sq[k]/nb,df1=nb,df2=residual.df,lower.tail=FALSE)
      } ## for (i in 1:nt)
    } ## if (nt>0)
  }

  if (npt) {
    attr(pTerms.pv,"names") <- term.labels
    if (!est.disp) {
      pTerms.table <- cbind(pTerms.df, pTerms.chi.sq, pTerms.pv)
      dimnames(pTerms.table) <- list(term.labels, c("df", "Chi.sq", "p-value"))
    } else {
      pTerms.table <- cbind(pTerms.df, pTerms.chi.sq/pTerms.df, pTerms.pv)
      dimnames(pTerms.table) <- list(term.labels, c("df", "F", "p-value"))
    }
  } else { pTerms.df<-pTerms.chi.sq<-pTerms.pv<-array(0,0)}

  ## Now deal with the smooth terms....

  m <- length(object$basis_info) # number of smooth terms
  #
  # if (p.type < 0 ) {
  #   kmax <- 0
  #   for (i in 1:m) {
  #     start <- object$smooth[[i]]$first.para
  #     stop <- object$smooth[[i]]$last.para
  #     k <- stop-start+1
  #     if (k>kmax) kmax <- k
  #   }
  # }
  #
  # df <- edf1 <- edf <- s.pv <- chi.sq <- array(0, m)
  # if (m>0) # form test statistics for each smooth
  # { if (p.type < 5) { ## Bayesian p-values required
  #   if (useR)  X <- object$R else {
  #     sub.samp <- max(1000,2*length(object$coefficients))
  #     if (nrow(object$model)>sub.samp) { ## subsample to get X for p-values calc.
  #       seed <- try(get(".Random.seed",envir=.GlobalEnv),silent=TRUE) ## store RNG seed
  #       if (inherits(seed,"try-error")) {
  #         runif(1)
  #         seed <- get(".Random.seed",envir=.GlobalEnv)
  #       }
  #       kind <- RNGkind(NULL)
  #       RNGkind("default","default")
  #       set.seed(11) ## ensure repeatability
  #       ind <- sample(1:nrow(object$model),sub.samp,replace=FALSE)  ## sample these rows from X
  #       X <- predict(object,object$model[ind,],type="lpmatrix")
  #       RNGkind(kind[1],kind[2])
  #       assign(".Random.seed",seed,envir=.GlobalEnv) ## RNG behaves as if it had not been used
  #     } else { ## don't need to subsample
  #       X <- model.matrix(object)
  #     }
  #     X <- X[!is.na(rowSums(X)),] ## exclude NA's (possible under na.exclude)
  #   }
  # } ## end if (p.type<5)
  #
  #   for (i in 1:m) { ## loop through smooths
  #
  #     start <- object$smooth[[i]]$first.para;stop <- object$smooth[[i]]$last.para
  #
  #     if (p.type==5) { ## use frequentist cov matrix
  #       V <- object$Ve[start:stop,start:stop,drop=FALSE]
  #     } else V <- object$Vp[start:stop,start:stop,drop=FALSE] ## Bayesian
  #
  #     p <- object$coefficients[start:stop]  # params for smooth
  #
  #     edf1[i] <- edf[i] <- sum(object$edf[start:stop]) # edf for this smooth
  #     ## extract alternative edf estimate for this smooth, if possible...
  #     if (!is.null(object$edf1)) edf1[i] <-  sum(object$edf1[start:stop])
  #
  #     if (p.type==5) { ## old style frequentist
  #       M1 <- object$smooth[[i]]$df
  #       M <- min(M1,ceiling(2*sum(object$edf[start:stop]))) ## upper limit of 2*edf on rank
  #       V <- pinv(V,M) # get rank M pseudoinverse of V
  #       chi.sq[i] <- t(p)%*%V%*%p
  #       df[i] <- attr(V, "rank")
  #     } else { ## Better founded alternatives...
  #       Xt <- X[,start:stop,drop=FALSE]
  #       if (object$smooth[[i]]$null.space.dim==0&&!is.null(object$R)) { ## random effect or fully penalized term
  #         res <- reTest(object,i)
  #       } else { ## Inverted Nychka interval statistics
  #         df[i] <- min(ncol(Xt),edf1[i])
  #         if (est.disp) rdf <- residual.df else rdf <- -1
  #         res <- testStat(p,Xt,V,df[i],type=p.type,res.df = rdf)
  #       }
  #       df[i] <- res$rank
  #       chi.sq[i] <- res$stat
  #       s.pv[i] <- res$pval
  #     }
  #     names(chi.sq)[i]<- object$smooth[[i]]$label
  #
  #     if (p.type == 5) {
  #       if (!est.disp)
  #         s.pv[i] <- pchisq(chi.sq[i], df = df[i], lower.tail = FALSE)
  #       else
  #         s.pv[i] <- pf(chi.sq[i]/df[i], df1 = df[i], df2 = residual.df, lower.tail = FALSE)
  #       ## p-values are meaningless for very small edf. Need to set to NA
  #       if (df[i] < 0.1) s.pv[i] <- NA
  #     }
  #   }
  #   if (!est.disp) {
  #     if (p.type==5) {
  #       s.table <- cbind(edf, df, chi.sq, s.pv)
  #       dimnames(s.table) <- list(names(chi.sq), c("edf", "Est.rank", "Chi.sq", "p-value"))
  #     } else {
  #       s.table <- cbind(edf, df, chi.sq, s.pv)
  #       dimnames(s.table) <- list(names(chi.sq), c("edf", "Ref.df", "Chi.sq", "p-value"))
  #     }
  #   } else {
  #     if (p.type==5) {
  #       s.table <- cbind(edf, df, chi.sq/df, s.pv)
  #       dimnames(s.table) <- list(names(chi.sq), c("edf", "Est.rank", "F", "p-value"))
  #     } else {
  #       s.table <- cbind(edf, df, chi.sq/df, s.pv)
  #       dimnames(s.table) <- list(names(chi.sq), c("edf", "Ref.df", "F", "p-value"))
  #     }
  #   }
  # }
  # check!
  w <- as.numeric(object$prior.weights)
  mean.y <- sum(w*object$y,na.rm = TRUE)/sum(w,na.rm = TRUE)
  w <- sqrt(w)
  #   nobs <- nrow(object$model)
  nobs <- sum(!is.na(object$y))
  # need to adjust for NA values
  r.sq <- if (inherits(object$family,"general.family")||!is.null(object$family$no.r.sq)) NULL else
    1 - var(w*(as.numeric(object$y)-object$fitted.values),na.rm = TRUE)*(nobs-1)/(var(w*(as.numeric(object$y)-mean.y),na.rm = TRUE)*residual.df)
  dev.expl<-(object$null.deviance-object$deviance)/object$null.deviance
  #if (object$method%in%c("REML","ML")) object$method <- paste("-",object$method,sep="")
  ret<-list(p.coeff=p.coeff,se=se,p.t=p.t,p.pv=p.pv,residual.df=residual.df,m=m,#chi.sq=chi.sq,
            #s.pv=s.pv,#scale=dispersion,
            r.sq=r.sq,family=object$family,formula=object$formula,n=nobs,
            dev.expl=dev.expl,#edf=edf,dispersion=dispersion,
            pTerms.pv=pTerms.pv,pTerms.chi.sq=pTerms.chi.sq,
            pTerms.df = pTerms.df, #cov.unscaled = covmat.unscaled,
            cov.scaled = covmat, p.table = p.table,
            pTerms.table = pTerms.table, #s.table = s.table,
            method=object$method,criterion=object$criterion,gcv_opt=object$gcv_opt,cv_opt=object$cv_opt,
            #rank=object$rank,
            np=length(object$coefficients),bands=bands,mhat=object$mhat)

  class(ret)<-"summary.plbpsm"
  ret
} ## end summary.plbpsm


#' @rdname summary.plbpsm
#' @importFrom stats printCoefmat
#' @method print summary.plbpsm
#' @export
print.summary.plbpsm <- function (x, digits = max(3, getOption("digits") - 3), signif.stars = getOption("show.signif.stars"),...)
{
  print(x$family)
  cat("Formula:\n")
  if (is.list(x$formula))
    for (i in 1:length(x$formula)) print(x$formula[[i]])
  else print(x$formula)
  if (length(x$p.coeff) > 0) {
    cat("\nParametric coefficients:\n")
    printCoefmat(x$p.table, digits = digits, signif.stars = signif.stars,
                 na.print = "NA")
  }
  cat("\n")
  # if (x$m > 0) {
  #   cat("Approximate significance of smooth terms:\n")
  #   printCoefmat(x$s.table, digits = digits, signif.stars = signif.stars,
  #                has.Pvalue = TRUE, na.print = "NA", cs.ind = 1)
  # }
  cat("\n")
  if (!is.null(x$rank) && x$rank < x$np)
    cat("Rank: ", x$rank, "/", x$np, "\n", sep = "")
  if (!is.null(x$r.sq))
    cat("R-sq.(adj) = ", formatC(x$r.sq, digits = 3, width = 5),
        "  ")
  if (length(x$dev.expl) > 0)
    cat("Deviance explained = ", formatC(x$dev.expl * 100,
                                         digits = 3, width = 4), "%", sep = "")
  cat("\n")
  if (!is.null(x$criterion) && (x$criterion %in% c("GCV")) && !is.null(x$gcv_opt))
    cat(x$criterion, " = ", formatC(x$gcv_opt, digits = 5),
        sep = "")
  if (!is.null(x$criterion) && (x$criterion %in% c("CV"))&& !is.null(x$cv_opt))
    cat(x$criterion, " = ", formatC(x$cv_opt, digits = 5),
        sep = "")
  if(x$family[1]!="gaussian"){
    cat(#"  Scale est. = ", formatC(x$sigma_2, digits = 5, width = 8,flag = "-"),
        "  n = ", x$n, "\n", sep = "")} else {
                                   cat("  n = ", x$n, "\n", sep = "")
                                 }


  ## Look into this!

  invisible(x)
}
