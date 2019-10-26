#' Partial Linear Bivariate Penalized Spline Model
#'
#' Fits a partial linear bivariate penalized spline model
#' @importFrom stats gaussian model.frame quantile reformulate residuals
#' @importFrom utils getFromNamespace
#' @param formula a PLBPSM formula. These are exactly like the formula for
#' a GLM except that univariable and bivariate smooth terms.
#' @param family This is a family object specifying the distribution
#' and link to use in fitting etc (see \code{\link{glm}} and \code{\link{family}}).
#' @param data A data frame or list containing the model response variable and
#' covariates required by the formula. By default the variables are taken from
#' environment(\code{formula}).
#' @param ind_c Defualt set to '\code{NULL}', if not, it is the arbitrary chosen variables from the parametric part
#' @param weights weights on the contribution of the data to the log likelihood.
#' Note that a weight of 2, for example, is equivalent to having made exactly the same observation twice.
#' If you want to reweight the contributions of each datum without changing
#' the overall magnitude of the log likelihood, then you should normalize the weights
#' (e.g. \code{weights <- weights/mean(weights)}).
#' @param na.action a function which indicates what should happen when the data contain ‘NA’s.
#'  The default is set by the 'na.action' setting of 'options',
#' and is 'na.fail' if that is unset. The 'factory-fresh' default is 'na.omit'.
#' @param offset Can be used to supply a model offset for use in fitting.
#' Note that this offset will always be completely ignored when predicting, unlike an offset
#' included in formula (this used to conform to the behaviour of lm and \code{glm}).
#' @param criterion The criterion to choose the penalty parameter lambda. \code{"GCV"} to use
#' generalized cross validation method and \code{"CV"} for cross validation
#' @param method The variable selection method for linear covariates. \code{"SCAD"} is the SCAD method in penalizing
#' the coefficients for linear covariates. \code{"ALASSO"} is the adaptive
#' LASSO method in penalizing the coeffiecient for linear covariates
#' @param scale scale parameter in exponential family
#' @param VS '\code{TRUE}' for using ALASSO/SCAD to select linear variables
#' @param fit  If this argument is \code{TRUE} then \code{gam} sets up the model and fits it, but if it is
#' \code{FALSE} then the model is set up and an object \code{G} containing what
#' would be required to fit is returned is returned. See argument \code{G}.
#' @param G Usually \code{NULL}, but may contain the object returned by a previous call to \code{gam} with
#' \code{fit=FALSE}, in which case all other arguments are ignored except for
#' \code{scale}, \code{control}, \code{method} and \code{fit}.
#' @param drop.unused.levels by default unused levels are dropped from factors before fitting. For some smooths
#' involving factor variables you might want to turn this off. Only do so if you know what you are doing.
#' @param drop.intercept Set to TRUE to force the model to really not have the a constant in the parametric model part,
#' even with factor variables present. Can be vector when formula is a list.
#' @param ... further arguments for passing on e.g. to \code{grplsfit}.
#' @param control A list of fit control parameters to replace defaults returned by \code{\link{plbpsm.control}}.
#' Any control parameters not supplied stay at their default values.
#' @param group A vector describing the grouping of the coefficients. For greatest efficiency and least ambiguity (see details),
#' it is best if group is a factor or vector of consecutive integers, although unordered groups and character vectors are also allowed. If there
#' are coefficients to be included in the model without being penalized, assign them to group 0 (or "0").
#' @param ecdf the choice of whether empirical conditional density function is used.
#' @param backfitting whether spline backfitted local polynomial estimation is applied.
#'@return If \code{fit=FALSE} the function returns a list \code{G} of items needed to fit a PLBPSM, but doesn't actually fit it.
#'Otherwise the function returns an object of class "\code{plbpsm}" as described in \code{\link{plbpsmObject}}.


#' @details A generalized geoadditive model (GgAM) is a generalized linear model (GLM) in which the linear
#' predictor is given by user specified univariate functions plus additive functions and
#' a bivariate function of the location covariates of the linear predictor. A simple example is:
#' \deqn{\log(E(y_i)) = u(z_1) + u(z_2) + f_{(x_{1i},x_{2i})}}{log(E(y_i))= u(z_1) + u(z_2) + f_{(x_{1i},x_{2i})}}
#' where the (independent) response variables \eqn{y_i \sim {\rm Poi }}{y_i~Poi}, \eqn{z_1}{z_1} and
#' \eqn{z_2}{z_2} are explantary variables,
#' \eqn{f}{f} is the bivariate smooth function of covariates \eqn{x_{1i}}{x_{1i}} and
#' \eqn{x_{2i}}{x_{2i}}. The log is an example of a link function.
#'
#'
#' A Generalized Geoadditive Model (GgAM) is a generalized linear model (GLM) in which the linear predictor is given by a user
#' specified bivariate functions of the location covariates plus a conventional parametric component of the linear predictor.  A simple example is:
#' \deqn{\log(E(y_i)) = z_1 + z_2 + u_{x_1} + u_{x_2} + b_{(s_{1i},s_{2i})}}{log(E(y_i))= z_1 + z_2 + u_{x_1} + u_{x_2} + b_{(s_{1i},s_{2i})}}
#' where the (independent) response variables \eqn{y_i \sim {\rm Poi }}{y_i~Poi}, \eqn{z_1}{z_1} and
#' \eqn{z_2}{z_2} are explantary variables,
#' \eqn{f}{f} is the bivariate smooth function of covariates \eqn{x_{1i}}{x_{1i}} and
#' \eqn{x_{2i}}{x_{2i}}. The log is an example of a link function.
#'
#' Model structure identification process is contained in \code{plbpsm} before model fitting to identify the linear
#' and nonlinear continuous variables by using group adaptive LASSO.
#'
#' We incorporate a variable selection
#' mechanism into the Partial linear Spatial Regression Model (PLSM) and propose a double penalized least squares approach based
#' on bivariate spline approximation over the spatial domain when the link is gaussian.
#'
#' \code{plbpsm} in \code{GgAM} chooses the smoothing penalty parameter by using the
#' Generalized Cross Validation (GCV) criterion.
#'
#' In terms of shrinkage penalty, Adaptive LASSO and SCAD penaltiy could be used to do variable selection under PLSM and GPLSM.
#' Broadly \code{plbpsm} works by first constructing basis functions and penalty
#' coefficient matrices for bivariate smooth term in the model formula, obtaining a model matrix for
#' the parametric part of the model formula. The algorithm chooses penalty parameter based on GCV/CV criterion
#' for the bivariate smoothing part and chooses significant linear covariates based on adaptive LASSO and
#' SCAD method using BIC criterion. And then, the refit is applied to the choosen model to get the final fit.
#'
#' @references
#' Wood, S., & Wood, M. S. (2015). Package ‘mgcv’. R package version, 1, 29.
#'
#' Yu, S., Wang, G., Wang, L., Liu, C., & Yang, L. (2019). Estimation and inference for generalized geoadditive models.
#' Journal of the American Statistical Association, 1-27.
#'
#' Wang, G., & Wang, J. (2019). On selection of semiparametric spatial regression models. Stat, 8(1), e221.
#'
#' Breheny P (2016).grpreg: Regularization Paths for Regression Models with Grouped Covari-ates.Rpackage
#' version 3.0-2, URLhttps://CRAN.R-project.org/packages=grpreg.
#' 
#' Wang L, Wang G, Li X, Mu J, Yu S, Wang Y, Kim M, Wang J (2019).BPST: Smoothing viaBivariate Spline over
#' Triangulation.Rpackage version 1.0,  URLhttps://GitHub.com/funstatpackages/BPST.
#' @examples
#' library(MASS)
#' library(grpreg)
#' library(mgcv)
#' library(Triangulation)
#' library(BPST)
#'
#' #################### Horseshoe Domain Example ######################
#' ###### VS in PLSM ######
#' data("eg1pop_dat")
#' eg1_V2 <- eg1pop_dat[['V2']]
#' eg1_T2 <- eg1pop_dat[['T2']]
#' eg1pop_rho03 <- eg1pop_dat[['rho03']]
#' n <- 100
#' Npop <- nrow(eg1pop_rho03)
#' # ind.pop <- (1:Npop)[!is.na(eg1pop_rho03[,1])]
#' ind.pop <- (1:Npop)
#' set.seed(1234)
#' sam.ind <- sort(sample(ind.pop,n,replace=FALSE))
#' sam <- eg1pop_rho03[sam.ind,]
#' lambda <- 10^(seq(-2,5,by=1))
#' data <- sam
#' data("horseshoe")
#' V <- TriMesh(horseshoe,4)$V
#' Tr <- TriMesh(horseshoe,4)$Tr
#'
#' formula <- Y~z1+z2+z3+z4+z5+z6+z7+z8+b(x1,x2,V=V,Tr=Tr,d=2,r=1)
#' # example 0:default lambda
#' res0 <- plbpsm(formula=formula,data=as.data.frame(data),VS=TRUE)
#' predict(res0,as.data.frame(data))
#' res0$fitted.values
#' summary(res0)
#'
#' ### Estimation ###
#' formula <- Y~z1+z2+z3+z4+z5+z6+z7+z8+b(x1,x2,V=eg1_V2,Tr=eg1_T2,lambda=lambda,d=2,r=1)
#' res <- plbpsm(formula=formula,data=as.data.frame(data),VS=TRUE)
#'
#' # ind_c is provided
#' res4 <- plbpsm(formula=formula,data=as.data.frame(sam),ind_c=c(1,2,4),VS=TRUE)
#'
#' ###### GPLSM ######
#' data(eg1pop_bin_rho00)
#' formula <- Y~z1+z2+z3+b(x1,x2,V=eg1_V2,Tr=eg1_T2,lambda=lambda,d=2,r=1)
#' n <- 100
#' Npop <- nrow(eg1pop_bin_rho00)
#' # ind.pop <- (1:Npop)[!is.na(eg1pop_bin_rho00[,1])]
#' ind.pop <- (1:Npop)
#' sam.ind <- sort(sample(ind.pop,n,replace=FALSE))
#' sam <- eg1pop_bin_rho00[sam.ind,]
#' data <- sam
#' res_eg1_bin <- plbpsm(formula=formula,data=as.data.frame(data),family=binomial())
#'
#' data(eg1pop_nb_rho00)
#' formula <- Y~z1+z2+z3+b(x1,x2,V=eg1_V2,Tr=eg1_T2,d=2,r=1)
#' n <- 100
#' Npop <- nrow(eg1pop_bin_rho00)
#' # ind.pop <- (1:Npop)[!is.na(eg1pop_bin_rho00[,1])]
#' ind.pop <- (1:Npop)
#' sam.ind <- sort(sample(ind.pop,n,replace=FALSE))
#' sam <- eg1pop_bin_rho00[sam.ind,]
#' res_nb <- plbpsm(formula=formula,data=as.data.frame(sam),family=negbin(c(2:8)))
#'
#' ### GgAM ####
#' # Poisson #
#' data(eg1pop_poi2)
#' n <- 1000
#' set.seed(2008)
#' Npop <- nrow(eg1pop_poi2)
#' # ind.pop <- (1:Npop)[!is.na(eg1pop_poi2[,1])]
#' ind.pop <- (1:Npop)
#' sam.ind <- sort(sample(ind.pop,n,replace=FALSE))
#' sam <- eg1pop_poi2[sam.ind,]
#' formula <- Y~u(z1)+u(z2)+u(z3)
#' res_eg1_poi_add0 <- plbpsm(formula=formula,data=as.data.frame(sam),family=quasipoisson(),
#' group=c(0,0,0))
#' summary(res_eg1_poi_add0)
#' formula <- Y~u(z1)+u(z2)+u(z3)+b(x1,x2,V=eg1_V2,Tr=eg1_T2,d=2,r=1)
#'
#' ### Estimation ###
#' res_eg1_poi_add <- plbpsm(formula=formula,data=as.data.frame(sam),family='poisson',group=c(0,0,0))
#' res_intercptonly <- plbpsm(formula=Y~1,data=as.data.frame(sam),family='poisson')
#' ### Inference & Plotting ###
#' summary(res_eg1_poi_add)
#' plot(res_eg1_poi_add)
#'
#' ### GgAM Prediction ###
#' pred <- predict(res_eg1_poi_add,as.data.frame(eg1pop_poi2[1:200,]))
#' origy <- eg1pop_poi2[1:200,'Y']
#' origy[which(is.na(pred))] <- NA
#' mean((pred-origy)^2,na.rm=TRUE)
#' # 2.8527
#'
#' ### GGAM-SMILE: Nonlinear detection ###
#' data('dat_poi_ggams')
#' formula <- y~u(x1)+u(x2)+u(x3)+b(s1,s2,V=V,Tr=Tr,d=2,r=1)
#' res_poi <- plbpsm(formula=formula,data=as.data.frame(dat_poi_ggams),family = "poisson")
#' # compare to the following Known GGAM model #
#' formula <- y~x1+u(x2)+u(x3)+b(s1,s2,V=V,Tr=Tr,d=2,r=1)
#' res_poi <- plbpsm(formula=formula,data=as.data.frame(dat_poi_ggams),family = "poisson",group=c(0,0))
#' # end of horseshoe domain example
#'
#' ################### Rectangular Domain Example #######################
#' ### VS in PLSM ###
#' data("eg2pop_dat")
#' eg2_V20=eg2pop_dat[['V20']]
#' eg2_T20=eg2pop_dat[['T20']]
#' eg2pop=eg2pop_dat[['pop']]
#' n=100
#' Npop=nrow(eg2pop)
#' # ind.pop=(1:Npop)[!is.na(eg2pop[,1])]
#' ind.pop <- (1:Npop)
#' sam.ind <- sort(sample(ind.pop,n,replace=FALSE))
#' sam=eg2pop[sam.ind,]
#' formula <- Y~z1+z2+z3+z4+z5+z6+z7+z8+b(x1,x2,V=eg2_V20,Tr=eg2_T20)
#' res0 <- plbpsm(formula=formula,data=as.data.frame(sam),VS=TRUE)
#'
#' ### More complicated
#' formula <- Y~z2+u(z1)+b(x1,x2,V=eg2_V20,Tr=eg2_T20)
#' res_ex <- plbpsm(formula=formula,data=as.data.frame(data),group=c(0,0))
#'
#' @export

plbpsm=function(formula,family = gaussian(), data = list(), ind_c=NULL, group=NULL,
                weights = NULL,na.action, offset = NULL, criterion = "GCV", method='SCAD',control=list(),
                scale = 1, VS = FALSE, fit = TRUE, G = NULL,
                drop.unused.levels = TRUE, drop.intercept = NULL,ecdf = FALSE, backfitting=TRUE,...){
  ### G is usually NULL, but may contain the object returned by a previous call to plbpsm with
  control <- do.call("plbpsm.control", control)
  if (is.null(G)) {
    gp <- interpret.plbpsm0(formula)
    cl <- match.call()
    mf<-match.call()
    mf <- match.call(expand.dots = FALSE)
    mf$formula <- gp$fake.formula
    mf$family <- mf$group <- mf$control <- mf$ind_c <- mf$scale <- mf$drop.intercept <-
      mf$method <- mf$criterion <- mf$se <- mf$VS <- mf$fit <- mf$G <- mf$ecdf <- mf$backfitting <- mf$... <- NULL
    mf$drop.unused.levels <- drop.unused.levels
    mf[[1]] <- quote(stats::model.frame)
    pmf <- mf
    mf <- eval(mf, parent.frame())
    #mf <- model.frame(mf,data)
    if (nrow(mf) < 2)
      stop("Not enough (non-NA) data to do anything meaningful")
    terms <- attr(mf, "terms")
    ###
    univ=sapply(gp$smooth.spec,function(x){class(x)})
    ind_univ=which(univ=="univariate.smooth")
    # if (length(ind_univ)>0){
    #   ecdf <- TRUE
    # }
    if (ecdf){
      ed_trans=function(x){
        Fn=ecdf(x)
        Fn(x)
      }
      data2=as.data.frame(data)
      names=attr(terms,"term.labels")
      loc=gp$smooth.spec[[length(gp$smooth.spec)]]$term
      names2=names[-which(names %in% loc )]
      X_ed=apply(data2[,names2],2,ed_trans)
      data22=data.frame(mf[,1],X_ed,mf[,loc])
      names(data22)<-names(attr(terms,"dataClasses"))
      data=data22
      mf22 <- model.frame(mf,data22)
      #mf <- eval(mf, parent.frame())
    } else {mf22=mf}
    ### all.vars1 in mgcv since deal with dat$
    # all_vars1 <- getFromNamespace("all.vars1", "mgcv")
    # vars <- all_vars1(gp$fake.formula[-2]) ## drop response here
    # inp <- parse(text = paste("list(", paste(vars, collapse = ","),")"))
    
    ## allow a bit of extra flexibility in what `data' is allowed to be (as model.frame actually does)
    if (!is.list(data)&&!is.data.frame(data)) data <- as.data.frame(data)
    
    # dl <- eval(inp, data, parent.frame())
    # names(dl) <- vars ## list of all variables needed
    # var.summary <- lapply(dl,
    #                       function(X){c(min(X,na.rm=TRUE),quantile(X,probs=0.5,type=3,na.rm=TRUE),
    #                                     max(X,na.rm=TRUE))})
    # variable_sum <- getFromNamespace("variable.summary", "mgcv")
    
    #var.summary<-variable_sum(gp$pf, dl, nrow(mf))
    #var.summary<- lapply(names(var.summary),function(X){names(X)<-c('min',"median","max")})
    # var.summary <- variable.summary(gp$pf,dl,nrow(mf)) ## summarize the input data
    # rm(dl) ## save space
    
    ## pterms are terms objects for the parametric model components used in
    ## model setup - don't try obtaining by evaluating pf in mf - doesn't
    ## work in general (e.g. with offset)...
    
    ## single linear predictor case
    pmf$formula <- gp$pf
    pmf <- eval(pmf, parent.frame())
    #pmf <- model.frame(pmf$formula, data)
    pterms <- attr(pmf, "terms") ## change
    
    if (is.character(family))
      family <- eval(parse(text = family))
    if (is.function(family))
      family <- family()
    if (is.null(family$family))
      stop("family not recognized")
    ## not sure useful or not
    if (family$family[1] == "gaussian" && family$link ==
        "identity")
      am <- TRUE else am <- FALSE
    
    # if (!control$keepData) rm(data) ## save space
    
    
    # Deal with intercept information
    if (is.null(family$drop.intercept)) {
      lengthf <- if (is.list(formula))
        length(formula) else 1
      if (is.null(drop.intercept))
        drop.intercept <- rep(FALSE, lengthf)
      else {
        drop.intercept <- rep(drop.intercept, length = lengthf)
        if (sum(drop.intercept))
          family$drop.intercept <- drop.intercept
      }
    } else drop.intercept <- as.logical(family$drop.intercept)
    
    # not sure whether we need this or not
    if (inherits(family, "general.family") && !is.null(family$presetup))
      eval(family$presetup)
    
    ### start to get all the model set up information
    gsname <- if (is.list(formula))
      ### have not done  "plbpsm.setup.list" ,not sure whether needed
      "plbpsm.setup.list" else "plbpsm.setup"
    G <- do.call(gsname, list(formula = gp, pterms = pterms,
                              data = mf22,  drop.intercept = drop.intercept,group=group))
    #
    #data = mf, data_ecdf = mf22
    #G$var.summary <- var.summary
    #rownames(G$var.summary) <- c('min',"median","max")
    G$family <- family
    if ((is.list(formula) && (is.null(family$nlp) || family$nlp !=
                              gp$nlp)) || (!is.list(formula) && !is.null(family$npl) &&
                                           (family$npl > 1)))
      stop("incorrect number of linear predictors for family")
    G$scale <- scale
    G$terms <- terms
    G$mf <- mf22
    G$cl <- cl
    G$am <- am
    G$backfitting <- backfitting
    
    if (is.null(G$offset))
      G$offset <- rep(0, G$n)
    # G$min.edf <- G$nsdf
    # if (G$m)
    #   for (i in 1:G$m) G$min.edf <- G$min.edf + G$basis_info[[i]]$null.space.dim
    ## check what is null space
    G$formula <- formula
    G$pred.formula <- gp$pred.formula
    environment(G$formula) <- environment(formula)
  }
  if (!fit)
    return(G)
  if (ncol(G$X) > nrow(G$X))
    stop("Model has more coefficients than data")
  ## main algorithm inside here
  object <- estimate.plbpsm(G,criterion,method,ind_c,VS,control,G$group0)
  
  # need to figure out which inputs we want to keep
  
  object$formula <- G$formula
  if (is.list(object$formula))
    attr(object$formula, "lpi") <- attr(G$X, "lpi")
  object$var.summary <- G$var.summary
  object$cmX <- G$cmX
  object$model <- mf
  object$model_ecdf <- mf22
  object$na.action <- attr(G$mf, "na.action")
  if(!G$am){
    object$control <- control
  }
  object$terms <- G$terms
  object$pred.formula <- G$pred.formula
  attr(object$pred.formula, "full") <- reformulate(all.vars(object$terms))
  object$pterms <- G$pterms
  object$assign <- G$assign
  object$contrasts <- G$contrasts
  object$xlevels <- G$xlevels
  object$offset <- G$offset
  if (!is.null(G$Xcentre))
    object$Xcentre <- G$Xcentre ###???
  object$ecdf=ecdf
  ## check
  object$df.residual <- nrow(G$X) - sum(object$min.edf)
  ##
  object$min.edf <- G$min.edf
  object$call <- G$cl
  class(object) <- c("plbpsm", "glm", "lm")
  if (is.null(object$deviance))
    object$deviance <- sum(residuals(object, "deviance")^2)
  if (!is.null(object$gcv_opt)){  names(object$gcv_opt) <- criterion}
  environment(object$formula) <- environment(object$pred.formula) <- environment(object$terms) <- environment(object$pterms) <- .GlobalEnv
  if (!is.null(object$model))
    environment(attr(object$model, "terms")) <- .GlobalEnv
  if (!is.null(attr(object$pred.formula, "full")))
    environment(attr(object$pred.formula, "full")) <- .GlobalEnv
  object
}
