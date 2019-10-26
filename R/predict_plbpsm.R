#' Prediction from fitted PLBPSM model
#'
#' Takes a fitted \code{plbpsm} object produced by \code{plbpsm()}
#' and produces predictions given a new set of values for the model covariates
#' or the original values used for the model fit. Predictions can be accompanied
#' by standard errors, based on the posterior distribution of the model
#' coefficients. The routine can optionally return the matrix by which the model
#' coefficients must be pre-multiplied in order to yield the values of the linear predictor at
#' the supplied covariate values.
#' @importFrom stats na.pass reformulate model.frame .checkMFClasses model.matrix model.offset napredict delete.response runif
#' @param object a fitted \code{plbpsm} object as produced by \code{plbpsm()}
#' @param newdata  A data frame or list containing the values of the model covariates at which predictions
#' are required. If this is not provided then predictions corresponding to the
#' original data are returned. If \code{newdata} is provided then
#' it should contain all the variables needed for prediction: a
#' warning is generated if not.
#' @param type the type of prediction required. The default is on the scale of the linear predictors;
#'  the alternative \code{"response"} is on the scale of the response variable. Thus for a default binomial
#'  model the default predictions are of log-odds (probabilities on logit scale) and
#'  \code{type = "response"} gives the predicted probabilities. The \code{"terms"} option returns
#'  a matrix giving the fitted values of each term in the model formula on the linear predictor scale.
#'  When this has the value \code{"link"} the linear predictor (possibly with associated standard errors) is returned.
#' @param se.fit when this is \code{TRUE} (not default) standard error estimates are returned for each prediction.
#' @param terms if \code{type=="terms"} then only results for the terms given in this array
#' will be returned.
#' @param exclude if \code{type=="terms"} or \code{type="iterms"} then terms (smooth or parametric) named in this
#'  array will not be returned. Otherwise any smooth terms named in this array will be set to zero.
#'   If \code{NULL} then no terms are excluded. Note that this is the term names as it appears in the model summary, see example.
#' @param block.size maximum number of predictions to process per call to underlying
#' code: larger is quicker, but more memory intensive. Set to < 1 to use total number
#' of predictions as this. If \code{NULL} then block size is 1000 if new data supplied,
#' and the number of rows in the model frame otherwise.
#' @param newdata.guaranteed Set to \code{TRUE} to turn off all checking of
#' \code{newdata}: this can speed things up
#' for large prediction tasks, but \code{newdata} must be complete, with no
#' \code{NA} values for predictors required in the model.
#' @param na.action what to do about \code{NA} values in \code{newdata}. With the
#' default \code{na.pass}, any row of \code{newdata} containing \code{NA} values
#' for required predictors, gives rise to \code{NA} predictions (even if the term concerned has no
#'                                                               \code{NA} predictors). \code{na.exclude} or \code{na.omit} result in the
#' dropping of \code{newdata} rows, if they contain any \code{NA} values for
#' required predictors. If \code{newdata} is missing then \code{NA} handling is
#' determined from \code{object$na.action}.
#' @param unconditional if TRUE then the smoothing parameter uncertainty corrected covariance matrix
#' is used to compute uncertainty bands, if available. Otherwise the bands treat the
#' smoothing parameters as fixed.
#' @param newB the given matrix of bivariate spline basis functions.
#' @param newind00 the given index of the data points in the triangulation.
#' @param backfitting whether SBL estimation is obtained.
#' @param ... other arguments.
#' @return if \code{se.fit} is \code{TRUE} then a 2 item list is returned with items (both arrays) \code{fit}
#' and \code{se.fit} containing predictions and associated standard error estimates, otherwise an
#' array of predictions is returned. The dimensions of the returned arrays depends on whether
#' \code{type} is \code{"terms"} or not: if it is then the array is 2 dimensional with each
#' term in the linear predictor separate, otherwise the array is 1 dimensional and contains the
#' linear predictor/predicted values (or corresponding s.e.s). The linear predictor returned termwise will
#' not include the offset or the intercept.
#' @details
#' See examples to see usages of different types.
#' @examples
#' library(MASS)
#' library(grpreg)
#' library(BPST)
#' data("eg1pop_dat")
#' eg1_V2 <- eg1pop_dat[['V2']]
#' eg1_T2 <- eg1pop_dat[['T2']]
#' eg1pop_rho03 <- eg1pop_dat[['rho03']]
#' sam <- eg1pop_rho03[sample(1:dim(eg1pop_rho03)[1],100),]
#'
#' ### Partial Linear Spatial Model ###
#' formula_d4 <- Y~z1+z2+z3+z4+z5+z6+z7+z8+b(x1,x2,d=4,r=1,V=eg1_V2,Tr=eg1_T2)
#' res <- plbpsm(formula=formula_d4,data=as.data.frame(sam),VS=TRUE)
#' # check deviance and the estimated error from predict function
#' newdata <- as.data.frame(sam[sample(1:dim(sam)[1],90),])
#' pred <- predict(res,newdata)
#'
#' # exclude one covariate
#' pred0 <- predict(res,newdata,exclude="z2",type='terms')
#' # check
#'
#' ### Generalized Partially Linear Spatial Model ###
#' # Poisson family
#' data("eg_poi_pl")
#' sam <- as.data.frame(eg_poi[['sam_poi']])
#' sam[1,]
#' V <- eg_poi[['V1']]
#' Tr <- eg_poi[['T1']]
#' formula <- y~c1+c2+c3+c4+c5+c6+c7+c8+c9+c10+c11+c12+c13+c14+c15+b(loc1,loc2,V=V,Tr=Tr,d=2,r=1)
#' res_poi <- plbpsm(formula=formula,family=poisson(),data=as.data.frame(sam))
#'
#' # with offset
#' res_poi2 <- plbpsm(formula=formula,family=poisson(),offset=rep(0.1,length(data$y)),
#' data <- as.data.frame(sam))
#' newdata <- as.data.frame(sam[sample(1:dim(sam)[1],90),])
#'
#' # return the predicted value for each term #
#' predict(res_poi,newdata,type='terms',se.fit = TRUE)
#' predict(res_poi2,newdata)
#'
#' data(eg1pop_bin_rho00)
#' formula <- Y~z1+z2+z3+b(x1,x2,V=eg1_V2,Tr=eg1_T2,d=2,r=1)
#' n <- 100
#' Npop <- nrow(eg1pop_bin_rho00)
#' # ind.pop <- (1:Npop)[!is.na(eg1pop_bin_rho00[,1])]
#' ind.pop <- (1:Npop)
#' sam.ind <- sort(sample(ind.pop,n,replace=FALSE))
#' sam <- eg1pop_bin_rho00[sam.ind,]
#' res_eg1_bin <- plbpsm(formula=formula,data=as.data.frame(sam),family=binomial())
#' eg1pop_bin_rho00=as.data.frame(eg1pop_bin_rho00)
#' newdata <- eg1pop_bin_rho00[1000:1100,]
#' predict(res_eg1_bin,as.data.frame(newdata),type='terms',se.fit = TRUE)
#'
#' ### Generalized Geoadditive Models with Model Identification ###
#' data(eg1pop_poi2)
#' formula <- Y~u(z1)+u(z2)+u(z3)+b(x1,x2,V=eg1_V2,Tr=eg1_T2,d=2,r=1)
#' n <- 100
#' Npop <- nrow(eg1pop_poi2)
#' # ind.pop <- (1:Npop)[!is.na(eg1pop_poi2[,1])]
#' ind.pop <- (1:Npop)
#' sam.ind <- sort(sample(ind.pop,n,replace=FALSE))
#' sam <- eg1pop_poi2[sam.ind,]
#' res_eg1_poi_add <- plbpsm(formula=formula,data=as.data.frame(sam),family='poisson')
#' newdata <- as.data.frame(sam[sample(1:dim(sam)[1],90),])
#'
#' ## backfitting = FALSE ##
#' pred_noBKF <- predict(res_eg1_poi_add,newdata,backfitting=FALSE)
#'
#' ## defualt backfitting for \code{res_eg1_poi_add}
#' pred_BKF <- predict(res_eg1_poi_add,newdata)

#'@export

predict.plbpsm <- function(object, newdata, type = "response", se.fit=FALSE,
                           terms = NULL, exclude = NULL,
                         block.size = NULL, newdata.guaranteed = FALSE, na.action = na.pass,
                         unconditional = FALSE, newB = NULL, newind00 = NULL,
                         backfitting = object$backfitting,...) {

  # This function is used for predicting from a PLBPSM 'object' is a plbpsm object, newdata a dataframe to
  # be used in prediction......
  #
  # Type == "link"     - for linear predictor (may be several for extended case)
  #      == "response" - for fitted values: may be several if several linear predictors,
  #                      and may return something other than inverse link of l.p. for some families
  #      == "terms"    - for individual terms on scale of linear predictor
  # Steps are:
  #  1. Set newdata to object$model if no newdata supplied
  #  2. split up newdata into manageable blocks if too large
  #  3. Obtain parametric model matrix (safely!)
  #  4. Work through nonparametric terms
  #  5. Work out required quantities
  #
  # if newdata.guaranteed == TRUE then the data.frame is assumed complete and
  # ready to go, so that only factor levels are checked for sanity.
  #
  # if `terms' is non null then it should be a list of terms to be returned
  # when type=="terms".
  # if `object' has an attribute `para.only' then only parametric terms of order
  # 1 are returned for type=="terms" : i.e. only what termplot can handle.
  #
  # if no new data is supplied then na.action does nothing, otherwise
  # if na.action == "na.pass" then NA predictors result in NA predictions (as lm
  #                   or glm)
  #              == "na.omit" or "na.exclude" then NA predictors result in
  #                       dropping

  # backfitting use its own algorithm
  # if (backfitting){
  #   H=predict_backfitting(object,as.data.frame(newdata),type=type,...)
  #   return(H)
  # } else {


  if (type!="link"&&type!="terms"&&type!="response")  {
    warning("Unknown type, reset to terms.")
    type<-"response"
  }
  if (!inherits(object,"plbpsm")) stop("predict.plbpsm can only be used to predict from plbpsm objects")

  ## to mimic behaviour of predict.lm, some resetting is required ...
  if (missing(newdata)) na.act <- object$na.action else {
    if (is.null(na.action)) na.act <- NULL
    else {
      #na.pass is a function
      na.txt <- "na.pass"
      if (is.character(na.action))
        na.txt <- substitute(na.action) else
          if (is.function(na.action)) na.txt <- deparse(substitute(na.action))
          if (na.txt=="na.pass") na.act <- "na.exclude" else
            if (na.txt=="na.exclude") na.act <- "na.omit" else
              na.act <- na.action
    }
  } ## ... done

  # get data from which to predict.....
  nd.is.mf <- FALSE # need to flag if supplied newdata is already a model frame
  ## get name of response...
  yname <- all.vars(object$terms)[attr(object$terms,"response")]
  if (newdata.guaranteed==FALSE) {
    if (missing(newdata)) { # then "fake" an object suitable for prediction
      newdata <- object$model
      new.data.ok <- FALSE
      nd.is.mf <- TRUE
      response <- newdata[[yname]]
    } else {  # do an R ``standard'' evaluation to pick up data
      new.data.ok <- TRUE
      if (is.data.frame(newdata)&&!is.null(attr(newdata,"terms"))) { # it's a model frame
        if (sum(!(names(object$model)%in%names(newdata)))) stop(
          "newdata is a model.frame: it should contain all required variables\n")
        nd.is.mf <- TRUE
      } else {
        yname <- all.vars(object$terms)[attr(object$terms,"response")]
        naresp <- FALSE
        if (!is.null(object$family$predict)&&!is.null(newdata[[yname]])) {
          ## response provided
          #if (!is.null(object$pred.formula)) object$pred.formula <- attr(object$pred.formula,"full")
          response <- TRUE
          Terms <- terms(object)
          resp <- newdata[[yname]]
          if (sum(is.na(resp))>0) {
            naresp <- TRUE ## there are NAs in supplied response
            ## replace them with a numeric code, so that rows are not dropped below
            rar <- range(resp,na.rm=TRUE)
            thresh <- rar[1]*1.01-rar[2]*.01
            resp[is.na(resp)] <- thresh
            newdata[[yname]] <- thresh
          }
        } else { ## response not provided
          response <- FALSE
          Terms <- delete.response(terms(object))
        }
        #allNames <- if (is.null(object$pred.formula)) all.vars(Terms) else all.vars(object$pred.formula)
        allNames <- all.vars(Terms)
        if (length(allNames) > 0) {
          ff <- if (is.null(object$pred.formula)) reformulate(allNames) else  object$pred.formula
          if (sum(!(allNames%in%names(newdata)))) {
            warning("not all required variables have been supplied in  newdata!\n")
          }
          newdata <- eval(model.frame(ff,data=newdata,na.action=na.act),parent.frame())
          if (naresp) newdata[[yname]][newdata[[yname]]<=thresh] <- NA ## reinstate as NA
        } ## otherwise it's intercept only and newdata can be left alone
        na.act <- attr(newdata,"na.action")
        response <- if (response) newdata[[yname]] else NULL
      }
    }
  } else { ## newdata.guaranteed == TRUE
    na.act <- NULL
    new.data.ok=TRUE ## it's guaranteed!
    if (!is.null(attr(newdata,"terms"))) nd.is.mf <- TRUE
    response <- newdata[[yname]]
  }

  ## now check the factor levels and split into blocks...

  if (new.data.ok) {
    # split prediction into blocks, to avoid running out of memory
    if (length(newdata)==1) newdata[[2]] <- newdata[[1]] # avoids data frame losing its labels and dimensions below!
    if (is.null(dim(newdata[[1]]))) np <- length(newdata[[1]])
    else np <- dim(newdata[[1]])[1]
    nb <- length(object$coefficients)+ as.numeric(object$intercept)
      #as.numeric(length(object$basis_info)>0 && object$intercept)
    #+ as.numeric(length(object$basis_info)==0 && object$intercept)
    if (is.null(block.size)) block.size <- 10000
    if (block.size < 1) block.size <- np
  } else { # no new data, just use object$model
    np <- nrow(object$model)
    nb <- length(object$coefficients)
  }

  ## split prediction into blocks, to avoid running out of memory
  if (is.null(block.size)) {
    ## use one block as predicting using model frame
    ## and no block size supplied...
    n.blocks <- 1
    b.size <- array(np,1)
  } else {
    n.blocks <- np %/% block.size
    b.size <- rep(block.size,n.blocks)
    last.block <- np-sum(b.size)
    if (last.block>0) {
      n.blocks <- n.blocks+1
      b.size[n.blocks] <- last.block
    }
  }

  term.smooth=sapply(object$basis_info, class)[1,]
  if (object$MI){
    smooth.univariate=object$ind.nl
  } else {
    smooth.univariate=which(term.smooth=="univariate.smooth")
  }
  smooth.bivariate=which(term.smooth=="bivariate.smooth")
  n.smooth1=length(smooth.univariate)
  n.smooth=length(smooth.univariate)+length(smooth.bivariate)
  # setup prediction arrays...
  ## in multi-linear predictor models, lpi[[i]][j] is the column of model matrix contributing the jth col to lp i
  #lpi <- if (is.list(object$formula)) attr(object$formula,"lpi") else NULL

  if (type=="terms") {
    term.labels <- attr(object$pterms,"term.labels")
    # if (is.null(attr(object,"para.only"))) para.only <-FALSE else
    #   para.only <- TRUE  # if true then only return information on parametric part
    n.pterms <- length(term.labels)
    if (object$MI){
      n.pterms=n.pterms+length(object$ind.l)
      term.labels=c(term.labels,object$linear.names)
    }
    fit <- array(NA,c(np,n.pterms+n.smooth))
    if (se.fit) se <- fit
    ColNames <- term.labels
  } else { ## "response" or "link"
    fit <- array(NA,np)
    if (se.fit) se <- fit
    #fit1 <- NULL ## "response" returned by fam$fv can be non-vector
  }
  ##se: 100 9
  stop <- 0
  # if (is.list(object$pterms)) { ## multiple linear predictors
  #   pstart <- attr(object$nsdf,"pstart") ## starts of parametric blocks in coef vector
  #   pind <- rep(0,0) ## index of parametric coefs
  #   Terms <- list();pterms <- object$pterms
  #   for (i in 1:length(object$nsdf)) {
  #     Terms[[i]] <- delete.response(object$pterms[[i]])
  #     if (object$nsdf[i]>0) pind <- c(pind,pstart[i]-1+1:object$nsdf[i])
  #   }
  # } else { ## normal single predictor case
    Terms <- list(delete.response(object$pterms)) ## make into a list anyway
    pterms <- list(object$pterms)
    #pstart <- 1
    #pind <- 1:object$nsdf ## index of parameteric coefficients
  #}

  ## check if extended family required intercept to be dropped...
  drop.intercept <- FALSE
  if (type=="terms") {
    if (object$MI){
      if (n.smooth1) for (k in 1:(n.smooth1)){
        #term.smooth=sapply(object$basis_info2, class)[1,]
        col.nonlinear=object$nonlinear.names
        k.names=object$basis_info[[k]]$label
        ColNames=c(ColNames,k.names)}
      if (length(smooth.bivariate)>0){
        ColNames=c(ColNames,basisinfo$label)
      }
    } else {
      if (n.smooth) for (k in 1:n.smooth){
        ColNames[n.pterms+k] <- object$basis_info[[k]]$label
      }
    }
  }

  ####################################
  ## Actual prediction starts here...
  ####################################
  #s.offset <- NULL # to accumulate any smooth term specific offset
  #any.soff <- FALSE # indicator of term specific offset existence
  if (n.blocks > 0) for (b in 1:n.blocks) { # work through prediction blocks
    start <- stop+1
    stop <- start + b.size[b] - 1
    if (n.blocks==1) data <- newdata else data <- newdata[start:stop,]
    X <- matrix(0,b.size[b],nb)
    Xoff <- matrix(0,b.size[b],n.smooth) ## term specific offsets
    for (i in 1:length(Terms)) { ## loop for parametric components (1 per lp)
      ## implements safe prediction for parametric part as described in
      ## http://developer.r-project.org/model-fitting-functions.txt
      if (new.data.ok) {
        if (nd.is.mf) mf <- model.frame(data) else {
          mf <- model.frame(Terms[[i]],data)
          if (!is.null(cl <- attr(pterms[[i]],"dataClasses"))) .checkMFClasses(cl,mf)
        }
        Xp <- model.matrix(Terms[[i]],mf,contrasts=object$contrasts)
      } else {
        Xp <- model.matrix(Terms[[i]],object$model)
        mf <- newdata # needed in case of offset, below
      }
      # if (drop.intercept) {
      #   xat <- attributes(Xp);ind <- xat$assign>0
      #   Xp <- Xp[,xat$assign>0,drop=FALSE] ## some extended families need to drop intercept
      #   xat$assign <- xat$assign[ind];xat$dimnames[[2]]<-xat$dimnames[[2]][ind];
      #   xat$dim[2] <- xat$dim[2]-1;attributes(Xp) <- xat
      # }
      if (object$MI){
        if (length(object$ind.l)>0){
            Xp=cbind(Xp,data[, sapply(object$basis_info_MI[object$ind.l],function(X){X$term})])
        }
        if (object$nsdf[i]>0) X[,1:(object$nsdf[i]+length(object$ind.l))] <- as.matrix(Xp)
      } else {if (object$nsdf[i]>0) X[,1:object$nsdf[i]] <- Xp}
    } ## end of parametric loop

    #if (!is.null(drop.ind)) X <- X[,-drop.ind]

    ### BPST Predictions ###

    if (backfitting){
      #newind=newind0
      res.fit <- predict_backfitting(object,data,type,ColNames,...)
      # newind <-res.fit$newind
      # newind2 <-(start:stop)[newind]
      # newind30 <-newind2[!is.na(newind2)]
      # newind3 <-newind30[apply(data, 1, function(x) !any(is.na(x))) ]
      # fit[newind3] <-res.fit$fit
      #print(res.fit$fit)
      fit <- res.fit$pred
    } else {
    if (length(smooth.bivariate)>0) {
      basisinfo <-object$basis_info[[smooth.bivariate]]
      # all the newdata below should be data
      loc <-data[,basisinfo$term]
      if (is.null(newB) | is.null(newind00)){
         B0 <-basis(as.matrix(basisinfo$V),as.matrix(basisinfo$Tr),basisinfo$d,basisinfo$r,as.matrix(loc))
         newB <-B0$B
         newind0 <-B0$Ind.inside
      } else {
        newind0 <-newind00
      }
   ## prediction of the nonlinear part
      # print(dim(newB))
      if (dim(newB)[2]>1){
        newbeta <-as.matrix(newB)%*%matrix(object$coefficients_bivariate,ncol=1)
        print(dim(newbeta))
      } else {
        newbeta <-matrix(newB,nrow=1)%*%matrix(object$coefficients_bivariate,ncol=1)
      }
      # print(newbeta)
      # Xfrag <- PredictMat(object$smooth[[k]],data)
      # X[,object$smooth[[k]]$first.para:object$smooth[[k]]$last.para] <- Xfrag
      #Xfrag.off <- attr(Xfrag,"offset") ## any term specific offsets?
      #if (!is.null(Xfrag.off)) { Xoff[,k] <- Xfrag.off; any.soff <- TRUE }
    } else {
      newind0 <-1:block.size
      newbeta <-0
    }

    ### univariate new basis
    first.para <-dim(Xp)[2]
    if (n.smooth1>0) for (k in 1:n.smooth1){
      basisinfo <-object$basis_info[[smooth.univariate[k]]]
      Xfrag <- Basis_generator(data[, basisinfo$term],basisinfo$N,
                               basisinfo$q,basisinfo$KnotsLocation,basisinfo$knots)
      Xb <- Xfrag$B
      last.para <-first.para+basisinfo$n.para
      X[, (first.para+1):last.para] <- as.matrix(Xb)
      first.para <-last.para
    }

    # Now have prediction matrix, X, for this block, need to do something with it...
    # Standardization for VS
    # NewX3 for SE
    # newind <-as.numeric(newind0-rep(block.size*(b-1),length(newind0)))
    newind <-newind0
    newind2 <-(start:stop)[newind]
    newX <-matrix(X[newind,,drop=FALSE],ncol=dim(X)[2])

    if (dim(newX)[2]>1){
      if (names(object$coefficients)[1]=='(Intercept)'){
        newX2 <-if(object$VS) {cbind(1,apply(newX[,-1], 2, scale))} else {newX}
        if (is.null(object$ind_c) | object$ind_c==0){
          newX3 <-matrix(0,ncol=1,nrow=length(newind))
        } else {
          newX3 <-newX2[,object$ind_c+1]
        }
      } else {
        # Bivariate term includes the intercept
        if (object$intercept){
          newX <-newX[,-1,drop=FALSE]
        }
        newX2 <-if(object$VS){apply(newX, 2, scale)} else {newX}
        newX3 <-newX2[,object$ind_c,drop=FALSE]
      }
    } else {
      newX2 <-matrix(0,ncol=1,nrow=length(newind))
      newX3 <-matrix(0,ncol=1,nrow=length(newind))
      object$coefficients <-0
    }
    #newind3 <-which(!is.na(newind2) && !is.na(newX2))
    newind30 <-newind2[!is.na(newind2)]
    newind3 <-newind30[apply(newX2, 1, function(x) !any(is.na(x))) ]
    newX2 <-newX2[apply(newX2, 1, function(x) !any(is.na(x))),,drop=FALSE]

    if (type=="terms") {      ## split results into terms
      ind_c=object$ind_c
      ### first deal with parametric terms
      #lass <- if (is.list(object$assign)) object$assign else list(object$assign)
      k <-0
        nptj <- max(length(ind_c))
        if (nptj>0) for (i in 1:length(ind_c)) { ## work through parametric part
          k <- k + 1 ## counts total number of parametric terms
          #fit <- fit[newind2,,drop=FALSE]
          fit[newind3,ind_c[k]] <- newX2[,ind_c[k],drop=FALSE]%*%object$coefficients[k+as.numeric(length(smooth.bivariate)==0)]

          ### need to check
          if (se.fit) {
            se[newind3,ind_c[k]] <-
              sqrt(pmax(0,rowSums((newX2[,k,drop=FALSE]%*%object$Ve[k,k])*newX2[,k,drop=FALSE])))
          }
        }
      #}
      ## assign list done, deal with univariate terms
      intecp=as.numeric(object$intercept && (length(smooth.bivariate)>0))

     # prediction for univariate terms
      if (n.smooth) {
        j=0
        if(n.smooth1>0){
          #first <- if (object$intercept) {dim(as.matrix(Xp[,-1]))[2]} else {dim(as.matrix(Xp))[2]}
          first=n.pterms
          for (j in 1:n.smooth1) # work through the smooth terms
          { #first <- object$basis_info[[j]]$first.para; last <- object$basis_info[[k]]$last.para
            n.para <- object$basis_info[[smooth.univariate[j]]]$n.para#-intecp
            last=first+n.para
            #fit <- fit[newind,,drop=FALSE]

            fit[newind3, n.pterms + j] <- newX2[, (first+1):last,
                                           drop = FALSE] %*%
              object$coefficients[(first+1):last] #+ Xoff[newind, j]
            first=last
          }
        }
        ### BPST Term
        #if (length(smooth.bivariate)>0) fit[newind,n.pterms+j+1] <- as.matrix(newbeta + Xoff[newind,j+1])
        if (length(smooth.bivariate)>0) fit[newind3,n.pterms+j+1] <- as.matrix(newbeta )
        colnames(fit) <- ColNames
        if (se.fit) colnames(se) <- ColNames
      } else {
        if (is.list(object$pterms)) {
          ## have to use term labels that match original data, or termplot fails
          ## to plot. This only applies for 'para.only' calls which are
          ## designed for use from termplot called from plot.gam
          term.labels <- unlist(lapply(object$pterms,attr,"term.labels"))
        }
        colnames(fit) <- term.labels
        if (se.fit) colnames(se) <- term.labels
        # retain only terms of order 1 - this is to make termplot work
        order <- if (is.list(object$pterms)) unlist(lapply(object$pterms,attr,"order")) else attr(object$pterms,"order")
        term.labels <- term.labels[order==1]
        ## fit <- as.matrix(as.matrix(fit)[,order==1])
        fit <- fit[,order==1,drop=FALSE]
        colnames(fit) <- term.labels
        if (se.fit) { ## se <- as.matrix(as.matrix(se)[,order==1])
          se <- se[,order==1,drop=FALSE]
          colnames(se) <- term.labels }
      }
    } else { ## "link" or "response" case
      fam <- object$family
      k <- attr(attr(object$model,"terms"),"offset")
      print(dim(newX2))
      print(length(object$coefficients))
      print(length(newbeta))
      fit[newind3] <- newX2%*%object$coefficients + as.matrix(newbeta)#+offs[newind]

      ### calculate se
        if (length(object$Ve)>0){
          if (se.fit) se[newind3] <- sqrt(pmax(0,rowSums(newX3%*%as.matrix(object$Ve))*newX3))
        }
        if (type=="response") { # transform
          linkinv <- fam$linkinv
          #if (is.null(fam$predict)) {
            dmu.deta <- fam$mu.eta
            if (se.fit) se[newind3]<-se[newind3]*abs(dmu.deta(fit))
            fit[newind3] <- linkinv(fit[newind3])
            # print(fit)
        }
    }
    }## end of link or response case
    rm(X)
  } ## end of prediction block loop


  if ((type == "terms") && (!is.null(terms) || !is.null(exclude))) {
    cnames <- colnames(fit)
    if (!is.null(terms)) {
      if (sum(!(terms %in% cnames)))
        warning("non-existent terms requested - ignoring")
      else {
        fit <- fit[, terms, drop = FALSE]
        if (se.fit) {
          se <- se[, terms, drop = FALSE]
        }
      }
    }
    if (!is.null(exclude)) {
      if (sum(!(exclude %in% cnames)))
        warning("non-existent exclude terms requested - ignoring")
      else {
        exclude <- which(cnames %in% exclude)
        fit <- fit[, -exclude, drop = FALSE]
        if (se.fit) {
          se <- se[, -exclude, drop = FALSE]
        }
      }
    }
  }

  if (type=="response") {
    #fit <- fit1
    #if (se.fit) se <- se1
    if (se.fit) se <-se
  }

  rn <- rownames(newdata)
  if (se.fit) {
    if (is.null(nrow(fit))) {
      names(fit) <- rn
      names(se) <- names(object$Ve)
      # Use missing value information to adjust residuals and predictions
      # fit <- napredict(na.act,fit)
      # se <- napredict(na.act,se)
    } else {
      rownames(fit)<-rn
      rownames(se)<-rn
      # fit <- napredict(na.act,fit)
      # se <- napredict(na.act,se)
    }
    H<-list(fit=fit,se.fit=se)
  } else {
    H <- fit
    # print(names(newdata))
    # print(fit)
    if (is.null(nrow(H))) names(H) <- rn else
      rownames(H)<-rn
    H <- napredict(na.act,H)
  }
  if ((type=="terms")&&attr(object$terms,"intercept")==1 && length(object$basis_info)==0)
  {attr(H,"constant") <- object$coefficients[1]} else{attr(H,"constant") <- NULL}
  H # ... and return
} ## end of predict.gam
