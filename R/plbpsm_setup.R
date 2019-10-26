#' @importFrom stats model.frame model.offset model.matrix .getXlevels
plbpsm.setup <- function(formula,pterms,
                       data=stop("No data supplied to plbpsm.setup"),
                       drop.intercept=FALSE,group)
  ## set up the model matrix, penalty matrices and auxilliary information about the Bersterin bases

{ # split the formula if the object being passed is a formula, otherwise it's already split

  if (inherits(formula,"split.plbpsm.formula")) split <- formula else
    if (inherits(formula,"formula")) split <- interpret.plbpsm0(formula) else stop("First argument is no sort of formula!")

    if (length(split$smooth.spec)==0) {
      if (split$pfok==0) stop("You've got no model....")
      m <- 0
    } else  m <- length(split$smooth.spec) # number of spline terms
    G <- list(m=m,pterms=pterms) ##
    if (is.null(attr(data,"terms"))) # then data is not a model frame
      mf <- model.frame(split$pf,data,drop.unused.levels=FALSE) else mf <- data # data is already a model frame
    G$intercept <-  attr(attr(mf,"terms"),"intercept")>0
    G$offset <- model.offset(mf)   # get the model offset (if any)
    if (!is.null(G$offset))  G$offset <- as.numeric(G$offset)

    # construct parametric model matrix
    if (drop.intercept) attr(pterms,"intercept") <- 1 ## ensure there is an intercept to drop
    X <- model.matrix(pterms,mf)
    G$Xp <- X
    # if (drop.intercept) { ## some extended families require intercept to be dropped
    #   xat <- attributes(X);ind <- xat$assign>0
    #   X <- X[,xat$assign>0,drop=FALSE] ## some extended families need to drop intercept
    #   xat$assign <- xat$assign[ind];xat$dimnames[[2]]<-xat$dimnames[[2]][ind];
    #   xat$dim[2] <- xat$dim[2]-1;attributes(X) <- xat
    #   G$intercept <- FALSE
    # }
    rownames(X) <- NULL ## save memory
    G$mf <- mf
    G$nsdf <- ncol(X)
    G$contrasts <- attr(X,"contrasts")
    G$xlevels <- .getXlevels(pterms,mf)
    G$assign <- attr(X,"assign") # used to tell which coeffs relate to which pterms
    ##
    term.smooth <- sapply(split$smooth.spec, class)
    smooth.bivariate <- which(term.smooth=="bivariate.smooth")
    smooth.univariate <- which(term.smooth=="univariate.smooth")
    m1 <- length(smooth.univariate)
    m2 <- length(smooth.bivariate)
    m_b <- 0
    if (m2>0){
      m_b <- m_b+1
      if (!is.null(split$smooth.spec[[smooth.bivariate]]$r))
        if (split$smooth.spec[[smooth.bivariate]]$r > split$smooth.spec[[smooth.bivariate]]$d) {
          stop ('r should be smaller than or equal to d')}
      G$fx <- split$smooth.spec[[smooth.bivariate]]$fixed
      G$lambda <- split$smooth.spec[[smooth.bivariate]]$lambda
      G$basis_info[[smooth.bivariate]] <- BasisCon(split$smooth.spec[[smooth.bivariate]],data)[[1]]
      data <- data[G$basis_info[[smooth.bivariate]]$ind,]
      X_B <- G$basis_info[[smooth.bivariate]]$B
      if (!is.null(dim(X)) && !is.null(dim(G$Xp))){
        X <- X[G$basis_info[[smooth.bivariate]]$ind,]
        G$Xp <- G$Xp[G$basis_info[[smooth.bivariate]]$ind,]
      } else{
        X <- X[G$basis_info[[smooth.bivariate]]$ind]
        G$Xp <- G$Xp[G$basis_info[[smooth.bivariate]]$ind]
      }

    } else {X_B <- NULL}

    if (m1>0){
      if (is.null(group)){
        vars <- sapply(split$smooth.spec[smooth.univariate],function(X){X$term})
        #group2 <- seq(1,length(vars))
        group2 <- as.factor(vars)
        group2 <- as.numeric(group2)
        #WHY?
        group2 <- sort(group2)
        #MI <- any(group %in% vars)
      } else {group2 <- group}
    } else{
      group2 <- group0 <- 0
    }

    split2 <- split # for constructing constant splines for model identification
    if (G$nsdf > 0) term.names <- colnames(X)[1:G$nsdf] else term.names<-array("",0)
    first.para<-G$nsdf+1
    if (m1>0) {
      Xorig <- c()
      group0 <- c()
      for (i in 1:m1){
          if(i <= length(group2) & group2[i]>0){
          split2$smooth.spec[[i]]$q <- 0
          split2$smooth.spec[[i]]$N <- split2$smooth.spec[[i]]$N_MI
          G$basis_info_MI[[i]] <- BasisCon(split2$smooth.spec[[i]],data)[[1]]
          Xorig2 <- data[,G$basis_info_MI[[i]]$term]
          Xorig <- cbind(Xorig,Xorig2)
          tem2 <- rep(group2[i],dim(G$basis_info_MI[[i]]$B)[2])
          group0 <- c(group0,tem2)
          }
        # no MI
        G$basis_info[[i]] <- BasisCon(split$smooth.spec[[i]],data)[[1]]
        if (!is.null(group) & i <= length(group2) & group2[i]==0){
          tem2 <- rep(group2[i],dim(G$basis_info[[i]]$B)[2])
          group0 <- c(group0,tem2)
        }

        n.para<-G$basis_info[[i]]$n.para
        # define which elements in the parameter vector this smooth relates to....
        G$basis_info[[i]]$first.para<-first.para
        first.para<-first.para+n.para
        G$basis_info[[i]]$last.para<-first.para-1
        if (class(split$smooth.spec[[i]])=='univariate.smooth'){
          X <- cbind(X,G$basis_info[[i]]$B)#;basis_info[[i]]$B<-NULL
          k<-1
          jj <- G$basis_info[[i]]$first.para:G$basis_info[[i]]$last.para
          for (j in jj) {
            term.names[j] <- paste(G$basis_info[[i]]$label,".",as.character(k),sep="")
            k <- k+1
          }
        }
      }
    } #else { warning('no smoothing functions in the model.')
    #   }
    G$m_b <- m_b
    G$group0 <- group0
    if (length(smooth.univariate)>0){
        if (length(Xorig)>0){
      XX <- sweep(Xorig,2,colMeans(Xorig,na.rm=T),"-")
      ## data/data_ecdf
      G$XX <- XX
      G$Xorig <- data[, sapply(G$basis_info_MI,function(X){X$term})]
    }
  }

    if (G$nsdf>0) G$cmX[-(1:G$nsdf)] <- 0 else G$cmX <- G$cmX * 0

    G$smooth.univariate <- smooth.univariate
    G$X <- X;rm(X)
    n.p <- ncol(G$X)
    G$y <- data[[split$response]]
    G$n <- nrow(data)

    if (is.null(data$"(weights)")) G$w <- rep(1,G$n) else G$w <- data$"(weights)"

    ## Create names for model coefficients...
    G$term.names <- term.names
    if (is.null(G$lambda)){
      G$lambda <- 10^(seq(-2,5,by <- 1))
    }
    return(G)
} ## plbpsm.setup
