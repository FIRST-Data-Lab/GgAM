#' Default PLBPSM plotting
#'
#' Takes a fitted \code{plbpsm} object produced by \code{plbpsm()} and plots the triangulation of
#' location data points, predicted surface of bivariate smooth function and optionally produces
#' histogram of residuals for the model.
#' @importFrom grDevices devAskNewPage
#' @importFrom graphics par plot hist
#'@param x a fitted \code{plbpsm} object as produced by \code{plbpsm}
#'@param residuals If \code{TRUE} then a histogram of standardized residuals will be added.
#'@param main Used as title for plots if supplied.
#'@param xlim If supplied then this pair of numbers are used as the x limits for each plot.
#'@param ylim If supplied then this pair of numbers are used as the y limits for each plot.
#'@param xlab If supplied then this will be used as the x label for all plots.
#'@param ylab If supplied then this will be used as the y label for all plots.
#'@param pages (default 0) the number of pages over which to spread the output. For example,
#'if \code{pages=1} then all terms will be plotted on one page with the layout performed
#'automatically. Set to 0 to have the routine leave all graphics settings as they are.
#'@param select Allows the plot for a single model term to be selected for printing. e.g. if you
#'just want the plot for the second smooth term set \code{select=2}.
#'@param n1 number of points used in x axis in each plot.
#'@param n2 number of points used in y axis in each plot.
#'@param ... other graphics parameters to pass on to plotting commands. See details for
#'smooth plot specific options.
#'@details Used R package \code{fdaPDE} and \code{plotly} to draw triangulation plot and predicted surfaces. See
#'\code{plbpsm:::plot.plbpsm.smooth}.
#' @examples
#' library(MASS)
#' library(grpreg)
#'
#' # irregular domain:
#' library(GgAM)
#' library(BPST)
#' data("eg1pop_dat")
#' eg1_V2=eg1pop_dat[['V2']]
#' eg1_T2=eg1pop_dat[['T2']]
#' eg1pop_rho03=eg1pop_dat[['rho03']]
#' n=1000
#' Npop=nrow(eg1pop_rho03)
#' # ind.pop=(1:Npop)[!is.na(eg1pop_rho03[,1])]
#' ind.pop=(1:Npop)
#' sam.ind=sort(sample(ind.pop,n,replace=FALSE))
#' sam=eg1pop_rho03[sam.ind,]
#' lambda=10^(seq(-2,5,by=1))
#' data=sam
#' formula=Y~z1+z2+z3+z4+z5+z6+z7+z8+b(x1,x2,V=eg1_V2,Tr=eg1_T2,d=2,r=1,lambda=lambda)
#' res=plbpsm(formula=formula,data=as.data.frame(data))
#' plot(res,residuals=TRUE,n1=80,n2=50)
#'
#' ### GGAM ###
#' data(dat_poi_ggams)
#' n=100
#' Npop=nrow(dat_poi_ggams)
#' # ind.pop=(1:Npop)[!is.na(eg1pop_poi2[,1])]
#' ind.pop=(1:Npop)
#' sam.ind=sort(sample(ind.pop,n,replace=FALSE))
#' sam=dat_poi_ggams[sam.ind,]
#' data=sam
#' formula=y~u(x1)+u(x2)+u(x3)+b(s1,s2,V=eg1_V2,Tr=eg1_T2,d=2,r=1)
#' res_eg1_poi_add=plbpsm(formula=formula,data=as.data.frame(data),family='poisson')
#' summary(res_eg1_poi_add)
#' plot(res_eg1_poi_add)
#'@export

plot.plbpsm=function(x, residuals = FALSE, pages = 0, select = NULL, xlab = NULL, ylab = NULL,
                     main = NULL, ylim = NULL, xlim = NULL, n1=40, n2=40,...) {
  term.smooth=sapply(x$basis_info, class)[1,]
  smooth.bivariate=which(term.smooth=="bivariate.smooth")
  
  if (x$MI){
    #term.smooth=sapply(x$basis_info2, class)[1,]
    smooth.univariate=x$ind.nl
    
  } else {
    smooth.univariate=which(term.smooth=="univariate.smooth")
  }
  m=length(term.smooth)
  m1=length(smooth.univariate)
  m2=length(smooth.bivariate)
  n.plots=0
  obj=summary(x,...)
  
  if (m > 0){
    ### Then plot the univariate smooth related plots, later, we may need to plot
    ### for several smoothing terms.
    if (x$backfitting){
      if (m1>0){ ### confidence band
        for (i in 1:m1) if ((is.null(select) ||i == select)) {
          plot(obj$bands[[i]], #xlab=NULL,ylab=NULL,main=NULL,
               xlim=NULL, ylim=NULL,...)
          # if (ask) {
          #   oask <- devAskNewPage(TRUE)
          #   on.exit(devAskNewPage(oask))
          #   ask <- FALSE
          oask<-devAskNewPage(ask = TRUE)
          on.exit(devAskNewPage(oask))
          # print(basisplot$p.beta)
          # oask<-devAskNewPage(ask = TRUE)
          # on.exit(devAskNewPage(oask))
          # plot(basisplot$triplot,main='Triangulation Plot')
        }
        #}
      }
    }
    if (m2>0){
      pd <- list()
      i <- 1
      for (i in 1:m2) {
        x$basis_info[[smooth.bivariate[i]]]$coefficients <- x$coefficients_bivariate
        P <- plot(x$basis_info[[smooth.bivariate[i]]], P = NULL, data = x$model,
                  n1 = n1, n2 = n2, xlab = xlab,
                  ylab = ylab, main = main, ylim = ylim,
                  xlim = xlim,...)
        pd[[i]] <- P
        rm(P)
      }
      ### Then plot the bivariate smooth related plots, later, we may need to plot
      ### for several smoothing terms.
      for (i in 1:m2) if ((is.null(select) ||i == select)) {
        basisplot=plot(x$basis_info[[smooth.bivariate[i]]], P = pd[[i]],xlab = xlab,
                       ylab = ylab, main = main, ylim = ylim,
                       xlim = xlim,n1=n1,n2=n2,...)
        # if (ask) {
        #   oask <- devAskNewPage(TRUE)
        #   on.exit(devAskNewPage(oask))
        #   ask <- FALSE
        oask<-devAskNewPage(ask = TRUE)
        on.exit(devAskNewPage(oask))
        print(basisplot$p.beta)
        # oask<-devAskNewPage(ask = TRUE)
        # on.exit(devAskNewPage(oask))
        #plot(basisplot$triplot,main='Triangulation Plot')
      }
      devAskNewPage(ask = FALSE)
    } else {
      basisplot=NULL}
  } else {
    basisplot=NULL
    warning('No smooth term, no confindence bands, predicted surface plot nor triangulation plot.')
  }
  if (m > 0)
    #for (i in 1:m)
    n.plots <- n.plots + m #as.numeric(pd[[i]]$plot.me)
  
  ### residue plots
  if (residuals==TRUE){
    n.plots=n.plots+1
    stdres=scale(x$residuals,...)
    oask<-devAskNewPage(ask = TRUE)
    on.exit(devAskNewPage(oask))
    hist(stdres,main='Histogram of Residuals')
  } else {stdres=NULL}
  
  
  if (n.plots == 0)
    stop("No terms to plot - nothing for plot.plbpsm() to do.")
  
  # We may need this later, since right now, plotly and other regular plots are not in the same ploting system,
  # we cannot plot them together.
  
  # if (pages > n.plots)
  #   pages <- n.plots
  # if (pages < 0)
  #   pages <- 0
  # if (pages != 0) {
  #   ppp <- n.plots%/%pages
  #   if (n.plots%%pages != 0) {
  #     ppp <- ppp + 1
  #     while (ppp * (pages - 1) >= n.plots) pages <- pages -
  #         1
  #   }
  #   c <- r <- trunc(sqrt(ppp))
  #   if (c < 1)
  #     r <- c <- 1
  #   if (c * r < ppp)
  #     c <- c + 1
  #   if (c * r < ppp)
  #     r <- r + 1
  #   oldpar <- par(mfrow = c(r, c))
  # } else {
  #   ppp <- 1
  #   oldpar <- par()
  # }
  # if ((pages == 0 && prod(par("mfcol")) < n.plots && dev.interactive()) ||
  #     pages > 1 && dev.interactive())
  #   ask <- TRUE else ask <- FALSE
  
  
  
  ### set up pages related to number of plots
  # if (pages > 0) {
  #   par(oldpar)
  # }
  allinfo=list(basisplot=basisplot,stdres=stdres)
  invisible(allinfo)
}
