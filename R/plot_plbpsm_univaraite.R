#' @importFrom graphics rug lines polygon abline
### The plots information for the univariate smooth term.
plot.band.univariate=function(x,xlab_u=NULL,main_u=NULL,xlim,ylim,shade=FALSE,shade.col="gray80",ab.line=FALSE,add.rug=FALSE,...){

  # if (x$dim==1) {
  # if (is.null(xlab)) xlabel <- xterm else xlabel <- xlab
  # yterm <- x$term[2]
  # if (is.null(ylab)) ylabel <- yterm else ylabel <- ylab
  if (is.null(ylim)) #ylim=c(-1.5,1.5)
  ylim <- c(range(x$lower)[1],range(x$upper)[2])
  if (is.null(xlim)) xlim <- range(x$x0)

  if (is.null(xlab_u)) xlab_u <- x$term
  if (is.null(main_u)) {
    main_u <- x$label
  }
  ## ?
  n <- 10
  # if (shade) {
    plot(x$x0,x$est,type='l',ylim=ylim,xlab = xlab_u,ylab='',main = main_u,xlim=xlim)
    #,
    polygon(c(x$x0, rev(x$x0)), c(x$upper, rev(x$lower)),
            col = "grey", border = NA)
    if(ab.line){
      abline(h=0,lty=2)
    }
    lines(x$x0,x$est,type='l',...)
    if (add.rug){
      rug(x$dat[,x$term], lwd = 1)
    }

    ####
    # plot(x$x0,x$est,type="n",xlab='',ylab = '',ylim=ylim,xlim=xlim,...)
    # polygon(c(x$x0,x$x0[n:1],x$x0[1]),
    #         c(x$upper,x$lower[n:1],x$upper[1]),col = shade.col,border = NA)
    # lines(P$x0,P$est,...)
  # } else {
  # plot(x$x0,x$est,type="l",ylab = '',ylim=ylim,xlim=xlim,xlab = xlab_u, main = main_u)
  # if (is.null(list(...)[["lty"]])) {
  #   lines(x$x0,x$lower,lty=2)
  #   lines(x$x0,x$upper,lty=2)
  #   } else {
  #   lines(x$x0,x$lower,...)
  #   lines(x$x0,x$upper,...)
  #   }
  # }
  ## end of plot production
  return(list(xlim=xlim,ylim=ylim))
  # }
  # else {
  #   return(NULL)
  # }
}
