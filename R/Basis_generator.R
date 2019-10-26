#' @importFrom splines bs
Basis_generator=function(x,N,q,KnotsLocation,knots){
  # Case I: x is a vector
  if (length(dim(x))==0) {
    # generate knots by supplied location
    if (is.null(knots)) {
      # knots generated in sample quantile
      if (KnotsLocation=="quantile")
        knots1=unique(quantile(x,probs=seq(0,1,by=1/(N+1))))
      # knots generated in uniform points
      if (KnotsLocation=="uniform")
        knots1=seq(min(x),max(x),by=(max(x)-min(x))/(N+1))
    } else knots1=knots
    N=length(knots1)-2
    # generate constant / polynomial B-spline basis
    if (q==0) {
      Bx0=cbs(x,knots1[-c(1,N+2)],Boundary.knots=knots1[c(1,N+2)])
    } else Bx0=bs(x,knots=knots1[-c(1,N+2)],degree=q,intercept=F,
                  Boundary.knots=knots1[c(1,N+2)])
    Knots=knots1
  } else { # Case II: x is a matrix
    d.x=ncol(x)
    block.size=floor(sqrt(d.x))
    nblock=ceiling(d.x/block.size)
    Bx0=NULL
    Knots=NULL

    # This routine reduces the memory needed in constructing spline
    #	bases and improve the computational efficiency
    for(nj in 1:nblock) {
      Bxnj=NULL
      Knotsj=NULL
      if(nj<nblock) block.ind=(block.size*(nj-1)+1):(block.size*nj)
      if(nj==nblock) block.ind=(block.size*(nj-1)+1):d.x
      for (j in block.ind) {
        # generate knots by supplied location
        if (is.null(knots)) {
          # knots generated in sample quantile
          if (KnotsLocation=="quantile") knots1=quantile(
            x[,j],probs=seq(0,1,by=1/(N+1)))
          # knots generated in uniform points
          if (KnotsLocation=="uniform") knots1=seq(min(x[,j]),
                                                   max(x[,j]),by=(max(x[,j])-min(x[,j]))/(N+1))
        } else knots1=knots[(N+2)*(j-1)+1:(N+2)]

        # generate constant / polynomial B-spline basis
        if (q==0) {
          bx0=cbs(x[,j],knots1[-c(1,N+2)],
                  Boundary.knots=knots1[c(1,N+2)])
        } else bx0=bs(x[,j],knots=knots1[-c(1,N+2)],degree=q,
                      intercept=F,Boundary.knots=knots1[c(1,N+2)])
        Bxnj=cbind(Bxnj,bx0)
        Knotsj=cbind(Knotsj,knots1)
      }
      Bx0=cbind(Bx0,Bxnj)
      Knots=cbind(Knots,Knotsj)
    }
  }
  BxMean=colMeans(Bx0)
  Bx=sweep(Bx0,2,BxMean,"-")
  list(Bx0=Bx0,B=Bx,BxMean=BxMean,prior.knots=knots,knots=Knots)
}
