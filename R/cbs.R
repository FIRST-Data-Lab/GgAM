cbs=function(x,knots,Boundary.knots=range(x)){
  n=length(x)
  N=length(knots)
  Ix=matrix(0L,nrow=n,ncol=N+1L)
  for (i in 1:n) {
    if (x[i]<knots[1L]) {Ix[i,1]=1L
    } else if (x[i]>=knots[N]) {Ix[i,N+1L]=1L
    } else {
      kl=max(which(x[i]>=knots))
      Ix[i,kl+1]=1L
    }
  }
  return(Ix)
}
