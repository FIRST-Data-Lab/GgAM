truncated_spline=function(x,knots){
  N=length(knots)-2
  XX=cbind(1,x,x^2,x^3)
  trun=function(x) max(x,0)
  for(i in 1:N){
    xx=sapply(x-knots[i+1],trun)^3
    XX=cbind(XX,xx)
  }

  XX1=cbind(0,1,2*x,3*x^2)
  for(i in 1:N){
    xx=3*sapply(x-knots[i+1],trun)^2
    XX1=cbind(XX1,xx)
  }

  XX2=cbind(0,0,2,6*x)
  for(i in 1:N){
    xx=6*sapply(x-knots[i+1],trun)
    XX2=cbind(XX2,xx)
  }
  list(XX=XX,XX1=XX1,XX2=XX2)
}
