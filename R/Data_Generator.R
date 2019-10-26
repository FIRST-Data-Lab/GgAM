#' Data Generating function
#' @param  n: Sample size.
#' @param p.x: Dimension of nonlinear components.
#' @param sig: Standard deviation of random error.
#' @param family: distribution family.
#' @param mu0: intercept
#' @param scale: scale parameter in gamma distribution
#' @param family: distribution family
#' @param m0: the value at grid points of bivariate function
#'
#' Output Arguments:
#'	(1) y: Generated response vector.
#'	(2) z: Generated design matrix of linear part; centralized for
#'		continuous z.
#'	(3) x: Generated design matrix of nonlinear part.
#' (4) fxTrue: Generated nonlinear functions values at generated x.
#' @export

Data_Generator=function(F1,F2,F3,n,p.x,mu0,sig,scale,
                        family,m0) {
  if (is.character(family))
    family <- eval(parse(text = family))
  if (is.function(family))
    family <- family()
  m=NA
  while (sum(is.na(m))>0){
    # Generate gridded bivariate function
    nmax=n*3
    fsb=list(fs.boundary())[[1]]
    names(fsb)=c("v","w")
    v=runif(nmax)*5-1
    w=runif(nmax)*2-1
    m=fs.test(v,w,b=1)
    ind=inSide(fsb,x=v,y=w) # remove outsiders
    v=v[ind]
    w=w[ind]
    m=m[ind]
    index=sample(1:sum(ind),size=n,replace=FALSE)
    s1=v[index]
    s2=w[index]
    s=cbind(s1,s2)
    m=m[index]
    m=m-mean(m)
    m=m/2
  }
  
  # Generate additive functions
  x=matrix(runif(n*p.x,-.5,.5),n,p.x)
  x=sweep(x,2,colMeans(x),"-")
  Fx=NULL
  for(l in 1:p.x) {Fx=cbind(Fx,do.call(paste0("F",l),list(x[,l])))}
  Fx=sweep(Fx,2,colMeans(Fx),"-")
  eta=as.numeric(mu0+rowSums(Fx)+m)
  
  # Generate response
  mu=family$linkinv(eta)
  if (family[[1]]=="gaussian") {y=rnorm(n,mu,sig)}
  if (family[[1]]=="binomial") {y=rbinom(n,1,mu)}
  if (family[[1]]=="poisson") {y=rpois(n,mu)}
  if(family[[1]]=="negative binomial") {y=rnbinom(n,mu=mu,size=5)}
  if (family[[1]]=='Gamma'){
    y=rgamma(n,shape=1/scale,scale=mu*scale)}
  list(y=y,eta=eta,mu=mu,x=x,fxTrue=Fx,s=s,psi=m,
       F1=F1,F2=F2,F3=F3,family=family)
}
