#' Setting PLBPSM fitting defaults
#'
#' This is an internal function of package \code{ggam} which allows control of the numerical options for
#' fitting a PLBPSM. Typically users will want to modify the defaults if model fitting fails to converge, or
#' if the warnings are generated which suggest a loss of numerical stability during fitting.
#' @param delta1 The convergence criterion in \code{\link{gplsfitGCV}}.
#' @param delta2 The convergence criterion in \code{\link{gplsfitGCV}}.
#' @param trace Set this to \code{TRUE} to turn on diagnostic output.
#' @param maxstep Maximum number of iterations to perform.
#' @param epsilon This is used for judging conversion of the loop in \code{\link{gplsfitGCV}}.
#' @export
plbpsm.control <- function(delta1=1,delta2=1,trace=FALSE,maxstep=10,epsilon=1e-7){
  list(trace=trace,maxstep=maxstep,epsilon=epsilon,delta1=delta1,delta2=delta2)
}
