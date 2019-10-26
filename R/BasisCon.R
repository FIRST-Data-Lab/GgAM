#' Basis Construction Function
#'
#' Smooth terms in a \code{plbpsm} formula are turned into smooth specification objects of class
#'  xx.smooth.spec during processing of the formula. Each of these objects is converted to
#'  a smooth object using an appropriate \code{BasisCon} function.
#'
#' @param object is a smooth specification object or a smooth object.
#' @param data A data frame, model frame or list containing the values of the (named) covariates
#' at which the smooth term is to be evaluated. If it’s a list then n must be supplied.
#'
#' @details It is the wrapper function which calls basis constructing method.
#' @return
#' a list of \code{smooth} objects are returned. Many of the outputs are from \code{b} function.
#' Other outputs include all the information related to berstein basis.
# \item{B}{The calculated basis matrix}
# \item{K}{The calculated energy matrix.}
# \item{H}{The calculated constraint matrix.}
# \item{Ind}{The index of points inside the triangulation.}
# \item{fx}{\code{TRUE} if the term is to be unpenalized, otherwise \code{FALSE}.}
#'
#' @examples
#' library(BPST)
#' data("eg1pop_dat")
#'eg1_V2=eg1pop_dat[['V2']]
#'eg1_T2=eg1pop_dat[['T2']]
#'eg1pop_rho03=eg1pop_dat[['rho03']]
#'sam=eg1pop_rho03[sample(1:dim(eg1pop_rho03)[1],70),]
#'
#'# bivariate spline
#'BI <- BasisCon(b(x1,x2,d=2,r=1,V=eg1_V2,Tr=eg1_T2),sam)
#'# univariate spline
#'BI <- BasisCon(u(z1),sam)
#' @export
######################################################
BasisCon <- function(object,data)
  UseMethod("BasisCon")
  ## wrapper function which calls basis constructing method
  ## Handles `by' variables, and summation convention.
  ## Note that `data' must be a data.frame or model.frame, unless n is provided explicitly,
  ## in which case a list will do.
  ## If present dataX specifies the data to be used to set up the model matrix, given the
  ## basis set up using data (but n same for both).
