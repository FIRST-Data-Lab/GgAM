#' Defining smooths in PLBPSM formulae
#'
#'Function used in definition of smooth terms within \code{plbpsm} model formulae. The function does not
#'evaluate a (spline) smooth - it exists purely to help set up a model using spline based smooths.
#' @importFrom  stats terms reformulate
#' @param ... a list of variables that are the covariates that this
#' smooth is a function of.
#' @param N Number of interior knots in generating spline matrix.
#' @param q Degree of polynomial spline. Default is 3.
#' @param KnotsLocation A character string naming the way for knots
#'			locations. Default is "quantile". The only alternative is
#'			"uniform".
#' @param knots An optional vector specifying the knots to be used in
#'			constructing spline bases.
#' @param N_MI Number of interior knots in generating spline matrix in the model identification process.
#' @param fx indicates whether the term is a fixed d.f. regression
#'   spline (\code{TRUE}) or a penalized regression spline (\code{FALSE}).
#' @param id An identifying label or number for the smooth, linking it to other smooths. Defaults to \code{NULL} for no linkage.
#' @return
#' These \code{smooth.spec} objects define uivariate smooths and are turned into bases and penalties by
#' \code{BasisCon} functions.
#'  The returned object contains the following items:
#' \item{term}{An array of text strings giving the names of the covariates that the term is a function of.}
#' \item{N}{Number of interior knots in generating spline matrix.}
#' \item{q}{Degree of polynomial spline. Default is 3.}
#' \item{knotsLocation}{A character string naming the way for knots
#'			locations. Default is "quantile". The only alternative is
#'		"uniform".}
#' \item{knots}{An optional vector specifying the knots to be used in
#'			constructing spline bases.}
#' \item{N_MI}{Number of interior knots in generating spline matrix in the model identification process.}
#' \item{fx}{\code{TRUE} if the term is to be treated as a pure regression
#' spline (with fixed degrees of freedom); FALSE if it is to be treated as a penalized regression spline}
#' \item{dim}{The dimension of the smoother - i.e. the number of covariates that it is a function of.}
#' \item{label}{A suitable text label for this smooth term.}
#' \item{id}{An identifying label or number for the smooth, linking it to other smooths. Defaults to \code{NULL} for no linkage. }
#' @examples
#' library(GgAM)
#' data(eg1pop_dat)
#'V=eg1pop_dat[['V2']]
#'Tr=eg1pop_dat[['T2']]
#' res <- u(x1,x2,N=2,q=3)
#' @export

u <- function (..., N=2,q=3,KnotsLocation='quantile',knots=NULL,N_MI=4,
               fx = FALSE, id = NULL)
{
  vars <- as.list(substitute(list(...)))[-1]
  dim <- length(vars)

  term <- deparse(vars[[1]], backtick = TRUE, width.cutoff = 500)
  if (term[1] == ".")
    stop("u(.) not supported.")
  if (dim > 1)
    for (i in 2:dim) {
      term[i] <- deparse(vars[[i]], backtick = TRUE, width.cutoff = 500)
      if (term[i] == ".")
        stop("f(.) not yet supported.")
    }
  for (i in 1:dim) term[i] <- attr(terms(reformulate(term[i])),
                                   "term.labels")
  d.new <- round(q)
  if (all.equal(d.new, q) != TRUE) {
    warning("argument q of u() should be integer and has been rounded")
  }
  d <- d.new
  if (length(unique(term)) != dim)
    stop("Repeated variables as arguments of a smooth are not permitted")
  full.call <- paste("u(", term[1], sep = "")
  if (dim > 1)
    for (i in 2:dim) full.call <- paste(full.call, ",", term[i],
                                        sep = "")
  label <- paste(full.call, ")", sep = "")
  if (!is.null(id)) {
    if (length(id) > 1) {
      id <- id[1]
      warning("only first element of `id' used")
    }
    id <- as.character(id)
  }
  ret <- list(term = term, N=N, q = d, KnotsLocation=KnotsLocation,knots=knots,N_MI=N_MI,fixed = fx, dim = dim,
               label = label, id = id)
  class(ret) <- "univariate.smooth"
  ret
}

