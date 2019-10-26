#' Defining smooths in PLBPSM formulae
#'
#'Function used in definition of smooth terms within \code{plbpsm} model formulae. The function does not
#'evaluate a (spline) smooth - it exists purely to help set up a model using spline based smooths.
#' @importFrom  stats terms reformulate
#' @param ... a list of variables that are the covariates that this
#' smooth is a function of.
#' @param d degree of polynomials.
#' @param r smoothness and r \eqn{\leq}{\leq} d
#' @param V an \code{N} by two matrix that lists vertices with the \code{i}th
#' row storing in Cartesian coordinates for the ith vertex. \code{N} is the number of vertices.
#' @param Tr a \code{K} by three matrix that each row represents one triangle.
#' All the elements are the integers that stand for the indices of vertices in \code{V}. \code{K} is the
#' number of triangles.
#' @param b Boundary of the domain of sample points.
#' @param nt A parameter controls the number of total triangles.
#' @param Holes Information of holes of polygon.
#' @param ind An ordering indices of observation points, in which the
#' cnt[j]+1th to cnt[j+1]th elements are indices of points in the jth triangle.
#' @param B Bernstein basis matrix
#' @param Q2 The \code{Q} matrix from QR decomposition of the constraint matrix.
#' @param K Energy matrix for constructing penalty matrix
#' @param lambda The default set of smoothing penalty parameter to be chosen from
#' @param fx indicates whether the term is a fixed d.f. regression
#'   spline (\code{TRUE}) or a penalized regression spline (\code{FALSE}).
#' @param id An identifying label or number for the smooth, linking it to other smooths. Defaults to \code{NULL} for no linkage.
#' @return
#' These \code{smooth.spec} objects define bivariate smooths and are turned into bases and penalties by
#' \code{BasisCon} functions.
#'  The returned object contains the following items:
#' \item{d}{degree of polynomials.}
#' \item{r}{smoothness and \code{r}\eqn{\leq}{\leq} \code{d}.}
#' \item{V}{an \code{N} by two matrix that lists vertices with the \code{i}th
#' row storing in Cartesian coordinates for the ith vertex. \code{N} is the number of vertices.}
#' \item{Tr}{a \code{K} by three matrix that each row represents one triangle.
#' All the elements are the integers that stand for the indices of vertices in \code{V}. \code{K} is the
#' number of triangles.}
#' \item{ind}{An ordering indices of observation points corresponding to index of triangles.}
#' \item{B}{Bernstein basis matrix}
#' \item{K}{Energy matrix for constructing penalty matrix}
#' \item{lambda}{The default set of smoothing penalty parameter to be chosen from}
#' \item{term}{An array of text strings giving the names of the covariates that the term is a function of.}
#' \item{fixed}{\code{TRUE} if the term is to be treated as a pure regression
#' spline (with fixed degrees of freedom); FALSE if it is to be treated as a penalized regression spline}
#' \item{dim}{The dimension of the smoother - i.e. the number of covariates that it is a function of.}
#' \item{label}{A suitable text label for this smooth term.}
#' \item{id}{An identifying label or number for the smooth, linking it to other smooths. Defaults to \code{NULL} for no linkage. }
#' @examples
#' library(GgAM)
#' data(eg1pop_dat)
#'V=eg1pop_dat[['V2']]
#'Tr=eg1pop_dat[['T2']]
#' res <- b(x1,x2,d=2,r=1,V=V,Tr=Tr,lambda=0)
#' @export

b <- function (..., d = NULL, r = NULL, V = NULL, Tr = NULL, b = NULL, nt=NULL , Holes = NULL,
               B = NULL, Q2 = NULL, K = NULL, ind = NULL, lambda = NULL,
               fx = FALSE, id = NULL)
{
  vars <- as.list(substitute(list(...)))[-1]
  dim <- length(vars)
  term <- deparse(vars[[1]], backtick = TRUE, width.cutoff = 500)
  if (term[1] == ".")
    stop("f(.) not supported.")
  if (dim > 1)
    for (i in 2:dim) {
      term[i] <- deparse(vars[[i]], backtick = TRUE, width.cutoff = 500)
      if (term[i] == ".")
        stop("f(.) not yet supported.")
    }
  for (i in 1:dim) term[i] <- attr(terms(reformulate(term[i])),
                                   "term.labels")
  if (length(unique(term)) != dim)
    stop("Repeated variables as arguments of a smooth are not permitted")
  full.call <- paste("b(", term[1], sep = "")
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
  ret <- list(term = term, d = d, r = r, V = V, Tr = Tr, b = b,nt=nt, Holes=Holes, B=B,Q2=Q2,K=K,ind=ind,lambda = lambda , fixed = fx, dim = dim,
              label = label, id = id)
  # if (!is.null(pc)) {
  #   if (length(pc) < d)
  #     stop("supply a value for each variable for a point constraint")
  #   if (!is.list(pc))
  #     pc <- as.list(pc)
  #   if (is.null(names(pc)))
  #     names(pc) <- unlist(lapply(vars, all.vars))
  #   ret$point.con <- pc
  # }
  class(ret) <- "bivariate.smooth"
  ret
}

