#' Basis Construction Function
#'
#' This is an internal function of package \code{ggam}.
#' @importFrom BPST basis qrH
#' @importFrom Matrix Matrix
#' @importFrom mgcv get.var
#' @param object is a smooth specification object or a smooth object.
# \item{term}{An array of text strings giving the names of the covariates that the term is a function of.}
# \item{N}{Number of interior knots in generating spline matrix.}
# \item{q}{Degree of polynomial spline. Default is 3.}
# \item{knotsLocation}{A character string naming the way for knots
# 		locations. Default is "quantile". The alternative is
# 	"uniform".}
# \item{knots}{An optional vector specifying the knots to be used in
# 		constructing spline bases.}
# \item{N_MI}{Number of interior knots in generating spline matrix in the model identification process.}
# \item{fx}{\code{TRUE} if the term is to be treated as a pure regression
# spline (with fixed degrees of freedom); FALSE if it is to be treated as a penalized regression spline}
# \item{dim}{The dimension of the smoother - i.e. the number of covariates that it is a function of.}
# \item{label}{A suitable text label for this smooth term.}
# \item{id}{An identifying label or number for the smooth, linking it to other smooths. Defaults to \code{NULL} for no linkage. }
#' @param data A data frame, model frame or list containing the values of the (named) covariates
#' at which the smooth term is to be evaluated. If it’s a list then \code{n} must be supplied.
#' @details It is the wrapper function which calls basis constructing method.
#' @return
#' a list of \code{smooth} objects are returned. Many of the outputs are from \code{f} function.
#' Other outputs include all the information related to berstein basis.
# \item{B}{The calculated basis matrix}
# \item{K}{The calculated energy matrix.}
# \item{H}{The calculated constraint matrix.}
# \item{Ind}{The index of points inside the triangulation.}
# \item{fx}{\code{TRUE} if the term is to be unpenalized, otherwise \code{FALSE}.}
#'
#' @examples
#' library(BPST)
#' library(GgAM)
#' data("eg1pop_dat")
#'eg1_V2=eg1pop_dat[['V2']]
#'eg1_T2=eg1pop_dat[['T2']]
#'eg1pop_rho03=eg1pop_dat[['rho03']]
#'sam=eg1pop_rho03[sample(1:dim(eg1pop_rho03)[1],70),]
#'aa=BasisCon(u(z1),sam)
#' @export
######################################################
BasisCon.univariate.smooth <- function(object,data)
  ## wrapper function which calls basis constructing method
  ## Handles `by' variables, and summation convention.
  ## Note that `data' must be a data.frame or model.frame, unless n is provided explicitly,
  ## in which case a list will do.
  ## If present dataX specifies the data to be used to set up the model matrix, given the
  ## basis set up using data (but n same for both).
{
  ###
  # x=data[, object$term]
  # N=object$N
  # q=object$q
  object2 <- Basis_generator(data[, object$term],object$N,object$q,object$KnotsLocation,object$knots)
  object$Bx0=object2$Bx0
  object$B=object2$B
  object$BxMean=object2$BxMean
  object$knots=object2$knots
  object$prior.knots=object2$prior.knots
  # null_sp_dim <- getFromNamespace("null.space.dimension", "mgcv")
  #
  # object$null.space.dimension=null_sp_dim(object$dim,0)
  ### should be the X outside since BasisCon does not deal with para part?

  #object$datanew=data[object$ind,]

  ## add plotting indicator if not present.
  ## plot.me tells `plot.plbpsm' whether or not to plot the term
  plot.me <- TRUE
  object$n.para=ncol(object$B)
  # may change in univariate case

  #offs <- NULL
  ## pick up "by variables" now, and handle summation convention ...
  # if (object$by!="NA") {
  #   if (is.null(dataX)) by <- get.var(object$by,datanew) else by <- get.var(object$by,dataX)
  #   if (is.null(by))
  #     if (is.null(object$IND)) { ## then by to be taken as sequence of 1s
  #       by <- rep(1,nrow(object$datanew))} else by <- object$IND
  #       if (is.null(by)) stop("Can't find by variable")
  #       #offs <- attr(sm$X,"offset")
  #       if (is.factor(by)) { ## generates smooth for each level of by
  #         objectl <- list()
  #         lev <- levels(by)
  #         ## if by variable is an ordered factor then first level is taken as a
  #         ## reference level, and smooths are only generated for the other levels
  #         ## this can help to ensure identifiability in complex models.
  #         if (is.ordered(by)&&length(lev)>1) lev <- lev[-1]
  #         for (j in 1:length(lev)) {
  #           objectl[[j]] <- object  ## replicate smooth for each factor level
  #           by.dum <- as.numeric(lev[j]==by)
  #           objectl[[j]]$X <- by.dum*sm$X  ## multiply model matrix by dummy for level
  #           objectl[[j]]$by.level <- lev[j] ## store level
  #           objectl[[j]]$label <- paste(object$label,":",object$by,lev[j],sep="")
  #           if (!is.null(offs)) {
  #             attr(objectl[[j]]$X,"offset") <- offs*by.dum
  #           }
  #         }
  #       } else { ## not a factor by variable
  #         objectl <- list(object)
  #         if ((is.null(object$IND)&&length(by)!=nrow(object$X)))
  #           stop("`by' variable must be same dimension as smooth arguments")
  #         if (object$by == "NA") objectl[[1]]$label <- object$label else
  #           objectl[[1]]$label <- paste(object$label,":",object$by,sep="")
  #       }
  # } else { ## no by variables
    class(object)<-c(class(object),"ggam.smooth")
    objectl <- list(object)
  # }
  objectl
}
