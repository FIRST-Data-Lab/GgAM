#' Basis Construction Function
#'
#' This is an internal function of package \code{ggam}.
#' @importFrom BPST basis qrH
#' @importFrom Matrix Matrix
#' @importFrom mgcv get.var
#' @importFrom geometry delaunayn
#' @importFrom Triangulation TriMesh
#' @param object is a bivariate smooth object.
# \item{term}{The names of the covariates for this smooth, in an array.}
# \item{fixed}{\code{TRUE} if the term is to be unpenalized, otherwise \code{FALSE}.}
# \item{dim}{the number covariates of which this smooth is a function.}
# \item{by}{the name of any \code{by} variable to multiply this term as supplied as an argument to \code{s}.
#   \code{"NA"} if there is no such term.}
# \item{label}{A suitable label for use with this term.}
# \item{id}{Any identity associated with this term --- used for linking bases
#   and smoothing parameters. \code{NULL} by default, indicating no linkage.}
# \item{r}{Smoothing parameter for the bivariate term.}
# \item{d}{degree of polynomial parameter for the term.}
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
#'BI <- BasisCon(b(x1,x2,d=2,r=1,V=eg1_V2,Tr=eg1_T2),eg1pop_rho03[sample(1:dim(eg1pop_rho03)[1],70),])
#'xx=seq(-0.8,3.8,length=100)
#'yy=seq(-0.8,0.8,length=100)
#'sam=data.frame(x1=rep(xx,length(xx)),x2=rep(yy,each=length(yy)))
#' BI <- BasisCon(b(x1,x2,d=2,r=1,V=eg1_V2,Tr=eg1_T2),sam)

#' @export
######################################################
BasisCon.bivariate.smooth <- function(object,data)
  ## wrapper function which calls basis constructing method
  ## Handles `by' variables, and summation convention.
  ## Note that `data' must be a data.frame or model.frame, unless n is provided explicitly,
  ## in which case a list will do.
  ## If present dataX specifies the data to be used to set up the model matrix, given the
  ## basis set up using data (but n same for both).
{
  ###
  loc=data[, object$term]
  if (is.null(object$V) & is.null(object$Tr)){
    object$Tr=TriMesh(object$b,n=object$nt,object$Holes)$Tr
    object$V=TriMesh(object$b,n=object$nt,object$Holes)$VT
  }
  if (is.null(object$B) || is.null(object$Q2) || is.null(object$K) || is.null(object$ind)){
    if ((!is.null(object$d)) & (!is.null(object$r))){
      B0=basis(as.matrix(object$V),as.matrix(object$Tr),object$d,object$r,as.matrix(loc))
    } else {
      B0=basis(as.matrix(object$V),as.matrix(object$Tr),Z=as.matrix(loc))
    }
    object$d=B0$d
    object$r=B0$r
    object$B=B0$B
    object$ind=B0$Ind.inside
    object$K=B0$K
    object$H=B0$H
    object$Q2=B0$Q2
  }
  object$n.para=dim(object$Q2)[1]
  object$loc=loc[object$ind,]
  null_sp_dim <- getFromNamespace("null.space.dimension", "mgcv")

  object$null.space.dimension=null_sp_dim(object$dim,0)
  ### should be the X outside since BasisCon does not deal with para part?

  #object$datanew=data[object$ind,]

  ## add plotting indicator if not present.
  ## plot.me tells `plot.plbpsm' whether or not to plot the term
  plot.me <- FALSE
  object$n.para=ncol(object$B)

  # may change in univariate case
  object$IND<-NULL
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
