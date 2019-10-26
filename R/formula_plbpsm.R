#' PLBPSM formula
#'
#' Description of \link{plbpsm} formula (see Details), and how to extract it from a fitted \code{plbpsm} object.
#' @param x fitted model objects of class gam (see \code{\link{plbpsmObject}}) as produced by \code{plbpsm()}.
#' @param ... un-used in this case
#' @method formula plbpsm
#' @details
#' The formulae supplied to \code{\link{plbpsm}} are exactly like those supplied to \code{\link{glm}}
#' except that univataible and bivariate smooth terms,\code{\link{u}} and \code{\link{b}} can be added to the right hand side (and . is not supported
#'  in \code{plbpsm} formulae).
#' Smooth terms are specified by expressions of the form:
#' \code{f(x1,x2,...,r=1,d=2,fx=FALSE)}.
#' If \code{d} is not specified then basis specific defaults are
#' used. Note that these defaults are essentially arbitrary, and it is important to check that they are not
#' so big that they cause oversmoothing.
#' fx is used to indicate whether or not this term should be unpenalized.
#' Formulae can involve nested or “overlapping” terms such as
#' \code{y~u(x)+u(z)+f(x,z)}.
#' @return
#' Returns the model formula, \code{x$formula}. Provided so that anova methods print an appropriate description of the model.
#' @seealso \code{\link{plbpsm}}
#' @export

formula.plbpsm <- function(x, ...)
  # formula.lm and formula.glm reconstruct the formula from x$terms, this is
  # problematic because of the way ggam handles f() terms.
{ x$formula
}
