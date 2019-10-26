#' Print a Bivariate Penalized Spline based on Triangulation object.
#'
#' The default print method for a \code{plbpsm} object.
#' @param x,... fitted model objects of class plbpsm as produced by \code{plbpsm()}.
#' @details
#' Prints out the family, model formula, and etc. (need to be decided) See \code{\link{plbpsmObject}} (or \code{names(x)}) for a listing
#'  of what the object contains. \code{\link{summary.plbpsm}} provides more detail.
#'  Note that the optimized smoothing penalty parameter selection criterion reported is one of GCV, CV.
#'@export
print.plbpsm<- function (x,...) {
  print(x$family)
  cat("Formula:\n")
  if (is.list(x$formula))
    for (i in 1:length(x$formula)) print(x$formula[[i]]) else print(x$formula)
  n.smooth <- length(x$basis_info)
  # if (n.smooth == 0)
  #   cat("Total model degrees of freedom", sum(x$edf), "\n") else {
  #   edf <- 0
  # cat("\nEstimated degrees of freedom:\n")
  ## edf from different parts of smooth
  # for (i in 1:n.smooth) edf[i] <- x$edf[x$basis_info[[i]]$first.para]
  # edf.str <- format(round(edf, digits = 4), digits = 3,
  #                   scientific = FALSE)
  # for (i in 1:n.smooth) {
  #   cat(edf.str[i], " ", sep = "")
  #   if (i%%7 == 0)
  #     cat("\n")
  # }
  cat(" total edf =", round(sum(x$edf), digits = 2), "\n")
  #}
  if (!is.null(x$criterion) && (x$criterion %in% c("GCV"
  )))
    cat("\n", x$criterion, " score: ", x$gcv_opt, "     ",
        sep = "")

  if (!is.null(x$criterion) && (x$criterion %in% c("CV"
  )))
    cat("\n", x$criterion, " score: ", x$cv_opt, "     ",
        sep = "")
  ## rank count non-para?
  if (!is.null(x$rank) && x$rank < length(x$coefficients_linear)+
      length(x$coefficients_nonlinear))
    cat("rank: ", x$rank, "/", length(x$coefficients_linear)+length(x$coefficients_nonlinear), sep = "")
  cat("\n")
  invisible(x)
}
