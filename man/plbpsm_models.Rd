\name{plbpsm.models}
\alias{plbpsm.models}
%- Also NEED an `\alias' for EACH other topic documented here.
\title{Specifying generalized geospatial additive models}
\description{ This page is intended to provide some more information on
 how to specify GGAMs. A GGAM is a GLM in which the linear predictor depends,
in part, on a sum of smooth functions of predictors and (possibly) linear
functionals of smooth functions of (possibly dummy) predictors.
Specifically let \eqn{y_i}{y_i} denote an independent random variable
with mean \eqn{\mu_i}{mu_i} and an exponential family distribution, or failing
that a known mean variance relationship suitable for use of quasi-likelihood methods.
Then the the linear predictor of a GGAM has a structure something like

\deqn{g(\mu_i) = {\bf X}_i{\beta} + u_1(x_{1i}) + u_2(x_{3i}) + L_i f_3(x_4) +
\ldots}{g(mu_i)=X_i b + f_1(x_1i,x_2i) + f_2(x_3i) + L_i f_3(x_4) + ...}
where \eqn{g}{g} is a known smooth monotonic `link' function, \eqn{{\bf X}_i\beta}{X_i b}
is the parametric part of the linear predictor, the \eqn{x_j}{x_j} are predictor variables,
the \eqn{u_j}{u_j} and \eqn{f_j}{f_j} are smooth functions. There may of course be multiple linear functional terms, or none.

The key idea here is that the
dependence of the response on the predictors can be represented as a
parametric sub-model plus the sum of some (functionals of) smooth functions of one or
more of the predictor variables. Thus the model is quite flexible
relative to strictly parametric linear or generalized linear models,
but still has much more structure than the completely general model
that says that the response is just some smooth function of all the
covariates.

The following sections discuss specifying model structures for \code{gam}.
Specification of the distribution and link function is done using the \code{\link{family}}
 argument to \code{\link{plbpsm}} and works in the same way as for \code{\link{glm}}.
This page therefore concentrates on the model formula for \code{plbpsm}.

}
\section{Models with simple smooth terms}{
Consider the example model.
\deqn{g(\mu_i) = \beta_0 + \beta_1 x_{1i} + \beta_2 x_{2i} + u(x_{3i}) +
f(x_{4i},x_{5i})}{ g(mu_i) = b_0 + b_1 x_1i + b_2 x_2i + u(x_3i) + f(x_4i,x_5i)}
where the response variables \eqn{y_i}{y_i} has expectation \eqn{\mu_i}{mu_i}
and \eqn{g}{g} is a link function.

The \code{plbpsm} formula for this would be \cr
  \code{y ~ x1 + x2 + u(x3) + f(x4,x5)}. \cr
  This would use the default basis for the smooths (a thin plate
  regression spline basis for each), with default choices of degree of polynomials for both smooths. Full details of how to control smooths are given in
  \code{\link{s}} and \code{\link{u}}.

The above assumes that \eqn{x_{4}}{x_4} and \eqn{x_5}{x_5} are naturally on
similar scales (e.g. they might be co-ordinates), so that isotropic smoothing
is appropriate. \cr
}
%All compulsory elements of \code{"glm"} and \code{"lm"} objects are present,
%but the fitting method for a PLBPSM is different to a linear model or GLM, so
%that the elements relating to the QR decomposition of the model matrix are
%absent.
\value{ A \code{plbpsm} object has the following elements:
\item{aic}{AIC of the fitted model: bear in mind that the degrees of freedom
used to calculate this are the effective degrees of freedom of the model, and
the likelihood is evaluated at the maximum of the penalized likelihood in most
cases, not at the MLE.}
\item{assign}{Array whose elements indicate which model term (listed in
\code{pterms}) each parameter relates to: applies only to non-smooth terms.}
\item{call}{the matched call (allows \code{update} to be used with \code{plbpsm} objects, for example). }
\item{cmX}{column means of the model matrix (with elements corresponding to smooths set to zero )
--- useful for componentwise CI calculation.}
\item{coefficients_linear}{the linear (parametric) coefficients of the fitted model.}
\item{coefficients_nonlinear}{the nonlinear (bivariate spline) coefficients of the fitted model.}
\item{data}{the original supplied data argument (for class \code{"gam"} compatibility).}
\item{deviance}{model deviance (not penalized deviance).}
\item{df.null}{null degrees of freedom.}
\item{df.residual}{effective residual degrees of freedom of the model.}
\item{edf}{estimated degrees of freedom.}
\item{family}{family object specifying distribution and link used.}
\item{fitted.values}{fitted model predictions of expected value for each
datum.}
\item{formula}{the model formula.}
\item{gcv_opt}{The minimized smoothing penalty parameter selection score: GCV.}
\item{cv_opt}{The minimized smoothing penalty parameter selection score: CV.}
\item{linear.predictors}{fitted model prediction of link function of
expected value for  each datum.}
\item{criterion}{One of \code{"GCV"} or \code{CV}, depending on the fitting
criterion used.}
\item{min.edf}{Minimum possible degrees of freedom for whole model.}
\item{model}{model frame containing all variables needed in original model fit.}
\item{na.action}{The \code{na.action} used in fitting.}
\item{nsdf}{number of parametric, non-smooth, model terms including the intercept.}
\item{null.deviance}{deviance for single parameter model.}
\item{method}{One of \code{ALASSO}, \code{SCAD}, depending on the variable selection method.}
\item{pred.formula}{one sided formula containing variables needed for prediction, used by \code{predict.plbpsm}}
\item{pterms}{\code{terms} object for strictly parametric part of model.}
\item{rank}{apparent rank of fitted model.}
\item{residuals}{the working residuals for the fitted model.}
\item{scale}{defualt=1.}
\item{se_beta}{estimated standard error for parametric covariates.}
\item{sse}{Sum of squared error of the estimation.}
\item{basis_info}{list of smooth objects, containing the basis information for each term in the
model formula in the order in which they appear. These smooth objects are what gets returned by
the \code{\link{BasisCon}} objects.}
\item{terms}{\code{terms} object of \code{model} model frame.}
\item{var.summary}{A named list of summary information on the predictor variables}
\item{Ve}{estimated covariance matrix for the parameter
estimators. Particularly useful for testing whether terms are zero. Not so
useful for CI's as smooths are usually biased.}
\item{weights}{final weights used in grplsfit.}
\item{y}{response data.}

}
