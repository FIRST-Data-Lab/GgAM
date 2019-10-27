\name{plbpsmObject}
\alias{plbpsmObject}
%- Also NEED an `\alias' for EACH other topic documented here.
\title{Fitted plbpsm object}
\description{A fitted PLBPSM object returned by function \code{plbpsm} and of class
\code{"plbpsm"} inheriting from classes \code{"gam"}, \code{"glm"} and \code{"lm"}. Method
functions \code{anova}, \code{logLik}, \code{influence}, \code{plot},
\code{predict}, \code{print}, \code{residuals} and \code{summary} exist for
this class.

}
%All compulsory elements of \code{"glm"} and \code{"lm"} objects are present,
%but the fitting method for a PLBPSM is different to a linear model or GLM, so
%that the elements relating to the QR decomposition of the model matrix are
%absent.
\value{ A \code{plbpsm} object has the following elements:
\item{backfitting}{Whether Spline-backfitted local polynomial (SBL) estimate is applied.}
\item{boundary}{did parameters end up at boundary of parameter space?}
\item{call}{the matched call (allows \code{update} to be used with \code{plbpsm} objects, for example). }
\item{coefficients_linear}{the linear (parametric) coefficients of the fitted model.}
\item{coefficients_nonlinear}{the nonlinear (bivariate spline) coefficients of the fitted model.}
\item{converged}{indicates whether or not the iterative fitting method converged.}
\item{data}{the original supplied data argument (for class \code{"gam"} compatibility).}
\item{deviance}{model deviance (not penalized deviance).}
\item{dev_sbl}{model deviance with SBL estimates.}
\item{df.null}{null degrees of freedom.}
\item{df.residual}{effective residual degrees of freedom of the model.}
\item{edf}{estimated degrees of freedom.}
\item{family}{family object specifying distribution and link used.}
\item{fitted.values}{fitted model predictions of expected value for each
datum.}
\item{fitted.values.sbl}{fitted model predictionsof expected value by using SBL estimation for each
datum.}
\item{formula}{the model formula.}
\item{gcv_opt}{The minimized smoothing penalty parameter selection score: GCV.}
\item{cv_opt}{The minimized smoothing penalty parameter selection score: CV.}
\item{linear.predictors}{fitted model prediction of link function of
expected value for  each datum.}
\item{criterion}{One of \code{"GCV"} or \code{CV}, depending on the fitting
criterion used.}
\item{iter}{number of iterations taken to get convergence.}
\item{model}{model frame containing all variables needed in original model fit.}
\item{na.action}{The \code{na.action} used in fitting.}
\item{null.deviance}{deviance for single parameter model.}
\item{method}{One of \code{ALASSO}, \code{SCAD}, depending on the variable selection method.}
\item{pred.formula}{one sided formula containing variables needed for prediction, used by \code{predict.plbpsm}}
\item{pterms}{\code{terms} object for strictly parametric part of model.}
\item{residuals}{the working residuals for the fitted model.}
\item{scale}{defualt=1.}
\item{se_beta}{estimated standard error for parametric covariates.}
\item{sigma_2}{estimated scale parameter.}
\item{basis_info}{list of smooth objects, containing the basis information for each term in the
model formula in the order in which they appear. These smooth objects are what gets returned by
the \code{\link{BasisCon}} objects.}
\item{terms}{\code{terms} object of \code{model} model frame.}
\item{var.summary}{A named list of summary information on the predictor variables}
\item{Ve}{estimated covariance matrix for the parameter
estimators. Particularly useful for testing whether terms are zero. Not so
useful for CI's as smooths are usually biased.}
\item{prior.weights}{prior weights of observations.}
\item{X2}{The design matrix with linear covariates and univariate functions.}
\item{mhat}{estimated component functions in the evaluated in the sample points.}
\item{y}{response data.}
\item{est_theta}{estimated theta in the negative binomial family.}
\item{ind.l}{The index for linear covariates in model identification.}
\item{ind.nl}{The index for nonlinear covariates in model identification.}
\item{ind_c}{The index for (selected) linear covariates.}
\item{VS}{Varaible selection conducted or not.}
\item{h_opt_all}{optimal bandwidth for each component functions.}
\item{backfitting}{Whether backfitting agorithm applied and SBL estimator is and used.}
\item{weights}{final weights used in the iteration.}
\item{Xp}{The design matrix for the parametric part of the model.}
\item{nsdf}{number of parametric, non-smooth, model terms including the intercept.}
}
