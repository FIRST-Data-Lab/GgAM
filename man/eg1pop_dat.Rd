\name{eg1pop_dat_gauss}
\alias{eg1pop_dat_gauss}
\docType{data}
%- Also NEED an `\alias' for EACH other topic documented here.
\title{Horseshore Domain Simulation Data}
\description{This dataset includes the simulated data from horseshore domain with three different triangulations and three different correlation structure of the error. Useful for simulation in geospacial models. }
\usage{
data("eg1pop_dat")}
%- maybe also `usage' for other objects documented here.

\format{ \code{eg1pop_dat} is a list with 9 data frames: \code{rho03}, \code{rho05} and \code{rho07} are the population datasets; \code{V1}, \code{V2}, \code{V3} are the verticies datasets corresponding to different triangulations; \code{T}s are corresponding indices of vertices in \code{V}s.
\describe{
\code{rho03}, \code{rho05} and \code{rho07}:
\item{Y}{Response variable.}
\item{m}{The true values for nonlinear part.}
\item{x1}{The location parameter (longitude).}
\item{x2}{The location parameter (latitude).}
\item{z1, z2, z3, z4,z5, z6, z7, z8}{The parametric covariates.}
}}


\details{The finite area test function based on one is proposed by Tim Ramsay (2002).
To study the performance of variable selection at different correlation
levels, we generate the covariates as follows:
\eqn{Z _ { i 1 } = - \frac { 2 } { 3 } \left\{ \arctan \pi \left( \rho \frac { X _ { i 1 } } { X _ { i 2 } } + ( 1 - \rho ) U _ { i } \right) \right\}}{Z _ { i 1 } = - \frac { 2 } { 3 } \left\{ \arctan \pi \left( \rho \frac { X _ { i 1 } } { X _ { i 2 } } + ( 1 - \rho ) U _ { i } \right) \right\}}, \eqn{Z _ { i 3 } = \cos \pi \left( \rho \frac { X _ { i 1 } } { X _ { i 2 } } + ( 1 - \rho ) U _ { i } \right)}{Z _ { i 3 } = \cos \pi \left( \rho \frac { X _ { i 1 } } { X _ { i 2 } } + ( 1 - \rho ) U _ { i } \right)},\eqn{Z _ { i j } \sim Uniform  ( - 1,1 )}{Z _ { i j } \sim  Uniform ( - 1,1 )}, j = 2, 4,...,8, \eqn{U _ { i } \sim  { Uniform } ( - 1,1 )}{U _ { i } \sim { Uniform } ( - 1,1 )}.
In particular, we consider the following three cases: (i) low correlation (\eqn{\rho = 0.3}{\rho = 0.3}); (ii) medium correlation (\eqn{\rho = 0.5}{\rho = 0.5}); and (iii) high correlation (\eqn{\rho = 0.7}{\rho = 0.7}).
}

\references{Wang et al. (2018)}



\examples{
## see ?plbpsm help files
}
