\name{eg2pop_dat}
\alias{eg2pop_dat}
\docType{data}
%- Also NEED an `\alias' for EACH other topic documented here.
\title{Rectangular Domain Simulation Data}
\description{This dataset includes the simulated data from rectangular domain with  20×20 lattice grid points and error term from stationary gaussian process with mean zero.}
\usage{
data("eg2pop_dat")}
%- maybe also `usage' for other objects documented here.

\format{ \code{eg2pop_dat} is a list with three data frames: \code{pop}, \code{V20} and \code{T20}. \code{V20} is the verticies dataset and \code{T20} is the corresponding indices of vertices in \code{V}s.
\describe{
\code{pop}:
\item{Y}{Response variable.}
\item{m}{The true values for nonlinear part.}
\item{x1}{The location parameter (longitude).}
\item{x2}{The location parameter (latitude).}
\item{z1, z2, z3, z4,z5, z6, z7, z8}{The parametric covariates.}
}}
\details{
Here is the data generating process:
\eqn{Y _ { i } = \mathbf { Z } _ { i } ^ { \top } \boldsymbol { \beta } + \varepsilon _ { i } , i = 1 , \ldots , n}{Y _ { i } = \mathbf { Z } _ { i } ^ { \top } \boldsymbol { \beta } + \varepsilon _ { i } , i = 1 , \ldots , n}, where \eqn{\boldsymbol { \beta } = ( 1 , - 1,0,0,0,0,0,0 ) ^ { \top }}{\boldsymbol { \beta } = ( 1 , - 1,0,0,0,0,0,0 ) ^ { \top }}. All covariates are generated from \eqn{N(0,1)}{N(0,1)}.
}

\references{Wang et al. (2018)}



\examples{
## see ?plbpsm help files
}
