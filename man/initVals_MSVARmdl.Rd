% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RcppExports.R
\name{initVals_MSVARmdl}
\alias{initVals_MSVARmdl}
\title{Initial values for Markov-switching vector autoregressive model}
\usage{
initVals_MSVARmdl(mdl, k)
}
\arguments{
\item{\code{mdl}}{List with parameter values of simple (one-regime) vector autoregressive model. This includes:
  \itemize{
   \item{phi: }{Matrix autoregressive coefficients.}
   \item{mu: }{Vector of means.}
   \item{sigma: }{Covariance matrix.}
   \item{msmu: }{Boolean indicator. If \code{TRUE}, mean is function of markov process. If \code{FALSE}, mean is constant across regimes.}
   \item{msvar: }{Boolean indicator. If \code{TRUE}, standard deviation is function of markov process. If \code{FALSE}, standard deviation is constant across regimes.}
}}

\item{\code{k}}{Number of regimes.}
}
\value{
Vector of initial parameter values.
}
\description{
This function generates a random parameter vector to be used as initial values for a Markov-switching vector autoregressive model.
}
\keyword{internal}
