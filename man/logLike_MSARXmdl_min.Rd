% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RcppExports.R
\name{logLike_MSARXmdl_min}
\alias{logLike_MSARXmdl_min}
\title{Markov-switching ARX log-likelihood objective function (minimization version)}
\usage{
logLike_MSARXmdl_min(theta, mdl, k)
}
\arguments{
\item{theta}{Vector of model parameters.}

\item{mdl}{List with model attributes.}

\item{k}{integer determining the number of regimes.}
}
\value{
Negative log-likelihood value.
}
\description{
This function computes the (negative) log-likelihood for a markov-switching ARX model.
}
\keyword{internal}
