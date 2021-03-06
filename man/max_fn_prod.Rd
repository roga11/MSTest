% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/moment_test.R
\name{max_fn_prod}
\alias{max_fn_prod}
\title{Product of p-value maximizing function}
\usage{
max_fn_prod(v)
}
\arguments{
\item{v}{nuissance parameter values}
}
\value{
p-value for given parameter value.
}
\description{
maximization function for nuissance parameters. This version uses min method 
of combining p-values as in  Fisher (1932) and Pearson (1933).
}
\references{
Dufour, J. M., & Luger, R. (2017). Identification-robust moment-based 
tests for Markov switching in autoregressive models. Econometric Reviews, 36(6-9), 713-727.

Pearson, K. (1933). On a method of determining whether a sample of size n
 supposed to have been drawn from a parent population having a known probability integral has probably been drawn at random. Biometrika 25:379–410.

Fisher, R. (1932). Statistical Methods for Research Workers. Edinburgh: 
Oliver and Boyd.
}
