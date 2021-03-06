% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Methods.R
\name{calc_mmcpval}
\alias{calc_mmcpval}
\title{Calculate MMC p-value}
\usage{
calc_mmcpval(search_type = "randSearch_paramCI", optimOptions = NULL)
}
\arguments{
\item{search_type}{Type of optimization algorithm when searching nuissance 
parameter space. Avaiable options are: GenSA, GA, PSO, randSearch_paramCI and 
gridSearch_paramCI. Default is set to randSearch_paramCI to match results in paper.}

\item{optimOptions}{List containing upper "upp" and lower "low" bounds of 
nuissance parameter space.}
}
\value{
MMC p-value of test
}
\description{
This function calculates a Maximized Monte-Carlo p-value
}
\references{
Dufour, J. M. (2006). Monte Carlo tests with nuisance parameters: 
A general approach to finite-sample inference and nonstandard asymptotics. 
Journal of Econometrics, 133(2), 443-477.

Dufour, J. M., & Luger, R. (2017). Identification-robust moment-based 
tests for Markov switching in autoregressive models. Econometric Reviews, 36(6-9), 713-727.
}
