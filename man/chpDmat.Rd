% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/CHP.R
\name{chpDmat}
\alias{chpDmat}
\title{Derivative matrix}
\usage{
chpDmat(Mdl, var_switch = 1)
}
\arguments{
\item{Mdl}{is a list containing AR model components
Specifically, it containg y, x, X (x plus constant), residuals, coefficients, stdev,
logLike}

\item{var_switch}{is an indicator = 1 if there is a switch in both Mean and Variance 
and = 0 if there is only a switch in the Mean. Less Second-Order dervatives are 
calculated if only the Mean is subject to regime switch.}
}
\value{
List containing relevant first and second derivatves of log likelihood function.
}
\description{
This function organizes the first and second derivatives of the log Likelihoood.
}
