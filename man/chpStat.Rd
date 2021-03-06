% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RcppExports.R
\name{chpStat}
\alias{chpStat}
\title{CHP Test Statistic}
\usage{
chpStat(mdl, rho_b, ltmt, var_switch)
}
\arguments{
\item{mdl}{List containing model information}

\item{rho_b}{bound for rho (nuisance param space)}

\item{ltmt}{List containing  relevant first and second derivatives of log likelihood function.}

\item{var_switch}{variance switch indicator}
}
\value{
Test Statistic
}
\description{
Calculate supTS and expTS test-statistics from CHP (2014).
}
\references{
Carrasco, Marine, Liang Hu, and Werner Ploberger. 2014. 
“Optimal test for Markov switch- ing parameters.” \emph{Econometrica} \bold{82 (2)}: 765–784.
}
