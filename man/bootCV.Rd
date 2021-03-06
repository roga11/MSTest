% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RcppExports.R
\name{bootCV}
\alias{bootCV}
\title{Bootstrap Critival Values CHP Test}
\usage{
bootCV(mdl, rho_b, N, var_switch)
}
\arguments{
\item{mdl}{List containing model information}

\item{rho_b}{bound for rho (nuisance param space)}

\item{N}{number of simulations}

\item{var_switch}{variance switch indicator}
}
\value{
Bootstrap critical values
}
\description{
This bootstrap procedure is described on page 771 of CHP (2014)
}
\references{
Carrasco, Marine, Liang Hu, and Werner Ploberger. 2014. 
“Optimal test for Markov switch- ing parameters.” \emph{Econometrica} \bold{82 (2)}: 765–784.
}
