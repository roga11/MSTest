% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/models.R
\name{VARXmdl}
\alias{VARXmdl}
\title{Vector X autoregressive model}
\usage{
VARXmdl(Y, p, Z, control = list())
}
\arguments{
\item{Y}{a \code{(T x q)} matrix of observations.}

\item{p}{integer determining the number of autoregressive lags.}

\item{Z}{a \code{(T x qz)} matrix of exogenous regressors.}

\item{control}{List with model options including:
\itemize{
 \item const: Boolean determining whether to estimate model with constant if \code{TRUE} or not if \code{FALSE}. Default is \code{TRUE}.
 \item getSE: Boolean determining whether to compute standard errors of parameters if \code{TRUE} or not if \code{FALSE}. Default is \code{TRUE}.
}}
}
\value{
List of class \code{VARmdl} (\code{S3} object) with model attributes including:
\itemize{
  \item y: a \code{(T-p x q)} matrix of observations.
  \item X: a \code{(T-p x p*q + const)} matrix of lagged observations with a leading column of \code{1}s if \code{const=TRUE} or not if \code{const=FALSE}.
  \item x: a \code{(T-p x p*q)} matrix of lagged observations.
  \item fitted: a \code{(T-p x q)} matrix of fitted values.
  \item resid: a \code{(T-p x q)} matrix of residuals.
  \item mu: a \code{(1 x q)} vector of estimated means of each process.
  \item beta: a \code{((1 + p + qz) x q)} matrix of estimated coefficients. 
  \item betaZ: a \code{(qz x q)} matrix of estimated exogenous regressor coefficients.
  \item intercept: estimate of intercepts.
  \item phi: a \code{(q x p*q)} matrix of estimated autoregressive coefficients.
  \item Fmat: Companion matrix containing autoregressive coefficients.
  \item stdev: a \code{(q x 1)} vector of estimated standard deviation of each process.
  \item sigma: a \code{(q x q)} estimated covariance matrix.
  \item theta: vector containing: \code{mu}, \code{vech(sigma)}, and \code{vec(t(phi))}.
  \item theta_mu_ind: vector indicating location of mean with \code{1} and \code{0} otherwise.
  \item theta_sig_ind: vector indicating location of variance and covariances with \code{1} and \code{0} otherwise.
  \item theta_var_ind: vector indicating location of variances with \code{1} and \code{0} otherwise.
  \item theta_phi_ind: vector indicating location of autoregressive coefficients with \code{1} and \code{0} otherwise.
  \item stationary: Boolean indicating if process is stationary if \code{TRUE} or non-stationary if \code{FALSE}.
  \item n: number of observations after lag transformation (i.e., \code{n = T-p}).
  \item p: number of autoregressive lags.
  \item q: number of series.
  \item k: number of regimes. This is always \code{1} in \code{VARmdl}.
  \item Fmat: matrix from companion form. Used to determine is process is stationary.
  \item control: List with model options used.
  \item logLike: log-likelihood.
  \item AIC: Akaike information criterion.
  \item BIC: Bayesian (Schwarz) information criterion.
  \item Hess: Hessian matrix. Approximated using \code{\link[numDeriv]{hessian}} and only returned if \code{getSE=TRUE}.
  \item info_mat: Information matrix. Computed as the inverse of \code{-Hess}. If matrix is not PD then nearest PD matrix is obtained using \code{\link[pracma]{nearest_spd}}. Only returned if \code{getSE=TRUE}.
  \item nearPD_used: Boolean determining whether \code{nearPD} function was used on \code{info_mat} if \code{TRUE} or not if \code{FALSE}. Only returned if \code{getSE=TRUE}.
  \item theta_se: standard errors of parameters in \code{theta}.  Only returned if \code{getSE=TRUE}.
}
}
\description{
This function estimates a vector autoregresive model with \code{p} lags. This can be used for the null hypothesis of a linear model against an alternative hypothesis of a Markov switching vector autoregressive model with \code{k} regimes.
}
\examples{
# ----- Bivariate VAR(1) process ----- #
set.seed(1234)
# Define DGP of VAR process
mdl_var <- list(n     = 1000, 
                p     = 1,
                q     = 2,
                mu    = c(5,-2),
                sigma = rbind(c(5.0, 1.5),
                              c(1.5, 1.0)),
                phi   = rbind(c(0.50, 0.30),
                              c(0.20, 0.70)))

# Simulate process using simuVAR() function
y_simu <- simuVAR(mdl_var)

# Set options for model estimation
control <- list(const  = TRUE, 
                getSE  = TRUE)

# Estimate model
y_var_mdl <- VARmdl(y_simu$y, p = 2, control = control)
summary(y_var_mdl)
}
\seealso{
\code{\link{MSVARmdl}}
}
