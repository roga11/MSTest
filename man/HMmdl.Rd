% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/models.R
\name{HMmdl}
\alias{HMmdl}
\title{Hidden Markov model}
\usage{
HMmdl(Y, k, Z = NULL, control = list())
}
\arguments{
\item{Y}{a \code{(T x q)} matrix of observations.}

\item{k}{integer determining the number of regimes to use in estimation. Must be greater than or equal to \code{2}.}

\item{Z}{an otpional  \code{(T x qz)} matrix of exogenous regressors. Default is NULL.}

\item{control}{List with model options including:
\itemize{
 \item getSE: Boolean. If \code{TRUE} standard errors are computed and returned. If \code{FALSE} standard errors are not computed. Default is \code{TRUE}.
 \item msmu: Boolean. If \code{TRUE} model is estimated with switch in mean. If \code{FALSE} model is estimated with constant mean. Default is \code{TRUE}.
 \item msvar: Boolean. If \code{TRUE} model is estimated with switch in variance. If \code{FALSE} model is estimated with constant variance. Default is \code{TRUE}.
 \item init_theta: vector of initial values. vector must contain \code{(1 x q)} vector \code{mu}, \code{vech(sigma)}, and \code{vec(P)} where sigma is a \code{(q x q)} covariance matrix.This is optional. Default is \code{NULL}, in which case \code{\link{initVals_MSARmdl}} is used to generate initial values.
 \item method: string determining which method to use. Options are \code{'EM'} for EM algorithm or \code{'MLE'} for Maximum Likelihood Estimation. Default is \code{'EM'}.
 \item maxit: integer determining the maximum number of EM iterations.
 \item thtol: double determining the convergence criterion for the absolute difference in parameter estimates \code{theta} between iterations. Default is \code{1e-6}.
 \item maxit_converge: integer determining the maximum number of initial values attempted until solution is finite. For example, if parameters in \code{theta} or \code{logLike} are \code{NaN} another set of initial values (up to \code{maxit_converge}) is attempted until finite values are returned. This does not occur frequently for most types of data but may be useful in some cases. Once finite values are obtained, this counts as one iteration towards \code{use_diff_init}. Default is \code{500}.
 \item use_diff_init: integer determining how many different initial values to try (that do not return \code{NaN}; see \code{maxit_converge}). Default is \code{1}.
 \item mle_variance_constraint: double used to determine the lower bound on the smallest eigenvalue for the covariance matrix of each regime. Default is \code{1e-3}.
 \item mle_theta_low: Vector with lower bounds on parameters (Used only if method = "MLE"). Default is \code{NULL}.
 \item mle_theta_upp: Vector with upper bounds on parameters (Used only if method = "MLE"). Default is \code{NULL}.
}}
}
\value{
List of class \code{HMmdl} (\code{S3} object) with model attributes including:
\itemize{
  \item y: a \code{(T x q)} matrix of observations.
  \item fitted: a \code{(T x q)} matrix of fitted values.
  \item resid: a \code{(T x q)} matrix of residuals.
  \item mu: a \code{(k x q)} matrix of estimated means of each process.
  \item beta: if \code{q=1}, this is a \code{((1 + qz) x k)} matrix of estimated coefficients. If \code{q>1}, this is a list containing \code{k} separate \code{((1 + qz) x q)} matrix of estimated coefficients for each regime.  
  \item betaZ: a \code{(qz x q)} matrix of estimated exogenous regressor coefficients.
  \item intercept: a \code{(k x q)} matrix of estimated intercept of each process. If Z is Null, this is the same as mu.
  \item stdev: If \code{q=1}, this is a \code{(k x 1)} matrix with estimated standard. If \code{q>1}, this is a List with \code{k} \code{(q x q)} matrices with estimated standard deviation on the diagonal.
  \item sigma: If \code{q=1}, this is a \code{(k x 1)} matrix with variances. If \code{q>1}, this is a List with \code{k} \code{(q x q)} estimated covariance matrix.
  \item theta: vector containing: \code{mu} and \code{vech(sigma)}.
  \item theta_mu_ind: vector indicating location of mean with \code{1} and \code{0} otherwise.
  \item theta_sig_ind: vector indicating location of variance and covariances with \code{1} and \code{0} otherwise.
  \item theta_var_ind: vector indicating location of variances with \code{1} and \code{0} otherwise.
  \item theta_P_ind: vector indicating location of transition matrix elements with \code{1} and \code{0} otherwise.
  \item n: number of observations (same as \code{T}).
  \item q: number of series.
  \item k: number of regimes in estimated model.
  \item P: a \code{(k x k)} transition matrix.
  \item pinf: a \code{(k x 1)} vector with limiting probabilities of each regime.
  \item St: a \code{(T x k)} vector with smoothed probabilities of each regime at each time \code{t}.
  \item deltath: double with maximum absolute difference in vector \code{theta} between last iteration.
  \item iterations: number of EM iterations performed to achieve convergence (if less than \code{maxit}).
  \item theta_0: vector of initial values used.
  \item init_used: number of different initial values used to get a finite solution. See description of input \code{maxit_converge}.
  \item msmu: Boolean. If \code{TRUE} model was estimated with switch in mean. If \code{FALSE} model was estimated with constant mean.
  \item msvar: Boolean. If \code{TRUE} model was estimated with switch in variance. If \code{FALSE} model was estimated with constant variance.
  \item control: List with model options used.
  \item logLike: log-likelihood.
  \item AIC: Akaike information criterion.
  \item BIC: Bayesian (Schwarz) information criterion.
  \item Hess: Hessian matrix. Approximated using \code{\link[numDeriv]{hessian}} and only returned if \code{getSE=TRUE}.
  \item info_mat: Information matrix. Computed as the inverse of \code{-Hess}. If matrix is not PD then nearest PD matrix is obtained using \code{\link[pracma]{nearest_spd}}. Only returned if \code{getSE=TRUE}.
  \item nearPD_used: Boolean determining whether \code{nearPD} function was used on \code{info_mat} if \code{TRUE} or not if \code{FALSE}. Only returned if \code{getSE=TRUE}.
  \item theta_se: standard errors of parameters in \code{theta}.  Only returned if \code{getSE=TRUE}.
  \item trace: List with Lists of estimation output for each initial value used due to \code{use_diff_init > 1}.
}
}
\description{
This function estimates a Hidden Markov model with \code{k} regimes.
}
\references{
Dempster, A. P., N. M. Laird, and D. B. Rubin. 1977. “Maximum Likelihood from Incomplete Data via the EM Algorithm.” \emph{Journal of the Royal Statistical Society}. Series B 39 (1): 1–38..

Hamilton, James D. 1990. “Analysis of time series subject to changes in regime.” \emph{Journal of econometrics}, 45 (1-2): 39–70.

Krolzig, Hans-Martin. 1997. “The markov-switching vector autoregressive model.”. Springer.
}
\seealso{
\code{\link{Nmdl}}
}
