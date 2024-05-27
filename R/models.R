#' @title Simulate normally distributed process
#' 
#' @description This function simulates a normally distributed process.
#' 
#' @param mdl_h0 List containing the following DGP parameters
#' \itemize{
#'   \item{n: }{Length of series.}
#'   \item{mu: }{A (\code{q x 1}) vector of means.}
#'   \item{sigma: }{A (\code{q x q}) covariance matrix.}
#'   \item{q: }{Number of series.}
#'   \item{eps: }{An optional (\code{T+burnin x q}) matrix with standard normal errors to be used. Errors will be generated if not provided.}
#' }
#' @param burnin Number of simulated observations to remove from beginning. Default is \code{100}.
#' 
#' @return List with simulated series and its DGP parameters.
#' 
#' @example /inst/examples/simuNorm_examples.R
#' @export
simuNorm <- function(mdl_h0, burnin = 0){
  simu_output <- simuNorm_cpp(mdl_h0, burnin)
  # Define class
  class(simu_output) <- "simuNorm"
  return(simu_output)
}


#' @title Simulate autoregressive process
#' 
#' @description This function simulates an autoregresive process.
#' 
#' @param mdl_h0 List containing the following DGP parameters
#' \itemize{
#'   \item{n: }{Length of series.}
#'   \item{mu: }{Mean of process.}
#'   \item{sigma: }{Standard deviation of process.}
#'   \item{phi: }{Vector of autoregressive coefficients.}
#'   \item{eps: }{An optional (\code{T+burnin x q}) matrix with standard normal errors to be used. Errors will be generated if not provided.}
#' }
#' @param burnin Number of simulated observations to remove from beginning. Default is \code{100}.
#' 
#' @return List with simulated autoregressive series and its DGP parameters.
#' 
#' @example /inst/examples/simuAR_examples.R
#' @export
simuAR <- function(mdl_h0, burnin = 100){
  simu_output <- simuAR_cpp(mdl_h0, burnin)  
  # Define class
  class(simu_output) <- "simuAR"
  return(simu_output)
}


#' @title Simulate VAR process
#' 
#' @description This function simulates a vector autoregresive process.
#' 
#' @param mdl_h0 List containing the following DGP parameters
#' \itemize{
#'   \item{n: }{Length of series.}
#'   \item{mu: }{A (\code{q x 1}) vector of means.}
#'   \item{sigma: }{A (\code{q x q}) covariance matrix.}
#'   \item{phi: }{ A (\code{q x qp}) matrix of autoregressive coefficients.}
#'   \item{p: }{Number of autoregressive lags.}
#'   \item{q: }{Number of series.}
#'   \item{eps: }{An optional (\code{T+burnin x q}) matrix with standard normal errors to be used. Errors will be generated if not provided.}
#' }
#' @param burnin Number of simulated observations to remove from beginning. Default is \code{100}.
#' 
#' @return List with simulated vector autoregressive series and its DGP parameters.
#' 
#' @export
simuVAR <- function(mdl_h0, burnin = 100){
  simu_output <- simuVAR_cpp(mdl_h0, burnin)
  # Define class
  class(simu_output) <- "simuVAR"
  return(simu_output)
}


#' @title Simulate Hidden Markov model with normally distributed errors
#' 
#' @description This function simulates a Hidden Markov Model process.
#' 
#' @param mdl_h0 List containing the following DGP parameters
#' \itemize{
#'   \item{n: }{Length of series.}
#'   \item{k: }{Number of regimes.}
#'   \item{mu: }{A (\code{k x q}) vector of means.}
#'   \item{sigma: }{A (\code{q x q}) covariance matrix.}
#'   \item{q: }{Number of series.}
#'   \item{P: }{A (\code{k x k}) transition matrix (columns must sum to one).}
#'   \item{eps: }{An optional (\code{T+burnin x q}) matrix with standard normal errors to be used. Errors will be generated if not provided.}
#' }
#' @param burnin Number of simulated observations to remove from beginning. Default is \code{100}.
#' 
#' @return List with simulated series and its DGP parameters.
#' 
#' @example /inst/examples/simuHMM_examples.R
#' @export
simuHMM <- function(mdl_h0, burnin = 100){
  simu_output <- simuHMM_cpp(mdl_h0, burnin)
  # Define class
  class(simu_output) <- "simuHMM"
  return(simu_output)
}


#' @title Simulate Markov-switching autoregressive process
#' 
#' @description This function simulates a Markov-switching autoregressive process.
#' 
#' @param mdl_h0 List containing the following DGP parameters
#' \itemize{
#'   \item{n: }{Length of series.}
#'   \item{k: }{Number of regimes.}
#'   \item{mu: }{A (\code{k x 1}) vector with mean of process in each regime.}
#'   \item{sigma: }{A (\code{k x 1}) vector with standard deviation of process in each regime.}
#'   \item{phi: }{Vector of autoregressive coefficients.}
#'   \item{P: }{A (\code{k x k}) transition matrix (columns must sum to one).}
#'   \item{eps: }{An optional (\code{T+burnin x q}) matrix with standard normal errors to be used. Errors will be generated if not provided.}
#' }
#' @param burnin Number of simulated observations to remove from beginning. Default is \code{100}.
#' 
#' @return List with simulated Markov-switching autoregressive process and its DGP properties.
#' 
#' @example /inst/examples/simuMSAR_examples.R
#' @export
simuMSAR <- function(mdl_h0, burnin = 100){
  simu_output <- simuMSAR_cpp(mdl_h0, burnin)
  # Define class
  class(simu_output) <- "simuMSAR"
  return(simu_output)
}


#' @title Simulate Markov-switching vector autoregressive process
#' 
#' @description This function simulates a Markov-switching vector autoregressive process.
#' 
#' @param mdl_h0 List containing the following DGP parameters
#' \itemize{
#'   \item{n: }{Length of series.}
#'   \item{k: }{Number of regimes.}
#'   \item{mu: }{A (\code{k x q}) matrix of means.}
#'   \item{sigma: }{List with \code{k} (\code{q x q}) covariance matrices.}
#'   \item{phi: }{A (\code{q x qp}) matrix of autoregressive coefficients.}
#'   \item{p: }{Number of autoregressive lags.}
#'   \item{q: }{Number of series.}
#'   \item{P: }{A (\code{k x k}) transition matrix (columns must sum to one).}
#'   \item{eps: }{An optional (\code{T+burnin x q}) matrix with standard normal errors to be used. Errors will be generated if not provided.}
#' }
#' @param burnin Number of simulated observations to remove from beginning. Default is \code{100}.
#' 
#' @return List with simulated vector autoregressive series and its DGP parameters.
#' 
#' @example /inst/examples/simuMSVAR_examples.R
#' @export
simuMSVAR <- function(mdl_h0, burnin = 100){
  simu_output <- simuMSVAR_cpp(mdl_h0, burnin)
  # Define class
  class(simu_output) <- "simuMSVAR"
  return(simu_output)
}


#' @title Normal distribution model
#' 
#' @description This function estimates a univariate or multivariate normally distributed model. This can be used for the null hypothesis of a linear model against an alternative hypothesis of a HMM with \code{k} regimes. 
#' 
#' @param Y a \code{(T x q)} matrix of observations. 
#' @param control List with model options including:
#' \itemize{
#'   \item{const: }{Boolean determining whether to estimate model with constant if \code{TRUE} or not if \code{FALSE}. Default is \code{TRUE}.}
#'   \item{getSE: }{Boolean determining whether to compute standard errors of parameters if \code{TRUE} or not if \code{FALSE}. Default is \code{TRUE}.}
#' }
#' 
#' @return List of class \code{Nmdl} (\code{S3} object) with model attributes including:
#' \itemize{
#'   \item{y: }{a \code{(T x q)} matrix of observations.}
#'   \item{fitted: }{a \code{(T x q)} matrix of fitted values.}
#'   \item{resid: }{a \code{(T x q)} matrix of residuals.}
#'   \item{mu: }{a \code{(1 x q)} vector of estimated means of each process.}
#'   \item{stdev: }{a \code{(q x 1)} vector of estimated standard deviation of each process.}
#'   \item{sigma: }{a \code{(q x q)} estimated covariance matrix.}
#'   \item{theta: }{vector containing: \code{mu} and \code{vech(sigma)}.}
#'   \item{theta_mu_ind: }{vector indicating location of mean with \code{1} and \code{0} otherwise.}
#'   \item{theta_sig_ind: }{vector indicating location of variance and covariances with \code{1} and \code{0} otherwise.}
#'   \item{theta_var_ind: }{vector indicating location of variances with \code{1} and \code{0} otherwise.}
#'   \item{n: }{number of observations (same as \code{T}).}
#'   \item{q: }{number of series.}
#'   \item{k: }{number of regimes. This is always \code{1} in \code{Nmdl}.}
#'   \item{control: }{List with model options used.}
#'   \item{logLike: }{log-likelihood.}
#'   \item{AIC: }{Akaike information criterion.}
#'   \item{BIC: }{Bayesian (Schwarz) information criterion.}
#'   \item{Hess: }{Hessian matrix. Approximated using \code{\link[numDeriv]{hessian}} and only returned if \code{getSE=TRUE}.}
#'   \item{info_mat: }{Information matrix. Computed as the inverse of \code{-Hess}. If matrix is not PD then nearest PD matrix is obtained using \code{\link[pracma]{nearest_spd}}. Only returned if \code{getSE=TRUE}.}
#'   \item{nearPD_used: }{Boolean determining whether \code{nearPD} function was used on \code{info_mat} if \code{TRUE} or not if \code{FALSE}. Only returned if \code{getSE=TRUE}.}
#'   \item{theta_se: }{standard errors of parameters in \code{theta}.  Only returned if \code{getSE=TRUE}.}
#' }
#' 
#' @example /inst/examples/Nmdl_examples.R
#' @export
Nmdl <- function(Y, control = list()){
  # ----- Set control values
  con <- list(const = TRUE, 
              getSE = TRUE)
  # Perform some checks for controls
  nmsC <- names(con)
  con[(namc <- names(control))] <- control
  if(length(noNms <- namc[!namc %in% nmsC])){
    warning("unknown names in control: ", paste(noNms,collapse=", ")) 
  }
  # ----- Get process dimensions
  n <- nrow(Y)
  q <- ncol(Y)
  # ----- estimate model
  if (con$const==TRUE){
    mu <- colMeans(Y)
  }else{
    mu <- rep(0, q)
  }
  # ----- Obtain variables of interest
  fitted      <- matrix(1, n, 1)%*%t(as.matrix(mu))
  resid       <- Y - fitted
  sigma       <- crossprod(resid)/(n-1)
  stdev       <- sqrt(diag(sigma))
  theta       <- c(mu,covar_vech(sigma))
  # theta indicators
  theta_mu_ind  <- c(rep(1, q), rep(0,length(theta)-q))
  theta_sig_ind <- c(rep(0, q), rep(1,q*(q+1)/2))
  theta_var_ind <- c(rep(0, q), t(covar_vech(diag(q))))
  # ----- Output
  out     <- list(y = Y, X = matrix(1, n, 1), fitted = fitted, resid = resid, mu = mu, stdev = stdev, sigma = sigma, theta = theta, 
                  theta_mu_ind = theta_mu_ind, theta_sig_ind = theta_sig_ind, theta_var_ind = theta_var_ind,
                  n = n, q = q, p = 0, k = 1, control = con)
  # Define class
  class(out) <- "Nmdl"
  # get log-likelihood
  out$logLike <- logLik(out)
  # get information criterion
  out$AIC <- stats::AIC(out)
  out$BIC <- stats::BIC(out)
  # names
  names(out$theta) <- c(paste0("mu_",(1:q)),
                        paste0("sig_",covar_vech(t(matrix(as.double(sapply((1:q),  function(x) paste0(x, (1:q)))), q,q)))))
  # get standard errors if 'con$getSE' is 'TRUE'
  if (con$getSE==TRUE){
    out <- thetaSE(out)
  }
  return(out)
}

#' @title Autoregressive Model
#' 
#' @description This function estimates an autoregresive model with \code{p} lags. This can be used for the null hypothesis of a linear model against an alternative hypothesis of a Markov switching autoregressive model with \code{k} regimes. 
#' 
#' @param Y A \code{(T x 1)} matrix of observations.  
#' @param p Integer determining the number of autoregressive lags.
#' @param control List with model options including: 
#' \itemize{
#'   \item{const: }{Boolean determining whether to estimate model with constant if \code{TRUE} or not if \code{FALSE}. Default is \code{TRUE}.}
#'   \item{getSE: }{Boolean determining whether to compute standard errors of parameters if \code{TRUE} or not if \code{FALSE}. Default is \code{TRUE}.}
#' }
#' 
#' @return List of class \code{ARmdl} (\code{S3} object) with model attributes including:
#' \itemize{
#'   \item{y: }{a \code{(T-p x 1)} matrix of observations.}
#'   \item{X: }{a \code{(T-p x p + const)} matrix of lagged observations with a leading column of \code{1}s if \code{const=TRUE} or not if \code{const=FALSE}.}
#'   \item{x: }{a \code{(T-p x p)} matrix of lagged observations.}
#'   \item{fitted: }{a \code{(T-p x 1)} matrix of fitted values.}
#'   \item{resid: }{a \code{(T-p x 1)} matrix of residuals.}
#'   \item{inter: }{estimated intercept of the process.}
#'   \item{mu: }{estimated mean of the process.}
#'   \item{coef: }{coefficient estimates. First value is the intercept (i.e., not \code{mu}) if \code{const=TRUE}. This is the same as \code{phi} if \code{const=FALSE}.}
#'   \item{intercept: }{estimate of intercept.}
#'   \item{phi: }{estimates of autoregressive coefficients.}
#'   \item{stdev: }{estimated standard deviation of the process.}
#'   \item{sigma: }{estimated variance of the process.}
#'   \item{theta: }{vector containing: \code{mu}, \code{sigma}, and \code{phi}.}
#'   \item{theta_mu_ind: }{vector indicating location of mean with \code{1} and \code{0} otherwise.}
#'   \item{theta_sig_ind: }{vector indicating location of variance with \code{1} and \code{0} otherwise.}
#'   \item{theta_var_ind: }{vector indicating location of variance with \code{1} and \code{0} otherwise. This is the same as \code{theta_sig_ind} in \code{ARmdl}.}
#'   \item{theta_phi_ind: }{vector indicating location of autoregressive coefficients with \code{1} and \code{0} otherwise.}
#'   \item{stationary: }{Boolean indicating if process is stationary if \code{TRUE} or non-stationary if \code{FALSE}.}
#'   \item{n: }{number of observations after lag transformation (i.e., \code{n = T-p}).}
#'   \item{p: }{number of autoregressive lags.}
#'   \item{q: }{number of series. This is always \code{1} in \code{ARmdl}.}
#'   \item{k: }{number of regimes. This is always \code{1} in \code{ARmdl}.}
#'   \item{control: }{List with model options used.}
#'   \item{logLike: }{log-likelihood.}
#'   \item{AIC: }{Akaike information criterion.}
#'   \item{BIC: }{Bayesian (Schwarz) information criterion.}
#'   \item{Hess: }{Hessian matrix. Approximated using \code{\link[numDeriv]{hessian}} and only returned if \code{getSE=TRUE}.}
#'   \item{info_mat: }{Information matrix. Computed as the inverse of \code{-Hess}. If matrix is not PD then nearest PD matrix is obtained using \code{\link[pracma]{nearest_spd}}. Only returned if \code{getSE=TRUE}.}
#'   \item{nearPD_used: }{Boolean determining whether \code{nearPD} function was used on \code{info_mat} if \code{TRUE} or not if \code{FALSE}. Only returned if \code{getSE=TRUE}.}
#'   \item{theta_se: }{standard errors of parameters in \code{theta}.  Only returned if \code{getSE=TRUE}.}
#' }
#' 
#' @seealso \code{\link{MSARmdl}}
#' @example /inst/examples/ARmdl_examples.R
#' @export
ARmdl <- function(Y, p, control = list()){
  # ----- Set control values
  con <- list(const = TRUE, 
              getSE = TRUE)
  # Perform some checks for controls
  nmsC <- names(con)
  con[(namc <- names(control))] <- control
  if(length(noNms <- namc[!namc %in% nmsC])){
    warning("unknown names in control: ", paste(noNms,collapse=", ")) 
  }
  # ----- transform data
  lagged_vals <- ts_lagged(Y, p)
  y <- lagged_vals$y
  x <- lagged_vals$X
  b0 <- matrix(0, p + con$const,1)
  phi <- matrix(0, p , 1)
  n = nrow(lagged_vals$y)
  npar = p
  # ----- estimate model
  if (con$const==TRUE){
    X     <- cbind(1, x)
    npar  <- npar + 1
    b0    <- crossprod(solve(crossprod(X)),crossprod(X,y))
    phi   <- b0[(2:npar),1]
    inter <- b0[1,1]
  }else{
    X     <- x
    b0    <- crossprod(solve(crossprod(X)),crossprod(X,y))
    phi   <- b0
    inter <- 0
  }
  # ----- Obtain variables of interest
  phisum  <- sum(phi)
  mu      <- inter/(1-phisum)
  fitted  <- X%*%b0
  resid   <- y - fitted
  stdev   <- sqrt((crossprod(resid))/(n-1))
  sigma   <- stdev^2
  theta   <- as.matrix(c(mu,sigma,phi))
  roots   <- all(Mod(as.complex(polyroot(c(1,-phi))))>1)
  #theta   <- as.matrix(c(mu,stdev,phi))
  # theta indicators
  theta_mu_ind  <- c(1,rep(0,length(theta)-1))
  theta_sig_ind <- c(0,1,rep(0,p))
  theta_var_ind <- c(0,1,rep(0,p))
  theta_phi_ind <- c(0,0,rep(1,p))
  # ----- Output
  out     <- list(y = y, X = X, x = x, fitted = fitted, resid = resid, mu = mu, inter = inter, coef = b0, intercept = inter, phi = phi,
                  stdev = stdev, sigma = sigma, theta = theta, theta_mu_ind = theta_mu_ind, 
                  theta_sig_ind = theta_sig_ind, theta_var_ind = theta_var_ind, theta_phi_ind = theta_phi_ind, 
                  stationary = roots, n = n, p = p, q = 1, k = 1, control = con)
  # Define class
  class(out) <- "ARmdl"
  # get log-likelihood
  out$logLike <- logLik(out)
  # get information criterion
  out$AIC <- stats::AIC(out)
  out$BIC <- stats::BIC(out)
  # names 
  names(out$theta) <- c(c("mu","sig"), paste0("phi_",(1:p)))
  # get standard errors if 'con$getSE' is 'TRUE'
  if (con$getSE==TRUE){
    out <- thetaSE(out)
  }
  return(out)
}

#' @title Vector autoregressive model
#' 
#' @description This function estimates a vector autoregresive model with \code{p} lags. This can be used for the null hypothesis of a linear model against an alternative hypothesis of a Markov switching vector autoregressive model with \code{k} regimes. 
#' 
#' @param Y a \code{(T x q)} matrix of observations.
#' @param p integer determining the number of autoregressive lags.
#' @param control List with model options including:
#' \itemize{
#' \item{const: }{Boolean determining whether to estimate model with constant if \code{TRUE} or not if \code{FALSE}. Default is \code{TRUE}.}
#' \item{getSE: }{Boolean determining whether to compute standard errors of parameters if \code{TRUE} or not if \code{FALSE}. Default is \code{TRUE}.}
#' }
#' 
#' @return List of class \code{VARmdl} (\code{S3} object) with model attributes including:
#' \itemize{
#'   \item{y: }{a \code{(T-p x q)} matrix of observations.}
#'   \item{X: }{a \code{(T-p x p*q + const)} matrix of lagged observations with a leading column of \code{1}s if \code{const=TRUE} or not if \code{const=FALSE}.}
#'   \item{x: }{a \code{(T-p x p*q)} matrix of lagged observations.}
#'   \item{fitted: }{a \code{(T-p x q)} matrix of fitted values}
#'   \item{resid: }{a \code{(T-p x q)} matrix of residuals.}
#'   \item{inter: }{a \code{(1 x q)} vector of estimated intercepts of each process.}
#'   \item{mu: }{a \code{(1 x q)} vector of estimated means of each process.}
#'   \item{coef: }{coefficient estimates. First row are the intercept (i.e., not \code{mu}) if \code{const=TRUE}. This is the same as \code{t(phi)} if \code{const=FALSE}.}
#'   \item{intercept: }{estimate of intercepts.}
#'   \item{phi: }{a \code{(q x p*q)} matrix of estimated autoregressive coefficients.}
#'   \item{Fmat: }{Companion matrix containing autoregressive coefficients.}
#'   \item{stdev: }{a \code{(q x 1)} vector of estimated standard deviation of each process.}
#'   \item{sigma: }{a \code{(q x q)} estimated covariance matrix.}
#'   \item{theta: }{vector containing: \code{mu}, \code{vech(sigma)}, and \code{vec(t(phi))}.}
#'   \item{theta_mu_ind: }{vector indicating location of mean with \code{1} and \code{0} otherwise.}
#'   \item{theta_sig_ind: }{vector indicating location of variance and covariances with \code{1} and \code{0} otherwise.}
#'   \item{theta_var_ind: }{vector indicating location of variances with \code{1} and \code{0} otherwise.}
#'   \item{theta_phi_ind: }{vector indicating location of autoregressive coefficients with \code{1} and \code{0} otherwise.}
#'   \item{stationary: }{Boolean indicating if process is stationary if \code{TRUE} or non-stationary if \code{FALSE}.}
#'   \item{n: }{number of observations after lag transformation (i.e., \code{n = T-p}).}
#'   \item{p: }{number of autoregressive lags.}
#'   \item{q: }{number of series.}
#'   \item{k: }{number of regimes. This is always \code{1} in \code{VARmdl}.}
#'   \item{Fmat: }{matrix from companion form. Used to determine is process is stationary.}
#'   \item{control: }{List with model options used.}
#'   \item{logLike: }{log-likelihood.}
#'   \item{AIC: }{Akaike information criterion.}
#'   \item{BIC: }{Bayesian (Schwarz) information criterion.}
#'   \item{Hess: }{Hessian matrix. Approximated using \code{\link[numDeriv]{hessian}} and only returned if \code{getSE=TRUE}.}
#'   \item{info_mat: }{Information matrix. Computed as the inverse of \code{-Hess}. If matrix is not PD then nearest PD matrix is obtained using \code{\link[pracma]{nearest_spd}}. Only returned if \code{getSE=TRUE}.}
#'   \item{nearPD_used: }{Boolean determining whether \code{nearPD} function was used on \code{info_mat} if \code{TRUE} or not if \code{FALSE}. Only returned if \code{getSE=TRUE}.}
#'   \item{theta_se: }{standard errors of parameters in \code{theta}.  Only returned if \code{getSE=TRUE}.}
#' }
#' 
#' @seealso \code{\link{MSVARmdl}}
#' @example /inst/examples/VARmdl_examples.R
#' @export
VARmdl <- function(Y, p, control = list()){
  # ----- Set control values
  con <- list(const = TRUE, 
              getSE = TRUE)
  # Perform some checks for controls
  nmsC <- names(con)
  con[(namc <- names(control))] <- control
  if(length(noNms <- namc[!namc %in% nmsC])){
    warning("unknown names in control: ", paste(noNms,collapse=", ")) 
  }
  # ----- transform data
  lagged_vals <- ts_lagged(Y, p)
  y <- lagged_vals$y
  x <- lagged_vals$X
  # ----- Get process dimensions
  n <- nrow(y)
  q <- ncol(y)
  # ----- estimate model
  if (con$const==TRUE){
    X     <- cbind(1,x)
    b0    <- crossprod(solve(crossprod(X)),crossprod(X,y))
    inter <- t(b0[1,])
    phi   <- t(b0[-1,])
  }else{
    X     <- x
    b0    <- crossprod(solve(crossprod(X)),crossprod(X,y))
    inter <- matrix(0,1,q)
    phi   <- t(b0)
  }
  # ----- Obtain variables of interest
  Fmat        <- companionMat(phi,p,q)
  stationary  <- all(abs(eigen(Fmat)[[1]])<1)
  mu_tmp      <- solve(diag(q*p)-Fmat)%*%as.matrix(c(inter,rep(0,q*(p-1))))
  mu          <- mu_tmp[(1:(q))]
  fitted      <- X%*%b0
  resid       <- y - fitted
  sigma       <- (crossprod(resid))/(n-1)
  stdev       <- sqrt(diag(sigma))
  theta       <- c(mu,covar_vech(sigma),c(t(phi)))
  # theta indicators
  theta_mu_ind  <- c(rep(1,q),rep(0,length(theta)-q))
  theta_sig_ind <- c(rep(0,q),rep(1,q*(q+1)/2),rep(0,q*q*p))
  theta_var_ind <- c(rep(0, q), t(covar_vech(diag(q))),rep(0,q*q*p))
  theta_phi_ind <- c(rep(0,length(theta)-q*q*p),rep(1,q*q*p))
  # ----- Output
  out     <- list(y = y, X = X, x = x, fitted = fitted, resid = resid, inter = inter, mu = mu, coef = b0, intercept = inter, phi = phi,
                  stdev = stdev, sigma = sigma, theta = theta, theta_mu_ind = theta_mu_ind, 
                  theta_sig_ind = theta_sig_ind, theta_var_ind = theta_var_ind, theta_phi_ind = theta_phi_ind, 
                  stationary = stationary, n = n, p = p, q = q, k = 1, Fmat = Fmat, control = con)
  out$Fmat    <- Fmat
  # Define class
  class(out) <- "VARmdl"
  # get log-likelihood
  out$logLike <- logLik(out)
  # get information criterion
  out$AIC <- stats::AIC(out)
  out$BIC <- stats::BIC(out)
  # names
  phi_n_tmp <- expand.grid((1:q),(1:p),(1:q))
  names(out$theta) <- c(paste0("mu_",(1:q)),
                        paste0("sig_",covar_vech(t(matrix(as.double(sapply((1:q),  function(x) paste0(x, (1:q)))), q,q)))),
                        paste0("phi_",paste0(phi_n_tmp[,2],",",phi_n_tmp[,3],phi_n_tmp[,1])))
  # get standard errors if 'con$getSE' is 'TRUE'
  if (con$getSE==TRUE){
    out <- thetaSE(out)
  }
  return(out)
}


#' @title Hidden Markov model 
#' 
#' @description This function estimates a Hidden Markov model with \code{k} regimes.
#' 
#' @param Y a \code{(T x q)} matrix of observations.
#' @param k integer determining the number of regimes to use in estimation. Must be greater than or equal to \code{2}. 
#' @param control List with model options including:
#' \itemize{
#'  \item{getSE: }{Boolean. If \code{TRUE} standard errors are computed and returned. If \code{FALSE} standard errors are not computed. Default is \code{TRUE}.}
#'  \item{msmu: }{Boolean. If \code{TRUE} model is estimated with switch in mean. If \code{FALSE} model is estimated with constant mean. Default is \code{TRUE}.}
#'  \item{msvar: }{Boolean. If \code{TRUE} model is estimated with switch in variance. If \code{FALSE} model is estimated with constant variance. Default is \code{TRUE}.}
#'  \item{init_theta: }{vector of initial values. vector must contain \code{(1 x q)} vector \code{mu}, \code{vech(sigma)}, and \code{vec(P)} where sigma is a \code{(q x q)} covariance matrix.This is optional. Default is \code{NULL}, in which case \code{\link{initVals_MSARmdl}} is used to generate initial values.}
#'  \item{method: }{string determining which method to use. Options are \code{'EM'} for EM algorithm or \code{'MLE'} for Maximum Likelihood Estimation.}
#'  \item{maxit: }{integer determining the maximum number of EM iterations.}
#'  \item{thtol: }{double determining the convergence criterion for the absolute difference in parameter estimates \code{theta} between iterations. Default is \code{1e-6}.}
#'  \item{maxit_converge: }{integer determining the maximum number of initial values attempted until solution is finite. For example, if parameters in \code{theta} or \code{logLike} are \code{NaN} another set of initial values (up to \code{maxit_converge}) is attempted until finite values are returned. This does not occur frequently for most types of data but may be useful in some cases. Once finite values are obtained, this counts as one iteration towards \code{use_diff_init}. Default is \code{500}.}
#'  \item{use_diff_init: }{integer determining how many different initial values to try (that do not return \code{NaN}; see \code{maxit_converge}). Default is \code{1}.}
#'  \item{mle_variance_constraint: }{double used to determine the lower bound on the smallest eigenvalue for the covariance matrix of each regime. Default is \code{1e-3}.}
#'  \item{mle_theta_low: }{Vector with lower bounds on parameters (Used only if method = "MLE"). Default is \code{NULL}.}
#'  \item{mle_theta_upp: }{Vector with upper bounds on parameters (Used only if method = "MLE"). Default is \code{NULL}.}
#' }
#' 
#' @return List of class \code{HMmdl} (\code{S3} object) with model attributes including:
#' \itemize{
#'   \item{y: }{a \code{(T x q)} matrix of observations.}
#'   \item{fitted: }{a \code{(T x q)} matrix of fitted values.}
#'   \item{resid: }{a \code{(T x q)} matrix of residuals.}
#'   \item{mu: }{a \code{(k x q)} matrix of estimated means of each process.}
#'   \item{stdev: }{List with \code{k} \code{(q x 1)} vector of estimated standard deviation of each process.}
#'   \item{sigma: }{List with \code{k} \code{(q x q)} estimated covariance matrix.}
#'   \item{theta: }{vector containing: \code{mu} and \code{vech(sigma)}.}
#'   \item{theta_mu_ind: }{vector indicating location of mean with \code{1} and \code{0} otherwise.}
#'   \item{theta_sig_ind: }{vector indicating location of variance and covariances with \code{1} and \code{0} otherwise.}
#'   \item{theta_var_ind: }{vector indicating location of variances with \code{1} and \code{0} otherwise.}
#'   \item{theta_var_ind: }{vector indicating location of transition matrix elements with \code{1} and \code{0} otherwise.}
#'   \item{n: }{number of observations (same as \code{T}).}
#'   \item{q: }{number of series.}
#'   \item{k: }{number of regimes in estimated model.}
#'   \item{P: }{a \code{(k x k)} transition matrix.}
#'   \item{pinf: }{a \code{(k x 1)} vector with limiting probabilities of each regime.}
#'   \item{St: }{a \code{(T x k)} vector with smoothed probabilities of each regime at each time \code{t}.}
#'   \item{deltath: }{double with maximum absolute difference in vector \code{theta} between last iteration.}
#'   \item{iterations: }{number of EM iterations performed to achieve convergence (if less than \code{maxit}).}
#'   \item{theta_0: }{vector of initial values used.}
#'   \item{init_used: }{number of different initial values used to get a finite solution. See description of input \code{maxit_converge}.}
#'   \item{msmu: }{Boolean. If \code{TRUE} model was estimated with switch in mean. If \code{FALSE} model was estimated with constant mean.}
#'   \item{msvar: }{Boolean. If \code{TRUE} model was estimated with switch in variance. If \code{FALSE} model was estimated with constant variance.}
#'   \item{control: }{List with model options used.}
#'   \item{logLike: }{log-likelihood.}
#'   \item{AIC: }{Akaike information criterion.}
#'   \item{BIC: }{Bayesian (Schwarz) information criterion.}
#'   \item{Hess: }{Hessian matrix. Approximated using \code{\link[numDeriv]{hessian}} and only returned if \code{getSE=TRUE}.}
#'   \item{info_mat: }{Information matrix. Computed as the inverse of \code{-Hess}. If matrix is not PD then nearest PD matrix is obtained using \code{\link[pracma]{nearest_spd}}. Only returned if \code{getSE=TRUE}.}
#'   \item{nearPD_used: }{Boolean determining whether \code{nearPD} function was used on \code{info_mat} if \code{TRUE} or not if \code{FALSE}. Only returned if \code{getSE=TRUE}.}
#'   \item{theta_se: }{standard errors of parameters in \code{theta}.  Only returned if \code{getSE=TRUE}.}
#'   \item{trace: }{List with Lists of estimation output for each initial value used due to \code{use_diff_init > 1}.}
#' }
#' 
#' @references Dempster, A. P., N. M. Laird, and D. B. Rubin. 1977. “Maximum Likelihood from Incomplete Data via the EM Algorithm.” \emph{Journal of the Royal Statistical Society}. Series B 39 (1): 1–38..
#' @references Hamilton, James D. 1990. “Analysis of time series subject to changes in regime.” \emph{Journal of econometrics}, 45 (1-2): 39–70.
#' @references Krolzig, Hans-Martin. 1997. “The markov-switching vector autoregressive model.”. Springer.
#' 
#' @seealso \code{\link{Nmdl}}
#' @export
HMmdl <- function(Y, k, control = list()){
  # ----- Set control values
  con <- list(getSE = TRUE,
              msmu = TRUE, 
              msvar = TRUE,
              init_theta = NULL,
              method = "EM",
              maxit = 1000,
              thtol = 1.e-6, 
              maxit_converge = 500, 
              use_diff_init = 1,
              mle_variance_constraint = 1e-3,
              mle_theta_low = NULL,
              mle_theta_upp = NULL)
  # ----- Perform some checks for controls
  nmsC <- names(con)
  con[(namc <- names(control))] <- control
  if(length(noNms <- namc[!namc %in% nmsC])){
    warning("unknown names in control: ", paste(noNms,collapse=", ")) 
  }
  if (k<2){
    stop("value for 'k' must be greater than or equal to 2.")
  }
  # ---------- Optimization options
  optim_options <- list(maxit = con$maxit, thtol = con$thtol)
  # pre-define list and matrix length for results
  output_all <- list(con$use_diff_init)
  max_loglik <- matrix(0, con$use_diff_init, 1)
  # ---------- Estimate linear model to use for initial values & transformed series
  init_control <- list(const = TRUE, getSE = FALSE)
  init_mdl <- Nmdl(Y, init_control)
  init_mdl$msmu <- con$msmu
  init_mdl$msvar <- con$msvar
  if (is.null(con$init_theta)==TRUE){
    if (con$method=="EM"){
      # ----- Estimate using 'use_diff_init' different initial values
      for (xi in 1:con$use_diff_init){
        init_used <- 0
        converge_check <- FALSE
        while ((converge_check==FALSE) & (init_used<con$maxit_converge)){
          # ----- Initial values
          theta_0 <- initVals_HMmdl(init_mdl, k)
          # ----- Estimate using EM algorithm and initial values provided
          output_tmp <- HMmdl_em(theta_0, init_mdl, k, optim_options)
          # ----- Convergence check
          logLike_tmp <- output_tmp$logLike
          theta_tmp <- output_tmp$theta
          converge_check <- ((is.finite(output_tmp$logLike)) & (all(is.finite(output_tmp$theta))))
          init_used <- init_used + 1
        }
        output_tmp$theta_0 <- theta_0
        max_loglik[xi] <- output_tmp$logLike
        output_tmp$init_used <- init_used
        output_all[[xi]] <- output_tmp
      }
    }else if (con$method=="MLE"){
      optim_options$mle_theta_low <- con$mle_theta_low
      optim_options$mle_theta_upp <- con$mle_theta_upp
      init_mdl$mle_variance_constraint <- con$mle_variance_constraint
      # ----- Estimate using 'use_diff_init' different initial values
      for (xi in 1:con$use_diff_init){
        init_used <- 0
        converge_check <- FALSE
        while ((converge_check==FALSE) & (init_used<con$maxit_converge)){
          # ----- Initial values
          theta_0 <- initVals_HMmdl(init_mdl, k)
          # ----- Estimate using roptim and initial values provided 
          output_tmp <- NULL
          try(
            output_tmp <- HMmdl_mle(theta_0, init_mdl, k, optim_options)  
          )
          # ----- Convergence check
          if (is.null(output_tmp)==FALSE){
            converge_check <- TRUE
          }else{
            converge_check <- FALSE
          }
          init_used = init_used + 1
        }
        output_tmp$theta_0 <- theta_0
        max_loglik[xi] <- output_tmp$logLike
        output_tmp$init_used <- init_used
        output_all[[xi]] <- output_tmp
      }
    }
    if (con$use_diff_init==1){
      output = output_tmp
    }else{
      xl = which.max(max_loglik)
      if (length(xl)==0){
        warning("Model(s) did not converge. Use higher 'use_diff_init' or 'maxit_converge'.")
        output <- output_all[[1]] 
      }else{
        output <- output_all[[xl]] 
      }
    }
  }else{
    if (con$method=="EM"){
      # ----- Estimate using EM algorithm and initial values provided
      output <- HMmdl_em(con$init_theta, init_mdl, k, optim_options)
      output$theta_0 <- con$init_theta
      output$init_used <- 1  
    }else if (con$method=="MLE"){
      init_mdl$mle_variance_constraint <- con$mle_variance_constraint
      # ----- Estimate using roptim and initial values provided 
      output <- HMmdl_mle(con$init_theta, init_mdl, k, optim_options)  
      output$theta_0 <- con$init_theta
      output$init_used <- 1  
    }
  }
  # ----- Obtain variables of interest
  q <- ncol(Y)
  Nsig <- (q*(q+1))/2
  theta_mu_ind <- c(rep(1, q + q*(k-1)*con$msmu), rep(0, Nsig + Nsig*(k-1)*con$msvar + k*k))
  theta_sig_ind <- c(rep(0, q + q*(k-1)*con$msmu), rep(1, Nsig + Nsig*(k-1)*con$msvar), rep(0, k*k))
  theta_var_ind <- c(rep(0, q + q*(k-1)*con$msmu), rep(t(covar_vech(diag(q))), 1+(k-1)*con$msvar), rep(0, k*k))
  theta_P_ind <- c(rep(0, q + q*(k-1)*con$msmu + Nsig + Nsig*(k-1)*con$msvar), rep(1, k*k))
  if (con$msmu==FALSE){
    output$mu <- matrix(output$mu, nrow=k, ncol=q, byrow=TRUE) 
  }
  out <- list(y = init_mdl$y, resid = output$resid, mu = output$mu, sigma = output$sigma, theta = output$theta, 
              theta_mu_ind = theta_mu_ind, theta_sig_ind = theta_sig_ind, theta_var_ind = theta_var_ind, theta_P_ind = theta_P_ind,
              n = init_mdl$n, q = q, p = 0, k = k, logLike = output$logLike, P = output$P, pinf = output$pinf, St = output$St,
              deltath = output$deltath,  iterations = output$iterations, theta_0 = output$theta_0,
              init_used = output$init_used, msmu = con$msmu, msvar = con$msvar, control = con)
  # Define class
  class(out) <- "HMmdl"
  # other output
  out$fitted  <- stats::fitted(out)
  out$AIC     <- stats::AIC(out)
  out$BIC     <- stats::BIC(out)
  stdev <- list(k)
  for (xk in 1:k){
    stdev[[xk]] <- diag(sqrt(diag(out$sigma[[xk]])))
  }
  out$stdev <- stdev
  # names
  if (q>1){
    mu_n_tmp <- expand.grid((1:q),(1:k))
    sig_n_tmp <- expand.grid(covar_vech(t(matrix(as.double(sapply((1:q),  function(x) paste0(x, (1:q)))), q,q))),(1:k))
    names(out$theta) <- c(if (con$msmu==TRUE) paste0("mu_",mu_n_tmp[,1],",",mu_n_tmp[,2]) else  paste0("mu_",(1:q)),
                          if (con$msvar==TRUE) paste0("sig_",sig_n_tmp[,1],",",sig_n_tmp[,2]) else paste0("sig_",covar_vech(t(matrix(as.double(sapply((1:q),  function(x) paste0(x, (1:q)))), q,q)))),
                          paste0("p_",c(sapply((1:k),  function(x) paste0(x, (1:k)) ))))  
  }else if (q==1){
    names(out$theta) <- c(if (con$msmu==TRUE) paste0("mu_", (1:k)) else "mu",
                          if (con$msvar==TRUE) paste0("sig_", (1:k)) else "sig",
                          paste0("p_",c(sapply((1:k),  function(x) paste0(x, (1:k)) ))))
  }
  if (con$getSE==TRUE){
    out <- thetaSE(out)
  }
  if (is.null(con$init_theta)){
    out$trace <- output_all
  }
  return(out)
}



#' @title Markov-switching autoregressive model 
#' 
#' @description This function estimates a Markov-switching autoregressive model
#' 
#' @param Y (T x 1) vector with observational data. 
#' @param p integer for the number of lags to use in estimation. Must be greater than or equal to \code{1}. 
#' @param k integer for the number of regimes to use in estimation. Must be greater than or equal to \code{2}.
#' @param control List with model options including:
#' \itemize{
#'  \item{getSE: }{Boolean. If \code{TRUE} standard errors are computed and returned. If \code{FALSE} standard errors are not computed. Default is \code{TRUE}.}
#'  \item{msmu: }{Boolean. If \code{TRUE} model is estimated with switch in mean. If \code{FALSE} model is estimated with constant mean. Default is \code{TRUE}.}
#'  \item{msvar: }{Boolean. If \code{TRUE} model is estimated with switch in variance. If \code{FALSE} model is estimated with constant variance. Default is \code{TRUE}.}
#'  \item{init_theta: }{vector of initial values. vector must contain \code{(1 x q)} vector \code{mu}, \code{vech(sigma)}, and \code{vec(P)} where sigma is a \code{(q x q)} covariance matrix.This is optional. Default is \code{NULL}, in which case \code{\link{initVals_MSARmdl}} is used to generate initial values.}
#'  \item{method: }{string determining which method to use. Options are \code{'EM'} for EM algorithm or \code{'MLE'} for Maximum Likelihood Estimation.}
#'  \item{maxit: }{integer determining the maximum number of EM iterations.}
#'  \item{thtol: }{double determining the convergence criterion for the absolute difference in parameter estimates \code{theta} between iterations. Default is \code{1e-6}.}
#'  \item{maxit_converge: }{integer determining the maximum number of initial values attempted until solution is finite. For example, if parameters in \code{theta} or \code{logLike} are \code{NaN} another set of initial values (up to \code{maxit_converge}) is attempted until finite values are returned. This does not occur frequently for most types of data but may be useful in some cases. Once finite values are obtained, this counts as one iteration towards \code{use_diff_init}. Default is \code{500}.}
#'  \item{use_diff_init: }{integer determining how many different initial values to try (that do not return \code{NaN}; see \code{maxit_converge}). Default is \code{1}.}
#'  \item{mle_stationary_constraint: }{Boolean determining if only stationary solutions are considered (if \code{TRUE}) or not (if \code{FALSE}). Default is \code{TRUE}.}
#'  \item{mle_variance_constraint: }{Double used to determine the lower bound for variance in each regime. Value should be between \code{0} and \code{1} as it is multiplied by single regime variance. Default is \code{0.01} (i.e., \code{1\%} of single regime variance.}
#'  \item{mle_theta_low: }{Vector with lower bounds on parameters (Used only if method = "MLE"). Default is \code{NULL}.}
#'  \item{mle_theta_upp: }{Vector with upper bounds on parameters (Used only if method = "MLE"). Default is \code{NULL}.}
#' }
#' 
#' @return List of class \code{MSARmdl} (\code{S3} object) with model attributes including:
#' \itemize{
#'   \item{y: }{a \code{(T x 1)} matrix of observations.}
#'   \item{X: }{a \code{(T-p x p + const)} matrix of lagged observations with a leading column of \code{1}s.}
#'   \item{x: }{a \code{(T-p x p)} matrix of lagged observations.}
#'   \item{fitted: }{a \code{(T x 1)} matrix of fitted values.}
#'   \item{resid: }{a \code{(T x 1)} matrix of residuals.}
#'   \item{inter: }{a \code{(k x 1)} vector of estimated intercepts of each process.}
#'   \item{mu: }{a \code{(k x 1)} vector of estimated means of each process.}
#'   \item{phi: }{estimates of autoregressive coefficients.}
#'   \item{stdev: }{a \code{(k x 1)} vector of estimated standard deviation of each process.}
#'   \item{sigma: }{a \code{(k x 1)} estimated covariance matrix.}
#'   \item{theta: }{vector containing: \code{mu} and \code{vech(sigma)}.}
#'   \item{theta_mu_ind: }{vector indicating location of mean with \code{1} and \code{0} otherwise.}
#'   \item{theta_sig_ind: }{vector indicating location of variances with \code{1} and \code{0} otherwise.}
#'   \item{theta_var_ind: }{vector indicating location of variances with \code{1} and \code{0} otherwise. This is the same as \code{theta_sig_ind} in \code{MSARmdl}}
#'   \item{theta_P_ind: }{vector indicating location of transition matrix elements with \code{1} and \code{0} otherwise.}
#'   \item{stationary: }{Boolean indicating if process is stationary if \code{TRUE} or non-stationary if \code{FALSE}.}
#'   \item{n: }{number of observations (same as \code{T}).}
#'   \item{p: }{number of autoregressive lags.}
#'   \item{q: }{number of series. This is always \code{1} in \code{MSARmdl}.}
#'   \item{k: }{number of regimes in estimated model.}
#'   \item{P: }{a \code{(k x k)} transition matrix.}
#'   \item{pinf: }{a \code{(k x 1)} vector with limiting probabilities of each regime.}
#'   \item{St: }{a \code{(T x k)} vector with smoothed probabilities of each regime at each time \code{t}.}
#'   \item{deltath: }{double with maximum absolute difference in vector \code{theta} between last iteration.}
#'   \item{iterations: }{number of EM iterations performed to achieve convergence (if less than \code{maxit}).}
#'   \item{theta_0: }{vector of initial values used.}
#'   \item{init_used: }{number of different initial values used to get a finite solution. See description of input \code{maxit_converge}.}
#'   \item{msmu: }{Boolean. If \code{TRUE} model was estimated with switch in mean. If \code{FALSE} model was estimated with constant mean.}
#'   \item{msvar: }{Boolean. If \code{TRUE} model was estimated with switch in variance. If \code{FALSE} model was estimated with constant variance.}
#'   \item{control: }{List with model options used.}
#'   \item{logLike: }{log-likelihood.}
#'   \item{AIC: }{Akaike information criterion.}
#'   \item{BIC: }{Bayesian (Schwarz) information criterion.}
#'   \item{Hess: }{Hessian matrix. Approximated using \code{\link[numDeriv]{hessian}} and only returned if \code{getSE=TRUE}.}
#'   \item{info_mat: }{Information matrix. Computed as the inverse of \code{-Hess}. If matrix is not PD then nearest PD matrix is obtained using \code{\link[pracma]{nearest_spd}}. Only returned if \code{getSE=TRUE}.}
#'   \item{nearPD_used: }{Boolean determining whether \code{nearPD} function was used on \code{info_mat} if \code{TRUE} or not if \code{FALSE}. Only returned if \code{getSE=TRUE}.}
#'   \item{theta_se: }{standard errors of parameters in \code{theta}.  Only returned if \code{getSE=TRUE}.}
#'   \item{trace: }{List with Lists of estimation output for each initial value used due to \code{use_diff_init > 1}.}
#' }
#' 
#' @references Dempster, A. P., N. M. Laird, and D. B. Rubin. 1977. “Maximum Likelihood from Incomplete Data via the EM Algorithm.” \emph{Journal of the Royal Statistical Society}. Series B 39 (1): 1–38..
#' @references Hamilton, James D. 1990. “Analysis of time series subject to changes in regime.” \emph{Journal of econometrics}, 45 (1-2): 39–70.
#' 
#' @seealso \code{\link{ARmdl}}
#' @example /inst/examples/MSARmdl_examples.R
#' @export
MSARmdl <- function(Y, p, k, control = list()){
  # ----- Set control values
  con <- list(getSE = TRUE,
              msmu = TRUE, 
              msvar = TRUE,
              init_theta = NULL,
              method = "EM",
              maxit = 1000, 
              thtol = 1.e-6, 
              maxit_converge = 500, 
              use_diff_init = 1, 
              mle_stationary_constraint = TRUE,
              mle_variance_constraint = 0.01,
              mle_theta_low = NULL,
              mle_theta_upp = NULL)
  # ----- Perform some checks for controls
  nmsC <- names(con)
  con[(namc <- names(control))] <- control
  if(length(noNms <- namc[!namc %in% nmsC])){
    warning("unknown names in control: ", paste(noNms,collapse=", ")) 
  }
  if (k<2){
    stop("value for 'k' must be greater than or equal to 2.")
  }
  # ---------- Optimization options
  optim_options <- list(maxit = con$maxit, thtol = con$thtol)
  # pre-define list and matrix length for results
  output_all <- list(con$use_diff_init)
  max_loglik <- matrix(0, con$use_diff_init, 1)
  # ---------- Estimate linear model to use for initial values & transformed series
  init_control <- list(const = TRUE, getSE = FALSE)
  init_mdl <- ARmdl(Y, p = p, init_control)
  init_mdl$msmu <- con$msmu
  init_mdl$msvar <- con$msvar
  if (is.null(con$init_theta)==TRUE){
    if (con$method=="EM"){
      # ----- Estimate using 'use_diff_init' different initial values
      for (xi in 1:con$use_diff_init){
        init_used <- 0
        converge_check <- FALSE
        while ((converge_check==FALSE) & (init_used<con$maxit_converge)){
          # ----- Initial values
          theta_0 <- initVals_MSARmdl(init_mdl, k)  
          # ----- Estimate using EM algorithm and initial values provided
          output_tmp <- MSARmdl_em(theta_0, init_mdl, k, optim_options)
          # ----- Convergence check
          logLike_tmp <- output_tmp$logLike
          theta_tmp <- output_tmp$theta
          converge_check <- ((is.finite(output_tmp$logLike)) & (all(is.finite(output_tmp$theta))))
          init_used <- init_used + 1
        }
        output_tmp$theta_0 <- theta_0
        max_loglik[xi] <- output_tmp$logLike
        output_tmp$init_used <- init_used
        output_all[[xi]] <- output_tmp
      }
    }else if (con$method=="MLE"){
      optim_options$mle_theta_low <- con$mle_theta_low
      optim_options$mle_theta_upp <- con$mle_theta_upp
      init_mdl$mle_stationary_constraint <- con$mle_stationary_constraint
      init_mdl$mle_variance_constraint <- con$mle_variance_constraint
      # ----- Estimate using 'use_diff_init' different initial values
      for (xi in 1:con$use_diff_init){
        init_used <- 0
        converge_check <- FALSE
        while ((converge_check==FALSE) & (init_used<con$maxit_converge)){
          # ----- Initial values
          theta_0 <- initVals_MSARmdl(init_mdl, k)  
          # ----- Estimate using roptim and initial values provided 
          output_tmp <- NULL
          try(
            output_tmp <- MSARmdl_mle(theta_0, init_mdl, k, optim_options)  
          )
          # ----- Convergence check
          if (is.null(output_tmp)==FALSE){
            converge_check <- TRUE
          }else{
            converge_check <- FALSE
          }
          init_used = init_used + 1 
        }
        output_tmp$theta_0 <- theta_0
        max_loglik[xi] <- output_tmp$logLike
        output_tmp$init_used <- init_used
        output_all[[xi]] <- output_tmp
      }
    }
    if (con$use_diff_init==1){
      output <- output_tmp
    }else{
      xl = which.max(max_loglik)
      if (length(xl)==0){
        warning("Model(s) did not converge. Use higher 'use_diff_init' or 'maxit_converge'.")
        output <- output_all[[1]] 
      }else{
        output <- output_all[[xl]] 
      }
    }
  }else{
    if (con$method=="EM"){
      # ----- Estimate using EM algorithm and initial values provided
      output <- MSARmdl_em(con$init_theta, init_mdl, k, optim_options)
      output$theta_0 <- con$init_theta
      output$init_used <- 1  
    }else if (con$method=="MLE"){
      init_mdl$mle_stationary_constraint <- con$mle_stationary_constraint
      init_mdl$mle_variance_constraint <- con$mle_variance_constraint
      # ----- Estimate using roptim and initial values provided 
      output <- MSARmdl_mle(con$init_theta, init_mdl, k, optim_options)  
      output$theta_0 <- con$init_theta
      output$init_used <- 1  
    }
  }
  roots <- NULL
  try(
    roots   <- all(Mod(as.complex(polyroot(c(1,-output$phi))))>1)  
  )
  # ----- Obtain variables of interest
  theta_mu_ind <- c(rep(1, 1 + (k-1)*con$msmu), rep(0, 1 + (k-1)*con$msvar + p + k*k))
  theta_sig_ind <- c(rep(0, 1 + (k-1)*con$msmu), rep(1, 1 + (k-1)*con$msvar), rep(0, p + k*k))
  theta_var_ind <- theta_sig_ind
  theta_phi_ind <- c(rep(0, 2 + (k-1)*con$msmu + (k-1)*con$msvar), rep(1, p), rep(0, k*k))
  theta_P_ind <- c(rep(0, 2 + (k-1)*con$msmu + (k-1)*con$msvar + p), rep(1, k*k))
  # ----- Output
  if (con$msmu==FALSE){
    output$mu <- matrix(output$mu, nrow=k, ncol=1, byrow=TRUE) 
  }
  if (con$msvar==FALSE){
    output$sigma <- matrix(output$sigma, nrow=k, ncol=1, byrow=TRUE) 
  }
  out <- list(y = init_mdl$y, X = init_mdl$X, x = init_mdl$x, resid = output$resid, mu = output$mu, 
              phi = output$phi, stdev = sqrt(output$sigma), sigma = output$sigma, 
              theta = output$theta, theta_mu_ind = theta_mu_ind, theta_sig_ind = theta_sig_ind, theta_var_ind = theta_var_ind, 
              theta_phi_ind = theta_phi_ind, theta_P_ind = theta_P_ind, stationary = roots, 
              n = init_mdl$n, p = p, q = 1, k = k, logLike = output$logLike, P = output$P, pinf = output$pinf, 
              St = output$St, deltath = output$deltath, iterations = output$iterations, theta_0 = output$theta_0, 
              init_used = output$init_used, msmu = con$msmu, msvar = con$msvar, control = con)
  # Define class
  class(out) <- "MSARmdl"
  # other output
  out$inter   <- interMSARmdl(out)
  out$fitted  <- stats::fitted(out)
  out$AIC     <- stats::AIC(out)
  out$BIC     <- stats::BIC(out)
  # names 
  names(out$theta) <- c(if (con$msmu==TRUE) paste0("mu_", (1:k)) else "mu",
                        if (con$msvar==TRUE) paste0("sig_", (1:k)) else "sig",
                        paste0("phi_",(1:p)),
                        paste0("p_",c(sapply((1:k),  function(x) paste0(x, (1:k)) ))))
  if (con$getSE==TRUE){
    out <- thetaSE(out)
  }
  if (is.null(con$init_theta)){
    out$trace <- output_all
  }
  return(out)
}



#' @title Markov-switching vector autoregressive model 
#' 
#' @description This function estimates a Markov-switching vector autoregressive model 
#' 
#' @param Y (\code{T x q}) vector with observational data.
#' @param p integer for the number of lags to use in estimation. Must be greater than or equal to \code{0}.
#' @param k integer for the number of regimes to use in estimation. Must be greater than or equal to \code{2}.
#' @param control List with optimization options including:
#' \itemize{
#'  \item{getSE: }{Boolean. If \code{TRUE} standard errors are computed and returned. If \code{FALSE} standard errors are not computed. Default is \code{TRUE}.}
#'  \item{msmu: }{Boolean. If \code{TRUE} model is estimated with switch in mean. If \code{FALSE} model is estimated with constant mean. Default is \code{TRUE}.}
#'  \item{msvar: }{Boolean. If \code{TRUE} model is estimated with switch in variance. If \code{FALSE} model is estimated with constant variance. Default is \code{TRUE}.}
#'  \item{init_theta: }{vector of initial values. vector must contain \code{(1 x q)} vector \code{mu}, \code{vech(sigma)}, and \code{vec(P)} where sigma is a \code{(q x q)} covariance matrix. This is optional. Default is \code{NULL}, in which case \code{\link{initVals_MSARmdl}} is used to generate initial values.}
#'  \item{method: }{string determining which method to use. Options are \code{'EM'} for EM algorithm or \code{'MLE'} for Maximum Likelihood Estimation.}
#'  \item{maxit: }{integer determining the maximum number of EM iterations.}
#'  \item{thtol: }{double determining the convergence criterion for the absolute difference in parameter estimates \code{theta} between iterations. Default is \code{1e-6}.}
#'  \item{maxit_converge: }{integer determining the maximum number of initial values attempted until solution is finite. For example, if parameters in \code{theta} or \code{logLike} are \code{NaN} another set of initial values (up to \code{maxit_converge}) is attempted until finite values are returned. This does not occur frequently for most types of data but may be useful in some cases. Once finite values are obtained, this counts as one iteration towards \code{use_diff_init}. Default is \code{500}.}
#'  \item{use_diff_init: }{integer determining how many different initial values to try (that do not return \code{NaN}; see \code{maxit_converge}). Default is \code{1}.}
#'  \item{mle_stationary_constraint: }{Boolean determining if only stationary solutions are considered (if \code{TRUE}) or not (if \code{FALSE}). Default is \code{TRUE}.}
#'  \item{mle_variance_constraint: }{double used to determine the lower bound on the smallest eigenvalue for the covariance matrix of each regime. Default is \code{1e-3}.}
#'  \item{mle_theta_low: }{Vector with lower bounds on parameters (Used only if method = "MLE"). Default is \code{NULL}.}
#'  \item{mle_theta_upp: }{Vector with upper bounds on parameters (Used only if method = "MLE"). Default is \code{NULL}.}
#' }
#' 
#' @return List of class \code{MSVARmdl} (\code{S3} object) with model attributes including:
#' \itemize{
#'   \item{y: }{a \code{(T-p x q)} matrix of observations.}
#'   \item{X: }{a \code{(T-p x p*q + const)} matrix of lagged observations with a leading column of \code{1}s.}
#'   \item{x: }{a \code{(T-p x p*q)} matrix of lagged observations.}
#'   \item{resid: }{a \code{(T-p x q)} matrix of residuals.}
#'   \item{inter: }{a \code{(k x q)} matrix of estimated intercepts of each process.}
#'   \item{mu: }{a \code{(k x q)} matrix of estimated means of each process.}
#'   \item{phi: }{estimates of autoregressive coefficients.}
#'   \item{Fmat: }{Companion matrix containing autoregressive coefficients.}
#'   \item{stdev: }{List with \code{k} \code{(q x q)} matrices with estimated standard deviation on the diagonal.}
#'   \item{sigma: }{List with \code{k} \code{(q x q)} matrices with estimated covariance matrix.}
#'   \item{theta: }{vector containing: \code{mu} and \code{vech(sigma)}.}
#'   \item{theta_mu_ind: }{vector indicating location of mean with \code{1} and \code{0} otherwise.}
#'   \item{theta_sig_ind: }{vector indicating location of variance and covariances with \code{1} and \code{0} otherwise.}
#'   \item{theta_var_ind: }{vector indicating location of variances with \code{1} and \code{0} otherwise.}
#'   \item{theta_P_ind: }{vector indicating location of transition matrix elements with \code{1} and \code{0} otherwise.}
#'   \item{stationary: }{Boolean indicating if process is stationary if \code{TRUE} or non-stationary if \code{FALSE}.}
#'   \item{n: }{number of observations (same as \code{T}).}
#'   \item{p: }{number of autoregressive lags.}
#'   \item{q: }{number of series.}
#'   \item{k: }{number of regimes in estimated model.}
#'   \item{P: }{a \code{(k x k)} transition matrix.}
#'   \item{pinf: }{a \code{(k x 1)} vector with limiting probabilities of each regime.}
#'   \item{St: }{a \code{(T x k)} vector with smoothed probabilities of each regime at each time \code{t}.}
#'   \item{deltath: }{double with maximum absolute difference in vector \code{theta} between last iteration.}
#'   \item{iterations: }{number of EM iterations performed to achieve convergence (if less than \code{maxit}).}
#'   \item{theta_0: }{vector of initial values used.}
#'   \item{init_used: }{number of different initial values used to get a finite solution. See description of input \code{maxit_converge}.}
#'   \item{msmu: }{Boolean. If \code{TRUE} model was estimated with switch in mean. If \code{FALSE} model was estimated with constant mean.}
#'   \item{msvar: }{Boolean. If \code{TRUE} model was estimated with switch in variance. If \code{FALSE} model was estimated with constant variance.}
#'   \item{control: }{List with model options used.}
#'   \item{logLike: }{log-likelihood.}
#'   \item{AIC: }{Akaike information criterion.}
#'   \item{BIC: }{Bayesian (Schwarz) information criterion.}
#'   \item{Hess: }{Hessian matrix. Approximated using \code{\link[numDeriv]{hessian}} and only returned if \code{getSE=TRUE}.}
#'   \item{info_mat: }{Information matrix. Computed as the inverse of \code{-Hess}. If matrix is not PD then nearest PD matrix is obtained using \code{\link[pracma]{nearest_spd}}. Only returned if \code{getSE=TRUE}.}
#'   \item{nearPD_used: }{Boolean determining whether \code{nearPD} function was used on \code{info_mat} if \code{TRUE} or not if \code{FALSE}. Only returned if \code{getSE=TRUE}.}
#'   \item{theta_se: }{standard errors of parameters in \code{theta}.  Only returned if \code{getSE=TRUE}.}
#'   \item{trace: }{List with Lists of estimation output for each initial value used due to \code{use_diff_init > 1}.}
#' }
#' 
#' @return List with model characteristics
#' 
#' @references Dempster, A. P., N. M. Laird, and D. B. Rubin. 1977. “Maximum Likelihood from Incomplete Data via the EM Algorithm.” \emph{Journal of the Royal Statistical Society}. Series B 39 (1): 1–38..
#' @references Krolzig, Hans-Martin. 1997. “The markov-switching vector autoregressive model.”. Springer.
#'  
#' @seealso \code{\link{VARmdl}}
#' @example /inst/examples/MSVARmdl_examples.R
#' @export
MSVARmdl <- function(Y, p, k, control = list()){
  # ----- Set control values
  con <- list(getSE = TRUE,
              msmu = TRUE, 
              msvar = TRUE,
              init_theta = NULL,
              method = "EM",
              maxit = 1000,
              thtol = 1.e-6, 
              maxit_converge = 500, 
              use_diff_init = 1,
              mle_stationary_constraint = TRUE,
              mle_variance_constraint = 1e-3,
              mle_theta_low = NULL,
              mle_theta_upp = NULL)
  # ----- Perform some checks for controls
  nmsC <- names(con)
  con[(namc <- names(control))] <- control
  if(length(noNms <- namc[!namc %in% nmsC])){
    warning("unknown names in control: ", paste(noNms,collapse=", ")) 
  }
  if (k<2){
    stop("value for 'k' must be greater than or equal to 2.")
  }
  # ---------- Optimization options
  optim_options <- list(maxit = con$maxit, thtol = con$thtol)
  # pre-define list and matrix length for results
  output_all <- list(con$use_diff_init)
  max_loglik <- matrix(0, con$use_diff_init, 1)
  # ---------- Estimate linear model to use for initial values & transformed series
  init_control <- list(const = TRUE, getSE = FALSE)
  init_mdl <- VARmdl(Y, p = p, init_control)
  init_mdl$msmu <- con$msmu
  init_mdl$msvar <- con$msvar
  if (is.null(con$init_theta)==TRUE){
    if (con$method=="EM"){
      # ----- Estimate using 'use_diff_init' different initial values
      for (xi in 1:con$use_diff_init){
        init_used <- 0
        converge_check <- FALSE
        while ((converge_check==FALSE) & (init_used<con$maxit_converge)){
          # ----- Initial values
          theta_0 <- initVals_MSVARmdl(init_mdl, k)
          # ----- Estimate using EM algorithm and initial values provided
          output_tmp <- MSVARmdl_em(theta_0, init_mdl, k, optim_options)
          # ----- Convergence check
          logLike_tmp <- output_tmp$logLike
          theta_tmp <- output_tmp$theta
          converge_check <- ((is.finite(output_tmp$logLike)) & (all(is.finite(output_tmp$theta))))
          init_used <- init_used + 1
        }
        output_tmp$theta_0 <- theta_0
        max_loglik[xi] <- output_tmp$logLike
        output_tmp$init_used <- init_used
        output_all[[xi]] <- output_tmp
      }
    }else if (con$method=="MLE"){
      optim_options$mle_theta_low <- con$mle_theta_low
      optim_options$mle_theta_upp <- con$mle_theta_upp
      init_mdl$mle_stationary_constraint <- con$mle_stationary_constraint
      init_mdl$mle_variance_constraint <- con$mle_variance_constraint
      # ----- Estimate using 'use_diff_init' different initial values
      for (xi in 1:con$use_diff_init){
        init_used <- 0
        converge_check <- FALSE
        while ((converge_check==FALSE) & (init_used<con$maxit_converge)){
          # ----- Initial values
          theta_0 <- initVals_MSVARmdl(init_mdl, k)
          # ----- Estimate using roptim and initial values provided 
          output_tmp <- NULL
          try(
            output_tmp <- MSVARmdl_mle(theta_0, init_mdl, k, optim_options)  
          )
          # ----- Convergence check
          if (is.null(output_tmp)==FALSE){
            converge_check <- TRUE
          }else{
            converge_check <- FALSE
          }
          init_used = init_used + 1
        }
        output_tmp$theta_0 <- theta_0
        max_loglik[xi] <- output_tmp$logLike
        output_tmp$init_used <- init_used
        output_all[[xi]] <- output_tmp
      }
    }
    if (con$use_diff_init==1){
      output = output_tmp
    }else{
      xl = which.max(max_loglik)
      if (length(xl)==0){
        warning("Model(s) did not converge. Use higher 'use_diff_init' or 'maxit_converge'.")
        output <- output_all[[1]] 
      }else{
        output <- output_all[[xl]] 
      }
    }
  }else{
    if (con$method=="EM"){
      # ----- Estimate using EM algorithm and initial values provided
      output <- MSVARmdl_em(con$init_theta, init_mdl, k, optim_options)
      output$theta_0 <- con$init_theta
      output$init_used <- 1  
    }else if (con$method=="MLE"){
      optim_options$lower <- con$lower
      optim_options$upper <- con$upper
      init_mdl$mle_stationary_constraint <- con$mle_stationary_constraint
      init_mdl$mle_variance_constraint <- con$mle_variance_constraint
      # ----- Estimate using roptim and initial values provided 
      output <- MSVARmdl_mle(con$init_theta, init_mdl, k, optim_options)  
      output$theta_0 <- con$init_theta
      output$init_used <- 1  
    }
  }
  # ----- Obtain variables of interest
  q <- ncol(Y)
  Nsig <- (q*(q+1))/2
  phi_len <- q*p*q 
  theta_mu_ind <- c(rep(1, q + q*(k-1)*con$msmu), rep(0, Nsig + Nsig*(k-1)*con$msvar + phi_len + k*k))
  theta_sig_ind <- c(rep(0, q + q*(k-1)*con$msmu), rep(1, Nsig + Nsig*(k-1)*con$msvar), rep(0, phi_len + k*k))
  theta_var_ind <- c(rep(0, q + q*(k-1)*con$msmu), rep(t(covar_vech(diag(q))), 1+(k-1)*con$msvar), rep(0, phi_len + k*k))
  theta_phi_ind <- c(rep(0, q + q*(k-1)*con$msmu + Nsig + Nsig*(k-1)*con$msvar), rep(1, phi_len), rep(0, k*k))
  theta_P_ind <- c(rep(0, q + q*(k-1)*con$msmu + Nsig + Nsig*(k-1)*con$msvar + phi_len), rep(1, k*k))
  if (con$msmu==FALSE){
    output$mu <- matrix(output$mu, nrow=k, ncol=q, byrow=TRUE) 
  }
  out <- list(y = init_mdl$y, X = init_mdl$X, x = init_mdl$x, resid = output$resid, mu = output$mu, phi = output$phi,
              sigma = output$sigma, theta = output$theta, theta_mu_ind = theta_mu_ind, theta_sig_ind = theta_sig_ind, theta_var_ind = theta_var_ind, 
              theta_phi_ind = theta_phi_ind, theta_P_ind = theta_P_ind, stationary = NULL, n = init_mdl$n, p = p, q = q, k = k, logLike = output$logLike, P = output$P, pinf = output$pinf, 
              St = output$St, deltath = output$deltath,  iterations = output$iterations, theta_0 = output$theta_0,
              init_used = output$init_used, msmu = con$msmu, msvar = con$msvar, control = con)
  # Define class
  class(out) <- "MSVARmdl"
  # other output
  out$Fmat    <- companionMat(out$phi, out$p, out$q)
  out$inter   <- interMSVARmdl(out)
  out$fitted  <- stats::fitted(out)
  out$AIC     <- stats::AIC(out)
  out$BIC     <- stats::BIC(out)
  stdev <- list(k)
  for (xk in 1:k){
    stdev[[xk]] <- diag(sqrt(diag(out$sigma[[xk]])))
  }
  out$stdev <- stdev
  # names
  phi_n_tmp <- expand.grid((1:q),(1:p),(1:q))
  mu_n_tmp <- expand.grid((1:q),(1:k))
  sig_n_tmp <- expand.grid(covar_vech(t(matrix(as.double(sapply((1:q),  function(x) paste0(x, (1:q)))), q,q))),(1:k))
  names(out$theta) <- c(if (con$msmu==TRUE) paste0("mu_",mu_n_tmp[,1],",",mu_n_tmp[,2]) else  paste0("mu_",(1:q)),
                        if (con$msvar==TRUE) paste0("sig_",sig_n_tmp[,1],",",sig_n_tmp[,2]) else paste0("sig_",covar_vech(t(matrix(as.double(sapply((1:q),  function(x) paste0(x, (1:q)))), q,q)))),
                        paste0("phi_",paste0(phi_n_tmp[,2],",",phi_n_tmp[,3],phi_n_tmp[,1])),
                        paste0("p_",c(sapply((1:k),  function(x) paste0(x, (1:k)) ))))
  if (con$getSE==TRUE){
    out <- thetaSE(out)
  }
  if (is.null(con$init_theta)){
    out$trace <- output_all
  }
  return(out)
}
