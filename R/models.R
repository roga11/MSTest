#' @title Normal dis model
#' 
#' @description This function estimates a normally distributed model
#' 
#' @param Y matrix of observations with dimension (n x q) 
#' @param control List with model options including:
#' \itemize{
#' \item{const - }{boolean determining whether to estimate model with constant, if 'TRUE', or not, if 'FALSE'.}
#' \item{getSE - }{boolean determining whether to compute standard errors of parameters, if 'TRUE', or not, if 'FALSE'.}
#' }
#' 
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
  resid       <- Y - matrix(1, n, 1)%*%t(as.matrix(mu))
  sigma       <- (t(resid)%*%resid)/(n-1)
  stdev       <- sqrt(diag(sigma))
  theta       <- c(mu,covar_vech(sigma))
  # theta indicators
  theta_mu_ind  <- c(rep(1, q),rep(0,length(theta)-q))
  theta_sig_ind <- c(rep(0, q),rep(1,q*(q+1)/2))
  # ----- Output
  out     <- list(y = Y, resid = resid, mu = mu, stdev = stdev, sigma = sigma, theta = theta, 
                  theta_mu_ind = theta_mu_ind, theta_sig_ind = theta_sig_ind, n = n, q = q, k = 1, control = con)
  # Define class
  class(out) <- "Nmdl"
  # get log-likelihood
  out$logLike <- logLikelihood(out)
  # get information criterion
  out$AIC <- aic(out$logLike, length(out$theta))
  out$BIC <- bic(out$logLike, out$n, length(out$theta))
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
#' @description This function estimates an autoregresive model
#' 
#' @param Y vector of observations with dimension (n x 1)
#' @param p integer determining the number of autoregressive lags
#' @param control List with model options including:
#' \itemize{
#' \item{const - }{boolean determining whether to estimate model with constant, if 'TRUE', or not, if 'FALSE'.}
#' \item{getSE - }{boolean determining whether to compute standard errors of parameters, if 'TRUE', or not, if 'FALSE'.}
#' }
#' 
#' @return List with model attributes which include:
#' \itemize{
#'   \item{y - }{vector of observations of dimension (n x 1).}
#'   \item{X - }{matrix of lagged observations (with or without vector of 1s depending on const='TRUE' or const='FALSE').}
#'   \item{x - }{matrix of lagged observations without vector of 1s.}
#'   \item{resid - }{vector of residuals.}
#'   \item{mu - }{mean of the process}
#'   \item{coef - }{coefficient estimates. This is the same as phi if const='FALSE'.}
#'   \item{intercept - }{coefficient estimate of intercept.}
#'   \item{phi - }{autoregressive coefficient estimates. This is the same as coef if const='FALSE'.}
#'   \item{stdev - }{standard deviations.}
#'   \item{sigma - }{variance.}
#'   \item{theta - }{vector containing: mu, sigma, and phi.}
#'   \item{theta_mu_ind - }{vector indicating location of mean.}
#'   \item{theta_sig_ind - }{vector indicating location of variance.}
#'   \item{theta_phi_ind - }{vector indicating location of autoregressive coefficients.}
#'   \item{stationary - }{bool indicating if process is stationary 'TRUE' or non-stationary 'FALSE'.}
#'   \item{n - }{number of observations after transofrmation due to lags (i.e., T-p observations).}
#'   \item{p - }{number of autoregressive parameters.}
#'   \item{q - }{number of serires. This is always 1 in ARmdl.}
#'   \item{k - }{number of regimes. This is always 1 in ARmdl.}
#'   \item{logLike - }{log-likelihood.}
#'   \item{Hess - }{Hessian matrix. Approximated using numDeriv package and only returned if getSE='TRUE'.}
#'   \item{info_mat - }{Information matrix. Computed as the inverse of -Hess which is approximated using numDeriv package. If matrix is not PD then nearest PD matrix is obtained using nearPD. Only returned if getSE='TRUE'.}
#'   \item{nearPD_used - }{Bool determining whether nearPD was used on infoMat 'TRUE' or not 'FALSE'. Only returned if getSE='TRUE'.}
#'   \item{theta_se - }{standard errors of parameters in theta. Only returned if getSE='TRUE'.}
#' }
#' 
#' @seealso \code{\link{MSARmdl}}
#' @example /examples/ARmdl_examples.R
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
    b0    <- solve(t(X)%*%X)%*%t(X)%*%y
    phi   <- b0[(2:npar),1]
    inter <- b0[1,1]
  }else{
    X     <- x
    b0    <- solve(t(X)%*%X)%*%t(X)%*%y
    phi   <- b0
    inter <- 0
  }
  # ----- Obtain variables of interest
  phisum  <- sum(phi)
  mu      <- inter/(1-phisum)
  resid   <- y - X%*%b0
  stdev   <- sqrt((t(resid)%*%resid)/(n-1))
  sigma   <- stdev^2
  theta   <- as.matrix(c(mu,sigma,phi))
  roots   <- all(Mod(as.complex(polyroot(c(1,-phi))))>1)
  #theta   <- as.matrix(c(mu,stdev,phi))
  # theta indicators
  theta_mu_ind  <- c(1,rep(0,length(theta)-1))
  theta_sig_ind <- c(0,1,rep(0,p))
  theta_phi_ind <- c(0,0,rep(1,p))
  # ----- Output
  out     <- list(y = y, X = X, x = x, resid = resid, mu = mu, coef = b0, intercept = inter, phi = phi,
                  stdev = stdev, sigma = sigma, theta = theta, theta_mu_ind = theta_mu_ind, theta_sig_ind = theta_sig_ind, 
                  theta_phi_ind = theta_phi_ind, stationary = roots, n = n, p = p, q = 1, k = 1, control = con)
  # Define class
  class(out) <- "ARmdl"
  # get log-likelihood
  out$logLike <- logLikelihood(out)
  # get information criterion
  out$AIC <- aic(out$logLike, length(out$theta))
  out$BIC <- bic(out$logLike, out$n, length(out$theta))
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
#' @description This function estimates a vector autoregresive model
#' 
#' @param Y matrix of observations with dimension (n x q) 
#' @param p integer determining the number of autoregressive lags
#' @param control List with model options including:
#' \itemize{
#' \item{const - }{boolean determining whether to estimate model with constant, if 'TRUE', or not, if 'FALSE'.}
#' \item{getSE - }{boolean determining whether to compute standard errors of parameters, if 'TRUE', or not, if 'FALSE'.}
#' }
#' 
#' @return List with model attributes which include:
#' \itemize{
#'   \item{y - }{vector of observations of dimension (n x q).}
#'   \item{X - }{matrix of lagged observations (with or without vector of 1s depending on const='TRUE' or const='FALSE').}
#'   \item{x - }{matrix of lagged observations without vector of 1s.}
#'   \item{resid - }{vector of residuals.}
#'   \item{mu - }{vector of dimension (q x 1) containing means of each process.}
#'   \item{coef - }{coefficient estimates. This is the same as phi if const='FALSE'.}
#'   \item{intercept - }{vector of dimension (q x 1) containing coefficient estimate of intercepts.}
#'   \item{phi - }{matrix of dimension (q x q*p) autoregressive coefficient estimates. This is the same as coef if const='FALSE'.}
#'   \item{stdev - }{standard deviations of each process (i.e., square root of diagonal of 'sigma'.)}
#'   \item{sigma - }{covariance matrix.}
#'   \item{theta - }{vector containing: mu, vech(sigma), and phi.}
#'   \item{theta_mu_ind - }{vector indicating location of mean.}
#'   \item{theta_sig_ind - }{vector indicating location of variance.}
#'   \item{theta_phi_ind - }{vector indicating location of autoregressive coefficients.}
#'   \item{stationary - }{bool indicating if process is stationary 'TRUE' or non-stationary 'FALSE'.}
#'   \item{n - }{number of observations after transofrmation due to lags (i.e., T-p observations).}
#'   \item{p - }{number of autoregressive lags.}
#'   \item{q - }{number of serires.}
#'   \item{k - }{number of regimes. This is always 1 in VARmdl.}
#'   \item{Fmat - }{matrix of dimension (qp x qp) of companion form.}
#'   \item{logLike - }{log-likelihood.}
#'   \item{Hess - }{Hessian matrix. Approximated using numDeriv package and only returned if getSE='TRUE'.}
#'   \item{info_mat - }{Information matrix. Computed as the inverse of -Hess which is approximated using numDeriv package. If matrix is not PD then nearest PD matrix is obtained using nearPD. Only returned if getSE='TRUE'.}
#'   \item{nearPD_used - }{Bool determining whether nearPD was used on infoMat 'TRUE' or not 'FALSE'. Only returned if getSE='TRUE'.}
#'   \item{theta_se - }{standard errors of parameters in theta. Only returned if getSE='TRUE'.}
#' }
#' 
#' @seealso \code{\link{MSVARmdl}}
#' @example /examples/VARmdl_examples.R
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
  # ----- Get process dimensions
  n <- nrow(Y)
  q <- ncol(Y)
  # ----- transform data
  lagged_vals <- ts_lagged(Y, p)
  y <- lagged_vals$y
  x <- lagged_vals$X
  # ----- estimate model
  if (con$const==TRUE){
    X     <- cbind(1,x)
    b0    <- solve(t(X)%*%X)%*%t(X)%*%y
    inter <- t(b0[1,])
    phi   <- t(b0[-1,])
  }else{
    X     <- x
    b0    <- solve(t(X)%*%X)%*%t(X)%*%y
    inter <- matrix(0,q,1)
    phi   <- t(b0)
  }
  # ----- Obtain variables of interest
  Fmat        <- companionMat(phi,p,q)
  stationary  <- all(abs(eigen(Fmat)[[1]])<1)
  mu_tmp      <- solve(diag(q*p)-Fmat)%*%as.matrix(c(inter,rep(0,q*(p-1))))
  mu          <- mu_tmp[(1:(q))]
  resid       <- y - X%*%b0
  sigma       <- (t(resid)%*%resid)/(n-1)
  stdev       <- sqrt(diag(sigma))
  theta       <- c(mu,covar_vech(sigma),c(t(phi)))
  # theta indicators
  theta_mu_ind  <- c(rep(1,q),rep(0,length(theta)-q))
  theta_sig_ind <- c(rep(0,q),rep(1,q*(q+1)/2),rep(0,q*q*p))
  theta_phi_ind <- c(rep(0,length(theta)-q*q*p),rep(1,q*q*p))
  # ----- Output
  out     <- list(y = y, X = X, x = x, resid = resid, mu = mu, coef = b0, intercept = inter, phi = phi,
                  stdev = stdev, sigma = sigma, theta = theta, theta_mu_ind = theta_mu_ind, theta_sig_ind = theta_sig_ind, theta_phi_ind = theta_phi_ind, 
                  stationary = stationary, n = n, p = p, q = q, k = 1, Fmat = Fmat, control = con)
  # Define class
  class(out) <- "VARmdl"
  # get log-likelihood
  out$logLike <- logLikelihood(out)
  # get information criterion
  out$AIC <- aic(out$logLike, length(out$theta))
  out$BIC <- bic(out$logLike, out$n, length(out$theta))
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
#' @description This function estimates a Hidden Markov model
#' 
#' @param Y (Txq) vector with observational data. Required argument.
#' @param p integer for the number of lags to use in estimation. Must be greater than or equal to 0. Default is 0.
#' @param k integer for the number of regimes to use in estimation. Must be greater than or equal to 2. Default is 2.
#' @param control List with optimization options including:
#' 
#' @return List with model characteristics
#' 
#' @references Hamilton, James D. 1990. “Analysis of time series subject to changes in regime.” Journal of econometrics, 45 (1-2): 39–70
#' @references Dempster, A. P., N. M. Laird, and D. B. Rubin. 1977. “Maximum Likelihood from Incomplete Data via the EM Algorithm.” Journal of the Royal Statistical Society. Series B 39 (1): 1–38
#' @references Krolzig, Hans-Martin. 1997. “The markov-switching vector autoregressive model.” In Markov-Switching Vector Autoregressions, 6–28. Springer
#' 
#' @seealso \code{\link{Nmdl}}
#' @example /examples/HMmdl_examples.R
#' @export
HMmdl <- function(Y, k, control = list()){
  # ----- Set control values
  con <- list(getSE = TRUE,
              msmu = TRUE, 
              msvar = TRUE,
              init_value = NULL,
              method = "EM",
              maxit = 10000,
              thtol = 1.e-6, 
              maxit_converge = 500, 
              use_diff_init = 1,
              mle_variance_constraint = 1e-6)
  # ----- Perform some checks for controls
  nmsC <- names(con)
  con[(namc <- names(control))] <- control
  if(length(noNms <- namc[!namc %in% nmsC])){
    warning("unknown names in control: ", paste(noNms,collapse=", ")) 
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
      output <- HMmdl_em(con$init_value, init_mdl, k, optim_options)
      output$theta_0 <- con$init_value
      output$init_used <- 1  
    }else if (con$method=="MLE"){
      init_mdl$mle_variance_constraint <- con$mle_variance_constraint
      # ----- Estimate using roptim and initial values provided 
      output_tmp <- HMmdl_mle(con$init_value, init_mdl, k, optim_options)  
      output$theta_0 <- con$init_value
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
  out <- list(y = init_mdl$y, resid = output$resid, mu = output$mu, sigma = output$sigma, theta = output$theta, 
              theta_mu_ind = theta_mu_ind, theta_sig_ind = theta_sig_ind, n = init_mdl$n, q = q, k = k, 
              logLike = output$logLike, P = output$P, pinf = output$pinf, St = output$St, eta = output$eta, 
              thl = output$thl, deltath = output$deltath,  iterations = output$iterations, theta_0 = output$theta_0,
              init_used = output$init_used, msmu = con$msmu, msvar = con$msvar, control = con)
  out$AIC <- aic(out$logLike, length(out$theta))
  out$BIC <- bic(out$logLike, out$n, length(out$theta))
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
  # Define class
  class(out) <- "HMmdl"
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
#' @param Y (T x 1) vector with observational data. Required argument.
#' @param p integer for the number of lags to use in estimation. Must be greater than or equal to 0. Default is 0.
#' @param k integer for the number of regimes to use in estimation. Must be greater than or equal to 2. Default is 2.
#' @param control List with optimization options including: 
#' \itemize{
#'  \item{getSE - }{bool if 'TRUE' standard errors are computed and returned. If 'FALSE' standard errors are not computed. Default is 'FALSE'.}
#'  \item{msmu - }{indicator for switch in mean (TRUE) or constant mean (FALSE). Default is 'TRUE.'}
#'  \item{msvar - }{bool indicator for switch in variance (TRUE) or constant variance (FALSE). Default is 'TRUE.'}
#'  \item{init_value - }{vector of initial values. This is optional. Default is 'NULL', in which case \code{initVals_MSARmdl} is used to generate initial values.}
#'  \item{method - }{string determining which method to use. Options are 'EM' for EM algorithm or 'MLE' for Maximum Likelihood Estimation.}
#'  \item{maxit - }{integer determining the maximum number of EM iterations.}
#'  \item{thtol - }{double determining the convergence criterion for the absolute difference in parameter estimates (theta) between iterations. Default is 1e-6.}
#'  \item{maxit_converge - }{integer determining the maximum number of initial values attempted until solution is finite. For example, if parameters in 'theta' or 'logLike' is NaN another set of initital values (up to 'maxit_converge') is attempted until finite values are returned. This does not occur frequently for most types of data but may be useful in some cases. Once finite values are obtained, this counts as one iteration in 'maxit'. Default is 500.}
#'  \item{use_diff_init - }{integer determining how many different initial values (that do not return NaN) to try. Default is 1.}
#'  \item{mle_stationary_constraint - }{bool indicator determining if only stationary solutions should be considered for autoregressive coefficients (if 'TRUE') or if non-stationary solutions are allowed (if 'FALSE'). This is only used when method='MLE'. Default is 'TRUE'.}
#'  \item{mle_variance_constraint - }{double used to determine the lower bound on variance. Specifically, this value is multiplied by the variance of the linear (one regime) model to determine the lower bound. For example, if '0.01' then lower bound is 1 \% of variance from linear model. This is only used when method='MLE' and should be between 0 and 1. Default is '0.01'.}
#' }
#' 
#' @return List with model characteristics
#' 
#' @references Dempster, A. P., N. M. Laird, and D. B. Rubin. 1977. “Maximum Likelihood from Incomplete Data via the EM Algorithm.” Journal of the Royal Statistical Society. Series B 39 (1): 1–38
#' @references Hamilton, James D. 1990. “Analysis of time series subject to changes in regime.” Journal of econometrics, 45 (1-2): 39–70
#' 
#' @seealso \code{\link{ARmdl}}
#' @example /examples/MSARmdl_examples.R
#' @export
MSARmdl <- function(Y, p, k, control = list()){
  # ----- Set control values
  con <- list(getSE = TRUE,
              msmu = TRUE, 
              msvar = TRUE,
              init_theta = NULL,
              method = "EM",
              maxit = 10000, 
              thtol = 1.e-6, 
              maxit_converge = 500, 
              use_diff_init = 1, 
              mle_stationary_constraint = TRUE,
              mle_variance_constraint = 0.01)
  # ----- Perform some checks for controls
  nmsC <- names(con)
  con[(namc <- names(control))] <- control
  if(length(noNms <- namc[!namc %in% nmsC])){
    warning("unknown names in control: ", paste(noNms,collapse=", ")) 
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
      output <- MSARmdl_em(con$init_value, init_mdl, k, optim_options)
      output$theta_0 <- con$init_value
      output$init_used <- 1  
    }else if (con$method=="MLE"){
      init_mdl$mle_stationary_constraint <- con$mle_stationary_constraint
      init_mdl$mle_variance_constraint <- con$mle_variance_constraint
      # ----- Estimate using roptim and initial values provided 
      output_tmp <- MSARmdl_mle(con$init_value, init_mdl, k, optim_options)  
      output$theta_0 <- con$init_value
      output$init_used <- 1  
    }
  }
  # ----- Obtain variables of interest
  theta_mu_ind <- c(rep(1, 1 + (k-1)*con$msmu), rep(0, 1 + (k-1)*con$msvar + p + k*k))
  theta_sig_ind <- c(rep(0, 1 + (k-1)*con$msmu), rep(1, 1 + (k-1)*con$msvar), rep(0, p + k*k))
  theta_phi_ind <- c(rep(0, 2 + (k-1)*con$msmu + (k-1)*con$msvar), rep(1, p), rep(0, k*k))
  theta_P_ind <- c(rep(0, 2 + (k-1)*con$msmu + (k-1)*con$msvar + p), rep(1, k*k))
  # ----- Output
  out <- list(y = init_mdl$y, X = init_mdl$X, x = init_mdl$x, resid = output$resid, mu = output$mu, 
              coef = NULL, intercept = NULL, phi = output$phi, stdev = sqrt(output$sigma), sigma = output$sigma, 
              theta = output$theta, theta_mu_ind = theta_mu_ind, theta_sig_ind = theta_sig_ind, 
              theta_phi_ind = theta_phi_ind, theta_P_ind = theta_P_ind, stationary = NULL, 
              n = init_mdl$n, p = p, q = 1, k = k, logLike = output$logLike, P = output$P, pinf = output$pinf, 
              St = output$St, eta = output$eta, thl = output$thl, deltath = output$deltath, 
              iterations = output$iterations, theta_0 = output$theta_0, init_used = output$init_used, 
              msmu = con$msmu, msvar = con$msvar, control = con)
  out$AIC <- aic(out$logLike, length(out$theta))
  out$BIC <- bic(out$logLike, out$n, length(out$theta))
  # names 
  names(out$theta) <- c(if (con$msmu==TRUE) paste0("mu_", (1:k)) else "mu",
                        if (con$msvar==TRUE) paste0("sig_", (1:k)) else "sig",
                        paste0("phi_",(1:p)),
                        paste0("p_",c(sapply((1:k),  function(x) paste0(x, (1:k)) ))))
  # Define class
  class(out) <- "MSARmdl"
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
#' @param Y (Txq) vector with observational data. Required argument.
#' @param p integer for the number of lags to use in estimation. Must be greater than or equal to 0. Default is 0.
#' @param k integer for the number of regimes to use in estimation. Must be greater than or equal to 2. Default is 2.
#' @param control List with optimization options including:
#' \itemize{
#'  \item{getSE - }{bool if 'TRUE' standard errors are computed and returned. If 'FALSE' standard errors are not computed. Default is 'FALSE'.}
#'  \item{msmu - }{indicator for switch in mean (TRUE) or constant mean (FALSE). Default is TRUE.}
#'  \item{msvar - }{bool indicator for switch in variance (TRUE) or constant variance (FALSE). Default is TRUE.}
#'  \item{init_value - }{vector of initial values. This is optional. Default is NULL, in which case \code{initVals_MSARmdl} is used to generate initial values.}
#'  \item{method - }{string determining which method to use. Options are 'EM' for EM algorithm or 'MLE' for Maximum Likelihood Estimation.}
#'  \item{maxit - }{integer determining the maximum number of EM iterations.}
#'  \item{thtol - }{double determining the convergence criterion for the absolute difference in parameter estimates (theta) between iterations. Default is 1e-6.}
#'  \item{maxit_converge - }{integer determining the maximum number of initial values attempted until solution is finite. For example, if parameters in 'theta' or 'logLike' is NaN another set of initital values (up to 'maxit_converge') is attempted until finite values are returned. This does not occur frequently for most types of data but may be useful in some cases. Once finite values are obtained, this counts as one iteration in 'maxit'. Default is 500.}
#'  \item{use_diff_init - }{integer determining how many different initial values (that do not return NaN) to try. Default is 1.}
#'  \item{mle_stationary_constraint - }{bool indicator determining if only stationary solutions should be considered for autoregressive coefficients (if 'TRUE') or if non-stationary solutions are allowed (if 'FALSE'). This is only used when method='MLE'. Default is 'TRUE'.}
#'  \item{mle_variance_constraint - }{double used to determine the lower bound on the smallest eigenvalue for the covariance matrix of each regime. Default is '1e-6'.}
#' }
#' 
#' @return List with model characteristics
#' 
#' @references Dempster, A. P., N. M. Laird, and D. B. Rubin. 1977. “Maximum Likelihood from Incomplete Data via the EM Algorithm.” Journal of the Royal Statistical Society. Series B 39 (1): 1–38
#' @references Krolzig, Hans-Martin. 1997. “The markov-switching vector autoregressive model.” In Markov-Switching Vector Autoregressions, 6–28. Springer
#' 
#' @seealso \code{\link{VARmdl}}
#' @example /examples/MSVARmdl_examples.R
#' @export
MSVARmdl <- function(Y, p, k, control = list()){
  # ----- Set control values
  con <- list(getSE = TRUE,
              msmu = TRUE, 
              msvar = TRUE,
              init_value = NULL,
              method = "EM",
              maxit = 10000,
              thtol = 1.e-6, 
              maxit_converge = 500, 
              use_diff_init = 1,
              mle_stationary_constraint = TRUE,
              mle_variance_constraint = 1e-6)
  # ----- Perform some checks for controls
  nmsC <- names(con)
  con[(namc <- names(control))] <- control
  if(length(noNms <- namc[!namc %in% nmsC])){
    warning("unknown names in control: ", paste(noNms,collapse=", ")) 
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
      output <- MSVARmdl_em(con$init_value, init_mdl, k, optim_options)
      output$theta_0 <- con$init_value
      output$init_used <- 1  
    }else if (con$method=="MLE"){
      init_mdl$mle_stationary_constraint <- con$mle_stationary_constraint
      init_mdl$mle_variance_constraint <- con$mle_variance_constraint
      # ----- Estimate using roptim and initial values provided 
      output_tmp <- MSVARmdl_mle(con$init_value, init_mdl, k, optim_options)  
      output$theta_0 <- con$init_value
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
  out <- list(y = init_mdl$y, X = init_mdl$X, x = init_mdl$x, resid = output$resid, mu = output$mu, coef = NULL, intercept = NULL, phi = output$phi,
              sigma = output$sigma, theta = output$theta, theta_mu_ind = theta_mu_ind, theta_sig_ind = theta_sig_ind, 
              theta_phi_ind = theta_phi_ind, stationary = NULL, n = init_mdl$n, p = p, q = q, k = k, logLike = output$logLike, P = output$P, pinf = output$pinf, 
              St = output$St, eta = output$eta, thl = output$thl, deltath = output$deltath,  iterations = output$iterations, theta_0 = output$theta_0,
              init_used = output$init_used, msmu = con$msmu, msvar = con$msvar, control = con)
  out$AIC <- aic(out$logLike, length(out$theta))
  out$BIC <- bic(out$logLike, out$n, length(out$theta))
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
  # Define class
  class(out) <- "MSVARmdl"
  if (con$getSE==TRUE){
    out <- thetaSE(out)
  }
  if (is.null(con$init_theta)){
    out$trace <- output_all
  }
  return(out)
}
