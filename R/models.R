

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
    b0    <- inv(t(X)%*%X)%*%t(X)%*%y
    phi   <- b0[(2:npar),1]
    inter <- b0[1,1]
  }else{
    X     <- x
    b0    <- inv(t(X)%*%X)%*%t(X)%*%y
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
                  theta_phi_ind = theta_phi_ind, stationary = roots, n = n, p = p, q = 1, k = 1)
  # Define class
  class(out) <- "ARmdl"
  # get log-likelihood
  out$logLike <- logLikelihood(out)
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
  npar <- p
  # ----- estimate model
  if (con$const==TRUE){
    X     <- cbind(1,x)
    npar  <- npar + 1
    b0    <- inv(t(X)%*%X)%*%t(X)%*%y
    inter <- t(b0[1,])
    phi   <- t(b0[-1,])
  }else{
    X     <- x
    b0    <- inv(t(X)%*%X)%*%t(X)%*%y
    inter <- matrix(0,q,1)
    phi   <- t(b0)
  }
  # ----- Obtain variables of interest
  Fmat        <- companionMat(phi,inter,p,q)
  stationary  <- all(abs(eigen(Fmat)[[1]])<1)
  mu_tmp      <- inv(diag(q*p)-Fmat)%*%as.matrix(c(inter,rep(0,q*(p-1))))
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
                  stationary = stationary, n = n, p = p, q = q, k = 1, Fmat = Fmat)
  # Define class
  class(out) <- "VARmdl"
  # get log-likelihood
  out$logLike <- logLikelihood(out)
  # get standard errors if 'con$getSE' is 'TRUE'
  if (con$getSE==TRUE){
    out <- thetaSE(out)
  }
  return(out)
}


# HMM


#' @title Markov-switching autoregressive model by Expectation Minimization algorithm
#' 
#' @description This function estimates a Markov-switching autoregressive model using the Expectation Minimization (EM) algorithm of Dempster, Laird and Rubin (1977) as explained in Hamilton (1990).
#' 
#' @param Y (Tx1) vector with observational data. Required argument.
#' @param p integer for the number of lags to use in estimation. Must be greater than or equal to 0. Default is 0.
#' @param k integer for the number of regimes to use in estimation. Must be greater than or equal to 2. Default is 2.
#' @param control List with optimization options including: 
#' \itemize{
#'  \item{msmu - }{indicator for switch in mean (TRUE) or constant mean (FALSE). Default is TRUE.}
#'  \item{msvar - }{bool indicator for switch in variance (TRUE) or constant variance (FALSE). Default is TRUE.}
#'  \item{maxit - }{integer determining the maximum number of EM iterations.}
#'  \item{thtol - }{double determining the convergence criterion for the absolute difference in parameter estimates (theta) between iterations. Default is 1e-6.}
#'  \item{getSE - }{bool if 'TRUE' standard errors are computed and returned. If 'FALSE' standard errors are not computed. Default is 'FALSE'.}
#'  \item{max_init - }{integer determining the maximum number of initial values to attempt in case solution is NaN. Default is 500.}
#'  \item{use_diff_init - }{integer determining how many different initial values (that do not return NaN) to try. Default is 1.}
#'  \item{init_value - }{vector of initial values. This is optional. Default is NUll, in which case 'initValsMS()' is used to generate initial values.}
#' }
#' 
#' @return List with model characteristics
#' 
#' @references Dempster, A. P., N. M. Laird, and D. B. Rubin. 1977. “Maximum Likelihood from Incomplete Data via the EM Algorithm.” Journal of the Royal Statistical Society. Series B 39 (1): 1–38
#' @references Hamilton, James D. 1990. “Analysis of time series subject to changes in regime.” Journal of econometrics, 45 (1-2): 39–70
#' 
#' @export
#' 
#' @examples 
#' # Define DGP of MS AR process
#' mdl_ms2 <- list(n     = 500, 
#'                 mu    = c(5,10),
#'                 sigma = c(1,2),
#'                 phi   = c(0.5,0.2),
#'                 k     = 2,
#'                 P     = rbind(c(0.1,0.9),
#'                               c(0.9,0.1)))
#'                               
#' # Simulate process using simuMSAR() function
#' y_ms_simu <- simuMSAR(mdl_ms2)
#' 
#' # Set options for model estimation
#' control <- list(msmu  = TRUE, 
#'                 msvar = TRUE, 
#'                 use_diff_init = 10)
#' 
#' # Estimate model
#' y_ms_mdl <- MSARmdl(y_ms_simu$y, ar = y_ms_simu$ar, k = y_ms_simu$k, control)
#' 
MSARmdl <- function(Y, p, k, control = list()){
  # ----- Set control values
  con <- list(const = TRUE,
              getSE = TRUE,
              msmu = TRUE, 
              msvar = TRUE,
              method = "EM",
              maxit = 10000, 
              thtol = 1.e-6, 
              max_init = 500, 
              use_diff_init = 1, 
              init_theta = NULL)
  # Perform some checks for controls
  nmsC <- names(con)
  con[(namc <- names(control))] <- control
  if(length(noNms <- namc[!namc %in% nmsC])){
    warning("unknown names in control: ", paste(noNms,collapse=", ")) 
  }
  # ---------- Optimization options
  optim_options <- list(maxit = con$maxit, thtol = con$thtol)
  # pre-define list and matrix length for results
  EM_output_all <- list(con$use_diff_init)
  max_loglik <- matrix(0, con$use_diff_init, 1)
  if (is.null(con$init_theta)==TRUE){
    # ---------- Estimate linear model to use for initial values
    init_control <- list(con$const, con$getSE)
    init_mdl <- ARmdl(Y, p = p, init_control)
    init_mdl$msmu = con$msmu
    init_mdl$msvar = con$msvar
    # ----- Estimate using 'use_diff_init' different initial values
    for (xi in 1:con$use_diff_init){
      init_used <- 0
      converge_check <- FALSE
      while ((converge_check==FALSE) & (init_used<con$max_init)){
        # ----- Initial values
        theta_0 <- initValsMS(init_mdl, k)
        if (con$method=="EM"){
          # ----- Estimate using EM algorithm and initial values provided
          EM_output_tmp <- MS_EMest(theta_0, init_mdl, k, optim_options)
          EM_output_tmp$theta_0 = theta_0
          # ----- Convergence check
          logLike_tmp = EM_output_tmp$logLike
          theta_tmp = EM_output_tmp$theta
          converge_check = ((is.finite(logLike_tmp)) & (all(is.finite(theta_tmp))))
        }else if (con$method=="NR"){
          # ----- Estimate using roptim and initial values provided 
          
          # ----- Convergence check
          
        }
        init_used = init_used + 1
      }
      max_loglik[xi] = logLike_tmp
      EM_output_tmp$init_used = init_used
      EM_output_all[[xi]] = EM_output_tmp
    }
    if (con[["use_diff_init"]]==1){
      EM_output = EM_output_tmp
    }else{
      xl = which.max(max_loglik)
      if (length(xl)==0){
        warning("Model(s) did not converge. Use higher 'use_diff_init' or 'max_init'.")
        EM_output = EM_output_all[[1]] 
      }else{
        EM_output = EM_output_all[[xl]] 
      }
    }
  }else{
    
    if (con$method=="EM"){
      # ----- Estimate using EM algorithm and initial values provided
      EM_output <- MS_EMest(con$init_value, mdl_out, k, optim_options)
      EM_output$theta_0 <- con$init_value
      EM_output$init_used <- 1  
    }else if (con$method=="NR"){
      # ----- Estimate using roptim and initial values provided 
      
    }
  }
  # ---------- organize output
  theta_mu_ind <- c(rep(1, 1 + (k-1)*con$msmu), rep(0, 1 + (k-1)*con$msvar + p + k*k))
  theta_sig_ind <- c(rep(0, 1 + (k-1)*con$msmu), rep(1, 1 + (k-1)*con$msvar), rep(0, p + k*k))
  theta_phi_ind <- c(rep(0, 2 + (k-1)*con$msmu + (k-1)*con$msvar), rep(1, p), rep(0, k*k))
  theta_P_ind <- c(rep(0, 2 + (k-1)*con$msmu + (k-1)*con$msvar + p), rep(1, k*k))
  # ----- output
  out <- list(y = init_mdl$y, X = init_mdl$X, x = init_mdl$x, resid = EM_output$residuals, mu = EM_output$mu, coef = NA, intercept = NA, phi = EM_output$phi,
              stdev = sqrt(EM_output$sigma), sigma = EM_output$sigma, theta = EM_output$theta, theta_mu_ind = theta_mu_ind, theta_sig_ind = theta_sig_ind, 
              theta_phi_ind = theta_phi_ind, stationary = NA, n = init_mdl$n, p = p, q = 1, k = k, logLike = EM_output$logLike, P = EM_output$P, pinf = EM_output$pinf, 
              St = EM_output$St, eta = EM_output$eta, thl = EM_output$thl, deltath = EM_output$deltath,  iterations = EM_output$iterations, theta_0 = EM_output$theta_0,
              init_used = EM_output$init_used, msmu = con$msmu, msvar = con$msvar, control = con)
  # Define class
  class(out) <- "MSARmdl"
  if (con$getSE==TRUE){
    out <- thetaSE(out)
  }
  if (is.null(con$init_theta)){
    out$trace <- EM_output_all
  }
  return(out)
}


#' @title Markov-switching vector autoregressive model by Expectation Minimization algorithm
#' 
#' @description This function estimates a Markov-switching vector autoregressive model using the Expectation Minimization (EM) algorithm of Dempster, Laird and Rubin (1977) as explained in Krolzig (1997).
#' 
#' @param Y (Tx1) vector with observational data. Required argument.
#' @param ar integer for the number of lags to use in estimation. Must be greater than or equal to 0. Default is 0.
#' @param k integer for the number of regimes to use in estimation. Must be greater than or equal to 2. Default is 2.
#' @param control List with optimization options including:
#' \itemize{
#'  \item{msmu - }{indicator for switch in mean (TRUE) or constant mean (FALSE). Default is TRUE.}
#'  \item{msvar - }{bool indicator for switch in variance (TRUE) or constant variance (FALSE). Default is TRUE.}
#'  \item{maxit - }{integer determining the maximum number of EM iterations.}
#'  \item{thtol - }{double determining the convergence criterion for the absolute difference in parameter estimates (theta) between iterations. Default is 1e-6.}
#'  \item{getSE - }{bool if 'TRUE' standard errors are computed and returned. If 'FALSE' standard errors are not computed. Default is 'FALSE'.}
#'  \item{max_init - }{integer determining the maximum number of initial values to attempt in case solution is NaN. Default is 500.}
#'  \item{use_diff_init - }{integer determining how many different initial values (that do not return NaN) to try. Default is 1.}
#'  \item{init_value - }{vector of initial values. This is optional. Default is NUll, in which case 'initValsMS()' is used to generate initial values.}
#' }
#' 
#' @return List with model characteristics
#' 
#' @references Dempster, A. P., N. M. Laird, and D. B. Rubin. 1977. “Maximum Likelihood from Incomplete Data via the EM Algorithm.” Journal of the Royal Statistical Society. Series B 39 (1): 1–38
#' @references Krolzig, Hans-Martin. 1997. “The markov-switching vector autoregressive model.” In Markov-Switching Vector Autoregressions, 6–28. Springer
#' @export
#' 
#' @examples 
#' # Define DGP of MS VAR process
#' mdl_msvar2 <- list(n     = 500, 
#'                    ar    = 1,
#'                    mu    = rbind(c(5,-2),
#'                                  c(10,2)),
#'                    sigma = list(rbind(c(5.0, 1.5),
#'                                       c(1.5, 1.0)),
#'                                 rbind(c(7.0, 3.0),
#'                                       c(3.0, 2.0))),
#'                    phi   = rbind(c(0.50, 0.30),
#'                                  c(0.20, 0.70)),
#'                    k     = 2,
#'                    P     = rbind(c(0.1,0.9),
#'                                  c(0.9,0.1)))
#'                                  
#' # Simulate process using simuMSVAR() function
#' y_msvar_simu <- simuMSVAR(mdl_msvar2)
#' 
#' # Set options for model estimation
#' control <- list(msmu = TRUE, 
#'                 msvar = TRUE, 
#'                 use_diff_init = 10)
#'                 
#' # Estimate model
#' y_msvar_mdl <- MSVARmdl_EM(y_msvar_simu$y, ar = y_msvar_simu$ar, k = y_msvar_simu$k, control)
#' 
MSVARmdl_EM <- function(Y, ar, k, control = list()){
  # ----- Set control values
  con <- list(msmu = TRUE, 
              msvar = TRUE, 
              maxit = 10000, 
              thtol = 1.e-6, 
              getSE = FALSE, 
              max_init = 500, 
              use_diff_init = 1, 
              init_value = NULL)
  # Perform some checks for controls
  nmsC <- names(con)
  con[(namc <- names(control))] <- control
  if(length(noNms <- namc[!namc %in% nmsC])){
    warning("unknown names in control: ", paste(noNms,collapse=", ")) 
  }
  # ---------- Optimization options
  optim_options <- list(maxit = con[["maxit"]], thtol = con[["thtol"]])
  # ---------- Estimate linear model to use for initial values
  mdl_out = VARmdl(Y, ar = ar, intercept = TRUE, getSE = con[["getSE"]])
  msmu <- con[["msmu"]]
  msvar <- con[["msvar"]]
  mdl_out[["msmu"]] = msmu
  mdl_out[["msvar"]] = msvar
  # ---------- Estimate model 
  EM_output_all <- list(con[["use_diff_init"]])
  max_loglik <- matrix(0, con[["use_diff_init"]], 1)
  if (is.null(con[["init_value"]])==FALSE){
    # ----- Estimate using initial values provided
    EM_output <- MSVAR_EMest(con[["init_value"]], mdl_out, k, optim_options)
    EM_output[["theta_0"]] <- con[["init_value"]]
    EM_output[["init_used"]] <- 1
  }else{
    for (xi in 1:con[["use_diff_init"]]){
      init_used <- 0
      converge_check <- FALSE
      while ((converge_check==FALSE) & (init_used<con[["max_init"]])){
        # ----- Initial values
        theta_0 <- initValsMSVAR(mdl_out, k)
        # ----- Estimate using EM algorithm 
        EM_output_tmp <- MSVAR_EMest(theta_0, mdl_out, k, optim_options)
        EM_output_tmp[["theta_0"]] = theta_0
        # ----- Convergence check
        logLike_tmp = EM_output_tmp[["logLike"]]
        theta_tmp = EM_output_tmp[["theta"]]
        converge_check = ((is.finite(logLike_tmp)) & (all(is.finite(theta_tmp))))
        init_used = init_used + 1
      }
      max_loglik[xi] = logLike_tmp
      EM_output_tmp[["init_used"]] = init_used
      EM_output_all[[xi]] = EM_output_tmp
    }
    if (con[["use_diff_init"]]==1){
      EM_output = EM_output_tmp
    }else{
      xl = which.max(max_loglik)
      if (length(xl)==0){
        warning("Model(s) did not converge. Use higher 'use_diff_init' or 'max_init'.")
        EM_output = EM_output_all[[1]] 
      }else{
        EM_output = EM_output_all[[xl]] 
      }
    }
  }
  # ---------- organize output
  q <- ncol(Y)
  Nsig <- (q*(q+1))/2
  phi_len <- length(c(mdl_out[["phi"]]))
  theta_mu_ind <- c(rep(1, q + q*(k-1)*msmu), rep(0, Nsig + Nsig*(k-1)*msvar + phi_len + k*k))
  theta_sig_ind <- c(rep(0, q + q*(k-1)*msmu), rep(1, Nsig + Nsig*(k-1)*msvar), rep(0, phi_len + k*k))
  theta_var_ind <- c(rep(0, q + q*(k-1)*msmu),rep(t(covar_vech(diag(q))), 1+(k-1)*msvar), rep(0, phi_len + k*k))
  theta_phi_ind <- c(rep(0, q + q*(k-1)*msmu + Nsig + Nsig*(k-1)*msvar), rep(1, phi_len), rep(0, k*k))
  theta_P_ind <- c(rep(0, q + q*(k-1)*msmu + Nsig + Nsig*(k-1)*msvar + phi_len), rep(1, k*k))
  MSVARmdl_output <- EM_output
  MSVARmdl_output[["theta_mu_ind"]] = theta_mu_ind
  MSVARmdl_output[["theta_sig_ind"]] = theta_sig_ind
  MSVARmdl_output[["theta_var_ind"]] = theta_var_ind
  MSVARmdl_output[["theta_phi_ind"]] = theta_phi_ind
  MSVARmdl_output[["theta_P_ind"]] = theta_P_ind
  MSVARmdl_output[["y"]] = mdl_out[["y"]]
  MSVARmdl_output[["ar"]] = ar
  MSVARmdl_output[["phi"]] = t(MSVARmdl_output[["phi"]])
  MSVARmdl_output[["n"]] = mdl_out[["n"]]
  MSVARmdl_output[["q"]] = mdl_out[["q"]]
  MSVARmdl_output[["k"]] = k
  MSVARmdl_output[["x"]] = mdl_out[["x"]]
  MSVARmdl_output[["X"]] = mdl_out[["X"]]
  MSVARmdl_output[["msmu"]] = con[["msmu"]]
  MSVARmdl_output[["msvar"]] = con[["msvar"]]
  MSVARmdl_output[["control"]] <- con
  if (con[["getSE"]]==TRUE){
    Hess <- getHess(MSVARmdl_output, k)
    info_mat <- solve(-Hess)
    nearPD_used <- FALSE
    if ((all(is.na(Hess)==FALSE)) & (any(diag(info_mat)<0))){
      info_mat <- nearPD(info_mat)
      nearPD_used <- TRUE
    }
    MSVARmdl_output[["Hess"]] <- Hess
    MSVARmdl_output[["theta_stderr"]] <- sqrt(diag(info_mat))
    MSVARmdl_output[["info_mat"]] <- info_mat
    MSVARmdl_output[["nearPD_used"]] <- nearPD_used
  }
  if (is.null(con[["init_value"]])){
    MSVARmdl_output[["trace"]] <- EM_output_all
  }
  return(MSVARmdl_output)
}


# MSVARmdl (Optimization)
