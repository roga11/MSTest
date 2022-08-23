#' @title Approximate CDF distribution 
#'
#' @description This function obtains the parameters in eq. 16 of the CDF distribution needed for combining moment-based test statistics. 
#'
#' @param Tsize Sample size.
#' @param simdist_N Number of simulations to approximate CDF distribution.
#' 
#' @return A (\code{2 x 4}) matrix with parameters of CDF distributions. The first row contains \eqn{\gamma_0} for each moment and the second row contains \eqn{\gamma_1} for each moment.
#' 
#' @references Dufour, J. M., & Luger, R. 2017. "Identification-robust moment-based 
#' tests for Markov switching in autoregressive models." \emph{Econometric Reviews}, 36(6-9), 713-727.
#' 
#' @export
approxDistDL <- function(Tsize, simdist_N){
  S_N2  <- sim_DLmoments(Tsize, simdist_N)
  x     <- apply(S_N2, 2, sort)
  Fx    <- approx_dist_loop(x)
  a_start <- 0.01
  b_start <- 0.01
  # initiate matrix
  a<-matrix(nrow=1,ncol=0)
  b<-matrix(nrow=1,ncol=0)
  # estimate params of ecdf for each moment statistic
  for (i in 1:4){
    mdl<-nls(Fx[,i]~exp(alpha+beta*x[,i])/(1+exp(alpha+beta*x[,i])),
             start=list(alpha=a_start,beta=b_start))
    params<-coef(mdl)
    a<-cbind(a,params[1])
    b<-cbind(b,params[2])
  }
  return(rbind(matrix(a,nrow=1,ncol=4),matrix(b,nrow=1,ncol=4)))
}


# ==============================================================================
#' @title  Monte Carlo moment-based test for Markov switching model
#' 
#' @description This function performs the Local Monte-Carlo moment-based test for
#' Markov-switching autoregressive models proposed in Dufour & Luger (2017).
#' 
#' @param Y Series to be tested
#' @param p Number of autoregressive lags.
#' @param control List with test procedure options including: 
#' \itemize{
#'   \item{\code{N}: }{Integer determining the number of Monte Carlo simulations. Default is set to \code{99} as in paper.}
#'   \item{\code{simdist_N}: }{Integer determining the number of simulations for CDF distribution approximation. Default is set to \code{10000}.}
#'   \item{\code{getSE}: }{Boolean indicator. If \code{TRUE}, standard errors for restricted model are estimated. If \code{FALSE} no standard errors are estimated. Default is \code{TRUE}.}
#' }
#' 
#' @return List of class \code{DLMCTest} (\code{S3} object) with model attributes including: 
#' \itemize{
#'   \item{\code{mdl_h0}: }{List with restricted model attributes. This will be of class \code{ARmdl} if \code{p>0} or \code{Nmdl} otherwise (\code{S3} objects). See \code{\link{ARmdl}} or \code{\link{Nmdl}}.}    
#'   \item{\code{Fmin}: }{test statistic value for min version of Local Monte Carlo moment-based test.}
#'   \item{\code{Fprod}: }{test statistic value for prod version of Local Monte Carlo moment-based test.}
#'   \item{\code{Fmin_N}: }{A (\code{N x 1}) vector with simulated test statistics for min version of Local Monte Carlo moment-based test under null hypothesis.}
#'   \item{\code{Fprod_N}: }{A (\code{N x 1}) vector with simulated test statistics for prod version of Local Monte Carlo moment-based test under null hypothesis.}
#'   \item{\code{pval_min}: }{P-value for min version of Local Monte Carlo moment-based test.}
#'   \item{\code{pval_prod}: }{P-value for prod version of Local Monte Carlo moment-based test.}
#'   \item{\code{Fmin_cv}: }{Vector with 90\%, 95\%, and 99\% Monte Carlo critical values for min version of Local Monte Carlo moment-based test.}
#'   \item{\code{Fprod_cv}: }{Vector with 90\%, 95\%, and 99\% Monte Carlo critical values for prod version of Local Monte Carlo moment-based test.}
#'   \item{\code{control}: }{List with test procedure options used.}
#' }
#'
#' @references Dufour, J. M., & Luger, R. 2017. "Identification-robust moment-based tests for 
#' Markov switching in autoregressive models." \emph{Econometric Reviews}, 36(6-9), 713-727.
#' 
#' @example /examples/DLMCTest_examples.R
#' @export
DLMCTest <- function(Y, p, control = list()){
  # ----- Set control values
  con <- list(N = 99,
              simdist_N = 10000,
              getSE = TRUE)
  # Perform some checks for controls
  nmsC <- names(con)
  con[(namc <- names(control))] <- control
  if(length(noNms <- namc[!namc %in% nmsC])){
    warning("unknown names in control: ", paste(noNms,collapse=", ")) 
  }
  # ----- Transform process (from eq. (20) to z_t(phi))
  Tsize <- length(Y)
  null_control <- list(const = TRUE, getSE = con$getSE)
  if(p>0){
    mdl_h0 <- ARmdl(Y, p, null_control)
  }else{
    mdl_h0 <- Nmdl(Y, null_control)
  }
  # --------- Get parameters from approximated distribution ----------
  params <- approxDistDL(Tsize-p, con$simdist_N)
  # -------------------------- Get P-Values --------------------------
  eps   <- mdl_h0$resid
  Fmin  <- calc_DLmcstat(eps, con$N, params, "min")
  Fprod <- calc_DLmcstat(eps, con$N, params, "prod")
  # ----- Obtain p-value
  Fmin0     <- Fmin[con$N+1]
  Fprod0    <- Fprod[con$N+1]
  FminSim   <- as.matrix(sort(Fmin[1:con$N]))
  FprodSim  <- as.matrix(sort(Fprod[1:con$N]))
  pval_min  <- MCpval(Fmin0, FminSim, "geq")
  pval_prod <- MCpval(Fprod0, FprodSim, "geq")
  Fmin_cv   <- FminSim[round(c(0.90,0.95,0.99)*nrow(FminSim)),]
  Fprod_cv  <- FprodSim[round(c(0.90,0.95,0.99)*nrow(FprodSim)),]
  names(Fmin_cv)  <- paste0(c("0.90","0.95","0.99"), "%")
  names(Fprod_cv) <- paste0(c("0.90","0.95","0.99"), "%")
  # ----- Organize output
  DLMCTest_output <- list(mdl_h0 = mdl_h0, Fmin = Fmin0, Fprod = Fprod0, Fmin_N = FminSim, Fprod_N = FprodSim, 
                          pval_min = pval_min, pval_prod = pval_prod, Fmin_cv = Fmin_cv, Fprod_cv = Fprod_cv, 
                          control = con)
  class(DLMCTest_output) <- "DLMCTest"
  return(DLMCTest_output)
}





# ==============================================================================
#' @title Maximized Monte Carlo moment-based test for Markov switching model
#' 
#' @description Performs the Maximized Monte Carlo (MMC) moment-based test described in Dufour & Luger (2017). 
#' 
#' @param Y Series to be tested
#' @param p Number of autoregressive lags.
#' @param pval_type Method to be used to combine p-values (see eq. (17) and (18) of Dufour & Luger (2017))
#' @param control List of control parameters. See “Details”.
#'
#' @return 
#' 
#' @references Dufour, J. M., & Luger, R. 2017. "Identification-robust moment-based tests for 
#' Markov switching in autoregressive models." \emph{Econometric Reviews}, 36(6-9), 713-727.
#' 
#' @export
DLMMCTest <- function(Y, p, control = list()){
  # ----- Set control values
  con <- list(N = 99,
              pval_type = "min",
              simdist_N = 10000,
              getSE = TRUE,
              eps = 0.1,
              CI_union = TRUE,
              lambda = 100,
              stationary_ind = TRUE,
              optim_type = "GenSA",
              silence = FALSE,
              threshold_stop = 1,
              type_control = list(maxit = 200))
  # Perform some checks for controls
  nmsC <- names(con)
  con[(namc <- names(control))] <- control
  if(length(noNms <- namc[!namc %in% nmsC])){
    warning("unknown names in control: ", paste(noNms,collapse=", ")) 
  }
  if ((con$CI_union==TRUE) & (con$getSE==FALSE)){
    con$getSE <- TRUE
    warning("getSE was changed to be 'TRUE' because CI_union is 'TRUE'.")
  }
  # ----- Set initial value as consistent estimate
  Tsize <- length(Y)
  null_control <- list(const = TRUE, getSE = con$getSE)
  if(p>0){
    mdl_h0 <- ARmdl(Y, p, null_control)
    y <- mdl_h0$y
    x <- mdl_h0$x
    theta_0 <- mdl_h0$phi
  }else{
    stop("No Nuisance parameters is model is not Autoregressive. Number of lags must be greater than 0.")
  }
  # ----- Get parameters from approximated distribution 
  params <- approxDistDL(Tsize-p, con$simdist_N)
  # ----- Simulted process eta_i is fixed (see pg. 721 of Dufour & Luger 2017) and so simulated statistics are fixed 
  sim_ms  <- sim_DLmoments(Tsize, con$N)
  Fsim    <-  as.matrix(sort(combine_stat(sim_ms, params, con$pval_type)))
  # ----- Define lower & upper bounds for MMC search
  theta_low <- theta_0 - con$eps
  theta_upp <- theta_0 + con$eps
  # if CI_union==TRUE use union of eps & 2*standard error to define bounds
  phiSE <- mdl_h0$theta_se[3:(3+p-1)]
  if ((con$CI_union==TRUE) & all(is.finite(phiSE))){
    theta_low <- apply(cbind(as.matrix(theta_0 - 2*phiSE), as.matrix(theta_low)), 1, FUN = min)
    theta_upp <- apply(cbind(as.matrix(theta_0 + 2*phiSE), as.matrix(theta_upp)), 1, FUN = max)
  }
  # ----- Search for Max p-value within bounds
  if(con$optim_type=="pso"){
    # Set PSO specific controls
    con$type_control$trace.stats <- TRUE
    con$type_control$trace <- as.numeric(con$silence==FALSE)
    con$type_control$abstol <- -con$threshold_stop
    # begin optimization
    mmc_out <- pso::psoptim(par = theta_0, fn = DLMMCpval_fun_min, lower = theta_low, upper = theta_upp, 
                            gr = NULL, control = con$type_control,
                            y = y, x = x, N = con$N, simdist_N = con$simdist_N, params = params, sim_stats = Fsim, 
                            pval_type = con$pval_type, stationary_ind = con$stationary_ind, lambda = con$lambda)
    theta <- mmc_out$par
    pval <- -mmc_out$value
  }else if(con$optim_type=="GenSA"){
    # Set GenSA specific controls
    con$type_control$verbose <- con$silence==FALSE
    con$type_control$threshold.stop <- -con$threshold_stop
    mmc_out <- GenSA::GenSA(par = theta_0, fn = DLMMCpval_fun_min, lower = theta_low, upper = theta_upp, 
                            control = con$type_control,
                            y = y, x = x, N = con$N, simdist_N = con$simdist_N, params = params, sim_stats = Fsim, 
                            pval_type = con$pval_type, stationary_ind = con$stationary_ind, lambda = con$lambda)
    theta <- mmc_out$par
    pval <- -mmc_out$value
  }else if(con$optim_type=="GA"){
    mmc_out <- GA::ga(type = "real-valued", fitness = DLMMCpval_fun, 
                      y = y, x = x, N = con$N, simdist_N = con$simdist_N, params = params, sim_stats = Fsim, 
                      pval_type = con$pval_type, stationary_ind = con$stationary_ind, lambda = con$lambda,
                      lower = theta_low, upper = theta_upp, 
                      maxiter = con$type_control$maxit, maxFitness = con$threshold_stop, 
                      monitor = (con$silence==FALSE), suggestions = t(theta_0))
    theta <- c(mmc_out@solution)
    pval <- mmc_out@fitnessValue
  }
  # ----- get test-stat using optimization output params
  z = y - x%*%theta
  # Compute test stats
  eps = z - mean(z)
  S0 = t(calc_DLmoments(eps))
  # ----- get critical values
  Fsim_cv   <- Fsim[round(c(0.90,0.95,0.99)*nrow(Fsim)),]
  names(Fsim_cv)  <- paste0(c("0.90","0.95","0.99"), "%")
  # combine moment stats
  F0 = combine_stat(S0, params, con$pval_type)
  # ----- organize remaining output
  DLMMCTest_output <- list(mdl_h0 = mdl_h0, moment_test_stats = S0, combined_test_stat = F0, combined_sim_test_stats = Fsim, 
                          pval = pval, theta_max = theta, theta_low = theta_low, theta_upp = theta_upp, sim_test_stats_cv = Fsim_cv, 
                          control = con, optim_output = mmc_out)
  class(DLMMCTest_output) <- "DLMMCTest"
  return(DLMMCTest_output)
}

