#' @title Approximate CDF distribution 
#'
#' @description This function obtains the parameters in eq. 16 of the CDF distribution needed for combining moment-based test statistics. 
#'
#' @param Tsize Sample size.
#' @param simdist_N Number of simulations to approximate CDF distribution.
#' 
#' @return A (\code{2 x 4}) matrix with parameters of CDF distributions. The first row contains \eqn{\gamma_0} for each moment and the second row contains \eqn{\gamma_1} for each moment.
#' 
#' @keywords internal
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
    mdl <- stats::nls(Fx[,i]~exp(alpha+beta*x[,i])/(1+exp(alpha+beta*x[,i])),
                      start=list(alpha=a_start,beta=b_start))
    params <- stats::coef(mdl)
    a<-cbind(a,params[1])
    b<-cbind(b,params[2])
  }
  return(rbind(matrix(a,nrow=1,ncol=4),matrix(b,nrow=1,ncol=4)))
}


#' @title  Monte Carlo moment-based test for Markov switching model
#' 
#' @description This function performs the Local Monte Carlo moment-based test for
#' Markov switching autoregressive models proposed in Dufour & Luger (2017).
#' 
#' @param Y Series to be tested
#' @param p Number of autoregressive lags.
#' @param control List with test procedure options including: 
#' \itemize{
#'   \item N: Integer determining the number of Monte Carlo simulations. Default is set to \code{99} as in paper.
#'   \item simdist_N: Integer determining the number of simulations for CDF distribution approximation. Default is set to \code{10000}.
#'   \item getSE: Boolean indicator. If \code{TRUE}, standard errors for restricted model are estimated. If \code{FALSE} no standard errors are estimated. Default is \code{TRUE}.
#' }
#' 
#' @return List of class \code{DLMCTest} (\code{S3} object) with attributes including: 
#' \itemize{
#'   \item mdl_h0: List with restricted model attributes. This will be of class \code{ARmdl} if \code{p>0} or \code{Nmdl} otherwise (\code{S3} objects). See \code{\link{ARmdl}} or \code{\link{Nmdl}}.    
#'   \item theta: Value of nuisance parameters. Specifically, these are the consistent estimates of nuisance parameters as discussed in Dufour & Luger (2017) LMC procedure.
#'   \item S0: A (\code{1 x 4})) matrix containing the four moment-based test statistics defined in (\code{11}) - (\code{14}) in Dufour & Luger (2017).
#'   \item F0_min: Test statistic value for min version of Local Monte Carlo moment-based test.
#'   \item F0_prod: Test statistic value for prod version of Local Monte Carlo moment-based test.
#'   \item FN_min: A (\code{N x 1}) vector with simulated test statistics for min version of Local Monte Carlo moment-based test under null hypothesis.
#'   \item FN_prod: A (\code{N x 1}) vector with simulated test statistics for prod version of Local Monte Carlo moment-based test under null hypothesis.
#'   \item pval_min: P-value for min version of Local Monte Carlo moment-based test.
#'   \item pval_prod: P-value for prod version of Local Monte Carlo moment-based test.
#'   \item FN_min_cv: Vector with 90\%, 95\%, and 99\% Monte Carlo critical values for min version of Local Monte Carlo moment-based test.
#'   \item FN_prod_cv: Vector with 90\%, 95\%, and 99\% Monte Carlo critical values for prod version of Local Monte Carlo moment-based test.
#'   \item control: List with test procedure options used.
#' }
#'
#' @references Dufour, J. M., & Luger, R. 2017. "Identification-robust moment-based tests for 
#' Markov switching in autoregressive models." \emph{Econometric Reviews}, 36(6-9), 713-727.
#' 
#' @example /inst/examples/DLMCTest_examples.R
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
  theta     <- as.matrix(mdl_h0$phi)
  rownames(theta) <- paste0("phi_", seq(1:p))
  # ----- Simulate distribution 
  params    <- approxDistDL(Tsize-p, con$simdist_N)
  sim_ms    <- sim_DLmoments(Tsize, con$N)
  Fmin_sim  <- as.matrix(sort(combine_stat(sim_ms, params, "min")))
  Fprd_sim  <- as.matrix(sort(combine_stat(sim_ms, params, "prod")))
  # ----- Compute test stat
  eps       <- mdl_h0$resid
  S0        <- t(calc_DLmoments(eps))
  colnames(S0) <- c("M(\U03B5)","V(\U03B5)","S(\U03B5)","K(\U03B5)")
  # ----- get critical values
  Fmin_sim_cv   <- Fmin_sim[round(c(0.90,0.95,0.99)*nrow(Fmin_sim)),]
  Fprd_sim_cv   <- Fprd_sim[round(c(0.90,0.95,0.99)*nrow(Fprd_sim)),]
  names(Fmin_sim_cv)  <- paste0(c("0.90","0.95","0.99"), "%")
  names(Fprd_sim_cv)  <- paste0(c("0.90","0.95","0.99"), "%")
  # combine moment stats
  F0_min  <- combine_stat(S0, params, "min")
  F0_prod <- combine_stat(S0, params, "prod")
  colnames(F0_min) <- "F(\U03B5)"
  colnames(F0_prod) <- "F(\U03B5)"
  # ----- Get p-values
  pval_min  <- MCpval(F0_min, Fmin_sim, "geq")
  pval_prod <- MCpval(F0_prod, Fprd_sim, "geq")
  # ----- Organize output
  DLMCTest_output <- list(mdl_h0 = mdl_h0, S0 = S0, F0_min = F0_min, F0_prod = F0_prod, FN_min = Fmin_sim, FN_prod = Fprd_sim, 
                          pval_min = pval_min, pval_prod = pval_prod, theta = theta, FN_min_cv = Fmin_sim_cv, FN_prod_cv = Fprd_sim_cv, 
                          control = con)
  class(DLMCTest_output) <- "DLMCTest"
  return(DLMCTest_output)
}


#' @title MMC nuisance parameter bounds for Moment-based test
#' 
#' @description This function is used to determine the lower and upper bounds for the MMC LRT parameter search.
#' 
#' @param mdl_h0 List with restricted model properties.
#' @param con List with control options provided to \code{DLMMCTest} procedure.
#' 
#' @return List with \code{theta_low}, vector of parameter lower bounds, and \code{theta_upp}, vector of parameter upper bounds.
#' 
#' @keywords internal
#' 
#' @export
DLMMC_bounds <- function(mdl_h0, con){
  theta_0 <- mdl_h0$phi
  # ----- Define lower & upper bounds for MMC search
  theta_low <- theta_0 - con$eps
  theta_upp <- theta_0 + con$eps
  # if CI_union==TRUE use union of eps & 2*standard error to define bounds
  phi_ind <- mdl_h0$theta_phi_ind
  phiSE <- mdl_h0$theta_se[phi_ind==1]
  if ((con$CI_union==TRUE) & all(is.finite(phiSE))){
    theta_low <- apply(cbind(as.matrix(theta_0 - 2*phiSE), as.matrix(theta_low)), 1, FUN = min)
    theta_upp <- apply(cbind(as.matrix(theta_0 + 2*phiSE), as.matrix(theta_upp)), 1, FUN = max)
  }
  if (is.null(con$phi_low)==FALSE){
    theta_low <- apply(cbind(as.matrix(theta_low),as.matrix(con$phi_low)), 1, function(x) max(x))
  }  
  if (is.null(con$phi_upp)==FALSE){
    theta_upp <- apply(cbind(as.matrix(theta_upp),as.matrix(con$phi_upp)), 1, function(x) min(x))
  }  
  # ----- output
  mmc_bounds <- list(theta_low = theta_low, theta_upp = theta_upp)
  return(mmc_bounds)
}






#' @title Maximized Monte Carlo moment-based test for Markov switching model
#' 
#' @description This function performs the maximized Monte Carlo moment-based test for
#' Markov switching autoregressive models proposed in Dufour & Luger (2017).
#' 
#' @param Y Series to be tested
#' @param p Number of autoregressive lags.
#' @param control List with test procedure options including: 
#' \itemize{
#'   \item N: Integer determining the number of Monte Carlo simulations. Default is set to \code{99} as in paper.
#'   \item simdist_N: Integer determining the number of simulations for CDF distribution approximation. Default is set to \code{10000}.
#'   \item getSE: Boolean indicator. If \code{TRUE}, standard errors for restricted model are estimated. If \code{FALSE} no standard errors are estimated. Default is \code{TRUE}.
#'   \item eps: Fixed positive constant that does not depend on \code{T} used to determine lower and upper bounds on consistent set considered for nuisance parameter space.
#'   \item CI_union: Boolean indicator determining if union between \code{eps} and confidence interval is used to determine lower and upper bound on consistent set considered for nuisance parameter space. If \code{TRUE} union is used and if \code{FALSE} only \code{eps} is used. Note that if standard errors obtained are not finite then only \code{eps} is used. Default is \code{FALSE}.       
#'   \item lambda: Numeric value for penalty on stationary constraint not being met. Default is \code{100}.
#'   \item stationary_ind: Boolean indicator determining if only stationary solutions should be considered if \code{TRUE} or any solution can be considered if \code{FALSE}. Default is \code{TRUE}.
#'   \item phi_low: Vector with lower bound for autoregressive parameters when optimizing. Default is \code{NULL}.
#'   \item phi_upp: Vector with upper bound for autoregressive parameters when optimizing. Default is \code{NULL}.
#'   \item optim_type: String determining type of numerical optimization algorithm to use. Available options are: "\code{\link{pso}}", ""\code{\link{GenSA}}", "\code{\link{GA}}". Default is "\code{\link{GenSA}}".
#'   \item silence: Boolean indicator determining if optimization updates should be silenced if \code{TRUE} or not if \code{FALSE}. Default is \code{FALSE}.
#'   \item threshold_stop: Numeric value determining the maximum possible p-value attainable. Default is \code{1}.
#'   \item type_control: List containing other optimization options specific to the numerical optimization algorithm used. This includes maximum number of iterations which is \code{200} b y default. For other options see documentation of numerical algorithm chosen.
#' }
#' 
#' @return List of class \code{DLMCTest} (\code{S3} object) with attributes including: 
#' \itemize{
#'   \item mdl_h0: List with restricted model attributes. This will be of class \code{ARmdl} if \code{p>0} or \code{Nmdl} otherwise (\code{S3} objects). See \code{\link{ARmdl}} or \code{\link{Nmdl}}.    
#'   \item theta_max_min: Value of nuisance parameters when min version of p-value is maximized as discussed in Dufour & Luger (2017) MMC procedure.
#'   \item theta_max_prod: Value of nuisance parameters when prod version of p-value is maximized as discussed in Dufour & Luger (2017) MMC procedure.
#'   \item theta_low: Lower bound on nuisance parameter values used when searching for maximum p-value.
#'   \item theta_upp: Upper bound on nuisance parameter values used when searching for maximum p-value.
#'   \item S0_min: A (\code{1 x 4})) matrix containing the four moment-based test statistics defined in (\code{11}) - (\code{14}) in Dufour & Luger (2017) when \code{theta_min} is used.
#'   \item S0_prod: A (\code{1 x 4})) matrix containing the four moment-based test statistics defined in (\code{11}) - (\code{14}) in Dufour & Luger (2017) when \code{theta_prod} is used.
#'   \item F0_min: Test statistic value for min version of Maximized Monte Carlo moment-based test.
#'   \item F0_prod: Test statistic value for prod version of Maximized Monte Carlo moment-based test.
#'   \item FN_min: A (\code{N x 1}) vector with simulated test statistics for min version of Maximized Monte Carlo moment-based test under null hypothesis.
#'   \item FN_prod: A (\code{N x 1}) vector with simulated test statistics for prod version of Maximized Monte Carlo moment-based test under null hypothesis.
#'   \item pval_min: Maximum p-value for min version of Maximized Monte Carlo moment-based test.
#'   \item pval_prod: Maximum p-value for prod version of Local Monte Carlo moment-based test.
#'   \item FN_min_cv: Vector with 90\%, 95\%, and 99\% Monte Carlo critical values for min version of Local Monte Carlo moment-based test.
#'   \item FN_prod_cv: Vector with 90\%, 95\%, and 99\% Monte Carlo critical values for prod version of Local Monte Carlo moment-based test.
#'   \item control: List with test procedure options used.
#'   \item optim_min_output: List with optimization output for min version of Maximized Monte Carlo moment-based test.
#'   \item optim_prod_output: List with optimization output for prod version of Maximized Monte Carlo moment-based test.
#' }
#'
#' @references Dufour, J. M., & Luger, R. 2017. "Identification-robust moment-based tests for 
#' Markov switching in autoregressive models." \emph{Econometric Reviews}, 36(6-9), 713-727.
#' 
#' @example /inst/examples/DLMMCTest_examples.R
#' @export
DLMMCTest <- function(Y, p, control = list()){
  # ----- Set control values
  con <- list(N = 99,
              simdist_N = 10000,
              getSE = TRUE,
              eps = 0.1,
              CI_union = FALSE,
              lambda = 100,
              stationary_ind = TRUE,
              phi_low = NULL,
              phi_upp = NULL,
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
  }else{
    stop("No Nuisance parameters is model is not Autoregressive. Number of lags must be greater than 0.")
  }
  # ----- Get parameters from approximated distribution 
  params <- approxDistDL(Tsize-p, con$simdist_N)
  # ----- Simulated process eta_i is fixed (see pg. 721 of Dufour & Luger 2017) and so simulated statistics are fixed 
  sim_ms    <- sim_DLmoments(Tsize, con$N)
  Fmin_sim  <-  as.matrix(sort(combine_stat(sim_ms, params, "min")))
  Fprd_sim  <-  as.matrix(sort(combine_stat(sim_ms, params, "prod")))
  # ----- Define starting values, lower & upper bounds for search
  theta_0 <- mdl_h0$phi
  mmc_bounds <- DLMMC_bounds(mdl_h0, con)
  theta_low <- mmc_bounds$theta_low
  theta_upp <- mmc_bounds$theta_upp
  # ----- Search for Max p-value within bounds
  if(con$optim_type=="pso"){
    # Set PSO specific controls
    con$type_control$trace.stats <- TRUE
    con$type_control$trace <- as.numeric(con$silence==FALSE)
    con$type_control$abstol <- -con$threshold_stop
    # begin optimization
    mmc_min_out <- pso::psoptim(par = theta_0, fn = DLMMCpval_fun_min, lower = theta_low, upper = theta_upp, 
                                gr = NULL, control = con$type_control,
                                y = y, x = x, params = params, sim_stats = Fmin_sim, 
                                pval_type = "min", stationary_ind = con$stationary_ind, lambda = con$lambda)
    mmc_prd_out <- pso::psoptim(par = theta_0, fn = DLMMCpval_fun_min, lower = theta_low, upper = theta_upp, 
                                gr = NULL, control = con$type_control,
                                y = y, x = x, params = params, sim_stats = Fprd_sim, 
                                pval_type = "prod", stationary_ind = con$stationary_ind, lambda = con$lambda)
    theta_min   <- as.matrix(mmc_min_out$par)
    theta_prod  <- as.matrix(mmc_prd_out$par)
    pval_min    <- -mmc_min_out$value
    pval_prod   <- -mmc_prd_out$value
  }else if(con$optim_type=="GenSA"){
    # Set GenSA specific controls
    con$type_control$verbose <- con$silence==FALSE
    con$type_control$threshold.stop <- -con$threshold_stop
    mmc_min_out <- GenSA::GenSA(par = theta_0, fn = DLMMCpval_fun_min, lower = theta_low, upper = theta_upp, 
                                control = con$type_control,
                                y = y, x = x, params = params, sim_stats = Fmin_sim, 
                                pval_type = "min", stationary_ind = con$stationary_ind, lambda = con$lambda)
    mmc_prd_out <- GenSA::GenSA(par = theta_0, fn = DLMMCpval_fun_min, lower = theta_low, upper = theta_upp, 
                                control = con$type_control,
                                y = y, x = x, params = params, sim_stats = Fprd_sim, 
                                pval_type = "prod", stationary_ind = con$stationary_ind, lambda = con$lambda)
    theta_min   <- as.matrix(mmc_min_out$par)
    theta_prod  <- as.matrix(mmc_prd_out$par)
    pval_min    <- -mmc_min_out$value
    pval_prod   <- -mmc_prd_out$value
  }else if(con$optim_type=="GA"){
    mmc_min_out <- GA::ga(type = "real-valued", fitness = DLMMCpval_fun, 
                          y = y, x = x, params = params, sim_stats = Fmin_sim, 
                          pval_type = "min", stationary_ind = con$stationary_ind, lambda = con$lambda,
                          lower = theta_low, upper = theta_upp, 
                          maxiter = con$type_control$maxit, maxFitness = con$threshold_stop, 
                          monitor = (con$silence==FALSE), suggestions = t(theta_0))
    mmc_prd_out <- GA::ga(type = "real-valued", fitness = DLMMCpval_fun, 
                          y = y, x = x, params = params, sim_stats = Fprd_sim, 
                          pval_type = "prod", stationary_ind = con$stationary_ind, lambda = con$lambda,
                          lower = theta_low, upper = theta_upp, 
                          maxiter = con$type_control$maxit, maxFitness = con$threshold_stop, 
                          monitor = (con$silence==FALSE), suggestions = t(theta_0))
    theta_min <- as.matrix(mmc_min_out@solution[1,])  # keeps on the first set that gives optim output
    theta_prod <- as.matrix(mmc_prd_out@solution[1,]) # keeps on the first set that gives optim output
    pval_min <- mmc_min_out@fitnessValue
    pval_prod <- mmc_prd_out@fitnessValue
    
  }
  # ----- get test output using optimization output params
  z_min <- y - x%*%theta_min
  z_prd <- y - x%*%theta_prod
  rownames(theta_min) <- paste0("phi_", seq(1:p))
  rownames(theta_prod) <- paste0("phi_", seq(1:p))
  # Compute test stats
  eps_min <- z_min - mean(z_min)
  eps_prd <- z_prd - mean(z_prd)
  S0_min  <- t(calc_DLmoments(eps_min))
  S0_prd  <- t(calc_DLmoments(eps_prd))
  colnames(S0_min) <- c("M(\U03B5)","V(\U03B5)","S(\U03B5)","K(\U03B5)")
  colnames(S0_prd) <- c("M(\U03B5)","V(\U03B5)","S(\U03B5)","K(\U03B5)")
  # get critical values
  Fmin_sim_cv   <- Fmin_sim[round(c(0.90,0.95,0.99)*nrow(Fmin_sim)),]
  Fprd_sim_cv   <- Fprd_sim[round(c(0.90,0.95,0.99)*nrow(Fprd_sim)),]
  names(Fmin_sim_cv)  <- paste0(c("0.90","0.95","0.99"), "%")
  names(Fprd_sim_cv)  <- paste0(c("0.90","0.95","0.99"), "%")
  # combine moment stats
  F0_min  <- combine_stat(S0_min, params, "min")
  F0_prod <- combine_stat(S0_prd, params, "prod")
  colnames(F0_min) <- "F(\U03B5)"
  colnames(F0_prod) <- "F(\U03B5)"
  # ----- organize test output
  DLMMCTest_output <- list(mdl_h0 = mdl_h0, S0_min = S0_min, S0_prod = S0_prd, F0_min = F0_min, F0_prod = F0_prod, 
                           FN_min = Fmin_sim, FN_prod = Fprd_sim, pval_min = pval_min, pval_prod = pval_prod, 
                           theta_max_min = theta_min, theta_max_prod = theta_prod, theta_low = theta_low, theta_upp = theta_upp, 
                           FN_min_cv = Fmin_sim_cv,  FN_prod_cv = Fprd_sim_cv, control = con, 
                           optim_min_output = mmc_min_out, optim_prod_output = mmc_prd_out)
  class(DLMMCTest_output) <- "DLMMCTest"
  return(DLMMCTest_output)
}

