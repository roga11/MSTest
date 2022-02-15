# ==============================================================================
#' @title  Monte-Carlo Moment-based test for Markov-switching model
#' 
#' @description This function performs the Local Monte-Carlo Moment-Based test for
#' Markov-switching autoregressive models proposed in Dufour & Luger (2017).
#' 
#' @param Y Series to be tested
#' @param ar Order of autoregressive components AR(p)
#' @param control List of control parameters. See “Details”.
#' 
#' @details The control argument is a list that can supply any of the following components:
#' \itemize{
#'   \item N: length of generate series fosample distribution. Default is N=99.
#'   \item simdist_N: length of simulation used to approximate distribution of p-value distribution. 
#'  Default is simdist_N=10000
#'  \item getSE: Boolean value which determines if standard errors should be returned with AR model 
#'  (only when AR model is used). If TRUE, standard errors are obtained. Default is getSE=FALSE.
#' }
#' 
#' @return List with model and test results
#'
#' @references Dufour, J. M., & Luger, R. (2017). Identification-robust moment-based tests for 
#' Markov switching in autoregressive models. Econometric Reviews, 36(6-9), 713-727.
#' 
#' @export
DLMCTest <- function(Y, ar, control = list()){
  # ----- Set control values
  con <- list(N = 99,
              simdist_N = 10000,
              getSE = FALSE)
  # Perform some checks for controls
  nmsC <- names(con)
  con[(namc <- names(control))] <- control
  if(length(noNms <- namc[!namc %in% nmsC])){
    warning("unknown names in control: ", paste(noNms,collapse=", ")) 
  }
  # ----- Transform process (from eq. (20) to z_t(phi))
  Tsize <- length(Y)
  if(ar>0){
    mdl_out <- ARmdl(Y, ar, TRUE, con[["getSE"]])
    y <- mdl_out[["y"]]
    x <- mdl_out[["x"]]
    phi <- mdl_out[["phi"]]
    z <- y - x%*%phi
  }else{
    z <- Y
  }
  # --------- Get parameters from approximated distribution ----------
  params <- approxDistDL(Tsize-ar, con[["simdist_N"]])
  # -------------------------- Get P-Values --------------------------
  N     <- con[["N"]]
  eps   <- z - mean(z)
  Fmin  <- calc_DLmcstat(eps, N, params, "min")
  Fprod <- calc_DLmcstat(eps, N, params, "prod")
  # ----- Obtain p-value
  Fmin0     <- Fmin[N+1]
  Fprod0    <- Fprod[N+1]
  FminSim   <- Fmin[1:N]
  FprodSim  <- Fprod[1:N]
  pval_min  <- MCpval(Fmin0, FminSim, "geq")
  pval_prod <- MCpval(Fprod0, FprodSim, "geq")
  LMC_N     <- cbind(Fmin,Fprod)
  # ----- Organize output
  DLMCMTest_output <- list()
  DLMCMTest_output[["eps"]] <- eps
  if (ar>0){
    DLMCMTest_output[["ARmdl"]] <- mdl_out
  }
  DLMCMTest_output[["dist_params"]] <- params
  DLMCMTest_output[["LMC_N"]] <- LMC_N
  DLMCMTest_output[["Fmin_0"]] <- Fmin0
  DLMCMTest_output[["Fprod_0"]] <- Fprod0
  DLMCMTest_output[["pval_min"]] <- pval_min
  DLMCMTest_output[["pval_prod"]] <- pval_prod
  return(DLMCMTest_output)
}





# ==============================================================================
#' @title Maximized Monte Carlo Likelihood Ratio Test
#' 
#' @description Performs the Maximized Monte Carlo (MMC) moment-based test described in Dufour & Luger (2017). 
#' 
#' @param Y Series to be tested
#' @param ar Order of autoregressive components AR(p)
#' @param pval_type Method to be used to combine p-values (see eq. (17) and (18) of Dufour & Luger (2017))
#' @param control List of control parameters. See “Details”.
#'
#' @return 
#' 
#' @export
DLMMCTest <- function(Y, ar, pval_type = "min", control = list()){
  # ----- Set control values
  con <- list(N = 99,
              simdist_N = 10000,
              getSE = TRUE,
              eps = 0.1,
              CI_union = TRUE,
              lambda = 100,
              stationary_ind = TRUE,
              type = "GenSA",
              silence = FALSE,
              threshold_stop = 1,
              type_control = list(maxit = 200))
  # Perform some checks for controls
  nmsC <- names(con)
  con[(namc <- names(control))] <- control
  if(length(noNms <- namc[!namc %in% nmsC])){
    arning("unknown names in control: ", paste(noNms,collapse=", ")) 
  }
  # ----- Set initial value as consistent estimate
  if(ar>0){
    mdl_h0 <- ARmdl(Y, ar, TRUE, con[["getSE"]])
    y <- mdl_h0[["y"]]
    x <- mdl_h0[["x"]]
    theta_0 <- mdl_h0[["phi"]]
  }else{
    stop("No Nuisance parameters is model is not Autoregressive. Number of lags must be greater than 0.")
  }
  # ----- Define lower & upper bounds for MMC search
  theta_low <- theta_0 - con[["eps"]]
  theta_upp <- theta_0 + con[["eps"]]
  # if CI_union==TRUE use union of eps & 2*standard error to define bounds
  phiSE <- mdl_h0[["theta_stderr"]][3:(3+ar-1)]
  if ((con[["CI_union"]]==TRUE) & all(is.finite(phiSE))){
    theta_low <- apply(cbind(as.matrix(theta_0 - 2*phiSE), as.matrix(theta_low)), 1, FUN = min)
    theta_upp <- apply(cbind(as.matrix(theta_0 + 2*phiSE), as.matrix(theta_upp)), 1, FUN = max)
  }
  # ----- Search for Max p-value within bounds
  MMCLRTest_output <- list()
  if(con[["type"]]=="pso"){
    # Set PSO specific controls
    con$type_control[["trace.stats"]] <- TRUE
    con$type_control[["trace"]] <- as.numeric(con[["silence"]]==FALSE)
    con$type_control[["abstol"]] <- -con[["threshold_stop"]]
    # begin optimization
    mmc_out <- pso::psoptim(par = theta_0, fn = DLMMCpval_fun, lower = theta_low, upper = theta_upp, 
                            gr = NULL, control = con[["type_control"]],
                            y = y, x = x, N = con[["N"]], simdist_N = con[["simdist_N"]], pval_type = pval_type,
                            stationary_ind = con[["stationary_ind"]], lambda = con[["lambda"]])
    MMCLRTest_output[["theta"]] <- mmc_out$par
    MMCLRTest_output[["pval"]] <- -mmc_out$value
  }else if(con[["type"]]=="GenSA"){
    # Set GenSA specific controls
    con$type_control[["verbose"]] <- con[["silence"]]==FALSE
    con$type_control[["threshold.stop"]] <- -con[["threshold_stop"]]
    mmc_out <- GenSA::GenSA(par = theta_0, fn = DLMMCpval_fun, lower = theta_low, upper = theta_upp, 
                            control = con[["type_control"]],
                            y = y, x = x, N = con[["N"]], simdist_N = con[["simdist_N"]], pval_type = pval_type,
                            stationary_ind = con[["stationary_ind"]], lambda = con[["lambda"]])
    MMCLRTest_output[["theta"]] <- mmc_out$par
    MMCLRTest_output[["pval"]] <- -mmc_out$value
  }else if(con[["type"]]=="GA"){
    mmc_out <- GA::ga(type = "real-valued", fitness = DLMMCpval_fun_max, 
                      y = y, x = x, N = con[["N"]], simdist_N = con[["simdist_N"]], pval_type = pval_type,
                      stationary_ind = con[["stationary_ind"]], lambda = con[["lambda"]],
                      lower = theta_low, upper = theta_upp, 
                      maxiter = con$type_control[["maxit"]], maxFitness = con[["threshold_stop"]], 
                      monitor = (con[["silence"]]==FALSE), suggestions = t(theta_0))
    MMCLRTest_output[["theta"]] <- c(mmc_out@solution)
    MMCLRTest_output[["pval"]] <- mmc_out@fitnessValue
  }else if(con[["type"]]=="gridSearch"){
    # Grid Search: not ready
  }
  MMCLRTest_output[["mdl_h0"]] <- mdl_h0
  MMCLRTest_output[["theta_0"]] <- theta_0
  MMCLRTest_output[["theta_low"]] <- theta_low
  MMCLRTest_output[["theta_upp"]] <- theta_upp
  MMCLRTest_output[["opt_output"]] <- mmc_out
  return(MMCLRTest_output)
}




# ==============================================================================
#' @title Monte-Carlo Moment-based test for MS AR model
#'
#' This function performs the Local Monte-Carlo Moment-Based test for
#' MS AR models presented in Dufour & Luger (2017) (i.e when no nuissance 
#' parameters are present). 
#'
#' @param Y Series to be tested
#' @param ar Order of autoregressive components AR(p)
#' @param k0 number of regime sunder the null
#' @param k1 number of regimes under the alternative
#' @param msmu bool indication if mean changes with regime (if k0>1)
#' @param msvar bool indication if variance changes with regime (if k0>1)
#' @param control List of control parameters. See “Details”.
#' 
#' @details The control argument is a list that can supply any of the following components:
#' \itemize{
#'   
#' }
#' 
#' @return List with model and test results.
#' 
#' @references Dufour, J. M., & Luger, R. (2017). Identification-robust moment-based 
#' tests for Markov switching in autoregressive models. Econometric Reviews, 36(6-9), 713-727.
#' 
#' @export
PMCTest <- function(Y, ar, k0, k1, msmu =  TRUE, msvar = TRUE, control = list()){
  # ----- Set control values
  con <- list(N = 99,
              simdist_N = 10000,
              maxit = 500,
              thtol = 1e-6,
              getSE = FALSE,
              converge_check = NULL)
  # Perform some checks for controls
  nmsC <- names(con)
  con[(namc <- names(control))] <- control
  if(length(noNms <- namc[!namc %in% nmsC])){
    warning("unknown names in control: ", paste(noNms,collapse=", ")) 
  }
  # ----- Transform process (from eq. (20) to z_t(phi))
  Tsize <- length(Y)
  if(ar>0){
    if (k0==1){
      mdl_h0  <- ARmdl(Y, ar, TRUE, con[["getSE"]])
      mdl_h0[["iterations"]] <- 1
      y       <- mdl_h0[["y"]]
      x       <- mdl_h0[["x"]]
      phi     <- mdl_h0[["phi"]]
      z       <- y - x%*%phi
    }else{
      mdl_h0  <- MSARmdl(Y, ar, k0, msmu, msvar, con[["maxit"]], con[["thtol"]], con[["getSE"]]) 
      y       <- mdl_h0[["y"]]
      x       <- mdl_h0[["x"]]
      phi     <- mdl_h0[["phi"]]
      z       <- y - x%*%phi 
      if (msmu == FALSE){
        mdl_h0[["mu"]] = rep(mdl_h0[["mu"]], k0)
      }
      if (msvar == FALSE){
        mdl_h0[["stdev"]] = rep(mdl_h0[["stdev"]], k0)
      }
    }
  }else{
    if (k0==1){
      z                       <- Y
      mdl_h0                  <- list()
      mdl_h0[["k"]]           <- k0
      mdl_h0[["ar"]]          <- 0
      mdl_h0[["mu"]]          <- mean(z)
      mdl_h0[["sigma"]]       <- sum((z-mean(z))*(z-mean(z)))/(Tsize-1)
      mdl_h0[["iterations"]]  <- 1
    }else{
      # **** USE MSmdl (i.e., HMM since no autoregressive params) to obtain mu and sigma for each regime
    }
  }
  # Optional model convergence checks
  if (is.null(con[["converge_check"]])==FALSE){
    if ((con[["converge_check"]]=="null") & (mdl_h0[["iterations"]]==con[["maxit"]])){
      stop("Model under null hypothesis did not converge. Run again to use different initial values and/or increase 'maxit'")
    }
  }
  # --------- Get parameters from approximated distribution ----------
  params <- approxDistQ(Tsize-ar, con[["simdist_N"]], mdl_h0, k1)
  # -------------------------- Get P-Values --------------------------
  N     <- con[["N"]]
  eps   <- z - mean(z)
  Fmin  <- calc_Qmcstat(eps, N, params, mdl_h0, k1, "min")
  Fprod <- calc_Qmcstat(eps, N, params, mdl_h0, k1, "prod")
  # ----- Obtain p-value
  Fmin0     <- Fmin[N+1]
  Fprod0    <- Fprod[N+1]
  FminSim   <- Fmin[1:N]
  FprodSim  <- Fprod[1:N]
  pval_min  <- MCpval(Fmin0, FminSim, "geq")
  pval_prod <- MCpval(Fprod0, FprodSim, "geq")
  LMC_N     <- cbind(Fmin,Fprod)
  # ----- Organize output
  Qtest_output <- list()
  Qtest_output[["eps"]] <- eps
  Qtest_output[["mdl_h0"]] <- mdl_h0
  Qtest_output[["dist_params"]] <- params
  Qtest_output[["LMC_N"]] <- LMC_N
  Qtest_output[["Fmin_0"]] <- Fmin0
  Qtest_output[["Fprod_0"]] <- Fprod0
  Qtest_output[["pval_min"]] <- pval_min
  Qtest_output[["pval_prod"]] <- pval_prod
  return(Qtest_output)
}