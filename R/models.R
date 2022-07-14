

#' @title Markov-switching autoregressive model by Expectation Minimization algorithm
#' 
#' @description This function estimates a Markov-switching autoregressive model using the Expectation Minimization (EM) algorithm of Dempster, Laird and Rubin (1977) as explained in Hamilton (1990).
#' 
#' @param Y (Tx1) vector with observational data. Required argument.
#' @param ar integer for the number of lags to use in estimation. Must be greater than or equal to 0. Default is 0.
#' @param k integer for the number of regimes to use in estimation. Must be greater than or equal to 2. Default is 2.
#' @param control List with optimization options including:
#'     - msmu indicator for switch in mean (TRUE) or constant mean (FALSE). Default is TRUE.
#'     - msvar bool indicator for switch in variance (TRUE) or constant variance (FALSE). Default is TRUE.
#'     - maxit integer determining the maximum number of EM iterations.
#'     - thtol double determining the convergence criterion for the absolute difference in parameter estimates (theta) between iterations. Default is 1e-6.
#'     - getSE bool if 'TRUE' standard errors are computed and returned. If 'FALSE' standard errors are not computed. Default is 'FALSE'.
#'     - max_init integer determining the maximum number of initial values to attempt in case solution is NaN. Default is 500.
#'     - use_diff_init integer determining how many different initital values (that do not return NaN) to try. Default is 1. 
#'     - init_value vector of initial values. This is optional. Default is NUll, in which case 'initValsMS()' is used to generate initial values.
#' 
#' @return List with model characteristics
#' 
#' @references Dempster, A. P., N. M. Laird, and D. B. Rubin. 1977. “Maximum Likelihood from Incomplete Data via the EM Algorithm.” Journal of the Royal Statistical Society. Series B 39 (1): 1–38
#' @references Hamilton, James D. 1990. “Analysis of time series subject to changes in regime.” Journal of econometrics, 45 (1-2): 39–70
#' 
#' @export
MSmdl_EM <- function(Y, ar, k, control = list()){
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
  mdl_out <- ARmdl(Y, ar = ar, intercept = TRUE, getSE = con[["getSE"]])
  msmu <- con[["msmu"]]
  msvar <- con[["msvar"]]
  mdl_out[["msmu"]] = msmu
  mdl_out[["msvar"]] = msvar
  # ---------- Estimate model 
  EM_output_all <- list(con[["use_diff_init"]])
  max_loglik <- matrix(0, con[["use_diff_init"]], 1)
  if (is.null(con[["init_value"]])==FALSE){
    # ----- Estimate using initial values provided
    EM_output <- MS_EMest(con[["init_value"]], mdl_out, k, optim_options)
    EM_output[["theta_0"]] <- con[["init_value"]]
    EM_output[["init_used"]] <- 1
  }else{
    # ----- Estimate using 'use_diff_init' different initial values
    for (xi in 1:con[["use_diff_init"]]){
      init_used <- 0
      converge_check <- FALSE
      while ((converge_check==FALSE) & (init_used<con[["max_init"]])){
        # ----- Initial values
        theta_0 <- initValsMS(mdl_out, k)
        # ----- Estimate using EM algorithm 
        EM_output_tmp <- MS_EMest(theta_0, mdl_out, k, optim_options)
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
  theta_mu_ind <- c(rep(1, 1 + (k-1)*msmu), rep(0, 1 + (k-1)*msvar + ar + k*k))
  theta_sig_ind <- c(rep(0, 1 + (k-1)*msmu), rep(1, 1 + (k-1)*msvar), rep(0, ar + k*k))
  theta_phi_ind <- c(rep(0, 2 + (k-1)*msmu + (k-1)*msvar), rep(1, ar), rep(0, k*k))
  theta_P_ind <- c(rep(0, 2 + (k-1)*msmu + (k-1)*msvar + ar), rep(1, k*k))
  MSARmdl_output <- EM_output
  MSARmdl_output[["theta_mu_ind"]] = theta_mu_ind
  MSARmdl_output[["theta_sig_ind"]] = theta_sig_ind
  MSARmdl_output[["theta_phi_ind"]] = theta_phi_ind
  MSARmdl_output[["theta_P_ind"]] = theta_P_ind
  MSARmdl_output[["stdev"]] <- sqrt(EM_output[["sigma"]])
  MSARmdl_output[["y"]] <- mdl_out[["y"]]
  MSARmdl_output[["ar"]] <- ar
  MSARmdl_output[["n"]] <- mdl_out[["n"]]
  MSARmdl_output[["q"]] <- 1
  MSARmdl_output[["k"]] <- k
  MSARmdl_output[["x"]] <- mdl_out[["x"]]
  MSARmdl_output[["X"]] <- mdl_out[["X"]]
  MSARmdl_output[["msmu"]] <- con[["msmu"]]
  MSARmdl_output[["msvar"]] <- con[["msvar"]]
  MSARmdl_output[["control"]] <- con
  if (con[["getSE"]]==TRUE){
    Hess <- getHess(MSARmdl_output, k)
    info_mat <- solve(-Hess)
    nearPD_used <- FALSE
    if ((all(is.na(Hess)==FALSE)) & (any(diag(info_mat)<0))){
      info_mat <- nearPD(info_mat)
      nearPD_used <- TRUE
    }
    MSARmdl_output[["Hess"]] <- Hess
    MSARmdl_output[["theta_stderr"]] <- sqrt(diag(info_mat))
    MSARmdl_output[["info_mat"]] <- info_mat
    MSARmdl_output[["nearPD_used"]] <- nearPD_used
  }
  if (is.null(con[["init_value"]])){
    MSARmdl_output[["trace"]] <- EM_output_all
  }
  return(MSARmdl_output)
}




#' @title Markov-switching vector autoregressive model by Expectation Minimization algorithm
#' 
#' @description This function estimates a Markov-switching vector autoregressive model using the Expectation Minimization (EM) algorithm of Dempster, Laird and Rubin (1977) as explained in Krolzig (1997).
#' 
#' @param Y (Tx1) vector with observational data. Required argument.
#' @param ar integer for the number of lags to use in estimation. Must be greater than or equal to 0. Default is 0.
#' @param k integer for the number of regimes to use in estimation. Must be greater than or equal to 2. Default is 2.
#' @param control List with optimization options including:
#'     - msmu indicator for switch in mean (TRUE) or constant mean (FALSE). Default is TRUE.
#'     - msvar bool indicator for switch in variance (TRUE) or constant variance (FALSE). Default is TRUE.
#'     - maxit integer determining the maximum number of EM iterations.
#'     - thtol double determining the convergence criterion for the absolute difference in parameter estimates (theta) between iterations. Default is 1e-6.
#'     - getSE bool if 'TRUE' standard errors are computed and returned. If 'FALSE' standard errors are not computed. Default is 'FALSE'.
#'     - max_init integer determining the maximum number of initial values to attempt in case solution is NaN. Default is 500.
#'     - use_diff_init integer determining how many different initital values (that do not return NaN) to try. Default is 1. 
#'     - init_value vector of initial values. This is optional. Default is NUll, in which case 'initValsMS()' is used to generate initial values.
#' 
#' @return List with model characteristics
#' 
#' @references Dempster, A. P., N. M. Laird, and D. B. Rubin. 1977. “Maximum Likelihood from Incomplete Data via the EM Algorithm.” Journal of the Royal Statistical Society. Series B 39 (1): 1–38
#' @references Krolzig, Hans-Martin. 1997. “The markov-switching vector autoregressive model.” In Markov-Switching Vector Autoregressions, 6–28. Springer
#' @export
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
  theta_var_ind <- c(rep(0, q + q*(k-1)*msmu),rep(t(sig_mattovec(diag(q),q)), 1+(k-1)*msvar), rep(0, phi_len + k*k))
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

