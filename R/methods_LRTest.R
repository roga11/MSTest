# ==============================================================================
#' @title Monte Carlo Likelihood Ratio Test sample disttribution (parallel version)
#' 
#' @description 
#' @param 
#'
#' @return 
#' 
#' @export
LR_samp_dist_par <- function(mdl_h0, k1, msmu, msvar, N, maxit, thtol, burnin, max_init, dist_converge_iter, workers){ 
  N_worker_i <- matrix(rep(floor(N/workers),workers),workers,1)
  N_worker_i[1:(N-floor(N/workers)*(workers))] <- N_worker_i[1:(N-floor(N/workers)*(workers))] +1 
  LRN_all <- matrix(0,N,1)
  `%dopar%` <- foreach::`%dopar%`
  LRN_all <- foreach::foreach(wi=1:workers, .inorder = FALSE, .packages = "MSTest") %dopar% {
    LRN <- LR_samp_dist(mdl_h0, k1, msmu, msvar, N_worker_i[wi], maxit, thtol, burnin, max_init, dist_converge_iter)
    LRN
  }
  return(unlist(LRN_all))
}
  
  

# ==============================================================================
#' @title Monte Carlo Likelihood Ratio Test
#' 
#' @description 
#' @param 
#'
#' @return 
#' 
#' @export
MCLRTest <- function(Y, ar, k0, k1, msmu = TRUE, msvar = TRUE, control = list()){
  # ----- Set control values
  con <- list(N = 99,
              maxit = 500,
              thtol = 1e-6,
              burnin = 200,
              getSE = FALSE,
              converge_check = "both",
              max_init = 100,
              dist_converge_iter = 100,
              workers = 0)
  # Perform some checks for controls
  nmsC <- names(con)
  con[(namc <- names(control))] <- control
  if(length(noNms <- namc[!namc %in% nmsC])){
    warning("unknown names in control: ", paste(noNms,collapse=", ")) 
  }
  # ----- Compute test statistic (using observed data)
  if (k0==1){
    mdl_h0 <- ARmdl(Y, ar, TRUE, con[["getSE"]])
    mdl_h0[["iterations"]] <- 1
    mdl_h1 <- MSARmdl(Y, ar, k1, msmu, msvar, con[["maxit"]], con[["thtol"]], con[["getSE"]], con[["max_init"]])
  }else if (k0>1){
    mdl_h0 <- MSARmdl(Y, ar, k0, msmu, msvar, con[["maxit"]], con[["thtol"]], con[["getSE"]], con[["max_init"]])  
    mdl_h1 <- MSARmdl(Y, ar, k1, msmu, msvar, con[["maxit"]], con[["thtol"]], con[["getSE"]], con[["max_init"]])  
  }
  # Optional model convergence checks
  if (is.null(con[["converge_check"]])==FALSE){
    if ((con[["converge_check"]]=="null") & (mdl_h0[["iterations"]]==con[["maxit"]])){
      stop("Model under null hypothesis did not converge. Run again to use different initial values and/or increase 'maxit'")
    }
    if((con[["converge_check"]]=="alt") & (mdl_h1[["iterations"]]==con[["maxit"]])){
      stop("Model under alternative hypothesis did not converge. Run again to use different initial values and/or increase 'maxit'")
    }
    if ((con[["converge_check"]]=="both") & ((mdl_h0[["iterations"]]==con[["maxit"]]) | (mdl_h1[["iterations"]]==con[["maxit"]]))){
      stop("Model did not converge. Run again to use different initial values and/or increase 'maxit'")
    }
  }
  logL0 <- mdl_h0[["logLike"]]
  logL1 <- mdl_h1[["logLike"]]
  theta_h0 <- mdl_h0[["theta"]]
  theta_h1 <- mdl_h1[["theta"]]
  LRT_0 <- -2*(logL0-logL1)
  # ----- Perform check of test stat and model parameters
  if ((is.finite(LRT_0)==FALSE) | (LRT_0<0) | (any(is.finite(theta_h0)==FALSE)) | (any(is.finite(theta_h1)==FALSE))){
    stop("LRT_0 or model parameters are not finite. Run again to use different initial values")
  }
  # ----- Simulate sample null distribution
  if(con[["workers"]]>0){
    LRN <- LR_samp_dist_par(mdl_h0, k1, msmu, msvar, con[["N"]], con[["maxit"]], con[["thtol"]], con[["burnin"]], con[["max_init"]], con[["dist_converge_iter"]], con[["workers"]])
  }else{
    LRN <- LR_samp_dist(mdl_h0, k1, msmu, msvar, con[["N"]], con[["maxit"]], con[["thtol"]], con[["burnin"]], con[["max_init"]], con[["dist_converge_iter"]]) 
  }
  
  # ----- Compute p-value
  pval <- MCpval(LRT_0, LRN, "geq")
  # ----- Organize output
  MCLRTest_output <- list()
  MCLRTest_output[["mdl_h0"]] <- mdl_h0
  MCLRTest_output[["mdl_h1"]] <- mdl_h1
  MCLRTest_output[["LRT_0"]] <- LRT_0
  MCLRTest_output[["LRN"]] <- LRN
  MCLRTest_output[["pval"]] <- pval
  return(MCLRTest_output)
}


# ==============================================================================
#' @title MMC nuisance parameter bounds for univariate models 
#' 
#' @description 
#' @param 
#'
#' @return 
#' 
#' @export
MMC_bounds_univariate <- function(theta_0, mdl_h0, mdl_h1, con, msmu, msvar){
  k0 <- mdl_h0[["k"]]
  k1 <- mdl_h1[["k"]]
  # ----- Define lower & upper bounds for search
  theta_low = theta_0 - con[["eps"]]
  theta_upp = theta_0 + con[["eps"]]
  # create ball around union of eps and 2*standard error (if set to true and SE are finite)
  if ((con[["CI_union"]]==TRUE) & all(is.finite(mdl_h0[["theta_stderr"]])) & all(is.finite(mdl_h1[["theta_stderr"]]))){
    theta_low <- apply(cbind(as.matrix(theta_0 - 2*c(mdl_h0[["theta_stderr"]],mdl_h1[["theta_stderr"]])),as.matrix(theta_low)), 1, FUN = min)
    theta_upp <- apply(cbind(as.matrix(theta_0 + 2*c(mdl_h0[["theta_stderr"]],mdl_h1[["theta_stderr"]])),as.matrix(theta_upp)), 1, FUN = max)
  }
  # check that bounds respect admissible regions
  sigma_h0_ind <- rep(FALSE, length(mdl_h0[["theta"]]))
  P_h0_ind <- rep(FALSE, length(mdl_h0[["theta"]]))
  sigma_h1_ind <- rep(FALSE, length(mdl_h1[["theta"]]))
  P_h1_ind <- rep(FALSE, length(mdl_h1[["theta"]]))
  sigma_h0_ind[(2+msmu*(k0-1)):(2+msmu*(k0-1)+msvar*(k0-1))] <- TRUE
  sigma_h1_ind[(2+msmu*(k1-1)):(2+msmu*(k1-1)+msvar*(k1-1))] <- TRUE
  P_h1_ind[(length(mdl_h1[["theta"]])-k1*k1+1):length(mdl_h1[["theta"]])] <- TRUE
  if (k0>1){
    P_h0_ind[(length(mdl_h0[["theta"]])-k0*k0+1):length(mdl_h0[["theta"]])] <- TRUE
  }
  sigma_ind <- c(sigma_h0_ind,sigma_h1_ind)
  P_ind <- c(P_h0_ind,P_h1_ind)
  # correct variances to be in admissible region
  theta_low[sigma_ind][theta_low[sigma_ind]<=0]=con[["variance_lower_bound"]]
  # correct transition probs to be in admissible region
  theta_low[P_ind][theta_low[P_ind]<0] <- 0
  theta_upp[P_ind][theta_upp[P_ind]>1] <- 1
  mmc_bounds <- list()
  mmc_bounds[["theta_low"]] <- theta_low
  mmc_bounds[["theta_upp"]] <- theta_upp
  return(mmc_bounds)
}




# ==============================================================================
#' @title Maximized Monte Carlo Likelihood Ratio Test
#' 
#' @description 
#' @param 
#'
#' @return 
#' 
#' @export
MMCLRTest <- function(Y, ar, k0, k1, msmu = TRUE, msvar = TRUE, control = list()){
  # ----- Set control values
  con <- list(N = 99,
              maxit = 500,
              thtol = 1e-6,
              burnin = 200, 
              getSE = TRUE,
              eps = 0.1,
              CI_union = TRUE,
              lambda = 100,
              variance_lower_bound = 0.1,
              stationary_ind = TRUE,
              type = "GenSA",
              silence = FALSE,
              threshold_stop = 1,
              converge_check = NULL,
              max_init = 100, 
              dist_converge_iter = 100,
              type_control = list(maxit = 200))
  # Perform some checks for controls
  nmsC <- names(con)
  con[(namc <- names(control))] <- control
  if(length(noNms <- namc[!namc %in% nmsC])){
    arning("unknown names in control: ", paste(noNms,collapse=", ")) 
  }
  # ----- Set initial value as consistent estimate
  if (k0==1){
    mdl_h0 <- ARmdl(Y, ar, TRUE, con[["getSE"]])
    mdl_h0[["iterations"]] <- 1
    mdl_h1 <- MSARmdl(Y, ar, k1, msmu, msvar, con[["maxit"]], con[["thtol"]], con[["getSE"]], con[["max_init"]])
  }else if (k0>1){
    mdl_h0 <- MSARmdl(Y, ar, k0, msmu, msvar, con[["maxit"]], con[["thtol"]], con[["getSE"]], con[["max_init"]])  
    mdl_h1 <- MSARmdl(Y, ar, k1, msmu, msvar, con[["maxit"]], con[["thtol"]], con[["getSE"]], con[["max_init"]])  
  }
  # Optional model convergence checks (model under null is only checked if k0>1 since EM is only used then)
  if (is.null(con[["converge_check"]])==FALSE){
    if ((con[["converge_check"]]=="null") & (mdl_h0[["iterations"]]==con[["maxit"]])){
      stop("Model under null hypothesis did not converge. Run again to use different initial values and/or increase 'maxit'")
    }
    if((con[["converge_check"]]=="alt") & (mdl_h1[["iterations"]]==con[["maxit"]])){
      stop("Model under alternative hypothesis did not converge. Run again to use different initial values and/or increase 'maxit'")
    }
    if ((con[["converge_check"]]=="both") & ((mdl_h0[["iterations"]]==con[["maxit"]]) | (mdl_h1[["iterations"]]==con[["maxit"]]))){
      stop("Model did not converge. Run again to use different initial values and/or increase 'maxit'")
    }
  }
  theta_0 <- c(mdl_h0[["theta"]], mdl_h1[["theta"]])
  # ----- Define lower & upper bounds for search
  mmc_bounds <- MMC_bounds_univariate(theta_0, mdl_h0, mdl_h1, con, msmu, msvar)
  theta_low <- mmc_bounds[["theta_low"]]
  theta_upp <- mmc_bounds[["theta_upp"]]
  # ----- Search for Max p-value within bounds
  MMCLRTest_output <- list()
  if(con[["type"]]=="pso"){
    # Set PSO specific controls
    con$type_control[["trace.stats"]] <- TRUE
    con$type_control[["trace"]] <- as.numeric(con[["silence"]]==FALSE)
    con$type_control[["abstol"]] <- -con[["threshold_stop"]]
    # begin optimization
    mmc_out <- pso::psoptim(par = theta_0, fn = MMCLRpval_fun, lower = theta_low, upper = theta_upp, 
                            gr = NULL, control = con[["type_control"]],
                            mdl_h0 = mdl_h0, mdl_h1 = mdl_h1, msmu = msmu, msvar = msvar, ar = ar, 
                            N = con[["N"]], maxit = con[["maxit"]], thtol = con[["thtol"]], burnin = con[["burnin"]],
                            stationary_ind = con[["stationary_ind"]], lambda = con[["lambda"]], max_init = con[["max_init"]], dist_converge_iter = con[["dist_converge_iter"]])
    MMCLRTest_output[["theta"]] <- mmc_out$par
    MMCLRTest_output[["pval"]] <- -mmc_out$value
  }else if(con[["type"]]=="GenSA"){
    # Set GenSA specific controls
    con$type_control[["trace.mat"]] <- TRUE
    con$type_control[["verbose"]] <- con[["silence"]]==FALSE
    con$type_control[["threshold.stop"]] <- -con[["threshold_stop"]]
    # begin optimization
    mmc_out <- GenSA::GenSA(par = theta_0, fn = MMCLRpval_fun, lower = theta_low, upper = theta_upp, 
                            control = con[["type_control"]],
                            mdl_h0 = mdl_h0, mdl_h1 = mdl_h1, msmu = msmu, msvar = msvar, ar = ar, 
                            N = con[["N"]], maxit = con[["maxit"]], thtol = con[["thtol"]], burnin = con[["burnin"]], 
                            stationary_ind = con[["stationary_ind"]], lambda = con[["lambda"]], max_init = con[["max_init"]], dist_converge_iter = con[["dist_converge_iter"]])
    MMCLRTest_output[["theta"]] <- mmc_out$par
    MMCLRTest_output[["pval"]] <- -mmc_out$value
  }else if(con[["type"]]=="GA"){
    # begin optimization
    mmc_out <- GA::ga(type = "real-valued", fitness = MMCLRpval_fun_max, 
                      mdl_h0 = mdl_h0, mdl_h1 = mdl_h1, msmu = msmu, msvar = msvar, ar = ar, 
                      N = con[["N"]], maxit = con[["maxit"]], thtol = con[["thtol"]], burnin = con[["burnin"]],
                      stationary_ind = con[["stationary_ind"]], lambda = con[["lambda"]], max_init = con[["max_init"]], dist_converge_iter = con[["dist_converge_iter"]], 
                      lower = theta_low, upper = theta_upp, 
                      maxiter = con$type_control[["maxit"]], maxFitness = con[["threshold_stop"]], 
                      monitor = (con[["silence"]]==FALSE), suggestions = t(theta_0))
    MMCLRTest_output[["theta"]] <- c(mmc_out@solution)
    MMCLRTest_output[["pval"]] <- mmc_out@fitnessValue
  }else if(con[["type"]]=="gridSearch"){
    # Grid Search: not ready
  }
  MMCLRTest_output[["mdl_h0"]] <- mdl_h0
  MMCLRTest_output[["mdl_h1"]] <- mdl_h1
  MMCLRTest_output[["theta_0"]] <- theta_0
  MMCLRTest_output[["theta_low"]] <- theta_low
  MMCLRTest_output[["theta_upp"]] <- theta_upp
  MMCLRTest_output[["opt_output"]] <- mmc_out
  return(MMCLRTest_output)
}



# ==============================================================================
#' @title Bootstrap Likelihood Ratio Test
#' 
#' @description 
#' @param 
#'
#' @return 
#' 
#' @export
BootLRTest <- function(Y, ar, k0, k1, msmu = TRUE, msvar = TRUE, control = list()){
  # ----- Set control values
  con <- list(B = 1000,
            maxit = 500,
            thtol = 1e-6,
            burnin = 200,
            getSE = FALSE,
            converge_check = "both",
            max_init = 100,
            dist_converge_iter = 100)
  # Perform some checks for controls
  nmsC <- names(con)
  con[(namc <- names(control))] <- control
  if(length(noNms <- namc[!namc %in% nmsC])){
    warning("unknown names in control: ", paste(noNms,collapse=", ")) 
  }
  # ----- Compute test statistic (using observed data)
  if (k0==1){
    mdl_h0 <- ARmdl(Y, ar, TRUE, con[["getSE"]])
    mdl_h0[["iterations"]] <- 1
    mdl_h1 <- MSARmdl(Y, ar, k1, msmu, msvar, con[["maxit"]], con[["thtol"]], con[["getSE"]], con[["max_init"]])
  }else if (k0>1){
    mdl_h0 <- MSARmdl(Y, ar, k0, msmu, msvar, con[["maxit"]], con[["thtol"]], con[["getSE"]], con[["max_init"]])  
    mdl_h1 <- MSARmdl(Y, ar, k1, msmu, msvar, con[["maxit"]], con[["thtol"]], con[["getSE"]], con[["max_init"]])  
  }
  # Optional model convergence checks (model under null is only checked if k0>1 since EM is only used then)
  if (is.null(con[["converge_check"]])==FALSE){
    if ((con[["converge_check"]]=="null") & (mdl_h0[["iterations"]]==con[["maxit"]])){
      stop("Model under null hypothesis did not converge. Run again to use different initial values and/or increase 'maxit'")
    }
    if((con[["converge_check"]]=="alt") & (mdl_h1[["iterations"]]==con[["maxit"]])){
      stop("Model under alternative hypothesis did not converge. Run again to use different initial values and/or increase 'maxit'")
    }
    if ((con[["converge_check"]]=="both") & ((mdl_h0[["iterations"]]==con[["maxit"]]) | (mdl_h1[["iterations"]]==con[["maxit"]]))){
      stop("Model did not converge. Run again to use different initial values and/or increase 'maxit'")
    }
  }
  logL0 <- mdl_h0[["logLike"]]
  logL1 <- mdl_h1[["logLike"]]
  theta_h0 <- mdl_h0[["theta"]]
  theta_h1 <- mdl_h1[["theta"]]
  LRT_0 <- -2*(logL0-logL1)
  # ----- Perform check of test stat and model parameters
  if ((is.finite(LRT_0)==FALSE) | (LRT_0<0) | (any(is.finite(theta_h0)==FALSE)) | (any(is.finite(theta_h1)==FALSE))){
    stop("LRT_0 or model parameters are not finite. Please check series")
  }
  # ----- Simulate sample null distribution
  LRN <- LR_samp_dist(mdl_h0, k1, msmu, msvar, con[["B"]], con[["maxit"]], con[["thtol"]], con[["burnin"]], con[["max_init"]], con[["dist_converge_iter"]])
  # ----- Compute p-value
  pval <- sum(LRN>LRT_0)/con[["B"]] # [eq. 4.62] (Davidson & MacKinnon, 2004)
  # ----- Organize output
  BootLRTest_output<-list()
  BootLRTest_output[["mdl_h0"]] <- mdl_h0
  BootLRTest_output[["mdl_h1"]] <- mdl_h1
  BootLRTest_output[["LRT_0"]] <- LRT_0
  BootLRTest_output[["LRN"]] <- LRN
  BootLRTest_output[["pval"]] <- pval
  return(BootLRTest_output)
}


