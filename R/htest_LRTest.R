
#' @title Monte Carlo Likelihood Ratio Test sample distribution (parallel version)
#' 
#' @description This function simulates the sample distribution under the null hypothesis using a parallel pool.
#'  
#' @param mdl_h0 List with restricted model properties.
#' @param k1 integer specifying the number of regimes under the alternative hypothesis.
#' @param N integer specifying the number of replications.
#' @param burnin integer specifying the number of observations to drop from beginning of simulation.
#' @param mdl_h0_control List with controls/options used to estimate restricted model.
#' @param mdl_h1_control List with controls/options used to estimate unrestricted model.
#' @param workers Integer determining the number of workers to use for parallel computing version of test. Note that parallel pool must already be open.
#'  
#' @return vector of simulated LRT statistics
#' 
#' @keywords internal
#' 
#' @references Rodriguez-Rondon, Gabriel and Jean-Marie Dufour. 2022. "Simulation-Based Inference for Markov Switching Models” \emph{JSM Proceedings, Business and Economic Statistics Section: American Statistical Association}.
#' @references Rodriguez-Rondon, Gabriel and Jean-Marie Dufour. 2023. “Monte Carlo Likelihood Ratio Tests for Markov Switching Models.” \emph{Unpublished manuscript}.
#'
#' @export
LR_samp_dist_par <- function(mdl_h0, k1, N, burnin, Z, mdl_h0_control, mdl_h1_control, workers){ 
  # ----- Set number of simulations per worker
  N_i_f <- floor(N/workers)
  N_worker_i <- matrix(rep(N_i_f, workers), workers, 1)
  if (sum(N_worker_i)<N){
    N_worker_i[1:(N-(N_i_f*workers))] <- N_worker_i[1:(N-(N_i_f*workers))] + 1  
  }
  # ----- Begin parallel simulations
  LRN_all <- matrix(0,N,1)
  `%dopar%` <- foreach::`%dopar%`
  wi <- NULL # needed to pass CMD check
  LRN_all <- foreach::foreach(wi = 1:workers, .inorder = FALSE, .packages = "MSTest") %dopar% {
    LRN <- LR_samp_dist(mdl_h0, k1, N_worker_i[wi], burnin, Z, mdl_h0_control, mdl_h1_control) 
    LRN
  }
  return(unlist(LRN_all))
}



#' @title Estimate model for likelihood ratio test
#' 
#' @description This function is used by the Monte Carlo testing procedures 
#' to estimate restricted and unrestricted models.
#' 
#' @param Y Series to be tested. Must be a (\code{T x q}) matrix.
#' @param p integer specifying the number of autoregressive lags.
#' @param q integer specifying the number of series.
#' @param k integer specifying the number of regimes.
#' @param Z exogeneous regressors. Defualt is NULL.
#' @param control List with control options for model estimation. For default values, see description of model being estimated.
#' 
#' @return List with estimated model properties. 
#' 
#' @keywords internal
#' 
#' @export
estimMdl <- function(Y, p, q, k, Z = NULL, control = list()){
  if ((k==1) & (p==0)){
    # Normally distributed model
    control$const <- TRUE # forced to be TRUE for hypothesis testing
    mdl <- Nmdl(Y, Z, control)
    mdl$converged = TRUE
  }else if ((k>1) & (p==0)){
    # Hidden Markov model
    mdl <- HMmdl(Y, k, Z, control)
    mdl$converged = (mdl$deltath <= mdl$control$thtol)
  }else if ((k==1) & (q==1) & (p>0)){
    # Autoregressive model
    control$const <- TRUE # forced to be TRUE for hypothesis testing
    if (is.null(Z) | (length(Z)==0)){
      mdl <- ARmdl(Y, p, control)
      mdl$converged = TRUE 
    }else{
      mdl <- ARXmdl(Y, p, Z, control)
      mdl$converged = TRUE
    }
  }else if ((k>1) & (q==1) & (p>0)){
    # Markov switching model
    if (is.null(Z) | (length(Z)==0)){
      mdl <- MSARmdl(Y, p, k, control)
      mdl$converged = (mdl$deltath <= mdl$control$thtol)  
    }else{
      mdl <- MSARXmdl(Y, p, k, Z, control)
      mdl$converged = (mdl$deltath <= mdl$control$thtol)
    }
  }else if ((k==1) & (q>1) & (p>0)){
    # Vector autoregressive model
    control$const <- TRUE # forced to be TRUE for hypothesis testing
    if (is.null(Z) | (length(Z)==0)){
      mdl <- VARmdl(Y, p, control)
      mdl$converged = TRUE 
    }else{
      mdl <- VARXmdl(Y, p, Z, control)
      mdl$converged = TRUE 
    }
  }else if ((k>1) & (q>1) & (p>0)){
    # Vector autoregressive Markov switching model
    if (is.null(Z) | (length(Z)==0)){
      mdl <- MSVARmdl(Y, p, k, control)
      mdl$converged = (mdl$deltath <= mdl$control$thtol)
    }else{
      mdl <- MSVARXmdl(Y, p, k, Z, control)
      mdl$converged = (mdl$deltath <= mdl$control$thtol)
    }
  }
  return(mdl)
}


#' @title Monte Carlo Likelihood Ratio Test
#' 
#' @description This function performs the Local Monte Carlo likelihood ratio 
#' test (LMC-LRT) proposed in Rodriguez-Rondon & Dufour (2024). As discussed in 
#' their work, this test can be applied in very general settings and can be used 
#' to compare varioous regimes under the null and under the alternative. 
#' 
#' @param Y  Series to be tested. Must be a (\code{T x q}) matrix where T is the number of time observations and q is the number of variables.
#' @param p  Number of autoregressive lags. Must be greater than or equal to 0. 
#' @param k0 Number of regimes under null hypothesis. Must be greater than or equal to 1.
#' @param k1 Number of regimes under alternative hypothesis. Must be greater than \code{k0}.
#' @param Z Exogenous regressors. Optional input and default is NULL. When used, it should be a (\code{T x qz}) matrix where T is the number of time observations and q is the number of exogenous variables.
#' @param control List with test procedure options including: 
#' \itemize{
#'   \item N: Integer determining the number of Monte Carlo simulations. Default is set to \code{99} as in paper.
#'   \item burnin: Number of simulated observations to remove from beginning. Default is \code{100}.
#'   \item converge_check: String or NULL determining if convergence of model(s) should be verified. Allowed inputs are: "null", "alt", "both", or \code{NULL}. If \code{NULL} (default) no model convergence is verified.
#'   \item workers: Integer determining the number of workers to use for parallel computing version of test. Note that parallel pool must already be open. Default is \code{0}.
#'   \item mdl_h0_control: List with restricted model options. See \code{\link{Nmdl}}, \code{\link{ARmdl}}, \code{\link{VARmdl}}, \code{\link{HMmdl}}, \code{\link{MSARmdl}}, or \code{\link{MSVARmdl}} documentation for available and default values.
#'   \item mdl_h1_control: List with unrestricted model options. See \code{\link{HMmdl}}, \code{\link{MSARmdl}}, or \code{\link{MSVARmdl}} documentation for available and default values.
#'   \item use_diff_init_sim: Value which determines the number of initial values to use when estimating models for null distribution. Default is set to use the same as specified in \code{mdl_h0_control} and \code{mdl_h1_control}.
#' }
#'
#' @return List of class \code{LMCLRTest} (\code{S3} object) with attributes including: 
#' \itemize{
#'   \item mdl_h0: List with restricted model attributes. 
#'   \item mdl_h1: List with unrestricted model attributes. 
#'   \item LRT_0: Value of test statistic from observed data.
#'   \item LRN: A (\code{N x 1}) vector of test statistics from data simulated under the null hypothesis.
#'   \item pval: P-value of Local Monte Carlo Likelihood Ratio Test.
#'   \item LRN_cv: Vector with 90\%, 95\%, and 99\% Monte Carlo simulated critical values (from vector \code{LRN}). These are not asymptotic critical values. 
#'   \item control: List with test procedure options used.
#' }
#'
#' @references Rodriguez-Rondon, Gabriel and Jean-Marie Dufour. 2022. "Simulation-Based Inference for Markov Switching Models” \emph{JSM Proceedings, Business and Economic Statistics Section: American Statistical Association}.
#' @references Rodriguez-Rondon, Gabriel and Jean-Marie Dufour. 2025. “Monte Carlo Likelihood Ratio Tests for Markov Switching Models.” \emph{Unpublished manuscript}.
#' @example /inst/examples/LMCLRTest_examples.R
#' @export
LMCLRTest <- function(Y, p, k0, k1, Z = NULL, control = list()){
  # ----- Set control values
  con <- list(N = 99,
              burnin = 100,
              converge_check = NULL,
              workers = 0,
              mdl_h0_control = list(),
              mdl_h1_control = list(),
              use_diff_init_sim = NULL)
  # ----- Perform some checks for controls
  nmsC <- names(con)
  con[(namc <- names(control))] <- control
  if (length(noNms <- namc[!namc %in% nmsC])){
    warning("unknown names in control: ", paste(noNms,collapse=", ")) 
  }
  # ----- Perform other checks
  if (is.matrix(Y)){
    q <- ncol(Y)
  }else{
    stop("Observations Y must be a (T x q) matrix.") 
  }
  # ----- Estimate models using observed data
  mdl_h0 <- estimMdl(Y, p, q, k0, Z, con$mdl_h0_control)
  mdl_h1 <- estimMdl(Y, p, q, k1, Z, con$mdl_h1_control)
  con$mdl_h0_control <- mdl_h0$control
  con$mdl_h1_control <- mdl_h1$control
  # ----- Optional model convergence checks
  if (is.null(con$converge_check)==FALSE){
    if ((con$converge_check=="null") & (mdl_h0$converged==FALSE)){
      stop("Model under null hypothesis did not converge. Run again to use different initial values and/or increase 'maxit' for restricted model.")
    }
    if ((con$converge_check=="alt") & (mdl_h1$converged==FALSE)){
      stop("Model under alternative hypothesis did not converge. Run again to use different initial values and/or increase 'maxit' for unrestricted model.")
    }
    if ((con$converge_check=="both") & ((mdl_h0$converged==FALSE) | (mdl_h1$converged==FALSE))){
      stop("Model did not converge. Run again to use different initial values and/or increase 'maxit' for each models.")
    }
  }
  # ----- Compute test statistic (LRT_0)
  logL0 <- mdl_h0$logLike
  logL1 <- mdl_h1$logLike
  theta_h0 <- mdl_h0$theta
  theta_h1 <- mdl_h1$theta
  LRT_0 <- -2*(logL0-logL1)
  # ----- Perform check of test stat and model parameters
  if ((is.finite(LRT_0)==FALSE) | (any(is.finite(theta_h0)==FALSE)) | (any(is.finite(theta_h1)==FALSE))){
    stop("LRT_0 or model parameters are not finite. Run again to use different initial values") 
  }
  if (LRT_0<0){
    stop("LRT_0 is negative. Run again to use different initial values")
  }
  names(LRT_0) <- c("LRT_0")
  # ----- Simulate sample null distribution
  if (is.null(Z)==FALSE){
    Zsim <- Z[(p+1):nrow(Z),,drop=F]
  }else{
    Zsim <- Z
  }
  mdl_h0_null_cont <- con$mdl_h0_control
  mdl_h1_null_cont <- con$mdl_h1_control
  if (is.null(con$use_diff_init_sim)==FALSE){
    mdl_h0_null_cont$use_diff_init <- con$use_diff_init_sim
    mdl_h1_null_cont$use_diff_init <- con$use_diff_init_sim
  }
  if (con$workers>0){
    LRN <- LR_samp_dist_par(mdl_h0, k1, con$N, con$burnin, Zsim, mdl_h0_null_cont, mdl_h1_null_cont, con$workers)
  }else{
    LRN <- LR_samp_dist(mdl_h0, k1, con$N, con$burnin, Zsim, mdl_h0_null_cont, mdl_h1_null_cont) 
  }
  # ----- get critical values
  LRN     <- as.matrix(sort(LRN))
  LRN_cv  <- LRN[round(c(0.90,0.95,0.99)*nrow(LRN)),]
  names(LRN_cv)  <- paste0(c("0.90","0.95","0.99"), "%")
  # ----- Compute p-value
  pval <- MCpval(LRT_0, LRN, "geq")
  # ----- Organize output
  MCLRTest_output <- list(mdl_h0 = mdl_h0, mdl_h1 = mdl_h1, LRT_0 = LRT_0, LRN = LRN,
                          pval = pval, LRN_cv = LRN_cv, control = con)
  class(MCLRTest_output) <- "LMCLRTest"
  return(MCLRTest_output)
}


#' @title MMC nuisance parameter bounds 
#' 
#' @description This function is used to determine the lower and upper bounds for the MMC LRT parameter search.
#' 
#' @param mdl_h0 List with restricted model properties.
#' @param con List with control options provided to MMC LRT procedure.
#' 
#' @return List with \code{theta_low}, vector of parameter lower bounds, and \code{theta_upp}, vector of parameter upper bounds.
#' 
#' @keywords internal
#' 
#' @references Rodriguez-Rondon, Gabriel and Jean-Marie Dufour. 2022. "Simulation-Based Inference for Markov Switching Models” \emph{JSM Proceedings, Business and Economic Statistics Section: American Statistical Association}.
#' @references Rodriguez-Rondon, Gabriel and Jean-Marie Dufour. 2025. “Monte Carlo Likelihood Ratio Tests for Markov Switching Models.” \emph{Unpublished manuscript}.
#' 
#' @export
MMC_bounds <- function(mdl_h0, con){
  theta_0   <- mdl_h0$theta
  k0        <- mdl_h0$k
  # ----- Define lower & upper bounds for search
  theta_low = theta_0 - con$eps
  theta_upp = theta_0 + con$eps
  # create ball around union of eps and 2*standard error (if set to true and SE are finite)
  if ((con$CI_union==TRUE) & all(is.finite(mdl_h0$theta_se))){
    theta_low <- apply(cbind(as.matrix(theta_0 - 2*c(mdl_h0$theta_se)), as.matrix(theta_low)), 1, FUN = min)
    theta_upp <- apply(cbind(as.matrix(theta_0 + 2*c(mdl_h0$theta_se)), as.matrix(theta_upp)), 1, FUN = max)
  }
  # ----- Check that bounds respect admissible regions
  # correct lower bound of variances to be in admissible region
  sigma_ind <- mdl_h0$theta_var_ind
  if (any(theta_low[sigma_ind==1]<=0)==TRUE){
    theta_low[sigma_ind==1][theta_low[sigma_ind==1]<=0] = theta_0[sigma_ind==1][theta_low[sigma_ind==1]<=0]*con$variance_constraint  
  }
  # correct transition probability bounds to be in admissible region
  if (k0>1){
    P_h0_ind <- mdl_h0$theta_P_ind
    theta_low[P_h0_ind==1][theta_low[P_h0_ind==1]<con$P_low] <- con$P_low
    theta_upp[P_h0_ind==1][theta_upp[P_h0_ind==1]>con$P_upp] <- con$P_upp
  }
  if (mdl_h0$p>0){
    # correct lower and upper bounds for autoregressive parameters, if any
    phi_ind <- mdl_h0$theta_phi_ind
    if (is.null(con$phi_low)==FALSE){
      theta_low[phi_ind==1] <- apply(cbind(as.matrix(theta_low[phi_ind==1]),as.matrix(con$phi_low)), 1, function(x) max(x))
    }  
    if (is.null(con$phi_upp)==FALSE){
      theta_upp[phi_ind==1] <- apply(cbind(as.matrix(theta_upp[phi_ind==1]),as.matrix(con$phi_upp)), 1, function(x) min(x))
    }  
  }
  # ----- output
  mmc_bounds <- list(theta_low = theta_low, theta_upp = theta_upp)
  return(mmc_bounds)
}




#' @title Maximized Monte Carlo Likelihood Ratio Test
#'
#' @description This function performs the Maximized Monte Carlo likelihood ratio 
#' test (MMC-LRT) proposed in Rodriguez-Rondon & Dufour (2024).
#' 
#' @param Y  Series to be tested. Must be a (\code{T x q}) matrix  where T is the number of time observations and q is the number of variables.
#' @param p  Number of autoregressive lags. Must be greater than or equal to 0. 
#' @param k0 Number of regimes under null hypothesis. Must be greater than or equal to 1.
#' @param k1 Number of regimes under alternative hypothesis. Must be greater than \code{k0}.
#' @param Z  Exogenous regressors. Optional input and default is NULL. When used, it should be a (\code{T x qz}) matrix where T is the number of time observations and q is the number of exogenous variables.
#' @param control List with test procedure options including: 
#' \itemize{
#'   \item N: Integer determining the number of Monte Carlo simulations. Default is set to \code{99} as in paper.
#'   \item burnin: Number of simulated observations to remove from beginning. Default is \code{100}.
#'   \item converge_check: String of NULL determining if convergence of model(s) should be verified. Allowed inputs are: "null", "alt", "both", or \code{NULL}. If \code{NULL} (default) no model convergence is verified.
#'   \item workers: Integer determining the number of workers to use for parallel computing version of test. Note that parallel pool must already be open. Default is \code{0}.
#'   \item type: String that determines the type of optimization algorithm used. Arguments allowed are: \code{"pso"}, \code{"GenSA"}, and \code{"GA"}. Default is \code{"pso"}.
#'   \item eps: Double determining the constant value that defines a consistent set for search. Default is \code{0.1}.
#'   \item CI_union: Boolean determining if union of set determined by \code{eps} and confidence set should be used to define consistent set for search. Default is \code{TRUE}.
#'   \item lambda: Double determining penalty on nonlinear constraint. Default is \code{100}.
#'   \item stationary_constraint: Boolean determining if only stationary solutions are considered (if \code{TRUE}) or not (if \code{FALSE}). Default is \code{TRUE}.
#'   \item phi_low: Vector with lower bound for autoregressive parameters when optimizing. Default is \code{NULL}.
#'   \item phi_upp: Vector with upper bound for autoregressive parameters when optimizing. Default is \code{NULL}.
#'   \item P_low: Value with lower bound for transition probabilities when optimizing. Default is \code{0}.
#'   \item P_upp: Value with upper bound for transition probabilities when optimizing. Default is \code{1}.
#'   \item variance_constraint: Double used to determine the lower bound for variance in parameter set for search. Value should be between \code{0} and \code{1} as it is multiplied by consistent point estimates of variances. Default is \code{0.01} (i.e., \code{1\%} of consistent point estimates.
#'   \item silence: Boolean determining if optimization steps should be silenced (if \code{TRUE}) or not (if \code{FALSE}). Default is \code{FALSE}.
#'   \item threshold_stop: Double determining the global optimum of function. Default is \code{1}.
#'   \item mdl_h0_control: List with restricted model options. See \code{\link{Nmdl}}, \code{\link{ARmdl}}, \code{\link{VARmdl}}, \code{\link{HMmdl}}, \code{\link{MSARmdl}}, or \code{\link{MSVARmdl}} documentation for available and default values.
#'   \item mdl_h1_control: List with unrestricted model options. See \code{\link{HMmdl}}, \code{\link{MSARmdl}}, or \code{\link{MSVARmdl}} documentation for available and default values.
#'   \item use_diff_init_sim: Value which determines the number of initial values to use when estimating models for null distribution. Default is set to use the same as specified in \code{mdl_h0_control} and \code{mdl_h1_control}.
#'   \item optim_control: List with optimization algorithm options. See \code{\link[pso]{psoptim}}, \code{\link[GenSA]{GenSA}}, \code{\link[GA]{ga}}. Default is to set \code{list(maxit = 200)} so that maximum number of iterations is \code{200}.
#' }
#'
#' @return List of class \code{LMCLRTest} (\code{S3} object) with attributes including: 
#' \itemize{
#'   \item mdl_h0: List with restricted model attributes. 
#'   \item mdl_h1: List with unrestricted model attributes. 
#'   \item LRT_0: Value of test statistic from observed data.
#'   \item LRN: A (\code{N x 1}) vector of test statistics from data simulated under the null hypothesis.
#'   \item pval: P-value of Local Monte Carlo Likelihood Ratio Test.
#'   \item LRN_cv: Vector with 90\%, 95\%, and 99\% Monte Carlo simulated critical values (from vector \code{LRN}). These are not asymptotic critical values. 
#'   \item control: List with test procedure options used.
#' }
#'
#' @references Rodriguez-Rondon, Gabriel and Jean-Marie Dufour. 2022. "Simulation-Based Inference for Markov Switching Models” \emph{JSM Proceedings, Business and Economic Statistics Section: American Statistical Association}.
#' @references Rodriguez-Rondon, Gabriel and Jean-Marie Dufour. 2025. “Monte Carlo Likelihood Ratio Tests for Markov Switching Models.” \emph{Unpublished manuscript}.
#' @example /inst/examples/MMCLRTest_examples.R
#' @export
MMCLRTest <- function(Y, p, k0, k1, Z = NULL, control = list()){
  # ----- Set control values
  con <- list(N = 99,
              burnin = 100,
              converge_check = NULL,
              workers = 0,
              type = "pso",
              eps = 0.1,
              CI_union = TRUE,
              lambda = 100,
              stationary_constraint = TRUE,
              phi_low = NULL,
              phi_upp = NULL,
              P_low = 0,
              P_upp = 1,
              variance_constraint = 0.01,
              silence = FALSE,
              threshold_stop = 1,
              mdl_h0_control = list(getSE = TRUE),
              mdl_h1_control = list(getSE = TRUE),
              use_diff_init_sim = NULL,
              maxit = 50,
              optim_control = list())
  # ----- Perform some checks for controls
  nmsC <- names(con)
  con[(namc <- names(control))] <- control
  if(length(noNms <- namc[!namc %in% nmsC])){
    warning("unknown names in control: ", paste(noNms,collapse=", ")) 
  }
  # ----- Perform other checks
  if (is.matrix(Y)){
    q <- ncol(Y)
  }else{
    stop("Observations Y must be a (T x q) matrix.") 
  }
  if ((con$CI_union==TRUE) & ((con$mdl_h0_control$getSE==FALSE) | (con$mdl_h1_control$getSE==FALSE))){
    con$mdl_h0_control$getSE <- TRUE
    con$mdl_h1_control$getSE <- TRUE
    warning("getSE was changed to be 'TRUE' because CI_union is 'TRUE'.")
  }
  # ----- Estimate models using observed data
  mdl_h0 <- estimMdl(Y, p, q, k0, Z, con$mdl_h0_control)
  mdl_h1 <- estimMdl(Y, p, q, k1, Z, con$mdl_h1_control)
  con$mdl_h0_control <- mdl_h0$control
  con$mdl_h1_control <- mdl_h1$control
  # ----- Optional model convergence checks
  if (is.null(con$converge_check)==FALSE){
    if ((con$converge_check=="null") & (mdl_h0$converged==FALSE)){
      stop("Model under null hypothesis did not converge. Run again to use different initial values and/or increase 'maxit' for restricted model.")
    }
    if ((con$converge_check=="alt") & (mdl_h1$converged==FALSE)){
      stop("Model under alternative hypothesis did not converge. Run again to use different initial values and/or increase 'maxit' for unrestricted model.")
    }
    if ((con$converge_check=="both") & ((mdl_h0$converged==FALSE) | (mdl_h1$converged==FALSE))){
      stop("Model did not converge. Run again to use different initial values and/or increase 'maxit' for each models.")
    }
  }
  theta_0 <- mdl_h0$theta
  # ----- Define lower & upper bounds for search
  mmc_bounds <- MMC_bounds(mdl_h0, con)
  theta_low <- mmc_bounds$theta_low
  theta_upp <- mmc_bounds$theta_upp
  # ----- Search for Max p-value within bounds
  mdl_h0_null_cont <- con$mdl_h0_control
  mdl_h1_null_cont <- con$mdl_h1_control
  if (is.null(con$use_diff_init_sim)==FALSE){
    mdl_h0_null_cont$use_diff_init <- con$use_diff_init_sim
    mdl_h1_null_cont$use_diff_init <- con$use_diff_init_sim
  }
  if (is.null(Z)==FALSE){
    Zsim <- Z[(p+1):nrow(Z),,drop=F]
    exog <- TRUE
  }else{
    Zsim <- Z
    exog <- FALSE
  }
  if (con$type=="pso"){
    # Set PSO specific controls
    con$optim_control$trace.stats <- TRUE
    con$optim_control$trace <- as.numeric(con$silence==FALSE)
    con$optim_control$abstol <- -con$threshold_stop
    con$optim_control$maxf <- con$maxit
    con$optim_control$REPORT <- 1
    # begin optimization
    mmc_out   <- pso::psoptim(par = theta_0, fn = MMCLRpval_fun_min, lower = theta_low, upper = theta_upp, 
                              gr = NULL, control = con$optim_control,
                              mdl_h0 = mdl_h0, k1 = k1, LT_h1 = mdl_h1$logLike, N = con$N, burnin = con$burnin, workers = con$workers,
                              lambda = con$lambda, stationary_constraint = con$stationary_constraint, 
                              thtol = mdl_h1$control$thtol, Z = Zsim, exog = exog, mdl_h0_control = mdl_h0_null_cont, 
                              mdl_h1_control = mdl_h1_null_cont)
    theta     <- mmc_out$par
    pval      <- -mmc_out$value
  }else if(con$type=="GenSA"){
    # Set GenSA specific controls
    con$optim_control$trace.mat <- TRUE
    con$optim_control$verbose <- con$silence==FALSE
    con$optim_control$threshold.stop <- -con$threshold_stop
    con$optim_control$max.call <- con$maxit
    # begin optimization
    mmc_out   <- GenSA::GenSA(par = theta_0, fn = MMCLRpval_fun_min, lower = theta_low, upper = theta_upp, 
                              control = con$optim_control,
                              mdl_h0 = mdl_h0, k1 = k1, LT_h1 = mdl_h1$logLike, N = con$N, burnin = con$burnin, workers = con$workers,
                              lambda = con$lambda, stationary_constraint = con$stationary_constraint, 
                              thtol = mdl_h1$control$thtol, Z = Zsim, exog = exog, mdl_h0_control = mdl_h0_null_cont, 
                              mdl_h1_control = mdl_h1_null_cont)
    theta     <- mmc_out$par
    pval      <- -mmc_out$value
  }else if(con$type=="GA"){
    # begin optimization
    mmc_out   <- GA::ga(type = "real-valued", fitness = MMCLRpval_fun, 
                      mdl_h0 = mdl_h0, k1 = k1, LT_h1 = mdl_h1$logLike, N = con$N, burnin = con$burnin, workers = con$workers,
                      lambda = con$lambda, stationary_constraint = con$stationary_constraint, 
                      thtol = mdl_h1$control$thtol, Z = Zsim,exog = exog,  mdl_h0_control = mdl_h0_null_cont, 
                      mdl_h1_control = mdl_h1_null_cont,
                      lower = theta_low, upper = theta_upp, 
                      maxiter = con$maxit, maxFitness = con$threshold_stop, 
                      monitor = (con$silence==FALSE), suggestions = t(theta_0))
    theta     <- as.matrix(mmc_out@solution[1,])
    pval      <- mmc_out@fitnessValue
  }else if(con$type=="gridSearch"){
    stop("Optim method 'gridSearch' is not available yet. Please use 'pso', 'GenSA', or 'GA' for 'type' in control List. ")
    # LT_h1 <- mdl_h1$logLike
    # LRT_0s <- matrix(0,con$maxit,1)
    # mmc_params_h0 <- matrix(0,con$maxit,length(theta_0))
    # for (xp in 1:length(theta_0)){
    #   mmc_params_h0[,xp] <- runif(con$maxit,min = theta_low[xp], max = theta_upp[xp])  
    # }
    # # Need to write soemthing that will make sure process is stationary, P has columns that sum to 1
    # mmc_pval_mat <- matrix(0,con$maxit,1)
    # LRN_ls <- list()
    # for (xs in 1:nrow(mmc_params_h0)){
    #    mdl_h0_tmp <- mdledit(mdl_h0,mmc_params_h0[xs,],p,q,k0,exog)
    #    LRT_0s[xs,]  <- compu_tstat(mmc_params_h0[xs,], mdl_h0_tmp, LT_h1, p, q, k0, exog)
    #    if (con$workers>0){
    #      LRN <- LR_samp_dist_par(mdl_h0_tmp, k1, con$N, con$burnin, Zsim, mdl_h0_null_cont, mdl_h1_null_cont, con$workers)
    #    }else{
    #      LRN <- LR_samp_dist(mdl_h0, k1, con$N, con$burnin, Zsim, mdl_h0_null_cont, mdl_h1_null_cont) 
    #    }
    #    LRN_ls[[xs]] <- LRN
    #    mmc_pval_mat[xs,] <- MCpval(LRT_0s[xs,],LRN)
    #    if (mmc_pval_mat[xs,]>con$threshold_stop){
    #      break
    #    }
    # }
    # pval <- mmc_pval_mat[which.max(mmc_pval_mat)[1],]
    # theta <- mmc_params_h0[which.max(mmc_pval_mat)[1],]
    # LRT_0 <- LRT_0s[which.max(mmc_pval_mat)[1],]
  }
  # ----- get test output using optimization output params
  theta_h0 <- theta
  theta_h1 <- mdl_h1$theta
  names(theta_h0) <- names(mdl_h0$theta)
  names(theta_h1) <- names(mdl_h1$theta)
  mdl_h0_mmc <- mdledit(mdl_h0, theta_h0, p, q, k0, exog)
  mdl_h0_mmc$logLike <- logLik(mdl_h0_mmc)
  mdl_h0_mmc$AIC <- stats::AIC(mdl_h0_mmc)
  mdl_h0_mmc$BIC <- stats::BIC(mdl_h0_mmc)
  if (mdl_h0$control$getSE==TRUE){
    mdl_h0_mmc <- thetaSE(mdl_h0_mmc)
  }
  # Compute test stats
  LRT_0 = compu_tstat(theta_h0, mdl_h0, mdl_h1$logLike, p, q, k0, exog)
  names(LRT_0) <- c("LRT_0")
  # ----- organize test output
  MMCLRTest_output <- list(mdl_h0 = mdl_h0, mdl_h1 = mdl_h1, mdl_h0_mmc = mdl_h0_mmc, mdl_h1_mmc = mdl_h1, 
                           LRT_0 = LRT_0, pval = pval,
                           theta_h0 = theta_h0, theta_h1 = theta_h1, control = con, 
                           mmc_optimout = mmc_out)
  class(MMCLRTest_output) <- "MMCLRTest"
  return(MMCLRTest_output)
}



