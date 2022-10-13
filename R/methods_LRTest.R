
#' @title Monte Carlo Likelihood Ratio Test sample distribution (parallel version)
#' 
#' @keywords internal
#'
#' @export
LR_samp_dist_par <- function(mdl_h0, k1, N, burnin, mdl_h0_control, mdl_h1_control, workers){ 
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
    LRN <- LR_samp_dist(mdl_h0, k1, N_worker_i[wi], burnin, mdl_h0_control, mdl_h1_control) 
    LRN
  }
  return(unlist(LRN_all))
}



#' @title Estimate model for likelihood ratio test
#' 
#' @keywords internal
#' 
#' @export
estimMdl <- function(Y, p, q, k, control = list()){
  if ((k==1) & (p==0)){
    # Normally distributed model
    control$const <- TRUE # forced to be TRUE for testing
    mdl <- Nmdl(Y, control)
    mdl$p <- 0
    mdl$converged = TRUE
  }else if ((k>1) & (p==0)){
    # Hidden Markov model
    mdl <- HMmdl(Y, k, control)
    mdl$p <- 0
    mdl$converged = (mdl$deltath <= mdl$control$thtol)
  }else if ((k==1) & (q==1) & (p>0)){
    # Autoregressive model
    control$const <- TRUE # forced to be TRUE for testing
    mdl <- ARmdl(Y, p, control)
    mdl$converged = TRUE
  }else if ((k>1) & (q==1) & (p>0)){
    # Markov switching model
    mdl <- MSARmdl(Y, p, k, control)
    mdl$converged = (mdl$deltath <= mdl$control$thtol)
  }else if ((k==1) & (q>1) & (p>0)){
    # Vector autoregressive model
    control$const <- TRUE # forced to be TRUE for testing
    mdl <- VARmdl(Y, p, control)
    mdl$converged = TRUE
  }else if ((k>1) & (q>1) & (p>0)){
    # Vector autoregressive Markov switching model
    mdl <- MSVARmdl(Y, p, k, control)
    mdl$converged = (mdl$deltath <= mdl$control$thtol)
  }
  return(mdl)
}


#' @title Monte Carlo Likelihood Ratio Test
#' 
#' @description This function performs the Local Monte Carlo likelihood ratio 
#' test (LMC-LRT) proposed in Rodriguez Rondon & Dufour (2022).
#' 
#' @param Y  Series to be tested. Must be a (\code{T x q}) matrix.
#' @param p  Number of autoregressive lags. Must be greater than or equal to 0. 
#' @param k0 Number of regimes under null hypothesis. Must be greater than or equal to 1.
#' @param k1 Number of regimes under alternative hypothesis. Must be greater than \code{k0}.
#' @param control List with test procedure options including: 
#' \itemize{
#'   \item{\code{N}: }{Integer determining the number of Monte Carlo simulations. Default is set to \code{99} as in paper.}
#'   \item{\code{burnin}: }{Number of simulated observations to remove from beginning. Default is \code{100}.}
#'   \item{\code{converge_check}: }{String of NULL determining if convergence of model(s) should be verified. Allowed inputs are: "null", "alt", "both", or \code{NULL}. If \code{NULL} (default) no model convergence is verified.}
#'   \item{\code{workers}: }{Integer determining the number of workers to use for parallel computing version of test. Note that parallel pool must already be open. Default is \code{0}.}
#'   \item{\code{mdl_h0_control}: }{List with restricted model options. See \code{\link{Nmdl}}, \code{\link{ARmdl}}, \code{\link{VARmdl}}, \code{\link{HMmdl}}, \code{\link{MSARmdl}}, or \code{\link{MSVARmdl}} documentation for available and default values.}
#'   \item{\code{mdl_h1_control}: }{List with unrestricted model options. See \code{\link{HMmdl}}, \code{\link{MSARmdl}}, or \code{\link{MSVARmdl}} documentation for available and default values.}
#' }
#'
#' @return List of class \code{LMCLRTest} (\code{S3} object) with attributes including: 
#' \itemize{
#'   \item{\code{mdl_h0}: }{List with restricted model attributes. See \code{\link{Nmdl}}, \code{\link{ARmdl}}, \code{\link{VARmdl}}, \code{\link{HMmdl}}, \code{\link{MSARmdl}}, or \code{\link{MSVARmdl}} documentation for return values.}
#'   \item{\code{mdl_h0}: }{List with unrestricted model attributes. See \code{\link{HMmdl}}, \code{\link{MSARmdl}}, or \code{\link{MSVARmdl}} documentation for return values.}
#'   \item{\code{LRT_0}: }{Value of test statistic from observed data.}
#'   \item{\code{LRN}: }{A (\code{N x 1}) vector of test statistics from data simulated under the null hypothesis.}
#'   \item{\code{pval}: }{P-value of Local Monte Carlo Likelihood Ratio Test.}
#'   \item{\code{LRN_cv}: }{Vector with 90\%, 95\%, and 99\% Monte Carlo critical values (from vector \code{LRN}).}
#'   \item{\code{control}: }{List with test procedure options used.}
#' }
#'
#' @references Rodriguez Rondon, Gabriel and Jean-Marie Dufour. 2022. “Monte Carlo Likelihood Ratio Tests for Markov Switching Models.” \emph{Unpublished manuscript}.
#'
#' @example /inst/examples/LMCLRTest_examples.R
#' @export
LMCLRTest <- function(Y, p, k0, k1, control = list()){
  # ----- Set control values
  con <- list(N = 99,
              burnin = 100,
              converge_check = NULL,
              workers = 0,
              mdl_h0_control = list(),
              mdl_h1_control = list())
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
  mdl_h0 <- estimMdl(Y, p, q, k0, con$mdl_h0_control)
  mdl_h1 <- estimMdl(Y, p, q, k1, con$mdl_h1_control)
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
  if (con$workers>0){
    LRN <- LR_samp_dist_par(mdl_h0, k1, con$N, con$burnin, con$mdl_h0_control, con$mdl_h1_control, con$workers)
  }else{
    LRN <- LR_samp_dist(mdl_h0, k1, con$N, con$burnin, con$mdl_h0_control, con$mdl_h1_control) 
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


#' @title MMC nuisance parameter bounds for univariate models 
#' 
#' @keywords internal
#' 
#' @references Rodriguez Rondon, Gabriel and Jean-Marie Dufour. 2022. “Monte Carlo Likelihood Ratio Tests for Markov Switching Models.” \emph{Unpublished manuscript}.
#' 
#' @export
MMC_bounds <- function(theta_0, mdl_h0, mdl_h1, con){
  k0 <- mdl_h0$k
  k1 <- mdl_h1$k
  # ----- Define lower & upper bounds for search
  theta_low = theta_0 - con$eps
  theta_upp = theta_0 + con$eps
  # create ball around union of eps and 2*standard error (if set to true and SE are finite)
  if ((con$CI_union==TRUE) & all(is.finite(mdl_h0$theta_se)) & all(is.finite(mdl_h1$theta_se))){
    theta_low <- apply(cbind(as.matrix(theta_0 - 2*c(mdl_h0$theta_se,mdl_h1$theta_se)),as.matrix(theta_low)), 1, FUN = min)
    theta_upp <- apply(cbind(as.matrix(theta_0 + 2*c(mdl_h0$theta_se,mdl_h1$theta_se)),as.matrix(theta_upp)), 1, FUN = max)
  }
  # ----- Check that bounds respect admissible regions
  # correct lower bound of variances to be in admissible region
  sigma_ind <- c(mdl_h0$theta_var_ind,mdl_h1$theta_var_ind)
  if (any(theta_low[sigma_ind==1]<=0)==TRUE){
    theta_low[sigma_ind==1][theta_low[sigma_ind==1]<=0] = theta_0[sigma_ind==1][theta_low[sigma_ind==1]<=0]*con$variance_constraint  
  }
  # correct transition probability bounds to be in admissible region
  if (k0==1){
    P_h0_ind <- rep(0,length(mdl_h0$theta))
  }else if (k0>1){
    P_h0_ind <- mdl_h0$theta_P_ind
  }
  P_ind <- c(P_h0_ind, mdl_h1$theta_P_ind)
  theta_low[P_ind==1][theta_low[P_ind==1]<0] <- 0
  theta_upp[P_ind==1][theta_upp[P_ind==1]>1] <- 1
  # ----- output
  mmc_bounds <- list(theta_low = theta_low, theta_upp = theta_upp)
  return(mmc_bounds)
}




#' @title Maximized Monte Carlo Likelihood Ratio Test
#'
#' @description This function performs the Maximized Monte Carlo likelihood ratio 
#' test (MMC-LRT) proposed in Rodriguez Rondon & Dufour (2022).
#' 
#' @param Y  Series to be tested. Must be a (\code{T x q}) matrix.
#' @param p  Number of autoregressive lags. Must be greater than or equal to 0. 
#' @param k0 Number of regimes under null hypothesis. Must be greater than or equal to 1.
#' @param k1 Number of regimes under alternative hypothesis. Must be greater than \code{k0}.
#' @param control List with test procedure options including: 
#' \itemize{
#'   \item{\code{N}: }{Integer determining the number of Monte Carlo simulations. Default is set to \code{99} as in paper.}
#'   \item{\code{burnin}: }{Number of simulated observations to remove from beginning. Default is \code{100}.}
#'   \item{\code{converge_check}: }{String of NULL determining if convergence of model(s) should be verified. Allowed inputs are: "null", "alt", "both", or \code{NULL}. If \code{NULL} (default) no model convergence is verified.}
#'   \item{\code{workers}: }{Integer determining the number of workers to use for parallel computing version of test. Note that parallel pool must already be open. Default is \code{0}.}
#'   \item{\code{type}: }{String that determines the type of optimization algorithm used. Arguments allowed are: \code{"pso"}, \code{"GenSA"}, and \code{"GA"}. Default is \code{"pso"}.}
#'   \item{\code{eps}: }{Double determining the constant value that defines a consistent set for search. Default is \code{0.1}.}
#'   \item{\code{CI_union}: }{Boolean determining if union of set determined by \code{eps} and confidence set should be used to define consistent set for search. Default is \code{TRUE}.}
#'   \item{\code{lambda}: }{Double determining penalty on nonlinear constraint. Default is \code{100}.}
#'   \item{\code{stationary_constraint}: }{Boolean determining if only stationary solutions are considered (if \code{TRUE}) or not (if \code{FALSE}). Default is \code{TRUE}.}
#'   \item{\code{variance_constraint}: }{Double used to determine the lower bound for variance in parameter set for search. Value should be between \code{0} and \code{1} as it is multiplied by consistent point estimates of variances. Default is \code{0.01} (i.e., \code{1\%} of consistent point estimates.}
#'   \item{\code{silence}: }{Boolean determining if optimization steps should be silenced (if \code{TRUE}) or not (if \code{FALSE}). Default is \code{FALSE}.}
#'   \item{\code{threshold_stop}: }{Default is \code{1}.}
#'   \item{\code{mdl_h0_control}: }{List with restricted model options. See \code{\link{Nmdl}}, \code{\link{ARmdl}}, \code{\link{VARmdl}}, \code{\link{HMmdl}}, \code{\link{MSARmdl}}, or \code{\link{MSVARmdl}} documentation for available and default values.}
#'   \item{\code{mdl_h1_control}: }{List with unrestricted model options. See \code{\link{HMmdl}}, \code{\link{MSARmdl}}, or \code{\link{MSVARmdl}} documentation for available and default values.}
#'   \item{\code{type_control}: }{List with optimization algorithm options. See \code{\link[pso]{psoptim}}, \code{\link[GenSA]{GenSA}}, \code{\link[GA]{ga}}. Default is to set \code{list(maxit = 200)} so that maximum number of iterations is \code{200}.}
#' }
#'
#' @return List of class \code{LMCLRTest} (\code{S3} object) with attributes including: 
#' \itemize{
#'   \item{\code{mdl_h0}: }{List with restricted model attributes. See \code{\link{Nmdl}}, \code{\link{ARmdl}}, \code{\link{VARmdl}}, \code{\link{HMmdl}}, \code{\link{MSARmdl}}, or \code{\link{MSVARmdl}} documentation for return values.}
#'   \item{\code{mdl_h0}: }{List with unrestricted model attributes. See \code{\link{HMmdl}}, \code{\link{MSARmdl}}, or \code{\link{MSVARmdl}} documentation for return values.}
#'   \item{\code{LRT_0}: }{Value of test statistic from observed data.}
#'   \item{\code{LRN}: }{A (\code{N x 1}) vector of test statistics from data simulated under the null hypothesis.}
#'   \item{\code{pval}: }{P-value of Local Monte Carlo Likelihood Ratio Test.}
#'   \item{\code{LRN_cv}: }{Vector with 90\%, 95\%, and 99\% Monte Carlo critical values (from vector \code{LRN}).}
#'   \item{\code{control}: }{List with test procedure options used.}
#' }
#'
#' @references Rodriguez Rondon, Gabriel and Jean-Marie Dufour. 2022. “Monte Carlo Likelihood Ratio Tests for Markov Switching Models.” \emph{Unpublished manuscript}.
#'
#' @example /inst/examples/MMCLRTest_examples.R
#' @export
MMCLRTest <- function(Y, p, k0, k1, control = list()){
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
              variance_constraint = 0.01,
              silence = FALSE,
              threshold_stop = 1,
              mdl_h0_control = list(getSE = TRUE),
              mdl_h1_control = list(getSE = TRUE),
              type_control = list(maxit = 200))
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
  mdl_h0 <- estimMdl(Y, p, q, k0, con$mdl_h0_control)
  mdl_h1 <- estimMdl(Y, p, q, k1, con$mdl_h1_control)
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
  theta_0 <- c(mdl_h0$theta, mdl_h1$theta)
  # ----- Define lower & upper bounds for search
  mmc_bounds <- MMC_bounds(theta_0, mdl_h0, mdl_h1, con)
  theta_low <- mmc_bounds$theta_low
  theta_upp <- mmc_bounds$theta_upp
  # ----- Search for Max p-value within bounds
  if (con$type=="pso"){
    # Set PSO specific controls
    con$type_control$trace.stats <- TRUE
    con$type_control$trace <- as.numeric(con$silence==FALSE)
    con$type_control$abstol <- -con$threshold_stop
    # begin optimization
    mmc_out   <- pso::psoptim(par = theta_0, fn = MMCLRpval_fun_min, lower = theta_low, upper = theta_upp, 
                              gr = NULL, control = con$type_control,
                              mdl_h0 = mdl_h0, mdl_h1 = mdl_h1, N = con$N, burnin = con$burnin, workers = con$workers,
                              lambda = con$lambda, stationary_constraint = con$stationary_constraint, 
                              thtol = mdl_h1$control$thtol, mdl_h0_control = con$mdl_h0_control, 
                              mdl_h1_control = con$mdl_h1_control)
    theta     <- mmc_out$par
    pval      <- -mmc_out$value
  }else if(con$type=="GenSA"){
    # Set GenSA specific controls
    con$type_control$trace.mat <- TRUE
    con$type_control$verbose <- con$silence==FALSE
    con$type_control$threshold.stop <- -con$threshold_stop
    # begin optimization
    mmc_out   <- GenSA::GenSA(par = theta_0, fn = MMCLRpval_fun_min, lower = theta_low, upper = theta_upp, 
                              control = con$type_control,
                              mdl_h0 = mdl_h0, mdl_h1 = mdl_h1, N = con$N, burnin = con$burnin, workers = con$workers,
                              lambda = con$lambda, stationary_constraint = con$stationary_constraint, 
                              thtol = mdl_h1$control$thtol, mdl_h0_control = con$mdl_h0_control, 
                              mdl_h1_control = con$mdl_h1_control)
    theta     <- mmc_out$par
    pval      <- -mmc_out$value
  }else if(con$type=="GA"){
    # begin optimization
    mmc_out   <- GA::ga(type = "real-valued", fitness = MMCLRpval_fun, 
                      mdl_h0 = mdl_h0, mdl_h1 = mdl_h1, N = con$N, burnin = con$burnin, workers = con$workers,
                      lambda = con$lambda, stationary_constraint = con$stationary_constraint, 
                      thtol = mdl_h1$control$thtol, mdl_h0_control = con$mdl_h0_control, 
                      mdl_h1_control = con$mdl_h1_control,
                      lower = theta_low, upper = theta_upp, 
                      maxiter = con$type_control$maxit, maxFitness = con$threshold_stop, 
                      monitor = (con$silence==FALSE), suggestions = t(theta_0))
    theta     <- as.matrix(mmc_out@solution[1,])
    pval      <- mmc_out@fitnessValue
  }else if(con$type=="gridSearch"){
    # Grid Search: not ready
  }
  # ----- get test output using optimization output params
  theta_h0 <- theta[1:length(mdl_h0[["theta"]])]
  theta_h1 <- theta[(length(mdl_h0[["theta"]])+1):length(theta)]
  names(theta_h0) <- names(mdl_h0$theta)
  names(theta_h1) <- names(mdl_h1$theta)
  mdl_h0_mmc <- mdledit(mdl_h0, theta_h0, p, q, k0)
  mdl_h0_mmc$logLike <- logLikelihood(mdl_h0_mmc)
  mdl_h0_mmc$AIC <- aic(mdl_h0_mmc$logLike, length(theta_h0))
  mdl_h0_mmc$BIC <- bic(mdl_h0_mmc$logLike, mdl_h0_mmc$n, length(theta_h0))
  mdl_h1_mmc <- mdledit(mdl_h0, theta_h0, p, q, k0)
  mdl_h1_mmc$logLike <- logLikelihood(mdl_h1_mmc)
  mdl_h1_mmc$AIC <- aic(mdl_h1_mmc$logLike, length(theta_h1))
  mdl_h1_mmc$BIC <- bic(mdl_h1_mmc$logLike, mdl_h1_mmc$n, length(theta_h1))
  if (mdl_h0$control$getSE==TRUE){
    mdl_h0_mmc <- thetaSE(mdl_h0_mmc)
  }
  if (mdl_h1$control$getSE==TRUE){
    mdl_h1_mmc <- thetaSE(mdl_h1_mmc)
  }
  # Compute test stats
  LRT_0 = compu_tstat(theta_h0, theta_h1, mdl_h0, mdl_h1, p, q, k0, k1)
  names(LRT_0) <- c("LRT_0")
  # ----- organize test output
  MMCLRTest_output <- list(mdl_h0 = mdl_h0, mdl_h1 = mdl_h1, mdl_h0_mmc = mdl_h0_mmc, mdl_h1_mmc = mdl_h1_mmc, 
                           LRT_0 = LRT_0, 
                           pval = pval,
                           theta_h0 = theta_h0, theta_h1 = theta_h1, control = con)
  class(MMCLRTest_output) <- "MMCLRTest"
  return(MMCLRTest_output)
}



#' @title Bootstrap Likelihood Ratio Test
#' 
#' @description This function performs the bootstrap likelihood ratio 
#' test discussed in Qu & Zhuo (2021) and Kasahara & Shimotsu (2018).
#' 
#' @param Y  Series to be tested. Must be a (\code{T x q}) matrix.
#' @param p  Number of autoregressive lags. Must be greater than or equal to 0. 
#' @param k0 Number of regimes under null hypothesis. Must be greater than or equal to 1.
#' @param k1 Number of regimes under alternative hypothesis. Must be greater than \code{k0}.
#' @param control List with test procedure options including: 
#' \itemize{
#'   \item{\code{B}: }{Integer determining the number of bootstrap simulations. Default is set to \code{999}.}
#'   \item{\code{burnin}: }{Number of simulated observations to remove from beginning. Default is \code{100}.}
#'   \item{\code{converge_check}: }{String of NULL determining if convergence of model(s) should be verified. Allowed inputs are: "null", "alt", "both", or \code{NULL}. If \code{NULL} (default) no model convergence is verified.}
#'   \item{\code{workers}: }{Integer determining the number of workers to use for parallel computing version of test. Note that parallel pool must already be open. Default is \code{0}.}
#'   \item{\code{mdl_h0_control}: }{List with restricted model options. See \code{\link{Nmdl}}, \code{\link{ARmdl}}, \code{\link{VARmdl}}, \code{\link{HMmdl}}, \code{\link{MSARmdl}}, or \code{\link{MSVARmdl}} documentation for available and default values.}
#'   \item{\code{mdl_h1_control}: }{List with unrestricted model options. See \code{\link{HMmdl}}, \code{\link{MSARmdl}}, or \code{\link{MSVARmdl}} documentation for available and default values.}
#' }
#'
#' @return List of class \code{LMCLRTest} (\code{S3} object) with attributes including: 
#' \itemize{
#'   \item{\code{mdl_h0}: }{List with restricted model attributes. See \code{\link{Nmdl}}, \code{\link{ARmdl}}, \code{\link{VARmdl}}, \code{\link{HMmdl}}, \code{\link{MSARmdl}}, or \code{\link{MSVARmdl}} documentation for return values.}
#'   \item{\code{mdl_h0}: }{List with unrestricted model attributes. See \code{\link{HMmdl}}, \code{\link{MSARmdl}}, or \code{\link{MSVARmdl}} documentation for return values.}
#'   \item{\code{LRT_0}: }{Value of test statistic from observed data.}
#'   \item{\code{LRN}: }{A (\code{N x 1}) vector of test statistics from data simulated under the null hypothesis.}
#'   \item{\code{pval}: }{P-value of Local Monte Carlo Likelihood Ratio Test.}
#'   \item{\code{LRN_cv}: }{Vector with 90\%, 95\%, and 99\% Monte Carlo critical values (from vector \code{LRN}).}
#'   \item{\code{control}: }{List with test procedure options used.}
#' }
#' 
#' @references Qu, Zhongjun, and Fan Zhuo. 2021. “Likelihood Ratio-Based Tests for Markov Regime Switching.” \emph{The Review of Economic Studies} 88 (2): 937–968.
#' @references Kasahara, Hiroyuk, and Katsum Shimotsu. 2018. “Testing the number of regimes in Markov regime switching models.” \emph{arXiv preprint arXiv:1801.06862}.
#' 
#' @example /inst/examples/BootLRTest_examples.R
#' @export
BootLRTest <- function(Y, p, k0, k1, control = list()){
  # ----- Set control values
  con <- list(B = 999,
              burnin = 100,
              converge_check = NULL,
              workers = 0,
              mdl_h0_control = list(),
              mdl_h1_control = list())
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
  mdl_h0 <- estimMdl(Y, p, q, k0, con$mdl_h0_control)
  mdl_h1 <- estimMdl(Y, p, q, k1, con$mdl_h1_control)
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
  if (con$workers>0){
    LRN <- LR_samp_dist_par(mdl_h0, k1, con$B, con$burnin, con$mdl_h0_control, con$mdl_h1_control, con$workers)
  }else{
    LRN <- LR_samp_dist(mdl_h0, k1, con$B, con$burnin, con$mdl_h0_control, con$mdl_h1_control) 
  }
  # ----- Compute p-value
  pval <- sum(LRN>LRT_0)/con$B # [eq. 4.62] (Davidson & MacKinnon, 2004)
  # ----- get critical values
  LRN     <- as.matrix(sort(LRN))
  LRN_cv  <- LRN[round(c(0.90,0.95,0.99)*nrow(LRN)),]
  names(LRN_cv)  <- paste0(c("0.90","0.95","0.99"), "%")
  # ----- Organize output
  BootLRTest_output <- list(mdl_h0 = mdl_h0, mdl_h1 = mdl_h1, LRT_0 = LRT_0, LRN = LRN,
                            pval = pval, LRN_cv = LRN_cv, control = con)
  class(BootLRTest_output) <- "BootLRTest"
  return(BootLRTest_output)
}


