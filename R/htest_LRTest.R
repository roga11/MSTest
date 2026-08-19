# Package-level environment for cluster caching (avoids repeated cluster creation in MMC)
.MSTest_cluster_env <- new.env(parent = emptyenv())

#' @title Monte Carlo Likelihood Ratio Test sample distribution (parallel version)
#'
#' @description This function simulates the sample distribution under the null hypothesis using parallel workers.
#'
#' @param mdl_h0 List with restricted model properties.
#' @param k1 integer specifying the number of regimes under the alternative hypothesis.
#' @param N integer specifying the number of replications.
#' @param burnin integer specifying the number of observations to drop from beginning of simulation.
#' @param Z Exogenous regressors matrix or \code{NULL}.
#' @param mdl_h0_control List with controls/options used to estimate restricted model.
#' @param mdl_h1_control List with controls/options used to estimate unrestricted model.
#' @param workers Integer determining the number of parallel workers. If \code{workers > N}, the effective count is silently capped at \code{N} (more workers than simulations provides no speedup and can thin the pre-drawn buffer below its retry budget).
#' @param seed Optional integer seed for reproducible parallel RNG (L'Ecuyer-CMRG). Default is \code{NULL} (random but independent streams).
#' @param predrawn_eps Optional list of pre-drawn standard normal matrices for fixed-error simulation (Dufour 2006, Prop. 4.2). Each element is a \code{(T+burnin) x q} matrix. Default is \code{NULL} (draw fresh innovations).
#' @param predrawn_state_rand Optional matrix of pre-drawn U[0,1] for state transitions. Columns correspond to simulation draws. Default is \code{NULL} (only needed when null model has multiple regimes).
#'
#' @return vector of simulated LRT statistics
#'
#' @keywords internal
#'
#' @references Rodriguez-Rondon, G., & Dufour, J.-M. 2026a. "Monte Carlo Likelihood-Ratio Tests for Markov Switching Models." \emph{Bank of Canada Staff Working Paper}, No. 2026-23. doi: 10.34989/swp-2026-23.
#'
#' @export
LR_samp_dist_par <- function(mdl_h0, k1, N, burnin, Z, mdl_h0_control, mdl_h1_control,
                              workers, seed = NULL,
                              predrawn_eps = NULL, predrawn_state_rand = NULL){
  # ----- Cap workers at N (no benefit to more workers than simulations;
  # also prevents per-worker pre-drawn buffer slice from falling below retry budget)
  if (workers > N) workers <- N
  # ----- Get or create cluster
  cl <- .MSTest_cluster_env$cl
  own_cluster <- FALSE
  if (is.null(cl) || length(cl) != workers) {
    cl <- parallel::makePSOCKcluster(workers)
    parallel::clusterEvalQ(cl, library(MSTest))
    own_cluster <- TRUE
  }
  if (own_cluster) on.exit(parallel::stopCluster(cl), add = TRUE)
  # ----- Set reproducible parallel RNG (L'Ecuyer-CMRG)
  # When pre-drawn innovations are provided AND seed is NULL (MMC via C++),
  # skip RNG reset — workers keep their deterministic state from cluster setup
  # (seeded once in MMCLRTest). This ensures EM starting values within workers
  # are reproducible across optimizer iterations.
  # When seed is provided (LMCLRTest) or no predrawn innovations, always seed.
  if (is.null(predrawn_eps) || !is.null(seed)) {
    if (!is.null(seed)) {
      parallel::clusterSetRNGStream(cl, iseed = seed)
    } else {
      parallel::clusterSetRNGStream(cl)
    }
  }
  if (!is.null(predrawn_eps)) {
    # ----- Pre-drawn mode: slice innovations per worker (Dufour 2006, Prop. 4.2)
    N_buffer <- length(predrawn_eps)
    N_per <- diff(round(seq(0, N, length.out = workers + 1)))
    N_buf_per <- diff(round(seq(0, N_buffer, length.out = workers + 1)))
    buf_ends <- cumsum(N_buf_per)
    buf_starts <- c(0, buf_ends[-length(buf_ends)])
    # Per-worker slices (only send what each worker needs)
    eps_slices <- lapply(seq_len(workers), function(w) {
      predrawn_eps[(buf_starts[w] + 1):buf_ends[w]]
    })
    sr_slices <- if (!is.null(predrawn_state_rand)) {
      lapply(seq_len(workers), function(w) {
        predrawn_state_rand[, (buf_starts[w] + 1):buf_ends[w], drop = FALSE]
      })
    } else {
      replicate(workers, NULL, simplify = FALSE)
    }
    LRN_all <- parallel::clusterMap(cl,
      function(n_i, eps_slice, sr_slice,
               mdl_h0, k1, burnin, Z, mdl_h0_control, mdl_h1_control) {
        LR_samp_dist(mdl_h0, k1, n_i, burnin, Z, mdl_h0_control, mdl_h1_control,
                     eps_slice, sr_slice)
      },
      N_per, eps_slices, sr_slices,
      MoreArgs = list(mdl_h0 = mdl_h0, k1 = k1, burnin = burnin, Z = Z,
                      mdl_h0_control = mdl_h0_control, mdl_h1_control = mdl_h1_control))
  } else {
    # ----- Fresh-draw mode: existing parLapply with L'Ecuyer-CMRG streams
    N_per <- diff(round(seq(0, N, length.out = workers + 1)))
    LRN_all <- parallel::parLapply(cl, N_per,
      function(n_i, mdl_h0, k1, burnin, Z, mdl_h0_control, mdl_h1_control) {
        LR_samp_dist(mdl_h0, k1, n_i, burnin, Z, mdl_h0_control, mdl_h1_control)
      }, mdl_h0, k1, burnin, Z, mdl_h0_control, mdl_h1_control)
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
    # EM fits report convergence per the criterion chosen by control$conv; MLE fits
    # do not set the flag, so fall back to the parameter-change test for them.
    if (is.null(mdl$converged)) mdl$converged = (mdl$deltath <= mdl$control$thtol)
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
      # EM fits report convergence per the criterion chosen by control$conv; MLE fits
    # do not set the flag, so fall back to the parameter-change test for them.
    if (is.null(mdl$converged)) mdl$converged = (mdl$deltath <= mdl$control$thtol)  
    }else{
      mdl <- MSARXmdl(Y, p, k, Z, control)
      # EM fits report convergence per the criterion chosen by control$conv; MLE fits
    # do not set the flag, so fall back to the parameter-change test for them.
    if (is.null(mdl$converged)) mdl$converged = (mdl$deltath <= mdl$control$thtol)
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
      # EM fits report convergence per the criterion chosen by control$conv; MLE fits
    # do not set the flag, so fall back to the parameter-change test for them.
    if (is.null(mdl$converged)) mdl$converged = (mdl$deltath <= mdl$control$thtol)
    }else{
      mdl <- MSVARXmdl(Y, p, k, Z, control)
      # EM fits report convergence per the criterion chosen by control$conv; MLE fits
    # do not set the flag, so fall back to the parameter-change test for them.
    if (is.null(mdl$converged)) mdl$converged = (mdl$deltath <= mdl$control$thtol)
    }
  }
  return(mdl)
}


#' @title Monte Carlo Likelihood Ratio Test
#' 
#' @description This function performs the Local Monte Carlo likelihood ratio 
#' test (LMC-LRT) proposed in Rodriguez-Rondon & Dufour (2026a). As discussed in 
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
#'   \item workers: Integer determining the number of workers to use for parallel computing. Default is \code{0} (sequential). If \code{workers > N}, the effective count is capped at \code{N} and an informational message is printed.
#'   \item mc_seed: Integer seed for reproducible Monte Carlo simulations. When set, seeds the RNG before observed-data estimation and pre-draws innovations once, ensuring fully reproducible results (including the observed LRT statistic). When \code{workers > 0}, also sets L'Ecuyer-CMRG parallel RNG streams. Default is \code{NULL}, in which case all internal randomness is drawn from the ambient RNG state, so a script-level \code{set.seed()} placed before the call makes the procedure reproducible as well (sequential and parallel).
#'   \item init_method: String determining how the alternative (\code{k1}-regime) model is initialized. \code{"warmstart"} (the default) adds deterministic warm starts built by embedding the null-model fit into the \code{k1}-regime space (see \code{\link{warmstart_theta}}) to the \code{use_diff_init} random starts, applied identically to the observed data and to every simulated draw; combined with the LR >= 0 policy this makes negative LR statistics essentially impossible. \code{"random"} reproduces the legacy behavior (random starts only).
#'   \item mdl_h0_control: List with restricted model options. See \code{\link{Nmdl}}, \code{\link{ARmdl}}, \code{\link{VARmdl}}, \code{\link{HMmdl}}, \code{\link{MSARmdl}}, or \code{\link{MSVARmdl}} documentation for available and default values.
#'   \item mdl_h1_control: List with unrestricted model options. See \code{\link{HMmdl}}, \code{\link{MSARmdl}}, or \code{\link{MSVARmdl}} documentation for available and default values.
#'   \item use_diff_init_sim: Value which determines the number of initial values to use when estimating models for null distribution. Default is set to use the same as specified in \code{mdl_h0_control} and \code{mdl_h1_control}.
#' }
#'
#' @section Numerical policy for the LR statistic: For nested models the LR statistic is non-negative by construction (the null fit embeds into the alternative parameter space with equal likelihood), so a negative computed value can only be a numerical optimization artifact. Both the observed statistic and every simulated statistic are therefore floored at 0 (identically on both sides, which preserves the exchangeability underlying the Monte Carlo p-value); a warning reports the magnitude when it is non-negligible. For \code{k0 > 1}, the null log-likelihood is additionally floored at the one-regime (linear) fit on both sides. Ties at 0 are handled by the randomized tie-breaking rule of Dufour (2006) already implemented in \code{\link{MCpval}}; ties at the bottom of the support can only make the test conservative. The regime variances are bounded away from zero during estimation (control \code{em_variance_constraint} of the model constructors, default \code{0.01} of the one-regime fit variance), because the Markov-switching likelihood is otherwise unbounded and degenerate solutions produce spurious likelihood-ratio spikes; the statistic is therefore a constrained-likelihood ratio, computed by the same rule for the observed and the simulated samples, which is all the Monte Carlo argument requires. Transition probabilities are not constrained.
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
#' @references Rodriguez-Rondon, G., & Dufour, J.-M. 2026a. "Monte Carlo Likelihood-Ratio Tests for Markov Switching Models." \emph{Bank of Canada Staff Working Paper}, No. 2026-23. doi: 10.34989/swp-2026-23.
#' @example /inst/examples/LMCLRTest_examples.R
#' @export
LMCLRTest <- function(Y, p, k0, k1, Z = NULL, control = list()){
  # ----- Set control values
  con <- list(N = 99,
              burnin = 100,
              converge_check = NULL,
              workers = 0,
              mc_seed = NULL,
              init_method = "warmstart",
              mdl_h0_control = list(),
              mdl_h1_control = list(),
              use_diff_init_sim = NULL)
  # ----- Perform some checks for controls
  nmsC <- names(con)
  con[(namc <- names(control))] <- control
  if (length(noNms <- namc[!namc %in% nmsC])){
    warning("unknown names in control: ", paste(noNms,collapse=", "))
  }
  con$init_method <- match.arg(con$init_method, c("warmstart", "random"))
  # ----- Perform other checks
  if (is.matrix(Y)){
    q <- ncol(Y)
  }else{
    stop("Observations Y must be a (T x q) matrix.")
  }
  # ----- Seed RNG for reproducible estimation (random initial values in EM)
  if (!is.null(con$mc_seed)) set.seed(con$mc_seed)
  # ----- Estimate models using observed data
  mdl_h0 <- estimMdl(Y, p, q, k0, Z, con$mdl_h0_control)
  if ((con$init_method == "warmstart") && (k1 > 1)){
    # Warm-start the alternative model from the null fit (regime-duplication
    # embedding, see warmstart_theta); deterministic, so the RNG stream and the
    # random starts are unaffected. The same rule is applied per draw when
    # simulating the null distribution (see LR_samp_dist).
    con$mdl_h1_control$init_theta_extra <- warmstart_theta(mdl_h0, k1,
                                                           msmu = con$mdl_h1_control$msmu,
                                                           msvar = con$mdl_h1_control$msvar)
  }
  mdl_h1 <- estimMdl(Y, p, q, k1, Z, con$mdl_h1_control)
  con$mdl_h0_control <- mdl_h0$control
  con$mdl_h1_control <- mdl_h1$control
  # Strip the observed-data warm starts from the stored control: they are specific
  # to the observed null fit; the simulated side rebuilds them per draw.
  con$mdl_h1_control$init_theta_extra <- NULL
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
  if (k0 > 1){
    # H0 linear floor: the one-regime model is a feasible point of the null space,
    # so its (closed-form) log-likelihood bounds logL0 from below. Applied
    # identically per draw when simulating the null (see LR_samp_dist).
    mdl_lin <- estimMdl(Y, p, q, 1, Z, list(getSE = FALSE))
    if (is.finite(mdl_lin$logLike) && (mdl_lin$logLike > logL0)) logL0 <- mdl_lin$logLike
  }
  LRT_0 <- -2*(logL0-logL1)
  # ----- Perform check of test stat and model parameters
  if ((is.finite(LRT_0)==FALSE) | (any(is.finite(theta_h0)==FALSE)) | (any(is.finite(theta_h1)==FALSE))){
    stop("LRT_0 or model parameters are not finite. Run again to use different initial values")
  }
  if (LRT_0 < 0){
    # LR >= 0 by construction for nested models (the null fit embeds into the
    # alternative space with equal likelihood); a negative value is a numerical
    # optimization artifact. Cap at 0, identically to the simulated statistics.
    if (LRT_0 < -1e-6){
      warning(sprintf("LMCLRTest: the observed LR statistic was negative (%.3g; numerical optimization artifact) and was set to 0. Consider a higher 'use_diff_init'.", LRT_0))
    }
    LRT_0 <- 0
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
    # k0 = 1 nulls are linear models fitted in closed form and do not take
    # 'use_diff_init'; passing it there would warn once per simulated draw.
    if (k0>1) mdl_h0_null_cont$use_diff_init <- con$use_diff_init_sim
    mdl_h1_null_cont$use_diff_init <- con$use_diff_init_sim
  }
  # Tell LR_samp_dist to warm-start each draw's alternative fit from that draw's
  # own null fit (identical procedure to the observed-data side).
  mdl_h1_null_cont$init_method <- con$init_method
  # ----- Pre-draw innovations for reproducibility when mc_seed is set
  predrawn_eps <- NULL
  predrawn_state_rand <- NULL
  if (!is.null(con$mc_seed)) {
    set.seed(con$mc_seed)
    N_buffer <- ceiling(con$N * 1.5) + 10
    Teps <- mdl_h0$n + con$burnin
    q_dim <- mdl_h0$q
    predrawn_eps <- lapply(seq_len(N_buffer), function(i) matrix(rnorm(Teps * q_dim), Teps, q_dim))
    if (k0 > 1) {
      predrawn_state_rand <- matrix(runif(Teps * N_buffer), Teps, N_buffer)
    }
  }
  # ----- Cap workers at N (no benefit to more workers than simulations)
  if (con$workers > 0 && con$workers > con$N) {
    message("LMCLRTest: capping workers at N = ", con$N,
            " (", con$workers, " requested).")
    con$workers <- con$N
  }
  if (con$workers>0){
    LRN <- LR_samp_dist_par(mdl_h0, k1, con$N, con$burnin, Zsim, mdl_h0_null_cont, mdl_h1_null_cont,
                              con$workers, seed = con$mc_seed,
                              predrawn_eps = predrawn_eps, predrawn_state_rand = predrawn_state_rand)
  }else{
    LRN <- LR_samp_dist(mdl_h0, k1, con$N, con$burnin, Zsim, mdl_h0_null_cont, mdl_h1_null_cont,
                          predrawn_eps, predrawn_state_rand)
  }
  # ----- Drop non-finite (failed) draws and check the null distribution is non-empty
  LRN <- LRN[is.finite(LRN)]
  if (length(LRN)==0){
    stop("LMCLRTest: the null distribution could not be simulated (all draws failed even after the re-draw safety); the p-value cannot be computed. Increase 'use_diff_init' or inspect the observed model fit.")
  }
  # ----- Compute p-value (before sorting, so the randomized tie-break in MCpval
  # pairs uniforms with draws in simulation order)
  pval <- MCpval(LRT_0, LRN, "geq")
  # ----- get critical values
  LRN     <- as.matrix(sort(LRN))
  LRN_cv  <- LRN[round(c(0.90,0.95,0.99)*nrow(LRN)),]
  names(LRN_cv)  <- paste0(c("0.90","0.95","0.99"), "%")
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
#' @references Rodriguez-Rondon, G., & Dufour, J.-M. 2026a. "Monte Carlo Likelihood-Ratio Tests for Markov Switching Models." \emph{Bank of Canada Staff Working Paper}, No. 2026-23. doi: 10.34989/swp-2026-23.
#' 
#' @export
MMC_bounds <- function(mdl_h0, con){
  theta_0   <- mdl_h0$theta
  k0        <- mdl_h0$k
  # ----- Define lower & upper bounds for search
  theta_low = theta_0 - con$eps
  theta_upp = theta_0 + con$eps
  # create ball around the union of eps and 2*standard error, per parameter; a parameter
  # without a finite standard error (e.g. a variance at the EM floor) contributes only its
  # eps band, without shrinking the other parameters' bounds
  if (isTRUE(con$CI_union) && !is.null(mdl_h0$theta_se) &&
      length(c(mdl_h0$theta_se)) == length(c(theta_0))){
    se <- c(mdl_h0$theta_se)
    se[!is.finite(se)] <- 0
    theta_low <- pmin(c(theta_low), c(theta_0) - 2*se)
    theta_upp <- pmax(c(theta_upp), c(theta_0) + 2*se)
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
#' test (MMC-LRT) proposed in Rodriguez-Rondon & Dufour (2026a).
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
#'   \item workers: Integer determining the number of workers to use for parallel computing. Default is \code{0} (sequential). If \code{workers > N}, the effective count is capped at \code{N} and an informational message is printed (unless \code{silence = TRUE}).
#'   \item mc_seed: Integer seed for reproducible Monte Carlo simulations and fixed-error MMC (Dufour 2006, Prop. 4.2). When set, seeds the RNG before observed-data estimation, pre-draws innovations once and holds them fixed across all optimizer evaluations, ensuring full reproducibility (including the observed LRT statistic) and theoretical size control. Internal seeds (optimizer trajectory, worker streams, CRN scheme) are drawn from the \code{mc_seed} stream rather than derived by arithmetic, so consecutive \code{mc_seed} values across replications of an outer simulation study do not produce colliding RNG streams; for studies with more than about 10^4 replications, well-separated seeds are still recommended. Default is \code{NULL}, in which case internal seeds are drawn from the ambient RNG state, so a script-level \code{set.seed()} placed before the call makes the procedure reproducible as well (sequential and parallel). Note: GenSA is deterministic given the fixed starting value used here, so its search trajectory does not depend on any seed; pso and GA use the optimizer seed.
#'   \item init_method: String determining how the alternative (\code{k1}-regime) model is initialized. \code{"warmstart"} (the default) adds deterministic warm starts built by embedding the null-model fit into the \code{k1}-regime space (see \code{\link{warmstart_theta}}) to the \code{use_diff_init} random starts, applied identically to the observed data and to every simulated draw. \code{"random"} reproduces the legacy behavior.
#'   \item crn: Boolean determining whether common random numbers are used across the nuisance-parameter search (default \code{TRUE}): an identical estimation-RNG stream (EM starting values; master and worker streams) is replayed at every candidate theta evaluation, so that together with the pre-drawn innovations the MC p-value is a deterministic function of theta (Dufour 2006, Prop. 4.2, with the estimation randomness folded into the fixed disturbance vector). This removes optimizer objective noise and guarantees \code{pval >= pval_0} exactly. Set to \code{FALSE} for the legacy behavior (estimation randomness evolves across evaluations; still valid, but the objective is noisy and the reported maximum tends to be conservative). Note: with \code{init_method = "random"}, \code{crn = FALSE}, \code{workers = 0} and no \code{mc_seed}, a script-seeded call reproduces version 0.1.9 bit-for-bit.
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
#'   \item threshold_stop: Double determining the p-value level at which the nuisance-parameter search stops early (the search seeks the maximum p-value, so once a candidate reaches this level the test decision at any smaller significance level is settled). If the p-value at the initial values already satisfies \code{pval_0 >= threshold_stop}, the optimizer is skipped entirely and \code{pval = pval_0} is returned (with \code{mmc_optimout} set to an early-stop marker). Default is \code{1} (no early stopping).
#'   \item mdl_h0_control: List with restricted model options. See \code{\link{Nmdl}}, \code{\link{ARmdl}}, \code{\link{VARmdl}}, \code{\link{HMmdl}}, \code{\link{MSARmdl}}, or \code{\link{MSVARmdl}} documentation for available and default values.
#'   \item mdl_h1_control: List with unrestricted model options. See \code{\link{HMmdl}}, \code{\link{MSARmdl}}, or \code{\link{MSVARmdl}} documentation for available and default values.
#'   \item use_diff_init_sim: Value which determines the number of initial values to use when estimating models for null distribution. Default is set to use the same as specified in \code{mdl_h0_control} and \code{mdl_h1_control}.
#'   \item maxit: Int determining the budget of the nuisance-parameter search, as a cap on objective-function evaluations: passed as \code{maxf} to \code{pso} and as \code{max.call} to \code{GenSA} (both are evaluation counts, not iterations), and as \code{maxiter} to \code{GA} (generations, so roughly \code{popSize} times as many evaluations). Note that \code{GenSA} treats its cap as a soft bound and can substantially exceed it for higher-dimensional searches. Default is \code{50}.
#'   \item optim_control: List with optimization algorithm options. See \code{\link[pso]{psoptim}}, \code{\link[GenSA]{GenSA}}, \code{\link[GA]{ga}}. Default is \code{list()} (an empty list); the search budget is instead controlled by the separate \code{maxit} control above.
#' }
#'
#' @return List of class \code{MMCLRTest} (\code{S3} object) with attributes including:
#' \itemize{
#'   \item mdl_h0: List with restricted model attributes.
#'   \item mdl_h1: List with unrestricted model attributes.
#'   \item mdl_h0_mmc: List with restricted model attributes at the maximizing nuisance parameter values.
#'   \item mdl_h1_mmc: List with unrestricted model attributes (same as \code{mdl_h1}).
#'   \item LRT_0: Value of test statistic from observed data.
#'   \item pval: Maximized p-value of the Maximized Monte Carlo Likelihood Ratio Test.
#'   \item theta_h0: Maximizing nuisance parameter vector for the restricted model.
#'   \item theta_h1: Parameter vector of the unrestricted model.
#'   \item control: List with test procedure options used.
#'   \item mmc_optimout: Optimization output object returned by the selected optimizer (\code{pso}, \code{GenSA}, or \code{GA}).
#'   \item pval_0: Monte Carlo p-value at the initial (estimated) nuisance parameter values \code{theta_0} (the Local MC point of the search). Under \code{crn = TRUE} the reported \code{pval} satisfies \code{pval >= pval_0} by construction.
#' }
#'
#' @section Numerical policy for the LR statistic: For nested models the LR statistic is non-negative by construction (the null fit embeds into the alternative parameter space with equal likelihood), so a negative computed value can only be a numerical optimization artifact. Both the observed statistic and every simulated statistic are therefore floored at 0 (identically on both sides, which preserves the exchangeability underlying the Monte Carlo p-value); a warning reports the magnitude when it is non-negligible. For \code{k0 > 1}, the null log-likelihood is additionally floored at the one-regime (linear) fit on both sides. Ties at 0 are handled by the randomized tie-breaking rule of Dufour (2006) already implemented in \code{\link{MCpval}}; ties at the bottom of the support can only make the test conservative. The regime variances are bounded away from zero during estimation (control \code{em_variance_constraint} of the model constructors, default \code{0.01} of the one-regime fit variance), because the Markov-switching likelihood is otherwise unbounded and degenerate solutions produce spurious likelihood-ratio spikes; the statistic is therefore a constrained-likelihood ratio, computed by the same rule for the observed and the simulated samples, which is all the Monte Carlo argument requires. Transition probabilities are not constrained.
#'
#' @references Rodriguez-Rondon, G., & Dufour, J.-M. 2026a. "Monte Carlo Likelihood-Ratio Tests for Markov Switching Models." \emph{Bank of Canada Staff Working Paper}, No. 2026-23. doi: 10.34989/swp-2026-23.
#' @example /inst/examples/MMCLRTest_examples.R
#' @export
MMCLRTest <- function(Y, p, k0, k1, Z = NULL, control = list()){
  # ----- Set control values
  con <- list(N = 99,
              burnin = 100,
              converge_check = NULL,
              workers = 0,
              mc_seed = NULL,
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
              init_method = "warmstart",
              crn = TRUE,
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
  con$init_method <- match.arg(con$init_method, c("warmstart", "random"))
  # Merge (rather than replace) the model-control sublists with their defaults so a
  # user-supplied sublist that omits e.g. 'getSE' keeps the default instead of
  # dropping it (a missing getSE previously crashed the CI_union check below).
  con$mdl_h0_control <- utils::modifyList(list(getSE = TRUE), con$mdl_h0_control)
  con$mdl_h1_control <- utils::modifyList(list(getSE = TRUE), con$mdl_h1_control)
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
  # ----- Seed RNG for reproducible estimation (random initial values in EM)
  if (!is.null(con$mc_seed)) set.seed(con$mc_seed)
  # ----- Estimate models using observed data
  mdl_h0 <- estimMdl(Y, p, q, k0, Z, con$mdl_h0_control)
  if ((con$init_method == "warmstart") && (k1 > 1)){
    # Warm-start the alternative model from the null fit (regime-duplication
    # embedding, see warmstart_theta); deterministic, so the RNG stream and the
    # random starts are unaffected. The same rule is applied per draw when
    # simulating the null distribution (see LR_samp_dist).
    con$mdl_h1_control$init_theta_extra <- warmstart_theta(mdl_h0, k1,
                                                           msmu = con$mdl_h1_control$msmu,
                                                           msvar = con$mdl_h1_control$msvar)
  }
  mdl_h1 <- estimMdl(Y, p, q, k1, Z, con$mdl_h1_control)
  con$mdl_h0_control <- mdl_h0$control
  con$mdl_h1_control <- mdl_h1$control
  # Strip the observed-data warm starts from the stored control: they are specific
  # to the observed null fit; the simulated side rebuilds them per draw.
  con$mdl_h1_control$init_theta_extra <- NULL
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
  # ----- Compute observed test statistic (LRT_0), held FIXED across the optimizer
  #       search per Dufour (2006, eq. 4.22). Mirrors LMCLRTest; only the simulated
  #       null distribution varies with the candidate nuisance parameters.
  logL0 <- mdl_h0$logLike
  logL1 <- mdl_h1$logLike
  if (k0 > 1){
    # H0 linear floor: the one-regime model is a feasible point of the null space,
    # so its (closed-form) log-likelihood bounds logL0 from below. Applied
    # identically per draw when simulating the null (see LR_samp_dist).
    mdl_lin <- estimMdl(Y, p, q, 1, Z, list(getSE = FALSE))
    if (is.finite(mdl_lin$logLike) && (mdl_lin$logLike > logL0)) logL0 <- mdl_lin$logLike
  }
  LRT_0 <- -2*(logL0-logL1)
  if ((is.finite(LRT_0)==FALSE) | (any(is.finite(mdl_h0$theta)==FALSE)) | (any(is.finite(mdl_h1$theta)==FALSE))){
    stop("LRT_0 or model parameters are not finite. Run again to use different initial values")
  }
  if (LRT_0 < 0){
    # LR >= 0 by construction for nested models (the null fit embeds into the
    # alternative space with equal likelihood); a negative value is a numerical
    # optimization artifact. Cap at 0, identically to the simulated statistics.
    if (LRT_0 < -1e-6){
      warning(sprintf("MMCLRTest: the observed LR statistic was negative (%.3g; numerical optimization artifact) and was set to 0. Consider a higher 'use_diff_init'.", LRT_0))
    }
    LRT_0 <- 0
  }
  names(LRT_0) <- c("LRT_0")
  # ----- Define lower & upper bounds for search
  mmc_bounds <- MMC_bounds(mdl_h0, con)
  theta_low <- mmc_bounds$theta_low
  theta_upp <- mmc_bounds$theta_upp
  # ----- Search for Max p-value within bounds
  mdl_h0_null_cont <- con$mdl_h0_control
  mdl_h1_null_cont <- con$mdl_h1_control
  if (is.null(con$use_diff_init_sim)==FALSE){
    # k0 = 1 nulls are linear models fitted in closed form and do not take
    # 'use_diff_init'; passing it there would warn once per simulated draw.
    if (k0>1) mdl_h0_null_cont$use_diff_init <- con$use_diff_init_sim
    mdl_h1_null_cont$use_diff_init <- con$use_diff_init_sim
  }
  # Tell LR_samp_dist to warm-start each draw's alternative fit from that draw's
  # own null fit (identical procedure to the observed-data side).
  mdl_h1_null_cont$init_method <- con$init_method
  if (is.null(Z)==FALSE){
    Zsim <- Z[(p+1):nrow(Z),,drop=F]
    exog <- TRUE
  }else{
    Zsim <- Z
    exog <- FALSE
  }
  # ----- Internal seeds (optimizer trajectory, worker streams, CRN scheme).
  # Drawn from the mc_seed stream and then re-seeding (draw-then-reseed) so the
  # pre-drawn innovations below are bit-identical for a given mc_seed. Drawn
  # seeds -- rather than mc_seed+k arithmetic -- avoid systematic RNG-stream
  # collisions when consecutive parent seeds are used across replications of an
  # outer Monte Carlo study. When mc_seed is NULL the seeds are drawn from the
  # ambient RNG state, so a script-level set.seed() upstream still makes the
  # whole procedure reproducible; in the sequential legacy configuration
  # (mc_seed = NULL, crn = FALSE, workers = 0) the draw is DEFERRED until after
  # the p_seed evaluation so that a script-seeded run reproduces version 0.1.9
  # bit-for-bit through the pre-drawn innovations and p_seed (only the optimizer
  # seed differs, which GenSA ignores).
  need_early_seeds <- (!is.null(con$mc_seed)) || isTRUE(con$crn) || (con$workers > 0)
  if (need_early_seeds){
    if (!is.null(con$mc_seed)) set.seed(con$mc_seed)
    internal_seeds <- sample.int(.Machine$integer.max, 3L)
  }
  # ----- Pre-draw innovations for fixed-error MMC (Dufour 2006, Prop. 4.2)
  N_buffer <- ceiling(con$N * 1.5) + 10
  Teps <- mdl_h0$n + con$burnin   # T + burnin
  q_dim <- mdl_h0$q
  # Set RNG seed for reproducibility (consumed by rnorm/runif below)
  if (!is.null(con$mc_seed)) set.seed(con$mc_seed)
  # Pre-draw Gaussian innovations: list of N_buffer matrices, each (T+burnin) x q
  predrawn_eps <- lapply(seq_len(N_buffer), function(i) matrix(rnorm(Teps * q_dim), Teps, q_dim))
  # Pre-draw state transition uniforms (only needed if null has multiple regimes)
  predrawn_state_rand <- NULL
  if (k0 > 1) {
    predrawn_state_rand <- matrix(runif(Teps * N_buffer), Teps, N_buffer)
  }
  # ----- Cap workers at N (no benefit to more workers than simulations)
  if (con$workers > 0 && con$workers > con$N) {
    if (!isTRUE(con$silence))
      message("MMCLRTest: capping workers at N = ", con$N,
              " (", con$workers, " requested).")
    con$workers <- con$N
  }
  # ----- Create and cache parallel cluster for MMC (reused across optimizer
  # iterations and, under crn = TRUE, re-seeded identically at each objective
  # evaluation). Created before the p_seed evaluation below so that p_seed uses
  # the same worker configuration as the optimizer's evaluations.
  if (con$workers > 0) {
    .MSTest_cluster_env$cl <- parallel::makePSOCKcluster(con$workers)
    parallel::clusterEvalQ(.MSTest_cluster_env$cl, library(MSTest))
    # Initial worker RNG streams (L'Ecuyer-CMRG). With crn = TRUE these are
    # re-seeded at every objective evaluation; this initialization governs the
    # draws only when crn = FALSE (legacy behavior: streams evolve across
    # evaluations since LR_samp_dist_par skips the reset when predrawn).
    parallel::clusterSetRNGStream(.MSTest_cluster_env$cl, iseed = internal_seeds[2L])
    on.exit({
      parallel::stopCluster(.MSTest_cluster_env$cl)
      .MSTest_cluster_env$cl <- NULL
    }, add = TRUE)
  }
  # ----- CRN (common random numbers) objective wrapper: replays an identical
  # estimation-RNG stream (master and, when parallel, worker streams) at every
  # candidate theta evaluation so that -- together with the pre-drawn
  # innovations -- the MC p-value is a deterministic function of theta during
  # the nuisance-parameter search (Dufour 2006, Prop. 4.2, with the EM
  # starting-value randomness folded into the fixed disturbance vector u). The
  # master RNG state is saved and restored around each evaluation (via on.exit,
  # so errors cannot corrupt it), leaving the optimizer's own search randomness
  # undisturbed.
  crn_seed <- if (isTRUE(con$crn)) internal_seeds[3L] else NA_integer_
  crn_wrap <- function(objfn){
    if (!isTRUE(con$crn)) return(objfn)
    function(...){
      if (!exists(".Random.seed", envir = globalenv(), inherits = FALSE)) stats::runif(1)
      old_seed <- get(".Random.seed", envir = globalenv())
      on.exit(assign(".Random.seed", old_seed, envir = globalenv()), add = TRUE)
      set.seed(crn_seed)
      cl <- .MSTest_cluster_env$cl   # looked up at call time (may be NULL)
      if (!is.null(cl)) parallel::clusterSetRNGStream(cl, iseed = crn_seed)
      objfn(...)
    }
  }
  obj_min <- crn_wrap(MMCLRpval_fun_min)   # pso / GenSA (minimize the negative p-value)
  obj_max <- crn_wrap(MMCLRpval_fun)       # GA (maximize the p-value) and p_seed below
  # ----- Verify the seed theta_0 (= the LMC point) yields a computable null distribution.
  #       theta_0 is exempt from the degenerate-candidate penalty inside MMCLRpval_fun, so if its
  #       null is fully degenerate (all draws fail even after the buffer + retry safety) the p-value
  #       is non-finite and would crash/garble the optimizer. Stop here instead. Under crn = TRUE
  #       this evaluation uses the optimizer's own worker configuration and CRN stream, so it is
  #       bit-identical to the optimizer's evaluation at theta_0: the guard is exact and the
  #       reported MMC p-value satisfies pval >= pval_0 by construction. Done BEFORE the optimizer
  #       seed is set below so it does not disturb the reproducible search trajectory.
  p_seed_workers <- if (isTRUE(con$crn)) con$workers else 0L
  p_seed <- obj_max(theta_0, mdl_h0, k1, LRT_0, con$N, con$burnin, p_seed_workers, con$lambda,
                    con$stationary_constraint, mdl_h1$control$thtol, Zsim, exog,
                    mdl_h0_null_cont, mdl_h1_null_cont, predrawn_eps, predrawn_state_rand)
  if (!is.finite(p_seed)) {
    stop("MMCLRTest: the null distribution could not be simulated at the initial parameter values (theta_0); the MMC p-value cannot be computed. Increase 'use_diff_init' or inspect the observed model fit.")
  }
  # Set optimizer seed for a reproducible search trajectory (PSO/GA internal
  # randomness; GenSA with a fixed 'par' is deterministic and unaffected). A drawn
  # seed is used instead of mc_seed+1 to avoid cross-replication stream collisions.
  # In the sequential legacy configuration (mc_seed = NULL, crn = FALSE,
  # workers = 0) no seed is set at all, exactly as in version 0.1.9: the
  # optimizer and estimation randomness continue the ambient (script-seeded)
  # stream, so script-seeded legacy runs reproduce 0.1.9 bit-for-bit.
  if (need_early_seeds) set.seed(internal_seeds[1L])
  if (p_seed >= con$threshold_stop){
    # ----- Early stop: the p-value at theta_0 already meets threshold_stop, so the
    # search cannot change the test decision (pval >= pval_0 >= threshold_stop by
    # construction). Skip the optimizer entirely -- exactly the threshold-stopping
    # semantics, at the cost of a single (LMC-type) evaluation. Note the reported
    # pval is then pval_0 itself (an optimizer run crossing the threshold could
    # report a weakly larger value; both are valid and imply the same decision).
    theta   <- theta_0
    pval    <- p_seed
    mmc_out <- list(early_stop = TRUE,
                    message = "threshold_stop reached at theta_0 (pval_0 >= threshold_stop); nuisance-parameter search skipped.")
    if (!isTRUE(con$silence)) message("MMCLRTest: threshold_stop reached at theta_0; skipping the nuisance-parameter search.")
  }else if (con$type=="pso"){
    # Set PSO specific controls
    con$optim_control$trace.stats <- TRUE
    con$optim_control$trace <- as.numeric(con$silence==FALSE)
    con$optim_control$abstol <- -con$threshold_stop
    con$optim_control$maxf <- con$maxit
    con$optim_control$REPORT <- 1
    # begin optimization
    mmc_out   <- pso::psoptim(par = theta_0, fn = obj_min, lower = theta_low, upper = theta_upp,
                              gr = NULL, control = con$optim_control,
                              mdl_h0 = mdl_h0, k1 = k1, LRT_0 = LRT_0, N = con$N, burnin = con$burnin, workers = con$workers,
                              lambda = con$lambda, stationary_constraint = con$stationary_constraint,
                              thtol = mdl_h1$control$thtol, Z = Zsim, exog = exog, mdl_h0_control = mdl_h0_null_cont,
                              mdl_h1_control = mdl_h1_null_cont,
                              predrawn_eps = predrawn_eps, predrawn_state_rand = predrawn_state_rand)
    theta     <- mmc_out$par
    pval      <- -mmc_out$value
  }else if(con$type=="GenSA"){
    # Set GenSA specific controls
    con$optim_control$trace.mat <- TRUE
    con$optim_control$verbose <- con$silence==FALSE
    con$optim_control$threshold.stop <- -con$threshold_stop
    con$optim_control$max.call <- con$maxit
    # begin optimization
    mmc_out   <- GenSA::GenSA(par = theta_0, fn = obj_min, lower = theta_low, upper = theta_upp,
                              control = con$optim_control,
                              mdl_h0 = mdl_h0, k1 = k1, LRT_0 = LRT_0, N = con$N, burnin = con$burnin, workers = con$workers,
                              lambda = con$lambda, stationary_constraint = con$stationary_constraint,
                              thtol = mdl_h1$control$thtol, Z = Zsim, exog = exog, mdl_h0_control = mdl_h0_null_cont,
                              mdl_h1_control = mdl_h1_null_cont,
                              predrawn_eps = predrawn_eps, predrawn_state_rand = predrawn_state_rand)
    theta     <- mmc_out$par
    pval      <- -mmc_out$value
  }else if(con$type=="GA"){
    # ----- GA controls: popSize defaults to 10 (GA's own default of 50 would make MMC perform
    #       popSize*maxiter expensive evaluations). Any extra GA control the user supplies via
    #       'optim_control' (e.g. popSize, pmutation, pcrossover, elitism, run) is passed through;
    #       keys that clash with the fixed MMC arguments are ignored.
    ga_ctrl <- con$optim_control
    if (is.null(ga_ctrl)) ga_ctrl <- list()
    if (is.null(ga_ctrl$popSize)) ga_ctrl$popSize <- 10
    ga_fixed <- list(type = "real-valued", fitness = obj_max,
                     mdl_h0 = mdl_h0, k1 = k1, LRT_0 = LRT_0, N = con$N, burnin = con$burnin, workers = con$workers,
                     lambda = con$lambda, stationary_constraint = con$stationary_constraint,
                     thtol = mdl_h1$control$thtol, Z = Zsim, exog = exog, mdl_h0_control = mdl_h0_null_cont,
                     mdl_h1_control = mdl_h1_null_cont,
                     predrawn_eps = predrawn_eps, predrawn_state_rand = predrawn_state_rand,
                     lower = theta_low, upper = theta_upp,
                     maxiter = con$maxit, maxFitness = con$threshold_stop,
                     monitor = (con$silence==FALSE), suggestions = t(theta_0))
    ga_ctrl   <- ga_ctrl[setdiff(names(ga_ctrl), names(ga_fixed))]
    mmc_out   <- do.call(GA::ga, c(ga_fixed, ga_ctrl))
    theta     <- as.matrix(mmc_out@solution[1,])
    pval      <- mmc_out@fitnessValue
  }else{
    stop("'type' must be one of 'pso', 'GenSA', or 'GA'.")
  }
  # ----- get test output using optimization output params
  if (pval < 0){
    stop("MMCLRTest: the nuisance-parameter search evaluated no admissible candidate (every candidate violated the transition-matrix or stationarity constraint, or could not be simulated), so no p-value is defined. Check that the search bounds admit the null model's parameter values.")
  }
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
  # ----- Reported test statistic: the FIXED observed LRT_0 computed above from the
  #       observed-data fits (Dufour 2006), NOT recomputed at the optimized nuisance value.
  # ----- organize test output
  MMCLRTest_output <- list(mdl_h0 = mdl_h0, mdl_h1 = mdl_h1, mdl_h0_mmc = mdl_h0_mmc, mdl_h1_mmc = mdl_h1,
                           LRT_0 = LRT_0, pval = pval,
                           theta_h0 = theta_h0, theta_h1 = theta_h1, control = con,
                           mmc_optimout = mmc_out,
                           pval_0 = p_seed)
  class(MMCLRTest_output) <- "MMCLRTest"
  return(MMCLRTest_output)
}



