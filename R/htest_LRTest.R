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


# Estimation for the SIMULATED side of the Monte Carlo tests, wrapped so a
# failed fit returns NULL instead of signaling an R error: the per-draw retry
# machinery in LR_samp_dist treats failure as "skip/retry", but an error
# raised through the C++ callback would still run every calling handler up
# the caller's stack before being absorbed, so condition-recording callers
# (testthat, logging harnesses built on withCallingHandlers) would record
# hundreds of phantom errors per penalized candidate. tryCatch's exiting
# handler pre-empts outer calling handlers, so nothing escapes.
.estimMdl_quiet <- function(Y, p, q, k, Z = NULL, control = list()){
  tryCatch(estimMdl(Y, p, q, k, Z, control), error = function(e) NULL)
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
#' @section Numerical policy for the LR statistic: For nested models the LR statistic is non-negative by construction (the null fit embeds into the alternative parameter space with equal likelihood), so a negative computed value can only be a numerical optimization artifact. Both the observed statistic and every simulated statistic are therefore floored at 0 (identically on both sides, which preserves the exchangeability underlying the Monte Carlo p-value); a warning reports the magnitude when it is non-negligible. For \code{k0 > 1}, the null log-likelihood is additionally floored at the one-regime (linear) fit on both sides. Ties at 0 are handled by the randomized tie-breaking rule of Dufour (2006) already implemented in \code{\link{MCpval}}; ties at the bottom of the support can only make the test conservative. The regime variances are bounded away from zero
#'   during estimation (control \code{em_variance_constraint} of the model constructors, default \code{0.01} of the one-regime fit variance), because the Markov-switching likelihood is otherwise unbounded and degenerate solutions produce spurious likelihood-ratio spikes; the statistic is therefore a constrained-likelihood ratio, computed by the same rule for the observed and the simulated samples, which is all the Monte Carlo argument requires. Transition probabilities are not constrained.
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
  if (k1 <= k0){
    stop("LMCLRTest: k1 must be greater than k0 (received k0 = ", k0, ", k1 = ", k1, ").")
  }
  # ----- Seed RNG for reproducible estimation (random initial values in EM)
  if (!is.null(con$mc_seed)) set.seed(con$mc_seed)
  # ----- Estimate models using observed data
  mdl_h0 <- estimMdl(Y, p, q, k0, Z, con$mdl_h0_control)
  # Resolve the alternative's switching flags before the warm starts are built and
  # before the alternative is estimated, so the two cannot disagree.
  con$mdl_h1_control <- .inherit_ms_flags(con$mdl_h1_control, mdl_h0)
  con$mdl_h1_control <- .inherit_stationarity_margin(con$mdl_h1_control, con$mdl_h0_control, "LMCLRTest", p)
  if ((con$init_method == "warmstart") && (k1 > 1)){
    # Warm-start the alternative model from the null fit (regime-duplication
    # embedding, see warmstart_theta); deterministic, so the RNG stream and the
    # random starts are unaffected. The same rule is applied per draw when
    # simulating the null distribution (see LR_samp_dist).
    # `[[` rather than `$`: the latter partial-matches, so a mistyped control name
    # would be read as the flag and silently build a warm start for the wrong space.
    con$mdl_h1_control$init_theta_extra <- warmstart_theta(mdl_h0, k1,
                                                           msmu = con$mdl_h1_control[["msmu"]],
                                                           msvar = con$mdl_h1_control[["msvar"]])
  }
  mdl_h1 <- estimMdl(Y, p, q, k1, Z, con$mdl_h1_control)
  # The flags are a proxy; what has to hold is that a warm start has the length the
  # alternative expects, which catches any future cause of a mismatch.
  ws_extra <- con$mdl_h1_control$init_theta_extra
  if (!is.null(ws_extra) && (length(ws_extra) > 0) &&
      any(vapply(ws_extra, length, 0L) != length(mdl_h1$theta))){
    warning("the warm starts did not match the alternative model's parameter vector ",
            "and were discarded.", call. = FALSE)
  }
  con$mdl_h0_control <- mdl_h0$control
  con$mdl_h1_control <- mdl_h1$control
  # The stationarity margin must be the same on both sides: with different
  # margins the alternative no longer nests the null's constrained space and
  # the LR statistic loses its constrained-ML meaning.
  d0 <- con[["mdl_h0_control"]][["stationarity_constraint"]]
  d1 <- con[["mdl_h1_control"]][["stationarity_constraint"]]
  if (!is.null(d0) && !is.null(d1) && !isTRUE(all.equal(as.numeric(d0), as.numeric(d1)))){
    warning("'stationarity_constraint' differs between mdl_h0_control (", d0, ") and mdl_h1_control (", d1, "); the same margin is required on both sides for the statistic to be a constrained-likelihood ratio.")
  }
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
    # The stationarity margin, when set on the null control, carries into the
    # floor fit: the floor's nesting rationale needs the floor point feasible
    # in the constrained space, by the same rule on both sides.
    mdl_lin_con <- list(getSE = FALSE)
    if (!is.null(con[["mdl_h0_control"]][["stationarity_constraint"]])){
      mdl_lin_con$stationarity_constraint <- con[["mdl_h0_control"]][["stationarity_constraint"]]
    }
    mdl_lin <- estimMdl(Y, p, q, 1, Z, mdl_lin_con)
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


# ---------------------------------------------------------------------------- #
# Column-stochastic parameterization of the MMC transition-matrix search.
#
# The optimizer would otherwise perturb the k0^2 entries of P independently, so
# a candidate satisfies "every column sums to 1" only by a measure-zero
# coincidence and nearly every candidate fails the colsum gate. Instead, only
# the free entries of each column (k0-1 of them) are searched; the remaining
# entry is derived as 1 - sum(free). theta itself, theta_P_ind, SEs, names, and
# every reported field stay in full k0^2 coordinates -- only the optimizer's
# search vector is reduced.

# Choose the entry of each P column to derive (never searched), fixed from
# theta_0 for the whole search. The largest entry in the column is derived,
# which maximizes the room the free entries have before it leaves [0,1].
# 'fixed_full_idx' excludes positions a caller has pinned (the boundary face
# search pins the absorbing column and the clean block) from being eligible.
# A column ALL of whose entries are pinned needs no derived entry -- its
# values come from the pins, and the caller asserts once that they sum to 1
# -- so it is marked NA and skipped by the reassembly and feasibility steps.
.mmc_derived_P_index <- function(theta_0, theta_P_ind, k0, fixed_full_idx = integer(0)){
  full_idx <- which(as.numeric(theta_P_ind) == 1)
  if (length(full_idx) != k0*k0){
    stop("MMC_bounds: theta_P_ind does not mark k0^2 entries of theta; cannot reduce the transition-matrix search.")
  }
  # vec(P) is stored column-major in theta, so filling a k0 x k0 matrix
  # column-major with these positions recovers exactly which theta index is
  # column j, row i of P.
  pos_mat <- matrix(full_idx, k0, k0)
  derived_full_idx <- rep(NA_integer_, k0)
  for (j in seq_len(k0)){
    col_positions <- pos_mat[, j]
    eligible <- setdiff(col_positions, fixed_full_idx)
    if (length(eligible) > 0L){
      derived_full_idx[j] <- eligible[which.max(theta_0[eligible])]
    }
  }
  list(pos_mat = pos_mat, derived_full_idx = derived_full_idx, full_idx = full_idx)
}

# Full-length theta from a reduced (free-coordinates-only) vector: write any
# pinned coordinates first, copy the free entries into place, then set each
# column's derived entry to 1 - sum(the OTHER entries of that column, pinned
# included -- the loop sums every non-derived position, so writing pins first
# is sufficient). A column all of whose non-derived entries are free (the
# ordinary case) reassembles exactly, to machine precision (verified at 0.0
# over 200k random draws at k0 in {3,4,5,8} during the spec's stage-0 review).
# A column with no derived entry (NA: every entry pinned, boundary face mode)
# is skipped -- its values are the pins, asserted to sum to 1 where they are
# constructed. Defaults reproduce the pre-pinning behavior exactly.
.mmc_reassemble_theta <- function(theta_reduced, reduced_idx, derived_full_idx, pos_mat, npar,
                                  pinned_idx = integer(0), pinned_vals = numeric(0)){
  theta_full <- numeric(npar)
  if (length(pinned_idx) > 0) theta_full[pinned_idx] <- pinned_vals
  theta_full[reduced_idx] <- as.numeric(theta_reduced)
  for (j in seq_len(ncol(pos_mat))){
    if (is.na(derived_full_idx[j])) next
    col_positions <- pos_mat[, j]
    free_in_col <- setdiff(col_positions, derived_full_idx[j])
    theta_full[derived_full_idx[j]] <- 1 - sum(theta_full[free_in_col])
  }
  theta_full
}

# Wraps an already-CRN-wrapped MMC objective (obj_min, which returns -pval, or
# obj_max, which returns pval) so it accepts the REDUCED search vector: reassemble
# to full theta, and if a derived entry falls outside the admissible region OR
# its own eps/CI box (the box MMC_bounds already built for that position, so no
# separate computation is needed), penalize the candidate WITHOUT calling the
# wrapped objective -- consuming no random numbers, at the ~10-flop cost the
# spec anticipated, exactly mirroring the colsum-gate rejection this replaces.
# A naive clamp instead of a penalty would not just move the candidate but break
# the column-sum identity itself (free entries can sum to more than 1, driving
# the clamped derived entry negative), which trips the simulators' own hard
# check and burns simulation draws before being penalized anyway; unreachable in
# practice under the largest-entry rule, so penalizing costs nothing.
.mmc_reduced_wrap <- function(objfn_full, reduced_idx, derived_full_idx, pos_mat, npar,
                              theta_low, theta_upp, P_low, P_upp, penalty_value,
                              pinned_idx = integer(0), pinned_vals = numeric(0)){
  force(objfn_full); force(reduced_idx); force(derived_full_idx); force(pos_mat)
  force(npar); force(theta_low); force(theta_upp); force(P_low); force(P_upp); force(penalty_value)
  force(pinned_idx); force(pinned_vals)
  # Columns with no derived entry (NA: fully pinned, boundary face mode) carry
  # no feasibility check here -- their values are constants asserted to sum to
  # 1 where they are constructed. Indexing with the NA left in would propagate
  # NA into the feasibility comparison and crash `if (infeasible)`.
  derived_chk <- derived_full_idx[!is.na(derived_full_idx)]
  function(theta_reduced, ...){
    theta_full <- .mmc_reassemble_theta(theta_reduced, reduced_idx, derived_full_idx, pos_mat, npar,
                                        pinned_idx, pinned_vals)
    dv <- theta_full[derived_chk]
    infeasible <- any(dv < P_low) || any(dv > P_upp) ||
      any(dv < theta_low[derived_chk]) || any(dv > theta_upp[derived_chk])
    if (infeasible) return(penalty_value)
    objfn_full(theta_full, ...)
  }
}

# A per-coordinate box cannot express positive definiteness, so for a
# multivariate model the searched set is the box intersected with the
# admissible parameter space: a candidate whose assembled regime covariance is
# not positive definite is excluded from the sup -- penalized like the other
# inadmissible-candidate gates (correct sign for the objective's direction,
# never a value that can win the max), before any simulation and consuming no
# random numbers. This wraps the FULL-coordinate objective, so the reduced-
# search layer reassembles first and this check sees complete vech segments.
# The reference evaluations (the fitted theta, the face projection) bypass the
# search wrappers entirely and are exempt, as they are for the other gates.
.mmc_admissible_wrap <- function(objfn, psd_blocks, q, penalty_value, count_env){
  force(objfn); force(psd_blocks); force(q); force(penalty_value); force(count_env)
  function(theta_full, ...){
    for (bidx in psd_blocks){
      Sig <- covar_unvech(as.numeric(theta_full)[bidx], q)
      lam <- eigen((Sig + t(Sig)) / 2, symmetric = TRUE, only.values = TRUE)$values
      if (min(lam) <= .mmc_psd_tol(max(lam))){
        count_env$pd_hits <- count_env$pd_hits + 1L
        return(penalty_value)
      }
    }
    objfn(theta_full, ...)
  }
}

# With an estimation stationarity margin in force (box_match_constraints), the
# searched set is the box intersected with {companion spectral radius <= 1 -
# delta}. At p = 1 the box clamp already enforces it; for p >= 2 that set is
# not a box, so a candidate outside it is excluded here before any simulation,
# consuming no random numbers, like the other inadmissible-candidate gates.
# Wraps the FULL-coordinate objective; the reference evaluation at the fitted
# theta bypasses it (a constrained fit satisfies the margin by construction).
.mmc_margin_wrap <- function(objfn, phi_idx, radius_bound, penalty_value, count_env){
  force(objfn); force(phi_idx); force(radius_bound); force(penalty_value); force(count_env)
  function(theta_full, ...){
    ph <- as.numeric(theta_full)[phi_idx]
    if (!all(is.finite(ph)) || .phi_companion_radius(ph) > radius_bound + 1e-12){
      count_env$margin_hits <- count_env$margin_hits + 1L
      return(penalty_value)
    }
    objfn(theta_full, ...)
  }
}

# Records every candidate the optimizer evaluates (reduced coordinates) and its
# objective value into 'rec_env' before delegating. A pure bookkeeping closure:
# consumes no random numbers and changes no numeric result, so it runs
# unconditionally and the harvest fields (feval, evaluated_range) are exact for
# every optimizer -- pso's own trace stores candidates only when trace > 0
# (i.e. never under silence = TRUE, exactly the batch configuration), and
# GenSA/GA never store them at all.
.mmc_record_wrap <- function(objfn, rec_env){
  force(objfn); force(rec_env)
  function(theta_reduced, ...){
    val <- objfn(theta_reduced, ...)
    n <- rec_env$n + 1L
    rec_env$n <- n
    rec_env$X[[n]] <- as.numeric(theta_reduced)
    rec_env$v[n] <- as.numeric(val)[1L]
    val
  }
}
# ---------------------------------------------------------------------------- #

# 'eps' may be a scalar (recycled to every coordinate) or a vector of
# length(theta_0), so the search half-width can differ per parameter. A length
# that merely divides length(theta_0) evenly would recycle silently with no
# warning, so only exactly 1 or length(theta_0) is accepted. eps = 0 is valid
# (typically relying entirely on CI_union to widen the box) but GenSA specifically
# treats a genuinely zero-width bound as invalid ("Lower and upper bounds are not
# consistent"), unlike pso/GA, which handle it fine; when the chosen optimizer is
# GenSA and CI_union (or a phi_low/P_low-type clamp) hasn't already widened a
# zero entry by the time the final box is built, that entry is nudged to 1e-8,
# which GenSA accepts. Returns the (possibly nudged) eps, so callers must use
# the return value, not the original input.
.validate_mmc_eps <- function(eps, n_theta, fn_name = "MMC_bounds", method = NULL){
  if (!is.numeric(eps) || !(length(eps) %in% c(1L, n_theta))){
    stop(fn_name, ": 'eps' must have length 1 or length(theta_0) = ", n_theta,
         " (received length ", length(eps), "). A length that merely divides ",
         n_theta, " evenly is not accepted: it would recycle silently rather than ",
         "apply to the parameter block a user likely intended.")
  }
  if (any(!is.finite(eps))){
    stop(fn_name, ": 'eps' must be finite (received ",
         paste(eps[!is.finite(eps)], collapse = ", "), ").")
  }
  if (any(eps < 0)){
    stop(fn_name, ": 'eps' must be non-negative (received ",
         paste(eps[eps < 0], collapse = ", "), ").")
  }
  if (isTRUE(method == "GenSA") && any(eps == 0)){
    eps[eps == 0] <- 1e-8
  }
  eps
}

# A box-shaping control (phi_low/phi_upp/P_low/P_upp/variance_constraint) that is
# non-finite -- an NA typo, most plausibly -- would otherwise flow silently into
# theta_low/theta_upp and only surface later as an opaque `if (NA)` crash in the
# box-shape checks below, far from the actual mistake. allow_null = FALSE is for
# controls that always carry a real default (e.g. variance_constraint = 0.01) and
# are unconditionally used, so NULL itself is also a mistake, not a "skip" signal.
.validate_box_control <- function(x, nm, fn_name, allow_null = TRUE){
  if (is.null(x)){
    if (allow_null) return(x)
    stop(fn_name, ": '", nm, "' must not be NULL.")
  }
  if (!is.numeric(x) || any(!is.finite(x))){
    stop(fn_name, ": '", nm, "' must be finite (received ",
         paste(x, collapse = ", "), ").")
  }
  x
}

# 'box_include' supplies per-parameter values the search box must contain (NA =
# leave that coordinate alone). Length must be exactly length(theta_0) -- no
# recycling of any kind, since a stretch aimed at the wrong coordinate is the
# failure mode -- and every non-NA entry finite.
.validate_box_include <- function(x, n_theta, fn_name = "MMC_bounds"){
  if (is.null(x)) return(x)
  # rep(NA, n) is logical, and the natural way to write an all-NA vector;
  # accept it as the numeric all-NA no-op it clearly means.
  if (is.logical(x) && all(is.na(x))) x <- as.numeric(x)
  if (!is.numeric(x)){
    stop(fn_name, ": 'box_include' must be a numeric vector (received ",
         class(x)[1], "), with NA for coordinates that are not stretched.")
  }
  if (length(x) != n_theta){
    stop(fn_name, ": 'box_include' must have length(theta_0) = ",
         n_theta, " (received length ", length(x), "), with NA for coordinates ",
         "that are not stretched.")
  }
  if (any(is.infinite(x) | is.nan(x))){
    stop(fn_name, ": 'box_include' entries must be finite or NA (received ",
         paste(x[is.infinite(x) | is.nan(x)], collapse = ", "), ").")
  }
  x
}

# The moment-anchored ranges R for the boundary search's rough block. Kept
# overridable but never auto-narrowed: the size bound of the boundary theory is
# alpha + Pr[theta_r0 outside R], so a too-narrow R errs anti-conservative.
.validate_boundary_R <- function(bR, fn_name = "MMCLRTest"){
  def <- list(mu_mult = 5, sigma_mult = 4, p_range = c(0.05, 0.95))
  if (is.null(bR)) bR <- list()
  if (!is.list(bR)) stop(fn_name, ": 'boundary_R' must be a list.")
  bad <- setdiff(names(bR), names(def))
  if (length(bad) > 0) stop(fn_name, ": unknown names in boundary_R: ", paste(bad, collapse = ", "))
  bR <- utils::modifyList(def, bR)
  if (!is.numeric(bR$mu_mult) || length(bR$mu_mult) != 1 || !is.finite(bR$mu_mult) || bR$mu_mult <= 0)
    stop(fn_name, ": boundary_R$mu_mult must be a single finite positive number.")
  if (!is.numeric(bR$sigma_mult) || length(bR$sigma_mult) != 1 || !is.finite(bR$sigma_mult) || bR$sigma_mult <= 0)
    stop(fn_name, ": boundary_R$sigma_mult must be a single finite positive number.")
  if (!is.numeric(bR$p_range) || length(bR$p_range) != 2 || any(!is.finite(bR$p_range)) ||
      bR$p_range[1] <= 0 || bR$p_range[2] >= 1 || bR$p_range[1] >= bR$p_range[2])
    stop(fn_name, ": boundary_R$p_range must be two increasing values strictly inside (0, 1).")
  bR
}

# Index sets for the boundary search's blocks, on a CANONICALIZED fit
# (canonicalize_regimes: transient regime first, absorbing regime last).
# The rough block is what only the transient era identifies: the SWITCHING
# parameters' transient entries plus the transient column of P. Parameters the
# whole sample identifies -- a non-switching common mean/variance, the AR
# coefficients, the exogenous-regressor coefficients (betaZ never switches),
# and the absorbing regime's own parameters -- are clean. For a multivariate
# fit the switching blocks are per-regime vech(Sigma) segments; the transient
# segment's internal structure (which entries are variances, and which series
# pair each covariance couples) is returned so the moment-anchored ranges can
# be built per coordinate.
.mmc_boundary_blocks <- function(mdl_h0){
  k0 <- mdl_h0$k
  q  <- mdl_h0$q
  Nsig <- q * (q + 1) / 2
  P_idx <- which(as.numeric(mdl_h0$theta_P_ind) == 1)
  pos_mat <- matrix(P_idx, k0, k0)
  mu_idx  <- which(as.numeric(mdl_h0$theta_mu_ind) == 1)
  sig_idx <- which(as.numeric(mdl_h0$theta_sig_ind) == 1)
  rough_mu <- integer(0); rough_sig <- integer(0); clean_idx <- integer(0)
  if (isTRUE(mdl_h0$msmu)){
    rough_mu  <- mu_idx[seq_len((k0 - 1) * q)]
    clean_idx <- c(clean_idx, mu_idx[((k0 - 1) * q + 1):(k0 * q)])
  } else {
    clean_idx <- c(clean_idx, mu_idx)
  }
  if (isTRUE(mdl_h0$msvar)){
    rough_sig <- sig_idx[seq_len((k0 - 1) * Nsig)]
    clean_idx <- c(clean_idx, sig_idx[((k0 - 1) * Nsig + 1):(k0 * Nsig)])
  } else {
    clean_idx <- c(clean_idx, sig_idx)
  }
  if (mdl_h0$p > 0){
    clean_idx <- c(clean_idx, which(as.numeric(mdl_h0$theta_phi_ind) == 1))
  }
  if (!is.null(mdl_h0$theta_x_ind)){
    clean_idx <- c(clean_idx, which(as.numeric(mdl_h0$theta_x_ind) == 1))
  }
  # vech(Sigma) stacks the lower triangle column by column (covar_vech), so
  # the r-th entry of a regime's segment couples the series pair vech_pairs[r, ].
  # theta_var_ind marks the variance (diagonal) entries within each segment.
  # For every rough off-diagonal entry, the theta POSITIONS of the two
  # variances it couples are resolved within that entry's OWN segment, so the
  # Cauchy-Schwarz ranges downstream are per-regime by construction, for any
  # number of transient regimes.
  vech_pairs <- which(lower.tri(diag(q), diag = TRUE), arr.ind = TRUE)
  vech_pairs <- vech_pairs[order(vech_pairs[, "col"], vech_pairs[, "row"]), , drop = FALSE]
  var_ind <- as.numeric(mdl_h0$theta_var_ind)
  rough_sig_diag <- rough_sig[var_ind[rough_sig] == 1]
  rough_sig_off  <- rough_sig[var_ind[rough_sig] == 0]
  diag_within <- which(vech_pairs[, 1] == vech_pairs[, 2])   # vech position of (i,i), i = 1..q
  off_within  <- which(vech_pairs[, 1] != vech_pairs[, 2])
  n_seg <- length(rough_sig) %/% Nsig
  off_diag_i <- integer(0); off_diag_j <- integer(0)
  for (m in seq_len(n_seg)){
    seg <- rough_sig[((m - 1) * Nsig + 1):(m * Nsig)]
    off_diag_i <- c(off_diag_i, seg[diag_within[vech_pairs[off_within, 1]]])
    off_diag_j <- c(off_diag_j, seg[diag_within[vech_pairs[off_within, 2]]])
  }
  list(rough_mu = rough_mu, rough_sig = rough_sig,
       rough_sig_diag = rough_sig_diag, rough_sig_off = rough_sig_off,
       off_diag_i = off_diag_i, off_diag_j = off_diag_j,
       Nsig = Nsig, n_transient_seg = n_seg, clean_idx = clean_idx,
       trans_P_idx = as.integer(pos_mat[, seq_len(k0 - 1)]),
       abs_diag_idx = pos_mat[k0, k0],
       abs_off_idx  = pos_mat[seq_len(k0 - 1), k0],
       pos_mat = pos_mat)
}

# The moment anchors themselves, from the observed estimation sample, per
# series. The variance lower end is the EM regime-variance floor (already the
# admissible lower bound of the estimation itself; for q > 1 the package
# anchors it at em_variance_constraint * tr(Sigma_lin)/q, one scalar); when
# the floor is disabled it is stored as 0, and the fallback is the same
# variance_constraint * theta_0 lower bound the ordinary box enforces (per
# coordinate, computed by the caller, so this returns NA there and the caller
# substitutes). Covariance off-diagonal ranges are NOT built here: they are
# +/- sqrt(U_i * U_j) with U_i the REALIZED diagonal upper ends, so they are
# derived by the caller after the diagonal ranges are final -- under a common
# multiplier that equals sigma_mult * s_i * s_j, but it stays the correct
# Cauchy-Schwarz object when the realized uppers differ per coordinate. Any
# positive semi-definite matrix whose diagonals lie within their ranges then
# has its off-diagonals within these automatically, so the maintained
# containment assumption for the covariance block reduces to the diagonal
# statement.
.mmc_R_ranges <- function(mdl_h0, bR, fn_name = "MMC_bounds"){
  ym <- as.matrix(mdl_h0$y)
  ybar <- colMeans(ym)
  s_y  <- apply(ym, 2, stats::sd)
  sig_low <- if (isTRUE(mdl_h0$sigma_floor > 0)) mdl_h0$sigma_floor else NA_real_
  sig_upp <- bR$sigma_mult * s_y^2
  if (is.finite(sig_low) && any(sig_low > sig_upp)){
    stop(fn_name, ": the regime-variance floor (", signif(sig_low, 6),
         ") exceeds a transient-variance range upper end (",
         paste(signif(sig_upp[sig_low > sig_upp], 6), collapse = ", "),
         "); the moment-anchored variance range is empty. This can only happen ",
         "on a degenerate sample or with a very small boundary_R$sigma_mult.")
  }
  list(mu_low = ybar - bR$mu_mult * s_y, mu_upp = ybar + bR$mu_mult * s_y,
       sig_low = sig_low, sig_upp = sig_upp,
       p_low = bR$p_range[1], p_upp = bR$p_range[2])
}

# Smallest-eigenvalue admissibility tolerance for an assembled covariance
# candidate: strictly positive and scaled to the matrix, sitting just above the
# Cholesky feasibility boundary -- anything admitted is simulable, and only
# nearly-singular candidates at rounding level are excluded. Deliberately NOT
# the EM eigenvalue floor: the floor is an estimation-space device the true
# covariance need not satisfy (a matrix with in-range diagonals can have a
# smallest eigenvalue far below it), and gating on it would silently
# strengthen the maintained containment assumption beyond the diagonal
# statement the size bound reports.
.mmc_psd_tol <- function(lam_max){
  1e-12 * max(lam_max, .Machine$double.xmin)
}

# Deterministic positive-definiteness repair by off-diagonal shrinkage:
# Sigma(s) = D + s * O with D the (positive, in-range) diagonal part and O the
# off-diagonal part. The repaired matrix must CLEAR THE SEARCH GATE'S OWN
# admissibility criterion (.mmc_admissible_wrap tests the smallest eigenvalue
# against .mmc_psd_tol of the largest), with an interior safety margin --
# landing merely on the positive side of singular would hand the optimizer a
# seed its first evaluation rejects. lambda_min(Sigma(s)) is concave in s and
# lambda_max convex (each an extremum over linear functions of s), so the
# margin criterion g(s) = lambda_min - margin * tol(lambda_max) is concave
# with g(0) > 0 for any positive diagonal: the admissible s form an interval
# containing 0, and the largest admissible s is found by bisection. Shrinking
# toward the diagonal keeps every coordinate inside its own range (the
# off-diagonal ranges are symmetric about 0 and the diagonals are untouched),
# which an eigenvalue clip would not.
.psd_repair_shrink <- function(Sig, iter = 60L, margin = 10){
  D <- diag(diag(Sig))
  O <- Sig - D
  clears_gate <- function(s){
    m <- D + s * O
    lam <- eigen((m + t(m)) / 2, symmetric = TRUE, only.values = TRUE)$values
    min(lam) > margin * .mmc_psd_tol(max(lam))
  }
  if (clears_gate(1)) return(list(Sigma = Sig, s = 1))
  lo <- 0; hi <- 1
  for (i in seq_len(iter)){
    mid <- (lo + hi) / 2
    if (clears_gate(mid)) lo <- mid else hi <- mid
  }
  list(Sigma = D + lo * O, s = lo)
}

# The boundary projection of a (canonicalized) fitted theta: the absorbing
# diagonal set to exactly 1 and its column's off-diagonals to exactly 0,
# everything else at its estimate. Shared by the face-mode search reference
# point and aBfixLMCTest's second component so the two can never drift apart.
.theta_boundary_projection <- function(theta_0, blocks){
  th <- as.numeric(theta_0)
  th[blocks$abs_diag_idx] <- 1
  th[blocks$abs_off_idx]  <- 0
  names(th) <- names(theta_0)
  th
}

# Face-mode box: the face set alone -- rough block over R, clean block pinned
# at its estimates, absorbing column pinned at (1, 0s). Does NOT go through
# MMC_bounds' eps/CI construction (none of it applies here); runs the same
# finiteness and inverted-box checks, with the in-box reference point being
# the CLAMPED projection (theta_hat itself legitimately lies outside the face).
.mmc_face_bounds <- function(mdl_h0, con){
  theta_0 <- mdl_h0$theta
  blocks  <- .mmc_boundary_blocks(mdl_h0)
  R <- .mmc_R_ranges(mdl_h0, con$boundary_R, fn_name = "MMCLRTest")
  theta_low <- as.numeric(theta_0)
  theta_upp <- as.numeric(theta_0)
  theta_low[blocks$rough_mu] <- R$mu_low
  theta_upp[blocks$rough_mu] <- R$mu_upp
  sig_low <- if (is.finite(R$sig_low)) rep(R$sig_low, length(blocks$rough_sig_diag))
             else con$variance_constraint * as.numeric(theta_0)[blocks$rough_sig_diag]
  theta_low[blocks$rough_sig_diag] <- sig_low
  theta_upp[blocks$rough_sig_diag] <- R$sig_upp
  # Covariance off-diagonals: +/- sqrt(U_i * U_j) from the coupled variances'
  # upper ends, resolved within each entry's own regime segment
  # (Cauchy-Schwarz; see .mmc_R_ranges).
  if (length(blocks$rough_sig_off) > 0){
    off_bound <- sqrt(theta_upp[blocks$off_diag_i] * theta_upp[blocks$off_diag_j])
    theta_low[blocks$rough_sig_off] <- -off_bound
    theta_upp[blocks$rough_sig_off] <-  off_bound
  }
  theta_low[blocks$trans_P_idx] <- R$p_low
  theta_upp[blocks$trans_P_idx] <- R$p_upp
  pinned_idx  <- c(blocks$clean_idx, blocks$abs_diag_idx, blocks$abs_off_idx)
  pinned_vals <- c(as.numeric(theta_0)[blocks$clean_idx], 1, rep(0, length(blocks$abs_off_idx)))
  theta_low[pinned_idx] <- pinned_vals
  theta_upp[pinned_idx] <- pinned_vals
  # The pinned absorbing column must sum to 1: the reassembly and feasibility
  # layers rely on it and no longer check it per evaluation, so it is asserted
  # once here, from the constructed values themselves (protects against any
  # future edit to the pin construction above).
  abs_col_sum <- sum(theta_low[c(blocks$abs_diag_idx, blocks$abs_off_idx)])
  if (abs(abs_col_sum - 1) > 1e-12){
    stop("MMCLRTest: internal error -- the pinned absorbing column does not sum to 1.")
  }
  proj_unclamped <- .theta_boundary_projection(theta_0, blocks)
  proj_clamped   <- pmin(pmax(proj_unclamped, theta_low), theta_upp)
  # Coordinate-wise clamping can hand the optimizer a non-positive-definite
  # transient covariance as its starting point (pso silently drops an
  # inadmissible seed), so the clamped seed's transient block is repaired
  # deterministically by off-diagonal shrinkage, which stays inside the box.
  # The unclamped projection is untouched: it is an EM estimate's boundary
  # projection, admissible by construction.
  seed_repair <- list(repaired = FALSE, shrink = 1)
  q <- mdl_h0$q
  if (q > 1 && length(blocks$rough_sig) > 0){
    for (m in seq_len(blocks$n_transient_seg)){
      seg <- blocks$rough_sig[((m - 1) * blocks$Nsig + 1):(m * blocks$Nsig)]
      Sig_seed <- covar_unvech(proj_clamped[seg], q)
      lam <- eigen((Sig_seed + t(Sig_seed)) / 2, symmetric = TRUE, only.values = TRUE)$values
      if (min(lam) <= .mmc_psd_tol(max(lam))){
        rep_out <- .psd_repair_shrink(Sig_seed)
        proj_clamped[seg] <- covar_vech(rep_out$Sigma)
        seed_repair <- list(repaired = TRUE,
                            shrink = min(seed_repair$shrink, rep_out$s))
      }
    }
  }
  names(proj_clamped) <- names(theta_0)
  if (any(!is.finite(theta_low)) || any(!is.finite(theta_upp))){
    bad <- which(!is.finite(theta_low) | !is.finite(theta_upp))
    nm  <- names(theta_0)[bad]; if (is.null(nm)) nm <- paste("position", bad)
    stop("MMCLRTest: the face search box is non-finite for ", paste(nm, collapse = ", "), ".")
  }
  if (any(theta_low > theta_upp)){
    bad <- which(theta_low > theta_upp)
    nm  <- names(theta_0)[bad]; if (is.null(nm)) nm <- paste("position", bad)
    stop("MMCLRTest: the face search box has a lower bound above its upper bound for ",
         paste(nm, collapse = ", "), "; check 'boundary_R'.")
  }
  if (any(proj_clamped < theta_low - 1e-12 | proj_clamped > theta_upp + 1e-12)){
    stop("MMCLRTest: internal error -- the clamped face projection falls outside the face box.")
  }
  list(theta_low = theta_low, theta_upp = theta_upp,
       pinned_idx = pinned_idx, pinned_vals = pinned_vals,
       proj_unclamped = proj_unclamped, proj_clamped = proj_clamped,
       seed_repair = seed_repair)
}

#' @title MMC nuisance parameter bounds
#'
#' @description This function is used to determine the lower and upper bounds for the MMC LRT parameter search.
#'
#' @param mdl_h0 List with restricted model properties.
#' @param con List with control options provided to MMC LRT procedure. \code{con$eps} may be a scalar, recycled to every coordinate, or a vector of length \code{length(mdl_h0$theta)} giving a different search half-width per parameter, in the same order as \code{names(mdl_h0$theta)}. Every entry must be finite and non-negative; \code{eps = 0} relies entirely on \code{CI_union} to widen the box, and when \code{con$type = "GenSA"} a zero entry is nudged to \code{1e-8} (GenSA treats a genuinely zero-width bound as invalid, unlike \code{pso}/\code{GA}). \code{con$phi_low}/\code{con$phi_upp}/\code{con$P_low}/\code{con$P_upp}, if not \code{NULL}, and \code{con$variance_constraint} must all be finite.
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
  eps <- .validate_mmc_eps(con$eps, length(theta_0), method = con$type)
  con$phi_low <- .validate_box_control(con$phi_low, "phi_low", "MMC_bounds")
  con$phi_upp <- .validate_box_control(con$phi_upp, "phi_upp", "MMC_bounds")
  con$P_low   <- .validate_box_control(con$P_low, "P_low", "MMC_bounds")
  con$P_upp   <- .validate_box_control(con$P_upp, "P_upp", "MMC_bounds")
  con$variance_constraint <- .validate_box_control(con$variance_constraint, "variance_constraint", "MMC_bounds", allow_null = FALSE)
  # ----- Define lower & upper bounds for search
  theta_low = theta_0 - eps
  theta_upp = theta_0 + eps
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
  # 'box_include': per-parameter values the box must contain (NA = untouched),
  # applied after the eps/CI union -- it only widens, and only where the union
  # does not already reach the value. The admissible-region clamps below keep
  # their existing semantics and can still remove an inadmissible value; that
  # contradiction is caught loudly after the clamps rather than silently
  # overridden (checked at the end of this function).
  box_inc <- .validate_box_include(con$box_include, length(theta_0))
  if (!is.null(box_inc)){
    inc_idx <- which(!is.na(box_inc))
    theta_low[inc_idx] <- pmin(c(theta_low)[inc_idx], box_inc[inc_idx])
    theta_upp[inc_idx] <- pmax(c(theta_upp)[inc_idx], box_inc[inc_idx])
    # With the search set aligned to an estimation stationarity margin, a
    # requested autoregressive point outside the margin can never be searched
    # (the margin gate excludes it), so the request is contradictory and errors
    # like a requested value excluded by an admissible-region clamp.
    if (isTRUE(con$box_match_constraints) && mdl_h0$p > 0 && isTRUE(mdl_h0$q == 1)){
      sc_inc <- .validate_stationarity_constraint(con[["mdl_h0_control"]][["stationarity_constraint"]], "MMC_bounds")
      if (sc_inc > 0){
        phi_pos_inc <- which(as.numeric(mdl_h0$theta_phi_ind) == 1)
        phi_inc <- box_inc[phi_pos_inc]
        phi_inc[is.na(phi_inc)] <- as.numeric(theta_0)[phi_pos_inc][is.na(phi_inc)]
        if (.phi_companion_radius(phi_inc) > 1 - sc_inc + 1e-12){
          stop("MMC_bounds: 'box_include' asks the search box to contain an autoregressive point with companion spectral radius ", signif(.phi_companion_radius(phi_inc), 6),
               ", outside the stationarity margin (radius <= ", 1 - sc_inc, ") that 'stationarity_constraint' imposes on every fit -- contradictory controls. Set box_match_constraints = FALSE to search beyond the estimation constraint set.")
        }
      }
    }
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
  # Boundary "union" stretch (C^bound_T-MMC): applied AFTER every clamp above,
  # deliberately -- the moment-anchored ranges R are load-bearing for the size
  # bound (alpha + Pr[theta_r0 outside R]) and must never be narrowed by the
  # conditional variance rewrite or the P/phi clamps. A user clamp that
  # contradicts the stretch (P_upp < 1, P_low > 0) is rejected by the caller's
  # scope checks before this point. The rough block's ranges are the hull of
  # the ordinary box and R; the clean block keeps its ordinary box; the
  # absorbing column is stretched so the face (p_AA = 1, off-diagonals 0) is
  # reachable whichever entry the reduced-P search derives.
  if (isTRUE(con$box_boundary == "union")){
    con$boundary_R <- .validate_boundary_R(con$boundary_R, "MMC_bounds")
    if (k0 < 2){
      stop("MMC_bounds: box_boundary = \"union\" requires k0 >= 2 (a one-regime null has no transition matrix, so the boundary face is undefined).")
    }
    bb <- .mmc_boundary_blocks(mdl_h0)
    # The block split assumes canonical regime order (absorbing = largest
    # fitted persistence = LAST). MMCLRTest canonicalizes before calling here;
    # a direct caller on a raw fit would otherwise silently stretch the wrong
    # column as "absorbing".
    pjj_chk <- diag(matrix(c(theta_0)[as.numeric(mdl_h0$theta_P_ind) == 1], k0, k0))
    if (max(pjj_chk) > pjj_chk[k0]){
      stop("MMC_bounds: box_boundary = \"union\" requires the fit in canonical regime order ",
           "(absorbing regime last -- largest fitted persistence in the final position); ",
           "apply canonicalize_regimes() first.")
    }
    R  <- .mmc_R_ranges(mdl_h0, con$boundary_R)
    theta_low[bb$rough_mu] <- pmin(theta_low[bb$rough_mu], R$mu_low)
    theta_upp[bb$rough_mu] <- pmax(theta_upp[bb$rough_mu], R$mu_upp)
    sig_low_R <- if (is.finite(R$sig_low)) rep(R$sig_low, length(bb$rough_sig_diag))
                 else con$variance_constraint * c(theta_0)[bb$rough_sig_diag]
    theta_low[bb$rough_sig_diag] <- pmin(theta_low[bb$rough_sig_diag], sig_low_R)
    theta_upp[bb$rough_sig_diag] <- pmax(theta_upp[bb$rough_sig_diag], R$sig_upp)
    # Covariance off-diagonals: +/- sqrt(U_i * U_j) from the coupled
    # variances' REALIZED upper ends, resolved within each entry's own regime
    # segment, so any positive semi-definite matrix whose diagonals fit the
    # box has its off-diagonals inside it (Cauchy-Schwarz); built after the
    # diagonal ranges are final so per-coordinate widening propagates.
    if (length(bb$rough_sig_off) > 0){
      off_bound <- sqrt(theta_upp[bb$off_diag_i] * theta_upp[bb$off_diag_j])
      theta_low[bb$rough_sig_off] <- pmin(theta_low[bb$rough_sig_off], -off_bound)
      theta_upp[bb$rough_sig_off] <- pmax(theta_upp[bb$rough_sig_off],  off_bound)
    }
    theta_low[bb$trans_P_idx] <- pmin(theta_low[bb$trans_P_idx], R$p_low)
    theta_upp[bb$trans_P_idx] <- pmax(theta_upp[bb$trans_P_idx], R$p_upp)
    theta_upp[bb$abs_diag_idx] <- 1
    theta_low[bb$abs_off_idx]  <- 0
    # Loud post-assertion: the final box must contain R on every rough
    # coordinate and reach the face on the absorbing column.
    reach_ok <- all(theta_low[bb$rough_mu] <= R$mu_low, theta_upp[bb$rough_mu] >= R$mu_upp,
                    theta_low[bb$rough_sig_diag] <= sig_low_R, theta_upp[bb$rough_sig_diag] >= R$sig_upp,
                    theta_low[bb$trans_P_idx] <= R$p_low, theta_upp[bb$trans_P_idx] >= R$p_upp,
                    theta_upp[bb$abs_diag_idx] >= 1, theta_low[bb$abs_off_idx] <= 0)
    if (!isTRUE(reach_ok)){
      stop("MMC_bounds: the boundary-stretched box failed to cover its own R ranges or the face; ",
           "this indicates contradictory box controls.")
    }
  }
  # 'box_include' containment check: if an admissible-region clamp removed a
  # requested value, the caller's controls are contradictory -- error naming
  # the coordinate rather than silently searching a box that excludes the
  # point the caller asked for.
  if (!is.null(box_inc)){
    inc_idx <- which(!is.na(box_inc))
    out_inc <- inc_idx[box_inc[inc_idx] < c(theta_low)[inc_idx] - 1e-12 |
                       box_inc[inc_idx] > c(theta_upp)[inc_idx] + 1e-12]
    if (length(out_inc) > 0){
      nm <- names(theta_0)[out_inc]
      if (is.null(nm)) nm <- paste("position", out_inc)
      stop("MMC_bounds: 'box_include' asks for ",
           paste(sprintf("%s = %.6g", nm, box_inc[out_inc]), collapse = ", "),
           " but an admissible-region control (variance_constraint / P_low / P_upp / ",
           "phi_low / phi_upp) excludes it from the box. The controls are contradictory; ",
           "relax the clamp or drop the box_include entry.")
    }
  }
  # ----- output
  # theta_low, theta_upp and theta_0 should all be finite by this point (eps and
  # every box-shaping control are validated above, and theta_0 comes from an
  # already-fitted model) -- but a pathological fit can still hand this function a
  # non-finite theta_0 having only triggered a warning upstream (the 'silent NaN
  # theta' check in the model constructors), and `if (any(NA > x))` throws R's own
  # generic "missing value where TRUE/FALSE needed" instead of an informative
  # message. Caught explicitly, once, so the failure names what's actually wrong;
  # the two checks below can then assume everything reaching them is finite.
  if (any(!is.finite(theta_low)) || any(!is.finite(theta_upp)) || any(!is.finite(theta_0))){
    bad <- which(!is.finite(theta_low) | !is.finite(theta_upp) | !is.finite(theta_0))
    nm  <- names(theta_0)[bad]
    if (is.null(nm)) nm <- paste("position", bad)
    stop("MMC_bounds: the search box is non-finite for ", paste(nm, collapse = ", "),
         " -- theta_0, theta_low, or theta_upp has a non-finite entry there. This ",
         "usually means the null model's own fit produced a non-finite parameter estimate.")
  }
  # A lower bound above its upper bound is unusable by any of the three optimizers:
  # GenSA errors, GA dies on a degenerate population, and pso silently searches the
  # inverted interval and returns a p-value for it. Reachable from 'eps' alone, and
  # more easily from an asymmetric 'phi_low'/'phi_upp'/'P_low'/'P_upp' clamp (only
  # one side of the pair is adjusted above), so this is caught here, once, right
  # after the box is built and before any simulation or search.
  if (any(theta_low > theta_upp)){
    bad <- which(theta_low > theta_upp)
    nm  <- names(theta_0)[bad]
    if (is.null(nm)) nm <- paste("position", bad)
    stop("MMC_bounds: the search box has a lower bound above its upper bound for ",
         paste(sprintf("%s (%.6g > %.6g)", nm, theta_low[bad], theta_upp[bad]), collapse = ", "),
         ". Widen 'eps', or check 'phi_low'/'phi_upp'/'P_low'/'P_upp'.")
  }
  # A fitted parameter outside its own search box means the observed fit is not a
  # point the search can return, so the reported p-value would not be a supremum
  # over a region containing theta_0 (and pso, whose 'par' seed must lie in the box,
  # would silently drop it). Reachable through 'phi_low'/'phi_upp' or
  # 'P_low'/'P_upp' set narrower than the fitted value, not just through 'eps'.
  if (any(theta_0 < theta_low | theta_0 > theta_upp)){
    bad <- which(theta_0 < theta_low | theta_0 > theta_upp)
    nm  <- names(theta_0)[bad]
    if (is.null(nm)) nm <- paste("position", bad)
    stop("MMC_bounds: the fitted parameter vector falls outside the search box for ",
         paste(sprintf("%s (%.6g, box [%.6g, %.6g])", nm, theta_0[bad], theta_low[bad], theta_upp[bad]), collapse = ", "),
         ". Widen 'eps', or check 'phi_low'/'phi_upp'/'P_low'/'P_upp'.")
  }
  mmc_bounds <- list(theta_low = as.numeric(theta_low), theta_upp = as.numeric(theta_upp))
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
#'   \item mc_seed: Integer seed for reproducible Monte Carlo simulations and fixed-error MMC (Dufour 2006, Prop. 4.2). When set, seeds the RNG before observed-data estimation, pre-draws innovations once and holds them fixed across all optimizer evaluations, ensuring full reproducibility (including the observed LRT statistic) and theoretical size control. Internal seeds (optimizer trajectory, worker streams, CRN scheme) are drawn from the \code{mc_seed} stream rather than derived by arithmetic, so consecutive \code{mc_seed} values across replications of an outer simulation study do not produce colliding RNG streams; for studies with more than about 10^4 replications, well-separated seeds are still recommended. Default is \code{NULL}, in which case internal seeds are drawn from the ambient RNG state, so a script-level \code{set.seed()} placed before the call makes the procedure
#'     reproducible as well (sequential and parallel). Note: GenSA is deterministic given the fixed starting value used here, so its search trajectory does not depend on any seed; pso and GA use the optimizer seed.
#'   \item init_method: String determining how the alternative (\code{k1}-regime) model is initialized. \code{"warmstart"} (the default) adds deterministic warm starts built by embedding the null-model fit into the \code{k1}-regime space (see \code{\link{warmstart_theta}}) to the \code{use_diff_init} random starts, applied identically to the observed data and to every simulated draw. \code{"random"} reproduces the legacy behavior.
#'   \item crn: Boolean determining whether common random numbers are used across the nuisance-parameter search (default \code{TRUE}): an identical estimation-RNG stream (EM starting values; master and worker streams) is replayed at every candidate theta evaluation, so that together with the pre-drawn innovations the MC p-value is a deterministic function of theta (Dufour 2006, Prop. 4.2, with the estimation randomness folded into the fixed disturbance vector). This removes optimizer objective noise and guarantees \code{pval >= pval_0} exactly. Set to \code{FALSE} for the legacy behavior (estimation randomness evolves across evaluations; still valid, but the objective is noisy and the reported maximum tends to be conservative). Note: with \code{init_method = "random"}, \code{crn = FALSE}, \code{workers = 0} and no \code{mc_seed}, a script-seeded call reproduces version 0.1.9 bit-for-bit.
#'   \item type: String that determines the type of optimization algorithm used. Arguments allowed are: \code{"pso"}, \code{"GenSA"}, and \code{"GA"}. Default is \code{"pso"}.
#'   \item eps: Double, or a vector of length \code{length(mdl_h0$theta)}, determining the half-width of the search box around \code{theta_0} for each parameter (a scalar is recycled to every coordinate; a vector lets the search be narrower for some parameters, e.g. transition probabilities, and wider for others, e.g. means -- entries are matched to \code{names(mdl_h0$theta)} positionally, i.e. in \code{theta}'s own order, so mind that a transition-probability entry named \code{p_12} is matrix entry \code{P[2,1]}). Every entry must be finite and non-negative; a variance coordinate whose lower bound would reach zero is overridden by \code{variance_constraint} regardless of \code{eps}. \code{eps = 0} is valid and typically relies on \code{CI_union} to widen the box; with \code{type = "GenSA"} a zero entry is nudged to \code{1e-8} (GenSA errors on a genuinely zero-width bound; \code{pso}/\code{GA} do not). Default is \code{0.1}.
#'   \item CI_union: Boolean determining if union of set determined by \code{eps} and confidence set should be used to define consistent set for search. Default is \code{TRUE}.
#'   \item lambda: Double determining penalty on nonlinear constraint. Default is \code{100}.
#'   \item stationary_constraint: Boolean determining if only stationary solutions are considered (if \code{TRUE}) or not (if \code{FALSE}). Default is \code{TRUE}.
#'   \item phi_low: Vector with lower bound for autoregressive parameters when optimizing. Default is \code{NULL}.
#'   \item phi_upp: Vector with upper bound for autoregressive parameters when optimizing. Default is \code{NULL}.
#'   \item P_low: Value with lower bound for transition probabilities when optimizing. Default is \code{0}.
#'   \item P_upp: Value with upper bound for transition probabilities when optimizing. Default is \code{1}.
#'   \item variance_constraint: Double used to determine the lower bound for variance in parameter set for search. Value should be between \code{0} and \code{1} as it is multiplied by consistent point estimates of variances. Default is \code{0.01} (i.e., \code{1\%} of consistent point estimates.
#'   \item box_include: Optional numeric vector of length \code{length(mdl_h0$theta)} giving, per parameter, a value the search box must contain (\code{NA} entries are left untouched). Applied after the \code{eps}/\code{CI_union} construction, so it only ever widens the box, and only where the ordinary construction does not already reach the value. The widening is TWO-SIDED: the coordinate's interval becomes the convex hull of the ordinary interval and the requested value, \code{[min(lower, v), max(upper, v)]}, so the value is covered even when the fitted estimate lies entirely on its far side -- under a true unit root the fitted autoregressive coefficient exceeds 1 in a non-negligible fraction of samples, and raising only the upper endpoint would then leave the unit root outside the box. The admissible-region controls (\code{variance_constraint}, \code{P_low}/\code{P_upp},
#'     \code{phi_low}/\code{phi_upp}) keep their usual meaning; if one of them excludes a requested value the call errors (the controls are contradictory) rather than silently searching a box that misses the requested point. The main use is unit-root testing: with \code{stationary_constraint = FALSE} and \code{phi_upp} not binding below 1, setting the autoregressive coordinate(s) of \code{box_include} to the unit-root point guarantees the search set contains it in every sample, per Rodriguez-Rondon & Dufour (2026a); for \code{p > 1} pass the coefficients obtained by moving the long-run coefficient \code{sum(phi)} to 1 while holding the short-run (augmented Dickey-Fuller) coefficients fixed, as computed e.g. by \code{\link{aURfixLMCTest}}'s reported \code{theta_fix}. In this configuration keep the box tight above 1 (a small \code{eps}, e.g. \code{0.02}, or
#'     \code{phi_upp} slightly above 1, e.g. \code{1.02}): a box reaching far into the explosive region admits candidates whose simulated series overflow, so every one of their draws fails estimation and the candidate is only penalized after the full per-draw retry budget -- valid, but a large waste of computation per such candidate. Default is \code{NULL}.
#'   \item box_boundary: Optional string enabling the boundary-setting search box for a null with \code{k0 >= 2} regimes, one of which may be absorbing. Requires \code{burnin = 0}: the boundary family simulates every candidate transient-start with no burn-in, and the regimes are relabeled canonically first via \code{\link{canonicalize_regimes}} so the absorbing regime is last. The family's initial condition is part of its definition: every candidate's chain is initialized in \strong{canonical state 1 of the fitted null} -- the least persistent transient regime at the fit -- and over the search box that initial state is inherited by index (candidates are never re-sorted, even where a candidate's persistence ordering differs from the fit's). With \code{k0 = 2} the transient regime is unique, so this pins the initial state with no choice made; with \code{k0 >= 3} it is a fixed-vertex
#'     family definition -- a maintained modeling statement about the data, on the same footing as transient-start itself -- and the procedure's finite-sample validity is conditional on it. At \code{p = 0} the chain is initialized in the transient regime at time 0 (the first observation's state is a draw from that state's transition distribution). Multivariate models are supported -- the transient block's ranges are built per series, and since a per-coordinate box cannot express positive definiteness, the search runs over the box intersected with the admissible parameter space (a candidate whose assembled regime covariance is not positive definite is excluded before simulation; see \code{pd_penalty_hits} in the output). Exogenous regressors are supported: their common coefficients belong to the clean block (identified by the whole sample, like the autoregressive
#'     coefficients), so they keep their ordinary range in \code{"union"} mode and are held at their estimates in \code{"face"} mode; \code{Z} must be fixed (deterministic) or strictly exogenous regressors -- dummy variables are deterministic and fully supported. One scope note for the clean-block ranges: the consistency argument behind them attaches per coefficient, to coefficients whose regressor's information grows with the sample; a coefficient identified only by a fixed observation window (e.g. a short-lived dummy) is instead covered through the maintained containment of its range, so widen that coordinate's range if it matters. \code{"union"} stretches the ordinary search box per coordinate: the transient (rough) block's ranges are widened to the moment-anchored ranges of \code{boundary_R}, and the absorbing regime's transition-matrix column is stretched so the
#'     boundary face (its diagonal equal to 1) is inside the box; one ordinary search runs over the stretched box. \code{"face"} searches the boundary face alone: the transient block over the \code{boundary_R} ranges with everything else fixed (the absorbing diagonal at exactly 1, the remaining parameters at their estimates); the reported \code{pval_proj} is then the evaluation at the (unclamped) boundary projection of the fitted parameters, and the search seed is that projection clamped into the box (for a multivariate model the clamped seed's transient covariance is repaired to positive definiteness by deterministic off-diagonal shrinkage when needed, reported in \code{seed_repair}). Default is \code{NULL} (off; the ordinary search).
#'   \item boundary_R: List of moment-anchored ranges for the transient block when \code{box_boundary} is used: \code{mu_mult} (default \code{5}; each transient mean searched over its own series' \code{mean(y_i) +/- mu_mult * sd(y_i)}), \code{sigma_mult} (default \code{4}; each transient variance up to \code{sigma_mult * var(y_i)} and down to the EM regime-variance floor; transient covariances range over \code{+/- sqrt(U_i * U_j)} with \code{U_i} the realized variance upper ends, so any positive semi-definite covariance matrix whose variances fit the box has its covariances inside it automatically), and \code{p_range} (default \code{c(0.05, 0.95)}; transient-column transition probabilities). These ranges do the theory's coverage work for the parameters the pre-absorption era identifies, so a too-narrow choice loosens the size control; widen rather than narrow them.
#'   \item box_match_constraints: Boolean. If \code{TRUE} (the default), the nuisance search set is aligned with the estimation constraints carried by \code{mdl_h0_control}: with a stationarity margin (\code{stationarity_constraint = delta}) the search stays inside the margin region -- the autoregressive box is clamped to \code{[-(1 - delta), 1 - delta]} at \code{p = 1}, and for any \code{p} a candidate whose companion spectral radius exceeds \code{1 - delta} is excluded before simulation (counted in \code{margin_penalty_hits}; the margin itself is defined for univariate models only); with a transition-probability bound (\code{em_transition_constraint = eps}, meaningful for \code{k0 > 1} and for EM nulls -- the bound is not part of the MLE path, so it is not aligned when \code{mdl_h0_control} sets \code{method = "MLE"}) \code{P_low} is raised to \code{eps} and \code{P_upp} lowered to \code{1 -
#'     (k0 - 1) * eps}. A \code{box_include} autoregressive point outside the margin is contradictory and errors. The constrained fits obey these constraints on the observed and simulated sides, so this keeps the maximized test's search set equal to the maintained null parameter space, and the local test remains the search's seed point. It has no effect when no estimation constraint is set. A stationarity margin together with \code{stationary_constraint = FALSE} is contradictory and errors. If \code{FALSE}, the box is built as usual and a warning notes that the search set may extend beyond the constrained space. Default is \code{TRUE}.
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
#'   \item mmc_optimout: Optimization output object returned by the selected optimizer (\code{pso}, \code{GenSA}, or \code{GA}). For \code{k0 >= 2} the search itself runs over a reduced parameter vector (the transition matrix contributes only its free entries), and the optimizer's reported solution (\code{$par} for \code{pso}/\code{GenSA}, \code{@solution} for \code{GA}) is mapped back to the full \code{theta} coordinates used everywhere else in the return value; other diagnostic fields inside \code{mmc_optimout} (e.g. \code{pso}'s \code{trace.stats$x}, \code{GenSA}'s \code{trace.mat}) are left in the reduced search coordinates.
#'   \item pval_0: Monte Carlo p-value at the initial (estimated) nuisance parameter values \code{theta_0} -- the Local MC point of the search. Present for every estimate-centered search; the reported \code{pval} satisfies \code{pval >= pval_0} by construction (enforced directly, not merely relied on from the optimizer's own search). With \code{box_boundary = "face"} the fitted point lies outside the searched set, so this field is replaced by \code{pval_proj}: the evaluation at the unclamped boundary projection of the fitted parameters (reported as \code{theta_face_proj}), with \code{pval >= pval_proj} enforced the same way.
#'   \item theta_low, theta_upp: The realized search box (full parameter coordinates, named) after every widening and clamp -- stored so it can be verified from saved output whether the box contained a target value.
#'   \item optim_diagnostics: List with \code{type} (the optimizer, or \code{"early_stop"}), \code{convergence} (the optimizer's own exit code where it reports one), \code{feval} (exact number of objective evaluations, counted directly), \code{optimizer_feval} (the optimizer's own reported call count where one exists, as a cross-check), \code{pd_penalty_hits} (number of candidates excluded because an assembled regime covariance was not positive definite; \code{NA} when that gate is not active, i.e. for univariate models or when neither \code{box_boundary} nor \code{box_include} is in use), \code{margin_penalty_hits} (number of candidates excluded by the stationarity-margin gate of \code{box_match_constraints}; \code{NA} when no margin is in force), and \code{evaluated_range} (per-coordinate min/max of the candidates the optimizer evaluated, in the reduced search coordinates, named by the searched positions; \code{NULL} when the search was skipped by \code{threshold_stop}).
#'   \item boundary_mode, absorbing_regime, theta_face_proj, seed_repair: Present only when \code{box_boundary} is used -- the mode, the absorbing regime's index after canonical labeling (always \code{k0}), and (face mode) the unclamped boundary projection whose evaluation is \code{pval_proj} together with \code{seed_repair} (whether the clamped search seed's transient covariance needed the positive-definiteness repair, and the off-diagonal shrinkage factor applied). With \code{box_boundary} set, \code{mdl_h0} is returned in canonical regime order (see \code{\link{canonicalize_regimes}}).
#' }
#'
#' @section Numerical policy for the LR statistic: For nested models the LR statistic is non-negative by construction (the null fit embeds into the alternative parameter space with equal likelihood), so a negative computed value can only be a numerical optimization artifact. Both the observed statistic and every simulated statistic are therefore floored at 0 (identically on both sides, which preserves the exchangeability underlying the Monte Carlo p-value); a warning reports the magnitude when it is non-negligible. For \code{k0 > 1}, the null log-likelihood is additionally floored at the one-regime (linear) fit on both sides. Ties at 0 are handled by the randomized tie-breaking rule of Dufour (2006) already implemented in \code{\link{MCpval}}; ties at the bottom of the support can only make the test conservative. The regime variances are bounded away from zero
#'   during estimation (control \code{em_variance_constraint} of the model constructors, default \code{0.01} of the one-regime fit variance), because the Markov-switching likelihood is otherwise unbounded and degenerate solutions produce spurious likelihood-ratio spikes; the statistic is therefore a constrained-likelihood ratio, computed by the same rule for the observed and the simulated samples, which is all the Monte Carlo argument requires. Transition probabilities are not constrained.
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
              box_include = NULL,
              box_boundary = NULL,
              boundary_R = NULL,
              box_match_constraints = TRUE,
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
  if (k1 <= k0){
    stop("MMCLRTest: k1 must be greater than k0 (received k0 = ", k0, ", k1 = ", k1, ").")
  }
  if ((con$CI_union==TRUE) & ((con$mdl_h0_control$getSE==FALSE) | (con$mdl_h1_control$getSE==FALSE))){
    con$mdl_h0_control$getSE <- TRUE
    con$mdl_h1_control$getSE <- TRUE
    warning("getSE was changed to be 'TRUE' because CI_union is 'TRUE'.")
  }
  # ----- Align the nuisance search set with the estimation constraints
  # (control box_match_constraints). An estimation constraint in the model
  # controls (a stationarity margin on the autoregressive block, a bound on
  # the transition probabilities) defines the maintained null parameter space;
  # the constrained fits obey it on the observed and simulated sides alike, so
  # by default the search runs over that same space. Where the constraint set
  # is a box it is imposed through the existing clamps (phi_low/phi_upp at
  # p = 1, P_low/P_upp); the margin for p >= 2 is a set on the companion
  # spectral radius, not a box, and is imposed as a penalty in the objective
  # wrapper below, next to the plain stationarity penalty.
  if (!is.logical(con$box_match_constraints) || length(con$box_match_constraints) != 1L || is.na(con$box_match_constraints)){
    stop("MMCLRTest: 'box_match_constraints' must be TRUE or FALSE.")
  }
  sc_delta <- .validate_stationarity_constraint(con$mdl_h0_control[["stationarity_constraint"]], "MMCLRTest")
  P_eps <- con$mdl_h0_control[["em_transition_constraint"]]
  # The transition-probability bound is an EM-path constraint (MLE fits do not
  # impose it), so it is aligned for EM nulls only.
  h0_is_mle <- isTRUE(toupper(as.character(con$mdl_h0_control[["method"]])[1]) == "MLE")
  P_eps <- if (k0 > 1 && !h0_is_mle && is.numeric(P_eps) && length(P_eps) == 1L && isTRUE(P_eps > 0)) as.numeric(P_eps) else 0
  margin_active <- con$box_match_constraints && sc_delta > 0 && p > 0 && q == 1
  if (con$box_match_constraints){
    if (sc_delta > 0 && p > 0){
      if (!isTRUE(con$stationary_constraint)){
        stop("MMCLRTest: 'stationarity_constraint' in mdl_h0_control keeps every fit inside the stationarity margin, but 'stationary_constraint = FALSE' asks the search to leave the stationary region -- contradictory controls. Set box_match_constraints = FALSE to search beyond the estimation constraint set (the search set then no longer matches the maintained null space).")
      }
      if (p == 1 && q == 1){
        con$phi_upp <- if (is.null(con$phi_upp)) 1 - sc_delta else pmin(con$phi_upp, 1 - sc_delta)
        con$phi_low <- if (is.null(con$phi_low)) -(1 - sc_delta) else pmax(con$phi_low, -(1 - sc_delta))
      }
    }
    if (P_eps > 0){
      con$P_low <- max(con$P_low, P_eps)
      con$P_upp <- min(con$P_upp, 1 - (k0 - 1) * P_eps)
    }
  } else if ((sc_delta > 0 && p > 0) || P_eps > 0){
    warning("MMCLRTest: box_match_constraints = FALSE with an estimation constraint in mdl_h0_control (stationarity_constraint / em_transition_constraint): the nuisance search set may extend beyond the constrained parameter space that the fits obey.")
  }
  # ----- Boundary box toggle: scope checks before any estimation. The boundary
  # family convention (transient-start simulation, burnin 0, canonical regime
  # labeling) is a different maintained model of the data, so a call that
  # contradicts it errors rather than being silently adjusted.
  if (!is.null(con$box_boundary)){
    con$box_boundary <- match.arg(con$box_boundary, c("union", "face"))
    if (k0 < 2){
      stop("MMCLRTest: 'box_boundary' requires a null with k0 >= 2 regimes (a one-regime null has no transition matrix, so the boundary face is undefined).")
    }
    if (con$burnin != 0){
      stop("MMCLRTest: 'box_boundary' requires burnin = 0 -- the boundary family simulates every candidate transient-start (the chain initialized in the transient regime), and a burnin would consume the transient era and compare against the wrong process. Set control$burnin = 0 explicitly.")
    }
    if (con$P_upp < 1 || con$P_low > 0){
      stop("MMCLRTest: 'box_boundary' requires P_low = 0 and P_upp = 1 -- the stretched box must reach the absorbing face (p = 1 and 0 entries), and a tighter clamp contradicts it (an 'em_transition_constraint' in mdl_h0_control imposes such a clamp through box_match_constraints).")
    }
    if (identical(con$box_boundary, "face") && !is.null(con$box_include)){
      stop("MMCLRTest: 'box_include' is not meaningful with box_boundary = \"face\" -- the face set is defined by the boundary_R ranges and the pinned coordinates, not by the eps/CI box that box_include stretches.")
    }
    con$boundary_R <- .validate_boundary_R(con$boundary_R)
  }
  # ----- Seed RNG for reproducible estimation (random initial values in EM)
  if (!is.null(con$mc_seed)) set.seed(con$mc_seed)
  # ----- Estimate models using observed data
  mdl_h0 <- estimMdl(Y, p, q, k0, Z, con$mdl_h0_control)
  # ----- Boundary modes work in canonical regime coordinates: transient
  # regimes first, absorbing regime last (ascending fitted persistence, with
  # canonicalize_regimes' deterministic tie-break). Deterministic and RNG-free,
  # so the stream layout below is unchanged; everything downstream (warm
  # starts, bounds, search, simulation, the returned mdl_h0) uses this
  # labeling.
  if (!is.null(con$box_boundary)){
    mdl_h0 <- canonicalize_regimes(mdl_h0, by = "persistence")
  }
  # Resolve the alternative's switching flags before the warm starts are built and
  # before the alternative is estimated, so the two cannot disagree.
  con$mdl_h1_control <- .inherit_ms_flags(con$mdl_h1_control, mdl_h0)
  con$mdl_h1_control <- .inherit_stationarity_margin(con$mdl_h1_control, con$mdl_h0_control, "MMCLRTest", p)
  if ((con$init_method == "warmstart") && (k1 > 1)){
    # Warm-start the alternative model from the null fit (regime-duplication
    # embedding, see warmstart_theta); deterministic, so the RNG stream and the
    # random starts are unaffected. The same rule is applied per draw when
    # simulating the null distribution (see LR_samp_dist).
    # `[[` rather than `$`: the latter partial-matches, so a mistyped control name
    # would be read as the flag and silently build a warm start for the wrong space.
    con$mdl_h1_control$init_theta_extra <- warmstart_theta(mdl_h0, k1,
                                                           msmu = con$mdl_h1_control[["msmu"]],
                                                           msvar = con$mdl_h1_control[["msvar"]])
  }
  mdl_h1 <- estimMdl(Y, p, q, k1, Z, con$mdl_h1_control)
  # The flags are a proxy; what has to hold is that a warm start has the length the
  # alternative expects, which catches any future cause of a mismatch.
  ws_extra <- con$mdl_h1_control$init_theta_extra
  if (!is.null(ws_extra) && (length(ws_extra) > 0) &&
      any(vapply(ws_extra, length, 0L) != length(mdl_h1$theta))){
    warning("the warm starts did not match the alternative model's parameter vector ",
            "and were discarded.", call. = FALSE)
  }
  con$mdl_h0_control <- mdl_h0$control
  con$mdl_h1_control <- mdl_h1$control
  # The stationarity margin must be the same on both sides: with different
  # margins the alternative no longer nests the null's constrained space and
  # the LR statistic loses its constrained-ML meaning.
  d0 <- con[["mdl_h0_control"]][["stationarity_constraint"]]
  d1 <- con[["mdl_h1_control"]][["stationarity_constraint"]]
  if (!is.null(d0) && !is.null(d1) && !isTRUE(all.equal(as.numeric(d0), as.numeric(d1)))){
    warning("'stationarity_constraint' differs between mdl_h0_control (", d0, ") and mdl_h1_control (", d1, "); the same margin is required on both sides for the statistic to be a constrained-likelihood ratio.")
  }
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
    # The stationarity margin, when set on the null control, carries into the
    # floor fit: the floor's nesting rationale needs the floor point feasible
    # in the constrained space, by the same rule on both sides.
    mdl_lin_con <- list(getSE = FALSE)
    if (!is.null(con[["mdl_h0_control"]][["stationarity_constraint"]])){
      mdl_lin_con$stationarity_constraint <- con[["mdl_h0_control"]][["stationarity_constraint"]]
    }
    mdl_lin <- estimMdl(Y, p, q, 1, Z, mdl_lin_con)
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
  # ----- Define lower & upper bounds for search. "face" mode builds the face
  # set directly (rough block over R, everything else pinned) and does not go
  # through MMC_bounds' eps/CI construction; its in-box reference point is the
  # CLAMPED boundary projection, since theta_hat itself legitimately lies
  # outside the face (p_hat_AA < 1).
  is_face <- isTRUE(con$box_boundary == "face")
  if (is_face){
    fb <- .mmc_face_bounds(mdl_h0, con)
    theta_low <- fb$theta_low
    theta_upp <- fb$theta_upp
    pinned_idx  <- fb$pinned_idx
    pinned_vals <- fb$pinned_vals
  } else {
    mmc_bounds <- MMC_bounds(mdl_h0, con)
    theta_low <- mmc_bounds$theta_low
    theta_upp <- mmc_bounds$theta_upp
    fb <- NULL
    pinned_idx  <- integer(0)
    pinned_vals <- numeric(0)
  }
  # ----- Reduce the search vector for k0 >= 2 to the free entries of P (see the
  # block above MMC_bounds), minus any pinned coordinates (face mode: the clean
  # block and the absorbing column are constants and are dropped from the
  # optimizer's dimensions entirely -- no zero-width bounds reach GenSA). k0 = 1
  # has no P block; 'reduced_idx' is then every position, so the '_search'
  # vectors below equal theta_low/theta_upp/theta_0.
  if (k0 > 1){
    mmc_rp <- .mmc_derived_P_index(theta_0, mdl_h0$theta_P_ind, k0, fixed_full_idx = pinned_idx)
    derived_nonNA <- mmc_rp$derived_full_idx[!is.na(mmc_rp$derived_full_idx)]
    reduced_idx <- setdiff(seq_along(theta_0), c(derived_nonNA, pinned_idx))
  } else {
    mmc_rp <- NULL
    reduced_idx <- seq_along(theta_0)
  }
  # The search starts at the fitted theta_0 (in-box by MMC_bounds' pre-flight
  # checks), except in face mode, where it starts at the clamped projection
  # (in-box by construction; pso silently drops an out-of-box 'par').
  theta_start      <- if (is_face) fb$proj_clamped else theta_0
  theta_0_search   <- theta_start[reduced_idx]
  theta_low_search <- theta_low[reduced_idx]
  theta_upp_search <- theta_upp[reduced_idx]
  if (length(theta_0_search) != length(theta_low_search) ||
      length(theta_0_search) != length(theta_upp_search)){
    # Defensive: pso accepts mismatched par/lower/upper lengths silently
    # (recycling or truncating); GenSA and GA already error loudly on this.
    stop("MMCLRTest: internal error building the nuisance-parameter search vector -- ",
         "lengths differ (par ", length(theta_0_search), ", lower ", length(theta_low_search),
         ", upper ", length(theta_upp_search), ").")
  }
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
  # ----- The optimizer sees the reduced search vector, so its objective must
  # reassemble to full theta before calling the underlying (already CRN-wrapped)
  # objective; an infeasible derived entry is penalized without that call, so it
  # costs no random draws (see .mmc_reduced_wrap above). k0 = 1 has no P block,
  # so the objective is used unwrapped.
  # ----- Search bookkeeping environment: the candidate recorder (exact
  # evaluation counts and evaluated ranges for the harvest fields) and the
  # positive-definiteness exclusion counter both write here. Pure bookkeeping
  # -- no RNG, no numeric change.
  mmc_rec <- new.env(parent = emptyenv())
  mmc_rec$n <- 0L; mmc_rec$X <- list(); mmc_rec$v <- numeric(0); mmc_rec$pd_hits <- 0L; mmc_rec$margin_hits <- 0L
  # ----- Positive-definiteness gate for multivariate covariance candidates:
  # a per-coordinate box (the boundary stretch, or a box_include hull on vech
  # coordinates) cannot express positive definiteness, so the searched set is
  # box-intersect-admissible-space, enforced by the wrapper. Enabled only when
  # one of those box-shaping controls is in use: the plain search keeps its
  # long-standing behavior (an inadmissible candidate fails in simulation and
  # is penalized there).
  mmc_rec$pd_active <- q > 1 && (!is.null(con$box_boundary) || !is.null(con$box_include))
  obj_min_full <- obj_min
  obj_max_full <- obj_max
  if (mmc_rec$pd_active){
    sig_pos <- which(as.numeric(mdl_h0$theta_sig_ind) == 1)
    Nsig <- q * (q + 1) / 2
    psd_blocks <- split(sig_pos, rep(seq_len(length(sig_pos) / Nsig), each = Nsig))
    obj_min_full <- .mmc_admissible_wrap(obj_min, psd_blocks, q, con$lambda, mmc_rec)
    obj_max_full <- .mmc_admissible_wrap(obj_max, psd_blocks, q, -(con$lambda), mmc_rec)
  }
  mmc_rec$margin_active <- margin_active
  if (margin_active){
    phi_pos <- which(as.numeric(mdl_h0$theta_phi_ind) == 1)
    obj_min_full <- .mmc_margin_wrap(obj_min_full, phi_pos, 1 - sc_delta, con$lambda, mmc_rec)
    obj_max_full <- .mmc_margin_wrap(obj_max_full, phi_pos, 1 - sc_delta, -(con$lambda), mmc_rec)
  }
  if (k0 > 1){
    obj_min_search <- .mmc_reduced_wrap(obj_min_full, reduced_idx, mmc_rp$derived_full_idx,
                                        mmc_rp$pos_mat, length(theta_0), theta_low, theta_upp,
                                        con$P_low, con$P_upp, penalty_value = con$lambda,
                                        pinned_idx = pinned_idx, pinned_vals = pinned_vals)
    obj_max_search <- .mmc_reduced_wrap(obj_max_full, reduced_idx, mmc_rp$derived_full_idx,
                                        mmc_rp$pos_mat, length(theta_0), theta_low, theta_upp,
                                        con$P_low, con$P_upp, penalty_value = -(con$lambda),
                                        pinned_idx = pinned_idx, pinned_vals = pinned_vals)
  } else {
    obj_min_search <- obj_min_full
    obj_max_search <- obj_max_full
  }
  # ----- Candidate recorder for the harvest fields: runs unconditionally for
  # every optimizer, recording each evaluated candidate into mmc_rec above.
  obj_min_search <- .mmc_record_wrap(obj_min_search, mmc_rec)
  obj_max_search <- .mmc_record_wrap(obj_max_search, mmc_rec)
  # ----- Verify the seed theta_0 (= the LMC point) yields a computable null distribution.
  #       theta_0 is exempt from the degenerate-candidate penalty inside MMCLRpval_fun, so if its
  #       null is fully degenerate (all draws fail even after the buffer + retry safety) the p-value
  #       is non-finite and would crash/garble the optimizer. Stop here instead. Under crn = TRUE
  #       this evaluation uses the optimizer's own worker configuration and CRN stream, so it is
  #       bit-identical to the optimizer's evaluation at theta_0: the guard is exact and the
  #       reported MMC p-value satisfies pval >= pval_0 by construction. Done BEFORE the optimizer
  #       seed is set below so it does not disturb the reproducible search trajectory.
  p_seed_workers <- if (isTRUE(con$crn)) con$workers else 0L
  # In face mode the reference evaluation is at the UNCLAMPED boundary
  # projection (theta_hat with the absorbing column at exactly (1, 0s)): an
  # evaluation needs no box membership, and this point is the same one
  # aBfixLMCTest evaluates its second component at (the two p-values differ
  # through their estimation-stream positions, but the point never silently
  # moves to the clamped one). Bypasses the reduced/recorded wrappers by
  # design (it is not a search evaluation).
  p_seed_point <- if (is_face) fb$proj_unclamped else theta_0
  p_seed <- obj_max(p_seed_point, mdl_h0, k1, LRT_0, con$N, con$burnin, p_seed_workers, con$lambda,
                    con$stationary_constraint, mdl_h1$control$thtol, Zsim, exog,
                    mdl_h0_null_cont, mdl_h1_null_cont, predrawn_eps, predrawn_state_rand)
  if (is_face && is.finite(p_seed) && p_seed < 0){
    # The projection is not the exempt seed inside MMCLRpval_fun, so an
    # unsimulable projection returns the raw penalty (-lambda). A classed
    # condition lets an outer simulation harness catch this and RECORD the
    # replication rather than redraw it -- redrawing on this error would
    # condition the sample on {projection simulable}, a silent selection
    # effect.
    stop(errorCondition(
      "MMCLRTest: the evaluation at the boundary face projection was penalized (its null distribution failed even after the re-draw safety, or it violated a constraint), so the face-mode p-value is not defined for this fit. Record this outcome rather than redrawing the sample.",
      class = "MSTest_face_unsimulable"))
  }
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
  # Tracked separately from mmc_out, whose type differs by branch (a plain list
  # here and for pso/GenSA, an S4 'ga' object for GA, which does not support $).
  is_early_stop <- FALSE
  if (p_seed >= con$threshold_stop){
    # ----- Early stop: the p-value at the seed evaluation point (theta_0, or the
    # boundary face projection in face mode) already meets threshold_stop, so the
    # search cannot change the test decision (pval >= p_seed >= threshold_stop by
    # construction). Skip the optimizer entirely -- exactly the threshold-stopping
    # semantics, at the cost of a single (LMC-type) evaluation. Note the reported
    # pval is then the seed p-value itself (an optimizer run crossing the threshold
    # could report a weakly larger value; both are valid and imply the same decision).
    # The reported theta must be the point the p-value belongs to: in face mode
    # that is the unclamped projection p_seed was evaluated at, not theta_0.
    theta   <- p_seed_point
    pval    <- p_seed
    seed_desc <- if (is_face) "the boundary face projection (pval_proj >= threshold_stop)" else "theta_0 (pval_0 >= threshold_stop)"
    mmc_out <- list(early_stop = TRUE,
                    message = paste0("threshold_stop reached at ", seed_desc, "; nuisance-parameter search skipped."))
    is_early_stop <- TRUE
    if (!isTRUE(con$silence)) message("MMCLRTest: threshold_stop reached at ", if (is_face) "the boundary face projection" else "theta_0", "; skipping the nuisance-parameter search.")
  }else if (con$type=="pso"){
    # Set PSO specific controls
    con$optim_control$trace.stats <- TRUE
    con$optim_control$trace <- as.numeric(con$silence==FALSE)
    con$optim_control$abstol <- -con$threshold_stop
    con$optim_control$maxf <- con$maxit
    con$optim_control$REPORT <- 1
    # begin optimization
    mmc_out   <- pso::psoptim(par = theta_0_search, fn = obj_min_search, lower = theta_low_search, upper = theta_upp_search,
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
    mmc_out   <- GenSA::GenSA(par = theta_0_search, fn = obj_min_search, lower = theta_low_search, upper = theta_upp_search,
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
    ga_fixed <- list(type = "real-valued", fitness = obj_max_search,
                     mdl_h0 = mdl_h0, k1 = k1, LRT_0 = LRT_0, N = con$N, burnin = con$burnin, workers = con$workers,
                     lambda = con$lambda, stationary_constraint = con$stationary_constraint,
                     thtol = mdl_h1$control$thtol, Z = Zsim, exog = exog, mdl_h0_control = mdl_h0_null_cont,
                     mdl_h1_control = mdl_h1_null_cont,
                     predrawn_eps = predrawn_eps, predrawn_state_rand = predrawn_state_rand,
                     lower = theta_low_search, upper = theta_upp_search,
                     maxiter = con$maxit, maxFitness = con$threshold_stop,
                     monitor = (con$silence==FALSE),
                     suggestions = if (all(theta_0_search >= theta_low_search) && all(theta_0_search <= theta_upp_search)) t(theta_0_search) else NULL)
    ga_ctrl   <- ga_ctrl[setdiff(names(ga_ctrl), names(ga_fixed))]
    mmc_out   <- do.call(GA::ga, c(ga_fixed, ga_ctrl))
    theta     <- as.matrix(mmc_out@solution[1,])
    pval      <- mmc_out@fitnessValue
  }else{
    stop("'type' must be one of 'pso', 'GenSA', or 'GA'.")
  }
  # ----- Map the optimizer's answer back to full theta coordinates. The
  # optimizer itself (mmc_out$par / @solution) still reports in the reduced
  # search space for k0 >= 2; iteration-level diagnostics inside mmc_out (e.g.
  # pso's trace.stats$x, GenSA's trace.mat) are not remapped and remain in
  # reduced coordinates, documented in the return value below.
  if (k0 > 1 && !is_early_stop){
    theta <- .mmc_reassemble_theta(theta, reduced_idx, mmc_rp$derived_full_idx, mmc_rp$pos_mat, length(theta_0),
                                   pinned_idx, pinned_vals)
    if (con$type %in% c("pso", "GenSA")) mmc_out$par <- theta
    if (con$type == "GA") mmc_out@solution <- matrix(theta, nrow = 1)
  } else {
    theta <- as.numeric(theta)
  }
  # ----- get test output using optimization output params
  # The search's starting point is always inside the box (theta_0 by
  # MMC_bounds' pre-flight checks; in face mode the clamped projection by the
  # face builder's own assertion), so the search always evaluates at least
  # that point (pso and GenSA both evaluate their 'par' starting value; GA is
  # seeded there when in-box, which it now always is). Under CRN the reported
  # pval is therefore never below p_seed by construction (in face mode p_seed
  # is the unclamped projection's evaluation, an extra member of the evaluated
  # set); enforced explicitly rather than relied on, since an optimizer's
  # internal bookkeeping is not a contract.
  # In face mode the evaluated set is the face box plus the unclamped
  # projection, and the reported theta must be the point the reported p-value
  # belongs to: when the projection's evaluation strictly beats the search, the
  # projection is the maximizer (ties go to the search point).
  if (is_face && !is_early_stop && p_seed > pval){
    theta <- p_seed_point
  }
  pval <- max(pval, p_seed)
  if (pval < 0){
    stop("MMCLRTest: the nuisance-parameter search evaluated no admissible candidate (every candidate violated the transition-matrix or stationarity constraint, or could not be simulated), so no p-value is defined. Check that the search bounds admit the null model's parameter values.")
  }
  # The reported p-value must belong to the reported parameter vector: an optimizer
  # can return a maximizer that was never the point its reported value came from.
  if (!all(is.finite(theta))){
    stop("MMCLRTest: the nuisance-parameter search returned a non-finite parameter vector, so the reported p-value would not correspond to it.")
  }
  if (k0>1){
    P_returned <- matrix(as.numeric(theta)[as.numeric(mdl_h0$theta_P_ind)==1], k0, k0)
    if (any(abs(colSums(P_returned) - 1) > min(mdl_h1$control$thtol, 1e-8))){
      stop("MMCLRTest: the nuisance-parameter search returned a transition matrix whose columns do not sum to 1, so the reported p-value would not correspond to it.")
    }
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
  # ----- Search diagnostics (harvest): the realized box, an exact evaluation
  # count, and the per-coordinate range of evaluated candidates (reduced search
  # coordinates, named by the searched positions), from the recorder -- exact
  # for every optimizer and independent of 'silence'.
  search_names <- names(mdl_h0$theta)[reduced_idx]
  eval_range <- NULL
  if (mmc_rec$n > 0){
    Xm <- do.call(rbind, mmc_rec$X)
    eval_range <- rbind(min = apply(Xm, 2, min), max = apply(Xm, 2, max))
    colnames(eval_range) <- search_names
  }
  optim_convergence <- if (is_early_stop) NA_integer_
    else if (con$type == "pso") mmc_out$convergence
    else if (con$type == "GA") as.integer(mmc_out@iter)
    else NA_integer_
  # The optimizer's own reported call count, where one exists, as a
  # cross-check against the recorder's exact count (GA reports iterations,
  # not calls, so it has none).
  optimizer_feval <- if (is_early_stop) 0L
    else if (con$type == "pso") as.integer(mmc_out$counts[["function"]])
    else if (con$type == "GenSA") as.integer(mmc_out$counts)
    else NA_integer_
  optim_diagnostics <- list(type = if (is_early_stop) "early_stop" else con$type,
                            convergence = optim_convergence,
                            feval = mmc_rec$n,
                            optimizer_feval = optimizer_feval,
                            pd_penalty_hits = if (isTRUE(mmc_rec$pd_active)) mmc_rec$pd_hits else NA_integer_,
                            margin_penalty_hits = if (isTRUE(mmc_rec$margin_active)) mmc_rec$margin_hits else NA_integer_,
                            evaluated_range = eval_range)
  names(theta_low) <- names(mdl_h0$theta)
  names(theta_upp) <- names(mdl_h0$theta)
  # ----- organize test output
  MMCLRTest_output <- list(mdl_h0 = mdl_h0, mdl_h1 = mdl_h1, mdl_h0_mmc = mdl_h0_mmc, mdl_h1_mmc = mdl_h1,
                           LRT_0 = LRT_0, pval = pval,
                           theta_h0 = theta_h0, theta_h1 = theta_h1, control = con,
                           mmc_optimout = mmc_out,
                           theta_low = theta_low, theta_upp = theta_upp,
                           optim_diagnostics = optim_diagnostics)
  # The reference evaluation is reported under the name that says what it is:
  # pval_0 is the Local MC p-value at the fitted theta_0, which exists for
  # every estimate-centered search; the face set does not contain theta_0, so
  # its reference evaluation -- at the boundary projection -- is pval_proj.
  if (is_face){
    MMCLRTest_output$pval_proj <- p_seed
  } else {
    MMCLRTest_output$pval_0 <- p_seed
  }
  if (!is.null(con$box_boundary)){
    MMCLRTest_output$boundary_mode <- con$box_boundary
    MMCLRTest_output$absorbing_regime <- k0
    if (is_face){
      MMCLRTest_output$theta_face_proj <- fb$proj_unclamped
      MMCLRTest_output$seed_repair <- fb$seed_repair
    }
  }
  class(MMCLRTest_output) <- "MMCLRTest"
  return(MMCLRTest_output)
}



