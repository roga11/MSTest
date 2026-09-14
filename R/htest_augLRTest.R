# Augmented Monte Carlo LR tests for the unit-root and boundary settings.
#
# Each test is an intersection-union composition (Berger 1982): two component
# p-values, reject only if BOTH are <= alpha, reported as pval = max(p1, p2).
# The composition is valid at any null point where at least one component is
# valid; no independence between the components is required, so they share
# the same pre-drawn simulation innovations and the same estimation-stream
# randomness ("only the parameter changes").

# Map an AR(p) coefficient vector to its unit-root counterpart in augmented
# Dickey-Fuller coordinates: the long-run coefficient rho = sum(phi) is moved
# to 1 while the short-run coefficients gamma_j = -sum(phi[(j+1):p]) are held
# fixed. Linear arithmetic on the fitted vector -- nothing is re-estimated.
# At p = 1 this is phi -> 1; for every p, sum of the result is exactly 1.
.adf_unit_root <- function(phi){
  p <- length(phi)
  if (p == 0) stop("no autoregressive coefficients to fix at the unit root.")
  if (p == 1) return(1)
  gam <- -rev(cumsum(rev(phi)))[-1]          # gamma_j = -sum_{i>j} phi_i, j = 1..p-1
  phi_new <- numeric(p)
  phi_new[1] <- 1 + gam[1]
  if (p > 2) phi_new[2:(p - 1)] <- gam[2:(p - 1)] - gam[1:(p - 2)]
  phi_new[p] <- -gam[p - 1]
  phi_new
}

# Strip the levels-specific search-box controls from a control list before it
# is reused for the DIFFERENCED component of an augmented test: 'box_include'
# (and a vector-valued 'eps') have the levels parameter vector's length, and
# the differenced model's is shorter; 'phi_low'/'phi_upp' and the boundary
# controls refer to the levels parameterization; and the differenced component
# is the stationary world by construction, so the stationarity constraint is
# restored to its default.
.aug_diff_control <- function(control){
  cd <- control
  cd$box_include <- NULL
  cd$phi_low <- NULL
  cd$phi_upp <- NULL
  cd$box_boundary <- NULL
  cd$boundary_R <- NULL
  if (!is.null(cd$eps) && length(cd$eps) > 1) cd$eps <- NULL
  cd$stationary_constraint <- TRUE
  cd
}

# Shared engine for the two "fix" compositions. Replicates LMCLRTest's front
# end exactly (same control surface, same estimation order, same RNG stream
# layout), then simulates TWO null distributions from the same pre-drawn
# innovations: one at the fitted theta (component 1) and one at the supplied
# fixed theta (component 2), with the estimation-stream state restored in
# between so the EM starting-value draws are repeated and the two components
# differ only through the simulation parameter. For "urfix", component 1 is
# bit-identical to LMCLRTest's p-value when mc_seed is set; for "bfix" it is
# the LMC of the transient-start family (canonical regime coordinates,
# burnin 0), which coincides with LMCLRTest only when the fitted labels are
# already canonical.
.aug_fix_lmc <- function(Y, p, k0, k1, Z, control, kind){
  fn_name <- if (kind == "urfix") "aURfixLMCTest" else "aBfixLMCTest"
  con <- list(N = 99,
              burnin = if (kind == "bfix") 0 else 100,
              converge_check = NULL,
              workers = 0,
              mc_seed = NULL,
              init_method = "warmstart",
              boundary_R = NULL,
              mdl_h0_control = list(),
              mdl_h1_control = list(),
              use_diff_init_sim = NULL)
  nmsC <- names(con)
  con[(namc <- names(control))] <- control
  if (length(noNms <- namc[!namc %in% nmsC])){
    warning("unknown names in control: ", paste(noNms, collapse = ", "))
  }
  con$init_method <- match.arg(con$init_method, c("warmstart", "random"))
  if (is.matrix(Y)){
    q <- ncol(Y)
  } else {
    stop("Observations Y must be a (T x q) matrix.")
  }
  if (k1 <= k0){
    stop(fn_name, ": k1 must be greater than k0 (received k0 = ", k0, ", k1 = ", k1, ").")
  }
  if (kind == "urfix"){
    if (q > 1){
      stop(fn_name, ": univariate models only (q = 1). The unit-root fixed point is built by the augmented Dickey-Fuller map on the fitted AR polynomial, which has no single multivariate counterpart (a VAR's unit-root null spans cointegration ranks 0 to q - 1, so a fixed point would impose an unstated rank choice); for a multivariate series use aURdiffLMCTest, whose differencing construction needs no fixed point.")
    }
    if (p < 1) stop(fn_name, ": requires p >= 1 (no autoregressive coefficient to fix at the unit root).")
    # boundary_R is a recognized control name (the shared engine reads it on
    # the boundary branch), so the unknown-names warning above cannot catch
    # it here; without this, a caller passing it would get no signal that it
    # has no effect on this test.
    if (!is.null(con$boundary_R)){
      warning(fn_name, ": 'boundary_R' belongs to the boundary family and is not used by this test; it was ignored.", call. = FALSE)
    }
  } else {
    if (k0 < 2) stop(fn_name, ": requires a null with k0 >= 2 regimes (a one-regime null has no transition matrix, so the boundary face is undefined).")
    if (con$burnin != 0){
      stop(fn_name, ": requires burnin = 0 -- the boundary family simulates every candidate transient-start (the chain initialized in the transient regime), and a burnin would consume the transient era and compare against the wrong process.")
    }
    con$boundary_R <- .validate_boundary_R(con$boundary_R, fn_name)
  }
  # ----- Observed-data estimation: exactly LMCLRTest's stream layout.
  if (!is.null(con$mc_seed)) set.seed(con$mc_seed)
  mdl_h0 <- estimMdl(Y, p, q, k0, Z, con$mdl_h0_control)
  if (kind == "bfix"){
    # The boundary family works in canonical regime coordinates (transient
    # first, absorbing last); deterministic and RNG-free, so the stream layout
    # is unchanged. At p = 0 the simulated chain is initialized in the
    # transient regime but the FIRST observation's state is already a draw
    # from P(.|transient) (the HMM simulator applies a transition before the
    # first emission).
    mdl_h0 <- canonicalize_regimes(mdl_h0, by = "persistence")
  }
  con$mdl_h1_control <- .inherit_ms_flags(con$mdl_h1_control, mdl_h0)
  con$mdl_h1_control <- .inherit_stationarity_margin(con$mdl_h1_control, con$mdl_h0_control, fn_name, p)
  if ((con$init_method == "warmstart") && (k1 > 1)){
    con$mdl_h1_control$init_theta_extra <- warmstart_theta(mdl_h0, k1,
                                                           msmu = con$mdl_h1_control[["msmu"]],
                                                           msvar = con$mdl_h1_control[["msvar"]])
  }
  mdl_h1 <- estimMdl(Y, p, q, k1, Z, con$mdl_h1_control)
  ws_extra <- con$mdl_h1_control$init_theta_extra
  if (!is.null(ws_extra) && (length(ws_extra) > 0) &&
      any(vapply(ws_extra, length, 0L) != length(mdl_h1$theta))){
    warning("the warm starts did not match the alternative model's parameter vector ",
            "and were discarded.", call. = FALSE)
  }
  con$mdl_h0_control <- mdl_h0$control
  con$mdl_h1_control <- mdl_h1$control
  con$mdl_h1_control$init_theta_extra <- NULL
  # The stationarity margin must be the same on both sides: with different
  # margins the alternative no longer nests the null's constrained space and
  # the LR statistic loses its constrained-ML meaning.
  d0 <- con[["mdl_h0_control"]][["stationarity_constraint"]]
  d1 <- con[["mdl_h1_control"]][["stationarity_constraint"]]
  if (!is.null(d0) && !is.null(d1) && !isTRUE(all.equal(as.numeric(d0), as.numeric(d1)))){
    warning("'stationarity_constraint' differs between mdl_h0_control (", d0, ") and mdl_h1_control (", d1, "); the same margin is required on both sides for the statistic to be a constrained-likelihood ratio.")
  }
  if (is.null(con$converge_check) == FALSE){
    if ((con$converge_check == "null") & (mdl_h0$converged == FALSE)){
      stop("Model under null hypothesis did not converge. Run again to use different initial values and/or increase 'maxit' for restricted model.")
    }
    if ((con$converge_check == "alt") & (mdl_h1$converged == FALSE)){
      stop("Model under alternative hypothesis did not converge. Run again to use different initial values and/or increase 'maxit' for unrestricted model.")
    }
    if ((con$converge_check == "both") & ((mdl_h0$converged == FALSE) | (mdl_h1$converged == FALSE))){
      stop("Model did not converge. Run again to use different initial values and/or increase 'maxit' for each models.")
    }
  }
  # ----- Observed statistic: fixed, shared by both components.
  logL0 <- mdl_h0$logLike
  logL1 <- mdl_h1$logLike
  if (k0 > 1){
    # H0 linear floor, by the same rule as LMCLRTest and the per-draw path:
    # the stationarity margin on the null control carries into the floor fit,
    # so the floor point is feasible in the constrained space on both sides.
    mdl_lin_con <- list(getSE = FALSE)
    if (!is.null(con[["mdl_h0_control"]][["stationarity_constraint"]])){
      mdl_lin_con$stationarity_constraint <- con[["mdl_h0_control"]][["stationarity_constraint"]]
    }
    mdl_lin <- estimMdl(Y, p, q, 1, Z, mdl_lin_con)
    if (is.finite(mdl_lin$logLike) && (mdl_lin$logLike > logL0)) logL0 <- mdl_lin$logLike
  }
  LRT_0 <- -2*(logL0 - logL1)
  if ((is.finite(LRT_0) == FALSE) | (any(is.finite(mdl_h0$theta) == FALSE)) | (any(is.finite(mdl_h1$theta) == FALSE))){
    stop("LRT_0 or model parameters are not finite. Run again to use different initial values")
  }
  if (LRT_0 < 0){
    if (LRT_0 < -1e-6){
      warning(sprintf("%s: the observed LR statistic was negative (%.3g; numerical optimization artifact) and was set to 0. Consider a higher 'use_diff_init'.", fn_name, LRT_0))
    }
    LRT_0 <- 0
  }
  names(LRT_0) <- c("LRT_0")
  # ----- The fixed second-component parameter: pure arithmetic on the fitted
  # theta, no RNG consumed, so the stream layout below matches LMCLRTest's.
  theta_fix <- as.numeric(mdl_h0$theta)
  if (kind == "urfix"){
    phi_pos <- which(as.numeric(mdl_h0$theta_phi_ind) == 1)
    theta_fix[phi_pos] <- .adf_unit_root(theta_fix[phi_pos])
  } else {
    theta_fix <- .theta_boundary_projection(mdl_h0$theta, .mmc_boundary_blocks(mdl_h0))
  }
  names(theta_fix) <- names(mdl_h0$theta)
  # ----- Simulation controls (mirrors LMCLRTest).
  if (is.null(Z) == FALSE){
    Zsim <- Z[(p + 1):nrow(Z), , drop = FALSE]
  } else {
    Zsim <- Z
  }
  mdl_h0_null_cont <- con$mdl_h0_control
  mdl_h1_null_cont <- con$mdl_h1_control
  if (is.null(con$use_diff_init_sim) == FALSE){
    if (k0 > 1) mdl_h0_null_cont$use_diff_init <- con$use_diff_init_sim
    mdl_h1_null_cont$use_diff_init <- con$use_diff_init_sim
  }
  mdl_h1_null_cont$init_method <- con$init_method
  if (con$workers > 0 && con$workers > con$N){
    message(fn_name, ": capping workers at N = ", con$N, " (", con$workers, " requested).")
    con$workers <- con$N
  }
  # ----- Worker seed: when mc_seed is set, both components pass seed = mc_seed
  # to the parallel path, exactly as LMCLRTest does (which both preserves
  # bit-identity of component 1 and re-seeds the worker streams identically
  # for component 2 -- the replay in parallel). Only when mc_seed is NULL with
  # workers > 0 is a seed drawn (from the ambient stream, BEFORE the predraw):
  # there LMCLRTest passes seed = NULL, but with pre-drawn innovations the
  # parallel path skips the worker RNG reset entirely, so without a real seed
  # component 2's worker streams would not replay component 1's.
  wseed <- con$mc_seed
  if (is.null(con$mc_seed) && con$workers > 0){
    wseed <- sample.int(.Machine$integer.max, 1L)
  }
  # ----- Pre-draw innovations (always -- both components must share them).
  if (!is.null(con$mc_seed)) set.seed(con$mc_seed)
  N_buffer <- ceiling(con$N * 1.5) + 10
  Teps <- mdl_h0$n + con$burnin
  predrawn_eps <- lapply(seq_len(N_buffer), function(i) matrix(rnorm(Teps * q), Teps, q))
  predrawn_state_rand <- NULL
  if (k0 > 1){
    predrawn_state_rand <- matrix(runif(Teps * N_buffer), Teps, N_buffer)
  }
  # ----- Component 1 (the LMC), then the estimation-stream replay, then
  # component 2 at theta_fix. Each component's p-value is computed IMMEDIATELY
  # after its own null draws (before any restore), matching LMCLRTest's stream
  # position for the randomized tie-break inside MCpval -- ties at 0 are
  # common under the LR >= 0 cap, and MCpval consumes random numbers on every
  # call. In parallel mode the master stream barely advances during the
  # parallel simulation, so the two components' tie-break draws coincide;
  # harmless, since the intersection-union composition needs no independence
  # between its components.
  S_anchor <- NULL
  if (exists(".Random.seed", envir = globalenv(), inherits = FALSE)){
    S_anchor <- get(".Random.seed", envir = globalenv())
  }
  sim_one <- function(mdl_sim){
    if (con$workers > 0){
      LR_samp_dist_par(mdl_sim, k1, con$N, con$burnin, Zsim, mdl_h0_null_cont, mdl_h1_null_cont,
                       con$workers, seed = wseed,
                       predrawn_eps = predrawn_eps, predrawn_state_rand = predrawn_state_rand)
    } else {
      LR_samp_dist(mdl_sim, k1, con$N, con$burnin, Zsim, mdl_h0_null_cont, mdl_h1_null_cont,
                   predrawn_eps, predrawn_state_rand)
    }
  }
  LRN_1 <- sim_one(mdl_h0)
  LRN_1 <- LRN_1[is.finite(LRN_1)]
  if (length(LRN_1) == 0){
    stop(fn_name, ": the null distribution could not be simulated at the fitted parameter values; the p-value cannot be computed. Increase 'use_diff_init' or inspect the observed model fit.")
  }
  pval_lmc <- MCpval(LRT_0, LRN_1, "geq")
  if (!is.null(S_anchor)) assign(".Random.seed", S_anchor, envir = globalenv())
  mdl_fix <- mdledit(mdl_h0, theta_fix, p, q, k0, !is.null(Z))
  LRN_2 <- sim_one(mdl_fix)
  LRN_2 <- LRN_2[is.finite(LRN_2)]
  if (length(LRN_2) == 0){
    stop(fn_name, ": the null distribution could not be simulated at the fixed (", if (kind == "urfix") "unit-root" else "boundary", ") parameter values; the second component's p-value cannot be computed.")
  }
  pval_fix <- MCpval(LRT_0, LRN_2, "geq")
  # ----- Composition and output.
  pval <- max(pval_lmc, pval_fix)
  LRN_1 <- as.matrix(sort(LRN_1))
  LRN_2 <- as.matrix(sort(LRN_2))
  LRN_cv <- LRN_1[round(c(0.90, 0.95, 0.99) * nrow(LRN_1)), ]
  LRN_fix_cv <- LRN_2[round(c(0.90, 0.95, 0.99) * nrow(LRN_2)), ]
  names(LRN_cv) <- paste0(c("0.90", "0.95", "0.99"), "%")
  names(LRN_fix_cv) <- paste0(c("0.90", "0.95", "0.99"), "%")
  out <- list(mdl_h0 = mdl_h0, mdl_h1 = mdl_h1, LRT_0 = LRT_0,
              pval_lmc = pval_lmc, pval_fix = pval_fix, pval = pval,
              LRN = LRN_1, LRN_fix = LRN_2, LRN_cv = LRN_cv, LRN_fix_cv = LRN_fix_cv,
              theta_fix = theta_fix, control = con)
  class(out) <- if (kind == "urfix") "aURfixLMCTest" else "aBfixLMCTest"
  out
}

#' @title Unit-root-augmented (fix) Local Monte Carlo Likelihood Ratio Test
#'
#' @description Runs the Local Monte Carlo likelihood ratio test of
#' Rodriguez-Rondon & Dufour (2026a) twice on the same observed series and the
#' same pre-drawn simulation draws: once with the null distribution simulated
#' at the fitted null-model parameters (the ordinary LMC-LRT), and once with
#' the null-simulation parameter's long-run autoregressive coefficient moved
#' to exactly 1 (everything else unchanged), so the second component's
#' simulated samples have an exact unit root -- the Monte Carlo counterpart of
#' imposing the unit root in the resampling scheme (Basawa, Mallik, McCormick,
#' Reeves & Taylor, 1991). The test rejects only if BOTH component p-values
#' are below the chosen level (an intersection-union composition, reported as
#' \code{pval = max(pval_lmc, pval_fix)}), which protects the level when the
#' data may carry a unit root. For \code{p > 1} the long-run coefficient
#' \code{sum(phi)} is moved to 1 while the short-run (augmented Dickey-Fuller)
#' coefficients are held fixed -- linear arithmetic on the fitted vector,
#' nothing re-estimated; the resulting vector is reported as \code{theta_fix}.
#' Exogenous regressors, when present, keep their estimated coefficients in
#' the unit-root-fixed simulation parameter. Estimation on the simulated
#' samples stays free (unrestricted) in both components, identically to the
#' observed side. At a one-regime null with \code{p = 1} and no exogenous
#' regressors the imposed-root component's guarantee at the unit root is
#' exact in finite samples: the fixed simulation parameter carries no
#' estimated autoregressive coefficient, and the remaining location-scale
#' nuisance is absorbed by invariance (the level cancels at the unit root,
#' the scale by equivariance). For \code{p >= 2} the short-run coefficients
#' enter the fixed parameter at their estimates and the component is
#' asymptotically valid at the unit root, as it is with regressors present;
#' it is the differencing construction that retains finite-sample exactness
#' at \code{p = 1} with regressors -- see \code{\link{aURdiffLMCTest}} when
#' exactness with \code{Z} is needed.
#' Univariate only: the augmented Dickey-Fuller map has no single multivariate
#' counterpart (a VAR's unit-root null spans cointegration ranks 0 to
#' \code{q - 1}, so a fixed point would impose an unstated rank choice); for a
#' multivariate series use \code{\link{aURdiffLMCTest}}, whose differencing
#' construction needs no fixed point.
#'
#' @param Y  Series to be tested. Must be a (\code{T x 1}) matrix (univariate only).
#' @param p  Number of autoregressive lags. Must be greater than or equal to 1.
#' @param k0 Number of regimes under null hypothesis. Must be greater than or equal to 1.
#' @param k1 Number of regimes under alternative hypothesis. Must be greater than \code{k0}.
#' @param Z Optional (\code{T x qz}) matrix of fixed (deterministic) or strictly exogenous regressors, or \code{NULL}. Dummy variables are deterministic, hence fixed regressors, and are fully supported. The regressors' common coefficients never switch and are identified by the whole sample, so both components treat them like the other common parameters.
#' @param control List with test procedure options: \code{N}, \code{burnin},
#'   \code{converge_check}, \code{workers}, \code{mc_seed}, \code{init_method},
#'   \code{mdl_h0_control}, \code{mdl_h1_control}, \code{use_diff_init_sim},
#'   with the same meanings and defaults as in \code{\link{LMCLRTest}}. With
#'   \code{mc_seed} set, \code{pval_lmc} is identical to \code{LMCLRTest}'s
#'   p-value for the same inputs, and the second component reuses the same
#'   pre-drawn innovations and the same estimation-stream randomness, so the
#'   two components differ only through the simulation parameter.
#'
#' @return List of class \code{aURfixLMCTest} with \code{mdl_h0},
#'   \code{mdl_h1}, \code{LRT_0} (the shared observed statistic),
#'   \code{pval_lmc} and \code{pval_fix} (the component p-values), \code{pval}
#'   (their maximum; reject only if it is below the level), \code{LRN} and
#'   \code{LRN_fix} (the two simulated null distributions), \code{LRN_cv} and
#'   \code{LRN_fix_cv}, \code{theta_fix} (the unit-root-fixed parameter
#'   vector), and \code{control}.
#'
#' @references Rodriguez-Rondon, G., & Dufour, J.-M. 2026a. "Monte Carlo Likelihood-Ratio Tests for Markov Switching Models." \emph{Bank of Canada Staff Working Paper}, No. 2026-23. doi: 10.34989/swp-2026-23.
#' @references Basawa, I. V., Mallik, A. K., McCormick, W. P., Reeves, J. H., & Taylor, R. L. 1991. "Bootstrapping Unstable First-Order Autoregressive Processes." \emph{The Annals of Statistics} 19(2): 1098-1101.
#' @references Basawa, I. V., Mallik, A. K., McCormick, W. P., & Taylor, R. L. 1991. "Bootstrap Test of Significance and Sequential Bootstrap Estimation for Unstable First Order Autoregressive Processes." \emph{Communications in Statistics - Theory and Methods} 20(3): 1015-1026.
#' @export
aURfixLMCTest <- function(Y, p, k0, k1, Z = NULL, control = list()){
  .aug_fix_lmc(Y, p, k0, k1, Z, control, kind = "urfix")
}

#' @title Boundary-augmented (fix) Local Monte Carlo Likelihood Ratio Test
#'
#' @description Runs the Local Monte Carlo likelihood ratio test of
#' Rodriguez-Rondon & Dufour (2026a) twice on the same observed series and the
#' same pre-drawn simulation draws, for a \code{k0 >= 2} null in the boundary
#' (possibly absorbing regime) setting: once at the fitted null-model
#' parameters, and once at their boundary projection -- the absorbing regime's
#' self-transition probability set to exactly 1 (its column's remaining
#' entries to 0), everything else at its estimate. The regimes are first
#' relabeled canonically (\code{\link{canonicalize_regimes}}: transient
#' regimes first, absorbing regime last), and every simulated sample starts
#' with the chain initialized in the transient regime with no burn-in (the
#' change-point-like family of Chib, 1998) -- which is why \code{burnin = 0}
#' is required. The test rejects only if BOTH component p-values are below the
#' chosen level (reported as \code{pval = max(pval_lmc, pval_fix)}). At
#' \code{p = 0} the first observation's regime is a draw from the transient
#' regime's transition distribution (the hidden Markov simulator applies a
#' transition before the first emission). Multivariate models are supported:
#' the boundary projection edits only the absorbing regime's transition-matrix
#' column, so the fixed simulation parameter keeps the fitted regime
#' covariances -- estimates of an admissible model, positive definite by
#' construction.
#'
#' @param Y  Series to be tested. Must be a (\code{T x q}) matrix.
#' @param p  Number of autoregressive lags. Must be greater than or equal to 0.
#' @param k0 Number of regimes under null hypothesis. Must be greater than or equal to 2. With \code{k0 >= 3} the boundary family's initial state is canonical state 1 of the fitted null (the least persistent transient regime at the fit, inherited by index for the simulation parameter) -- a fixed-vertex family definition on which the finite-sample validity is conditional.
#' @param k1 Number of regimes under alternative hypothesis. Must be greater than \code{k0}.
#' @param Z Optional (\code{T x qz}) matrix of fixed (deterministic) or strictly exogenous regressors, or \code{NULL}. Dummy variables are deterministic, hence fixed regressors, and are fully supported. The regressors' common coefficients never switch and are identified by the whole sample, so both components treat them like the other common parameters.
#' @param control List with test procedure options as in
#'   \code{\link{aURfixLMCTest}}, except that \code{burnin} defaults to 0 and
#'   must be 0. \code{boundary_R} is accepted for interface symmetry with
#'   \code{\link{MMCLRTest}}'s boundary modes but is not used here (the second
#'   component is a single fixed evaluation, not a search).
#'
#' @return List of class \code{aBfixLMCTest} with the same structure as
#'   \code{\link{aURfixLMCTest}} (component p-values \code{pval_lmc},
#'   \code{pval_fix}, composition \code{pval}, the two null distributions, and
#'   \code{theta_fix} -- here the boundary projection). \code{mdl_h0} is
#'   returned in canonical regime order.
#'
#' @references Rodriguez-Rondon, G., & Dufour, J.-M. 2026a. "Monte Carlo Likelihood-Ratio Tests for Markov Switching Models." \emph{Bank of Canada Staff Working Paper}, No. 2026-23. doi: 10.34989/swp-2026-23.
#' @references Chib, S. 1998. "Estimation and Comparison of Multiple Change-Point Models." \emph{Journal of Econometrics} 86(2): 221-241.
#' @export
aBfixLMCTest <- function(Y, p, k0, k1, Z = NULL, control = list()){
  .aug_fix_lmc(Y, p, k0, k1, Z, control, kind = "bfix")
}

#' @title Unit-root-augmented (differencing) Local Monte Carlo Likelihood Ratio Test
#'
#' @description Runs the Local Monte Carlo likelihood ratio test of
#' Rodriguez-Rondon & Dufour (2026a) twice: on the observed series in levels
#' with \code{p} autoregressive lags, and on its first difference with
#' \code{p - 1} lags (at \code{p = 1} the differenced null has no lagged term).
#' The test rejects only if BOTH component p-values are below the chosen level
#' (an intersection-union composition, reported as
#' \code{pval = max(pval_levels, pval_diff)}). Differencing imposes the null
#' unit root without estimating it: at a one-regime null with \code{p = 1} the
#' differenced component's statistic does not depend on any estimated
#' autoregressive parameter, and its unit-root guarantee is exact in finite
#' samples -- including with exogenous regressors present. For \code{p >= 2}
#' the differenced model keeps lagged differences, and the component is
#' asymptotically valid at the unit root rather than exact. Exogenous
#' regressors are kept IN LEVELS for the differenced component (rows 2 to T of
#' \code{Z}) -- they enter the observation equation additively, outside the
#' mean adjustment, so differencing the series does not difference them. Both
#' components are full \code{\link{LMCLRTest}} runs sharing the same
#' \code{control} (including \code{mc_seed}, so their simulation streams start
#' from the same seed position); each holds its own pre-drawn innovations
#' fixed internally. Multivariate series run the same way -- the difference is
#' taken per series and both components are ordinary multivariate tests -- but
#' the unit-root validity argument is a univariate construction: a
#' multivariate non-stationary null spans cointegration ranks 0 to q - 1, and
#' first-differencing all series imposes the no-cointegration case, so at
#' \code{q > 1} the composition's unit-root guarantee is conditional on that
#' maintained assumption.
#'
#' @param Y  Series to be tested. Must be a (\code{T x q}) matrix.
#' @param p  Number of autoregressive lags for the levels component. Must be greater than or equal to 1.
#' @param k0 Number of regimes under null hypothesis. Must be greater than or equal to 1.
#' @param k1 Number of regimes under alternative hypothesis. Must be greater than \code{k0}.
#' @param Z Optional (\code{T x qz}) matrix of fixed (deterministic) or strictly exogenous regressors, or \code{NULL}. Dummy variables are deterministic, hence fixed regressors, and are fully supported. The regressors' common coefficients never switch and are identified by the whole sample, so both components treat them like the other common parameters.
#' @param control List with test procedure options, passed to BOTH component
#'   \code{\link{LMCLRTest}} calls unchanged.
#'
#' @return List of class \code{aURdiffLMCTest} with \code{pval_levels},
#'   \code{pval_diff} (the component p-values), \code{pval} (their maximum;
#'   reject only if it is below the level), \code{lmc_levels} and
#'   \code{lmc_diff} (the two full \code{LMCLRTest} objects), and
#'   \code{control}.
#'
#' @references Rodriguez-Rondon, G., & Dufour, J.-M. 2026a. "Monte Carlo Likelihood-Ratio Tests for Markov Switching Models." \emph{Bank of Canada Staff Working Paper}, No. 2026-23. doi: 10.34989/swp-2026-23.
#' @export
aURdiffLMCTest <- function(Y, p, k0, k1, Z = NULL, control = list()){
  if (!is.matrix(Y)) stop("Observations Y must be a (T x q) matrix.")
  if (p < 1) stop("aURdiffLMCTest: requires p >= 1 (the differenced component uses p - 1 lags).")
  # The differenced observation y_t - y_{t-1} exists for t = 2..T, so the
  # differenced component keeps rows 2..T of Z -- in LEVELS. Z is never
  # differenced: it enters the observation equation additively, outside the
  # mean adjustment, so differencing y does not difference Z.
  Z_diff <- if (is.null(Z)) NULL else Z[-1, , drop = FALSE]
  lmc_levels <- LMCLRTest(Y, p, k0, k1, Z = Z, control = control)
  lmc_diff   <- LMCLRTest(diff(Y), p - 1L, k0, k1, Z = Z_diff, control = control)
  out <- list(pval_levels = lmc_levels$pval, pval_diff = lmc_diff$pval,
              pval = max(lmc_levels$pval, lmc_diff$pval),
              lmc_levels = lmc_levels, lmc_diff = lmc_diff,
              control = lmc_levels$control)
  class(out) <- "aURdiffLMCTest"
  out
}

#' @title Unit-root-augmented (differencing) Maximized Monte Carlo Likelihood Ratio Test
#'
#' @description Runs the Maximized Monte Carlo likelihood ratio test of
#' Rodriguez-Rondon & Dufour (2026a) twice: on the observed series in levels
#' with \code{p} autoregressive lags (typically with the unit-root search box
#' of \code{\link{MMCLRTest}}'s \code{box_include} control, so the search set
#' contains the unit root), and on its first difference with \code{p - 1}
#' lags, each with its own consistent-set search box around its own fitted
#' null model. The test rejects only if BOTH component p-values are below the
#' chosen level (reported as \code{pval = max(pval_levels, pval_diff)}).
#' Multivariate series are handled as in \code{\link{aURdiffLMCTest}}, with
#' the same scope note: first-differencing all series imposes the
#' no-cointegration case, so at \code{q > 1} the unit-root guarantee is
#' conditional on that maintained assumption.
#'
#' @param Y  Series to be tested. Must be a (\code{T x q}) matrix.
#' @param p  Number of autoregressive lags for the levels component. Must be greater than or equal to 1.
#' @param k0 Number of regimes under null hypothesis. Must be greater than or equal to 1.
#' @param k1 Number of regimes under alternative hypothesis. Must be greater than \code{k0}.
#' @param Z Optional (\code{T x qz}) matrix of fixed (deterministic) or strictly exogenous regressors, or \code{NULL}. Dummy variables are deterministic, hence fixed regressors, and are fully supported. The regressors' common coefficients never switch and are identified by the whole sample, so both components treat them like the other common parameters.
#' @param control List with test procedure options for the LEVELS component,
#'   passed to \code{\link{MMCLRTest}} unchanged (for a unit-root setting,
#'   pass \code{stationary_constraint = FALSE}, leave \code{phi_upp} unset,
#'   and put the unit-root point in \code{box_include}).
#' @param control_diff List with test procedure options for the DIFFERENCED
#'   component. Default \code{NULL} derives it from \code{control} by removing
#'   the levels-specific box settings -- \code{box_include} (whose length is
#'   the levels parameter vector's), \code{phi_low}/\code{phi_upp}, the
#'   boundary controls, and any vector-valued \code{eps} -- and restoring
#'   \code{stationary_constraint = TRUE}, since the differenced component is
#'   the stationary world by construction. Pass it explicitly to override.
#'
#' @return List of class \code{aURdiffMMCTest} with \code{pval_levels},
#'   \code{pval_diff}, \code{pval} (their maximum), \code{mmc_levels} and
#'   \code{mmc_diff} (the two full \code{MMCLRTest} objects), \code{control},
#'   and \code{control_diff}.
#'
#' @references Rodriguez-Rondon, G., & Dufour, J.-M. 2026a. "Monte Carlo Likelihood-Ratio Tests for Markov Switching Models." \emph{Bank of Canada Staff Working Paper}, No. 2026-23. doi: 10.34989/swp-2026-23.
#' @export
aURdiffMMCTest <- function(Y, p, k0, k1, Z = NULL, control = list(), control_diff = NULL){
  if (!is.matrix(Y)) stop("Observations Y must be a (T x q) matrix.")
  if (p < 1) stop("aURdiffMMCTest: requires p >= 1 (the differenced component uses p - 1 lags).")
  if (is.null(control_diff)) control_diff <- .aug_diff_control(control)
  # Z handling as in aURdiffLMCTest: rows 2..T, in levels, never differenced.
  Z_diff <- if (is.null(Z)) NULL else Z[-1, , drop = FALSE]
  mmc_levels <- MMCLRTest(Y, p, k0, k1, Z = Z, control = control)
  mmc_diff   <- MMCLRTest(diff(Y), p - 1L, k0, k1, Z = Z_diff, control = control_diff)
  out <- list(pval_levels = mmc_levels$pval, pval_diff = mmc_diff$pval,
              pval = max(mmc_levels$pval, mmc_diff$pval),
              mmc_levels = mmc_levels, mmc_diff = mmc_diff,
              control = mmc_levels$control, control_diff = mmc_diff$control)
  class(out) <- "aURdiffMMCTest"
  out
}
