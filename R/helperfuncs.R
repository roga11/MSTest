#' @title Companion Matrix
#'
#' @description This function converts the (\code{q x 1})  vector of constants and (\code{q x qp}) matrix of autoregressive coefficients into (\code{qp x qp}) matrix belonging to the companion form
#'
#' @param phi matrix of dimension (\code{q x qp}) containing autoregressive coefficients.
#' @param p integer for number of autoregressive lags.
#' @param q integer for number of series.
#'
#' @return matrix of dimension (\code{qp x qp}) of companion form.
#' 
#' @keywords internal
#' 
#' @export
companionMat <- function(phi, p, q){
  interzero   <- matrix(0,q*(p-1), 1)
  F_tmp       <- phi
  diagmat     <- diag(q*(p-1))
  diagzero    <- matrix(0, q*(p-1), q)
  Mn          <- cbind(diagmat,diagzero)
  compMat     <- rbind(F_tmp,Mn)
  return(compMat)
}


#' @title Autoregressive transition matrix
#'
#' @description This function converts a transition matrix to the transition matrix consistent with a Markov-switching autoregressive model.
#'
#' @param P original transition matrix.
#' @param k integer determining the number of regimes.
#' @param ar number of autoregressive lags.
#'
#' @return transformed transition matrix.
#' 
#' @keywords internal
#' 
#' @export
arP <- function(P, k, ar){
  if (ar>0){
    Pnew_tmp <- matrix(0,k*k,k)
    for (xk in 1:k){
      ind <- rep(FALSE,k)
      ind[xk] <- TRUE
      Pnew_tmp[(1+(xk-1)*k):(1+(xk-1)*k+(k-1)),xk] <- P[,ind]
    }
    repmat <- k^(ar-1)
    Pnew <- do.call(cbind,replicate(k,kronecker(diag(repmat), Pnew_tmp),simplify=FALSE)) 
  }else{
    Pnew <- P
  }
  return(Pnew)  
}

#' @title Autoregressive moment grid
#'
#' @description This function creates a grid of mean and variance consistent with a Markov-switching autoregressive model.
#'
#' @param mu vector (\code{k x 1}) of mu in each regime.
#' @param sig vector (\code{k x 1}) of sigma in each regime.
#' @param k integer determining the number of regimes.
#' @param ar number of autoregressive lags.
#' @param msmu Boolean indicator. If \code{TRUE} mean is subject to change. If \code{FALSE} mean is constant across regimes.
#' @param msvar Boolean indicator. If \code{TRUE} variance is subject to change. If \code{FALSE} variance is constant across regimes.
#'
#' @return List with (\code{M x ar+1}) matrix of means for each regime \code{M} (where \code{M = k^(ar+1)}) and each time \code{t,... t-ar}, vector with variance for each regime \code{M}, and vector indicating the corresponded \code{1,..., k} regime. 
#' 
#' @keywords internal
#'
#' @export
argrid_MSARmdl <- function(mu, sig, k, ar, msmu, msvar){
  if (msmu==FALSE){
    mu <- rep(mu, k)
  }
  if (msvar==FALSE){
    sig <- rep(sig, k)
  }
  mu <- as.vector(mu)
  sig <- as.vector(sig)
  # create grid of regimes
  lx <-list()
  for (xi in 1:(ar+1)) lx[[xi]] <- 1:k
  mu_stategrid <- as.matrix(expand.grid(lx))
  sig_stategrid <- as.matrix(expand.grid(lx))
  state_indicator = mu_stategrid[,1]
  # get matrix of mu and sigma
  mu_stategrid[] <- mu[mu_stategrid]
  sig_stategrid[] <- sig[sig_stategrid]
  sig_stategrid <- sig_stategrid[,1] # sig doesnt have AR seq
  musig_out <- list()
  musig_out[["mu"]] <- as.matrix(mu_stategrid)
  musig_out[["sig"]] <- as.matrix(sig_stategrid)
  musig_out[["state_ind"]] <- as.matrix(state_indicator)
  return(musig_out)
}


#' @title Vector autoregressive moment grid
#'
#' @description Creates grid of means and covariance matrices consistent with a Markov-switching vector autoregressive model.
#'
#' @param mu a (\code{k x q}) matrix of means in each regime (for \code{k} regimes and \code{q} time series).
#' @param sigma list with \code{k} regime specific (\code{q x q}) covariance matrices.
#' @param k integer determining the number of regimes.
#' @param ar number of autoregressive lags.
#' @param msmu Boolean indicator. If \code{TRUE} mean is subject to change. If \code{FALSE} mean is constant across regimes.
#' @param msvar Boolean indicator. If \code{TRUE} variance is subject to change. If \code{FALSE} variance is constant across regimes.
#'
#' @return List with M regime specific (\code{q x k}) matrices of means, List with \code{M} regime specific covariance matrices, and vector indicating the corresponded \code{1,..., k} regime. 
#' 
#' @keywords internal
#' 
#' @export
argrid_MSVARmdl <- function(mu, sigma, k, ar, msmu, msvar){
  # create grid of regimes
  lx <-list()
  for (xi in 1:(ar+1)) lx[[xi]] <- 1:k
  mu_stategrid <- as.matrix(expand.grid(lx))
  state_indicator <- mu_stategrid[,1]
  M <- k^(ar+1) 
  sig_stategrid <- list()
  muN_stategrid <- list()
  for (xm in 1:M){
    sig_stategrid[[xm]] <- sigma[[state_indicator[xm]]]
    muN_stategrid[[xm]] <- t(mu[mu_stategrid[xm,],])
  }
  musig_out <- list()
  musig_out[["mu"]] <- muN_stategrid
  musig_out[["sig"]] <- sig_stategrid
  musig_out[["state_ind"]] <- as.matrix(state_indicator)
  return(musig_out)
}




#' @title Theta standard errors
#'
#' @description This function computes the standard errors of the parameters in vector theta. This is done using an approximation of the Hessian matrix (using \code{\link[numDeriv]{hessian}} and \code{nearPD} if \code{info_mat} is not PD).
#'
#' For Markov switching models (k > 1), a reduced parameterization is used for the
#' transition matrix P to avoid rank deficiency. The k x k transition matrix has k
#' column-sum constraints, so only k(k-1) entries are free. The Hessian is computed
#' w.r.t. these free entries, and SEs for the constrained entries are obtained via the
#' delta method. See Krolzig (1997, Sections 6.3.2 and 6.6).
#'
#' @param mdl List with model properties
#'
#' @return List provided as input with additional attributes \code{HESS},\code{theta_se}, \code{info_mat}, and \code{nearPD_used}.
#'
#' @keywords internal
#'
#' @export
thetaSE <- function(mdl){
  k <- if (!is.null(mdl$k)) mdl$k else 1L
  if (k > 1L && inherits(mdl, c("HMmdl", "MSARmdl", "MSVARmdl"))) {
    return(.thetaSE_reduced(mdl))
  }
  Hess <- getHessian(mdl)
  nearPD_used <- FALSE
  if (any(!is.finite(Hess))) {
    warning("Hessian contains non-finite values; standard errors cannot be computed.")
    mdl$Hess <- Hess
    mdl$theta_se <- rep(NA_real_, length(mdl$theta))
    mdl$info_mat <- matrix(NA_real_, nrow(Hess), ncol(Hess))
    mdl$nearPD_used <- nearPD_used
    return(mdl)
  }
  cond_num <- tryCatch(kappa(-Hess, exact = FALSE), error = function(e) Inf)
  if (cond_num > 1e10)
    warning("Ill-conditioned Hessian (cond. number approx. ", format(cond_num, scientific = TRUE, digits = 2),
            "); standard errors may be unreliable.")
  info_mat <- tryCatch(solve(-Hess),
                       error = function(e) matrix(NA_real_, nrow(Hess), ncol(Hess)))
  if (any(!is.finite(info_mat)) || any(diag(info_mat) < 0)) {
    info_mat <- tryCatch(solve(pracma::nearest_spd(-Hess)),
                         error = function(e) matrix(NA_real_, nrow(Hess), ncol(Hess)))
    nearPD_used <- TRUE
  }
  mdl$Hess <- Hess
  mdl$theta_se <- sqrt(pmax(diag(info_mat), 0))
  mdl$info_mat <- info_mat
  mdl$nearPD_used <- nearPD_used
  return(mdl)
}

# ============================================================================== #
# Reduced parameterization helpers for transition matrix P (Krolzig 1997, Sec. 6.3.2)
# ============================================================================== #
#
# The k x k transition matrix P is column-stochastic (columns sum to 1). This means
# only k(k-1) of the k^2 entries are free. Computing the Hessian w.r.t. all k^2
# entries produces a rank-deficient matrix. The fix: project the full Hessian to the
# k(k-1)-dimensional free parameter space, invert there, then expand back.
#
# Since theta_full = J * theta_free + b is a LINEAR mapping, the chain rule gives
# H_reduced = J' H_full J (exact, no approximation). SEs for constrained entries
# are recovered via V_full = J V_free J' (delta method, exact for linear constraints).
# See Krolzig (1997, Sec. 6.6.1-6.6.2).
#
# Drop pattern (preserves all diagonal/persistence entries as free):
#   Columns 1 to k-1: drop P[k, j] (last row)
#   Column k: drop P[1, k] (first row)
# ============================================================================== #

# Compute the reduced parameterization setup: Jacobian J and index mappings
.reduced_P_setup <- function(k, npar) {
  P_start <- npar - k * k + 1L
  # P[i,j] at block position (j-1)*k + i (column-major)
  drop_block <- sort(c(
    seq_len(k - 1L) * k,     # P[k, 1], ..., P[k, k-1]
    (k - 1L) * k + 1L        # P[1, k]
  ))
  free_block <- setdiff(seq_len(k * k), drop_block)
  drop_theta <- P_start + drop_block - 1L
  free_theta <- setdiff(seq_len(npar), drop_theta)
  nfree      <- length(free_theta)
  # Per-column info for Jacobian construction
  drop_col_info <- lapply(seq_len(k), function(j) {
    col_block  <- ((j - 1L) * k + 1L):(j * k)
    col_theta  <- P_start + col_block - 1L
    list(
      constrained_pos   = intersect(drop_theta, col_theta),
      free_pos_in_theta = setdiff(col_theta, drop_theta)
    )
  })
  # Build Jacobian J (npar x nfree): theta_full = J * theta_free + b
  J <- matrix(0, npar, nfree)
  for (idx in seq_len(nfree))
    J[free_theta[idx], idx] <- 1
  for (j in seq_along(drop_col_info)) {
    ci <- drop_col_info[[j]]
    free_in_col_idx <- match(ci$free_pos_in_theta, free_theta)
    J[ci$constrained_pos, free_in_col_idx] <- -1
  }
  list(
    J          = J,
    free_theta = free_theta,
    drop_theta = drop_theta,
    nfree      = nfree,
    drop_col_info = drop_col_info
  )
}

# thetaSE with reduced parameterization for MS models
# Projects the full Hessian to the free-parameter subspace via J' H J
.thetaSE_reduced <- function(mdl) {
  theta <- mdl$theta
  npar  <- length(theta)
  k     <- mdl$k
  nearPD_used <- FALSE
  # 1. Compute full Hessian (same code path as linear models)
  Hess_full <- getHessian(mdl)
  if (any(!is.finite(Hess_full))) {
    warning("Hessian contains non-finite values; standard errors cannot be computed.")
    mdl$Hess      <- Hess_full
    mdl$theta_se  <- rep(NA_real_, npar)
    mdl$info_mat  <- matrix(NA_real_, npar, npar)
    mdl$nearPD_used <- FALSE
    return(mdl)
  }
  # 2. Project to reduced space: H_reduced = J' H_full J
  setup      <- .reduced_P_setup(k, npar)
  J          <- setup$J
  Hess_free  <- t(J) %*% Hess_full %*% J
  # 3. Check conditioning
  cond_num <- tryCatch(kappa(-Hess_free, exact = FALSE), error = function(e) Inf)
  if (cond_num > 1e10)
    warning("Ill-conditioned Hessian (cond. number approx. ",
            format(cond_num, scientific = TRUE, digits = 2),
            "); standard errors may be unreliable.")
  # 4. Invert in reduced space
  V_free <- tryCatch(solve(-Hess_free),
                     error = function(e) matrix(NA_real_, setup$nfree, setup$nfree))
  if (any(!is.finite(V_free)) || any(diag(V_free) < 0)) {
    V_free <- tryCatch(solve(pracma::nearest_spd(-Hess_free)),
                       error = function(e) matrix(NA_real_, setup$nfree, setup$nfree))
    nearPD_used <- TRUE
  }
  # 5. Expand back: V_full = J V_free J'
  V_full <- J %*% V_free %*% t(J)
  mdl$Hess        <- Hess_free
  mdl$theta_se    <- sqrt(pmax(diag(V_full), 0))
  mdl$info_mat    <- V_full
  mdl$nearPD_used <- nearPD_used
  return(mdl)
}

#' @title Louis (1982) standard errors for Markov-switching models
#'
#' @description Computes standard errors using the expected complete-data information
#' matrix \code{B} from Louis (1982). \code{B = -Hessian(Q)}, where \code{Q(theta)}
#' is the expected complete-data log-likelihood (the EM Q-function), evaluated at the
#' converged parameter values with smoothed state probabilities held fixed.
#'
#' The complete-data information \code{B} is typically better conditioned than the
#' observed-data Hessian. Since \code{I_obs = B - M} where \code{M >= 0} is the
#' missing information, \code{B >= I_obs} (PSD), so \code{B^{-1}} slightly
#' underestimates the true variance. The SEs are therefore slightly optimistic but
#' numerically more stable, especially for weakly identified models.
#'
#' @param mdl List with model properties (output from an MS model constructor).
#'   Must contain \code{theta}, \code{St} (smoothed probabilities), and model-specific fields.
#'
#' @return List provided as input with additional attributes \code{theta_se}, \code{info_mat},
#'   and \code{louis_used = TRUE}.
#'
#' @keywords internal
#'
#' @references Louis, T. A. (1982). "Finding the Observed Information Matrix
#'   When Using the EM Algorithm." \emph{Journal of the Royal Statistical Society,
#'   Series B}, 44(2), 226-233.
#'
#' @export
thetaSE_louis <- function(mdl) {
  theta  <- mdl$theta
  npar   <- length(theta)
  k      <- mdl$k

  # --- Recover extended-state smoothed probabilities ---
  # mdl$St is always T x k (marginal). For MSAR/MSVAR models with extended states

  # (M = k^(p+1)), we re-run the forward-backward pass to get T x M probabilities.
  xi   <- NULL
  exog <- FALSE

  if (inherits(mdl, "HMmdl")) {
    xi   <- mdl$St  # T x k (no extended state for HMM)
    exog <- isTRUE(mdl$exog)
  } else if (inherits(mdl, "MSARmdl")) {
    exog  <- !is.null(mdl$control$Z)
    em_fn <- if (exog) ExpectationM_MSARXmdl else ExpectationM_MSARmdl
    em_out <- tryCatch(em_fn(theta, mdl, k), error = function(e) NULL)
    if (!is.null(em_out)) xi <- em_out$xi_t_T_AR  # T x M
  } else if (inherits(mdl, "MSVARmdl")) {
    exog  <- !is.null(mdl$control$Z)
    em_fn <- if (exog) ExpectationM_MSVARXmdl else ExpectationM_MSVARmdl
    em_out <- tryCatch(em_fn(theta, mdl, k), error = function(e) NULL)
    if (!is.null(em_out)) xi <- em_out$xi_t_T_AR  # T x M
  }

  if (is.null(xi)) {
    warning("thetaSE_louis: could not recover smoothed probabilities; falling back to Hessian-based SEs.")
    return(thetaSE(mdl))
  }

  # --- Build Q function (expected complete-data log-likelihood) ---
  Q_fn <- tryCatch(.louis_build_Q(mdl, xi), error = function(e) {
    warning("thetaSE_louis: Q function build failed: ", e$message)
    NULL
  })
  if (is.null(Q_fn)) return(thetaSE(mdl))

  # --- Compute full Q-function Hessian, then project to reduced space ---
  # The Q function has the same P column-sum constraints as the observed-data
  # loglikelihood. Project via H_reduced = J' H_full J (exact for linear constraints).
  Q_Hess_full <- tryCatch(
    numDeriv::hessian(Q_fn, theta, method = "Richardson"),
    error = function(e) {
      warning("Louis (1982): Hessian computation failed: ", e$message)
      matrix(NA_real_, npar, npar)
    }
  )

  if (any(!is.finite(Q_Hess_full))) {
    warning("Louis (1982): Q-function Hessian contains non-finite values; standard errors cannot be computed.")
    mdl$theta_se    <- rep(NA_real_, npar)
    mdl$info_mat    <- matrix(NA_real_, npar, npar)
    mdl$louis_used  <- TRUE
    return(mdl)
  }

  # Project to reduced space
  setup  <- .reduced_P_setup(k, npar)
  J      <- setup$J
  nfree  <- setup$nfree
  Q_Hess <- t(J) %*% Q_Hess_full %*% J
  B      <- -Q_Hess  # expected complete-data information (positive semi-definite)

  # --- Eigenvalue decomposition for robust inversion ---
  # Standard solve() fails silently when B is near-singular (cond >> 1e16).
  # Use spectral decomposition: B = V D V', invert only well-conditioned directions.
  # Parameters in the null space of B are unidentifiable -- their SEs are set to NA.
  eig <- tryCatch(eigen(B, symmetric = TRUE), error = function(e) NULL)
  if (is.null(eig) || !all(is.finite(eig$values))) {
    warning("Louis (1982): Eigendecomposition of B failed; standard errors cannot be computed.")
    mdl$theta_se   <- rep(NA_real_, npar)
    mdl$info_mat   <- matrix(NA_real_, npar, npar)
    mdl$louis_used <- TRUE
    return(mdl)
  }

  eig_max <- max(eig$values)
  tol     <- nfree * .Machine$double.eps * eig_max
  good    <- eig$values > tol

  cond_num <- if (any(good)) eig_max / min(eig$values[good]) else Inf
  if (cond_num > 1e10)
    warning("Louis (1982): Ill-conditioned complete-data information (cond. number approx. ",
            format(cond_num, scientific = TRUE, digits = 2),
            "); standard errors may be unreliable.")

  if (!any(good)) {
    warning("Louis (1982): B matrix has no positive eigenvalues; standard errors cannot be computed.")
    mdl$theta_se   <- rep(NA_real_, npar)
    mdl$info_mat   <- matrix(NA_real_, npar, npar)
    mdl$louis_used <- TRUE
    return(mdl)
  }

  # Pseudoinverse using well-conditioned eigenvalues only
  V_eig    <- eig$vectors[, good, drop = FALSE]
  d_inv    <- 1 / eig$values[good]
  V_free   <- V_eig %*% (d_inv * t(V_eig))

  # Expand to full dimension via Jacobian (delta method for constrained P entries)
  V_full   <- J %*% V_free %*% t(J)

  # Compute SEs; set to NA for parameters whose variance is non-positive
  se_vec  <- rep(NA_real_, npar)
  d_ii    <- diag(V_full)
  ok      <- is.finite(d_ii) & (d_ii > 0)
  se_vec[ok] <- sqrt(d_ii[ok])

  mdl$Hess        <- Q_Hess
  mdl$theta_se    <- se_vec
  mdl$info_mat    <- V_full
  mdl$louis_used  <- TRUE
  return(mdl)
}

# ============================================================================== #
# Louis (1982) Q-function builders (internal)
# ============================================================================== #

# Dispatcher: build Q(theta) closure for the appropriate model class
.louis_build_Q <- function(mdl, xi) {
  exog <- if (inherits(mdl, "HMmdl")) isTRUE(mdl$exog) else !is.null(mdl$control$Z)
  if (inherits(mdl, "HMmdl"))    return(.louis_Q_HM(mdl, xi, exog))
  if (inherits(mdl, "MSARmdl"))  return(.louis_Q_MSAR(mdl, xi, exog))
  if (inherits(mdl, "MSVARmdl")) return(.louis_Q_MSVAR(mdl, xi, exog))
  stop("Unsupported model class: ", class(mdl)[1])
}

# Q-function for HMmdl (no AR dynamics, k regimes, optional exogenous)
# Q(theta) = sum_t sum_k xi_{t|T}(k) * log N(y_t - mu_k - Z*beta; 0, Sigma_k)
#           + sum_{i,j} n_hat(j,i) * log P(j,i)
.louis_Q_HM <- function(mdl, xi, exog) {
  y     <- as.matrix(mdl$y)
  k     <- mdl$k
  q     <- ncol(y)
  Tsize <- nrow(y)
  msmu  <- mdl$msmu
  msvar <- mdl$msvar
  Zdm   <- NULL; qz <- 0L
  if (exog) {
    Z   <- as.matrix(mdl$Z)
    qz  <- ncol(Z)
    Zdm <- sweep(Z, 2, colMeans(Z))
  }
  # Expected transition counts (exact at EM convergence)
  xi_sum <- colSums(xi[seq_len(Tsize - 1L), , drop = FALSE])
  n_hat  <- mdl$P * matrix(xi_sum, k, k, byrow = TRUE)
  sigN   <- as.integer(q * (q + 1) / 2)

  function(th) {
    mu_end <- q + q * msmu * (k - 1)
    mu_vec <- th[1:mu_end]
    mu_k   <- matrix(0, k, q)
    if (msmu) {
      for (xk in seq_len(k)) mu_k[xk, ] <- mu_vec[((xk-1)*q + 1):(xk*q)]
    } else {
      for (xk in seq_len(k)) mu_k[xk, ] <- mu_vec[1:q]
    }
    pos <- mu_end + 1L
    betaZ_mat <- NULL
    if (exog) {
      betaZ_mat <- matrix(th[pos:(pos + qz*q - 1L)], qz, q)
      pos <- pos + qz*q
    }
    n_sig   <- sigN + sigN * msvar * (k - 1)
    sig_vec <- th[pos:(pos + n_sig - 1L)]
    pos     <- pos + n_sig
    P       <- matrix(th[pos:(pos + k*k - 1L)], k, k)

    Q_obs <- 0
    for (xk in seq_len(k)) {
      eps <- sweep(y, 2, mu_k[xk, ])
      if (exog) eps <- eps - Zdm %*% betaZ_mat
      sig_idx <- if (msvar) ((xk-1)*sigN + 1):(xk*sigN) else 1:sigN
      sigma_k <- covar_unvech(sig_vec[sig_idx], q)
      # Floor eigenvalues to keep sigma PD for numDeriv perturbations
      eig_s <- eigen(sigma_k, symmetric = TRUE)
      eig_s$values <- pmax(eig_s$values, 1e-10)
      sigma_k <- eig_s$vectors %*% diag(eig_s$values, nrow = q) %*% t(eig_s$vectors)
      ch <- tryCatch(chol(sigma_k), error = function(e) NULL)
      if (is.null(ch)) return(-Inf)
      log_det <- 2 * sum(log(diag(ch)))
      eps_std <- backsolve(ch, t(eps), transpose = TRUE)
      maha    <- colSums(eps_std^2)
      Q_obs   <- Q_obs + sum(xi[, xk] * (-0.5 * (q * log(2*pi) + log_det + maha)))
    }
    Q_obs + sum(n_hat * log(pmax(P, 1e-300)))
  }
}

# Q-function for MSARmdl / MSARXmdl (univariate, extended states M = k^(p+1))
# Uses T x M smoothed probabilities from ExpectationM.
# For extended state n = (i_0, ..., i_p), the observation density is:
#   N(y_t; mu_{i_0} + sum_j phi_j*(y_{t-j} - mu_{i_j}) + Z*beta, sigma^2_{i_0})
.louis_Q_MSAR <- function(mdl, xi, exog) {
  y     <- as.vector(mdl$y)
  x     <- as.matrix(mdl$x)
  k     <- mdl$k
  p     <- mdl$p
  M     <- k^(p + 1)
  Tsize <- length(y)
  msmu  <- mdl$msmu
  msvar <- mdl$msvar
  # Extended state grid (expand.grid: first index = current regime, varies fastest)
  lx <- replicate(p + 1, seq_len(k), simplify = FALSE)
  sg <- as.matrix(expand.grid(lx))  # M x (p+1), 1-indexed
  Zdm <- NULL; qz <- 0L
  if (exog) {
    Z_orig <- as.matrix(mdl$control$Z)
    qz     <- ncol(Z_orig)
    Z      <- Z_orig[(p + 1):(Tsize + p), , drop = FALSE]
    Zdm    <- sweep(Z, 2, colMeans(Z))
  }
  # Expected transition counts: rows t = 2..T of the extended smoothed probabilities,
  # the same range as the EM transition-count update
  n_hat <- matrix(0, k, k)
  for (n in seq_len(M)) {
    n_hat[sg[n, 1], sg[n, 2]] <- n_hat[sg[n, 1], sg[n, 2]] + sum(xi[-1L, n])
  }

  function(th) {
    n_mu <- 1L + msmu * (k - 1L)
    mu   <- th[1:n_mu]
    phi  <- th[(n_mu + 1L):(n_mu + p)]
    pos  <- n_mu + p + 1L
    betaZ_v <- NULL
    if (exog) {
      betaZ_v <- th[pos:(pos + qz - 1L)]
      pos <- pos + qz
    }
    n_sig <- 1L + msvar * (k - 1L)
    sig   <- th[pos:(pos + n_sig - 1L)]
    pos   <- pos + n_sig
    P     <- matrix(th[pos:(pos + k*k - 1L)], k, k)
    if (!msmu)  mu  <- rep(mu, k)
    if (!msvar) sig <- rep(sig, k)
    zdm_eff <- if (exog) as.vector(Zdm %*% betaZ_v) else 0

    Q_obs <- 0
    for (n in seq_len(M)) {
      reg  <- sg[n, ]
      z    <- y - mu[reg[1]] - zdm_eff
      xz   <- rep(0, Tsize)
      for (j in seq_len(p)) xz <- xz + phi[j] * (x[, j] - mu[reg[j + 1]])
      eps_n <- z - xz
      s2    <- max(sig[reg[1]], 1e-10)
      Q_obs <- Q_obs + sum(xi[, n] * (-0.5 * (log(2*pi) + log(s2) + eps_n^2/s2)))
    }
    Q_obs + sum(n_hat * log(pmax(P, 1e-300)))
  }
}

# Q-function for MSVARmdl / MSVARXmdl (multivariate, extended states M = k^(ar+1))
# Uses T x M smoothed probabilities from ExpectationM.
# For extended state n = (i_0, ..., i_ar), the observation density is:
#   N(y_t; mu_{i_0} + sum_j A_j*(y_{t-j} - mu_{i_j}) + Z*beta, Sigma_{i_0})
.louis_Q_MSVAR <- function(mdl, xi, exog) {
  y     <- as.matrix(mdl$y)
  x     <- as.matrix(mdl$x)
  k     <- mdl$k
  ar    <- mdl$p
  q     <- ncol(y)
  M     <- k^(ar + 1)
  Tsize <- nrow(y)
  msmu  <- mdl$msmu
  msvar <- mdl$msvar
  lx <- replicate(ar + 1, seq_len(k), simplify = FALSE)
  sg <- as.matrix(expand.grid(lx))
  Zdm <- NULL; qz <- 0L
  if (exog) {
    Z_orig <- as.matrix(mdl$control$Z)
    qz     <- ncol(Z_orig)
    Z      <- Z_orig[(ar + 1):(Tsize + ar), , drop = FALSE]
    Zdm    <- sweep(Z, 2, colMeans(Z))
  }
  # Expected transition counts: rows t = 2..T of the extended smoothed probabilities,
  # the same range as the EM transition-count update
  n_hat <- matrix(0, k, k)
  for (n in seq_len(M)) {
    n_hat[sg[n, 1], sg[n, 2]] <- n_hat[sg[n, 1], sg[n, 2]] + sum(xi[-1L, n])
  }
  sigN <- as.integer(q * (q + 1) / 2)

  function(th) {
    n_mu   <- q + q * msmu * (k - 1)
    mu_vec <- th[1:n_mu]
    mu_k   <- matrix(0, k, q)
    if (msmu) {
      for (xk in seq_len(k)) mu_k[xk, ] <- mu_vec[((xk-1)*q + 1):(xk*q)]
    } else {
      for (xk in seq_len(k)) mu_k[xk, ] <- mu_vec[1:q]
    }
    pos     <- n_mu + 1L
    phi_vec <- th[pos:(pos + q*q*ar - 1L)]
    phi     <- t(matrix(phi_vec, q*ar, q))  # q x (q*ar), matches C++ convention
    pos     <- pos + q*q*ar
    betaZ_mat <- NULL
    if (exog) {
      betaZ_mat <- matrix(th[pos:(pos + qz*q - 1L)], qz, q)
      pos <- pos + qz*q
    }
    n_sig   <- sigN + sigN * msvar * (k - 1)
    sig_vec <- th[pos:(pos + n_sig - 1L)]
    pos     <- pos + n_sig
    P       <- matrix(th[pos:(pos + k*k - 1L)], k, k)
    # Pre-compute Cholesky per regime
    chol_list  <- vector("list", k)
    logdet_vec <- numeric(k)
    for (xk in seq_len(k)) {
      sig_idx <- if (msvar) ((xk-1)*sigN + 1):(xk*sigN) else 1:sigN
      sigma_k <- covar_unvech(sig_vec[sig_idx], q)
      # Floor eigenvalues to keep sigma PD for numDeriv perturbations
      eig_s <- eigen(sigma_k, symmetric = TRUE)
      eig_s$values <- pmax(eig_s$values, 1e-10)
      sigma_k <- eig_s$vectors %*% diag(eig_s$values, nrow = q) %*% t(eig_s$vectors)
      ch <- tryCatch(chol(sigma_k), error = function(e) NULL)
      if (is.null(ch)) return(-Inf)
      chol_list[[xk]]  <- ch
      logdet_vec[xk]   <- 2 * sum(log(diag(ch)))
    }
    zdm_eff <- if (exog) Zdm %*% betaZ_mat else matrix(0, Tsize, q)

    Q_obs <- 0
    for (n in seq_len(M)) {
      reg    <- sg[n, ]
      y_tild <- sweep(y, 2, mu_k[reg[1], ]) - zdm_eff
      xz     <- matrix(0, Tsize, q*ar)
      for (j in seq_len(ar)) {
        cols <- ((j-1)*q + 1):(j*q)
        xz[, cols] <- sweep(x[, cols], 2, mu_k[reg[j+1], ])
      }
      eps_n   <- y_tild - xz %*% t(phi)
      xk      <- reg[1]
      eps_std <- backsolve(chol_list[[xk]], t(eps_n), transpose = TRUE)
      maha    <- colSums(eps_std^2)
      Q_obs   <- Q_obs + sum(xi[, n] * (-0.5 * (q * log(2*pi) + logdet_vec[xk] + maha)))
    }
    Q_obs + sum(n_hat * log(pmax(P, 1e-300)))
  }
}

#' @title Intercept from mu for MSARmdl
#' 
#' @description This function computes the intercept for each regime k for an Markov switching AR model
#'
#' @param mdl List with model properties
#'
#' @return a \code{(k x 1)} vector of intercepts
#' 
#' @keywords internal
#' 
#' @export
interMSARmdl <- function(mdl){
  inter <- matrix(0,mdl$k,1)
  for(xk in 1:mdl$k){
    inter[xk,] <- t(as.matrix(mdl$mu[xk,]))*(1-sum(mdl$phi))
  }
  return(inter)
}

#' @title Intercept from mu for MSVARmdl
#' 
#' @description This function computes the intercept for each regime k for an Markov switching VAR model
#'
#' @param mdl List with model properties
#'
#' @return a \code{(k x q)} vector of intercepts
#' 
#' @keywords internal
#' 
#' @export
interMSVARmdl <- function(mdl){
  inter <- matrix(0,mdl$k,mdl$q)
  for(xk in 1:mdl$k){
    nu_tmp <- (diag(mdl$q*mdl$p)-mdl$Fmat)%*%as.matrix(c(mdl$mu[xk,],rep(0,mdl$q*(mdl$p-1))))
    inter[xk,] <- nu_tmp[(1:(mdl$q))]
  }
  return(inter)
}

# Project an MLE starting vector into the user-supplied box [lo, up].
# nloptr::slsqp() throws a hard error ("at least one element in x0 < lb") when the
# starting value violates the box, which burns a 'maxit_converge' retry per draw.
# The trailing k*k entries (transition matrix P, column-stochastic) are additionally
# projected onto {column sums = 1} within the box by iterative proportional slack
# redistribution (clip-and-renormalise alone can leave the box for k >= 3).
# Returns theta_0 unchanged when both bounds are NULL or malformed.
.mle_project_start <- function(theta_0, lo, up, k){
  n <- length(theta_0)
  if (is.null(lo) && is.null(up)) return(theta_0)
  if (is.null(lo)) lo <- rep(-Inf, n)
  if (is.null(up)) up <- rep(Inf, n)
  if (length(lo) != n || length(up) != n) return(theta_0)  # malformed: let slsqp report it
  th <- pmin(pmax(theta_0, lo), up)
  Pidx <- (n - k*k + 1L):n
  Plo <- matrix(lo[Pidx], k, k); Pup <- matrix(up[Pidx], k, k)
  Pm  <- matrix(th[Pidx], k, k)
  for (j in seq_len(k)){
    lj <- pmax(Plo[, j], 0); uj <- pmin(Pup[, j], 1)
    if (sum(lj) > 1 + 1e-12 || sum(uj) < 1 - 1e-12) next  # infeasible column box: keep clamped values
    pj <- pmin(pmax(Pm[, j], lj), uj)
    for (it in seq_len(3L*k)){
      d <- 1 - sum(pj)
      if (abs(d) < 1e-12) break
      slack <- if (d > 0) uj - pj else pj - lj
      tot <- sum(slack)
      if (tot <= 0) break
      pj <- pmin(pmax(pj + d*slack/tot, lj), uj)
    }
    Pm[, j] <- pj
  }
  th[Pidx] <- as.vector(Pm)
  th
}

#' @title Hidden Markov model maximum likelihood estimation
#'
#' @description This function computes estimate a Hidden Markov model using MLE.
#'
#' @param theta_0 vector containing initial values to use in optimization
#' @param mdl_in List with model properties (can be obtained from estimating linear model i.e., using \code{\link{ARmdl}})
#' @param k integer determining the number of regimes
#' @param optim_options List containing 
#' \itemize{
#'  \item maxit: maximum number of iterations.
#'  \item thtol: convergence criterion.
#'}
#'
#' @return List with model attributes
#' 
#' @keywords internal
#' 
#' @export
HMmdl_mle <- function(theta_0, mdl_in, k, optim_options){
  if (isTRUE(mdl_in$mle_variance_constraint >= 0) && is.null(mdl_in$sigma)) {
    stop("mdl_in$sigma (one-regime fit covariance) is required when mle_variance_constraint >= 0.")
  }
  # ----- Define function environment
  HMmdl_mle_env <- rlang::env()
  # set environment variables
  HMmdl_mle_env$mdl_in <- mdl_in
  HMmdl_mle_env$k <- k
  HMmdl_mle_env$q <- mdl_in$q
  HMmdl_mle_env$msmu <- mdl_in$msmu
  HMmdl_mle_env$msvar <- mdl_in$msvar
  HMmdl_mle_env$var_k1 <- mdl_in$sigma
  HMmdl_mle_env$betaZ <- mdl_in$betaZ
  HMmdl_mle_env$mle_variance_constraint <- mdl_in$mle_variance_constraint
  # ----------- HMmdl_mle equality constraint functions
  loglik_const_eq_HMmdl <- function(theta){
    # ----- Load equality constraint parameters
    k <- get("k", envir = HMmdl_mle_env)
    # ----- constraint
    P = matrix(theta[(length(theta)-k*k+1):(length(theta))],k,k)
    constraint = colSums(P)-1
    return(constraint)
  }
  # ---------- HMmdl_mle inequality constraint functions
  loglik_const_ineq_HMmdl <- function(theta){
    # ----- Load inequality constraint parameters
    k <- get("k", envir = HMmdl_mle_env)
    q <- get("q", envir = HMmdl_mle_env)
    msmu <- get("msmu", envir = HMmdl_mle_env)
    msvar <- get("msvar", envir = HMmdl_mle_env)
    var_k1 <- get("var_k1", envir = HMmdl_mle_env)
    betaZ <- get("betaZ", envir = HMmdl_mle_env)
    mle_variance_constraint <- get("mle_variance_constraint", envir = HMmdl_mle_env)
    # ----- constraint
    # transition probabilities
    Pvec <- theta[(length(theta)-k*k+1):(length(theta))]
    ineq_constraint <- c(Pvec, 1-Pvec)
    if (mle_variance_constraint>=0){
      Nsig <- (q*(q+1))/2
      eigen_vals <- c()
      # theta layout: [mu, betaZ (if exog), vech(sigma) x regimes, vec(P)]; the
      # selector must skip betaZ and cover the full theta length (a short logical
      # mask would be recycled by R into the P block)
      sig <- theta[c(rep(0, q + q*(k-1)*msmu + length(betaZ)), rep(1, Nsig + Nsig*(k-1)*msvar), rep(0, k*k))==1]
      if (msvar==TRUE){
        for (xk in  1:k){
          sig_tmp <- sig[(Nsig*(xk-1)+1):(Nsig*(xk-1)+Nsig)]
          sigma <- covar_unvech(sig_tmp, q)
          eigen_vals <- c(eigen_vals,eigen(sigma)[[1]])
        } 
      }else{
        sigma <- covar_unvech(sig, q)
        eigen_vals <- c(eigen_vals,eigen(sigma)[[1]])
      }
      # relative floor: smallest eigenvalue >= c * tr(Sigma_lin)/q, the same
      # anchoring as em_variance_constraint (var_k1 is the one-regime fit covariance)
      sig_floor <- mle_variance_constraint * sum(diag(as.matrix(var_k1)))/q
      ineq_constraint = c(eigen_vals - sig_floor, ineq_constraint)
    }
    return(-ineq_constraint)
  }
  loglike_wrapper <- function(theta) logLike_HMmdl_min(theta, mdl_in, k)
  # project starting values into the user box (no-op when bounds are NULL)
  theta_0 <- .mle_project_start(theta_0, optim_options$mle_theta_low, optim_options$mle_theta_upp, k)
  # use nloptr optimization to minimize (maximize) likelihood
  res <- nloptr::slsqp(x0 = theta_0,
                       fn = loglike_wrapper,
                       gr = NULL,
                       lower = optim_options$mle_theta_low,
                       upper = optim_options$mle_theta_upp,
                       hin = loglik_const_ineq_HMmdl,
                       heq = loglik_const_eq_HMmdl,
                       nl.info = FALSE,
                       control = list(maxeval = optim_options$maxit, xtol_rel = optim_options$thtol),
                       deprecatedBehavior = FALSE) 
  ok <- isTRUE(res$convergence > 0) && all(is.finite(res$par))
  if (!ok && all(is.finite(res$par))) {
    # slsqp reported failure but returned a finite iterate (e.g. NLOPT_ROUNDOFF_LIMITED):
    # keep it only if it improves on the starting value (objective is -logLike).
    v_res <- suppressWarnings(tryCatch(loglike_wrapper(res$par), error=function(e) Inf))
    v_0   <- suppressWarnings(tryCatch(loglike_wrapper(theta_0), error=function(e) Inf))
    feas  <- all(loglik_const_ineq_HMmdl(res$par) <= 1e-8) &&
             all(abs(loglik_const_eq_HMmdl(res$par)) <= 1e-6)
    if (is.finite(v_res) && v_res < v_0 && feas) ok <- TRUE
  }
  if (!ok) {
    warning("HMmdl MLE optimization failed (convergence=", res$convergence, "). Using initial values.")
    par_use <- theta_0
  } else {
    par_use <- res$par
    if (!isTRUE(res$convergence > 0)) {
      warning("HMmdl MLE optimization reported failure (convergence=", res$convergence,
              ") but its iterate improves on the starting value; keeping it.")
    }
  }
  # keep the transition block a valid column-stochastic matrix (slsqp honors its
  # constraints only to tolerance); conditional, so clean fits are untouched
  npar_use <- length(par_use); Pblk <- (npar_use - k*k + 1):npar_use
  Pm <- matrix(par_use[Pblk], k, k)
  if (any(Pm < 0) || any(Pm > 1) || any(abs(colSums(Pm) - 1) > 1e-8)) {
    Pm[Pm < 0] <- 0; Pm[Pm > 1] <- 1
    cs <- colSums(Pm); cs[cs <= 0] <- 1
    par_use[Pblk] <- as.numeric(sweep(Pm, 2, cs, "/"))
  }
  output <- tryCatch(ExpectationM_HMmdl(par_use, mdl_in, k), error=function(e) NULL)
  if (is.null(output)) return(NULL)
  output$iterations <- res$iter
  output$converged <- isTRUE(res$convergence %in% c(1, 2, 3, 4))
  output$convergence_code <- res$convergence
  output$St <- output$xi_t_T
  return(output)
}

#' @title Markov-switching autoregressive maximum likelihood estimation
#' 
#' @description This function computes estimate a Markov-switching autoregressive model using MLE.
#' 
#' @param theta_0 vector containing initial values to use in optimization
#' @param mdl_in List with model properties (can be obtained from estimating linear model i.e., using \code{\link{ARmdl}})
#' @param k integer determining the number of regimes
#' @param optim_options List containing 
#' \itemize{
#'  \item maxit: maximum number of iterations.
#'  \item thtol: convergence criterion.
#'}
#'
#' @return List with model attributes
#' 
#' @keywords internal
#' 
#' @export
MSARmdl_mle <- function(theta_0, mdl_in, k, optim_options){
  if (isTRUE(mdl_in$mle_variance_constraint >= 0) && is.null(mdl_in$sigma)) {
    stop("mdl_in$sigma (one-regime fit covariance) is required when mle_variance_constraint >= 0.")
  }
  # ----- Define function environment
  MSARmdl_mle_env         <- rlang::env()
  # set environment variables
  MSARmdl_mle_env$p       <- mdl_in$p
  MSARmdl_mle_env$k       <- k
  MSARmdl_mle_env$msmu    <- mdl_in$msmu
  MSARmdl_mle_env$msvar   <- mdl_in$msvar
  MSARmdl_mle_env$var_k1  <- mdl_in$sigma
  MSARmdl_mle_env$betaZ   <- mdl_in$betaZ
  MSARmdl_mle_env$mle_stationary_constraint <- mdl_in$mle_stationary_constraint
  MSARmdl_mle_env$mle_variance_constraint   <- mdl_in$mle_variance_constraint
  # ----------- MSARmdl_mle equality constraint functions
  loglik_const_eq_MSARmdl <- function(theta){
    # ----- Load equality constraint parameters
    k <- get("k", envir = MSARmdl_mle_env)
    # ----- constraint
    P <- matrix(theta[(length(theta)-k*k+1):(length(theta))],k,k)
    constraint <- colSums(P)-1
    return(constraint)
  }
  # ---------- MSARmdl_mle inequality constraint functions
  loglik_const_ineq_MSARmdl <- function(theta){
    # ----- Load inequality constraint parameters
    p       <- get("p", envir = MSARmdl_mle_env)
    k       <- get("k", envir = MSARmdl_mle_env)
    msmu    <- get("msmu", envir = MSARmdl_mle_env)
    msvar   <- get("msvar", envir = MSARmdl_mle_env)
    var_k1  <- get("var_k1", envir = MSARmdl_mle_env)
    betaZ   <- get("betaZ", envir = MSARmdl_mle_env)
    mle_stationary_constraint <- get("mle_stationary_constraint", envir = MSARmdl_mle_env)
    mle_variance_constraint   <- get("mle_variance_constraint", envir = MSARmdl_mle_env)
    # ----- constraint
    # transition probabilities
    Pvec <- theta[(length(theta)-k*k+1):(length(theta))]
    ineq_constraint <- c(Pvec, 1-Pvec)
    if (mle_stationary_constraint==TRUE){
      # roots of characteristic function
      n_p <- 2 + msmu*(k-1)
      poly_fun <- c(1,-theta[n_p:(n_p+p-1)])
      if (any(!is.finite(poly_fun))) return(rep(Inf, length(poly_fun) - 1))
      roots <- Mod(polyroot(poly_fun))
      ineq_constraint = c((roots-1), ineq_constraint)
    }
    if (mle_variance_constraint>=0){
      # variances (lower bound is 1% of single regime variance)
      vars <- theta[c(rep(0, 1+(k-1)*msmu+p+length(betaZ)), rep(1, 1 + (k-1)*msvar), rep(0, k*k))==1] - c(mle_variance_constraint*var_k1)
      ineq_constraint <- c(vars, ineq_constraint)
    }
    return(-ineq_constraint)
  }
  
  loglike_wrapper <- function(theta) {
    k <- get("k", envir = MSARmdl_mle_env)
    mdl_in <- get("mdl_in", envir = MSARmdl_mle_env)
    if (length(mdl_in$betaZ)>0){
      logLike_MSARXmdl_min(theta, mdl_in, k) 
    }else{
      logLike_MSARmdl_min(theta, mdl_in, k)
    }
  }
  # project starting values into the user box (no-op when bounds are NULL)
  theta_0 <- .mle_project_start(theta_0, optim_options$mle_theta_low, optim_options$mle_theta_upp, k)
  # use nloptr optimization to minimize (maximize) likelihood
  res <- nloptr::slsqp(x0 = theta_0,
                       fn = loglike_wrapper,
                       gr = NULL,
                       lower = optim_options$mle_theta_low,
                       upper = optim_options$mle_theta_upp,
                       hin = loglik_const_ineq_MSARmdl,
                       heq = loglik_const_eq_MSARmdl,
                       nl.info = FALSE,
                       control = list(maxeval = optim_options$maxit, xtol_rel = optim_options$thtol),
                       deprecatedBehavior = FALSE) 
  ok <- isTRUE(res$convergence > 0) && all(is.finite(res$par))
  if (!ok && all(is.finite(res$par))) {
    # slsqp reported failure but returned a finite iterate (e.g. NLOPT_ROUNDOFF_LIMITED):
    # keep it only if it improves on the starting value (objective is -logLike).
    v_res <- suppressWarnings(tryCatch(loglike_wrapper(res$par), error=function(e) Inf))
    v_0   <- suppressWarnings(tryCatch(loglike_wrapper(theta_0), error=function(e) Inf))
    feas  <- all(loglik_const_ineq_MSARmdl(res$par) <= 1e-8) &&
             all(abs(loglik_const_eq_MSARmdl(res$par)) <= 1e-6)
    if (is.finite(v_res) && v_res < v_0 && feas) ok <- TRUE
  }
  if (!ok) {
    warning("MSARmdl MLE optimization failed (convergence=", res$convergence, "). Using initial values.")
    par_use <- theta_0
  } else {
    par_use <- res$par
    if (!isTRUE(res$convergence > 0)) {
      warning("MSARmdl MLE optimization reported failure (convergence=", res$convergence,
              ") but its iterate improves on the starting value; keeping it.")
    }
  }
  # keep the transition block a valid column-stochastic matrix (slsqp honors its
  # constraints only to tolerance); conditional, so clean fits are untouched
  npar_use <- length(par_use); Pblk <- (npar_use - k*k + 1):npar_use
  Pm <- matrix(par_use[Pblk], k, k)
  if (any(Pm < 0) || any(Pm > 1) || any(abs(colSums(Pm) - 1) > 1e-8)) {
    Pm[Pm < 0] <- 0; Pm[Pm > 1] <- 1
    cs <- colSums(Pm); cs[cs <= 0] <- 1
    par_use[Pblk] <- as.numeric(sweep(Pm, 2, cs, "/"))
  }
  if (length(mdl_in$betaZ)>0){
    output <- tryCatch(ExpectationM_MSARXmdl(par_use, mdl_in, k),
                       error=function(e) NULL)
  }else{
    output <- tryCatch(ExpectationM_MSARmdl(par_use, mdl_in, k),
                       error=function(e) NULL)
  }
  if (is.null(output)) return(NULL)
  output$iterations <- res$iter
  output$converged <- isTRUE(res$convergence %in% c(1, 2, 3, 4))
  output$convergence_code <- res$convergence
  output$St <- output$xi_t_T
  return(output)
}


#' @title Markov-switching vector autoregressive maximum likelihood estimation
#' 
#' @description This function computes estimate a Markov-switching vector autoregressive model using MLE.
#' 
#' @param theta_0 vector containing initial values to use in optimization
#' @param mdl_in List with model properties (can be obtained from estimating linear model i.e., using \code{\link{VARmdl}})
#' @param k integer determining the number of regimes
#' @param optim_options List containing 
#' \itemize{
#'  \item maxit: maximum number of iterations.
#'  \item thtol: convergence criterion.
#'}
#'
#' @return List with model attributes
#' 
#' @keywords internal
#' 
#' @export
MSVARmdl_mle <- function(theta_0, mdl_in, k, optim_options){
  if (isTRUE(mdl_in$mle_variance_constraint >= 0) && is.null(mdl_in$sigma)) {
    stop("mdl_in$sigma (one-regime fit covariance) is required when mle_variance_constraint >= 0.")
  }
  # ---------- Define function environment
  MSVARmdl_mle_env <- rlang::env()
  # set environment variables
  MSVARmdl_mle_env$mdl_in <- mdl_in
  MSVARmdl_mle_env$p <- mdl_in$p
  MSVARmdl_mle_env$k <- k
  MSVARmdl_mle_env$q <- mdl_in$q
  MSVARmdl_mle_env$msmu <- mdl_in$msmu
  MSVARmdl_mle_env$msvar <- mdl_in$msvar
  MSVARmdl_mle_env$betaZ   <- mdl_in$betaZ
  MSVARmdl_mle_env$mle_stationary_constraint <- mdl_in$mle_stationary_constraint
  MSVARmdl_mle_env$mle_variance_constraint <- mdl_in$mle_variance_constraint
  MSVARmdl_mle_env$var_k1 <- mdl_in$sigma
  # ---------- MSVARmdl_mle equality constraint functions
  loglik_const_eq_MSVARmdl <- function(theta){
    # ----- Load equality constraint parameters
    k <- get("k", envir = MSVARmdl_mle_env)
    # ----- constraint
    P <- matrix(theta[(length(theta)-k*k+1):(length(theta))],k,k)
    constraint <- colSums(P)-1
    return(constraint)
  }
  # ---------- MSVARmdl_mle inequality constraint functions
  loglik_const_ineq_MSVARmdl <- function(theta){
    # ----- Load inequality constraint parameters
    p <- get("p", envir = MSVARmdl_mle_env)
    k <- get("k", envir = MSVARmdl_mle_env)
    q <- get("q", envir = MSVARmdl_mle_env)
    msmu <- get("msmu", envir = MSVARmdl_mle_env)
    msvar <- get("msvar", envir = MSVARmdl_mle_env)
    betaZ   <- get("betaZ", envir = MSVARmdl_mle_env)
    mle_stationary_constraint <- get("mle_stationary_constraint", envir = MSVARmdl_mle_env)
    mle_variance_constraint <- get("mle_variance_constraint", envir = MSVARmdl_mle_env)
    var_k1 <- get("var_k1", envir = MSVARmdl_mle_env)
    Nsig <- (q*(q+1))/2
    phi_len <- q*p*q
    # ----- constraint 
    Pvec = theta[(length(theta)-k*k+1):(length(theta))]
    ineq_constraint = c(Pvec, 1-Pvec)
    if (mle_stationary_constraint==TRUE){
      # eigen values of companion matrix: stationarity requires all |eigenvalue| <= 1.
      # (The previous form c(|m-1|, 1-|m-1|) only bounded |eigenvalue| to [0, 2],
      # admitting explosive solutions.)
      phi <- t(matrix(theta[c(rep(0, q + q*(k-1)*msmu), rep(1, phi_len), rep(0, length(theta) - (q + q*(k-1)*msmu)-phi_len))==1], q*p, q))
      Fmat <- companionMat(phi,p,q)
      eig_mod <- Mod(eigen(Fmat)[[1]])
      eig_mod[!is.finite(eig_mod)] <- 1e6  # finite penalty keeps slsqp's FD jacobian defined
      ineq_constraint = c(1-eig_mod, ineq_constraint)
    }
    if (mle_variance_constraint>=0){
      # eigen values of covariance matrices
      eigen_vals <- c()
      sig <- theta[c(rep(0, q + q*(k-1)*msmu + phi_len + length(betaZ)), rep(1, Nsig + Nsig*(k-1)*msvar), rep(0, k*k))==1]
      if (msvar==TRUE){
        for (xk in  1:k){
          sig_tmp <- sig[(Nsig*(xk-1)+1):(Nsig*(xk-1)+Nsig)]
          sigma <- covar_unvech(sig_tmp, q)
          eigen_vals <- c(eigen_vals,eigen(sigma)[[1]])
        } 
      }else{
        sigma <- covar_unvech(sig, q)
        eigen_vals <- c(eigen_vals,eigen(sigma)[[1]])
      }
      # relative floor: smallest eigenvalue >= c * tr(Sigma_lin)/q, the same
      # anchoring as em_variance_constraint (var_k1 is the one-regime fit covariance)
      sig_floor <- mle_variance_constraint * sum(diag(as.matrix(var_k1)))/q
      ineq_constraint = c(eigen_vals - sig_floor, ineq_constraint)
    }
    return(-ineq_constraint)
  }
  # Wrapper for log-likelihood with additional arguments
  loglike_wrapper <- function(theta) {
    k <- get("k", envir = MSVARmdl_mle_env)
    mdl_in <- get("mdl_in", envir = MSVARmdl_mle_env)
    if (length(mdl_in$betaZ)>0){
      logLike_MSVARXmdl_min(theta, mdl_in, k) 
    }else{
      logLike_MSVARmdl_min(theta, mdl_in, k)
    }
  }
  # project starting values into the user box (no-op when bounds are NULL)
  theta_0 <- .mle_project_start(theta_0, optim_options$mle_theta_low, optim_options$mle_theta_upp, k)
  # use nloptr optimization to minimize (maximize) likelihood
  res <- nloptr::slsqp(x0 = theta_0,
                       fn = loglike_wrapper,
                       gr = NULL,
                       lower = optim_options$mle_theta_low,
                       upper = optim_options$mle_theta_upp,
                       hin = loglik_const_ineq_MSVARmdl,
                       heq = loglik_const_eq_MSVARmdl,
                       nl.info = FALSE,
                       control = list(maxeval = optim_options$maxit, xtol_rel = optim_options$thtol),
                       deprecatedBehavior = FALSE) 
  ok <- isTRUE(res$convergence > 0) && all(is.finite(res$par))
  if (!ok && all(is.finite(res$par))) {
    # slsqp reported failure but returned a finite iterate (e.g. NLOPT_ROUNDOFF_LIMITED):
    # keep it only if it improves on the starting value (objective is -logLike).
    v_res <- suppressWarnings(tryCatch(loglike_wrapper(res$par), error=function(e) Inf))
    v_0   <- suppressWarnings(tryCatch(loglike_wrapper(theta_0), error=function(e) Inf))
    feas  <- all(loglik_const_ineq_MSVARmdl(res$par) <= 1e-8) &&
             all(abs(loglik_const_eq_MSVARmdl(res$par)) <= 1e-6)
    if (is.finite(v_res) && v_res < v_0 && feas) ok <- TRUE
  }
  if (!ok) {
    warning("MSVARmdl MLE optimization failed (convergence=", res$convergence, "). Using initial values.")
    par_use <- theta_0
  } else {
    par_use <- res$par
    if (!isTRUE(res$convergence > 0)) {
      warning("MSVARmdl MLE optimization reported failure (convergence=", res$convergence,
              ") but its iterate improves on the starting value; keeping it.")
    }
  }
  # keep the transition block a valid column-stochastic matrix (slsqp honors its
  # constraints only to tolerance); conditional, so clean fits are untouched
  npar_use <- length(par_use); Pblk <- (npar_use - k*k + 1):npar_use
  Pm <- matrix(par_use[Pblk], k, k)
  if (any(Pm < 0) || any(Pm > 1) || any(abs(colSums(Pm) - 1) > 1e-8)) {
    Pm[Pm < 0] <- 0; Pm[Pm > 1] <- 1
    cs <- colSums(Pm); cs[cs <= 0] <- 1
    par_use[Pblk] <- as.numeric(sweep(Pm, 2, cs, "/"))
  }
  if (length(mdl_in$betaZ)>0){
    output <- tryCatch(ExpectationM_MSVARXmdl(par_use, mdl_in, k), error=function(e) NULL)
  }else{
    output <- tryCatch(ExpectationM_MSVARmdl(par_use, mdl_in, k), error=function(e) NULL)
  }
  if (is.null(output)) return(NULL)
  output$iterations <- res$iter
  output$converged <- isTRUE(res$convergence %in% c(1, 2, 3, 4))
  output$convergence_code <- res$convergence
  output$St <- output$xi_t_T
  return(output)
}



# ============================================================================== #
# Warm-start embeddings for Monte Carlo LR tests (regime-duplication)
# ============================================================================== #

# Enumerate all ways to distribute `extra` clones among k0 regimes
# (weak compositions of `extra` into k0 non-negative parts).
.ws_compositions <- function(k0, extra){
  if (k0 == 1) return(list(extra))
  out <- list()
  for (e1 in 0:extra){
    for (rest in .ws_compositions(k0 - 1, extra - e1)){
      out[[length(out) + 1]] <- c(e1, rest)
    }
  }
  return(out)
}

# Half-vectorization matching the C++ covar_vech() convention
# (row-wise upper triangle == column-major lower triangle for symmetric M).
.ws_vech <- function(M){
  M <- as.matrix(M)
  return(M[lower.tri(M, diag = TRUE)])
}

# Extract regime-level parameters from a fitted model object (any estimMdl class)
# in normalized shapes: mu0 (k0 x q), sig0 (list of k0 q x q), phi0, betaZ0, P0.
.ws_extract <- function(mdl_h0){
  k0 <- mdl_h0$k
  if (is.null(k0)) k0 <- 1
  q  <- mdl_h0$q
  p  <- mdl_h0$p
  if (is.null(p)) p <- 0
  exog <- !is.null(mdl_h0$betaZ)
  # ----- means (k0 x q; MS outputs duplicate rows when msmu = FALSE)
  mu_m <- as.matrix(mdl_h0$mu)
  if (nrow(mu_m) != k0) mu_m <- matrix(as.numeric(mu_m), k0, q, byrow = TRUE)
  # ----- covariances (list of k0 q x q; variances for q = 1)
  sig_raw <- mdl_h0$sigma
  if (is.list(sig_raw)){
    sig0 <- lapply(sig_raw, as.matrix)
  }else{
    sig_m <- as.matrix(sig_raw)
    if (q == 1){
      sig0 <- lapply(seq_len(k0), function(j) matrix(sig_m[min(j, nrow(sig_m)), 1], 1, 1))
    }else{
      sig0 <- replicate(k0, sig_m, simplify = FALSE)
    }
  }
  # ----- autoregressive coefficients
  phi0 <- NULL
  if (p > 0){
    phi0 <- if (q == 1) as.numeric(mdl_h0$phi) else as.matrix(mdl_h0$phi)
  }
  # ----- exogenous regressor coefficients (qz x q)
  betaZ0 <- if (exog) as.matrix(mdl_h0$betaZ) else NULL
  # ----- transition matrix
  P0 <- if (k0 > 1) as.matrix(mdl_h0$P) else NULL
  return(list(k0 = k0, q = q, p = p, exog = exog, mu0 = mu_m, sig0 = sig0,
              phi0 = phi0, betaZ0 = betaZ0, P0 = P0))
}

# Build one embedded (and optionally perturbed) k1-regime theta vector from a
# fitted k0-regime model, for the clone-count composition `extras` (length k0,
# extras[j] = number of extra clones of regime j; sum(extras) = k1 - k0).
.ws_build_theta <- function(src, k1, msmu, msvar, c_pert, extras){
  k0 <- src$k0; q <- src$q; p <- src$p
  nj <- 1 + extras
  # ----- per-regime split weights (equal; asymmetric if nothing else can
  # break the clone symmetry because neither mean nor variance switches)
  wts <- lapply(seq_len(k0), function(j){
    n <- nj[j]
    if (!msmu && !msvar && n > 1) return((n:1) / sum(n:1))
    return(rep(1 / n, n))
  })
  # ----- clone-level mean and covariance with antisymmetric perturbation
  mu_new  <- matrix(0, k1, q)
  sig_new <- vector("list", k1)
  r <- 0
  for (j in seq_len(k0)){
    n  <- nj[j]
    sdj <- sqrt(pmax(diag(src$sig0[[j]]), 0))
    for (a in seq_len(n)){
      r <- r + 1
      s <- if (n == 1) 0 else c_pert * (2 * a - n - 1) / (n - 1)
      mu_new[r, ]   <- src$mu0[j, ] + if (msmu) s * sdj else 0
      sig_new[[r]]  <- src$sig0[[j]] * (if (msvar) (1 + c_pert)^(if (n == 1) 0 else (2 * a - n - 1) / (n - 1)) else 1)
    }
  }
  # ----- transition matrix embedding
  if (k0 == 1){
    # identical clone emissions make the likelihood invariant to P; use a
    # persistent chain as the EM starting point
    P_new <- matrix(0.2 / max(1, k1 - 1), k1, k1)
    diag(P_new) <- 0.8
  }else{
    # lumpability-exact expansion: P1[(i,a),(c,b)] = w_{i,a} * P0[i,c]
    Rm <- matrix(0, k1, k0)
    Cm <- matrix(0, k0, k1)
    r <- 0
    for (j in seq_len(k0)){
      for (a in seq_len(nj[j])){
        r <- r + 1
        Rm[r, j] <- wts[[j]][a]
        Cm[j, r] <- 1
      }
    }
    P_new <- Rm %*% src$P0 %*% Cm
  }
  # ----- assemble theta in the target class layout:
  # c(mu, phi, betaZ, sigma, vec(P)) with empty blocks dropped
  mu_block <- if (msmu){
    if (q == 1) mu_new[, 1] else as.vector(t(mu_new))
  }else{
    as.numeric(src$mu0[1, ])
  }
  phi_block <- if (p > 0){
    if (q == 1) src$phi0 else as.vector(t(src$phi0))
  }else{
    numeric(0)
  }
  x_block <- if (src$exog) as.vector(src$betaZ0) else numeric(0)
  sig_block <- if (q == 1){
    if (msvar) vapply(sig_new, function(m) m[1, 1], 0.0) else src$sig0[[1]][1, 1]
  }else{
    if (msvar) unlist(lapply(sig_new, .ws_vech)) else .ws_vech(src$sig0[[1]])
  }
  return(c(mu_block, phi_block, x_block, sig_block, as.vector(P_new)))
}

# Unperturbed embedding (likelihood exactly equal to the k0-regime fit's;
# all extra clones assigned to `split_regime`). Used by tests and as the
# theoretical basis of the LR >= 0 floor.
.embed_theta <- function(mdl_h0, k1, msmu = NULL, msvar = NULL, split_regime = 1){
  src <- .ws_extract(mdl_h0)
  if (is.null(msmu))  msmu  <- if (!is.null(mdl_h0$msmu))  mdl_h0$msmu  else TRUE
  if (is.null(msvar)) msvar <- if (!is.null(mdl_h0$msvar)) mdl_h0$msvar else TRUE
  extras <- rep(0, src$k0)
  extras[split_regime] <- k1 - src$k0
  return(.ws_build_theta(src, k1, msmu, msvar, c_pert = 0, extras = extras))
}

#' @title Warm-start values for the alternative model of a Monte Carlo LR test
#'
#' @description Builds deterministic starting values for the \code{k1}-regime
#' (alternative) model by embedding a fitted \code{k0}-regime (null) model into
#' the larger parameter space. Each embedding duplicates one or more regimes of
#' the null fit (with the transition-matrix column split so that the likelihood
#' at the exact embedding equals the null model's log-likelihood) and then
#' applies a small deterministic antisymmetric perturbation to the duplicated
#' regimes (means shifted by \code{+/- c_pert} regime standard deviations when
#' the mean switches; variances scaled by \code{(1 + c_pert)^{+/-1}} when the
#' variance switches) so that the EM algorithm can leave the symmetric fixed
#' point. One starting vector is returned per way of distributing the
#' \code{k1 - k0} extra regimes among the \code{k0} regimes of the null fit.
#' These starts are used in addition to the random starting values when
#' \code{init_method = "warmstart"} in \code{\link{LMCLRTest}} and
#' \code{\link{MMCLRTest}}.
#'
#' @param mdl_h0 Fitted model under the null hypothesis (any model object
#' returned by \code{\link{estimMdl}}; a one-regime linear model or a Markov
#' switching model with \code{k0 > 1} regimes).
#' @param k1 integer specifying the number of regimes of the alternative model
#' (must exceed the number of regimes of \code{mdl_h0}).
#' @param msmu Boolean indicating whether the mean switches in the alternative
#' model. Default is \code{NULL}, in which case it is taken from \code{mdl_h0}
#' when available and \code{TRUE} otherwise.
#' @param msvar Boolean indicating whether the variance switches in the
#' alternative model. Default is \code{NULL} (as for \code{msmu}).
#' @param c_pert Double determining the size of the antisymmetric perturbation
#' in regime standard-deviation units. Default is \code{0.5}.
#'
#' @return List of numeric vectors, each a valid starting value for the
#' \code{k1}-regime model (in the parameter ordering used by the corresponding
#' model constructor). Starting vectors containing non-finite values are
#' dropped; the list may be empty if the null fit is degenerate.
#'
#' @references Kasahara, H. & Shimotsu, K. 2018. "Testing the number of regimes in Markov regime switching models." arXiv:1801.06862.
#' @references Dufour, J.-M. 2006. "Monte Carlo tests with nuisance parameters: A general approach to finite-sample inference and nonstandard asymptotics." \emph{Journal of Econometrics}, 133(2), 443-477.
#' @references Rodriguez-Rondon, G. & Dufour, J.-M. 2026. "Monte Carlo Likelihood-Ratio Tests for Markov Switching Models." Bank of Canada Staff Working Paper 2026-23.
#'
#' @export
warmstart_theta <- function(mdl_h0, k1, msmu = NULL, msvar = NULL, c_pert = 0.5){
  src <- .ws_extract(mdl_h0)
  if (k1 <= src$k0) stop("'k1' must exceed the number of regimes of 'mdl_h0'.")
  if (is.null(msmu))  msmu  <- if (!is.null(mdl_h0$msmu))  mdl_h0$msmu  else TRUE
  if (is.null(msvar)) msvar <- if (!is.null(mdl_h0$msvar)) mdl_h0$msvar else TRUE
  comps <- .ws_compositions(src$k0, k1 - src$k0)
  out <- lapply(comps, function(extras) .ws_build_theta(src, k1, msmu, msvar, c_pert, extras))
  out <- out[vapply(out, function(th) all(is.finite(th)), TRUE)]
  return(out)
}


# ============================================================================== #
# EM regime-variance floor helpers (constrained ML; see em_variance_constraint)
# ============================================================================== #

# Per-regime indicator of whether the fitted variance (min eigenvalue for
# covariance matrices) sits at the floor (within relative tolerance 1e-6).
.sigma_floor_binding <- function(sigma, floor_val){
  k <- if (is.list(sigma)) length(sigma) else nrow(as.matrix(sigma))
  if (!isTRUE(floor_val > 0)) return(rep(FALSE, k))
  tol <- floor_val * (1 + 1e-6)
  if (is.list(sigma)){
    vapply(sigma, function(S){
      min(eigen(as.matrix(S), symmetric = TRUE, only.values = TRUE)$values) <= tol
    }, TRUE)
  }else{
    as.numeric(as.matrix(sigma)[, 1]) <= tol
  }
}

# Set the theta-SE entries of floored regimes' variance components to NA
# (boundary point: the asymptotic normal approximation is invalid there).
.floor_na_se <- function(mdl){
  if (is.null(mdl$theta_se) || is.null(mdl$sigma_floor_binding) ||
      !any(mdl$sigma_floor_binding, na.rm = TRUE)) return(mdl)
  ind <- which(as.numeric(mdl$theta_sig_ind) == 1)
  if (length(ind) == 0) return(mdl)
  if (isTRUE(mdl$msvar)){
    bl <- length(ind) / mdl$k
    for (j in which(mdl$sigma_floor_binding)){
      mdl$theta_se[ind[((j - 1) * bl + 1):(j * bl)]] <- NA
    }
  }else{
    # single shared variance block
    mdl$theta_se[ind] <- NA
  }
  return(mdl)
}
