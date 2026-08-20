# ---------------------------------------------------------------------------
# Phase 3: canonical regime labeling.
#
# A fitted model's regime labels are arbitrary: estimation is invariant to
# permuting them, and which permutation a given multistart run lands on depends
# on which local optimum it found. That is harmless whenever a simulation burns
# in (the chain forgets which state it started in), but the boundary-augmented
# test's transient-start convention (burnin = 0) is decisive: if the fitted
# state 1 happens to be the absorbing regime, every simulated null draw starts
# absorbed and never visits the transient regime. This canonicalizes state 1 to
# the LEAST persistent (transient) regime and the last state to the MOST
# persistent, so that convention is well defined. Never called automatically --
# an automatic relabeling would reorder St, theta and P for every existing user
# and every stored result, for no benefit outside this specific use.
# ---------------------------------------------------------------------------

# Index permutation of a FULL theta vector implementing a regime relabeling.
# theta's shared layout (see .theta_len_ms) is mu | phi | betaZ | sigma | vec(P):
# mu and sigma are permuted regime-block-wise (only if they switch at all -- a
# regime-invariant mu or sigma has no regime axis to permute); phi and betaZ are
# unchanged (never regime-indexed); P is permuted on BOTH the row and column
# index, since a transition probability names two regimes, not one.
.canon_theta_perm_idx <- function(k, q, p, qz, msmu, msvar, perm){
  msmu  <- isTRUE(as.logical(msmu))
  msvar <- isTRUE(as.logical(msvar))
  mu_len    <- q * (1 + msmu * (k - 1))
  phi_len   <- q * q * p
  betaZ_len <- qz * q
  Nsig      <- q * (q + 1) / 2
  sig_len   <- Nsig * (1 + msvar * (k - 1))
  P_len     <- k * k
  npar <- mu_len + phi_len + betaZ_len + sig_len + P_len
  idx <- seq_len(npar)
  off <- 0L
  if (msmu){
    for (j in seq_len(k)){
      new_pos <- (off + (j - 1L) * q + 1L):(off + j * q)
      old_pos <- (off + (perm[j] - 1L) * q + 1L):(off + perm[j] * q)
      idx[new_pos] <- old_pos
    }
  }
  off <- off + mu_len + phi_len + betaZ_len
  if (msvar){
    for (j in seq_len(k)){
      new_pos <- (off + (j - 1L) * Nsig + 1L):(off + j * Nsig)
      old_pos <- (off + (perm[j] - 1L) * Nsig + 1L):(off + perm[j] * Nsig)
      idx[new_pos] <- old_pos
    }
  }
  off <- off + sig_len
  for (j in seq_len(k)) for (i in seq_len(k)){
    new_pos <- off + (j - 1L) * k + i
    old_pos <- off + (perm[j] - 1L) * k + perm[i]
    idx[new_pos] <- old_pos
  }
  idx
}

# A scalar-per-regime proxy for the sigma tie-break in the ordering key (D3.1):
# the (1,1) entry of each regime's covariance, which exists identically whether
# sigma is stored as a list of q x q matrices (multivariate) or a k x 1 matrix
# of variances (univariate / MSAR-type classes).
.canon_sigma_key <- function(sigma, k){
  if (is.list(sigma)) vapply(sigma, function(s) s[1, 1], numeric(1)) else as.numeric(sigma)[seq_len(k)]
}

# Permute a regime-indexed field by 'perm', dispatching on its actual shape
# rather than on the model class: sigma/stdev/beta are a list of k matrices for
# every multivariate class but a plain k-row (sigma/stdev) or k-column (beta)
# matrix for the univariate ones, so branching on is.list()/is.matrix() is the
# robust rule (matches the corrected Phase 3 field checklist).
.canon_permute_field <- function(x, perm, along){
  if (is.null(x)) return(x)
  if (is.list(x) && !is.data.frame(x)) return(x[perm])
  if (along == "rows")    return(x[perm, , drop = FALSE])
  if (along == "cols")    return(x[, perm, drop = FALSE])
  if (along == "both")    return(x[perm, perm, drop = FALSE])
  if (along == "vector")  return(x[perm])
  stop("canonicalize_regimes: internal error, unknown 'along' = ", along)
}

#' @title Canonical regime labeling
#'
#' @description Relabels the regimes of a fitted Markov-switching model in
#'   ascending order of fitted persistence (\code{P[j,j]}), so that state 1 is
#'   the least persistent (transient) regime and the last state is the most
#'   persistent. Estimation is invariant to regime relabeling, so a fitted
#'   model's labels are otherwise arbitrary; this makes them deterministic.
#'   Never called automatically by any estimation or testing function.
#'
#' @param mdl A fitted \code{HMmdl}, \code{MSARmdl} or \code{MSVARmdl} object
#'   (the exogenous variants \code{MSARXmdl}/\code{MSVARXmdl} return one of
#'   these base classes) with \code{k > 1} regimes.
#' @param by Character, the ordering rule. Currently only \code{"persistence"}
#'   (ascending fitted \code{P[j,j]}) is implemented.
#'
#' @return An object of the same class as \code{mdl}, with every field on the
#'   regime axis permuted to the canonical order (\code{theta}, \code{mu},
#'   \code{sigma}, \code{stdev}, \code{P}, \code{pinf}, \code{St},
#'   \code{intercept}, \code{beta}, standard errors when present, and the
#'   constraint-binding diagnostics). \code{logLike}, \code{AIC}, \code{BIC},
#'   \code{resid} and \code{fitted} are carried unchanged, not recomputed: a
#'   relabeling does not change the likelihood (verified to 2.3e-13 by
#'   construction of the embedding this equality also underwrites). Fields with
#'   no regime axis (\code{phi}, \code{betaZ}, \code{y}, \code{Z}, \code{n},
#'   \code{p}, \code{q}, \code{k}, \code{control}, and the \code{theta_*_ind}
#'   position maps) are carried unchanged. \code{theta_0}, \code{init_used} and
#'   \code{trace} describe the ORIGINAL labels and are left as-is.
#'
#' @references Rodriguez-Rondon, G., & Dufour, J.-M. 2026a. "Monte Carlo Likelihood-Ratio Tests for Markov Switching Models." \emph{Bank of Canada Staff Working Paper}, No. 2026-23. doi: 10.34989/swp-2026-23.
#'
#' @export
canonicalize_regimes <- function(mdl, by = "persistence"){
  by <- match.arg(by, "persistence")
  cls <- class(mdl)[1]
  if (!(cls %in% c("HMmdl", "MSARmdl", "MSVARmdl")) || is.null(mdl$k) || mdl$k <= 1){
    stop("canonicalize_regimes: 'mdl' must be a fitted HMmdl, MSARmdl or MSVARmdl object (a k > 1 Markov-switching fit); calling it on a one-regime or non-switching-model object is not meaningful.")
  }
  k <- mdl$k; q <- mdl$q; p <- mdl$p
  qz <- if (is.null(mdl$Z)) 0L else NCOL(mdl$Z)
  npar <- length(mdl$theta)
  # Read the ordering key -- and every subsequent structural quantity -- from
  # theta, not mdl$P: mdledit() products (e.g. an MMC mdl_h0_mmc) can carry a
  # stale P that theta itself does not.
  P_ind  <- as.numeric(mdl$theta_P_ind)
  P_th   <- matrix(as.numeric(mdl$theta)[P_ind == 1], k, k)
  if (!is.null(mdl$P) && !isTRUE(all.equal(P_th, mdl$P, tolerance = 1e-8, check.attributes = FALSE))){
    warning("canonicalize_regimes: mdl$P disagrees with the transition matrix encoded in mdl$theta; using theta's, which is canonical.")
  }
  pjj <- diag(P_th)
  mu_mat  <- .theta_mu_kq(mdl$theta, k, q, mdl$msmu)
  sig_key <- .canon_sigma_key(mdl$sigma, k)
  # D3.1: a full deterministic order -- persistence, then each series' mean
  # (lexicographically), then a sigma proxy, then the original index -- so the
  # face search set stays deterministic even at exact ties.
  ord_args <- c(list(pjj), lapply(seq_len(ncol(mu_mat)), function(j) mu_mat[, j]),
                list(sig_key), list(seq_len(k)))
  perm <- do.call(order, ord_args)
  # No early return for an already-canonical fit (perm == identity): every
  # permutation below is a pure reindexing, so applying it with the identity
  # is itself an exact no-op, and running the same code path unconditionally
  # means the returned object's field SHAPES never depend on whether anything
  # actually moved (e.g. pinf is always returned as a plain vector, not a
  # matrix on some calls and a vector on others).
  theta_idx <- .canon_theta_perm_idx(k, q, p, qz, mdl$msmu, mdl$msvar, perm)
  theta_new <- as.numeric(mdl$theta)[theta_idx]
  names(theta_new) <- names(mdl$theta)   # names describe the SLOT (regime-1's mean, ...), not the value

  # theta_0, init_used and trace describe the pre-canonicalization fit and are
  # deliberately left as mdl's own values below (see @return).
  out <- mdl
  out$theta <- theta_new
  out$mu        <- .canon_permute_field(mdl$mu, perm, "rows")
  out$intercept <- .canon_permute_field(mdl$intercept, perm, "rows")
  # 'along' only governs the non-list branch; a list argument is always
  # detected and reordered as a list regardless of what is passed here.
  out$sigma     <- .canon_permute_field(mdl$sigma, perm, "rows")
  out$stdev     <- .canon_permute_field(mdl$stdev, perm, "rows")
  out$beta      <- .canon_permute_field(mdl$beta, perm, "cols")
  out$P    <- mdl$P[perm, perm, drop = FALSE]
  out$pinf <- as.numeric(mdl$pinf)[perm]
  out$St   <- mdl$St[, perm, drop = FALSE]
  if (!is.null(mdl$sigma_floor_binding)) out$sigma_floor_binding <- as.logical(mdl$sigma_floor_binding)[perm]
  if (!is.null(mdl$P_constraint_binding)) out$P_constraint_binding <- mdl$P_constraint_binding[perm, perm, drop = FALSE]
  if (!is.null(mdl$theta_se)){
    v <- as.numeric(mdl$theta_se)[theta_idx]; names(v) <- names(mdl$theta_se)
    out$theta_se <- v
  }
  if (!is.null(mdl$info_mat)){
    out$info_mat <- mdl$info_mat[theta_idx, theta_idx, drop = FALSE]
    dimnames(out$info_mat) <- dimnames(mdl$info_mat)
  }
  if (!is.null(mdl$Hess)){
    # M3-1: Hess is stored in the REDUCED free-P-entry coordinates of the
    # standard-error parameterization (Section 20), not in full theta
    # coordinates -- e.g. 7x7 against length(theta) = 9 for an MSARmdl k = 2
    # fit -- so "permute rows and columns by theta_idx" is dimensionally
    # impossible, and the drop set of that reduction is not itself
    # permutation-equivariant. The tangent space span(J) IS permutation
    # invariant, so Pi J = J A for A = solve(t(J)J) t(J) Pi J (the projection
    # of the permuted tangent directions back onto that span), and the
    # curvature transforms by the congruence A' Hess A -- independent of the
    # reduced parameterization's additive constant, since a Hessian is
    # insensitive to affine shifts of where it is evaluated.
    setup <- .reduced_P_setup(k, npar)
    J <- setup$J
    Pi_full <- diag(npar)[, theta_idx]   # (Pi_full %*% v)[i] == v[theta_idx[i]]
    A <- solve(t(J) %*% J) %*% t(J) %*% Pi_full %*% J
    out$Hess <- t(A) %*% mdl$Hess %*% A
  }
  out
}
