# Phase 3: canonicalize_regimes() relabels a fitted model's regimes in ascending
# persistence order. Estimation is label-permutation invariant, and the simulators
# always start a chain in state 1, so a burnin = 0 simulation (the transient-start
# convention the boundary-augmented test needs) is decisive: if the fitted state 1
# happens to be the absorbing regime, every simulated null draw starts absorbed and
# never visits the transient regime. This is a manual test, never automatic.

.canon_swap_fit <- function(seed, near_absorbing = FALSE){
  set.seed(seed)
  P_true <- if (near_absorbing) cbind(c(0.998, 0.002), c(0.2, 0.8)) else cbind(c(0.85, 0.15), c(0.8, 0.2))
  S <- simuHMM(list(n = 400, q = 2, k = 2, mu = rbind(c(0, 5), c(4, -3)),
                    sigma = list(diag(2), 2 * diag(2)), P = P_true))
  HMmdl(S$y, k = 2, control = list(msmu = TRUE, msvar = TRUE, getSE = FALSE, use_diff_init = 5))
}

test_that("k = 1 and non-Markov-switching objects are rejected", {
  y <- matrix(rnorm(100), 100, 1)
  expect_error(canonicalize_regimes(ARmdl(y, p = 1)), "HMmdl, MSARmdl or MSVARmdl")
  expect_error(canonicalize_regimes(list(k = 1)), "HMmdl, MSARmdl or MSVARmdl")
  expect_error(canonicalize_regimes(structure(list(k = 2), class = "not_a_real_class")),
              "HMmdl, MSARmdl or MSVARmdl")
})

test_that("an already-canonical fit has unchanged VALUES, and the input is not mutated", {
  # Values are unchanged (identity permutation is a no-op on the numbers). Every
  # field's SHAPE is also unconditional on whether anything actually moved --
  # pinf is a (k x 1) matrix on both an already-canonical and a genuinely
  # relabeled call, never a matrix on one path and a flattened vector on the
  # other -- see the code comment on the removed early-return in
  # canonicalize_regimes for why that uniformity matters.
  skip_on_cran()
  m <- .canon_swap_fit(1)   # confirmed already ascending-persistence at this seed
  expect_true(diff(diag(m$P)) >= 0)
  P_before <- m$P
  seed_before <- .Random.seed
  m2 <- canonicalize_regimes(m)
  expect_identical(.Random.seed, seed_before)   # no RNG side effect
  expect_equal(as.numeric(m2$theta), as.numeric(m$theta))
  expect_identical(m$P, P_before)   # m itself untouched (not mutated by reference)
  expect_true(is.matrix(m2$pinf))   # the exact shape the removed early-return got wrong
  expect_equal(dim(m2$pinf), dim(m$pinf))
})

test_that("every field is permuted (or invariant) correctly on a real, non-canonical fit", {
  skip_on_cran()
  m <- .canon_swap_fit(5)   # confirmed not already ascending-persistence at this seed
  expect_true(diag(m$P)[1] > diag(m$P)[2])
  m2 <- canonicalize_regimes(m)
  perm <- order(diag(m$P))   # the permutation this fit must have used

  expect_equal(diag(m2$P), sort(diag(m$P)))
  expect_equal(m2$P, m$P[perm, perm])
  expect_equal(m2$mu, m$mu[perm, , drop = FALSE])
  expect_equal(m2$sigma, m$sigma[perm])
  expect_equal(m2$stdev, m$stdev[perm])
  expect_equal(m2$intercept, m$intercept[perm, , drop = FALSE])
  expect_equal(m2$beta, m$beta[perm])
  expect_equal(m2$pinf, m$pinf[perm, , drop = FALSE])
  expect_true(is.matrix(m2$pinf))   # shape preserved, not flattened
  expect_equal(m2$St, m$St[, perm, drop = FALSE])

  # Carried unchanged: the likelihood-affecting and bookkeeping fields.
  for (f in c("logLike", "AIC", "BIC", "resid", "fitted", "y", "Z", "n", "p", "q", "k",
             "msmu", "msvar", "control", "theta_mu_ind", "theta_sig_ind", "theta_var_ind",
             "theta_P_ind", "theta_0", "init_used")){
    expect_identical(m2[[f]], m[[f]], label = f)
  }

  # Invariance: the permuted theta reproduces the SAME likelihood, recomputed from
  # scratch, not merely carried over.
  mdl_for_ll <- list(y = m2$y, q = 2L, msmu = TRUE, msvar = TRUE, exog = FALSE)
  ll <- MSTest:::logLike_HMmdl(matrix(as.numeric(m2$theta), ncol = 1), mdl_for_ll, 2L)
  expect_lt(abs(ll - m$logLike), 2.3e-13 + 1e-9)   # spec's own verified tolerance, plus slack
})

test_that("standard errors and the Hessian (M3-1) permute correctly, including idempotence", {
  skip_on_cran()
  m <- .canon_swap_fit(21)
  m$control$getSE <- TRUE
  m <- thetaSE(m)
  expect_true(diag(m$P)[1] > diag(m$P)[2])
  m2 <- canonicalize_regimes(m)
  perm <- order(diag(m$P))
  theta_idx <- MSTest:::.canon_theta_perm_idx(2L, 2L, 0L, 0L, TRUE, TRUE, perm)
  expect_equal(as.numeric(m2$theta_se), as.numeric(m$theta_se)[theta_idx])
  expect_equal(m2$info_mat, m$info_mat[theta_idx, theta_idx])
  expect_equal(dim(m2$Hess), dim(m$Hess))   # still in the SE-reduced space, same size

  # Idempotence: canonicalizing an already-canonical object changes nothing, bitwise.
  m3 <- canonicalize_regimes(m2)
  expect_identical(m3$theta, m2$theta)
  expect_identical(m3$P, m2$P)
  expect_identical(m3$Hess, m2$Hess)
  expect_identical(m3$theta_se, m2$theta_se)
})

test_that("a near-absorbing raw fit starts every burnin=0 draw absorbed; the canonicalized one does not", {
  # The actual failure mode this primitive exists to prevent, on a real EM fit --
  # not a constructed example.
  skip_on_cran()
  m <- .canon_swap_fit(21, near_absorbing = TRUE)
  expect_true(diag(m$P)[1] > diag(m$P)[2])   # raw state 1 is the near-absorbing one
  m2 <- canonicalize_regimes(m)
  expect_true(diag(m2$P)[1] < diag(m2$P)[2])

  never_left_raw <- 0L; never_left_canon <- 0L
  for (i in 1:25) {
    set.seed(4000 + i)
    d_raw <- simuHMM(list(n = 100, q = 2, k = 2, mu = m$mu, sigma = m$sigma, P = m$P), burnin = 0)
    set.seed(4000 + i)
    d_can <- simuHMM(list(n = 100, q = 2, k = 2, mu = m2$mu, sigma = m2$sigma, P = m2$P), burnin = 0)
    if (all(d_raw$St == 0)) never_left_raw <- never_left_raw + 1L
    if (all(d_can$St == 0)) never_left_canon <- never_left_canon + 1L
  }
  expect_gt(never_left_raw, 10L)     # the raw fit gets stuck routinely
  expect_equal(never_left_canon, 0L) # the canonicalized fit essentially never does
})

test_that("MSARmdl (univariate) permutes correctly, and phi/betaZ are left untouched", {
  skip_on_cran()
  set.seed(33)
  s <- simuMSAR(list(n = 300, k = 2, mu = c(3, -1), sigma = c(1, 2), phi = 0.3,
                     P = cbind(c(0.2, 0.8), c(0.85, 0.15))), burnin = 100)
  y <- matrix(s$y, ncol = 1)
  m <- MSARmdl(y, p = 1, k = 2, control = list(msmu = TRUE, msvar = TRUE, getSE = FALSE, use_diff_init = 5))
  skip_if(diag(m$P)[1] <= diag(m$P)[2], "fit happened to already be canonical at this seed")
  m2 <- canonicalize_regimes(m)
  perm <- order(diag(m$P))
  expect_equal(diag(m2$P), sort(diag(m$P)))
  expect_equal(m2$mu, m$mu[perm, , drop = FALSE])
  expect_equal(m2$sigma, m$sigma[perm, , drop = FALSE])
  expect_equal(m2$beta, m$beta[, perm, drop = FALSE])   # regime COLUMNS in the univariate case
  expect_identical(m2$phi, m$phi)                       # no regime axis: unchanged
  ll <- MSTest:::logLike_MSARmdl(matrix(as.numeric(m2$theta), ncol = 1),
                                 list(y = m2$y, x = m2$x, q = 1L, p = 1L, msmu = TRUE, msvar = TRUE, exog = FALSE), 2L)
  expect_lt(abs(ll - m$logLike), 1e-8)
})
