# Regression tests for the hidden-Markov EM mean update: the M-step for mu must be
# the exact conditional maximizer of Q given the previous betaZ and sigma. Covers the
# common-mean (precision-weighted) branch, the switching-mean branch, the exogenous
# adjustment, and the degenerate-regime fallback.

# Q_obs(mu) with the smoothed probabilities and (betaZ, sigma) held at the values the
# M-step conditions on. Written out here so the assertions do not depend on the code
# under test.
.hm_Qmu <- function(y, Zdm, xi, betaZ, Sig, mu_k) {
  Tn <- nrow(y); q <- ncol(y); tot <- 0
  for (kk in seq_len(nrow(mu_k))) {
    R <- y - matrix(mu_k[kk, ], Tn, q, byrow = TRUE) - Zdm %*% betaZ
    tot <- tot + sum(xi[, kk] * (-0.5 * rowSums((R %*% solve(Sig[[kk]])) * R)))
  }
  tot
}
.hm_vech <- function(S) S[lower.tri(S, diag = TRUE)]

test_that("the common-mean M-step is the sample mean when the variance does not switch", {
  # With a single common covariance the precisions cancel and the exact maximizer of
  # Q over mu is colMeans(y); the previous iterate must not enter the update.
  set.seed(4101)
  for (q in c(1, 3)) for (k in c(2, 3)) {
    Tn <- 150
    y <- matrix(rnorm(Tn * q, 2, 1.5), Tn, q)
    P <- matrix(0.2 / (k - 1), k, k); diag(P) <- 0.8
    mdl <- list(y = y, q = as.integer(q), msmu = FALSE, msvar = FALSE, exog = FALSE)
    # start the mean far from ybar so an accumulating update cannot coincide with it
    th <- c(rep(10, q), .hm_vech(diag(q)), as.numeric(P))
    E <- MSTest:::ExpectationM_HMmdl(matrix(th, ncol = 1), mdl, as.integer(k))
    M <- MSTest:::EMaximization_HMmdl(matrix(th, ncol = 1), mdl, E, as.integer(k))
    expect_equal(as.numeric(M$mu), colMeans(y), tolerance = 1e-12)
  }
})

test_that("the common-mean M-step maximizes Q under regime-specific covariances", {
  # msvar = TRUE is the case the sample-mean identity cannot discriminate: the exact
  # maximizer is precision weighted, mu = [sum_k w_k S_k^-1]^-1 sum_k S_k^-1 m_k.
  set.seed(4102)
  q <- 2; k <- 2; Tn <- 200
  y <- matrix(rnorm(Tn * q, 2, 1.5), Tn, q)
  Sg <- list(matrix(c(1, .3, .3, 1), 2), matrix(c(4, -1, -1, 2), 2))
  P <- cbind(c(.9, .1), c(.2, .8))
  mdl <- list(y = y, q = 2L, msmu = FALSE, msvar = TRUE, exog = FALSE)
  th <- c(c(1.2, 0.4), unlist(lapply(Sg, .hm_vech)), as.numeric(P))
  E <- MSTest:::ExpectationM_HMmdl(matrix(th, ncol = 1), mdl, 2L)
  M <- MSTest:::EMaximization_HMmdl(matrix(th, ncol = 1), mdl, E, 2L)
  xi <- E$xi_t_T
  w  <- colSums(xi)
  A  <- w[1] * solve(Sg[[1]]) + w[2] * solve(Sg[[2]])
  b  <- solve(Sg[[1]]) %*% (t(y) %*% xi[, 1]) + solve(Sg[[2]]) %*% (t(y) %*% xi[, 2])
  expect_equal(as.numeric(M$mu), as.numeric(solve(A, b)), tolerance = 1e-10)
  # and it is not the unweighted sample mean, so the test has content
  expect_gt(max(abs(as.numeric(M$mu) - colMeans(y))), 1e-4)
})

test_that("the mean update removes the exogenous contribution before maximizing", {
  # mu is a conditional maximizer given the PREVIOUS betaZ, so it must be computed on
  # y - Zdm*betaZ. Checked against a direct maximization of Q for both msmu settings.
  skip_on_cran()
  set.seed(4103)
  q <- 2; k <- 2; qz <- 2; Tn <- 200
  y <- matrix(rnorm(Tn * q, 2, 1.5), Tn, q)
  Z <- matrix(rnorm(Tn * qz), Tn, qz)
  Zdm <- scale(Z, TRUE, FALSE)
  bz <- matrix(c(0.8, -0.5, 0.3, 0.6), qz, q)
  Sg <- list(matrix(c(1, .3, .3, 1), 2), matrix(c(3, -.5, -.5, 2), 2))
  P <- cbind(c(.9, .1), c(.2, .8))
  for (msmu in c(TRUE, FALSE)) {
    mu0 <- if (msmu) matrix(c(1, 0, 3, -1), 2, 2, byrow = TRUE) else matrix(c(1.5, 0), 1, 2)
    mdl <- list(y = y, q = 2L, msmu = msmu, msvar = TRUE, exog = TRUE, Z = Z)
    th <- c(as.numeric(t(mu0)), as.numeric(bz), unlist(lapply(Sg, .hm_vech)), as.numeric(P))
    E <- MSTest:::ExpectationM_HMmdl(matrix(th, ncol = 1), mdl, 2L)
    M <- MSTest:::EMaximization_HMmdl(matrix(th, ncol = 1), mdl, E, 2L)
    mu_k <- if (msmu) M$mu else matrix(rep(as.numeric(M$mu), each = k), k, q)
    f <- function(v) -.hm_Qmu(y, Zdm, E$xi_t_T, bz, Sg,
                              if (msmu) matrix(v, k, q, byrow = TRUE)
                              else matrix(rep(v, each = k), k, q))
    opt <- optim(as.numeric(t(mu0)), f, method = "BFGS",
                 control = list(reltol = 1e-14, maxit = 5000))
    expect_lt(f(as.numeric(t(M$mu))) - opt$value, 1e-8)
  }
})

test_that("the mean update keeps the previous iterate when every regime is degenerate", {
  # Section 18 "Bug 3" fallback: an unusable weighting must leave mu at theta's value.
  set.seed(4104)
  Tn <- 120; q <- 2
  y <- matrix(rnorm(Tn * q, 3, 1), Tn, q)
  P <- cbind(c(.9, .1), c(.2, .8))
  mdl <- list(y = y, q = 2L, msmu = FALSE, msvar = TRUE, exog = FALSE)
  th <- c(c(2.5, -0.5), rep(0, 3), rep(0, 3), as.numeric(P))   # both covariances zero
  M <- MSTest:::EMaximization_HMmdl(matrix(th, ncol = 1), mdl,
         MSTest:::ExpectationM_HMmdl(matrix(th, ncol = 1), mdl, 2L), 2L)
  expect_equal(as.numeric(M$mu), c(2.5, -0.5))
})

test_that("the M-step never lowers the part of Q it controls", {
  # EM splits Q into an observation part, a transition part, and a term in the initial
  # distribution that the M-step does not optimize (xi_0 is pinned to the ergodic
  # distribution of the new P). Ascent of the observed log-likelihood is therefore
  # confounded by that third term; ascent of Q_obs + Q_trans is not.
  #
  # Two design points, both established by measurement rather than assumed. The series
  # must carry real regime structure, since with iid data the smoothed weights are near
  # uniform and a wrong mean update satisfies the invariant anyway. And the covariance
  # must NOT switch: with a covariance per regime, the covariance update absorbs a wrong
  # mean and total Q still rises, so the invariant holds on a broken build too.
  skip_on_cran()
  .hm_Qparts <- function(E, y, Zdm, k, q, mu_k, betaZ, Sig, Pnew) {
    Tn <- nrow(y)
    q_obs <- 0
    for (kk in seq_len(k)) {
      R  <- y - matrix(mu_k[kk, ], Tn, q, byrow = TRUE) - Zdm %*% betaZ
      Si <- solve(Sig[[kk]])
      dens <- -0.5 * (q * log(2 * pi) + log(det(Sig[[kk]])) + rowSums((R %*% Si) * R))
      q_obs <- q_obs + sum(E$xi_t_T[, kk] * dens)
    }
    # smoothed joint Pr(s_t = j, s_{t+1} = i | Y) from Kim's identity, exact at p = 0
    q_tr <- 0
    for (t in seq_len(Tn - 1)) for (i in seq_len(k)) {
      den <- max(E$xi_tp1_t[t, i], 1e-300)
      for (j in seq_len(k)) {
        jt <- E$xi_t_T[t + 1, i] * E$P[i, j] * E$xi_t_t[t, j] / den
        if (jt > 0) q_tr <- q_tr + jt * log(max(Pnew[i, j], 1e-300))
      }
    }
    q_obs + q_tr
  }
  set.seed(4106)
  S <- simuHMM(list(n = 180, q = 2, k = 2, mu = rbind(c(0, 0), c(3, -2)),
                    sigma = list(diag(2), matrix(c(2, .5, .5, 1.5), 2)),
                    P = cbind(c(.9, .1), c(.15, .85))))
  y <- S$y; Tn <- 180L; q <- 2L; k <- 2L
  Zfull <- matrix(rnorm(Tn * 2), Tn, 2)
  for (msmu in c(TRUE, FALSE)) for (exog in c(FALSE, TRUE)) {
    Z   <- if (exog) Zfull else NULL
    Zdm <- if (exog) scale(Z, TRUE, FALSE) else matrix(0, Tn, 1)
    mdl <- list(y = y, q = q, msmu = msmu, msvar = FALSE, exog = exog)
    if (exog) mdl$Z <- Z
    # start at the fitted optimum, where a correct M-step is stationary
    fit <- HMmdl(y, k, Z, control = list(msmu = msmu, msvar = FALSE, getSE = FALSE,
                                         use_diff_init = 3, maxit = 400))
    th <- as.numeric(fit$theta); nmu <- if (msmu) 4L else 2L
    worst <- 0
    for (it in 1:8) {
      E <- MSTest:::ExpectationM_HMmdl(matrix(th, ncol = 1), mdl, k)
      M <- MSTest:::EMaximization_HMmdl(matrix(th, ncol = 1), mdl, E, k)
      bz_old <- if (exog) matrix(th[(nmu + 1):(nmu + 4)], 2, 2) else matrix(0, 1, 2)
      mk_old <- if (msmu) matrix(th[1:4], 2, 2, byrow = TRUE) else matrix(rep(th[1:2], each = 2), 2, 2)
      mk_new <- if (msmu) M$mu else matrix(rep(as.numeric(M$mu), each = 2), 2, 2)
      bz_new <- if (exog) M$betaZ else matrix(0, 1, 2)
      worst <- min(worst,
                   .hm_Qparts(E, y, Zdm, k, q, mk_new, bz_new, M$sigma, M$P) -
                   .hm_Qparts(E, y, Zdm, k, q, mk_old, bz_old, E$sigma, E$P))
      th <- as.numeric(M$theta)
    }
    # round-off only; a broken mean update loses 1e-6 to 1e-3 here
    expect_gt(worst, -1e-9)
  }
})

test_that("a non-finite exogenous coefficient does not become an absorbing state", {
  # The mean update conditions on the previous betaZ, so a non-finite value there
  # reaches mu through the exogenous-adjusted series. Without a guard mu goes
  # non-finite, the betaZ solve then fails and keeps the non-finite value, and the
  # iteration is stuck. The mean must fall back to its previous value instead.
  skip_on_cran()
  set.seed(303)
  Tn <- 160
  Zx <- matrix(rnorm(Tn), Tn, 1)
  Y  <- matrix(rnorm(Tn), Tn, 1) + 0.7 * Zx
  for (bad in c(NaN, Inf, -Inf)) {
    th <- c(-1, 1, bad, 1, 2, 0.9, 0.1, 0.2, 0.8)   # mu1, mu2, betaZ, sig1, sig2, vec(P)
    m <- HMmdl(Y, k = 2, Z = Zx,
               control = list(msmu = TRUE, msvar = TRUE, getSE = FALSE, init_theta = th))
    expect_true(all(is.finite(m$theta)))
    expect_true(is.finite(m$logLike))
  }
})

test_that("a per-regime covariance list shorter than k is rejected", {
  # The covariance list is indexed per regime in two places; Rcpp::List does not
  # bounds check, so a short list reads past the end and aborts the session.
  skip_on_cran()
  set.seed(77)
  Tn <- 120; q <- 2; k <- 2
  y   <- matrix(rnorm(Tn * q), Tn, q)
  mdl <- list(y = y, q = q, msmu = FALSE, msvar = TRUE, exog = FALSE)
  th  <- c(0, 0, 1, 0, 1, 1, 0, 1, 0.9, 0.1, 0.1, 0.9)
  E   <- MSTest::ExpectationM_HMmdl(th, mdl, k)
  expect_length(E$sigma, k)                       # every public path supplies k
  for (short in list(E$sigma[1], list())) {
    E_bad <- E; E_bad$sigma <- short
    expect_error(MSTest::EMaximization_HMmdl(th, mdl, E_bad, k),
                 "one covariance matrix per regime")
  }
})
