test_that("EM does not descend the likelihood (MS-AR, switching mean and variance)", {
  skip_on_cran()
  set.seed(41)
  simu <- simuMSAR(list(n = 250, k = 2, mu = c(0, 2), sigma = c(1, 4), phi = 0.5,
                        P = cbind(c(0.9, 0.1), c(0.15, 0.85))), burnin = 100)
  Y <- matrix(simu$y, ncol = 1)
  ctl <- list(getSE = FALSE, msmu = TRUE, msvar = TRUE, use_diff_init = 1)
  set.seed(7)
  th <- as.numeric(MSARmdl(Y, p = 1, k = 2, control = c(ctl, list(maxit = 3)))$theta)
  # A tighter stopping tolerance must not attain a lower likelihood than a looser
  # one from the same start (allowing the small ergodic-initialization effect).
  fit <- function(ltol) MSARmdl(Y, p = 1, k = 2,
    control = c(ctl, list(maxit = 1000, conv = "loglik", ltol = ltol, init_theta = th)))
  expect_gte(logLik(fit(1e-7)), logLik(fit(1e-5)) - 1e-2)
  # Chained single iterations: no step may lose more than the documented
  # ergodic-initialization scale at interior parameters.
  ll_of <- function(m) MSTest:::logLike_MSARmdl(matrix(as.numeric(m$theta), ncol = 1), m, 2L)
  cur <- th; ll <- -Inf; worst <- 0
  for (i in 1:15) {
    m <- MSARmdl(Y, p = 1, k = 2, control = c(ctl,
           list(maxit = 1, conv = "theta", thtol = 0, init_theta = cur)))
    ll_new <- ll_of(m)
    if (is.finite(ll) && is.finite(ll_new)) worst <- min(worst, ll_new - ll)
    cur <- as.numeric(m$theta); ll <- ll_new
  }
  expect_gt(worst, -1e-2)
})

test_that("limP returns a valid probability vector at boundary transition matrices", {
  cases <- list(diag(2), diag(3),
                cbind(c(0.9, 0.1), c(0, 1)),                        # absorbing column
                cbind(c(1, 0, 0), c(0, 1, 0), c(0.2, 0.3, 0.5)),    # two absorbing states
                cbind(c(0, 1), c(1, 0)))                            # 2-cycle
  for (P in cases) {
    v <- limP(P)
    expect_true(all(is.finite(v)))
    expect_true(all(v >= 0))
    expect_lt(abs(sum(v) - 1), 1e-12)
  }
})

test_that("smoothed probabilities equal the extended-state aggregation", {
  skip_on_cran()
  set.seed(42)
  simu <- simuMSAR(list(n = 200, k = 2, mu = c(0, 2), sigma = c(1, 4), phi = 0.5,
                        P = cbind(c(0.9, 0.1), c(0.15, 0.85))), burnin = 100)
  m <- MSARmdl(matrix(simu$y, ncol = 1), p = 1, k = 2,
               control = list(getSE = FALSE, msmu = TRUE, msvar = TRUE,
                              use_diff_init = 3, maxit = 300))
  E <- MSTest:::ExpectationM_MSARmdl(matrix(as.numeric(m$theta), ncol = 1), m, 2L)
  agg <- cbind(E$xi_t_T_AR[, 1] + E$xi_t_T_AR[, 3],
               E$xi_t_T_AR[, 2] + E$xi_t_T_AR[, 4])
  expect_lt(max(abs(E$xi_t_T - agg)), 1e-12)
  expect_lt(max(abs(rowSums(m$St) - 1)), 1e-8)
})

test_that("EM does not descend the likelihood for any switching model or variant", {
  # Full coverage of the five Markov-switching classes against the three switching
  # specifications that matter. The mean update, the variance update and the
  # transition update are each exercised alone (one flag off) and together.
  #
  # The tolerance is not arbitrary: the M-step maximizes the observation and
  # transition parts of Q but not the term in the initial distribution, which is
  # reset to the ergodic distribution of the new P. That channel can lose a small
  # amount of likelihood per iteration at interior parameters, which is why the
  # designs below keep the transition matrix away from its boundary.
  skip_on_cran()
  set.seed(8801)
  Tn <- 200L
  P2 <- cbind(c(0.90, 0.10), c(0.15, 0.85))
  uni <- simuMSAR(list(n = Tn, k = 2, mu = c(0, 3), sigma = c(1, 4), phi = 0.4, P = P2),
                  burnin = 100)
  biv <- simuMSVAR(list(n = Tn, p = 1, q = 2, k = 2, mu = rbind(c(0, 0), c(3, -2)),
                        sigma = list(diag(2), 2 * diag(2)),
                        phi = matrix(c(0.3, 0, 0, 0.3), 2, 2), P = P2), burnin = 100)
  hm  <- simuHMM(list(n = Tn, q = 2, k = 2, mu = rbind(c(0, 0), c(3, -2)),
                      sigma = list(diag(2), 2 * diag(2)), P = P2), burnin = 100)
  y1 <- matrix(uni$y, ncol = 1)
  y2 <- biv$y
  yh <- hm$y
  Z1 <- matrix(rnorm(Tn), Tn, 1)
  Z2 <- matrix(rnorm(Tn * 2), Tn, 2)
  y1x <- y1 + Z1 * 0.8
  y2x <- y2 + Z2 %*% matrix(c(0.6, -0.4, 0.3, 0.5), 2, 2)
  yhx <- yh + Z2 %*% matrix(c(0.7, -0.5, 0.4, 0.6), 2, 2)

  # HMmdl carries the exogenous case itself rather than through a separate class, so
  # it appears twice: without regressors, and with them.
  fits <- list(
    HMmdl     = function(ct) HMmdl(yh, k = 2, control = ct),
    HMmdlZ    = function(ct) HMmdl(yhx, k = 2, Z = Z2, control = ct),
    MSARmdl   = function(ct) MSARmdl(y1, p = 1, k = 2, control = ct),
    MSARXmdl  = function(ct) MSARXmdl(y1x, p = 1, k = 2, Z = Z1, control = ct),
    MSVARmdl  = function(ct) MSVARmdl(y2, p = 1, k = 2, control = ct),
    MSVARXmdl = function(ct) MSVARXmdl(y2x, p = 1, k = 2, Z = Z2, control = ct)
  )
  lls <- list(
    HMmdl     = MSTest:::logLike_HMmdl,     HMmdlZ    = MSTest:::logLike_HMmdl,
    MSARmdl   = MSTest:::logLike_MSARmdl,
    MSARXmdl  = MSTest:::logLike_MSARXmdl,  MSVARmdl  = MSTest:::logLike_MSVARmdl,
    MSVARXmdl = MSTest:::logLike_MSVARXmdl
  )
  variants <- list(c(msmu = TRUE,  msvar = FALSE),
                   c(msmu = FALSE, msvar = TRUE),
                   c(msmu = TRUE,  msvar = TRUE))

  for (cls in names(fits)) for (v in variants) {
    base_ct <- list(getSE = FALSE, msmu = unname(v["msmu"]), msvar = unname(v["msvar"]),
                    use_diff_init = 1, em_best_iterate = FALSE)
    m <- fits[[cls]](c(base_ct, list(maxit = 100)))
    cur <- as.numeric(m$theta); ll <- -Inf; worst <- 0; n_fin <- 0L
    for (i in 1:6) {
      f <- fits[[cls]](c(base_ct, list(maxit = 1, conv = "theta", thtol = 0,
                                       init_theta = cur)))
      ll_new <- lls[[cls]](matrix(as.numeric(f$theta), ncol = 1), f, 2L)
      if (is.finite(ll) && is.finite(ll_new)) {
        worst <- min(worst, ll_new - ll); n_fin <- n_fin + 1L
      }
      cur <- as.numeric(f$theta); ll <- ll_new
    }
    lab <- paste(cls, paste(v, collapse = "/"))
    # without this the assertion below would hold vacuously on an all-NaN build
    expect_gt(n_fin, 4L, label = lab)
    expect_gt(worst, -1e-2, label = lab)
  }
})
