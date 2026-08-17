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
