# Tests for the EM transition-probability bound (em_transition_constraint):
# the exact water-filling M-step, the off switch, binding behavior, start
# projection, and input validation.

set.seed(7100)
pb_sim <- simuMSAR(list(n = 300, k = 2, mu = c(-1, 1), sigma = c(1, 1), phi = 0.3,
                        P = cbind(c(0.90, 0.10), c(0.10, 0.90))), burnin = 100)
pb_y <- matrix(pb_sim$y, ncol = 1)
set.seed(7200)
pb_sim_persist <- simuMSAR(list(n = 400, k = 2, mu = c(0, 2), sigma = c(1, 1), phi = 0.3,
                                P = cbind(c(0.99, 0.01), c(0.01, 0.99))), burnin = 100)
pb_y_persist <- matrix(pb_sim_persist$y, ncol = 1)

test_that("constrain_P_col: exactness, pass-through, and defensive fallbacks", {
  cpc <- MSTest:::constrain_P_col
  f3 <- rep(1/3, 3)
  expect_equal(as.numeric(cpc(c(0, 0, 10), 0.05, f3)), c(0.05, 0.05, 0.90))
  expect_identical(as.numeric(cpc(c(3, 4, 5), 0.05, f3)), c(3, 4, 5) / sum(c(3, 4, 5)))
  # scale invariance
  x <- c(40, 30, 20, 0.2); f4 <- rep(0.25, 4)
  expect_identical(as.numeric(cpc(x, 0.1, f4)), as.numeric(cpc(1e6 * x, 0.1, f4)))
  # fallback on unusable inputs
  expect_identical(as.numeric(cpc(c(-1, 2), 0.05, c(0.5, 0.5))), c(0.5, 0.5))
  expect_identical(as.numeric(cpc(c(NaN, 1), 0.05, c(0.5, 0.5))), c(0.5, 0.5))
  expect_identical(as.numeric(cpc(c(0, 0), 0.05, c(0.5, 0.5))), c(0.5, 0.5))
  expect_identical(as.numeric(cpc(c(9, 1), 0.5, c(0.5, 0.5))), c(0.5, 0.5))   # eps = 1/k
  expect_identical(as.numeric(cpc(c(1e308, 1e308), 0.02, c(0.6, 0.4))), c(0.6, 0.4))  # sum overflow
  # feasible and exact over a randomized grid
  set.seed(7300)
  for (i in 1:25) {
    k <- sample(c(2, 3, 6), 1); e <- runif(1, 0.005, 0.9 / k)
    n <- 10^runif(k, -3, 3)
    p <- as.numeric(cpc(n, e, rep(1 / k, k)))
    expect_true(all(p >= e - 1e-12))
    expect_lt(abs(sum(p) - 1), 1e-12)
  }
})

test_that("off switch: default control is bit-identical to em_transition_constraint = 0", {
  skip_on_cran()
  ctl0 <- list(getSE = FALSE, msmu = TRUE, msvar = FALSE, use_diff_init = 2)
  set.seed(7401)
  f1 <- MSARmdl(pb_y, p = 1, k = 2, control = ctl0)
  set.seed(7401)
  f2 <- MSARmdl(pb_y, p = 1, k = 2, control = c(ctl0, list(em_transition_constraint = 0)))
  expect_identical(f1$theta, f2$theta)
  expect_identical(f1$logLike, f2$logLike)
  expect_null(f1$P_constraint_binding)
  expect_null(f2$P_constraint_binding)
})

test_that("binding fit pins at eps with exact flags, NA SEs, and valid P", {
  skip_on_cran()
  # informed start at the projected truth: tests the pinning mechanics
  # deterministically (random multistart can settle on interior local modes)
  set.seed(7500)
  base <- MSARmdl(pb_y_persist, p = 1, k = 2,
                  control = list(getSE = FALSE, msmu = TRUE, msvar = FALSE, use_diff_init = 2))
  th <- base$theta; np <- length(th)
  th[1:2] <- c(0, 2); th[3] <- 0.3; th[4] <- 1
  th[(np - 3):np] <- c(0.95, 0.05, 0.05, 0.95)
  set.seed(7501)
  fit <- suppressWarnings(MSARmdl(pb_y_persist, p = 1, k = 2,
           control = list(getSE = TRUE, msmu = TRUE, msvar = FALSE, use_diff_init = 1,
                          init_theta = th, em_transition_constraint = 0.05)))
  P <- as.matrix(fit$P)
  expect_equal(min(P), 0.05, tolerance = 1e-10)
  expect_true(all(abs(colSums(P) - 1) < 1e-12))
  expect_identical(as.vector(fit$P_constraint_binding), as.vector(P <= 0.05 + 5e-11))
  expect_equal(fit$em_transition_constraint, 0.05)
  ind <- which(as.numeric(fit$theta_P_ind) == 1)
  expect_true(all(is.na(fit$theta_se[ind[as.vector(fit$P_constraint_binding)]])))
  expect_true(all(is.finite(fit$theta_se[ind[!as.vector(fit$P_constraint_binding)]])))
})

test_that("start projection: infeasible or malformed init P becomes feasible", {
  skip_on_cran()
  set.seed(7601)
  base <- MSARmdl(pb_y, p = 1, k = 2,
                  control = list(getSE = FALSE, msmu = TRUE, msvar = FALSE, use_diff_init = 2))
  th <- base$theta; np <- length(th)
  th[(np - 3):np] <- c(0.999, 0.001, 0.001, 0.999)   # below the bound
  set.seed(7602)
  f <- suppressWarnings(MSARmdl(pb_y, p = 1, k = 2,
         control = list(getSE = FALSE, msmu = TRUE, msvar = FALSE, use_diff_init = 1,
                        init_theta = th, em_transition_constraint = 0.05)))
  expect_true(min(f$P) >= 0.05 - 1e-12)
})

test_that("init_theta P-block validation errors on malformed input, all paths", {
  set.seed(7605)
  base <- MSARmdl(pb_y, p = 1, k = 2,
                  control = list(getSE = FALSE, msmu = TRUE, msvar = FALSE, use_diff_init = 2))
  th <- base$theta; np <- length(th)
  bad_sum <- th; bad_sum[(np - 3):np] <- rep(0.9, 4)        # column sums 1.8
  bad_rng <- th; bad_rng[(np - 3):np] <- c(1.2, -0.2, 0.1, 0.9)  # outside [0,1]
  bad_nan <- th; bad_nan[np] <- NaN
  for (ctl_extra in list(list(), list(em_transition_constraint = 0.05))) {
    ctl <- c(list(getSE = FALSE, msmu = TRUE, msvar = FALSE, use_diff_init = 1), ctl_extra)
    expect_error(MSARmdl(pb_y, p = 1, k = 2, control = c(ctl, list(init_theta = bad_sum))),
                 "sum to 1")
    expect_error(MSARmdl(pb_y, p = 1, k = 2, control = c(ctl, list(init_theta = bad_rng))),
                 "in \\[0,1\\]")
    expect_error(MSARmdl(pb_y, p = 1, k = 2, control = c(ctl, list(init_theta = bad_nan))),
                 "finite")
  }
  # a valid stored start (columns sum to 1 up to ulp) passes
  expect_silent(suppressWarnings(MSARmdl(pb_y, p = 1, k = 2,
    control = list(getSE = FALSE, msmu = TRUE, msvar = FALSE, use_diff_init = 1,
                   init_theta = th))))
})

test_that("validation: eps >= 1/k errors for EM and MLE; MLE fits have NULL flags", {
  expect_error(MSARmdl(pb_y, p = 1, k = 2, control = list(em_transition_constraint = 0.5)),
               "less than 1/k")
  expect_error(MSARmdl(pb_y, p = 1, k = 2,
                       control = list(method = "MLE", em_transition_constraint = 0.5)),
               "less than 1/k")
  skip_on_cran()
  set.seed(7701)
  m <- suppressWarnings(MSARmdl(pb_y, p = 1, k = 2,
         control = list(getSE = FALSE, msmu = TRUE, msvar = FALSE, use_diff_init = 2,
                        method = "MLE", em_transition_constraint = 0.05)))
  expect_null(m$P_constraint_binding)
})

test_that("k = 3: bound respected, diagonal cap honored, ascent under binding", {
  skip_on_cran()
  set.seed(7801)
  s3 <- simuMSAR(list(n = 500, k = 3, mu = c(-3, 0, 3), sigma = c(1, 1, 1), phi = 0.2,
                      P = cbind(c(0.96, 0.02, 0.02), c(0.02, 0.96, 0.02),
                                c(0.02, 0.02, 0.96))), burnin = 100)
  y3 <- matrix(s3$y, ncol = 1)
  # informed start at the (projected) truth
  set.seed(7802)
  b3 <- MSARmdl(y3, p = 1, k = 3,
                control = list(getSE = FALSE, msmu = TRUE, msvar = FALSE, use_diff_init = 2))
  th <- b3$theta; np <- length(th)
  th[1:3] <- c(-3, 0, 3); th[4] <- 0.2; th[5] <- 1
  th[(np - 8):np] <- c(0.96, 0.02, 0.02, 0.02, 0.96, 0.02, 0.02, 0.02, 0.96)
  set.seed(7803)
  f3 <- suppressWarnings(MSARmdl(y3, p = 1, k = 3,
          control = list(getSE = FALSE, msmu = TRUE, msvar = FALSE, use_diff_init = 1,
                         init_theta = th, em_transition_constraint = 0.05)))
  expect_true(min(f3$P) >= 0.05 - 1e-12)
  expect_true(all(diag(as.matrix(f3$P)) <= 1 - 2 * 0.05 + 1e-10))
  expect_true(all(abs(colSums(as.matrix(f3$P)) - 1) < 1e-10))
  # chained single-iteration ascent trace on a binding k = 2 fit
  ctl1 <- list(getSE = FALSE, msmu = TRUE, msvar = FALSE, use_diff_init = 1, maxit = 1,
               em_transition_constraint = 0.05, em_best_iterate = FALSE)
  set.seed(7804)
  m <- suppressWarnings(MSARmdl(pb_y_persist, p = 1, k = 2, control = ctl1))
  lls <- m$logLike; thc <- m$theta
  for (it in 1:40) {
    ctl1$init_theta <- thc
    m <- suppressWarnings(MSARmdl(pb_y_persist, p = 1, k = 2, control = ctl1))
    lls <- c(lls, m$logLike); thc <- m$theta
  }
  expect_gt(min(diff(lls)), -1e-2)
})
