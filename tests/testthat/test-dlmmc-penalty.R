# Regression test for the DL-MMC stationarity-penalty sign.
#
# Background: DLMMCpval_fun applies a penalty (lambda) to non-stationary
# candidate parameters so the optimizer avoids them. The penalty must be
# signed so that it REPELS the optimizer (it is a constraint, not an optimum),
# matching the convention used in MMCLRpval_fun. A previous sign error made the
# penalty an ATTRACTOR: when the nuisance search box extended past the
# stationarity boundary (phi + eps > 1, which happens when a linear AR model is
# fit to mean-switching data and the OLS estimate is close to 1), the optimizer
# selected the non-stationary candidate and returned the raw penalty constant
# (e.g. 100) as the "p-value".
#
# This test constructs exactly that situation and asserts the returned p-values
# are valid probabilities in [0, 1] for every available optimizer.

make_trigger_data <- function(seed = 12345) {
  # Mean-switching MS-AR data; the linear AR(1) null fit exhibits spurious
  # persistence (phi_hat ~ 0.75), so phi_hat + eps pushes the search box past 1.
  set.seed(seed)
  simuMSAR(list(n = 500, mu = c(5, 10), sigma = c(1, 2), phi = c(0.5),
                k = 2, P = rbind(c(0.9, 0.1), c(0.1, 0.9))))$y
}

test_that("DLMMCTest never returns a p-value outside [0, 1] (penalty repels)", {
  y <- make_trigger_data()
  for (ot in c("GenSA", "pso")) {
    set.seed(12345)
    out <- DLMMCTest(y, p = 1,
                     control = list(N = 19, eps = 0.3, CI_union = FALSE,
                                    optim_type = ot, threshold_stop = 1, maxit = 50,
                                    silence = TRUE))
    expect_true(is.finite(out$pval_min),  info = ot)
    expect_true(is.finite(out$pval_prod), info = ot)
    expect_gte(as.numeric(out$pval_min),  0)
    expect_lte(as.numeric(out$pval_min),  1)
    expect_gte(as.numeric(out$pval_prod), 0)
    expect_lte(as.numeric(out$pval_prod), 1)
    # the maximizing nuisance value must itself be stationary (|phi| < 1)
    tmax <- as.numeric(out$theta_max_min)
    expect_true(all(abs(tmax) < 1), info = ot)
  }
})

test_that("DLMMCpval_fun penalty is negative for non-stationary theta", {
  # Directly exercise the internal objective: a non-stationary phi (>= 1) must
  # score WORSE (more negative) than any stationary candidate, so that a
  # maximizer of the p-value never selects it.
  skip_if_not(exists("DLMMCpval_fun"))
  y <- as.matrix(make_trigger_data())
  x <- as.matrix(y[-length(y)])          # single lag
  yt <- as.matrix(y[-1])
  params <- matrix(0, nrow = 2, ncol = 4)
  sim_stats <- rep(0, 20)
  pen <- DLMMCpval_fun(theta = c(1.2), y = yt, x = x, params = params,
                       sim_stats = sim_stats, pval_type = "min",
                       stationary_ind = TRUE, lambda = 100)
  expect_lt(pen, 0)   # non-stationary -> negative penalty (repeller)
})
