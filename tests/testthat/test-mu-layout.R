# Regime means are reported as a (k x q) matrix with regimes in rows, and the
# reported standard deviation is a standard deviation. Neither was covered before:
# the packaged tests exercised only configurations where the layout defect cannot
# appear (one series, or a switching mean).

test_that("the regime-mean helper inverts the packing the estimator uses", {
  # Fit-free: build theta with unmistakable values and compare against the C++
  # decoder itself. k == q is the case the helper exists for, because a shape test
  # cannot tell rows from columns there.
  set.seed(4210)
  for (q in 2:4) for (k in 2:4) for (msmu in c(TRUE, FALSE)) {
    mu_len <- q * (1 + msmu * (k - 1))
    mu_blk <- seq_len(mu_len) * 100
    Tn  <- 60
    y   <- matrix(rnorm(Tn * q), Tn, q)
    P   <- matrix(0.2 / (k - 1), k, k); diag(P) <- 0.8
    sig <- diag(q)[lower.tri(diag(q), diag = TRUE)]
    th  <- c(mu_blk, rep(sig, k), as.numeric(P))
    mdl <- list(y = y, q = as.integer(q), msmu = msmu, msvar = FALSE, exog = FALSE)
    E   <- MSTest:::ExpectationM_HMmdl(matrix(th, ncol = 1), mdl, as.integer(k))
    expect_equal(MSTest:::.theta_mu_kq(th, k, q, msmu), E$mu, tolerance = 0)
  }
})

test_that("the regime-mean helper reads the flag the way the estimator does", {
  # The estimator converts this flag to a C++ bool, so a numeric 1 is a switching
  # mean there; the helper must not read it as a common mean.
  th <- c(1, 2, 3, 4, rep(0, 3), 0.9, 0.1, 0.1, 0.9)
  expect_equal(MSTest:::.theta_mu_kq(th, 2, 2, 1),
               MSTest:::.theta_mu_kq(th, 2, 2, TRUE))
  expect_equal(MSTest:::.theta_mu_kq(th, 2, 2, 0),
               MSTest:::.theta_mu_kq(th, 2, 2, FALSE))
  # the two branches must differ, or the comparisons above hold vacuously
  expect_false(identical(MSTest:::.theta_mu_kq(th, 2, 2, TRUE),
                         MSTest:::.theta_mu_kq(th, 2, 2, FALSE)))
})

test_that("a non-switching mean is reported identically across regimes", {
  skip_on_cran()
  set.seed(5150)
  Tn <- 250
  # per-series distinct levels: a common mean that is the same in every series is a
  # fixed point of the layout defect and cannot detect it
  S <- simuHMM(list(n = Tn, q = 2, k = 2, mu = rbind(c(4, 11), c(4, 11)),
                    sigma = list(diag(2), 2 * diag(2)),
                    P = cbind(c(.9, .1), c(.15, .85))))
  ctl <- list(msmu = FALSE, msvar = TRUE, getSE = FALSE, use_diff_init = 2)
  fits <- list(HMmdl(S$y, k = 2, control = ctl),
               MSVARmdl(S$y, p = 1, k = 2, control = ctl))
  for (m in fits) expect_equal(apply(m$mu, 2, function(z) diff(range(z))), rep(0, 2))
})

test_that("HMmdl reports an intercept equal to the mean when there are no regressors", {
  skip_on_cran()
  set.seed(5151)
  Tn <- 250
  S <- simuHMM(list(n = Tn, q = 2, k = 2, mu = rbind(c(0, 5), c(4, -3)),
                    sigma = list(diag(2), 2 * diag(2)),
                    P = cbind(c(.9, .1), c(.15, .85))))
  m <- HMmdl(S$y, k = 2, control = list(msmu = TRUE, msvar = TRUE, getSE = FALSE,
                                        use_diff_init = 2))
  expect_equal(m$intercept, m$mu)   # man/HMmdl.Rd: with no Z the two coincide
})

test_that("HMmdl reports a standard deviation, not a variance", {
  skip_on_cran()
  set.seed(5152)
  y <- matrix(c(rnorm(120, 0, 1), rnorm(120, 3, 2)), ncol = 1)
  m <- HMmdl(y, k = 2, control = list(msvar = TRUE, getSE = FALSE, use_diff_init = 2))
  expect_equal(as.numeric(m$stdev), sqrt(as.numeric(unlist(m$sigma))))
})
