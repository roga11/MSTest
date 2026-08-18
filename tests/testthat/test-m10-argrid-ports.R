# Tests for the C++ ports of argrid_MSARmdl, argrid_MSVARmdl, and arP (the R
# helpers remain the exported reference implementations) and for namespace-only
# estimation (no library() attach).

norm_m <- function(x){ x <- as.matrix(x); dimnames(x) <- NULL; storage.mode(x) <- "double"; x }

test_that("argrid_MSARmdl_cpp matches argrid_MSARmdl including non-switching recycling", {
  set.seed(9101)
  for (k in c(2, 3, 4)) for (ar in c(0, 1, 2, 4)) for (msmu in c(TRUE, FALSE)) for (msvar in c(TRUE, FALSE)) {
    mu  <- if (msmu) rnorm(k, 0, 5) else rnorm(1)
    sig <- if (msvar) 10^runif(k, -6, 6) else 10^runif(1, -6, 6)
    rR <- argrid_MSARmdl(mu, sig, k, ar, msmu, msvar)
    rC <- MSTest:::argrid_MSARmdl_cpp(mu, sig, k, ar)
    expect_identical(norm_m(rR$mu), norm_m(rC$mu))
    expect_identical(norm_m(rR$sig), norm_m(rC$sig))
    expect_identical(as.numeric(rR$state_ind), as.numeric(rC$state_ind))
  }
})

test_that("argrid_MSVARmdl_cpp matches argrid_MSVARmdl for ar >= 1", {
  # ar = 0 excluded: the R version returns a 1 x q dropped matrix there and no
  # live C++ path calls it with ar = 0
  set.seed(9102)
  for (k in c(2, 3)) for (ar in c(1, 2, 3)) for (q in c(1, 2, 3)) {
    mu <- matrix(rnorm(k * q, 0, 3), k, q)
    sigma <- lapply(1:k, function(j){ A <- matrix(rnorm(q * q), q, q); crossprod(A) + diag(q) })
    rR <- argrid_MSVARmdl(mu, sigma, k, ar, TRUE, TRUE)
    rC <- MSTest:::argrid_MSVARmdl_cpp(mu, sigma, k, ar)
    M <- k^(ar + 1)
    for (m in 1:M) {
      expect_identical(norm_m(rR$mu[[m]]), norm_m(rC$mu[[m]]))
      expect_identical(norm_m(rR$sig[[m]]), norm_m(rC$sig[[m]]))
    }
    expect_identical(as.numeric(rR$state_ind), as.numeric(rC$state_ind))
  }
})

test_that("arP_cpp matches arP; ar = 0 returns P; columns sum to 1", {
  set.seed(9103)
  for (k in c(2, 3, 4)) for (ar in 0:3) {
    P <- matrix(abs(rnorm(k * k)), k, k); P <- sweep(P, 2, colSums(P), "/")
    aR <- arP(P, k, ar)
    aC <- MSTest:::arP_cpp(P, k, ar)
    expect_identical(norm_m(aR), norm_m(aC))
    expect_true(max(abs(colSums(norm_m(aC)) - 1)) < 1e-12)
  }
})

test_that("malformed inputs error instead of crashing", {
  expect_error(MSTest:::argrid_MSVARmdl_cpp(matrix(rnorm(4), 2, 2), list(diag(2)), 2L, 1L),
               "at least k")
  expect_error(MSTest:::argrid_MSARmdl_cpp(c(0, 1), c(1, 2), 2L, 31L), "too large")
  expect_error(MSTest:::arP_cpp(diag(2), 2L, 40L), "too large")
})

test_that("estimation works namespace-only (no attach)", {
  skip_on_cran()
  skip_if_not_installed("callr")
  set.seed(9104)
  y <- as.numeric(arima.sim(list(ar = 0.3), n = 120)) + rep(c(0, 2), each = 60)
  fit_child <- callr::r(function(y){
    stopifnot(!"package:MSTest" %in% search())
    set.seed(9105)
    m <- MSTest::MSARmdl(matrix(y, ncol = 1), p = 1, k = 2,
                         control = list(getSE = FALSE, use_diff_init = 2, maxit = 200))
    stopifnot(!"package:MSTest" %in% search())
    list(theta = m$theta, logLike = m$logLike)
  }, args = list(y = y), libpath = .libPaths())
  set.seed(9105)
  m2 <- MSARmdl(matrix(y, ncol = 1), p = 1, k = 2,
                control = list(getSE = FALSE, use_diff_init = 2, maxit = 200))
  expect_identical(fit_child$theta, m2$theta)
  expect_identical(fit_child$logLike, m2$logLike)
})
