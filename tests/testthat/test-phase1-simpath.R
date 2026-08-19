# Regression cells for the Phase-1 simulation-path fixes that test-mdledit-P.R does not reach:
# the column-sum gate (2), MSVARXmdl's $Z (3), the univariate hidden-Markov sigma shape (4),
# exogenous regressors in simulated hidden-Markov draws (5), the non-finite-P simulator guard (6),
# the control-name leak (7) and the "no admissible candidate" guard.

P2 <- cbind(c(0.90, 0.10), c(0.15, 0.85))

# ---------------------------------------------------------------- item 6
test_that("every Markov-switching simulator rejects a non-finite transition matrix", {
  mk <- list(
    MSAR   = list(n = 60, k = 2, mu = c(-1, 2), sigma = c(1, 2), phi = 0.4, P = P2),
    MSARX  = list(n = 60, k = 2, mu = c(-1, 2), sigma = c(1, 2), phi = 0.4, P = P2,
                  betaZ = as.matrix(-0.7), Z = matrix(seq_len(60) / 60, ncol = 1)),
    MSVAR  = list(n = 60, p = 1, q = 2, k = 2, mu = rbind(c(-1, -1), c(2, 2)),
                  sigma = list(diag(2), 2 * diag(2)),
                  phi = rbind(c(0.4, 0.1), c(0.05, 0.3)), P = P2),
    MSVARX = list(n = 60, p = 1, q = 2, k = 2, mu = rbind(c(-1, -1), c(2, 2)),
                  sigma = list(diag(2), 2 * diag(2)),
                  phi = rbind(c(0.4, 0.1), c(0.05, 0.3)), P = P2,
                  betaZ = matrix(c(-0.8, 0.5), 1, 2), Z = matrix(seq_len(60) / 60, ncol = 1)),
    HMM    = list(n = 60, k = 2, q = 1, mu = matrix(c(0, 2.5), 2, 1),
                  sigma = list(matrix(1, 1, 1), matrix(1, 1, 1)), P = P2))
  fn <- list(MSAR = simuMSAR_cpp, MSARX = simuMSARX_cpp, MSVAR = simuMSVAR_cpp,
             MSVARX = simuMSVARX_cpp, HMM = function(m, b) simuHMM_cpp(m, b, FALSE))
  bad <- list(one_NaN_column = { M <- P2; M[, 1] <- NaN; M },
              all_NaN        = matrix(NaN, 2, 2),
              positive_Inf   = { M <- P2; M[1, 1] <- Inf; M },
              mixed_Inf      = { M <- P2; M[1, 1] <- -Inf; M[2, 1] <- Inf; M })
  for (nm in names(fn)) {
    # a NaN column sums to NaN, and NaN > tol is false, so the column-sum test alone
    # lets it through: the finiteness test is what stops it
    for (b in names(bad)) {
      m <- mk[[nm]]; m$P <- bad[[b]]
      expect_error(fn[[nm]](m, 20), "finite", info = paste(nm, b))
    }
    # and a valid P is still accepted, including one a whisker inside the 1e-8 tolerance
    expect_silent(fn[[nm]](mk[[nm]], 20))
    ok <- mk[[nm]]; ok$P[1, 1] <- ok$P[1, 1] + 1e-9
    expect_silent(fn[[nm]](ok, 20))
  }
})

# ---------------------------------------------------------------- item 2
test_that("the MMC column-sum gate is at least as strict as the simulators", {
  skip_on_cran()
  set.seed(9200)
  s <- simuMSAR(list(n = 250, k = 2, mu = c(-1.5, 1.5), sigma = c(1, 1), phi = 0.3, P = P2),
                burnin = 100)
  set.seed(9201)
  m <- MSARmdl(matrix(s$y, ncol = 1), p = 1, k = 2,
               control = list(getSE = FALSE, use_diff_init = 1))
  ind <- which(as.numeric(m$theta_P_ind) == 1)
  h0c <- list(getSE = FALSE, use_diff_init = 1)
  gate <- function(d, thtol = 1e-6) {
    P <- as.matrix(m$P); P[1, 1] <- P[1, 1] + d
    th <- m$theta; th[ind] <- c(P)
    suppressWarnings(MMCLRpval_fun(th, m, 3L, 8, 5L, 100L, 0L, 100, TRUE, thtol,
                                   NULL, FALSE, h0c, h0c))
  }
  # the simulators hard-stop above 1e-8; every such candidate must be rejected by the
  # gate (value = -lambda) rather than sent to a simulation that can only fail
  for (d in c(1.1e-8, 5e-8, 1e-7, 5e-7, 1e-6, 5e-6, 1e-3)) {
    expect_equal(gate(d), -100, info = sprintf("deviation %.1e", d))
  }
  # a deviation the simulators tolerate must not be rejected by the gate
  expect_silent(simuMSAR_cpp(list(n = 40, k = 2, mu = c(-1, 2), sigma = c(1, 1), phi = 0.3,
                                  P = { P <- P2; P[1, 1] <- P[1, 1] + 1e-9; P }), 20))
  # the gate never becomes looser than it was: it is min(thtol, 1e-8), so a user thtol
  # below 1e-8 still binds
  expect_equal(gate(5e-10, thtol = 1e-12), -100)
  # a fitted P is exactly column-stochastic, so the seed is never caught by the gate
  expect_lt(max(abs(colSums(as.matrix(m$P)) - 1)), 1e-12)
})

# ---------------------------------------------------------------- item 3
test_that("MSVARXmdl returns the row-trimmed Z its simulator needs", {
  skip_on_cran()
  set.seed(9210)
  Z <- matrix(rnorm(200), ncol = 1)
  s <- simuMSVARX(list(n = 200, k = 2, q = 2, p = 1, mu = rbind(c(0, 0), c(2, 1)),
                       sigma = list(diag(2), 2 * diag(2)),
                       phi = matrix(c(.3, .05, .05, .25), 2, 2), P = P2,
                       betaZ = matrix(c(-0.8, 0.5), 1, 2), Z = Z), burnin = 100)
  set.seed(9211)
  m <- MSVARXmdl(as.matrix(s$y), p = 1, k = 2, Z = Z,
                 control = list(getSE = FALSE, use_diff_init = 1))
  expect_false(is.null(m$Z))
  # same convention as the other exogenous classes, and the length the simulator assumes
  expect_equal(as.numeric(m$Z), as.numeric(Z[2:200, , drop = FALSE]))
  expect_identical(nrow(m$Z), m$n)
  expect_identical(nrow(m$Z), nrow(m$y))
  # the Z term enters the simulated series at the right rows: kill the dynamics and the
  # noise and what is left is the demeaned Z times betaZ
  e <- m
  e$mu <- rbind(c(0, 0), c(0, 0)); e$phi <- matrix(0, 2, 2); e$P <- P2
  e$sigma <- list(diag(2) * 1e-12, diag(2) * 1e-12)
  e$betaZ <- matrix(c(3, -2), 1, 2)
  y <- as.matrix(simuMSVARX_cpp(e, 100)$y)
  Zdm <- m$Z - matrix(colMeans(m$Z), nrow(m$Z), ncol(m$Z), byrow = TRUE)
  expect_equal(y[-1, ], (Zdm %*% e$betaZ)[-1, ], tolerance = 1e-4)
  # and the whole k0 >= 2 exogenous null path runs
  a <- suppressWarnings(suppressMessages(
    LMCLRTest(as.matrix(s$y), p = 1, k0 = 2, k1 = 3, Z = Z,
              control = list(N = 5, mc_seed = 11, workers = 0,
                             mdl_h0_control = list(getSE = FALSE, use_diff_init = 1),
                             mdl_h1_control = list(getSE = FALSE, use_diff_init = 1)))))
  expect_true(is.finite(a$LRT_0))
  expect_true(a$pval >= 0 && a$pval <= 1)
})

# ---------------------------------------------------------------- item 4
test_that("a univariate hidden-Markov null is edited into the shape simuHMM_cpp reads", {
  skip_on_cran()
  set.seed(9220)
  s <- simuHMM(list(n = 200, k = 2, q = 1, mu = matrix(c(0, 2.5), 2, 1),
                    sigma = list(matrix(1, 1, 1), matrix(1, 1, 1)), P = P2), burnin = 100)
  y <- matrix(s$y, ncol = 1)
  for (msv in c(TRUE, FALSE)) {
    set.seed(9221)
    m <- HMmdl(y, k = 2, control = list(getSE = FALSE, use_diff_init = 1, msvar = msv))
    ed <- mdledit(m, m$theta, p = 0, q = 1, k0 = 2, exog = FALSE)
    # simuHMM_cpp reads sigma as a List of (q x q) matrices for every q, not as a vector
    expect_type(ed$sigma, "list")
    expect_length(ed$sigma, 2L)
    expect_equal(dim(as.matrix(ed$sigma[[1]])), c(1L, 1L))
    expect_equal(as.numeric(unlist(ed$sigma)),
                 as.numeric(if (msv) m$theta[as.numeric(m$theta_sig_ind) == 1]
                            else rep(m$theta[as.numeric(m$theta_sig_ind) == 1], 2)))
    e <- ed; e$n <- 100L
    expect_silent(simuHMM_cpp(e, 20, FALSE))
    # simuMdl routes (k > 1, p == 0) here for every q, so the same list must work there
    expect_silent(simuMdl(e, 0L, 1L, 2L, 20L, FALSE))
  }
  # the p > 0 branch keeps the vector convention simuMSAR_cpp reads
  set.seed(9222)
  sa <- simuMSAR(list(n = 200, k = 2, mu = c(-1, 2), sigma = c(1, 2), phi = 0.4, P = P2),
                 burnin = 100)
  set.seed(9223)
  ma <- MSARmdl(matrix(sa$y, ncol = 1), p = 1, k = 2,
                control = list(getSE = FALSE, use_diff_init = 1))
  eda <- mdledit(ma, ma$theta, p = 1, q = 1, k0 = 2, exog = FALSE)
  expect_false(is.list(eda$sigma))
  expect_length(as.numeric(eda$sigma), 2L)
  # the multivariate p == 0 branch was already a list and stays one
  set.seed(9224)
  s2 <- simuHMM(list(n = 250, k = 2, q = 2, mu = rbind(c(-1, -1), c(2, 2)),
                     sigma = list(diag(2), 2 * diag(2)), P = P2), burnin = 100)
  set.seed(9225)
  m2 <- HMmdl(s2$y, k = 2, control = list(getSE = FALSE, use_diff_init = 1))
  expect_type(mdledit(m2, m2$theta, p = 0, q = 2, k0 = 2, exog = FALSE)$sigma, "list")
  # and the whole univariate hidden-Markov MMC path runs
  a <- suppressWarnings(suppressMessages(
    MMCLRTest(y, p = 0, k0 = 2, k1 = 3,
              control = list(N = 5, mc_seed = 12, workers = 0, maxit = 2, type = "GenSA",
                             threshold_stop = 1, silence = TRUE,
                             mdl_h0_control = list(getSE = TRUE, use_diff_init = 1),
                             mdl_h1_control = list(getSE = FALSE, use_diff_init = 1)))))
  expect_true(a$pval >= 0 && a$pval <= 1)
})

# ---------------------------------------------------------------- item 5
test_that("simulated hidden-Markov draws carry the exogenous regressors", {
  skip_on_cran()
  set.seed(9230)
  s <- simuHMM(list(n = 200, k = 2, q = 1, mu = matrix(c(0, 2.5), 2, 1),
                    sigma = list(matrix(1, 1, 1), matrix(1, 1, 1)), P = P2), burnin = 100)
  Z <- matrix(rnorm(200), ncol = 1)
  set.seed(9231)
  m <- HMmdl(matrix(s$y, ncol = 1) + Z * (-0.8), k = 2, Z = Z,
             control = list(getSE = FALSE, use_diff_init = 2))
  set.seed(77); with_z <- simuMdl(m, 0L, 1L, 2L, 100L, TRUE)
  set.seed(77); no_z   <- simuMdl(m, 0L, 1L, 2L, 100L, FALSE)
  d <- as.numeric(with_z$y) - as.numeric(no_z$y)
  Zdm <- m$Z - matrix(colMeans(m$Z), nrow(m$Z), ncol(m$Z), byrow = TRUE)
  # exog = TRUE must reach simuHMM_cpp; the difference is exactly the demeaned Z term
  expect_gt(max(abs(d)), 1e-6)
  expect_equal(d, as.numeric(Zdm %*% m$betaZ))
})

# ---------------------------------------------------------------- item 7
test_that("the simulated-draw controls do not leak use_diff_init into a linear null", {
  skip_on_cran()
  set.seed(9240)
  s <- simuMSAR(list(n = 200, k = 2, mu = c(-1, 1.5), sigma = c(1, 1), phi = 0.3, P = P2),
                burnin = 100)
  y <- matrix(s$y, ncol = 1)
  count <- function(k0, k1) {
    n <- 0
    withCallingHandlers(
      suppressMessages(LMCLRTest(y, p = 1, k0 = k0, k1 = k1,
        control = list(N = 5, mc_seed = 13, workers = 0, use_diff_init_sim = 2,
                       mdl_h0_control = list(getSE = FALSE),
                       mdl_h1_control = list(getSE = FALSE)))),
      warning = function(w) {
        if (grepl("unknown names in control", conditionMessage(w))) n <<- n + 1
        invokeRestart("muffleWarning")
      })
    n
  }
  expect_identical(count(1, 2), 0)   # linear null: the option does not apply
  expect_identical(count(2, 3), 0)   # switching null: the option does apply, and is accepted
})

test_that("the retired BootLRTest methods are gone", {
  expect_length(as.character(utils::methods(class = "BootLRTest")), 0L)
  expect_false(exists("print.BootLRTest", where = asNamespace("MSTest"), inherits = FALSE))
  expect_false(exists("summary.BootLRTest", where = asNamespace("MSTest"), inherits = FALSE))
})

# ---------------------------------------------------------------- the guard
test_that("MMCLRTest reports when the search box admits no candidate", {
  skip_on_cran()
  set.seed(9250)
  s <- simuMSAR(list(n = 200, k = 2, mu = c(-1, 1.5), sigma = c(1, 1), phi = 0.3, P = P2),
                burnin = 100)
  y <- matrix(s$y, ncol = 1)
  # bounds that admit only non-stationary autoregressive parameters: no candidate the
  # optimizer can see is admissible, so there is no p-value to report
  for (ty in c("pso", "GenSA")) {
    expect_error(
      suppressWarnings(suppressMessages(MMCLRTest(y, p = 1, k0 = 1, k1 = 2,
        control = list(N = 5, mc_seed = 21, workers = 0, maxit = 8, type = ty,
                       threshold_stop = 1, eps = 2, phi_low = 1.05, phi_upp = 1.5,
                       silence = TRUE,
                       mdl_h0_control = list(getSE = FALSE),
                       mdl_h1_control = list(getSE = FALSE, use_diff_init = 1))))),
      "no admissible candidate", info = ty)
  }
  # and it does not fire under the default bounds, where theta_0 is always inside the box
  # (k0 = 2 is the tightest case: only theta_0 itself satisfies the column-sum constraint)
  for (ty in c("pso", "GenSA")) {
    a <- suppressWarnings(suppressMessages(MMCLRTest(y, p = 1, k0 = 2, k1 = 3,
      control = list(N = 5, mc_seed = 22, workers = 0, maxit = 12, type = ty,
                     threshold_stop = 1, silence = TRUE,
                     mdl_h0_control = list(getSE = TRUE, use_diff_init = 1),
                     mdl_h1_control = list(getSE = FALSE, use_diff_init = 1)))))
    expect_true(a$pval >= 0 && a$pval <= 1, info = ty)
    expect_gte(a$pval, a$pval_0)
  }
})
