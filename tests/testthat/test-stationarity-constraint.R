# Tests for the stationarity margin on estimation (control
# stationarity_constraint) and the linear constructors' inert acceptance of
# forwarded Markov-switching controls.

rho_max <- function(phi){
  p <- length(phi)
  if (p == 1) return(abs(phi))
  F <- matrix(0, p, p); F[1, ] <- phi; F[2:p, 1:(p-1)] <- diag(p-1)
  max(Mod(eigen(F, only.values = TRUE)$values))
}

set.seed(61)
y_ur <- as.matrix(cumsum(rnorm(300)))   # unit-root series: OLS phi_hat near 1

test_that("delta = 0 reproduces the unconstrained fits bitwise", {
  f0 <- ARmdl(y_ur, p = 1, control = list(getSE = TRUE))
  f1 <- ARmdl(y_ur, p = 1, control = list(getSE = TRUE, stationarity_constraint = 0))
  expect_identical(f0$theta, f1$theta)
  expect_identical(f0$logLike, f1$logLike)
  expect_identical(f0$theta_se, f1$theta_se)
})

test_that("ARmdl p = 1: exact constrained fit (clip + exact intercept re-solve)", {
  skip_on_cran()
  f0  <- ARmdl(y_ur, p = 1, control = list(getSE = FALSE))
  del <- 1 - as.numeric(f0$phi) + 0.02      # bound below phi_hat
  f   <- ARmdl(y_ur, p = 1, control = list(getSE = TRUE, stationarity_constraint = del))
  expect_identical(as.numeric(f$phi), 1 - del)     # the clip, exactly
  expect_true(isTRUE(f$phi_constraint_binding))
  # the intercept is the exact conditional least-squares solution at phi*
  yv <- as.numeric(f$y); x1 <- as.numeric(f$x)
  expect_equal(as.numeric(f$intercept), mean(yv - (1 - del) * x1), tolerance = 1e-12)
  # exactness against a fine grid of the concentrated sum of squares
  ssr <- function(phi){ r <- yv - phi * x1; r <- r - mean(r); sum(r^2) }
  grid <- seq(-(1 - del), 1 - del, length.out = 20001)
  expect_lte(ssr(as.numeric(f$phi)), min(vapply(grid, ssr, 0)) + 1e-9)
  # boundary SEs are NA for the AR block, finite elsewhere
  phi_ind <- which(as.numeric(f$theta_phi_ind) == 1)
  expect_true(all(is.na(f$theta_se[phi_ind])))
  expect_true(any(is.finite(f$theta_se[-phi_ind])))
  expect_identical(f$stationarity_constraint, del)
  # the constraint is scale-invariant: the same series at a large level
  fL <- ARmdl(y_ur + 5000, p = 1, control = list(getSE = FALSE, stationarity_constraint = del))
  expect_identical(as.numeric(fL$phi), 1 - del)
})

test_that("ARXmdl p = 1: constrained fit re-solves intercept AND betaZ exactly", {
  skip_on_cran()
  set.seed(62)
  Zx <- as.matrix(rnorm(300))
  yx <- y_ur + 0.5 * Zx
  f0 <- ARXmdl(yx, p = 1, Z = Zx, control = list(getSE = FALSE))
  del <- max(0.02, 1 - as.numeric(f0$phi) + 0.02)
  f  <- ARXmdl(yx, p = 1, Z = Zx, control = list(getSE = FALSE, stationarity_constraint = del))
  expect_identical(as.numeric(f$phi), 1 - del)
  yv <- as.numeric(f$y); x1 <- as.numeric(f$x)
  W  <- cbind(1, f$Z)
  gam <- solve(crossprod(W), crossprod(W, yv - (1 - del) * x1))
  expect_equal(as.numeric(f$intercept), gam[1], tolerance = 1e-12)
  expect_equal(as.numeric(f$betaZ), gam[-1], tolerance = 1e-12)
})

test_that("unidentified AR block: constraint declined with a warning, binding flag NA", {
  skip_on_cran()
  set.seed(63)
  yy <- as.matrix(cumsum(rnorm(200)))
  Zlag <- as.matrix(c(0, yy[-nrow(yy)])) + 1e-6 * rnorm(200)   # the lagged series itself
  expect_warning(
    fz <- ARXmdl(yy, p = 1, Z = Zlag, control = list(getSE = FALSE, stationarity_constraint = 0.5)),
    "collinear")
  expect_true(is.na(fz$phi_constraint_binding))
})

test_that("ARmdl p = 3: projection lands inside the margin", {
  skip_on_cran()
  set.seed(64)
  y3 <- as.matrix(cumsum(rnorm(400)) + 0.2 * as.numeric(arima.sim(list(ar = c(0.3, 0.2)), n = 400)))
  f0 <- ARmdl(y3, p = 3, control = list(getSE = FALSE))
  del <- max(0.05, 1 - rho_max(as.numeric(f0$phi)) + 0.03)
  f  <- ARmdl(y3, p = 3, control = list(getSE = FALSE, stationarity_constraint = del))
  expect_lte(rho_max(as.numeric(f$phi)), (1 - del) + 1e-5)
  expect_true(isTRUE(f$phi_constraint_binding))
})

test_that("projection helper: real rebuild, only offending roots moved, degenerate clusters", {
  phi_A <- c(2.85, -2.7075, 0.857375)     # triple root 0.95
  phi_B <- c(-2.85, -2.7075, -0.857375)   # triple root -0.95
  mid   <- (phi_A + phi_B) / 2
  expect_gt(rho_max(mid), 1)              # nonconvex: midpoint explosive
  pr <- as.numeric(phi_margin_project_cpp(mid, 0.05))
  expect_lte(rho_max(pr), 0.95 + 1e-5)
  phi_in <- c(0.4, 0.2, 0.1)
  expect_identical(as.numeric(phi_margin_project_cpp(phi_in, 1 - rho_max(phi_in) - 0.01)), phi_in)
  expect_identical(as.numeric(phi_margin_project_cpp(1.05, 0.02)), 0.98)
  # an m-fold root at 1.02, m = 4..6: the coefficient-scaling fallback holds
  for (m in 4:6){
    cf <- Re(-Reduce(function(acc, l) c(acc, 0) - c(0, acc) * l, rep(1.02, m), init = 1)[-1])
    pm <- as.numeric(phi_margin_project_cpp(cf, 0.05))
    expect_lte(rho_max(pm), 0.95 + 1e-5)
  }
})

test_that("MSARmdl EM p = 1: constrained fit respects the margin, binding flagged, SEs NA", {
  skip_on_cran()
  set.seed(65)
  f <- suppressWarnings(MSARmdl(y_ur, p = 1, k = 2,
        control = list(getSE = TRUE, use_diff_init = 2, stationarity_constraint = 0.03)))
  expect_lte(abs(as.numeric(f$phi)), 0.97 + 1e-12)
  expect_true(is.finite(f$logLike))
  expect_identical(f$stationarity_constraint, 0.03)
  expect_true(isTRUE(f$phi_constraint_binding))
  phi_ind <- which(as.numeric(f$theta_phi_ind) == 1)
  expect_true(all(is.na(f$theta_se[phi_ind])))
})

test_that("MSARmdl EM p = 1: the M-step clip is the exact constrained maximizer", {
  skip_on_cran()
  set.seed(66)
  ysim <- simuMSAR(list(n = 300, k = 2, mu = c(0, 1.5), sigma = c(1, 1.5),
                        phi = 0.97, P = matrix(c(0.95, 0.05, 0.05, 0.95), 2, 2)))
  y <- as.matrix(ysim$y)
  delta <- 0.05
  f <- suppressWarnings(MSARmdl(y, p = 1, k = 2, control = list(getSE = FALSE, use_diff_init = 2,
                                                                stationarity_constraint = delta)))
  mdl <- f; mdl$phi_margin <- delta
  th <- as.numeric(f$theta)
  E <- ExpectationM_MSARmdl(th, mdl, 2)
  M <- EMaximization_MSARmdl(th, mdl, E, 2)
  g <- MSTest:::argrid_MSARmdl_cpp(as.numeric(M$mu), as.numeric(E$sigma), 2, 1)
  yv <- as.numeric(f$y); X <- as.matrix(f$x); xi <- as.matrix(E$xi_t_T_AR)
  den <- 0; num <- 0
  for (m in seq_len(nrow(g$mu))){
    if (g$sig[m] < 1e-14) next
    yt <- yv - g$mu[m, 1]; xt <- X - g$mu[m, 2]; w <- xi[, m] / g$sig[m]
    den <- den + sum(w * xt^2); num <- num + sum(w * xt * yt)
  }
  phi_clip <- min(max(num / den, -(1 - delta)), 1 - delta)
  expect_equal(as.numeric(M$phi), phi_clip, tolerance = 1e-12)
  Q <- function(phi) num * phi - 0.5 * den * phi^2
  grid <- seq(-(1 - delta), 1 - delta, length.out = 200001)
  expect_gte(Q(as.numeric(M$phi)), max(Q(grid)) - 1e-9)
})

test_that("MSARmdl EM: a nonstationary init_theta is projected under the margin, preserved at delta = 0", {
  skip_on_cran()
  set.seed(67)
  yb <- as.matrix(cumsum(rnorm(200)))
  f0 <- suppressWarnings(MSARmdl(yb, p = 1, k = 2, control = list(getSE = FALSE, use_diff_init = 1)))
  th <- as.numeric(f0$theta); th[which(as.numeric(f0$theta_phi_ind) == 1)] <- 1.05
  fL0a <- suppressWarnings(MSARmdl(yb, p = 1, k = 2,
           control = list(getSE = FALSE, use_diff_init = 1, init_theta = th)))
  fL0b <- suppressWarnings(MSARmdl(yb, p = 1, k = 2,
           control = list(getSE = FALSE, use_diff_init = 1, init_theta = th,
                          stationarity_constraint = 0)))
  expect_identical(fL0a$theta, fL0b$theta)
  expect_identical(fL0a$logLike, fL0b$logLike)
  fC <- suppressWarnings(MSARmdl(yb, p = 1, k = 2,
         control = list(getSE = FALSE, use_diff_init = 1, init_theta = th,
                        stationarity_constraint = 0.03)))
  expect_lte(abs(as.numeric(fC$phi)), 0.97 + 1e-12)
})

test_that("MSARXmdl EM p = 1: the constrained M-step re-solves betaZ at the clipped phi", {
  skip_on_cran()
  set.seed(85)                            # the constrained fit binds on this draw
  Tn <- 300
  Zx <- as.matrix(rnorm(Tn))
  yx <- as.matrix(cumsum(rnorm(Tn))) + 0.6 * Zx
  delta <- 0.05
  f <- suppressWarnings(MSARXmdl(yx, p = 1, k = 2, Z = Zx,
         control = list(getSE = FALSE, use_diff_init = 3, stationarity_constraint = delta)))
  expect_true(isTRUE(f$phi_constraint_binding))
  # one E/M step from the fitted theta, with the M-step's own objects rebuilt
  mdl <- f; mdl$phi_margin <- delta
  th <- as.numeric(f$theta)
  E <- ExpectationM_MSARXmdl(th, mdl, 2)
  M <- EMaximization_MSARXmdl(th, mdl, E, 2)
  yv <- as.numeric(f$y); X <- as.matrix(f$x)
  Zsl <- as.matrix(f$control$Z)[2:Tn, , drop = FALSE]
  zdm <- sweep(Zsl, 2, colMeans(Zsl))
  g <- MSTest:::argrid_MSARmdl_cpp(as.numeric(M$mu), as.numeric(E$sigma), 2, 1)
  xi <- as.matrix(E$xi_t_T_AR)
  denom <- matrix(0, 2, 2); num <- rep(0, 2)
  for (m in seq_len(nrow(g$mu))){
    if (g$sig[m] < 1e-14) next
    yt <- yv - g$mu[m, 1]; xt <- X - g$mu[m, 2]
    xz <- cbind(xt, zdm); w <- xi[, m] / g$sig[m]
    denom <- denom + t(xz) %*% (w * xz); num <- num + as.numeric(t(xz) %*% (w * yt))
  }
  b_unc <- solve(denom, num)
  phi_star <- min(max(b_unc[1], -(1 - delta)), 1 - delta)
  betaZ_star <- (num[2] - denom[2, 1] * phi_star) / denom[2, 2]
  expect_equal(as.numeric(M$phi), phi_star, tolerance = 1e-10)
  expect_equal(as.numeric(M$betaZ), betaZ_star, tolerance = 1e-10)
  # keeping the jointly-solved betaZ at the clipped phi would lower the
  # M-step objective (the ascent violation the re-solve prevents)
  Qj <- function(b) sum(num * b) - 0.5 * as.numeric(t(b) %*% denom %*% b)
  expect_gt(Qj(c(phi_star, betaZ_star)), Qj(c(phi_star, b_unc[2])) - 1e-12)
})

test_that("MLE path: the margin imposes the root bound, with or without mle_stationary_constraint", {
  skip_on_cran()
  set.seed(69)
  ym <- as.matrix(cumsum(rnorm(200)))
  for (msc in c(TRUE, FALSE)){
    fm <- suppressWarnings(try(MSARmdl(ym, p = 1, k = 2,
           control = list(getSE = FALSE, use_diff_init = 2, method = "MLE",
                          mle_stationary_constraint = msc, stationarity_constraint = 0.05)), silent = TRUE))
    skip_if(inherits(fm, "try-error"), "MLE did not converge on this draw")
    expect_lte(abs(as.numeric(fm$phi)), 0.95 + 1e-6)
  }
})

test_that("multivariate classes error loudly and validate; p = 0 classes accept the name silently", {
  set.seed(70)
  y2 <- matrix(rnorm(200), 100, 2)
  expect_error(VARmdl(y2, p = 1, control = list(stationarity_constraint = 0.05)),
               "not implemented for multivariate")
  expect_error(VARXmdl(y2, p = 1, Z = as.matrix(rnorm(100)), control = list(stationarity_constraint = 0.05)),
               "not implemented for multivariate")
  expect_error(suppressWarnings(MSVARmdl(y2, p = 1, k = 2,
               control = list(use_diff_init = 1, stationarity_constraint = 0.05))),
               "not implemented for multivariate")
  expect_error(VARmdl(y2, p = 1, control = list(stationarity_constraint = c(0.1, 0.1))),
               "single finite number")
  expect_silent(VARmdl(y2, p = 1, control = list(getSE = FALSE, stationarity_constraint = 0)))
  expect_silent(Nmdl(as.matrix(rnorm(100)), control = list(getSE = FALSE, stationarity_constraint = 0.05,
                                                            use_diff_init = 2)))
  expect_silent(suppressMessages(HMmdl(as.matrix(rnorm(100)), k = 2,
                control = list(getSE = FALSE, use_diff_init = 1, stationarity_constraint = 0.05))))
})

test_that("validation: malformed delta errors by name", {
  expect_error(ARmdl(y_ur, p = 1, control = list(stationarity_constraint = 1)),
               "less than 1")
  expect_error(ARmdl(y_ur, p = 1, control = list(stationarity_constraint = NA)),
               "finite")
})

test_that("linear constructors accept forwarded MS controls silently; typos still warn", {
  expect_silent(ARmdl(y_ur, p = 1, control = list(getSE = FALSE, use_diff_init = 5,
                       em_variance_constraint = 0.01, em_transition_constraint = 0.02)))
  expect_warning(ARmdl(y_ur, p = 1, control = list(getSE = FALSE, use_dif_init = 5)),
                 "unknown names")
  set.seed(71)
  y2 <- matrix(rnorm(200), 100, 2)
  expect_silent(VARmdl(y2, p = 1, control = list(getSE = FALSE, use_diff_init = 5)))
})

test_that("LMCLRTest end-to-end: no unknown-name warnings, mismatch warned, floor threaded", {
  skip_on_cran()
  set.seed(72)
  yb <- as.matrix(cumsum(rnorm(200)))
  ctl <- list(N = 19, mc_seed = 73,
              mdl_h0_control = list(stationarity_constraint = 0.02, use_diff_init = 2),
              mdl_h1_control = list(stationarity_constraint = 0.02, use_diff_init = 2))
  warns <- character(0)
  l <- withCallingHandlers(
    LMCLRTest(yb, p = 1, k0 = 1, k1 = 2, control = ctl),
    warning = function(w){ warns <<- c(warns, conditionMessage(w)); invokeRestart("muffleWarning") })
  expect_false(any(grepl("unknown names in control", warns)))   # the old storm is gone
  expect_true(is.finite(l$pval) && l$pval >= 0 && l$pval <= 1)
  expect_lte(abs(as.numeric(l$mdl_h0$phi)), 0.98 + 1e-12)
  # p = 0 null with the control on both sides: no warning storm either
  warns0 <- character(0)
  l0 <- withCallingHandlers(
    LMCLRTest(as.matrix(rnorm(120)), p = 0, k0 = 1, k1 = 2,
              control = list(N = 9, mc_seed = 74,
                             mdl_h0_control = list(stationarity_constraint = 0.05),
                             mdl_h1_control = list(stationarity_constraint = 0.05, use_diff_init = 1))),
    warning = function(w){ warns0 <<- c(warns0, conditionMessage(w)); invokeRestart("muffleWarning") })
  expect_false(any(grepl("unknown names in control", warns0)))
  # mismatched margins on the two sides warn
  expect_warning(
    suppressMessages(LMCLRTest(yb, p = 1, k0 = 1, k1 = 2,
      control = list(N = 9, mc_seed = 75,
                     mdl_h0_control = list(stationarity_constraint = 0.02, use_diff_init = 1),
                     mdl_h1_control = list(stationarity_constraint = 0.10, use_diff_init = 1)))),
    "differs between")
  # floor threading (k0 > 1), observed side: the null log-likelihood implied
  # by the reported statistic equals max(MS fit, CONSTRAINED linear fit)
  # exactly on a draw where the floor binds; an unconstrained floor would sit
  # strictly above the constrained linear value
  set.seed(103)                          # OLS phi 1.002: the floor binds here
  yb2 <- as.matrix(cumsum(rnorm(150)))
  l2 <- suppressWarnings(LMCLRTest(yb2, p = 1, k0 = 2, k1 = 3,
         control = list(N = 19, mc_seed = 104,
                        mdl_h0_control = list(stationarity_constraint = 0.02, use_diff_init = 1),
                        mdl_h1_control = list(stationarity_constraint = 0.02, use_diff_init = 1))))
  expect_true(is.finite(l2$pval))
  lin_c <- ARmdl(yb2, p = 1, control = list(getSE = FALSE, stationarity_constraint = 0.02))
  lin_u <- ARmdl(yb2, p = 1, control = list(getSE = FALSE))
  skip_if(lin_u$logLike <= lin_c$logLike + 1e-8, "floor does not discriminate on this draw")
  # On this seeded draw the constrained floor keeps the statistic positive; an
  # unconstrained floor would cap it at 0, so this is a hard assertion, not a skip.
  expect_gt(as.numeric(l2$LRT_0), 0)
  implied_logL0 <- l2$mdl_h1$logLike - as.numeric(l2$LRT_0) / 2
  expect_equal(implied_logL0, max(l2$mdl_h0$logLike, lin_c$logLike), tolerance = 1e-8)
  expect_lt(implied_logL0, lin_u$logLike)
})

test_that("parallel path with the constraint: deterministic under mc_seed, no warnings", {
  skip_on_cran()
  set.seed(76)
  yb <- as.matrix(cumsum(rnorm(150)))
  ctl <- list(N = 19, mc_seed = 77, workers = 2,
              mdl_h0_control = list(stationarity_constraint = 0.02, use_diff_init = 1),
              mdl_h1_control = list(stationarity_constraint = 0.02, use_diff_init = 1))
  warns <- character(0)
  a <- withCallingHandlers(suppressMessages(LMCLRTest(yb, p = 1, k0 = 1, k1 = 2, control = ctl)),
         warning = function(w){ warns <<- c(warns, conditionMessage(w)); invokeRestart("muffleWarning") })
  b <- suppressWarnings(suppressMessages(LMCLRTest(yb, p = 1, k0 = 1, k1 = 2, control = ctl)))
  expect_false(any(grepl("unknown names in control", warns)))
  expect_identical(a$pval, b$pval)
  expect_identical(a$LRN, b$LRN)
  expect_lte(abs(as.numeric(a$mdl_h0$phi)), 0.98 + 1e-12)
})

test_that("the alternative inherits the null's margin when its control is silent; an explicit mismatch warns without inheriting", {
  skip_on_cran()
  set.seed(21); y <- as.matrix(cumsum(rnorm(200)))
  w <- character(0)
  r1 <- withCallingHandlers(
    LMCLRTest(y, p = 1, k0 = 1, k1 = 2, control = list(N = 19, mc_seed = 77,
              mdl_h0_control = list(stationarity_constraint = 0.05, getSE = FALSE),
              mdl_h1_control = list(getSE = FALSE, use_diff_init = 2))),
    warning = function(cnd){ w <<- c(w, conditionMessage(cnd)); invokeRestart("muffleWarning") })
  expect_equal(r1$mdl_h1$stationarity_constraint, 0.05)
  expect_true(any(grepl("has been applied to the alternative", w)))
  expect_false(any(grepl("differs between", w)))
  expect_lte(MSTest:::.phi_companion_radius(as.numeric(r1$mdl_h1$phi)), 0.95 + 1e-12)
  w <- character(0)
  r2 <- withCallingHandlers(
    LMCLRTest(y, p = 1, k0 = 1, k1 = 2, control = list(N = 19, mc_seed = 77,
              mdl_h0_control = list(stationarity_constraint = 0.05, getSE = FALSE),
              mdl_h1_control = list(stationarity_constraint = 0, getSE = FALSE, use_diff_init = 2))),
    warning = function(cnd){ w <<- c(w, conditionMessage(cnd)); invokeRestart("muffleWarning") })
  expect_equal(r2$mdl_h1$stationarity_constraint, 0)
  expect_true(any(grepl("differs between", w)))
  expect_false(any(grepl("has been applied to the alternative", w)))
})

test_that("margin gate helper: candidates beyond the radius bound are penalized before the objective runs", {
  calls <- 0L
  f <- function(theta, ...){ calls <<- calls + 1L; 1 }
  env <- new.env(); env$margin_hits <- 0L
  g <- MSTest:::.mmc_margin_wrap(f, phi_idx = 2:3, radius_bound = 0.95, penalty_value = 100, count_env = env)
  expect_equal(g(c(0, 0.5, 0.3, 1)), 1)               # radius 0.86 < 0.95: passes through
  expect_equal(g(c(0, 0.9, 0.09, 1)), 100)            # radius 0.9877 > 0.95: penalized
  expect_equal(g(c(0, 1.9, -0.9, 1)), 100)            # roots 1.0 and 0.9: radius 1.0 > 0.95
  expect_equal(g(c(0, NaN, 0.1, 1)), 100)             # non-finite candidate: penalized, never evaluated
  expect_equal(calls, 1L); expect_equal(env$margin_hits, 3L)
})

test_that("MMCLRTest: box_match_constraints clamps the search to the margin at p = 1; FALSE warns and leaves the box", {
  skip_on_cran()
  set.seed(21); y <- as.matrix(cumsum(rnorm(250)))   # LRT_0 > 0 on this draw, so the search runs
  base <- list(N = 19, mc_seed = 77, maxit = 3, type = "GenSA", eps = 0.3, CI_union = FALSE, threshold_stop = 1, silence = TRUE,
               mdl_h0_control = list(stationarity_constraint = 0.05, getSE = FALSE),
               mdl_h1_control = list(stationarity_constraint = 0.05, getSE = FALSE, use_diff_init = 2))
  m <- suppressWarnings(MMCLRTest(y, p = 1, k0 = 1, k1 = 2, control = base))
  pp <- which(as.numeric(m$mdl_h0$theta_phi_ind) == 1)
  expect_true(isTRUE(m$mdl_h0$phi_constraint_binding))
  expect_equal(as.numeric(m$theta_upp[pp]), 1 - 0.05, tolerance = 0)
  expect_gte(as.numeric(m$theta_low[pp]), -(1 - 0.05))
  expect_false(is.null(m$optim_diagnostics$evaluated_range))
  expect_lte(as.numeric(m$optim_diagnostics$evaluated_range["max", pp]), 1 - 0.05 + 1e-12)
  expect_true(is.integer(m$optim_diagnostics$margin_penalty_hits))
  w <- character(0)
  m0 <- withCallingHandlers(MMCLRTest(y, p = 1, k0 = 1, k1 = 2, control = c(base, list(box_match_constraints = FALSE))),
          warning = function(cnd){ w <<- c(w, conditionMessage(cnd)); invokeRestart("muffleWarning") })
  expect_true(any(grepl("may extend beyond the constrained parameter space", w)))
  expect_gt(as.numeric(m0$theta_upp[pp]), 1 - 0.05)
  expect_true(is.na(m0$optim_diagnostics$margin_penalty_hits))
  # without any estimation constraint the control is inert: no warning, ordinary box
  m1 <- suppressWarnings(MMCLRTest(y, p = 1, k0 = 1, k1 = 2, control = list(N = 19, mc_seed = 77, maxit = 3, type = "GenSA", eps = 0.3, CI_union = FALSE, threshold_stop = 1, silence = TRUE,
                                   mdl_h0_control = list(getSE = FALSE), mdl_h1_control = list(getSE = FALSE, use_diff_init = 2))))
  expect_equal(as.numeric(m1$theta_upp[pp]), as.numeric(m1$mdl_h0$phi) + 0.3)
  expect_true(is.na(m1$optim_diagnostics$margin_penalty_hits))
})

test_that("MMCLRTest: a transition-probability bound clamps the P box to the bounded simplex; contradictions error", {
  skip_on_cran()
  set.seed(31)
  sim <- simuMSAR(list(n = 150, k = 2, mu = c(-1, 1), sigma = c(1, 1), phi = 0.3, P = matrix(c(0.9, 0.1, 0.1, 0.9), 2, 2), burnin = 50))
  y <- as.matrix(sim$y)
  ctl <- list(N = 19, mc_seed = 5, maxit = 2, type = "GenSA", eps = 0.3, CI_union = FALSE, threshold_stop = 1, silence = TRUE,
              mdl_h0_control = list(em_transition_constraint = 0.05, getSE = FALSE, use_diff_init = 2),
              mdl_h1_control = list(em_transition_constraint = 0.05, getSE = FALSE, use_diff_init = 2))
  m <- suppressWarnings(MMCLRTest(y, p = 1, k0 = 2, k1 = 3, control = ctl))
  Pp <- which(as.numeric(m$mdl_h0$theta_P_ind) == 1)
  expect_true(all(as.numeric(m$theta_low[Pp]) >= 0.05))
  expect_true(all(as.numeric(m$theta_upp[Pp]) <= 1 - 0.05))
  w <- character(0)
  m0 <- withCallingHandlers(MMCLRTest(y, p = 1, k0 = 2, k1 = 3, control = c(ctl, list(box_match_constraints = FALSE))),
          warning = function(cnd){ w <<- c(w, conditionMessage(cnd)); invokeRestart("muffleWarning") })
  expect_true(any(grepl("may extend beyond the constrained parameter space", w)))
  expect_true(any(as.numeric(m0$theta_low[Pp]) < 0.05))
  expect_error(MMCLRTest(y, p = 1, k0 = 2, k1 = 3, control = list(box_boundary = "union", burnin = 0,
                         mdl_h0_control = list(em_transition_constraint = 0.05))), "P_low = 0")
  expect_error(MMCLRTest(y, p = 1, k0 = 1, k1 = 2, control = list(stationary_constraint = FALSE,
                         mdl_h0_control = list(stationarity_constraint = 0.05))), "contradictory")
  expect_error(MMCLRTest(y, p = 1, k0 = 1, k1 = 2, control = list(box_match_constraints = NA)), "TRUE or FALSE")
})

test_that("augmented tests: the alternative inherits the null's margin, mismatches warn, and the observed floor is threaded", {
  skip_on_cran()
  set.seed(103)                          # OLS phi 1.002: the floor binds on this draw
  yb2 <- as.matrix(cumsum(rnorm(150)))
  w <- character(0)
  a <- withCallingHandlers(
    aURfixLMCTest(yb2, p = 1, k0 = 2, k1 = 3,
                  control = list(N = 19, mc_seed = 104,
                                 mdl_h0_control = list(stationarity_constraint = 0.02, use_diff_init = 1),
                                 mdl_h1_control = list(use_diff_init = 1))),
    warning = function(cnd){ w <<- c(w, conditionMessage(cnd)); invokeRestart("muffleWarning") })
  expect_equal(a$mdl_h1$stationarity_constraint, 0.02)
  expect_true(any(grepl("has been applied to the alternative", w)))
  expect_false(any(grepl("differs between", w)))
  lin_c <- ARmdl(yb2, p = 1, control = list(getSE = FALSE, stationarity_constraint = 0.02))
  lin_u <- ARmdl(yb2, p = 1, control = list(getSE = FALSE))
  skip_if(lin_u$logLike <= lin_c$logLike + 1e-8, "floor does not discriminate on this draw")
  expect_gt(as.numeric(a$LRT_0), 0)   # hard: an unconstrained floor would cap the statistic at 0 on this draw
  implied_logL0 <- a$mdl_h1$logLike - as.numeric(a$LRT_0) / 2
  expect_equal(implied_logL0, max(a$mdl_h0$logLike, lin_c$logLike), tolerance = 1e-8)
  expect_lt(implied_logL0, lin_u$logLike)
  # component 1 stays LMCLRTest under the same controls (inheritance included)
  l <- suppressWarnings(LMCLRTest(yb2, p = 1, k0 = 2, k1 = 3,
         control = list(N = 19, mc_seed = 104,
                        mdl_h0_control = list(stationarity_constraint = 0.02, use_diff_init = 1),
                        mdl_h1_control = list(use_diff_init = 1))))
  expect_identical(a$pval_lmc, l$pval)
  expect_identical(as.numeric(a$LRT_0), as.numeric(l$LRT_0))
  w <- character(0)
  b <- withCallingHandlers(
    aURfixLMCTest(yb2, p = 1, k0 = 2, k1 = 3,
                  control = list(N = 19, mc_seed = 104,
                                 mdl_h0_control = list(stationarity_constraint = 0.02, use_diff_init = 1),
                                 mdl_h1_control = list(stationarity_constraint = 0, use_diff_init = 1))),
    warning = function(cnd){ w <<- c(w, conditionMessage(cnd)); invokeRestart("muffleWarning") })
  expect_equal(b$mdl_h1$stationarity_constraint, 0)
  expect_true(any(grepl("differs between", w)))
})

test_that("MMCLRTest: at p = 2 the margin gate excludes candidates beyond the margin and the winner stays inside", {
  skip_on_cran()
  set.seed(12)   # AR(2) with a mean shift: the unconstrained fit has radius 0.956, so the constrained fit binds at 0.95 and the eps box straddles the margin
  y <- as.matrix(as.numeric(arima.sim(list(ar = c(0.75, 0.22)), n = 300)) + rep(c(0, 3), each = 150))
  ctl <- list(N = 19, mc_seed = 42, maxit = 8, type = "GenSA", eps = 0.3, CI_union = FALSE, threshold_stop = 1, silence = TRUE,
              mdl_h0_control = list(stationarity_constraint = 0.05, getSE = FALSE),
              mdl_h1_control = list(stationarity_constraint = 0.05, getSE = FALSE, use_diff_init = 2))
  m <- suppressWarnings(MMCLRTest(y, p = 2, k0 = 1, k1 = 2, control = ctl))
  pp <- which(as.numeric(m$mdl_h0$theta_phi_ind) == 1)
  expect_lte(MSTest:::.phi_companion_radius(as.numeric(m$mdl_h0$theta)[pp]), 0.95 + 1e-12)   # seed inside
  expect_true(is.integer(m$optim_diagnostics$margin_penalty_hits))
  expect_gt(m$optim_diagnostics$margin_penalty_hits, 0L)                                     # the box reaches beyond the margin
  expect_lte(MSTest:::.phi_companion_radius(as.numeric(m$theta_h0)[pp]), 0.95 + 1e-12)      # a penalized candidate never wins
  expect_gte(m$pval, m$pval_0)
  # at p = 1 the clamp makes the gate unreachable
  set.seed(21); y1 <- as.matrix(cumsum(rnorm(250)))
  m1 <- suppressWarnings(MMCLRTest(y1, p = 1, k0 = 1, k1 = 2, control = c(ctl[setdiff(names(ctl), "maxit")], list(maxit = 3))))
  expect_identical(m1$optim_diagnostics$margin_penalty_hits, 0L)
  # a box_include point outside the margin is contradictory
  binc <- rep(NA_real_, length(m$mdl_h0$theta)); binc[pp] <- c(0.5, 0.45)   # radius 0.966 > 0.95
  expect_error(MMCLRTest(y, p = 2, k0 = 1, k1 = 2, control = c(ctl, list(box_include = binc))), "contradictory")
})

test_that("MLE null with a binding margin: the fit sits inside the margin and MMCLRTest's aligned box accepts it; the P bound is not aligned for MLE", {
  skip_on_cran()
  set.seed(21); y <- as.matrix(cumsum(rnorm(150)))
  f <- suppressWarnings(MSARmdl(y, p = 1, k = 2, control = list(method = "MLE", stationarity_constraint = 0.05, getSE = FALSE, use_diff_init = 1)))
  expect_lte(MSTest:::.phi_companion_radius(as.numeric(f$phi)), 0.95)
  ctl <- list(N = 19, mc_seed = 7, maxit = 2, type = "GenSA", eps = 0.3, CI_union = FALSE, threshold_stop = 1, silence = TRUE,
              mdl_h0_control = list(method = "MLE", stationarity_constraint = 0.05, em_transition_constraint = 0.05, getSE = FALSE, use_diff_init = 1),
              mdl_h1_control = list(method = "MLE", stationarity_constraint = 0.05, getSE = FALSE, use_diff_init = 1))
  m <- suppressWarnings(MMCLRTest(y, p = 1, k0 = 2, k1 = 3, control = ctl))
  pp <- which(as.numeric(m$mdl_h0$theta_phi_ind) == 1)
  expect_equal(as.numeric(m$theta_upp[pp]), 1 - 0.05, tolerance = 0)
  Pp <- which(as.numeric(m$mdl_h0$theta_P_ind) == 1)
  expect_true(any(as.numeric(m$theta_low[Pp]) < 0.05))     # MLE ignores em_transition_constraint, so the P box is not clamped
})
