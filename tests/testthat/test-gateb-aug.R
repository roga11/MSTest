# Tests for the augmented (intersection-union) Monte Carlo LR tests:
# aURfixLMCTest, aBfixLMCTest, aURdiffLMCTest, aURdiffMMCTest, plus the ADF
# unit-root arithmetic and the p = 0 testing path they rely on.

valid_pval <- function(p) is.finite(p) && p >= 0 && p <= 1

set.seed(31)
y_ur <- as.matrix(cumsum(rnorm(400)))
set.seed(32)
y_b <- as.matrix(simuMSAR(list(n = 200, k = 2, mu = c(0, 2), sigma = c(1, 2),
                               phi = 0.3, P = matrix(c(0.9, 0.1, 0, 1), 2, 2)),
                          burnin = 0)$y)
ctl2 <- list(use_diff_init = 2)

test_that("ADF unit-root map: exact identities at p = 1..5", {
  adf <- MSTest:::.adf_unit_root
  expect_identical(adf(0.7), 1)
  for (p in 2:5){
    set.seed(p); phi <- runif(p, -0.2, 0.25)
    phin <- adf(phi)
    # long-run coefficient moved to exactly 1
    expect_lt(abs(sum(phin) - 1), 1e-12)
    # short-run (ADF) coefficients held fixed
    gam0 <- -rev(cumsum(rev(phi)))[-1]
    gam1 <- -rev(cumsum(rev(phin)))[-1]
    expect_lt(max(abs(gam0 - gam1)), 1e-12)
    # the AR polynomial has a unit root exactly: 1 - sum(phi~) = 0
    expect_lt(abs(1 - sum(phin)), 1e-12)
  }
  # the map moves ONLY the long-run coordinate: applying it to an already
  # unit-root vector is the identity
  phi_u <- adf(c(0.4, 0.3, 0.1))
  expect_equal(adf(phi_u), phi_u, tolerance = 1e-14)
})

test_that("aURfixLMCTest: component 1 is bitwise LMCLRTest; component 2 at the unit root", {
  skip_on_cran()
  ctl <- list(N = 19, mc_seed = 74, mdl_h1_control = ctl2)
  a <- suppressWarnings(aURfixLMCTest(y_ur, p = 1, k0 = 1, k1 = 2, control = ctl))
  l <- suppressWarnings(LMCLRTest(y_ur, p = 1, k0 = 1, k1 = 2, control = ctl))
  expect_identical(a$pval_lmc, l$pval)
  expect_identical(unname(a$LRT_0), unname(l$LRT_0))
  expect_identical(a$LRN, l$LRN)
  phi_pos <- which(as.numeric(a$mdl_h0$theta_phi_ind) == 1)
  expect_identical(a$theta_fix[[phi_pos]], 1)
  other <- setdiff(seq_along(a$theta_fix), phi_pos)
  expect_identical(unname(a$theta_fix[other]), unname(as.numeric(a$mdl_h0$theta)[other]))
  expect_true(valid_pval(a$pval_fix))
  expect_identical(a$pval, max(a$pval_lmc, a$pval_fix))
  # deterministic repeat
  a2 <- suppressWarnings(aURfixLMCTest(y_ur, p = 1, k0 = 1, k1 = 2, control = ctl))
  expect_identical(a$pval_fix, a2$pval_fix)
  # scope errors
  expect_error(aURfixLMCTest(y_ur, 0, 1, 2), "p >= 1")
})

test_that("aURfix component 2 replays the estimation stream (the shared-randomness contract)", {
  skip_on_cran()
  # Independent reimplementation of the documented choreography: seed ->
  # estimate -> re-seed -> predraw -> anchor -> component 1 (+ its MCpval) ->
  # restore anchor -> component 2 at theta_fix (+ its MCpval). If the
  # implementation ever stops restoring the anchor (so the components no
  # longer share estimation randomness) or moves the MCpval calls relative to
  # the stream, this test fails.
  ctl <- list(N = 19, mc_seed = 74, mdl_h1_control = ctl2)
  a <- suppressWarnings(aURfixLMCTest(y_ur, p = 1, k0 = 1, k1 = 2, control = ctl))
  seed_after_pkg <- get(".Random.seed", envir = globalenv())
  set.seed(74)
  mdl_h0 <- estimMdl(y_ur, 1, 1, 1, NULL, list())
  con_h1 <- MSTest:::.inherit_ms_flags(ctl2, mdl_h0)
  con_h1$init_theta_extra <- warmstart_theta(mdl_h0, 2, msmu = con_h1[["msmu"]],
                                             msvar = con_h1[["msvar"]])
  mdl_h1 <- estimMdl(y_ur, 1, 1, 2, NULL, con_h1)
  LRT_0 <- -2 * (mdl_h0$logLike - mdl_h1$logLike)
  if (LRT_0 < 0) LRT_0 <- 0
  theta_fix <- as.numeric(mdl_h0$theta)
  theta_fix[which(as.numeric(mdl_h0$theta_phi_ind) == 1)] <- 1
  set.seed(74)
  Nb <- ceiling(19 * 1.5) + 10; Teps <- mdl_h0$n + 100
  pe <- lapply(seq_len(Nb), function(i) matrix(rnorm(Teps), Teps, 1))
  h0c <- mdl_h0$control
  h1c <- mdl_h1$control; h1c$init_theta_extra <- NULL; h1c$init_method <- "warmstart"
  S <- get(".Random.seed", envir = globalenv())
  LRN1 <- suppressWarnings(LR_samp_dist(mdl_h0, 2, 19, 100, NULL, h0c, h1c, pe, NULL))
  LRN1 <- LRN1[is.finite(LRN1)]
  p1m <- MCpval(unname(LRT_0), LRN1, "geq")
  assign(".Random.seed", S, envir = globalenv())
  LRN2 <- suppressWarnings(LR_samp_dist(mdledit(mdl_h0, theta_fix, 1, 1, 1, FALSE),
                                        2, 19, 100, NULL, h0c, h1c, pe, NULL))
  LRN2 <- LRN2[is.finite(LRN2)]
  p2m <- MCpval(unname(LRT_0), LRN2, "geq")
  seed_after_manual <- get(".Random.seed", envir = globalenv())
  expect_identical(p1m, a$pval_lmc)
  expect_identical(p2m, a$pval_fix)
  # The stream-state pin is what actually discriminates a dropped restore: on
  # a clean configuration EM reaches the same optimum from any random start,
  # so the p-values alone can coincide even when the components no longer
  # share their estimation randomness. The implementation's post-call RNG
  # state must equal this choreography's.
  expect_identical(seed_after_pkg, seed_after_manual)
})

test_that("aBfixLMCTest: canonical labeling, exact boundary projection, composition", {
  skip_on_cran()
  a <- suppressWarnings(aBfixLMCTest(y_b, p = 1, k0 = 2, k1 = 3,
        control = list(N = 19, mc_seed = 75, burnin = 0,
                       mdl_h0_control = ctl2, mdl_h1_control = ctl2)))
  # canonical: absorbing regime last
  P_hat <- matrix(as.numeric(a$mdl_h0$theta)[as.numeric(a$mdl_h0$theta_P_ind) == 1], 2, 2)
  expect_lte(P_hat[1, 1], P_hat[2, 2])
  P_ind <- which(as.numeric(a$mdl_h0$theta_P_ind) == 1)
  pos <- matrix(P_ind, 2, 2)
  expect_identical(a$theta_fix[[pos[2, 2]]], 1)
  expect_identical(a$theta_fix[[pos[1, 2]]], 0)
  other <- setdiff(seq_along(a$theta_fix), c(pos[2, 2], pos[1, 2]))
  expect_identical(unname(a$theta_fix[other]), unname(as.numeric(a$mdl_h0$theta)[other]))
  expect_true(valid_pval(a$pval_lmc)); expect_true(valid_pval(a$pval_fix))
  expect_identical(a$pval, max(a$pval_lmc, a$pval_fix))
  # scope errors
  expect_error(aBfixLMCTest(y_b, 1, 2, 3, control = list(burnin = 100)), "burnin = 0")
  expect_error(aBfixLMCTest(y_ur, 1, 1, 2), "k0 >= 2")
})

test_that("aURdiffLMCTest: components equal standalone runs; p = 1 exercises the p = 0 path", {
  skip_on_cran()
  ctl <- list(N = 19, mc_seed = 76, mdl_h1_control = ctl2)
  a <- suppressWarnings(aURdiffLMCTest(y_ur, p = 1, k0 = 1, k1 = 2, control = ctl))
  l_lev <- suppressWarnings(LMCLRTest(y_ur, p = 1, k0 = 1, k1 = 2, control = ctl))
  l_dif <- suppressWarnings(LMCLRTest(diff(y_ur), p = 0, k0 = 1, k1 = 2, control = ctl))
  expect_identical(a$pval_levels, l_lev$pval)
  expect_identical(a$pval_diff, l_dif$pval)
  expect_identical(a$lmc_diff$mdl_h0$p, 0)   # the differenced null is p = 0 (Nmdl)
  expect_identical(a$pval, max(a$pval_levels, a$pval_diff))
  expect_error(aURdiffLMCTest(y_ur, 0, 1, 2), "p >= 1")
})

test_that("exogenous regressors flow through every augmented construction", {
  skip_on_cran()
  set.seed(37)
  Zx <- matrix(rnorm(400), 400, 1)
  y_urz <- as.matrix(cumsum(rnorm(400)) + 2 * Zx)
  ctlz <- list(N = 19, mc_seed = 88, mdl_h1_control = ctl2)
  # aURfix: component 1 bitwise-identical to LMCLRTest with the same Z; the
  # unit-root edit touches only the AR coefficient -- the exogenous
  # coefficient stays at its estimate in theta_fix
  a <- suppressWarnings(aURfixLMCTest(y_urz, p = 1, k0 = 1, k1 = 2, Z = Zx, control = ctlz))
  l <- suppressWarnings(LMCLRTest(y_urz, p = 1, k0 = 1, k1 = 2, Z = Zx, control = ctlz))
  expect_identical(a$pval_lmc, l$pval)
  expect_identical(a$LRN, l$LRN)
  x_pos <- which(as.numeric(a$mdl_h0$theta_x_ind) == 1)
  expect_identical(unname(a$theta_fix[x_pos]), unname(as.numeric(a$mdl_h0$theta)[x_pos]))
  expect_true(valid_pval(a$pval_fix))
  # aURdiff: components equal standalone runs; the differenced component keeps
  # Z in LEVELS, rows 2..T (never differenced)
  d <- suppressWarnings(aURdiffLMCTest(y_urz, p = 1, k0 = 1, k1 = 2, Z = Zx, control = ctlz))
  l_lev <- suppressWarnings(LMCLRTest(y_urz, p = 1, k0 = 1, k1 = 2, Z = Zx, control = ctlz))
  l_dif <- suppressWarnings(LMCLRTest(diff(y_urz), p = 0, k0 = 1, k1 = 2,
                                      Z = Zx[-1, , drop = FALSE], control = ctlz))
  expect_identical(d$pval_levels, l_lev$pval)
  expect_identical(d$pval_diff, l_dif$pval)
  # the variance floor's anchor on the differenced sample includes Z: the
  # alternative fit's floor equals the floor from the one-regime fit WITH Z,
  # and differs from the Z-free anchor whenever Z carries real signal
  y_d <- diff(y_urz); Z_d <- Zx[-1, , drop = FALSE]
  anchor_withZ <- Nmdl(y_d, Z = Z_d, control = list(getSE = FALSE, const = TRUE))
  anchor_noZ   <- Nmdl(y_d, control = list(getSE = FALSE, const = TRUE))
  expect_equal(unname(d$lmc_diff$mdl_h1$sigma_floor),
               0.01 * mean(diag(as.matrix(anchor_withZ$sigma))), tolerance = 1e-10)
  expect_gt(abs(0.01 * mean(diag(as.matrix(anchor_noZ$sigma))) -
                unname(d$lmc_diff$mdl_h1$sigma_floor)), 1e-8)
  # aBfix: boundary DGP with a Z effect; the projection leaves the exogenous
  # coefficient at its estimate
  set.seed(38)
  Zb <- matrix(rnorm(200), 200, 1)
  y_bz <- as.matrix(simuMSAR(list(n = 200, k = 2, mu = c(0, 2), sigma = c(1, 2), phi = 0.3,
                                  P = matrix(c(0.9, 0.1, 0, 1), 2, 2)), burnin = 0)$y + 1.5 * Zb)
  b <- suppressWarnings(aBfixLMCTest(y_bz, p = 1, k0 = 2, k1 = 3, Z = Zb,
        control = list(N = 19, mc_seed = 89, burnin = 0,
                       mdl_h0_control = ctl2, mdl_h1_control = ctl2)))
  xb_pos <- which(as.numeric(b$mdl_h0$theta_x_ind) == 1)
  expect_identical(unname(b$theta_fix[xb_pos]), unname(as.numeric(b$mdl_h0$theta)[xb_pos]))
  P_ind <- which(as.numeric(b$mdl_h0$theta_P_ind) == 1)
  expect_identical(b$theta_fix[[matrix(P_ind, 2, 2)[2, 2]]], 1)
  expect_true(valid_pval(b$pval_lmc)); expect_true(valid_pval(b$pval_fix))
})

test_that("aBfixLMCTest generalizes over transient regimes (k0 = 3)", {
  skip_on_cran()
  set.seed(39)
  y3 <- as.matrix(simuHMM(list(n = 150, k = 3, mu = as.matrix(c(-2, 0, 2)),
                               sigma = list(as.matrix(1), as.matrix(1.5), as.matrix(0.5)), q = 1,
                               P = matrix(c(0.8, 0.1, 0.1,
                                            0.15, 0.8, 0.05,
                                            0, 0, 1), 3, 3)), burnin = 0)$y)
  a <- suppressWarnings(aBfixLMCTest(y3, p = 0, k0 = 3, k1 = 4,
        control = list(N = 19, mc_seed = 93, burnin = 0,
                       mdl_h0_control = ctl2, mdl_h1_control = ctl2)))
  pos3 <- matrix(which(as.numeric(a$mdl_h0$theta_P_ind) == 1), 3, 3)
  # projection: the absorbing (third, canonical) column at exactly (0, 0, 1)
  expect_identical(a$theta_fix[[pos3[3, 3]]], 1)
  expect_true(all(a$theta_fix[pos3[1:2, 3]] == 0))
  other <- setdiff(seq_along(a$theta_fix), as.integer(pos3[, 3]))
  expect_identical(unname(a$theta_fix[other]), unname(as.numeric(a$mdl_h0$theta)[other]))
  expect_true(valid_pval(a$pval_lmc)); expect_true(valid_pval(a$pval_fix))
  expect_identical(a$pval, max(a$pval_lmc, a$pval_fix))
})

test_that("p = 0 testing path: k0 = 2 (hidden Markov) null runs through LMCLRTest", {
  skip_on_cran()
  # the regression test the Gate B work requested: the p = 0 path with a
  # k0 >= 2 null (differenced boundary-style data) must run end to end
  set.seed(33)
  y_hm <- as.matrix(simuHMM(list(n = 150, k = 2, mu = as.matrix(c(-1, 1.5)),
                                 sigma = list(as.matrix(2), as.matrix(1)), q = 1,
                                 P = matrix(c(0.9, 0.1, 0.15, 0.85), 2, 2)))$y)
  l <- suppressWarnings(LMCLRTest(y_hm, p = 0, k0 = 2, k1 = 3,
        control = list(N = 19, mc_seed = 78,
                       mdl_h0_control = ctl2, mdl_h1_control = ctl2)))
  expect_true(valid_pval(l$pval))
  expect_identical(l$mdl_h0$p, 0)
  expect_true(is.finite(l$mdl_h0$logLike))
})

test_that("aURdiffMMCTest: control_diff derivation strips the levels-specific box settings", {
  skip_on_cran()
  fit <- ARmdl(y_ur, p = 1, control = list(getSE = TRUE, const = TRUE))
  phi_pos <- which(as.numeric(fit$theta_phi_ind) == 1)
  binc <- rep(NA_real_, length(fit$theta)); binc[phi_pos] <- 1
  # phi_upp = 1.02 keeps the levels box out of the explosive region (documented
  # cost guidance); the vector eps still exercises the control_diff stripping
  ctl <- list(N = 19, mc_seed = 77, maxit = 4, silence = TRUE,
              stationary_constraint = FALSE, box_include = binc,
              phi_upp = 1.02,
              eps = rep(0.1, length(fit$theta)),
              mdl_h1_control = ctl2)
  a <- suppressWarnings(aURdiffMMCTest(y_ur, p = 1, k0 = 1, k1 = 2, control = ctl))
  expect_identical(a$mmc_diff$mdl_h0$p, 0)
  # derived control: box_include and vector eps gone, stationarity restored
  expect_null(a$control_diff$box_include)
  expect_true(isTRUE(a$control_diff$stationary_constraint))
  expect_true(length(a$control_diff$eps) == 1)
  # the levels component kept its own settings
  expect_false(isTRUE(a$control$stationary_constraint))
  expect_identical(a$pval, max(a$pval_levels, a$pval_diff))
  # the derivation itself, unit level
  cd <- MSTest:::.aug_diff_control(ctl)
  expect_null(cd$box_include); expect_null(cd$phi_upp)
  expect_true(isTRUE(cd$stationary_constraint))
  expect_null(cd$eps)   # vector eps dropped entirely (defaults reapply)
})

test_that("aURfixLMCTest parallel: component 1 matches parallel LMCLRTest, workers replay", {
  skip_on_cran()
  ctl <- list(N = 19, mc_seed = 87, workers = 2, mdl_h1_control = ctl2)
  a <- suppressWarnings(aURfixLMCTest(y_ur, p = 1, k0 = 1, k1 = 2, control = ctl))
  l <- suppressWarnings(LMCLRTest(y_ur, p = 1, k0 = 1, k1 = 2, control = ctl))
  expect_identical(a$pval_lmc, l$pval)
  expect_identical(a$LRN, l$LRN)
  expect_true(valid_pval(a$pval_fix))
  a2 <- suppressWarnings(aURfixLMCTest(y_ur, p = 1, k0 = 1, k1 = 2, control = ctl))
  expect_identical(a$pval_fix, a2$pval_fix)
})

test_that("aBfixLMCTest runs at p = 0 (hidden Markov path)", {
  skip_on_cran()
  set.seed(34)
  y_hb <- as.matrix(simuHMM(list(n = 150, k = 2, mu = as.matrix(c(0, 2)),
                                 sigma = list(as.matrix(1), as.matrix(2)), q = 1,
                                 P = matrix(c(0.9, 0.1, 0, 1), 2, 2)), burnin = 0)$y)
  a <- suppressWarnings(aBfixLMCTest(y_hb, p = 0, k0 = 2, k1 = 3,
        control = list(N = 19, mc_seed = 86, burnin = 0,
                       mdl_h0_control = ctl2, mdl_h1_control = ctl2)))
  P_ind <- which(as.numeric(a$mdl_h0$theta_P_ind) == 1)
  pos <- matrix(P_ind, 2, 2)
  expect_identical(a$theta_fix[[pos[2, 2]]], 1)
  expect_true(valid_pval(a$pval_lmc)); expect_true(valid_pval(a$pval_fix))
})

test_that("contained per-draw estimation failures do not leak error conditions", {
  skip_on_cran()
  # A strongly explosive parameter makes every simulated draw's estimation
  # fail; LR_samp_dist must contain those failures completely: the returned
  # vector is all-NaN (correct: the candidate is unsimulable) AND no error
  # condition reaches a calling handler upstack -- condition-recording
  # harnesses (testthat itself, logging wrappers) would otherwise record
  # hundreds of phantom errors per penalized candidate.
  set.seed(35)
  y_s <- as.matrix(arima.sim(list(ar = 0.5), n = 100, sd = 1))
  fit_s <- ARmdl(y_s, p = 1, control = list(getSE = FALSE, const = TRUE))
  theta_x <- as.numeric(fit_s$theta)
  theta_x[which(as.numeric(fit_s$theta_phi_ind) == 1)] <- 2
  mdl_x <- mdledit(fit_s, theta_x, 1, 1, 1, FALSE)
  set.seed(36)
  pe <- lapply(1:4, function(i) matrix(rnorm(200), 200, 1))
  n_seen <- 0L
  out <- withCallingHandlers(
    suppressWarnings(LR_samp_dist(mdl_x, 2, 2, 100, NULL,
                                  fit_s$control, list(use_diff_init = 1),
                                  pe, NULL)),
    error = function(e) n_seen <<- n_seen + 1L)
  expect_identical(n_seen, 0L)
  expect_true(all(!is.finite(out)))
})

test_that("print and summary methods for the augmented classes run", {
  skip_on_cran()
  ctl <- list(N = 19, mc_seed = 74, mdl_h1_control = ctl2)
  a <- suppressWarnings(aURfixLMCTest(y_ur, p = 1, k0 = 1, k1 = 2, control = ctl))
  expect_output(print(a), "Unit-Root-Augmented")
  expect_output(summary(a), "composition")
  d <- suppressWarnings(aURdiffLMCTest(y_ur, p = 1, k0 = 1, k1 = 2, control = ctl))
  expect_output(print(d), "differencing")
  expect_output(summary(d), "Component models")
})

test_that("aBfixLMCTest multivariate: projection edits only the absorbing column, covariances untouched", {
  skip_on_cran()
  set.seed(52)
  y2b <- as.matrix(simuMSVAR(list(n = 150, p = 1, q = 2,
                                  mu = rbind(c(0, 0), c(2, -1)),
                                  sigma = list(diag(2), diag(2) * 2),
                                  phi = diag(2) * 0.2, k = 2,
                                  P = matrix(c(0.9, 0.1, 0, 1), 2, 2)),
                             burnin = 0)$y)
  ctl_mv <- list(N = 19, mc_seed = 91, burnin = 0,
                 mdl_h0_control = list(use_diff_init = 1),
                 mdl_h1_control = list(use_diff_init = 1))
  a <- suppressWarnings(aBfixLMCTest(y2b, p = 1, k0 = 2, k1 = 3, control = ctl_mv))
  expect_identical(a$mdl_h0$q, 2L)
  pos2 <- matrix(which(as.numeric(a$mdl_h0$theta_P_ind) == 1), 2, 2)
  # canonical regime order: persistence ascending, absorbing last
  expect_lte(a$mdl_h0$theta[[pos2[1, 1]]], a$mdl_h0$theta[[pos2[2, 2]]])
  # projection: the absorbing (second, canonical) column at exactly (0, 1)
  expect_identical(a$theta_fix[[pos2[2, 2]]], 1)
  expect_identical(a$theta_fix[[pos2[1, 2]]], 0)
  # everything else -- means, vech(Sigma) including the covariances, phi, the
  # transient column -- is the fitted (canonicalized) theta, bitwise
  other <- setdiff(seq_along(a$theta_fix), as.integer(pos2[, 2]))
  expect_identical(unname(a$theta_fix[other]), unname(as.numeric(a$mdl_h0$theta)[other]))
  # the fitted covariances carried into the fixed parameter assemble to
  # positive definite matrices (estimates of an admissible model)
  sig_pos <- which(as.numeric(a$mdl_h0$theta_sig_ind) == 1)
  for (m in 1:2){
    Sig <- covar_unvech(a$theta_fix[sig_pos[((m - 1) * 3 + 1):(m * 3)]], 2)
    expect_gt(min(eigen(Sig, symmetric = TRUE, only.values = TRUE)$values), 0)
  }
  expect_true(valid_pval(a$pval_lmc)); expect_true(valid_pval(a$pval_fix))
  expect_identical(a$pval, max(a$pval_lmc, a$pval_fix))
  # deterministic repeat, both components
  a2 <- suppressWarnings(aBfixLMCTest(y2b, p = 1, k0 = 2, k1 = 3, control = ctl_mv))
  expect_identical(a$pval_lmc, a2$pval_lmc)
  expect_identical(a$pval_fix, a2$pval_fix)
})

test_that("aURdiff pair multivariate: per-series difference, components equal standalone runs", {
  skip_on_cran()
  set.seed(53)
  y2u <- apply(matrix(rnorm(300), 150, 2), 2, cumsum)
  ctl_mv <- list(N = 19, mc_seed = 95,
                 mdl_h0_control = list(use_diff_init = 1),
                 mdl_h1_control = list(use_diff_init = 1))
  d <- suppressWarnings(aURdiffLMCTest(y2u, p = 1, k0 = 1, k1 = 2, control = ctl_mv))
  l_lev <- suppressWarnings(LMCLRTest(y2u, 1, 1, 2, control = ctl_mv))
  l_dif <- suppressWarnings(LMCLRTest(diff(y2u), 0, 1, 2, control = ctl_mv))
  expect_identical(d$pval_levels, l_lev$pval)
  expect_identical(d$pval_diff, l_dif$pval)
  expect_identical(d$pval, max(d$pval_levels, d$pval_diff))
})

test_that("aURfixLMCTest stays univariate: the fixed point needs a rank choice at q > 1", {
  set.seed(54)
  y2x <- matrix(rnorm(120), 60, 2)
  expect_error(aURfixLMCTest(y2x, 1, 1, 2), "univariate models only")
})

test_that("aURfixLMCTest warns that boundary_R is ignored", {
  skip_on_cran()
  ctl <- list(N = 19, mc_seed = 74, mdl_h1_control = ctl2,
              boundary_R = list(mu_mult = 5))
  expect_warning(aURfixLMCTest(y_ur, p = 1, k0 = 1, k1 = 2, control = ctl),
                 "boundary_R.*not used")
})

test_that("aBfix multivariate replays the stream with the (T+burnin) x q predraw layout", {
  skip_on_cran()
  # The q > 1 analogue of the univariate shared-randomness test: an
  # independent reimplementation of the engine's choreography, including the
  # multivariate predraw layout (one (T+burnin) x q normal matrix per buffer
  # slot, column-major fill), must reproduce both component p-values AND the
  # post-call RNG state. A byrow or per-column predraw yields a different
  # stream layout and fails the pin even where the p-values coincide.
  set.seed(52)
  y2b <- as.matrix(simuMSVAR(list(n = 150, p = 1, q = 2,
                                  mu = rbind(c(0, 0), c(2, -1)),
                                  sigma = list(diag(2), diag(2) * 2),
                                  phi = diag(2) * 0.2, k = 2,
                                  P = matrix(c(0.9, 0.1, 0, 1), 2, 2)),
                             burnin = 0)$y)
  ctl_mv <- list(N = 19, mc_seed = 91, burnin = 0,
                 mdl_h0_control = list(use_diff_init = 1),
                 mdl_h1_control = list(use_diff_init = 1))
  a <- suppressWarnings(aBfixLMCTest(y2b, p = 1, k0 = 2, k1 = 3, control = ctl_mv))
  seed_after_pkg <- get(".Random.seed", envir = globalenv())
  set.seed(91)
  mdl_h0 <- estimMdl(y2b, 1, 2, 2, NULL, list(use_diff_init = 1))
  mdl_h0 <- canonicalize_regimes(mdl_h0, by = "persistence")
  con_h1 <- MSTest:::.inherit_ms_flags(list(use_diff_init = 1), mdl_h0)
  con_h1$init_theta_extra <- warmstart_theta(mdl_h0, 3, msmu = con_h1[["msmu"]],
                                             msvar = con_h1[["msvar"]])
  mdl_h1 <- estimMdl(y2b, 1, 2, 3, NULL, con_h1)
  logL0 <- mdl_h0$logLike
  mdl_lin <- estimMdl(y2b, 1, 2, 1, NULL, list(getSE = FALSE))
  if (is.finite(mdl_lin$logLike) && mdl_lin$logLike > logL0) logL0 <- mdl_lin$logLike
  LRT_0 <- -2 * (logL0 - mdl_h1$logLike)
  if (LRT_0 < 0) LRT_0 <- 0
  theta_fix <- MSTest:::.theta_boundary_projection(mdl_h0$theta,
                                                   MSTest:::.mmc_boundary_blocks(mdl_h0))
  set.seed(91)
  Nb <- ceiling(19 * 1.5) + 10; Teps <- mdl_h0$n + 0
  pe  <- lapply(seq_len(Nb), function(i) matrix(rnorm(Teps * 2), Teps, 2))
  psr <- matrix(runif(Teps * Nb), Teps, Nb)
  h0c <- mdl_h0$control
  h1c <- mdl_h1$control; h1c$init_theta_extra <- NULL; h1c$init_method <- "warmstart"
  S <- get(".Random.seed", envir = globalenv())
  LRN1 <- suppressWarnings(LR_samp_dist(mdl_h0, 3, 19, 0, NULL, h0c, h1c, pe, psr))
  p1m <- MCpval(unname(LRT_0), LRN1[is.finite(LRN1)], "geq")
  assign(".Random.seed", S, envir = globalenv())
  LRN2 <- suppressWarnings(LR_samp_dist(mdledit(mdl_h0, as.numeric(theta_fix), 1, 2, 2, FALSE),
                                        3, 19, 0, NULL, h0c, h1c, pe, psr))
  p2m <- MCpval(unname(LRT_0), LRN2[is.finite(LRN2)], "geq")
  seed_after_manual <- get(".Random.seed", envir = globalenv())
  expect_identical(p1m, a$pval_lmc)
  expect_identical(p2m, a$pval_fix)
  expect_identical(seed_after_pkg, seed_after_manual)
})
