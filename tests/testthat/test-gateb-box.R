# Tests for the unit-root and boundary search-box primitives:
# box_include, box_boundary ("union"/"face"), boundary_R, and the search
# diagnostics (realized box + candidate recorder).

valid_pval <- function(p) is.finite(p) && p >= 0 && p <= 1

# Shared datasets -------------------------------------------------------------
set.seed(31)
y_ur <- as.matrix(cumsum(rnorm(400)))          # exact unit root, T = 400
set.seed(32)
y_b <- as.matrix(simuMSAR(list(n = 200, k = 2, mu = c(0, 2), sigma = c(1, 2),
                               phi = 0.3, P = matrix(c(0.9, 0.1, 0, 1), 2, 2)),
                          burnin = 0)$y)       # boundary DGP: p22 = 1, transient start
ctl2 <- list(use_diff_init = 2)

test_that("box_include stretches only where needed and validates loudly", {
  skip_on_cran()
  fit <- ARmdl(y_ur, p = 1, control = list(getSE = TRUE, const = TRUE))
  phi_pos <- which(as.numeric(fit$theta_phi_ind) == 1)
  binc <- rep(NA_real_, length(fit$theta)); binc[phi_pos] <- 1
  con <- list(eps = 0.1, CI_union = TRUE, P_low = 0, P_upp = 1,
              phi_low = NULL, phi_upp = NULL, variance_constraint = 0.01,
              type = "pso", box_include = binc)
  b1 <- MMC_bounds(fit, con)
  # the box must contain the unit root at the phi coordinate, exactly
  expect_gte(b1$theta_upp[phi_pos], 1)
  expect_lte(b1$theta_low[phi_pos], 1)
  # a coordinate already covered is untouched: include the fitted value itself
  binc2 <- rep(NA_real_, length(fit$theta)); binc2[phi_pos] <- as.numeric(fit$theta)[phi_pos]
  con2 <- con; con2$box_include <- binc2
  b_ref <- MMC_bounds(fit, utils::modifyList(con, list(box_include = NULL)))
  b2 <- MMC_bounds(fit, con2)
  expect_identical(b2$theta_low, b_ref$theta_low)
  expect_identical(b2$theta_upp, b_ref$theta_upp)
  # validation: wrong length, non-finite entry, contradiction with a clamp
  expect_error(MMC_bounds(fit, utils::modifyList(con, list(box_include = c(1, 2)))),
               "box_include")
  binc_bad <- binc; binc_bad[phi_pos] <- Inf
  expect_error(MMC_bounds(fit, utils::modifyList(con, list(box_include = binc_bad))),
               "finite")
  expect_error(MMC_bounds(fit, utils::modifyList(con, list(phi_upp = 0.9))),
               "contradictory")
})

test_that("box_include is a two-sided hull: covers the point when the ordinary interval lies entirely above it", {
  skip_on_cran()
  # With a fitted phi above 1 + w, the ordinary (eps/CI) interval [lo, hi]
  # excludes 1 from BELOW. One-sided "raise the upper endpoint to reach v"
  # semantics would leave 1 outside the box, defeating the deterministic-
  # coverage requirement of the unit-root search set; the stretched interval
  # must be the convex hull [min(lo, v), max(hi, v)]. Fitted phi > 1 is a
  # real event under a true unit root (~4% of replications at T = 500), and
  # phi_hat - w > 1 is exactly the configuration where the two forms differ.
  fit <- ARmdl(y_ur, p = 1, control = list(getSE = TRUE, const = TRUE))
  phi_pos <- which(as.numeric(fit$theta_phi_ind) == 1)
  fit$theta[phi_pos] <- 1.15                   # mock an explosive fitted value
  binc <- rep(NA_real_, length(fit$theta)); binc[phi_pos] <- 1
  con <- list(eps = 0.02, CI_union = FALSE, P_low = 0, P_upp = 1,
              phi_low = NULL, phi_upp = NULL, variance_constraint = 0.01,
              type = "pso", box_include = binc)
  b <- MMC_bounds(fit, con)
  # the discriminating premise: the ordinary interval sits entirely above 1
  expect_gt(unname(fit$theta[phi_pos]) - 0.02, 1)
  # hull: the lower endpoint is pulled down to exactly the include point...
  expect_identical(b$theta_low[phi_pos], 1)
  # ...and the whole ordinary interval [phi_hat - eps, phi_hat + eps] survives
  expect_equal(b$theta_upp[phi_pos], unname(fit$theta[phi_pos]) + 0.02, tolerance = 1e-12)
})

test_that("unit-root MMC end-to-end: box_include + non-binding phi_upp, no E7, diagnostics", {
  skip_on_cran()
  fit <- ARmdl(y_ur, p = 1, control = list(getSE = TRUE, const = TRUE))
  phi_pos <- which(as.numeric(fit$theta_phi_ind) == 1)
  binc <- rep(NA_real_, length(fit$theta)); binc[phi_pos] <- 1
  # phi_upp = 1.02: the box must contain the unit root but should not reach
  # far into the explosive region, where every draw of a candidate fails
  # estimation and the candidate is only penalized after the full retry
  # budget -- the documented cost guidance for this configuration.
  ur_ctl <- list(N = 19, mc_seed = 71, maxit = 6, silence = TRUE,
                 stationary_constraint = FALSE, box_include = binc,
                 phi_upp = 1.02, mdl_h1_control = ctl2)
  m <- MMCLRTest(y_ur, p = 1, k0 = 1, k1 = 2, control = ur_ctl)
  expect_gte(m$theta_upp[phi_pos], 1)
  expect_lte(m$theta_upp[phi_pos], 1.02)
  expect_true(valid_pval(m$pval))
  expect_gte(m$pval, m$pval_0)
  # harvest: realized box named and consistent; recorder exact and populated,
  # and every evaluated candidate lies inside the realized box
  expect_identical(names(m$theta_low), names(m$mdl_h0$theta))
  expect_true(all(m$theta_low <= m$theta_upp))
  expect_gt(m$optim_diagnostics$feval, 0)
  er <- m$optim_diagnostics$evaluated_range
  expect_true(is.matrix(er))
  expect_true(all(er["min", ] <= er["max", ]))
  expect_true(all(er["min", ] >= m$theta_low[colnames(er)] - 1e-12))
  expect_true(all(er["max", ] <= m$theta_upp[colnames(er)] + 1e-12))
  # determinism: identical repeat
  m2 <- MMCLRTest(y_ur, p = 1, k0 = 1, k1 = 2, control = ur_ctl)
  expect_identical(m$pval, m2$pval)
  expect_identical(m$theta_h0, m2$theta_h0)
})

test_that("boundary 'union' box: canonical labeling, exact stretch, hull; scope errors", {
  skip_on_cran()
  set.seed(41)
  m <- MMCLRTest(y_b, p = 1, k0 = 2, k1 = 3,
                 control = list(N = 19, mc_seed = 72, maxit = 6, silence = TRUE, burnin = 0,
                                box_boundary = "union",
                                mdl_h0_control = ctl2, mdl_h1_control = ctl2))
  expect_identical(m$boundary_mode, "union")
  expect_identical(m$absorbing_regime, 2)
  # canonical: absorbing regime last (larger persistence in slot 2)
  P_hat <- matrix(as.numeric(m$mdl_h0$theta)[as.numeric(m$mdl_h0$theta_P_ind) == 1], 2, 2)
  expect_lte(P_hat[1, 1], P_hat[2, 2])
  # exact stretch on the absorbing column, in full coordinates
  P_ind <- which(as.numeric(m$mdl_h0$theta_P_ind) == 1)
  pos <- matrix(P_ind, 2, 2)
  expect_identical(m$theta_upp[[pos[2, 2]]], 1)
  expect_identical(m$theta_low[[pos[1, 2]]], 0)
  # rough block covers the moment-anchored ranges (defaults: 5, 4, [.05, .95])
  yv <- as.numeric(m$mdl_h0$y); ybar <- mean(yv); s_y <- stats::sd(yv)
  mu_pos  <- which(as.numeric(m$mdl_h0$theta_mu_ind) == 1)
  sig_pos <- which(as.numeric(m$mdl_h0$theta_sig_ind) == 1)
  expect_lte(m$theta_low[[mu_pos[1]]], ybar - 5 * s_y)
  expect_gte(m$theta_upp[[mu_pos[1]]], ybar + 5 * s_y)
  expect_gte(m$theta_upp[[sig_pos[1]]], 4 * s_y^2)
  expect_lte(m$theta_low[[pos[1, 1]]], 0.05)
  expect_gte(m$theta_upp[[pos[1, 1]]], 0.95)
  # clean block NOT stretched to the mu range (absorbing regime's mean)
  expect_gt(m$theta_low[[mu_pos[2]]], ybar - 5 * s_y)
  expect_true(valid_pval(m$pval)); expect_gte(m$pval, m$pval_0)
  # scope errors
  expect_error(MMCLRTest(y_ur, 1, 1, 2, control = list(N = 19, burnin = 0, box_boundary = "union")),
               "k0 >= 2")
  expect_error(MMCLRTest(y_b, 1, 2, 3, control = list(N = 19, box_boundary = "union")),
               "burnin = 0")
  expect_error(MMCLRTest(y_b, 1, 2, 3, control = list(N = 19, burnin = 0, box_boundary = "union",
                                                      P_upp = 0.99)),
               "P_low = 0 and P_upp = 1")
  expect_error(MMCLRTest(y_b, 1, 2, 3, control = list(N = 19, burnin = 0, box_boundary = "union",
                                                      boundary_R = list(p_range = c(0.5, 0.4)))),
               "p_range")
  expect_error(MMCLRTest(y_b, 1, 2, 3, control = list(N = 19, burnin = 0, box_boundary = "face",
                                                      box_include = rep(NA_real_, 9))),
               "not meaningful")
})

test_that("boundary 'union' stretch is load-bearing: reaches the face when the ordinary box does not", {
  skip_on_cran()
  # Interior-persistence data (p22_hat well below 1) with a deliberately tiny
  # ordinary box (eps = 0.01, no CI union): WITHOUT the stretch the box cannot
  # reach the face, so these assertions discriminate a dropped stretch --
  # unlike an absorbing DGP, where p22_hat ~ 1 and the P clamp alone already
  # lands the ordinary box on 1, making the same assertions vacuous there.
  set.seed(43)
  y_int <- as.matrix(simuMSAR(list(n = 200, k = 2, mu = c(-1, 1.5), sigma = c(1, 2),
                                   phi = 0.3, P = matrix(c(0.9, 0.1, 0.15, 0.85), 2, 2)),
                              burnin = 100)$y)
  fit <- MSARmdl(y_int, p = 1, k = 2, control = list(getSE = TRUE, use_diff_init = 2))
  fit_c <- canonicalize_regimes(fit, by = "persistence")
  P_ind <- which(as.numeric(fit_c$theta_P_ind) == 1)
  pos <- matrix(P_ind, 2, 2)
  con <- list(eps = 0.01, CI_union = FALSE, P_low = 0, P_upp = 1,
              phi_low = NULL, phi_upp = NULL, variance_constraint = 0.01,
              type = "pso", box_include = NULL,
              box_boundary = "union", boundary_R = NULL)
  b <- MMC_bounds(fit_c, con)
  b0 <- MMC_bounds(fit_c, utils::modifyList(con, list(box_boundary = NULL)))
  # the ordinary box does NOT reach the face here (this is what makes the
  # stretched assertions below non-vacuous)
  expect_lt(b0$theta_upp[pos[2, 2]], 1)
  expect_gt(b0$theta_low[pos[1, 2]], 0)
  # the stretched box reaches it exactly
  expect_identical(b$theta_upp[pos[2, 2]], 1)
  expect_identical(b$theta_low[pos[1, 2]], 0)
  # and the rough block covers R while the clean block stays ordinary
  yv <- as.numeric(fit_c$y); ybar <- mean(yv); s_y <- stats::sd(yv)
  mu_pos <- which(as.numeric(fit_c$theta_mu_ind) == 1)
  expect_lte(b$theta_low[mu_pos[1]], ybar - 5 * s_y)
  expect_gte(b$theta_upp[mu_pos[1]], ybar + 5 * s_y)
  expect_identical(b$theta_low[mu_pos[2]], b0$theta_low[mu_pos[2]])
  expect_identical(b$theta_upp[mu_pos[2]], b0$theta_upp[mu_pos[2]])
})

test_that("boundary 'union' with a common (non-switching) parameter keeps it clean", {
  skip_on_cran()
  set.seed(42)
  y_bv <- as.matrix(simuMSAR(list(n = 200, k = 2, mu = c(1, 1), sigma = c(1, 4),
                                  phi = 0.3, P = matrix(c(0.9, 0.1, 0, 1), 2, 2)),
                             burnin = 0)$y)
  m <- MMCLRTest(y_bv, p = 1, k0 = 2, k1 = 3,
                 control = list(N = 19, mc_seed = 79, maxit = 4, silence = TRUE, burnin = 0,
                                box_boundary = "union",
                                mdl_h0_control = list(use_diff_init = 2, msmu = FALSE),
                                mdl_h1_control = ctl2))
  # msmu = FALSE: the single common mean is identified by the whole sample ->
  # clean -> NOT stretched to the moment-anchored range
  yv <- as.numeric(m$mdl_h0$y); ybar <- mean(yv); s_y <- stats::sd(yv)
  mu_pos <- which(as.numeric(m$mdl_h0$theta_mu_ind) == 1)
  expect_length(mu_pos, 1)
  expect_gt(m$theta_low[[mu_pos]], ybar - 5 * s_y)
  expect_lt(m$theta_upp[[mu_pos]], ybar + 5 * s_y)
  # the switching variances' transient entry IS stretched
  sig_pos <- which(as.numeric(m$mdl_h0$theta_sig_ind) == 1)
  expect_gte(m$theta_upp[[sig_pos[1]]], 4 * s_y^2)
  expect_true(valid_pval(m$pval))
})

test_that("boundary 'face' box: pins, search dimension, projection, determinism", {
  skip_on_cran()
  m <- MMCLRTest(y_b, p = 1, k0 = 2, k1 = 3,
                 control = list(N = 19, mc_seed = 73, maxit = 6, silence = TRUE, burnin = 0,
                                box_boundary = "face",
                                mdl_h0_control = ctl2, mdl_h1_control = ctl2))
  expect_identical(m$boundary_mode, "face")
  P_ind <- which(as.numeric(m$mdl_h0$theta_P_ind) == 1)
  pos <- matrix(P_ind, 2, 2)
  # projection: absorbing column exactly (0, 1); everything else at the estimate
  expect_identical(m$theta_face_proj[[pos[2, 2]]], 1)
  expect_identical(m$theta_face_proj[[pos[1, 2]]], 0)
  other <- setdiff(seq_along(m$theta_face_proj), c(pos[2, 2], pos[1, 2]))
  expect_identical(unname(m$theta_face_proj[other]), unname(as.numeric(m$mdl_h0$theta)[other]))
  # pinned coordinates are zero-width in the realized box
  expect_identical(m$theta_low[[pos[2, 2]]], 1); expect_identical(m$theta_upp[[pos[2, 2]]], 1)
  phi_pos <- which(as.numeric(m$mdl_h0$theta_phi_ind) == 1)
  expect_identical(m$theta_low[[phi_pos]], m$theta_upp[[phi_pos]])
  # the search runs over the rough block only: msmu+msvar k0=2 -> 3 dimensions
  expect_identical(ncol(m$optim_diagnostics$evaluated_range), 3L)
  expect_true(valid_pval(m$pval)); expect_gte(m$pval, m$pval_proj)
  # the returned theta corresponds to the reported p-value (projection or
  # in-face point); its P column sums to 1 exactly either way
  P_ret <- matrix(as.numeric(m$theta_h0)[as.numeric(m$mdl_h0$theta_P_ind) == 1], 2, 2)
  expect_true(all(abs(colSums(P_ret) - 1) < 1e-8))
  m2 <- MMCLRTest(y_b, p = 1, k0 = 2, k1 = 3,
                  control = list(N = 19, mc_seed = 73, maxit = 6, silence = TRUE, burnin = 0,
                                 box_boundary = "face",
                                 mdl_h0_control = ctl2, mdl_h1_control = ctl2))
  expect_identical(m$pval, m2$pval)
  expect_identical(m$pval_proj, m2$pval_proj)
  expect_identical(m$theta_h0, m2$theta_h0)
})

test_that("reassembly with pinned coordinates is exact and NA-derived columns are skipped", {
  # synthetic k0 = 3 layout: theta = (mu1..3, sig1..3, vec(P) 9) -> npar = 15
  k0 <- 3; npar <- 15
  P_ind <- c(rep(0, 6), rep(1, 9))
  theta_P_pos <- which(P_ind == 1)
  pos_mat <- matrix(theta_P_pos, k0, k0)
  # pin the absorbing (3rd) column at (0, 0, 1) plus a clean coordinate (mu3)
  pinned_idx <- c(3L, pos_mat[, 3])
  pinned_vals <- c(0.5, 0, 0, 1)
  theta_0 <- c(0.1, 0.2, 0.5, 1, 2, 3, 0.8, 0.15, 0.05, 0.1, 0.7, 0.2, 0, 0, 1)
  rp <- MSTest:::.mmc_derived_P_index(theta_0, P_ind, k0, fixed_full_idx = pinned_idx)
  expect_true(is.na(rp$derived_full_idx[3]))          # fully pinned column skipped
  expect_false(anyNA(rp$derived_full_idx[1:2]))
  derived <- rp$derived_full_idx[!is.na(rp$derived_full_idx)]
  reduced_idx <- setdiff(seq_len(npar), c(derived, pinned_idx))
  x <- theta_0[reduced_idx]
  th <- MSTest:::.mmc_reassemble_theta(x, reduced_idx, rp$derived_full_idx, rp$pos_mat, npar,
                                       pinned_idx, pinned_vals)
  # pins written; searched entries written; every P column sums to exactly 1
  expect_identical(th[3], 0.5)
  expect_identical(th[pos_mat[3, 3]], 1)
  P_re <- matrix(th[theta_P_pos], k0, k0)
  # a + (1 - a) can differ from 1 by an ulp in floating point; the simulators'
  # own gate is 1e-8, so assert well inside that
  expect_true(all(abs(colSums(P_re) - 1) < 1e-14))
  # default (no pins) behavior unchanged: bitwise round-trip of theta_0
  rp0 <- MSTest:::.mmc_derived_P_index(theta_0, P_ind, k0)
  red0 <- setdiff(seq_len(npar), rp0$derived_full_idx)
  th0 <- MSTest:::.mmc_reassemble_theta(theta_0[red0], red0, rp0$derived_full_idx, rp0$pos_mat, npar)
  expect_equal(th0, theta_0, tolerance = 1e-15)
})

test_that("three-regime boundary: block split, per-segment structure, union box, absorbing index", {
  skip_on_cran()
  # The rough/clean split and the union stretch generalize over the number of
  # transient regimes: at k0 = 3 the first TWO regimes' switching parameters
  # and the first two P columns are rough, the absorbing regime is the third,
  # and the family's initial state is canonical state 1 of the fitted null.
  set.seed(44)
  y3 <- as.matrix(simuHMM(list(n = 200, k = 3, mu = as.matrix(c(-2, 0, 2)),
                               sigma = list(as.matrix(1), as.matrix(1.5), as.matrix(0.5)), q = 1,
                               P = matrix(c(0.8, 0.1, 0.1,
                                            0.15, 0.8, 0.05,
                                            0, 0, 1), 3, 3)), burnin = 0)$y)
  fit3 <- HMmdl(y3, k = 3, control = list(getSE = TRUE, use_diff_init = 2))
  fit3c <- canonicalize_regimes(fit3, by = "persistence")
  bb3 <- MSTest:::.mmc_boundary_blocks(fit3c)
  mu_pos  <- which(as.numeric(fit3c$theta_mu_ind) == 1)
  sig_pos <- which(as.numeric(fit3c$theta_sig_ind) == 1)
  pos3 <- matrix(which(as.numeric(fit3c$theta_P_ind) == 1), 3, 3)
  expect_identical(bb3$rough_mu, mu_pos[1:2])
  expect_identical(bb3$rough_sig, sig_pos[1:2])
  expect_true(all(c(mu_pos[3], sig_pos[3]) %in% bb3$clean_idx))
  expect_identical(sort(bb3$trans_P_idx), sort(as.integer(pos3[, 1:2])))
  expect_identical(bb3$abs_diag_idx, pos3[3, 3])
  expect_identical(bb3$abs_off_idx, pos3[1:2, 3])
  expect_identical(bb3$n_transient_seg, 2)
  # the per-segment off-diagonal map at q = 3, k0 = 3, on a mocked layout:
  # each off-diagonal must couple the two variances of its OWN segment
  qm <- 3L; Nsig <- 6L
  mock <- list(k = 3L, q = qm, p = 0L, msmu = TRUE, msvar = TRUE,
               theta_mu_ind = c(rep(1, 9), rep(0, 18 + 9)),
               theta_sig_ind = c(rep(0, 9), rep(1, 18), rep(0, 9)),
               theta_var_ind = c(rep(0, 9), rep(c(1, 0, 0, 1, 0, 1), 3), rep(0, 9)),
               theta_P_ind = c(rep(0, 27), rep(1, 9)),
               theta_phi_ind = NULL, theta_x_ind = NULL)
  bbm <- MSTest:::.mmc_boundary_blocks(mock)
  expect_identical(bbm$n_transient_seg, 2)
  seg1 <- bbm$rough_sig[1:6]; seg2 <- bbm$rough_sig[7:12]
  expect_true(all(bbm$off_diag_i[1:3] %in% seg1) && all(bbm$off_diag_j[1:3] %in% seg1))
  expect_true(all(bbm$off_diag_i[4:6] %in% seg2) && all(bbm$off_diag_j[4:6] %in% seg2))
  # vech order (1,1),(2,1),(3,1),(2,2),(3,2),(3,3): off entries couple
  # (2,1),(3,1),(3,2) -> diagonals at within-segment positions (4,1),(6,1),(6,4)
  expect_identical(bbm$off_diag_i, c(seg1[c(4, 6, 6)], seg2[c(4, 6, 6)]))
  expect_identical(bbm$off_diag_j, c(seg1[c(1, 1, 4)], seg2[c(1, 1, 4)]))
  # union end-to-end at k0 = 3 (early stop keeps it cheap): both transient
  # columns stretched, absorbing column reaches the face, index reported
  m3 <- MMCLRTest(y3, p = 0, k0 = 3, k1 = 4,
                  control = list(N = 19, mc_seed = 81, maxit = 3, silence = TRUE, burnin = 0,
                                 box_boundary = "union", threshold_stop = 1/(19 + 1),
                                 mdl_h0_control = ctl2, mdl_h1_control = ctl2))
  expect_identical(m3$absorbing_regime, 3)
  expect_identical(m3$theta_upp[[pos3[3, 3]]], 1)
  expect_true(all(m3$theta_low[pos3[1:2, 3]] == 0))
  expect_true(all(m3$theta_low[as.integer(pos3[, 1:2])] <= 0.05))
  expect_true(all(m3$theta_upp[as.integer(pos3[, 1:2])] >= 0.95))
  expect_true(valid_pval(m3$pval))
})

test_that("boundary modes: GA face, moved boundary_R, msmu = FALSE face, compositions", {
  skip_on_cran()
  # GA + face: the S4 @solution remap with pinned coordinates
  m_ga <- MMCLRTest(y_b, p = 1, k0 = 2, k1 = 3,
                    control = list(N = 19, mc_seed = 82, maxit = 1, type = "GA", silence = TRUE,
                                   burnin = 0, box_boundary = "face",
                                   optim_control = list(popSize = 4),
                                   mdl_h0_control = ctl2, mdl_h1_control = ctl2))
  P_ind <- which(as.numeric(m_ga$mdl_h0$theta_P_ind) == 1)
  pos <- matrix(P_ind, 2, 2)
  P_ret <- matrix(as.numeric(m_ga$theta_h0)[P_ind], 2, 2)
  expect_true(all(abs(colSums(P_ret) - 1) < 1e-8))
  expect_true(valid_pval(m_ga$pval))
  # a non-default boundary_R actually moves the realized box
  m_bR <- MMCLRTest(y_b, p = 1, k0 = 2, k1 = 3,
                    control = list(N = 19, mc_seed = 83, maxit = 2, silence = TRUE, burnin = 0,
                                   box_boundary = "union", threshold_stop = 1/(19 + 1),
                                   boundary_R = list(mu_mult = 3),
                                   mdl_h0_control = ctl2, mdl_h1_control = ctl2))
  yv <- as.numeric(m_bR$mdl_h0$y); ybar <- mean(yv); s_y <- stats::sd(yv)
  mu_pos <- which(as.numeric(m_bR$mdl_h0$theta_mu_ind) == 1)
  expect_lte(m_bR$theta_low[[mu_pos[1]]], ybar - 3 * s_y)
  expect_gt(m_bR$theta_low[[mu_pos[1]]], ybar - 5 * s_y)
  # face with msmu = FALSE: the rough block is variance-only (+ transient P)
  set.seed(45)
  y_bv <- as.matrix(simuMSAR(list(n = 200, k = 2, mu = c(1, 1), sigma = c(1, 4),
                                  phi = 0.3, P = matrix(c(0.9, 0.1, 0, 1), 2, 2)),
                             burnin = 0)$y)
  m_f0 <- MMCLRTest(y_bv, p = 1, k0 = 2, k1 = 3,
                    control = list(N = 19, mc_seed = 84, maxit = 2, silence = TRUE, burnin = 0,
                                   box_boundary = "face",
                                   mdl_h0_control = list(use_diff_init = 2, msmu = FALSE),
                                   mdl_h1_control = ctl2))
  expect_identical(ncol(m_f0$optim_diagnostics$evaluated_range), 2L)  # sig_1, p11
  expect_true(valid_pval(m_f0$pval))
  # union + box_include compose (both widen; no interaction hazard)
  fitb <- MSARmdl(y_b, p = 1, k = 2, control = list(getSE = TRUE, use_diff_init = 2))
  fitb_c <- canonicalize_regimes(fitb, by = "persistence")
  binc2 <- rep(NA_real_, length(fitb_c$theta))
  mu_pos_b <- which(as.numeric(fitb_c$theta_mu_ind) == 1)
  binc2[mu_pos_b[2]] <- 25   # stretch a CLEAN coordinate far beyond CI and R
  con_u <- list(eps = 0.1, CI_union = TRUE, P_low = 0, P_upp = 1,
                phi_low = NULL, phi_upp = NULL, variance_constraint = 0.01,
                type = "pso", box_include = binc2,
                box_boundary = "union", boundary_R = NULL)
  b_u <- MMC_bounds(fitb_c, con_u)
  expect_gte(b_u$theta_upp[mu_pos_b[2]], 25)
  expect_identical(b_u$theta_upp[matrix(which(as.numeric(fitb_c$theta_P_ind) == 1), 2, 2)[2, 2]], 1)
  # box_include on a P coordinate at k0 = 2: per-coordinate containment holds
  binc3 <- rep(NA_real_, length(fitb_c$theta))
  p11_pos <- matrix(which(as.numeric(fitb_c$theta_P_ind) == 1), 2, 2)[1, 1]
  binc3[p11_pos] <- 0.5
  con_p <- list(eps = 0.01, CI_union = FALSE, P_low = 0, P_upp = 1,
                phi_low = NULL, phi_upp = NULL, variance_constraint = 0.01,
                type = "pso", box_include = binc3)
  b_p <- MMC_bounds(fitb_c, con_p)
  expect_lte(b_p$theta_low[p11_pos], 0.5)
  expect_gte(b_p$theta_upp[p11_pos], 0.5)
  # exported-helper misuse guard: union bounds on an UNcanonicalized fit whose
  # absorbing regime is not last must error, not silently stretch the wrong column
  fit_swap <- fitb_c
  perm_theta <- fit_swap$theta
  # swap the two regimes' P columns/rows in theta to fake a raw fit with the
  # absorbing regime FIRST
  P_hat <- matrix(as.numeric(fitb_c$theta)[which(as.numeric(fitb_c$theta_P_ind) == 1)], 2, 2)
  if (P_hat[1, 1] < P_hat[2, 2]){
    P_sw <- P_hat[2:1, 2:1]
    perm_theta[which(as.numeric(fitb_c$theta_P_ind) == 1)] <- as.numeric(P_sw)
    fit_swap$theta <- perm_theta
    expect_error(MMC_bounds(fit_swap, con_u), "canonical")
  }
})

test_that("face mode early stop returns the projection as theta_h0", {
  skip_on_cran()
  m <- MMCLRTest(y_b, p = 1, k0 = 2, k1 = 3,
                 control = list(N = 19, mc_seed = 85, maxit = 4, silence = TRUE, burnin = 0,
                                box_boundary = "face", threshold_stop = 1/(19 + 1),
                                mdl_h0_control = ctl2, mdl_h1_control = ctl2))
  expect_identical(m$optim_diagnostics$type, "early_stop")
  expect_identical(unname(m$theta_h0), unname(m$theta_face_proj))
  expect_identical(m$pval, m$pval_proj)
})

test_that("multivariate boundary box: per-series ranges, Cauchy-Schwarz off-diagonals, PSD gate", {
  skip_on_cran()
  set.seed(46)
  y2 <- as.matrix(simuMSVAR(list(n = 200, p = 1, q = 2, mu = rbind(c(0, 0), c(2, -1)),
                                 sigma = list(diag(2), diag(2) * 2), phi = diag(2) * 0.2,
                                 k = 2, P = matrix(c(0.9, 0.1, 0, 1), 2, 2)), burnin = 0)$y)
  fit2 <- MSVARmdl(y2, p = 1, k = 2, control = list(getSE = TRUE, use_diff_init = 2))
  fit2c <- canonicalize_regimes(fit2, by = "persistence")
  con_u <- list(eps = 0.1, CI_union = TRUE, P_low = 0, P_upp = 1,
                phi_low = NULL, phi_upp = NULL, variance_constraint = 0.01,
                type = "pso", box_include = NULL,
                box_boundary = "union", boundary_R = NULL)
  b <- MMC_bounds(fit2c, con_u)
  bb <- MSTest:::.mmc_boundary_blocks(fit2c)
  ym <- as.matrix(fit2c$y); ybar <- colMeans(ym); s_y <- apply(ym, 2, stats::sd)
  # per-series ranges on the transient regime's means and variances
  expect_true(all(b$theta_low[bb$rough_mu] <= ybar - 5 * s_y))
  expect_true(all(b$theta_upp[bb$rough_mu] >= ybar + 5 * s_y))
  expect_true(all(b$theta_upp[bb$rough_sig_diag] >= 4 * s_y^2))
  expect_true(all(b$theta_low[bb$rough_sig_diag] <= fit2c$sigma_floor))
  # the covariance range is exactly +/- sqrt(U_i U_j) from the REALIZED
  # diagonal upper ends
  U <- b$theta_upp[bb$rough_sig_diag]
  ob <- sqrt(b$theta_upp[bb$off_diag_i] * b$theta_upp[bb$off_diag_j])
  expect_equal(unname(b$theta_upp[bb$rough_sig_off]), unname(ob), tolerance = 1e-12)
  expect_equal(unname(b$theta_low[bb$rough_sig_off]), unname(-ob), tolerance = 1e-12)
  # absorbing column reaches the face; the clean block keeps its ordinary box
  pos <- matrix(which(as.numeric(fit2c$theta_P_ind) == 1), 2, 2)
  expect_identical(b$theta_upp[pos[2, 2]], 1)
  b0 <- MMC_bounds(fit2c, utils::modifyList(con_u, list(box_boundary = NULL)))
  expect_identical(b$theta_low[bb$clean_idx], b0$theta_low[bb$clean_idx])
  expect_identical(b$theta_upp[bb$clean_idx], b0$theta_upp[bb$clean_idx])
  # a degenerate variance range (floor above the upper end) errors by name
  expect_error(MMC_bounds(fit2c, utils::modifyList(con_u, list(boundary_R = list(sigma_mult = 1e-9)))),
               "variance range is empty")
  # the covariance bound follows the REALIZED diagonal uppers, not the literal
  # sigma_mult * s_i * s_j: when the CI union pushes one diagonal's upper end
  # beyond its moment range, the covariance bound must widen with it -- the
  # two formulas differ exactly here
  fit_wide <- fit2c
  d1 <- bb$rough_sig_diag[1]
  fit_wide$theta_se[d1] <- 10 * (4 * s_y[1]^2)   # CI-dominant upper on series 1
  b_w <- MMC_bounds(fit_wide, con_u)
  U_w <- b_w$theta_upp[bb$rough_sig_diag]
  expect_gt(U_w[1], 4 * s_y[1]^2)                # the premise: CI dominates
  ob_w <- sqrt(b_w$theta_upp[bb$off_diag_i] * b_w$theta_upp[bb$off_diag_j])
  expect_equal(unname(b_w$theta_upp[bb$rough_sig_off]), unname(ob_w), tolerance = 1e-12)
  expect_gt(unname(b_w$theta_upp[bb$rough_sig_off])[1], 4 * s_y[1] * s_y[2])
  # the PSD gate: every coordinate inside its range, yet the assembled
  # transient covariance is not positive definite (small variances, a
  # covariance near the Cauchy-Schwarz bound of the UPPER ends) -- excluded
  # before the objective is reached, counted, no delegation
  calls <- 0L
  mock <- function(theta_full, ...){ calls <<- calls + 1L; 0.5 }
  cnt <- new.env(); cnt$pd_hits <- 0L
  sig_pos <- which(as.numeric(fit2c$theta_sig_ind) == 1)
  blocks_psd <- split(sig_pos, rep(1:2, each = 3L))
  wrapped <- MSTest:::.mmc_admissible_wrap(mock, blocks_psd, 2L, 100, cnt)
  th_bad <- as.numeric(fit2c$theta)
  th_bad[bb$rough_sig_diag] <- c(0.1, 0.1)
  th_bad[bb$rough_sig_off]  <- 0.9 * sqrt(U[1] * U[2])
  expect_true(all(th_bad[bb$rough_sig_off] <= b$theta_upp[bb$rough_sig_off]))
  expect_identical(wrapped(th_bad), 100)
  expect_identical(cnt$pd_hits, 1L)
  expect_identical(calls, 0L)
  expect_identical(wrapped(as.numeric(fit2c$theta)), 0.5)
  expect_identical(calls, 1L)
})

test_that("positive-definiteness repair: deterministic off-diagonal shrinkage, diagonals untouched, gate-admissible", {
  B <- matrix(c(0.1, 0.5, 0.5, 0.1), 2, 2)
  r1 <- MSTest:::.psd_repair_shrink(B)
  expect_lt(r1$s, 1)
  expect_identical(diag(r1$Sigma), diag(B))
  expect_identical(r1, MSTest:::.psd_repair_shrink(B))
  # The repaired matrix must clear the search gate's OWN admissibility
  # criterion -- being barely positive definite is not enough, since a seed at
  # the singular boundary would be penalized at its own first evaluation.
  lam <- eigen(r1$Sigma, symmetric = TRUE, only.values = TRUE)$values
  expect_gt(min(lam), MSTest:::.mmc_psd_tol(max(lam)))
  delegated <- 0L
  mock <- function(theta_full, ...){ delegated <<- delegated + 1L; 0.5 }
  cnt <- new.env(); cnt$pd_hits <- 0L
  wrapped <- MSTest:::.mmc_admissible_wrap(mock, list(1:3), 2L, 100, cnt)
  expect_identical(wrapped(covar_vech(r1$Sigma)), 0.5)
  expect_identical(delegated, 1L)
  expect_identical(cnt$pd_hits, 0L)
  # a harder case: a larger-scale non-PSD matrix repairs to gate-admissible too
  B2 <- matrix(c(1, 1.6, 1.6, 2), 2, 2)
  r2 <- MSTest:::.psd_repair_shrink(B2)
  lam2 <- eigen(r2$Sigma, symmetric = TRUE, only.values = TRUE)$values
  expect_gt(min(lam2), MSTest:::.mmc_psd_tol(max(lam2)))
  # an already-comfortably-PSD matrix passes through unchanged
  C <- matrix(c(1, 0.3, 0.3, 1), 2, 2)
  expect_identical(MSTest:::.psd_repair_shrink(C)$Sigma, C)
})

test_that("multivariate boundary end-to-end: union and face run with the harvest fields", {
  skip_on_cran()
  set.seed(48)
  y2 <- as.matrix(simuMSVAR(list(n = 150, p = 1, q = 2, mu = rbind(c(0, 0), c(2, -1)),
                                 sigma = list(diag(2), diag(2) * 2), phi = diag(2) * 0.2,
                                 k = 2, P = matrix(c(0.9, 0.1, 0, 1), 2, 2)), burnin = 0)$y)
  ctl1 <- list(use_diff_init = 1)
  m_u <- MMCLRTest(y2, p = 1, k0 = 2, k1 = 3,
                   control = list(N = 19, mc_seed = 90, maxit = 3, silence = TRUE, burnin = 0,
                                  box_boundary = "union",
                                  mdl_h0_control = ctl1, mdl_h1_control = ctl1))
  expect_true(valid_pval(m_u$pval))
  expect_gte(m_u$pval, m_u$pval_0)
  expect_true(is.integer(m_u$optim_diagnostics$pd_penalty_hits))
  expect_gte(m_u$optim_diagnostics$pd_penalty_hits, 0L)
  m_f <- MMCLRTest(y2, p = 1, k0 = 2, k1 = 3,
                   control = list(N = 19, mc_seed = 91, maxit = 2, silence = TRUE, burnin = 0,
                                  box_boundary = "face",
                                  mdl_h0_control = ctl1, mdl_h1_control = ctl1))
  expect_true(valid_pval(m_f$pval))
  expect_gte(m_f$pval, m_f$pval_proj)
  expect_true(is.list(m_f$seed_repair))
  expect_true(is.logical(m_f$seed_repair$repaired))
  # face searches the rough block only: per-series means (2) + vech (3) + p11 (1)
  expect_identical(ncol(m_f$optim_diagnostics$evaluated_range), 6L)
})

test_that("union mode with exogenous regressors: the common coefficients stay clean", {
  skip_on_cran()
  set.seed(49)
  Zb <- matrix(rnorm(200), 200, 1)
  y_bz <- as.matrix(simuMSAR(list(n = 200, k = 2, mu = c(0, 2), sigma = c(1, 2), phi = 0.3,
                                  P = matrix(c(0.9, 0.1, 0, 1), 2, 2)), burnin = 0)$y + 1.5 * Zb)
  fitx <- MSARXmdl(y_bz, p = 1, k = 2, Z = Zb, control = list(getSE = TRUE, use_diff_init = 2))
  fitx_c <- canonicalize_regimes(fitx, by = "persistence")
  con_u <- list(eps = 0.1, CI_union = TRUE, P_low = 0, P_upp = 1,
                phi_low = NULL, phi_upp = NULL, variance_constraint = 0.01,
                type = "pso", box_include = NULL,
                box_boundary = "union", boundary_R = NULL)
  b  <- MMC_bounds(fitx_c, con_u)
  b0 <- MMC_bounds(fitx_c, utils::modifyList(con_u, list(box_boundary = NULL)))
  x_pos <- which(as.numeric(fitx_c$theta_x_ind) == 1)
  expect_identical(b$theta_low[x_pos], b0$theta_low[x_pos])
  expect_identical(b$theta_upp[x_pos], b0$theta_upp[x_pos])
  m <- MMCLRTest(y_bz, p = 1, k0 = 2, k1 = 3, Z = Zb,
                 control = list(N = 19, mc_seed = 92, maxit = 2, silence = TRUE, burnin = 0,
                                box_boundary = "union",
                                mdl_h0_control = ctl2, mdl_h1_control = ctl2))
  expect_true(valid_pval(m$pval))
  expect_gte(m$pval, m$pval_0)
})

test_that("search diagnostics: recorder count matches pso's own, early stop marked", {
  skip_on_cran()
  fit_seed <- 71
  m <- MMCLRTest(y_ur, p = 1, k0 = 1, k1 = 2,
                 control = list(N = 19, mc_seed = fit_seed, maxit = 6, silence = TRUE,
                                mdl_h1_control = ctl2))
  expect_identical(m$optim_diagnostics$type, "pso")
  expect_identical(m$optim_diagnostics$feval,
                   as.integer(m$mmc_optimout$counts[["function"]]))
  expect_identical(m$optim_diagnostics$feval, m$optim_diagnostics$optimizer_feval)
  # early stop: threshold_stop below any p-value floor forces the marker path
  m_es <- MMCLRTest(y_ur, p = 1, k0 = 1, k1 = 2,
                    control = list(N = 19, mc_seed = fit_seed, maxit = 6, silence = TRUE,
                                   threshold_stop = 1/(19 + 1),
                                   mdl_h1_control = ctl2))
  expect_identical(m_es$optim_diagnostics$type, "early_stop")
  expect_identical(m_es$optim_diagnostics$feval, 0L)
})
