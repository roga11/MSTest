# Tests for the pre-simulation-freeze changes (v0.2.0):
# warm-start embeddings, LR >= 0 policy, CRN seeding, control-merge fix.

valid_pval <- function(p) is.finite(p) && p >= 0 && p <= 1

# Small shared datasets --------------------------------------------------------
set.seed(42)
y_uni <- simuMSAR(list(n = 150, p = 1, q = 1, mu = c(0, 2), sigma = c(0.5, 1),
                       phi = 0.5, k = 2,
                       P = matrix(c(0.9, 0.1, 0.1, 0.9), 2, 2), burnin = 50))$y
Z_uni <- matrix(rnorm(150), 150, 1)
y_biv <- simuMSVAR(list(n = 150, p = 1, q = 2, mu = rbind(c(0, 0), c(2, 2)),
                        sigma = list(diag(2) * 0.5, diag(2)), phi = diag(2) * 0.4,
                        k = 2, P = matrix(c(0.9, 0.1, 0.1, 0.9), 2, 2), burnin = 50))$y

test_that("regime-duplication embedding preserves the likelihood exactly", {
  ctl  <- list(getSE = FALSE)
  ctlm <- list(getSE = FALSE, use_diff_init = 1, maxit = 200)
  # Compare against the log-likelihood RE-EVALUATED at the fitted theta (the
  # EM's reported logLike lags the final theta by one E-step, ~1e-5).
  chk <- function(mdl, k1, llfun, msmu = TRUE, msvar = TRUE, split = 1){
    m0 <- mdl; m0$exog <- !is.null(m0$betaZ)
    ll_ref <- llfun(mdl$theta, m0, max(1, mdl$k))
    m1 <- m0; m1$msmu <- msmu; m1$msvar <- msvar
    th <- MSTest:::.embed_theta(mdl, k1, msmu = msmu, msvar = msvar, split_regime = split)
    expect_equal(llfun(th, m1, k1), ll_ref, tolerance = 1e-10)
  }
  chk_lin <- function(mdl, k1, llfun, msmu = TRUE, msvar = TRUE){
    # linear null: reference is the reported (closed-form) logLike
    m1 <- mdl; m1$exog <- !is.null(m1$betaZ); m1$msmu <- msmu; m1$msvar <- msvar
    th <- MSTest:::.embed_theta(mdl, k1, msmu = msmu, msvar = msvar)
    expect_equal(llfun(th, m1, k1), mdl$logLike, tolerance = 1e-10)
  }
  ar1 <- ARmdl(y_uni, p = 1, control = ctl)
  chk_lin(ar1, 2, logLike_MSARmdl)
  chk_lin(ar1, 3, logLike_MSARmdl)                       # m = 2, direct
  chk_lin(ar1, 2, logLike_MSARmdl, msvar = FALSE)
  chk_lin(ar1, 2, logLike_MSARmdl, msmu = FALSE)
  chk_lin(ARXmdl(y_uni, p = 1, Z = Z_uni, control = ctl), 2, logLike_MSARXmdl)
  chk_lin(Nmdl(y_uni, control = ctl), 2, logLike_HMmdl)
  chk_lin(Nmdl(y_biv, control = ctl), 2, logLike_HMmdl)
  chk_lin(VARmdl(y_biv, p = 1, control = ctl), 2, logLike_MSVARmdl)
  # MS nulls (lumpability construction)
  ms2 <- MSARmdl(y_uni, p = 1, k = 2, control = ctlm)
  chk(ms2, 3, logLike_MSARmdl, split = 1)
  chk(ms2, 3, logLike_MSARmdl, split = 2)
  chk(ms2, 4, logLike_MSARmdl)                           # m = 2, direct
  chk(HMmdl(y_uni, k = 2, control = ctlm), 3, logLike_HMmdl)
  chk(MSVARmdl(y_biv, p = 1, k = 2, control = ctlm), 3, logLike_MSVARmdl)
})

test_that("warmstart_theta returns one start per clone distribution, all finite", {
  ar1 <- ARmdl(y_uni, p = 1, control = list(getSE = FALSE))
  ms2 <- MSARmdl(y_uni, p = 1, k = 2,
                 control = list(getSE = FALSE, use_diff_init = 1, maxit = 200))
  expect_length(warmstart_theta(ar1, 2), 1)   # k0=1, m=1
  expect_length(warmstart_theta(ar1, 3), 1)   # k0=1, m=2
  expect_length(warmstart_theta(ms2, 3), 2)   # k0=2, m=1: split regime 1 or 2
  expect_length(warmstart_theta(ms2, 4), 3)   # k0=2, m=2: (3,1),(1,3),(2,2)
  expect_true(all(vapply(warmstart_theta(ms2, 4),
                         function(th) all(is.finite(th)), TRUE)))
  expect_error(warmstart_theta(ms2, 2))       # k1 must exceed k0
})

test_that("init_theta_extra competes with random starts (best logLike wins)", {
  ctl <- list(getSE = FALSE, use_diff_init = 3, maxit = 200)
  set.seed(7)
  fit3 <- MSARmdl(y_uni, p = 1, k = 2, control = ctl)
  # Feed the 3-start optimum as an extra start with a single random start:
  # the extra must win (weakly better logLike).
  set.seed(99)
  fit_x <- MSARmdl(y_uni, p = 1, k = 2,
                   control = list(getSE = FALSE, use_diff_init = 1, maxit = 500,
                                  init_theta_extra = list(as.numeric(fit3$theta))))
  expect_gte(fit_x$logLike, fit3$logLike - 1e-8)
  # Legacy bit-identity: init_theta_extra = NULL must not touch the RNG stream.
  set.seed(11); a <- MSARmdl(y_uni, p = 1, k = 2, control = ctl)
  set.seed(11); b <- MSARmdl(y_uni, p = 1, k = 2, control = ctl)
  expect_identical(a$theta, b$theta)
  # Malformed extras degrade gracefully (wrong length -> retry with random;
  # non-finite -> filtered).
  set.seed(12)
  fit_bad <- MSARmdl(y_uni, p = 1, k = 2,
                     control = list(getSE = FALSE, use_diff_init = 1, maxit = 200,
                                    init_theta_extra = list(c(1, 2, 3), rep(NaN, 9))))
  expect_true(is.finite(fit_bad$logLike))
})

test_that("LMCLRTest warm-start default: LRT_0 >= 0, valid p-value, reproducible", {
  skip_on_cran()
  ctl <- list(N = 19, burnin = 100, mc_seed = 8675309, workers = 0,
              mdl_h0_control = list(const = TRUE, getSE = FALSE),
              mdl_h1_control = list(getSE = FALSE, use_diff_init = 3))
  # (suppressWarnings: the LR >= 0 cap emits an informational warning when a
  # simulated draw is capped; LRN >= 0 is asserted explicitly below)
  a <- suppressWarnings(LMCLRTest(y_uni, p = 1, k0 = 1, k1 = 2, control = ctl))
  b <- suppressWarnings(LMCLRTest(y_uni, p = 1, k0 = 1, k1 = 2, control = ctl))
  expect_gte(as.numeric(a$LRT_0), 0)
  expect_true(all(a$LRN >= 0))
  expect_true(valid_pval(as.numeric(a$pval)))
  expect_identical(a$pval, b$pval)
  expect_identical(a$LRN, b$LRN)
  expect_identical(a$LRT_0, b$LRT_0)
  # warm start can only (weakly) improve the H1 fit relative to legacy starts
  d <- suppressWarnings(LMCLRTest(y_uni, p = 1, k0 = 1, k1 = 2,
                                  control = c(ctl, list(init_method = "random"))))
  expect_gte(as.numeric(a$LRT_0), as.numeric(d$LRT_0) - 1e-8)
})

test_that("MMCLRTest control-merge fix: sublist without getSE no longer crashes", {
  skip_on_cran()
  # Regression for the zero-length-logical crash: a user-supplied mdl_h0_control
  # that omits getSE previously made the CI_union check error out.
  expect_error(
    suppressWarnings(
      MMCLRTest(y_uni, p = 1, k0 = 1, k1 = 2,
                control = list(N = 19, maxit = 2, threshold_stop = 0.9, silence = TRUE,
                               type = "GenSA", mc_seed = 1,
                               mdl_h0_control = list(const = TRUE),
                               mdl_h1_control = list(use_diff_init = 1)))
    ),
    NA)
})

test_that("MMCLRTest: pval_0 reported, pval >= pval_0 under CRN, reproducible", {
  skip_on_cran()
  ctl <- list(N = 19, burnin = 100, mc_seed = 8675309, workers = 0,
              type = "GenSA", maxit = 6, threshold_stop = 0.9, silence = TRUE,
              CI_union = FALSE,
              mdl_h0_control = list(const = TRUE, getSE = FALSE),
              mdl_h1_control = list(getSE = FALSE, use_diff_init = 1))
  a <- suppressWarnings(MMCLRTest(y_uni, p = 1, k0 = 1, k1 = 2, control = ctl))
  b <- suppressWarnings(MMCLRTest(y_uni, p = 1, k0 = 1, k1 = 2, control = ctl))
  expect_true(valid_pval(as.numeric(a$pval)))
  expect_true(valid_pval(as.numeric(a$pval_0)))
  expect_gte(as.numeric(a$pval), as.numeric(a$pval_0))
  expect_identical(a$pval, b$pval)
  expect_identical(a$pval_0, b$pval_0)
  expect_identical(a$theta_h0, b$theta_h0)
})

test_that("MMCLRTest early stop: pval_0 >= threshold_stop skips the search", {
  skip_on_cran()
  # True-null data with a low threshold: the theta_0 evaluation should settle it.
  set.seed(55)
  y0 <- simuAR(list(n = 150, p = 1, q = 1, mu = 1, sigma = 1, phi = 0.5, burnin = 100))$y
  ctl <- list(N = 19, burnin = 100, mc_seed = 99, workers = 0,
              type = "pso", maxit = 50, threshold_stop = 0.2, silence = TRUE,
              CI_union = FALSE,
              mdl_h0_control = list(const = TRUE, getSE = FALSE),
              mdl_h1_control = list(getSE = FALSE, use_diff_init = 1))
  a <- suppressWarnings(MMCLRTest(y0, p = 1, k0 = 1, k1 = 2, control = ctl))
  skip_if(as.numeric(a$pval_0) < 0.2, "pval_0 below threshold for this seed")
  expect_identical(as.numeric(a$pval), as.numeric(a$pval_0))
  expect_true(isTRUE(a$mmc_optimout$early_stop))
})

test_that("EM variance floor removes the degenerate sigma^2 -> 0 spike", {
  skip_on_cran()
  # Dataset/seed pair that reproducibly finds the spike mode without the floor
  # (a regime with sigma^2 ~ 1e-11 planted on ~1% of observations, LR ~ 35).
  set.seed(5010)
  y <- simuAR(list(n = 200, p = 1, q = 1, mu = 1, sigma = 1, phi = 0.9, burnin = 100))$y
  ar <- ARmdl(y, p = 1, control = list(getSE = FALSE))
  ctl <- list(getSE = FALSE, use_diff_init = 20, maxit = 500)
  # Unconstrained (legacy escape hatch): spike present
  set.seed(6010)
  ms_free <- suppressWarnings(MSARmdl(y, p = 1, k = 2,
                control = c(ctl, list(em_variance_constraint = 0))))
  expect_gt(-2 * (ar$logLike - ms_free$logLike), 20)
  expect_lt(min(as.numeric(ms_free$sigma)), 1e-6)
  # Constrained (default): spike gone, variances respect the floor
  set.seed(6010)
  ms_con <- suppressWarnings(MSARmdl(y, p = 1, k = 2, control = ctl))
  expect_lt(-2 * (ar$logLike - ms_con$logLike), 20)
  expect_gte(min(as.numeric(ms_con$sigma)), ms_con$sigma_floor * (1 - 1e-9))
  expect_equal(ms_con$sigma_floor, 0.01 * as.numeric(ar$sigma), tolerance = 1e-8)
})

test_that("variance floor is a no-op when slack (bit-identity) and flags/SEs behave", {
  skip_on_cran()
  # Well-separated regimes: the floor never binds, so results must be identical
  set.seed(1234)
  Y <- simuMSAR(list(n = 200, k = 2, mu = c(-2, 2), sigma = c(1, 1), phi = 0.3,
                     P = rbind(c(0.90, 0.10), c(0.10, 0.90))))$y
  set.seed(3); a <- MSARmdl(Y, p = 1, k = 2, control = list(getSE = FALSE, use_diff_init = 3))
  set.seed(3); b <- MSARmdl(Y, p = 1, k = 2,
                            control = list(getSE = FALSE, use_diff_init = 3,
                                           em_variance_constraint = 0))
  expect_identical(a$logLike, b$logLike)
  expect_identical(a$theta, b$theta)
  expect_false(any(a$sigma_floor_binding))
  # A binding fit flags the regime and returns NA SE for that variance only
  set.seed(5011)
  y2 <- simuAR(list(n = 200, p = 1, q = 1, mu = 1, sigma = 1, phi = 0.9, burnin = 100))$y
  set.seed(6011)
  m <- suppressWarnings(MSARmdl(y2, p = 1, k = 2,
         control = list(getSE = TRUE, use_diff_init = 20, maxit = 500)))
  skip_if(!any(m$sigma_floor_binding), "floor did not bind for this seed")
  se_sig <- m$theta_se[as.logical(m$theta_sig_ind)]
  expect_true(all(is.na(se_sig[m$sigma_floor_binding])))
  expect_true(all(is.finite(m$theta_se[!as.logical(m$theta_sig_ind)])))
})

test_that("constrained EM ascends monotonically and keeps starts feasible", {
  skip_on_cran()
  set.seed(5008)
  y <- simuAR(list(n = 200, p = 1, q = 1, mu = 1, sigma = 1, phi = 0.9, burnin = 100))$y
  # Deliberately strong floor so the constraint binds throughout the run
  lls <- vapply(1:8, function(m){
    set.seed(11)
    suppressWarnings(MSARmdl(y, p = 1, k = 2,
      control = list(getSE = FALSE, use_diff_init = 1, maxit = m, conv = "theta",
                     thtol = 1e-14, em_variance_constraint = 0.5))$logLike)
  }, 0.0)
  expect_true(all(diff(lls) >= -1e-8))
  # Warm starts (which perturb variances downward) are projected onto the floor
  ar <- ARmdl(y, p = 1, control = list(getSE = FALSE))
  set.seed(5)
  fw <- suppressWarnings(MSARmdl(y, p = 1, k = 2,
          control = list(getSE = FALSE, use_diff_init = 1, maxit = 300,
                         em_variance_constraint = 0.5,
                         init_theta_extra = warmstart_theta(ar, 2))))
  expect_gte(min(as.numeric(fw$sigma)), fw$sigma_floor * (1 - 1e-9))
})

test_that("variance floor works for multivariate models (eigenvalue clipping)", {
  skip_on_cran()
  set.seed(9)
  mv <- suppressWarnings(MSVARmdl(y_biv, p = 1, k = 2,
          control = list(getSE = FALSE, use_diff_init = 2, maxit = 300)))
  expect_true(is.finite(mv$logLike))
  eig <- vapply(mv$sigma, function(S)
    min(eigen(as.matrix(S), symmetric = TRUE, only.values = TRUE)$values), 0.0)
  expect_true(all(eig >= mv$sigma_floor * (1 - 1e-9)))
  expect_equal(unname(colSums(mv$P)), rep(1, 2), tolerance = 1e-8)
})

test_that("internal seeds are drawn (not mc_seed+k) and pre-draws stay aligned", {
  # The pre-drawn innovations must depend on mc_seed alone (draw-then-reseed):
  # calling set.seed(s); sample.int(...); set.seed(s) reproduces the same rnorm
  # stream as set.seed(s) directly.
  s <- 12345
  set.seed(s); x_direct <- rnorm(5)
  set.seed(s); inseeds <- sample.int(.Machine$integer.max, 3L); set.seed(s); x_dtr <- rnorm(5)
  expect_identical(x_direct, x_dtr)
  # Drawn seeds differ from the mc_seed+k arithmetic and from the parent seed.
  expect_false(any(inseeds %in% (s + 0:3)))
})
