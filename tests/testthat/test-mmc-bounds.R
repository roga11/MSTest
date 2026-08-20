# MMC_bounds() must not hand an inverted search box to the optimizers: GenSA errors on
# one, GA dies on a degenerate population from another, and pso silently searches the
# inverted interval and reports a p-value for it. The reachable cause is not 'eps' alone
# (which only ever widens the box symmetrically) but an asymmetric 'phi_low'/'phi_upp' or
# 'P_low'/'P_upp' clamp, which adjusts only one side of the pair.

.mmc_bounds_mdl <- function(phi_hat = 0.422){
  # A minimal MSAR(1), k = 2 parameter vector -- mu(2) + phi(1) + sigma(2) + vec(P)(4) --
  # built by hand so the guard can be tested without paying for an actual EM fit.
  theta <- c(0, 3, phi_hat, 1, 1, 0.9, 0.1, 0.1, 0.9)
  names(theta) <- c("mu_1", "mu_2", "phi_1", "sig_1", "sig_2", "P_11", "P_21", "P_12", "P_22")
  list(theta = theta, k = 2L, p = 1L, theta_se = NULL,
       theta_var_ind = c(0, 0, 0, 1, 1, 0, 0, 0, 0),
       theta_phi_ind = c(0, 0, 1, 0, 0, 0, 0, 0, 0),
       theta_P_ind   = c(0, 0, 0, 0, 0, 1, 1, 1, 1))
}

test_that("an asymmetric phi clamp that inverts the box is rejected", {
  mdl <- .mmc_bounds_mdl(phi_hat = 0.422)
  # eps = 0.05 gives phi in [0.372, 0.472]; phi_low = 0.8 raises only the lower end.
  con <- list(eps = 0.05, CI_union = FALSE, variance_constraint = 0.01,
              P_low = 0, P_upp = 1, phi_low = 0.8, phi_upp = NULL)
  expect_error(MMC_bounds(mdl, con), "lower bound above its upper bound")
  expect_error(MMC_bounds(mdl, con), "phi_1")
})

test_that("an asymmetric P clamp that inverts the box is rejected", {
  mdl <- .mmc_bounds_mdl()
  # eps = 0.05 gives P_11 in [0.85, 0.95]; P_upp = 0.6 lowers only the upper end.
  con <- list(eps = 0.05, CI_union = FALSE, variance_constraint = 0.01,
              P_low = 0, P_upp = 0.6, phi_low = NULL, phi_upp = NULL)
  expect_error(MMC_bounds(mdl, con), "lower bound above its upper bound")
})

test_that("an ordinary box is returned silently, and contains the fitted parameters", {
  mdl <- .mmc_bounds_mdl()
  con <- list(eps = 0.1, CI_union = FALSE, variance_constraint = 0.01,
              P_low = 0, P_upp = 1, phi_low = NULL, phi_upp = NULL)
  out <- expect_silent(MMC_bounds(mdl, con))
  expect_true(all(out$theta_low <= out$theta_upp))
  expect_true(all(mdl$theta >= out$theta_low & mdl$theta <= out$theta_upp))
})

test_that("a P clamp that excludes theta_0, without inverting the box, is rejected", {
  # P_11 = 0.9; eps = 0.05 gives the natural interval [0.85, 0.95]; P_upp = 0.85 lowers
  # the upper end to exactly the natural lower end, so low == upp (not inverted) but
  # theta_0 = 0.9 sits above both.
  mdl <- .mmc_bounds_mdl()
  con <- list(eps = 0.05, CI_union = FALSE, variance_constraint = 0.01,
              P_low = 0, P_upp = 0.85, phi_low = NULL, phi_upp = NULL)
  out <- expect_error(MMC_bounds(mdl, con), "falls outside the search box")
})

test_that("a phi clamp that excludes theta_0, without inverting the box, is rejected", {
  # phi_1 = 0.422; eps = 0.05 gives [0.372, 0.472]; phi_upp = 0.4 lowers the upper end
  # below theta_0 while staying above the natural lower end -- not inverted, but
  # theta_0 is now outside.
  mdl <- .mmc_bounds_mdl(phi_hat = 0.422)
  con <- list(eps = 0.05, CI_union = FALSE, variance_constraint = 0.01,
              P_low = 0, P_upp = 1, phi_low = NULL, phi_upp = 0.4)
  expect_error(MMC_bounds(mdl, con), "falls outside the search box")
  expect_error(MMC_bounds(mdl, con), "phi_1")
})

test_that("MMCLRTest itself surfaces the theta_0-outside-box error, not just an inverted one", {
  skip_on_cran()
  set.seed(9002)
  y <- matrix(rnorm(150), 150, 1)
  # phi_upp below the fitted phi but above the natural eps-based lower bound: box is
  # not inverted, but the fitted phi is excluded from it. This is the D4 repro
  # recorded in CLAUDE.md's history (pso used to silently seed outside the box).
  fit0 <- ARmdl(y, p = 1, control = list(getSE = FALSE))
  phi_hat <- fit0$phi
  expect_error(
    MMCLRTest(y, p = 1, k0 = 1, k1 = 2, control = list(
      N = 9, mc_seed = 9002, CI_union = FALSE,
      phi_upp = phi_hat - 0.01, eps = 0.1,
      mdl_h0_control = list(getSE = FALSE), mdl_h1_control = list(getSE = FALSE))),
    "falls outside the search box")
})

test_that("MMCLRTest itself surfaces the same error when the box the user asked for inverts", {
  skip_on_cran()
  set.seed(9001)
  y <- matrix(rnorm(150), 150, 1)
  expect_error(
    MMCLRTest(y, p = 1, k0 = 1, k1 = 2, control = list(
      N = 9, mc_seed = 9001, phi_low = 0.99, CI_union = FALSE,
      mdl_h0_control = list(getSE = FALSE), mdl_h1_control = list(getSE = FALSE))),
    "lower bound above its upper bound")
})
