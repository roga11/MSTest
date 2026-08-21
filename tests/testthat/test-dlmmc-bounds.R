# DLMMC_bounds() mirrors MMC_bounds()'s guards: a non-finite standard error only
# falls back to eps for that ONE coefficient, not the whole box (previously all-or-
# nothing, via all(is.finite(phiSE))); 'eps' may be a per-coefficient vector; and the
# box must never invert or exclude the fitted phi.

.dlmmc_mdl <- function(phi = c(0.5, -0.2), se = c(0.05, 0.05)){
  list(phi = phi, theta_phi_ind = rep(1, length(phi)), theta_se = se)
}
.dlmmc_con <- function(eps, CI_union = FALSE, phi_low = NULL, phi_upp = NULL, optim_type = NULL){
  list(eps = eps, CI_union = CI_union, phi_low = phi_low, phi_upp = phi_upp, optim_type = optim_type)
}

test_that("a scalar eps reproduces the documented default behaviour", {
  mdl <- .dlmmc_mdl()
  out <- DLMMC_bounds(mdl, .dlmmc_con(0.1))
  expect_equal(out$theta_low, mdl$phi - 0.1)
  expect_equal(out$theta_upp, mdl$phi + 0.1)
})

test_that("a vector eps produces hand-computed per-coefficient bounds", {
  mdl <- .dlmmc_mdl()
  ev <- c(0.05, 0.2)
  out <- DLMMC_bounds(mdl, .dlmmc_con(ev))
  expect_equal(out$theta_low, mdl$phi - ev)
  expect_equal(out$theta_upp, mdl$phi + ev)
})

test_that("eps of the wrong length errors, naming both lengths", {
  mdl <- .dlmmc_mdl(phi = c(0.5, -0.2, 0.1))
  expect_error(DLMMC_bounds(mdl, .dlmmc_con(c(0.1, 0.2))),
              "length 1 or length\\(theta_0\\) = 3")
})

test_that("non-finite or negative eps entries are rejected", {
  mdl <- .dlmmc_mdl()
  expect_error(DLMMC_bounds(mdl, .dlmmc_con(NA_real_)), "finite")
  expect_error(DLMMC_bounds(mdl, .dlmmc_con(Inf)), "finite")
  expect_error(DLMMC_bounds(mdl, .dlmmc_con(-0.1)), "non-negative")
})

test_that("eps = 0 is valid: a zero-width box with no optimizer specified, nudged for GenSA", {
  mdl <- .dlmmc_mdl()
  out <- expect_silent(DLMMC_bounds(mdl, .dlmmc_con(0)))
  expect_equal(out$theta_low, mdl$phi)
  expect_equal(out$theta_upp, mdl$phi)
  # optim_type = "GenSA" (DLMMCTest's own default): a zero entry is nudged to
  # 1e-8, since GenSA errors on a genuinely zero-width bound.
  out2 <- expect_silent(DLMMC_bounds(mdl, .dlmmc_con(0, optim_type = "GenSA")))
  expect_equal(out2$theta_upp - out2$theta_low, rep(2e-8, length(mdl$phi)))
  # pso is not nudged -- it handles a zero-width bound fine.
  out3 <- expect_silent(DLMMC_bounds(mdl, .dlmmc_con(0, optim_type = "pso")))
  expect_equal(out3$theta_low, mdl$phi)
})

test_that("a non-finite SE on one coefficient only costs that coefficient the CI union", {
  # The actual E3 fix: pre-fix, all(is.finite(phiSE)) would skip the union for BOTH
  # coefficients because of coefficient 2's NA. Post-fix, coefficient 1 still gets
  # the wider of eps and 2*SE; only coefficient 2 falls back to eps alone.
  # se[1] = 0.2 so 2*se = 0.4 > eps = 0.1: the union must actually widen coefficient 1,
  # not merely tie eps (0.05 would leave 2*se == eps, unable to distinguish the fix).
  mdl <- .dlmmc_mdl(phi = c(0.5, -0.2), se = c(0.2, NA_real_))
  out <- DLMMC_bounds(mdl, .dlmmc_con(0.1, CI_union = TRUE))
  expect_equal(out$theta_low[1], 0.5 - 2 * 0.2)  # 2*se=0.4 > eps=0.1: union widens coefficient 1
  expect_equal(out$theta_upp[1], 0.5 + 2 * 0.2)
  expect_equal(out$theta_low[2], -0.2 - 0.1)  # NA SE -> eps alone, NOT collapsed for coef 1 too
  expect_equal(out$theta_upp[2], -0.2 + 0.1)
})

test_that("CI_union widens the box beyond eps alone when every SE is finite", {
  mdl <- .dlmmc_mdl(phi = c(0.5), se = c(0.2))
  out <- DLMMC_bounds(mdl, .dlmmc_con(0.1, CI_union = TRUE))
  expect_equal(out$theta_low, 0.5 - 2*0.2)
  expect_equal(out$theta_upp, 0.5 + 2*0.2)
})

test_that("a NULL theta_se (getSE = FALSE) degrades to eps alone under CI_union = TRUE", {
  mdl <- .dlmmc_mdl(phi = c(0.5), se = NULL)
  out <- expect_silent(DLMMC_bounds(mdl, .dlmmc_con(0.1, CI_union = TRUE)))
  expect_equal(out$theta_low, 0.5 - 0.1)
  expect_equal(out$theta_upp, 0.5 + 0.1)
})

test_that("an asymmetric phi clamp that inverts the box is rejected", {
  mdl <- .dlmmc_mdl(phi = c(0.422))
  expect_error(DLMMC_bounds(mdl, .dlmmc_con(0.05, phi_low = 0.8)),
              "lower bound above its upper bound")
})

test_that("a phi clamp that excludes theta_0, without inverting the box, is rejected", {
  mdl <- .dlmmc_mdl(phi = c(0.422))
  expect_error(DLMMC_bounds(mdl, .dlmmc_con(0.05, phi_upp = 0.4)),
              "falls outside the search box")
})

test_that("a non-finite phi_low/phi_upp errors by name, not with an opaque crash", {
  mdl <- .dlmmc_mdl()
  expect_error(DLMMC_bounds(mdl, .dlmmc_con(0.1, phi_low = NA_real_)),
               "'phi_low' must be finite")
  expect_error(DLMMC_bounds(mdl, .dlmmc_con(0.1, phi_upp = NA_real_)),
               "'phi_upp' must be finite")
})

test_that("a non-finite theta_0 (a pathological upstream fit) errors informatively instead of crashing on if(NA)", {
  mdl <- .dlmmc_mdl(phi = c(0.5, NaN))
  expect_error(DLMMC_bounds(mdl, .dlmmc_con(0.1)),
               "DLMMC_bounds: the search box is non-finite")
})

test_that("an ordinary box is returned silently, and contains the fitted coefficients", {
  mdl <- .dlmmc_mdl()
  out <- expect_silent(DLMMC_bounds(mdl, .dlmmc_con(0.1)))
  expect_true(all(out$theta_low <= out$theta_upp))
  expect_true(all(mdl$phi >= out$theta_low & mdl$phi <= out$theta_upp))
})

test_that("DLMMCTest itself rejects a bad eps before reporting a result", {
  skip_on_cran()
  set.seed(4242)
  y <- matrix(rnorm(150), 150, 1)
  expect_error(
    DLMMCTest(y, p = 1, control = list(N = 9, eps = c(0.1, 0.2))),
    "length 1 or length\\(theta_0\\)")
  expect_error(
    DLMMCTest(y, p = 1, control = list(N = 9, eps = -0.1)),
    "non-negative")
})
