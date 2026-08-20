# Phase 4: MMC_bounds()'s 'eps' control accepts a scalar (recycled to every
# coordinate, the original default behaviour, unchanged) or a vector of length
# length(theta_0), so the search half-width can differ per parameter. The bounds
# arithmetic already vectorizes by ordinary recycling either way, so this is mostly
# validation: reject a length that isn't exactly 1 or length(theta_0) -- a length
# that merely divides it evenly would otherwise recycle silently, with no warning,
# which is the actual failure mode a per-block user is likely to hit.

.eps_mdl <- function(k0 = 2L, other_len = 3, seed = 1){
  set.seed(seed)
  theta <- c(round(stats::runif(other_len), 3), 0.9, 0.1, 0.1, 0.9)
  names(theta) <- c(paste0("x_", seq_len(other_len)), "P_11", "P_21", "P_12", "P_22")
  list(theta = theta, k = k0, p = 0L,
       theta_var_ind = c(rep(0, other_len), 0, 0, 0, 0),
       theta_P_ind   = c(rep(0, other_len), 1, 1, 1, 1),
       theta_se = NULL)
}
.eps_con <- function(eps){
  list(eps = eps, CI_union = FALSE, variance_constraint = 0.01, P_low = 0, P_upp = 1,
      phi_low = NULL, phi_upp = NULL)
}

test_that("a scalar eps reproduces the documented default behaviour, unchanged", {
  mdl <- .eps_mdl()
  out <- MMC_bounds(mdl, .eps_con(0.1))
  expect_equal(out$theta_low, as.numeric(mdl$theta) - 0.1)
  expect_equal(out$theta_upp, as.numeric(mdl$theta) + 0.1)
})

test_that("a vector eps produces hand-computed per-coordinate bounds", {
  mdl <- .eps_mdl()
  # Small enough on the P block (theta = ...,0.9,0.1,0.1,0.9) that the pre-existing,
  # unrelated P_low/P_upp = [0,1] admissible-region clamp never fires -- that clamp
  # is exercised deliberately in its own test below.
  ev <- c(0.05, 0.2, 0.01, 0.03, 0.03, 0.03, 0.03)
  out <- MMC_bounds(mdl, .eps_con(ev))
  expect_equal(out$theta_low, as.numeric(mdl$theta) - ev)
  expect_equal(out$theta_upp, as.numeric(mdl$theta) + ev)
})

test_that("a vector of identical values equals the scalar case bitwise", {
  mdl <- .eps_mdl()
  out_scalar <- MMC_bounds(mdl, .eps_con(0.15))
  out_vector <- MMC_bounds(mdl, .eps_con(rep(0.15, length(mdl$theta))))
  expect_identical(out_scalar, out_vector)
})

test_that("eps of the wrong length errors, naming both lengths -- including lengths that merely divide length(theta_0)", {
  mdl <- .eps_mdl()   # length(theta) == 7
  expect_error(MMC_bounds(mdl, .eps_con(c(0.1, 0.2))), "length 1 or length\\(theta_0\\) = 7")
  # length 7 %% length 1 == 0 for any length, but note the DIVISOR case specifically:
  # length(theta) = 7 is prime here, so pick a model where a divisor actually exists.
  mdl9 <- .eps_mdl(other_len = 5)   # length(theta) == 9
  expect_error(MMC_bounds(mdl9, .eps_con(rep(0.1, 3))), "length 1 or length\\(theta_0\\) = 9")
})

test_that("non-finite or non-positive eps entries are rejected, not silently accepted", {
  mdl <- .eps_mdl()
  expect_error(MMC_bounds(mdl, .eps_con(NA_real_)), "finite")
  expect_error(MMC_bounds(mdl, .eps_con(Inf)), "finite")
  expect_error(MMC_bounds(mdl, .eps_con(-0.1)), "strictly positive")
  expect_error(MMC_bounds(mdl, .eps_con(0)), "strictly positive")
  ev <- rep(0.1, length(mdl$theta)); ev[3] <- 0
  expect_error(MMC_bounds(mdl, .eps_con(ev)), "strictly positive")
  ev2 <- rep(0.1, length(mdl$theta)); ev2[3] <- NA_real_
  expect_error(MMC_bounds(mdl, .eps_con(ev2)), "finite")
})

test_that("a vector eps composes correctly with the P clamp", {
  # P_low/P_upp are GLOBAL: they apply to every P entry, not just the one under
  # test. Since MMC_bounds also requires theta_0 to stay inside its own box (this
  # session's earlier pre-flight check), a P_low above any entry's own fitted value
  # would exclude that entry structurally, regardless of eps -- so this necessarily
  # targets an entry with a small fitted value (P_21 = 0.1) and a small P_low that
  # stays below every entry's fitted value (0.1 and 0.9), isolating eps's own
  # per-coordinate effect from that separate, expected constraint.
  mdl <- .eps_mdl()
  ev <- rep(0.2, length(mdl$theta)); ev[5] <- 0.05   # P_21's own eps: tight
  out <- MMC_bounds(mdl, utils::modifyList(.eps_con(ev), list(P_low = 0.08)))
  # P_21 = 0.1, eps = 0.05 gives [0.05, 0.15]; P_low = 0.08 raises only the lower end.
  P_21_pos <- which(names(mdl$theta) == "P_21")
  expect_equal(out$theta_low[P_21_pos], 0.08)
  expect_equal(out$theta_upp[P_21_pos], 0.1 + 0.05)
})

test_that("a vector eps composes correctly with the variance guard", {
  mdl <- .eps_mdl()
  mdl$theta["x_1"] <- 0.05
  names(mdl$theta_var_ind) <- names(mdl$theta)
  mdl$theta_var_ind[1] <- 1   # treat x_1 as a variance coordinate for this test
  ev <- rep(0.2, length(mdl$theta))
  ev[1] <- 0.2   # eps alone would take the lower bound to 0.05 - 0.2 = -0.15 (<=0)
  out <- MMC_bounds(mdl, .eps_con(ev))
  expect_equal(out$theta_low[1], 0.05 * 0.01)   # variance_constraint override, not eps
})

test_that("a vector eps composes correctly with the CI union and a per-coordinate NA SE", {
  mdl <- .eps_mdl()
  se <- rep(0.05, length(mdl$theta)); se[2] <- NA_real_   # one non-finite SE
  mdl$theta_se <- se
  ev <- rep(0.1, length(mdl$theta))
  con <- .eps_con(ev); con$CI_union <- TRUE
  out <- MMC_bounds(mdl, con)
  th <- as.numeric(mdl$theta)
  # coordinate 1: eps=0.1 vs 2*se=0.1 -- tied, union has no effect either way
  expect_equal(out$theta_low[1], th[1] - 0.1)
  # coordinate 2: NA SE contributes nothing (treated as 0), eps alone applies
  expect_equal(out$theta_low[2], th[2] - 0.1)
  expect_equal(out$theta_upp[2], th[2] + 0.1)
})

test_that("MMCLRTest itself rejects a bad eps before any simulation or search", {
  skip_on_cran()
  set.seed(4141)
  y <- matrix(rnorm(150), 150, 1)
  expect_error(
    MMCLRTest(y, p = 1, k0 = 1, k1 = 2, control = list(N = 9, mc_seed = 4141, eps = c(0.1, 0.2))),
    "length 1 or length\\(theta_0\\)")
  expect_error(
    MMCLRTest(y, p = 1, k0 = 1, k1 = 2, control = list(N = 9, mc_seed = 4141, eps = 0)),
    "strictly positive")
})
