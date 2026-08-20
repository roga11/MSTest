# Phase 2: the maximised Monte Carlo search perturbs the k0^2 entries of the null's
# transition matrix independently, so a candidate satisfies "every column sums to 1"
# only by a measure-zero coincidence; almost every candidate failed the colsum gate,
# so for k0 >= 2 the search evaluated essentially only theta_0 itself while still
# paying the full evaluation budget. The fix searches only the free entries of each
# column of P (k0-1 of them) and derives the remaining one as 1 - sum(free); theta
# itself, theta_P_ind, SEs, names, mdledit and every reported field stay in full
# k0^2 coordinates -- only the optimizer's search vector is reduced.

.make_theta_P <- function(k0, other_len = 3, seed = 1){
  set.seed(seed)
  other <- round(stats::runif(other_len), 3)
  P <- matrix(stats::runif(k0 * k0), k0, k0)
  P <- sweep(P, 2, colSums(P), "/")   # column-stochastic
  theta <- c(other, as.numeric(P))
  theta_P_ind <- c(rep(0, other_len), rep(1, k0 * k0))
  list(theta = theta, theta_P_ind = theta_P_ind, P = P, other_len = other_len)
}

test_that(".mmc_derived_P_index derives the largest entry of each column, one per column", {
  for (k0 in 2:5) {
    fx <- .make_theta_P(k0, seed = k0)
    rp <- MSTest:::.mmc_derived_P_index(fx$theta, fx$theta_P_ind, k0)
    expect_length(rp$derived_full_idx, k0)
    expect_length(unique(rp$derived_full_idx), k0)          # one per column, all distinct
    expect_equal(dim(rp$pos_mat), c(k0, k0))
    for (j in seq_len(k0)) {
      col_positions <- rp$pos_mat[, j]
      expect_equal(rp$derived_full_idx[j], col_positions[which.max(fx$theta[col_positions])],
                  label = paste("k0 =", k0, "col", j))
    }
  }
})

test_that("the derived choice is frozen: it does not change if theta changes after the call", {
  # D2.1: the map is computed ONCE from theta_0 and frozen for the whole search, not
  # recomputed per candidate. This test pins that the helper itself has no hidden state
  # that would make a second call see a different theta -- i.e. it is a pure function
  # of its arguments, so freezing it is the caller's responsibility (done once in
  # MMCLRTest) and is not silently defeated by the helper recomputing internally.
  k0 <- 3
  fx <- .make_theta_P(k0, seed = 42)
  rp1 <- MSTest:::.mmc_derived_P_index(fx$theta, fx$theta_P_ind, k0)
  theta_perturbed <- fx$theta
  theta_perturbed[fx$other_len + 1] <- theta_perturbed[fx$other_len + 1] + 10   # a non-P entry
  rp2 <- MSTest:::.mmc_derived_P_index(theta_perturbed, fx$theta_P_ind, k0)
  expect_identical(rp1$derived_full_idx, rp2$derived_full_idx)   # unaffected by non-P entries
})

test_that("a fixed coordinate is never chosen as the derived entry", {
  # M3 / Phase 5 compatibility: forward-looking hook, unused until a caller fixes a
  # coordinate, but the selection rule must already honor it correctly when it does.
  k0 <- 2
  fx <- .make_theta_P(k0, other_len = 0, seed = 9)
  # Force column 1's largest entry to be the one we then declare fixed.
  col1 <- MSTest:::.mmc_derived_P_index(fx$theta, fx$theta_P_ind, k0)$pos_mat[, 1]
  largest_pos <- col1[which.max(fx$theta[col1])]
  rp <- MSTest:::.mmc_derived_P_index(fx$theta, fx$theta_P_ind, k0, fixed_full_idx = largest_pos)
  expect_false(rp$derived_full_idx[1] == largest_pos)
  expect_true(rp$derived_full_idx[1] %in% setdiff(col1, largest_pos))
})

test_that("a column with every entry fixed errors instead of silently deriving nothing", {
  k0 <- 2
  fx <- .make_theta_P(k0, other_len = 0, seed = 11)
  col1 <- MSTest:::.mmc_derived_P_index(fx$theta, fx$theta_P_ind, k0)$pos_mat[, 1]
  expect_error(MSTest:::.mmc_derived_P_index(fx$theta, fx$theta_P_ind, k0, fixed_full_idx = col1),
              "fixed")
})

test_that("reassembly reproduces column sums of exactly 1, for random free vectors", {
  # The formula's own property (derived = 1 - sum(free)), independent of whether the
  # free vector came from any particular model: colSum = sum(free) + (1 - sum(free))
  # is exactly 1.0 in floating point for values in [0,1] (Sterbenz-style exactness of
  # the subtraction), which is the "reassembly is exact" property the spec measured
  # at 0.0 deviation over 200k random draws.
  set.seed(2026)
  for (k0 in c(2, 3, 4, 5, 8)) {
    fx <- .make_theta_P(k0, other_len = 2, seed = 100 + k0)
    rp <- MSTest:::.mmc_derived_P_index(fx$theta, fx$theta_P_ind, k0)
    reduced_idx <- setdiff(seq_along(fx$theta), rp$derived_full_idx)
    for (rep in 1:20) {
      theta_reduced <- fx$theta[reduced_idx]
      # randomize just the free P entries, in [0,1]
      free_P_reduced_pos <- which(reduced_idx %in% setdiff(rp$full_idx, rp$derived_full_idx))
      theta_reduced[free_P_reduced_pos] <- stats::runif(length(free_P_reduced_pos))
      theta_full <- MSTest:::.mmc_reassemble_theta(theta_reduced, reduced_idx, rp$derived_full_idx,
                                                    rp$pos_mat, length(fx$theta))
      P_full <- matrix(theta_full[rp$full_idx], k0, k0)
      expect_equal(colSums(P_full), rep(1, k0), tolerance = 0, label = paste("k0=", k0, "rep", rep))
    }
  }
})

test_that("round trip (full -> reduced -> full) reproduces a valid theta_0", {
  set.seed(303)
  for (k0 in 2:5) {
    fx <- .make_theta_P(k0, other_len = 3, seed = 200 + k0)
    rp <- MSTest:::.mmc_derived_P_index(fx$theta, fx$theta_P_ind, k0)
    reduced_idx <- setdiff(seq_along(fx$theta), rp$derived_full_idx)
    theta_full2 <- MSTest:::.mmc_reassemble_theta(fx$theta[reduced_idx], reduced_idx,
                                                  rp$derived_full_idx, rp$pos_mat, length(fx$theta))
    expect_equal(theta_full2, fx$theta, tolerance = 1e-10)
  }
})

test_that("the reduced search vector has exactly npar - k0 entries", {
  for (k0 in 2:5) {
    fx <- .make_theta_P(k0, other_len = 4, seed = 300 + k0)
    rp <- MSTest:::.mmc_derived_P_index(fx$theta, fx$theta_P_ind, k0)
    reduced_idx <- setdiff(seq_along(fx$theta), rp$derived_full_idx)
    expect_length(reduced_idx, length(fx$theta) - k0)
  }
})

test_that(".mmc_reduced_wrap calls the wrapped objective only when derived entries are admissible", {
  k0 <- 3
  fx <- .make_theta_P(k0, other_len = 2, seed = 555)
  rp <- MSTest:::.mmc_derived_P_index(fx$theta, fx$theta_P_ind, k0)
  reduced_idx <- setdiff(seq_along(fx$theta), rp$derived_full_idx)
  npar <- length(fx$theta)
  theta_low <- fx$theta - 0.3
  theta_upp <- fx$theta + 0.3
  calls <- 0
  mock_obj <- function(theta_full, ...) { calls <<- calls + 1; sum(theta_full) }
  wrapped <- MSTest:::.mmc_reduced_wrap(mock_obj, reduced_idx, rp$derived_full_idx, rp$pos_mat, npar,
                                        theta_low, theta_upp, P_low = 0, P_upp = 1, penalty_value = -999)
  # Feasible: the reduced vector taken directly from theta_0 reassembles to theta_0
  # itself, which is inside every one of its own bounds -> the wrapped objective runs.
  out_ok <- wrapped(fx$theta[reduced_idx])
  expect_equal(calls, 1L)
  expect_equal(out_ok, sum(fx$theta))

  # Infeasible via [P_low, P_upp]: push a free entry in column 1 far up, forcing that
  # column's derived value below P_low = 0 -- the penalty returns without calling the mock.
  bad_reduced <- fx$theta[reduced_idx]
  free_pos_col1 <- which(reduced_idx %in% setdiff(rp$pos_mat[, 1], rp$derived_full_idx[1]))
  bad_reduced[free_pos_col1[1]] <- 5
  out_bad <- wrapped(bad_reduced)
  expect_equal(calls, 1L)   # unchanged: mock NOT called again
  expect_equal(out_bad, -999)
})

test_that("M2: a derived entry outside its own eps/CI box is penalized even inside [P_low, P_upp]", {
  k0 <- 2
  fx <- .make_theta_P(k0, other_len = 2, seed = 777)
  rp <- MSTest:::.mmc_derived_P_index(fx$theta, fx$theta_P_ind, k0)
  reduced_idx <- setdiff(seq_along(fx$theta), rp$derived_full_idx)
  npar <- length(fx$theta)
  # A box pinned exactly at theta_0: any perturbation of a free P entry moves the
  # derived value away from its own tight bound, even though [0,1] admits it.
  theta_low <- fx$theta
  theta_upp <- fx$theta
  calls <- 0
  mock_obj <- function(theta_full, ...) { calls <<- calls + 1; 0 }
  wrapped <- MSTest:::.mmc_reduced_wrap(mock_obj, reduced_idx, rp$derived_full_idx, rp$pos_mat, npar,
                                        theta_low, theta_upp, P_low = 0, P_upp = 1, penalty_value = -777)
  bad_reduced <- fx$theta[reduced_idx]
  free_pos_col1 <- which(reduced_idx %in% setdiff(rp$pos_mat[, 1], rp$derived_full_idx[1]))
  bad_reduced[free_pos_col1[1]] <- min(1, bad_reduced[free_pos_col1[1]] + 0.01)
  out <- wrapped(bad_reduced)
  expect_equal(out, -777)
  expect_equal(calls, 0L)
})

# ---------------------------------------------------------------------------
# Integration: the real MMCLRTest, k0 = 1 unaffected, k0 >= 2 now finds
# admissible candidates.

test_that("k0 = 1 is bit-identical to a build with none of tonight's changes", {
  skip_on_cran()
  skip_if_not(nzchar(Sys.getenv("MSTEST_BASELINE_LIB")), "no baseline library configured")
  base_lib <- Sys.getenv("MSTEST_BASELINE_LIB")
  script <- tempfile(fileext = ".R")
  writeLines(c(
    sprintf('.libPaths(c("%s", .libPaths()))', base_lib),
    'library(MSTest)',
    'set.seed(31415)',
    's <- simuMSAR(list(n = 200, k = 2, mu = c(-1, 1.5), sigma = c(1, 1), phi = 0.3,',
    '                    P = cbind(c(.9,.1), c(.15,.85))), burnin = 100)',
    'y <- matrix(s$y, ncol = 1)',
    'out <- MMCLRTest(y, p = 1, k0 = 1, k1 = 2, control = list(N = 9, mc_seed = 707,',
    '                  maxit = 8, mdl_h0_control = list(getSE = FALSE),',
    '                  mdl_h1_control = list(getSE = FALSE, use_diff_init = 1)))',
    'saveRDS(list(pval = out$pval, pval_0 = out$pval_0, theta_h0 = as.numeric(out$theta_h0),',
    '             LRT_0 = out$LRT_0), commandArgs(trailingOnly = TRUE)[1])'
  ), script)
  out_base <- tempfile(fileext = ".rds")
  system2("Rscript", c(shQuote(script), shQuote(out_base)))
  base_res <- readRDS(out_base)

  set.seed(31415)
  s <- simuMSAR(list(n = 200, k = 2, mu = c(-1, 1.5), sigma = c(1, 1), phi = 0.3,
                     P = cbind(c(.9,.1), c(.15,.85))), burnin = 100)
  y <- matrix(s$y, ncol = 1)
  out <- MMCLRTest(y, p = 1, k0 = 1, k1 = 2, control = list(N = 9, mc_seed = 707,
                   maxit = 8, mdl_h0_control = list(getSE = FALSE),
                   mdl_h1_control = list(getSE = FALSE, use_diff_init = 1)))
  expect_identical(out$pval, base_res$pval)
  expect_identical(out$pval_0, base_res$pval_0)
  expect_identical(as.numeric(out$theta_h0), base_res$theta_h0)
  expect_identical(out$LRT_0, base_res$LRT_0)
})

test_that("k0 >= 2 now reaches candidates that pass the transition-matrix gate", {
  # Before Phase 2, essentially 0% of candidates satisfied the colsum constraint by
  # chance (measured 0 of 16,000-200,000 uniform draws in the spec's own review), so
  # the search evaluated only theta_0 itself. After the fix, a search vector's
  # reassembly always sums to 1 by construction, so every candidate that also clears
  # [P_low, P_upp] and its own eps box is admissible -- the gate should no longer be
  # the binding constraint. Evidence, not a fixed threshold: a wider search finds a
  # pval that is not simply the theta_0 evaluation, for at least one of a few seeds.
  skip_on_cran()
  found_improvement <- FALSE
  for (sd in c(101, 202)) {
    set.seed(sd)
    s <- simuMSAR(list(n = 150, k = 2, mu = c(-1.5, 1.5), sigma = c(1, 1.5), phi = 0.2,
                       P = cbind(c(.85, .15), c(.2, .8))), burnin = 100)
    y <- matrix(s$y, ncol = 1)
    out <- MMCLRTest(y, p = 1, k0 = 2, k1 = 3, control = list(
      N = 9, mc_seed = sd, maxit = 20, eps = 0.15, type = "pso", silence = TRUE, CI_union = FALSE,
      mdl_h0_control = list(getSE = FALSE, use_diff_init = 2),
      mdl_h1_control = list(getSE = FALSE, use_diff_init = 2)))
    expect_true(out$pval >= 0 && out$pval <= 1, label = paste("seed", sd))
    expect_gte(out$pval, out$pval_0)
    P_h0 <- matrix(as.numeric(out$theta_h0)[as.numeric(out$mdl_h0$theta_P_ind) == 1], 2, 2)
    expect_equal(colSums(P_h0), c(1, 1), tolerance = 1e-8, label = paste("seed", sd))
    if (out$pval > out$pval_0 + 1e-9) found_improvement <- TRUE
  }
  expect_true(found_improvement)
})

test_that("k0 >= 2 works under GenSA and GA too, and pval >= pval_0 in every case", {
  skip_on_cran()
  set.seed(9911)
  s <- simuMSAR(list(n = 150, k = 2, mu = c(-1, 1.2), sigma = c(1, 1), phi = 0.25,
                     P = cbind(c(.88, .12), c(.18, .82))), burnin = 100)
  y <- matrix(s$y, ncol = 1)
  for (ty in c("GenSA", "GA")) {
    out <- MMCLRTest(y, p = 1, k0 = 2, k1 = 3, control = list(
      N = 9, mc_seed = 9911, maxit = 8, eps = 0.1, type = ty, silence = TRUE, CI_union = FALSE,
      optim_control = if (ty == "GA") list(popSize = 6) else list(),
      mdl_h0_control = list(getSE = FALSE, use_diff_init = 2),
      mdl_h1_control = list(getSE = FALSE, use_diff_init = 2)))
    expect_true(out$pval >= 0 && out$pval <= 1, label = ty)
    expect_gte(out$pval, out$pval_0)
    P_h0 <- matrix(as.numeric(out$theta_h0)[as.numeric(out$mdl_h0$theta_P_ind) == 1], 2, 2)
    expect_equal(colSums(P_h0), c(1, 1), tolerance = 1e-8, label = ty)
  }
})

test_that("threshold_stop early-stop at k0 >= 2 still returns a full-coordinate theta_h0", {
  skip_on_cran()
  set.seed(1414)
  s <- simuMSAR(list(n = 150, k = 2, mu = c(-1, 1.2), sigma = c(1, 1), phi = 0.25,
                     P = cbind(c(.88, .12), c(.18, .82))), burnin = 100)
  y <- matrix(s$y, ncol = 1)
  out <- MMCLRTest(y, p = 1, k0 = 2, k1 = 3, control = list(
    N = 9, mc_seed = 1414, maxit = 5, threshold_stop = 1e-9, silence = TRUE, CI_union = FALSE,
    mdl_h0_control = list(getSE = FALSE, use_diff_init = 2),
    mdl_h1_control = list(getSE = FALSE, use_diff_init = 2)))
  expect_true(isTRUE(out$mmc_optimout$early_stop))
  expect_length(as.numeric(out$theta_h0), length(out$mdl_h0$theta))
  P_h0 <- matrix(as.numeric(out$theta_h0)[as.numeric(out$mdl_h0$theta_P_ind) == 1], 2, 2)
  expect_equal(colSums(P_h0), c(1, 1), tolerance = 1e-8)
})

test_that("a near-absorbing candidate P reaches the simulator without error (limP exposure)", {
  # Before Phase 2, a candidate P this close to absorbing almost never passed the
  # colsum gate, so this path through limP() was essentially unreachable from MMC.
  skip_on_cran()
  set.seed(5252)
  s <- simuMSAR(list(n = 150, k = 2, mu = c(-1, 1.2), sigma = c(1, 1), phi = 0.2,
                     P = cbind(c(0.999, 0.001), c(0.05, 0.95))), burnin = 100)
  y <- matrix(s$y, ncol = 1)
  out <- expect_no_error(MMCLRTest(y, p = 1, k0 = 2, k1 = 3, control = list(
    N = 9, mc_seed = 5252, maxit = 10, eps = 0.05, silence = TRUE, CI_union = FALSE,
    mdl_h0_control = list(getSE = FALSE, use_diff_init = 2),
    mdl_h1_control = list(getSE = FALSE, use_diff_init = 2))))
  expect_true(is.finite(out$pval))
})
