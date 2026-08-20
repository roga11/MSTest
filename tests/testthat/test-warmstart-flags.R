# The alternative model of a Monte Carlo LR test inherits a switching restriction
# stated for the null, so the test varies the regime count alone. Without this the
# alternative was fitted with a switching mean while its warm start was built for a
# common one, the warm start was discarded for the observed data only, and the
# observed statistic was then computed by a different rule from the simulated ones.

test_that("a restriction stated only for the null is inherited, with a warning", {
  h0 <- list(msmu = FALSE, msvar = TRUE)
  expect_warning(out <- MSTest:::.inherit_ms_flags(list(), h0), "applied to the alternative")
  expect_false(out$msmu)
  expect_false("msvar" %in% names(out))          # not restricted, so not inherited
})

test_that("a flag stated for the alternative is never overridden", {
  h0 <- list(msmu = FALSE, msvar = FALSE)
  expect_silent(out <- MSTest:::.inherit_ms_flags(list(msmu = TRUE, msvar = TRUE), h0))
  expect_true(out$msmu); expect_true(out$msvar)
})

test_that("nothing is inherited from an unrestricted or a linear null", {
  expect_silent(a <- MSTest:::.inherit_ms_flags(list(), list(msmu = TRUE, msvar = TRUE)))
  expect_false("msmu" %in% names(a))
  expect_silent(b <- MSTest:::.inherit_ms_flags(list(), list()))   # linear null
  expect_false("msmu" %in% names(b))
})

test_that("flags are matched by exact name", {
  # R's `$` partial-matches, so a mistyped control name would otherwise be read as
  # the flag and the restriction would be silently skipped.
  h0 <- list(msmu = FALSE, msvar = TRUE)
  expect_warning(out <- MSTest:::.inherit_ms_flags(list(msmuTYPO = TRUE), h0),
                 "applied to the alternative")
  expect_false(out$msmu)
  expect_true(out$msmuTYPO)
})

test_that("the alternative is estimated with the inherited restriction", {
  skip_on_cran()
  set.seed(6006)
  Tn <- 200
  S <- simuHMM(list(n = Tn, q = 2, k = 2, mu = rbind(c(2, 5), c(2, 5)),
                    sigma = list(diag(2), 3 * diag(2)),
                    P = cbind(c(.9, .1), c(.15, .85))))
  base   <- list(N = 19, mc_seed = 606)
  mdl_ct <- list(getSE = FALSE, use_diff_init = 2, maxit = 200)
  # Collect warnings rather than matching one: the call also emits unrelated ones
  # about control names, and only the inheritance warning is under test here.
  catch <- function(expr){
    seen <- character()
    val  <- withCallingHandlers(expr,
              warning = function(cond){ seen <<- c(seen, conditionMessage(cond))
                                        invokeRestart("muffleWarning") })
    list(value = val, warnings = seen)
  }
  got <- catch(LMCLRTest(S$y, p = 0, k0 = 2, k1 = 3,
                         control = c(base, list(mdl_h0_control = c(mdl_ct, list(msmu = FALSE, msvar = TRUE)),
                                                mdl_h1_control = mdl_ct))))
  tst <- got$value
  expect_true(any(grepl("applied to the alternative", got$warnings)))
  expect_false(tst$mdl_h1$msmu)          # the restriction reached the fitted model
  expect_true(is.finite(tst$LRT_0))
  expect_true(tst$pval >= 0 && tst$pval <= 1)

  # stating it on both sides is the same test, and says nothing about inheriting.
  # Only that warning is asserted absent: the call emits unrelated ones about
  # control names, which are a separate, pre-existing matter.
  set.seed(6006)
  got2 <- catch(LMCLRTest(S$y, p = 0, k0 = 2, k1 = 3,
                          control = c(base, list(mdl_h0_control = c(mdl_ct, list(msmu = FALSE, msvar = TRUE)),
                                                 mdl_h1_control = c(mdl_ct, list(msmu = FALSE, msvar = TRUE))))))
  expect_false(any(grepl("applied to the alternative", got2$warnings)))
  expect_equal(tst$LRT_0, got2$value$LRT_0)
})

test_that("warm start length matches the alternative under a restricted null", {
  # Pins warmstart_theta's own length contract for a restricted null with more than
  # one regime, which the packaged tests cover only for linear nulls. Both flags are
  # supplied here, so this does not exercise the resolution rule; that is tested
  # above.
  skip_on_cran()
  set.seed(6007)
  Tn <- 200
  S <- simuHMM(list(n = Tn, q = 2, k = 2, mu = rbind(c(1, 4), c(1, 4)),
                    sigma = list(diag(2), 2.5 * diag(2)),
                    P = cbind(c(.9, .1), c(.2, .8))))
  m0 <- HMmdl(S$y, k = 2, control = list(msmu = FALSE, msvar = TRUE, getSE = FALSE,
                                         use_diff_init = 2))
  m1 <- HMmdl(S$y, k = 3, control = list(msmu = FALSE, msvar = TRUE, getSE = FALSE,
                                         use_diff_init = 2))
  ws <- warmstart_theta(m0, 3, msmu = FALSE, msvar = TRUE)
  expect_gt(length(ws), 0)
  for (w in ws) expect_equal(length(w), length(m1$theta))
})

test_that("warm starts of the wrong length are reported, not silently dropped", {
  # The assertion that makes the fix durable. The flags are only a proxy; what has
  # to hold is that a warm start has the length the alternative expects. Simulated
  # here by making the warm-start builder return a vector for a different parameter
  # space, which is the shape any future drift in the resolution rule would take.
  skip_on_cran()
  skip_if_not_installed("testthat", "3.2.0")
  set.seed(6008)
  Tn <- 200
  S <- simuHMM(list(n = Tn, q = 2, k = 2, mu = rbind(c(0, 3), c(2, -1)),
                    sigma = list(diag(2), 2 * diag(2)),
                    P = cbind(c(.9, .1), c(.15, .85))))
  mdl_ct <- list(getSE = FALSE, use_diff_init = 1, maxit = 100)
  local_mocked_bindings(warmstart_theta = function(...) list(rep(0.5, 3)),
                        .package = "MSTest")
  seen <- character()
  withCallingHandlers(
    LMCLRTest(S$y, p = 0, k0 = 2, k1 = 3,
              control = list(N = 9, mc_seed = 608,
                             mdl_h0_control = mdl_ct, mdl_h1_control = mdl_ct)),
    warning = function(cond){ seen <<- c(seen, conditionMessage(cond))
                              invokeRestart("muffleWarning") })
  expect_true(any(grepl("did not match the alternative", seen)))
})

test_that("the embedded warm start reproduces the null fit's log-likelihood", {
  # The property the LR >= 0 policy rests on, at more than one regime with a
  # restricted flag: the packaged coverage reaches only linear nulls, where the
  # length guard cannot fire and the embedding cannot be wrong.
  skip_on_cran()
  set.seed(6009)
  Tn <- 220
  S <- simuHMM(list(n = Tn, q = 2, k = 2, mu = rbind(c(1, 4), c(1, 4)),
                    sigma = list(diag(2), 2.5 * diag(2)),
                    P = cbind(c(.9, .1), c(.2, .8))))
  for (flags in list(c(msmu = FALSE, msvar = TRUE), c(msmu = TRUE, msvar = FALSE))) {
    m0 <- HMmdl(S$y, k = 2, control = list(msmu = unname(flags["msmu"]),
                                           msvar = unname(flags["msvar"]),
                                           getSE = FALSE, use_diff_init = 2))
    for (k1 in 3:4) {
      emb <- MSTest:::.embed_theta(m0, k1, unname(flags["msmu"]), unname(flags["msvar"]))
      mdl <- list(y = m0$y, q = 2L, msmu = unname(flags["msmu"]),
                  msvar = unname(flags["msvar"]), exog = FALSE)
      ll  <- MSTest:::logLike_HMmdl(matrix(as.numeric(emb), ncol = 1), mdl, as.integer(k1))
      expect_lt(abs(ll - m0$logLike), 1e-8)
    }
  }
})
