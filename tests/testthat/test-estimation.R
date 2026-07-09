# Tests for the model estimation functions. These check that estimation
# returns a finite log-likelihood, a valid (column-stochastic) transition
# matrix for Markov switching models, and the documented output fields.

test_that("ARmdl estimates a linear AR model", {
  set.seed(201)
  y <- simuAR(list(n = 300, mu = 5, sigma = 1, phi = c(0.5), p = 1))$y
  mdl <- ARmdl(y, p = 1)
  expect_s3_class(mdl, "ARmdl")
  expect_true(is.finite(mdl$logLike))
  expect_length(as.numeric(mdl$phi), 1)
})

test_that("MSARmdl estimates a 2-regime MS-AR model with valid P", {
  set.seed(202)
  mdl_dgp <- list(n = 400, mu = c(5, 10), sigma = c(1, 2), phi = c(0.5),
                  k = 2, P = rbind(c(0.9, 0.1), c(0.1, 0.9)))
  y <- simuMSAR(mdl_dgp)$y
  mdl <- MSARmdl(y, p = 1, k = 2,
                 control = list(msmu = TRUE, msvar = TRUE, use_diff_init = 5))
  expect_s3_class(mdl, "MSARmdl")
  expect_true(is.finite(mdl$logLike))
  # transition matrix columns sum to 1 (column-stochastic convention)
  expect_equal(unname(colSums(mdl$P)), rep(1, 2), tolerance = 1e-6)
  # smoothed regime probabilities are a T x k matrix in [0, 1]
  expect_true(all(mdl$St >= -1e-8 & mdl$St <= 1 + 1e-8))
})

test_that("HMmdl estimates a 2-regime hidden Markov model", {
  set.seed(203)
  mdl_dgp <- list(n = 400, q = 1, mu = as.matrix(c(0, 5)),
                  sigma = list(as.matrix(1), as.matrix(1)), k = 2,
                  P = rbind(c(0.9, 0.1), c(0.1, 0.9)))
  y <- simuHMM(mdl_dgp)$y
  mdl <- HMmdl(y, k = 2, control = list(msmu = TRUE, msvar = FALSE, use_diff_init = 5))
  expect_s3_class(mdl, "HMmdl")
  expect_true(is.finite(mdl$logLike))
  expect_equal(unname(colSums(mdl$P)), rep(1, 2), tolerance = 1e-6)
})

test_that("S3 methods are available for fitted models", {
  set.seed(204)
  y <- simuAR(list(n = 200, mu = 5, sigma = 1, phi = c(0.5), p = 1))$y
  mdl <- ARmdl(y, p = 1)
  expect_true(is.finite(AIC(mdl)))
  expect_true(is.finite(BIC(mdl)))
  expect_length(as.numeric(fitted(mdl)), length(as.numeric(residuals(mdl))))
})
