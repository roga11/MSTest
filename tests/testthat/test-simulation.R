# Tests for the simulation functions: they should return objects with the
# documented structure and the requested dimensions.

test_that("simuNorm returns a series of the requested length", {
  set.seed(101)
  y <- simuNorm(list(n = 200, q = 1, mu = 0, sigma = as.matrix(1)))
  expect_type(y, "list")
  expect_true("y" %in% names(y))
  expect_equal(nrow(as.matrix(y$y)), 200)
  expect_true(all(is.finite(y$y)))
})

test_that("simuAR returns a univariate AR series", {
  set.seed(102)
  y <- simuAR(list(n = 300, mu = 5, sigma = 1, phi = c(0.5), p = 1))
  expect_true("y" %in% names(y))
  expect_equal(nrow(as.matrix(y$y)), 300)
  expect_true(all(is.finite(y$y)))
})

test_that("simuMSAR returns series and true regime states", {
  set.seed(103)
  mdl <- list(n = 300, mu = c(5, 10), sigma = c(1, 2), phi = c(0.5),
              k = 2, P = rbind(c(0.9, 0.1), c(0.1, 0.9)))
  y <- simuMSAR(mdl)
  expect_true(all(c("y", "St") %in% names(y)))
  expect_equal(nrow(as.matrix(y$y)), 300)
  # a true regime path is returned, one entry per simulated observation
  expect_equal(nrow(as.matrix(y$St)), 300)
})

test_that("simuMSVAR returns a multivariate series", {
  set.seed(104)
  mdl <- list(n = 200, p = 1, q = 2, k = 2,
              mu = rbind(c(5, -2), c(10, 2)),
              sigma = list(rbind(c(5.0, 1.5), c(1.5, 1.0)),
                           rbind(c(7.0, 3.0), c(3.0, 2.0))),
              phi = rbind(c(0.5, 0.3), c(0.2, 0.7)),
              P = rbind(c(0.9, 0.1), c(0.1, 0.9)))
  y <- simuMSVAR(mdl)
  expect_equal(ncol(as.matrix(y$y)), 2)
  expect_equal(nrow(as.matrix(y$y)), 200)
})
