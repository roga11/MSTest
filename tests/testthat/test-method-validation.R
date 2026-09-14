# control$method must be exactly "EM" or "MLE" on the Markov-switching
# constructors; anything else errors by name before any estimation.

test_that("an invalid 'method' errors by constructor name, before estimation", {
  set.seed(3); y <- as.matrix(rnorm(60)); Y2 <- cbind(y, as.matrix(rnorm(60))); Z <- as.matrix(rnorm(60))
  expect_error(HMmdl(y, k = 2, control = list(method = "mle")),        "HMmdl: 'method' must be \"EM\" or \"MLE\" \\(received \"mle\"\\)")
  expect_error(MSARmdl(y, p = 1, k = 2, control = list(method = "em")), "MSARmdl: 'method' must be")
  expect_error(MSARXmdl(y, p = 1, k = 2, Z = Z, control = list(method = "Mle")), "MSARXmdl: 'method' must be")
  expect_error(MSVARmdl(Y2, p = 1, k = 2, control = list(method = "ML")), "MSVARmdl: 'method' must be")
  expect_error(MSVARXmdl(Y2, p = 1, k = 2, Z = Z, control = list(method = NULL)), "MSVARXmdl: 'method' must be .*received NULL")
  expect_error(MSARmdl(y, p = 1, k = 2, control = list(method = c("EM", "MLE"))), "MSARmdl: 'method' must be")
  expect_error(MSARmdl(y, p = 1, k = 2, control = list(method = NA_character_)), "MSARmdl: 'method' must be")
  # through the dispatcher too
  expect_error(estimMdl(y, p = 1, q = 1, k = 2, Z = NULL, control = list(method = "mle")), "MSARmdl: 'method' must be")
})

test_that("valid 'method' values are untouched: the explicit default reproduces the implicit one bitwise", {
  skip_on_cran()
  set.seed(4); y <- as.matrix(rnorm(120) + rep(c(0, 2), each = 60))
  set.seed(5); a <- suppressWarnings(MSARmdl(y, p = 1, k = 2, control = list(getSE = FALSE, use_diff_init = 2)))
  set.seed(5); b <- suppressWarnings(MSARmdl(y, p = 1, k = 2, control = list(getSE = FALSE, use_diff_init = 2, method = "EM")))
  expect_identical(a$theta, b$theta); expect_identical(a$logLike, b$logLike)
  m <- suppressWarnings(MSARmdl(y, p = 1, k = 2, control = list(getSE = FALSE, use_diff_init = 1, method = "MLE")))
  expect_true(is.finite(m$logLike))
  # string-equivalent representations that the previous dispatch accepted keep working
  set.seed(5); f <- suppressWarnings(MSARmdl(y, p = 1, k = 2, control = list(getSE = FALSE, use_diff_init = 2, method = factor("EM"))))
  expect_identical(f$theta, a$theta); expect_identical(f$control$method, "EM")
})
