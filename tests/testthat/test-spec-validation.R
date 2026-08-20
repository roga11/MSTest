# A model specification has to describe a model that can be estimated. Both checks
# below reject a specification rather than let it reach the estimator, where the
# first fails with a message naming no parameter and the second, for the classes
# that read the transition block from the end of the parameter vector, is read as a
# different parameterization without complaint.

test_that("a switching model needs at least one regime-dependent parameter", {
  set.seed(770)
  y1 <- matrix(rnorm(200), 200, 1)
  y2 <- matrix(rnorm(400), 200, 2)
  Z1 <- matrix(rnorm(200), 200, 1)
  msg <- "regime-dependent parameter"
  none <- list(msmu = FALSE, msvar = FALSE, getSE = FALSE)
  expect_error(HMmdl(y1, k = 2, control = none), msg)
  expect_error(MSARmdl(y1, p = 1, k = 2, control = none), msg)
  expect_error(MSARXmdl(y1, p = 1, k = 2, Z = Z1, control = none), msg)
  expect_error(MSVARmdl(y2, p = 1, k = 2, control = none), msg)
  expect_error(MSVARXmdl(y2, p = 1, k = 2, Z = Z1, control = none), msg)
  # numeric flags are read as the estimator reads them, so these are the same case
  expect_error(HMmdl(y1, k = 3, control = list(msmu = 0, msvar = 0, getSE = FALSE)), msg)
  # one regime-dependent parameter is enough, and one regime needs none
  expect_silent(HMmdl(y1, k = 2, control = list(msmu = FALSE, msvar = TRUE,
                                                getSE = FALSE, use_diff_init = 1)))
  expect_silent(Nmdl(y1))
})

test_that("init_theta must have the length its model implies", {
  set.seed(771)
  y1 <- matrix(rnorm(200), 200, 1)
  y2 <- matrix(rnorm(400), 200, 2)
  ct <- list(msmu = FALSE, msvar = TRUE, getSE = FALSE)
  # The transition block is kept valid throughout, so only the length is wrong and
  # the pre-existing transition-matrix check cannot be what fires.
  vecP <- c(0.9, 0.1, 0.2, 0.8)
  # MSAR(1), k = 2, common mean: mu(1) + phi(1) + sigma(2) + vec(P)(4) = 8
  th9 <- c(rep(0.3, 5), vecP)
  expect_error(MSARmdl(y1, p = 1, k = 2, control = c(ct, list(init_theta = th9))),
               "8 parameters")
  # declaring a common mean while supplying one per regime is a mismatch, not a hint
  expect_error(MSARmdl(y1, p = 1, k = 2, control = c(ct, list(init_theta = th9))),
               "msmu = FALSE")
  # MSVAR(1), q = 2, k = 2, common mean: 2 + 4 + 3*2 + 4 = 16
  expect_error(MSVARmdl(y2, p = 1, k = 2, control = c(ct, list(init_theta = c(rep(0.3, 11), vecP)))),
               "16 parameters")
  # a correct length still reaches the existing transition-matrix validation
  bad_P <- c(0, 0.5, 1, 1, 0.9, 0.1, 0.3, 0.3)
  expect_error(MSARmdl(y1, p = 1, k = 2, control = c(ct, list(init_theta = bad_P))),
               "sum to 1")
})

test_that("the alternative must have more regimes than the null", {
  # Fires before any estimation, so garbage data is fine: neither function reaches
  # estimMdl() for k1 <= k0.
  y1 <- matrix(rnorm(50), 50, 1)
  msg <- "k1 must be greater than k0"
  expect_error(LMCLRTest(y1, p = 1, k0 = 2, k1 = 2), msg)
  expect_error(LMCLRTest(y1, p = 1, k0 = 2, k1 = 1), msg)
  expect_error(MMCLRTest(y1, p = 1, k0 = 2, k1 = 2), msg)
  expect_error(MMCLRTest(y1, p = 1, k0 = 2, k1 = 1), msg)
})
