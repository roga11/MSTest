# Tests for mdledit's transition-matrix block: the candidate P carried by theta must
# reach the named `P` field that every Markov-switching simulator reads, in the
# column-major orientation the package stores it in (P[i,j] = Pr(s_t = i | s_{t-1} = j)).
# Guards the MMC-LRT nuisance search, which otherwise simulates every candidate at the
# FITTED transition matrix.

set.seed(9100)
me_sim <- simuMSAR(list(n = 300, k = 2, mu = c(-1, 2), sigma = c(1, 2), phi = 0.4,
                        P = cbind(c(0.90, 0.10), c(0.15, 0.85))), burnin = 100)
me_y <- matrix(me_sim$y, ncol = 1)

test_that("mdledit writes the candidate transition matrix, column-major", {
  skip_on_cran()
  set.seed(9101)
  m <- MSARmdl(me_y, p = 1, k = 2, control = list(getSE = FALSE, use_diff_init = 1))
  ind <- which(as.numeric(m$theta_P_ind) == 1)
  # the P block is the LAST k*k entries and holds vec(P)
  expect_identical(ind, (length(m$theta) - 4L + 1L):length(m$theta))
  expect_equal(matrix(m$theta[ind], 2, 2), as.matrix(m$P), tolerance = 0)

  # a candidate P that is asymmetric (so a transposed read would be detected)
  Pc <- cbind(c(0.70, 0.30), c(0.40, 0.60))
  th <- m$theta; th[ind] <- c(Pc)
  ed <- mdledit(m, th, p = 1, q = 1, k0 = 2, exog = FALSE)
  expect_equal(as.matrix(ed$P), Pc, tolerance = 1e-12)
  expect_false(isTRUE(all.equal(as.matrix(ed$P), as.matrix(m$P))))
  # the other blocks still travel
  expect_equal(as.numeric(ed$mu),    as.numeric(th[which(as.numeric(m$theta_mu_ind)  == 1)]))
  expect_equal(as.numeric(ed$sigma), as.numeric(th[which(as.numeric(m$theta_sig_ind) == 1)]))
})

test_that("simulation from the edited model follows the candidate P, not the fitted P", {
  skip_on_cran()
  set.seed(9102)
  m <- MSARmdl(me_y, p = 1, k = 2, control = list(getSE = FALSE, use_diff_init = 1))
  ind <- which(as.numeric(m$theta_P_ind) == 1)
  Pc <- cbind(c(0.70, 0.30), c(0.40, 0.60))          # far from the fitted, persistent P
  th <- m$theta; th[ind] <- c(Pc)
  ed <- mdledit(m, th, p = 1, q = 1, k0 = 2, exog = FALSE)
  ed$n <- 20000
  set.seed(9103)
  st <- as.integer(round(as.numeric(simuMSAR(ed, burnin = 100)$St)))
  u  <- sort(unique(st))
  expect_length(u, 2L)
  Phat <- matrix(0, 2, 2)
  for (t in 2:length(st)) {
    i <- match(st[t], u); j <- match(st[t - 1], u)
    Phat[i, j] <- Phat[i, j] + 1
  }
  Phat <- sweep(Phat, 2, colSums(Phat), "/")
  expect_equal(Phat, Pc, tolerance = 0.03)           # matches the candidate
  expect_gt(max(abs(Phat - as.matrix(m$P))), 0.10)   # and NOT the fitted P
  expect_gt(max(abs(Phat - t(Pc))), 0.05)            # orientation is not transposed
})

test_that("mdledit is inert at theta_0 and for a one-regime null", {
  skip_on_cran()
  set.seed(9104)
  m <- MSARmdl(me_y, p = 1, k = 2, control = list(getSE = FALSE, use_diff_init = 1))
  # round-trip at the seed: theta stores vectorise(P), so this is exact
  ed0 <- mdledit(m, m$theta, p = 1, q = 1, k0 = 2, exog = FALSE)
  expect_identical(as.numeric(ed0$P), as.numeric(m$P))

  # k0 = 1: the null model has no transition matrix at all; mdledit must not touch it
  set.seed(9105)
  a <- ARmdl(me_y, p = 1, control = list(getSE = FALSE))
  expect_null(a$theta_P_ind)
  expect_null(a$P)
  th <- a$theta; th[1] <- th[1] + 0.25
  eda <- mdledit(a, th, p = 1, q = 1, k0 = 1, exog = FALSE)
  expect_null(eda$P)
  expect_setequal(names(eda), names(a))
  expect_equal(as.numeric(eda$mu), as.numeric(th[which(as.numeric(a$theta_mu_ind) == 1)]))
})

test_that("the edited model reproduces an independently built DGP list, all MS classes", {
  skip_on_cran()
  Tn <- 150; BI <- 100
  set.seed(9106)
  eps1 <- matrix(rnorm(Tn + BI), Tn + BI, 1)
  eps2 <- matrix(rnorm((Tn + BI) * 2), Tn + BI, 2)
  sr   <- runif(Tn + BI)
  Pc   <- cbind(c(0.65, 0.35), c(0.35, 0.65))

  # --- univariate MS-AR
  set.seed(9107)
  m <- MSARmdl(me_y, p = 1, k = 2, control = list(getSE = FALSE, use_diff_init = 1))
  th <- m$theta
  th[which(as.numeric(m$theta_mu_ind)  == 1)] <- c(-2, 3)
  th[which(as.numeric(m$theta_phi_ind) == 1)] <- 0.3
  th[which(as.numeric(m$theta_sig_ind) == 1)] <- c(0.5, 4)
  th[which(as.numeric(m$theta_P_ind)   == 1)] <- c(Pc)
  ed <- mdledit(m, th, p = 1, q = 1, k0 = 2, exog = FALSE)
  ed$n <- Tn; ed$eps <- eps1; ed$state_rand <- sr
  hand <- list(n = Tn, k = 2, q = 1, p = 1, mu = c(-2, 3), sigma = c(0.5, 4), phi = 0.3,
               P = Pc, eps = eps1, state_rand = sr)
  expect_identical(as.numeric(simuMSAR_cpp(ed, BI)$y), as.numeric(simuMSAR_cpp(hand, BI)$y))

  # --- multivariate hidden Markov (p = 0) and MS-VAR
  set.seed(9108)
  dh <- simuHMM(list(n = 300, k = 2, q = 2, mu = rbind(c(-1, -1), c(2, 2)),
                     sigma = list(diag(2), 2 * diag(2)),
                     P = cbind(c(0.9, 0.1), c(0.15, 0.85))), burnin = 100)
  mh <- HMmdl(dh$y, k = 2, control = list(getSE = FALSE, use_diff_init = 1))
  muC <- rbind(c(-2, 0), c(3, 1))
  S1  <- matrix(c(1, 0.3, 0.3, 1), 2, 2); S2 <- matrix(c(4, -0.5, -0.5, 2), 2, 2)
  th <- mh$theta
  th[which(as.numeric(mh$theta_mu_ind)  == 1)] <- c(t(muC))
  th[which(as.numeric(mh$theta_sig_ind) == 1)] <- c(S1[lower.tri(S1, diag = TRUE)],
                                                    S2[lower.tri(S2, diag = TRUE)])
  th[which(as.numeric(mh$theta_P_ind)   == 1)] <- c(Pc)
  ed <- mdledit(mh, th, p = 0, q = 2, k0 = 2, exog = FALSE)
  ed$n <- Tn; ed$eps <- eps2; ed$state_rand <- sr
  hand <- list(n = Tn, k = 2, q = 2, mu = muC, sigma = list(S1, S2), P = Pc,
               eps = eps2, state_rand = sr)
  expect_identical(as.numeric(simuHMM_cpp(ed, BI, FALSE)$y),
                   as.numeric(simuHMM_cpp(hand, BI, FALSE)$y))

  set.seed(9109)
  dv <- simuMSVAR(list(n = 300, p = 1, q = 2, k = 2, mu = rbind(c(-1, -1), c(2, 2)),
                       sigma = list(diag(2), 2 * diag(2)),
                       phi = rbind(c(0.4, 0.1), c(0.05, 0.3)),
                       P = cbind(c(0.9, 0.1), c(0.15, 0.85))), burnin = 100)
  mv <- MSVARmdl(dv$y, p = 1, k = 2, control = list(getSE = FALSE, use_diff_init = 1))
  phiC <- rbind(c(0.2, 0.05), c(0.1, 0.25))
  th <- mv$theta
  th[which(as.numeric(mv$theta_mu_ind)  == 1)] <- c(t(muC))
  th[which(as.numeric(mv$theta_phi_ind) == 1)] <- c(t(phiC))
  th[which(as.numeric(mv$theta_sig_ind) == 1)] <- c(S1[lower.tri(S1, diag = TRUE)],
                                                    S2[lower.tri(S2, diag = TRUE)])
  th[which(as.numeric(mv$theta_P_ind)   == 1)] <- c(Pc)
  ed <- mdledit(mv, th, p = 1, q = 2, k0 = 2, exog = FALSE)
  ed$n <- Tn; ed$eps <- eps2; ed$state_rand <- sr
  hand <- list(n = Tn, k = 2, q = 2, p = 1, mu = muC, sigma = list(S1, S2), phi = phiC,
               P = Pc, eps = eps2, state_rand = sr)
  expect_identical(as.numeric(simuMSVAR_cpp(ed, BI)$y), as.numeric(simuMSVAR_cpp(hand, BI)$y))
})
