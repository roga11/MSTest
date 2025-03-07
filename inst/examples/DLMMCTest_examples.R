set.seed(1234)
# load data used in Hamilton 1989 and extended data used in CHP 2014 
y84 <- as.matrix(hamilton84GNP$GNP_gr)
y10 <- as.matrix(chp10GNP$GNP_gr)

# Set test procedure options
mmc_control <- list(N = 99,
                    simdist_N = 10000,
                    getSE = TRUE,
                    eps = 0, 
                    CI_union = TRUE,
                    lambda = 100,
                    stationary_ind = TRUE,
                    optim_type = "GenSA",
                    silence = FALSE,
                    threshold_stop = 1,
                    maxit = 200)


# perform test on Hamilton 1989 data
\donttest{
  mmc_gnp84 <- DLMMCTest(y84, p = 4, control = mmc_control)
  summary(mmc_gnp84)
}

