set.seed(1234)
# load data used in Hamilton 1989 and extended data used in CHP 2014 
y84 <- as.matrix(hamilton84GNP$GNP_logdiff)
y10 <- as.matrix(carhuplo10GNP$GNP_logdiff)

# Set test procedure options
mmc_control <- list(N = 99,
                    simdist_N = 10000,
                    getSE = TRUE,
                    eps = 0.0000001, 
                    CI_union = TRUE,
                    lambda = 100,
                    stationary_ind = TRUE,
                    optim_type = "GenSA",
                    silence = FALSE,
                    threshold_stop = 1,
                    type_control = list(maxit = 200))


# perform test on Hamilton 1989 data
mmc_gnp84 <- DLMMCTest(y84, p = 4, control = mmc_control)
mmc_gnp84

# perform test on extended data used in Carrasco, Hu & Ploberger 2014 & Dufour & Luger 2017
mmc_gnp10 <- DLMMCTest(y10, p = 4, control = mmc_control)
mmc_gnp10