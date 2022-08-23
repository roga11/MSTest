# load data used in Hamilton 1989 and extended data used in CHP 2014 
y84 <- as.matrix(hamilton84GNP$GNP_logdiff)
y10 <- as.matrix(carhuplo10GNP$GNP_logdiff)

# Set test procedure options
lmc_control = list(N = 99,
                   simdist_N = 10000,
                   getSE = TRUE)

mmc_control <- list(N = 99,
                    pval_type = "min",
                    simdist_N = 10000,
                    getSE = TRUE,
                    eps = 0.0000001,
                    CI_union = TRUE,
                    lambda = 100,
                    stationary_ind = TRUE,
                    type = "GenSA",
                    silence = FALSE,
                    threshold_stop = 1,
                    type_control = list(maxit = 200))

  
# perform test with switch in mean only on Hamilton 1989 data
lmc_gnp84 <- DLMCTest(y84, p = 4, control = lmc_control)
lmc_gnp84

st <- Sys.time()
mmc_gnp84 <- DLMMCTest(y84, p = 4, control = mmc_control)
mmc_gnp84
end <- Sys.time() - st
end


# perform test with switch in mean only on extended data used in CHP 2014
lmc_gnp10 <- DLMCTest(y10, p = 4, control = lmc_control)
lmc_gnp10




st <- Sys.time()
mmc_gnp10 <- DLMMCTest(y10, p = 4, control = mmc_control)
mmc_gnp10
end <- Sys.time() - st
end




