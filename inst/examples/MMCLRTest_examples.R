set.seed(1234)
# load data used in Hamilton 1989 and extended data used in CHP 2014 
y84 <- as.matrix(hamilton84GNP$GNP_logdiff)
y10 <- as.matrix(carhuplo10GNP$GNP_logdiff)




# Set test procedure options
mmc_control = list(N = 19,
                   burnin = 100,
                   converge_check = NULL,
                   eps = 0.000000001,
                   CI_union = TRUE,
                   silence = FALSE,
                   threshold_stop = 1,
                   mdl_h0_control = list(const  = TRUE, 
                                         getSE  = TRUE),
                   mdl_h1_control = list(msmu   = TRUE, 
                                         msvar  = FALSE,
                                         getSE  = TRUE,
                                         method = "EM",
                                         maxit  = 500,
                                         use_diff_init = 1),
                   type_control = list(maxit = 50))




st <- Sys.time()
mdl <- MMCLRTest(y84, p = 4 , k0 = 1 , k1 = 2, mmc_control)
end <- Sys.time() - st
mdl
end

