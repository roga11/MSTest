library(foreach)
library(doParallel)


set.seed(1234)
# load data used in Hamilton 1989 and extended data used in CHP 2014 
y84 <- as.matrix(hamilton84GNP$GNP_logdiff)
y10 <- as.matrix(carhuplo10GNP$GNP_logdiff)

# Set test procedure options
lmc_control = list(N = 99,
                   burnin = 100,
                   converge_check = NULL,
                   dist_converge_iter = 1,
                   workers = 11)

mdl_h0_control <- list(const  = TRUE, 
                       getSE  = FALSE)


mdl_h1_control <- list(msmu   = TRUE, 
                       msvar  = FALSE,
                       getSE  = FALSE,
                       method = "MLE",
                       maxit  = 500,
                       use_diff_init = 3)


cl <- makeCluster(11)
registerDoParallel(cl)
st <- Sys.time()
mdl <- LMCLRTest(y84, p = 4 , k0 = 1 , k1 = 2, lmc_control, mdl_h0_control, mdl_h1_control)
end <- Sys.time() - st
stopCluster(cl)
mdl
end

