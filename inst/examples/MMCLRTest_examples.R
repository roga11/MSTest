
set.seed(1234)
# ----- Univariate ----- # 
# Define DGP 
mdl_norm <- list(n     = 1000, 
                 q     = 1,
                 mu    = as.matrix(5),
                 sigma = as.matrix(5.0))

# Simulate process using simuNorm() function
y_norm_simu <- simuNorm(mdl_norm)


# Set test procedure options
mmc_control = list(N = 19,
                   burnin = 100,
                   converge_check = NULL,
                   eps = 0.1,
                   CI_union = TRUE,
                   silence = FALSE,
                   threshold_stop = 0.05 + 1e-6,
                   type = "GenSA",
                   mdl_h0_control = list(const  = TRUE, 
                                         getSE  = TRUE),
                   mdl_h1_control = list(msmu   = TRUE, 
                                         msvar  = TRUE,
                                         getSE  = TRUE,
                                         method = "EM",
                                         maxit  = 300,
                                         use_diff_init = 1),
                   type_control = list(maxit = 10))

#mdl <- MMCLRTest(y_ms_simu$y, p = 0 , k0 = 1 , k1 = 2, mmc_control)
#mdl
