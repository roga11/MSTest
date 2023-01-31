# ------ MS-AR example ----- #
set.seed(123)
# Define DGP of MS AR process
mdl_ms2 <- list(n     = 200, 
                mu    = c(5,10),
                sigma = c(1,2),
                phi   = c(0.5),
                k     = 2,
                P     = rbind(c(0.90, 0.10),
                              c(0.10, 0.90)))

# Simulate process using simuMSAR() function
y_ms_simu <- simuMSAR(mdl_ms2)

# Set test procedure options
mmc_control = list(N = 19,
                   burnin = 100,
                   converge_check = NULL,
                   eps = 0.1,
                   CI_union = TRUE,
                   silence = FALSE,
                   threshold_stop = 0.05 + 1e-6,
                   mdl_h0_control = list(const  = TRUE, 
                                         getSE  = TRUE),
                   mdl_h1_control = list(msmu   = TRUE, 
                                         msvar  = TRUE,
                                         getSE  = TRUE,
                                         method = "EM",
                                         maxit  = 500,
                                         use_diff_init = 1),
                   type_control = list(maxit = 50))

st <- Sys.time()
mdl <- MMCLRTest(y_ms_simu$y, p = 1 , k0 = 1 , k1 = 2, mmc_control)
mdl
end <- Sys.time() - st


