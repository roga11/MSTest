set.seed(1234)
# Define DGP of MS AR process
mdl_ms2 <- list(n     = 200, 
                mu    = c(5,10),
                sigma = c(1,4),
                phi   = c(0.5),
                k     = 2,
                P     = rbind(c(0.90, 0.10),
                              c(0.10, 0.90)))

# Simulate process using simuMSAR() function
y_ms_simu <- simuMSAR(mdl_ms2)


# ------ MS-AR example ----- #
# Set test procedure options
lmc_control = list(N = 19,
                   burnin = 100,
                   converge_check = NULL,
                   mdl_h0_control = list(const  = TRUE, 
                                         getSE  = TRUE),
                   mdl_h1_control = list(msmu   = TRUE, 
                                         msvar  = TRUE,
                                         getSE  = TRUE,
                                         method = "EM",
                                         maxit  = 300,
                                         use_diff_init = 1))


#mdl <- LMCLRTest(y_ms_simu$y, p = 1, k0 = 1 , k1 = 2, lmc_control)
#mdl
