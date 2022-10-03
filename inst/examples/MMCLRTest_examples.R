# ------ MS-AR example ----- #
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

set.seed(1234)
# load data used in Hamilton 1989 and extended data used in CHP 2014 
y84 <- as.matrix(hamilton84GNP$GNP_logdiff)

mdl <- MMCLRTest(y84, p = 4 , k0 = 1 , k1 = 2, mmc_control)
mdl


# ------ MS-VAR example ----- #

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
                                         msvar  = TRUE,
                                         getSE  = TRUE,
                                         method = "EM",
                                         maxit  = 500,
                                         use_diff_init = 1),
                   type_control = list(maxit = 50))


set.seed(1234)
# Define DGP of MS VAR process
mdl_msvar2 <- list(n     = 1000, 
                   p     = 1,
                   q     = 2,
                   mu    = rbind(c(5, -2),
                                 c(10, 2)),
                   sigma = list(rbind(c(5.0, 1.5),
                                      c(1.5, 1.0)),
                                rbind(c(7.0, 3.0),
                                      c(3.0, 2.0))),
                   phi   = rbind(c(0.50, 0.30),
                                 c(0.20, 0.70)),
                   k     = 2,
                   P     = rbind(c(0.90, 0.10),
                                 c(0.10, 0.90)))

# Simulate process using simuMSVAR() function
y_msvar_simu <- simuMSVAR(mdl_msvar2)


st <- Sys.time()
mdl <- MMCLRTest(y_msvar_simu$y, p = 1 , k0 = 1 , k1 = 2, mmc_control)
end <- Sys.time() - st
mdl
end

