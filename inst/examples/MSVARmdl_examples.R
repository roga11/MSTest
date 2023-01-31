set.seed(123)
# Define DGP of MS VAR process
mdl_msvar2 <- list(n     = 200, 
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

# Set options for model estimation
control <- list(msmu   = TRUE, 
                msvar  = TRUE,
                method = "EM",
                use_diff_init = 1)
                
# Estimate model
y_msvar_mdl <- MSVARmdl(y_msvar_simu$y, p = 1, k = 2, control)
y_msvar_mdl
