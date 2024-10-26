# --------------------------- Use simulated process ----------------------------
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

# Set options for model estimation
control <- list(msmu   = TRUE, 
                msvar  = TRUE, 
                method = "EM",
                use_diff_init = 1)

# Estimate model
\donttest{
  ms_mdl <- MSARmdl(y_ms_simu$y, p = 1, k = 2, control = control)
  summary(ms_mdl)
}



