# --------------------------- Use simulated process ----------------------------
set.seed(1234)
# Define DGP of MS AR process
mdl_ms2 <- list(n     = 200, 
                mu    = c(5,1),
                sigma = c(1,1),
                phi   = c(0.5),
                k     = 2,
                P     = rbind(c(0.90, 0.10),
                              c(0.10, 0.90)))

# Simulate process using simuMSAR() function
y_ms_simu <- simuMSAR(mdl_ms2)

hlrt_control  <- list(ix          = 1, 
                      gridsize    = 5,
                      p_gridsize  = 9,
                      p_stepsize  = 0.1,
                      mugrid_from = 0,
                      mugrid_by   = 1)

\donttest{
  hlrt <- HLRTest(y_ms_simu$y, p = 1, control = hlrt_control)
  summary(hlrt)
}

