set.seed(1234)
# Define DGP 
mdl_norm <- list(n     = 1000, 
                 q     = 2,
                 mu    = c(5, -2),
                 sigma = rbind(c(5.0, 1.5),
                               c(1.5, 1.0)))

# Simulate process using simuNorm() function
y_norm_simu <- simuNorm(mdl_norm)

plot(y_norm_simu)