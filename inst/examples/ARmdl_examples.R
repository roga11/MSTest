set.seed(1234)
# Define DGP of AR process
mdl_ar <- list(n     = 500, 
               mu    = 5,
               sigma = 2,
               phi   = c(0.5,0.2))

# Simulate process using simuAR() function
y_simu <- simuAR(mdl_ar)

# Set options for model estimation
control <- list(const  = TRUE, 
                getSE  = TRUE)

# Estimate model
y_ar_mdl <- ARmdl(y_simu$y, p = 2, control)
y_ar_mdl

