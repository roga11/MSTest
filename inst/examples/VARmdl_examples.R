# ----- Bivariate VAR(1) process ----- #
set.seed(1234)
# Define DGP of VAR process
mdl_var <- list(n     = 1000, 
                p     = 1,
                q     = 2,
                mu    = c(5,-2),
                sigma = rbind(c(5.0, 1.5),
                              c(1.5, 1.0)),
                phi   = rbind(c(0.50, 0.30),
                              c(0.20, 0.70)))

# Simulate process using simuVAR() function
y_simu <- simuVAR(mdl_var)

# Set options for model estimation
control <- list(const  = TRUE, 
                getSE  = TRUE)

# Estimate model
y_var_mdl <- VARmdl(y_simu$y, p = 2, control = control)
summary(y_var_mdl)
