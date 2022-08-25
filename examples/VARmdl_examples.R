# -------------------------- Bivariate VAR(1) process --------------------------
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
y_var_mdl <- VARmdl(y_simu$y, p = 1, control)

y_var_mdl

# ----------------------- VAR(2) process with 3 variables ----------------------
set.seed(1234)
# Define DGP of VAR process
mdl_3var2 <- list(n     = 1000, 
                  p     = 2,
                  q     = 3,
                  mu    = c(5, -2, 1),
                  sigma = rbind(c(5.0, 1.5, 2.5),
                                c(1.5, 1.0, 1.5),
                                c(2.5, 1.5, 4.2)),
                  phi   = rbind(c(0.70, 0.30, 0.35,  -0.50, -0.20,   0.25),
                                c(0.20, 0.40, 0.35,  -0.30,  0.30,   0.25),
                                c(0.20, 0.30, 0.50,  -0.30, -0.20,  -0.40)))

# Simulate process using simuVAR() function
y3var2_simu <- simuVAR(mdl_3var2)

# Set options for model estimation
control <- list(const  = TRUE, 
                getSE  = TRUE)

# Estimate model
y_3var2_mdl <- VARmdl(y3var2_simu$y, p = 2, control)

y_3var2_mdl
