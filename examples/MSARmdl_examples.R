

# ----------------------- Use GNP growth 1951Q2 - 1984Q4 ----------------------- 
y_gnp_gw_84 <- hamilton84GNP$GNP_logdiff


# Set options for model estimation
control <- list(msmu   = TRUE, 
                msvar  = FALSE, 
                method = "MLE",
                use_diff_init = 10)

# Estimate model with p=4 and switch in mean only as in Hamilton (1989)
hamilton89_mdl <- MSARmdl(as.matrix(y_gnp_gw_84), p = 4, k = 2, control)

hamilton89_mdl

# plot smoothed probability of each state
plot(hamilton89_mdl$St[,1], type = 'l')
plot(hamilton89_mdl$St[,2], type = 'l')

# --------------------------- Use simulated process ----------------------------
set.seed(1234)
# Define DGP of MS AR process
mdl_ms2 <- list(n     = 500, 
                mu    = c(5,10),
                sigma = c(1,2),
                phi   = c(0.5,0.2),
                k     = 2,
                P     = rbind(c(0.90, 0.10),
                              c(0.10, 0.90)))

# Simulate process using simuMSAR() function
y_ms_simu <- simuMSAR(mdl_ms2)

# Set options for model estimation
control <- list(msmu   = TRUE, 
                msvar  = TRUE, 
                method = "EM",
                use_diff_init = 10)

# Estimate model
y_ms_mdl <- MSARmdl(y_ms_simu$y, p = 2, k = 2, control)

y_ms_mdl



