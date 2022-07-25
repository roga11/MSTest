

# Use GNP growth 1951Q2 - 1984Q4
y_gnp_gw_84 <- hamilton84GNP$GNP_logdiff


# Set options for model estimation
control <- list(msmu  = TRUE, 
                msvar = FALSE, 
                use_diff_init = 10)

# Estimate model with p=4 and switch in mean only as in Hamilton (1989)
hamilton89_mdl <- MSARmdl_EM(y_gnp_gw_84, ar = 4, k = 2, control)

# plot smoothed probability of each state
plot(hamilton89_mdl$St, type = 'l')

