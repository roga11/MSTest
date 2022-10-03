# load data used in Hamilton 1989 and extended data used in CHP 2014 
y84 <- as.matrix(hamilton84GNP$GNP_logdiff)
y10 <- as.matrix(chp10GNP$GNP_logdiff)

# Set test procedure options
control = list(N = 3000, 
               rho_b = 0.7, 
               msvar = FALSE)

# perform test with switch in mean only on Hamilton 1989 data
mdl_84_msmu <- CHPTest(y84, p = 4, control = control)

mdl_84_msmu

# perform test with switch in mean only on extended data used in CHP 2014
mdl_10_msmu <- CHPTest(y10, p = 4, control = control)

mdl_10_msmu