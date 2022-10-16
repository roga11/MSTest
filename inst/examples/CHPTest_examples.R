# load data used in Hamilton 1989 
y84 <- as.matrix(hamilton84GNP$GNP_logdiff)


# Set test procedure options
control = list(N = 1000, 
               rho_b = 0.7, 
               msvar = FALSE)

# perform test with switch in mean only on Hamilton 1989 data
mdl_84_msmu <- CHPTest(y84, p = 4, control = control)
mdl_84_msmu

