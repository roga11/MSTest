set.seed(1234)
# load data used in Hamilton 1989 and extended data used in CHP 2014 
y84 <- as.matrix(hamilton84GNP$GNP_gr)
y10 <- as.matrix(chp10GNP$GNP_gr)

# Set test procedure options
lmc_control = list(N = 99,
                   simdist_N = 10000,
                   getSE = TRUE)

# perform test on Hamilton 1989 data
lmc_gnp84 <- DLMCTest(y84, p = 4, control = lmc_control)
summary(lmc_gnp84)

