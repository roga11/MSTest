

y84 <- as.matrix(hamilton84GNP$GNP_logdiff)


Y <- y84


hlrt_control  <- list(ix = 1, 
                      p_gridsize = 17,
                      p_stepsize = 0.05)


hlrt <- HLRTest(Y, p = 4, control = hlrt_control)
hlrt



