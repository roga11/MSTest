set.seed(1234)

# ----- Univariate ----- #
# Define DGP 
mdl_hmm <- list(n     = 1000, 
                q     = 1,
                mu    = as.matrix(c(5,
                                   -2)),
                sigma = list(as.matrix(5.0),
                             as.matrix(7.0)),
                k     = 2,
                P     = rbind(c(0.90, 0.10),
                              c(0.10, 0.90)))

# Simulate process using simuHMM() function
y_hmm_simu <- simuHMM(mdl_hmm)

# ----- Multivariate ----- #
# Define DGP 
mdl_hmm <- list(n     = 1000, 
                q     = 2,
                mu    = rbind(c(5, -2),
                              c(10, 2)),
                sigma = list(rbind(c(5.0, 1.5),
                                   c(1.5, 1.0)),
                             rbind(c(7.0, 3.0),
                                   c(3.0, 2.0))),
                k     = 2,
                P     = rbind(c(0.90, 0.10),
                              c(0.10, 0.90)))

# Simulate process using simuHMM() function
y_hmm_simu <- simuHMM(mdl_hmm)

