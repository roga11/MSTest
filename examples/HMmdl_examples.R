# Define DGP of MS VAR process
mdl_hmm <- list(n     = 1000, 
                mu    = rbind(c(5,-2),
                              c(10,2)),
                sigma = list(rbind(c(5.0, 1.5),
                                   c(1.5, 1.0)),
                             rbind(c(7.0, 3.0),
                                   c(3.0, 2.0))),
                k     = 2,
                P     = rbind(c(0.98, 0.02),
                              c(0.02, 0.98)))

# Simulate process using simuHMM() function
y_hmm_simu <- simuHMM(mdl_hmm)

# Set options for model estimation
control <- list(msmu   = TRUE, 
                msvar  = TRUE,
                method = "EM",
                use_diff_init = 10)

# Estimate model
y_hmm_mdl <- HMmdl(y_hmm_simu$y, k = 2, control)

y_hmm_mdl

plot(y_hmm_mdl$St[,1], type = 'l')
lines(y_hmm_simu$St, col='red')


# init_mdl <- Nmdl(y_hmm_simu$y, control)
# theta_0 <- initVals_HMmdl(init_mdl,2)
# init_mdl$msmu <- control$msmu
# init_mdl$msvar <- control$msvar
# theta<- theta_0
# 
# 
# out_loglike <- ExpectationM_HMmdl(theta, init_mdl, 2)
# out_max <- EMaximization_HMmdl(theta, init_mdl, out_loglike, 2)
# out_loglike$logLike
# max(abs(out_max$theta - theta))
# theta<-out_max$theta
# 
# out<- HMmdl_em(theta_0, init_mdl, 2, list(maxit=10000, thtol=1e-6))