## Install Package

# CRAN version
#install.packages("MSTest")

# DEV version
#library(devtools)
#devtools::install_github("roga11/MSTest")

start <- proc.time()

## Package and options
library("MSTest")
library("foreach")
library("doParallel")
options(prompt = "R> ", continue = "+  ", width = 70,
        useFancyQuotes = FALSE)

seed <- 12345

# =========================================================== #
## ----- Data ----- 
# =========================================================== #
data("hamilton84GNP", package = "MSTest")
data("chp10GNP", package = "MSTest")
data("USGNP", package = "MSTest")
data("USRGDP", package = "MSTest")
# =========================================================== #
## ----- Simulation ----- 
# =========================================================== #

set.seed(seed)

### ----- Simulate Multivariate Normal process ----- 
# Define DGP of multivariate normal process
mdl_norm <- list(n     = 500, 
                 q     = 2,
                 mu    = c(5, -2),
                 sigma = rbind(c(5.0, 1.5),
                               c(1.5, 1.0)))
# Simulate process 
simu_norm <- simuNorm(mdl_norm)

### ----- Simulate Autoregressive process ----- 
# Define DGP of AR(2) process
mdl_ar <- list(n     = 500, 
               mu    = 5,
               sigma = 1,
               phi   = c(0.75))
# Simulate process 
simu_ar <- simuAR(mdl_ar)

### ----- Simulate Vector Autoregressive process ----- 
# Define DGP of VAR(2) process
mdl_var <- list(n     = 500, 
                p     = 1,
                q     = 2,
                mu    = c(5, -2),
                sigma = rbind(c(5.0, 1.5),
                              c(1.5, 1.0)),
                phi   = rbind(c(0.50, 0.30),
                              c(0.20, 0.70)))
# Simulate process 
simu_var <- simuVAR(mdl_var)


### ----- Simulate Hidden Markov process ----- 
# Define DGP of HMM
mdl_hmm <- list(n     = 500, 
                q     = 2,
                mu    = rbind(c(5, -2),
                              c(10, 2)),
                sigma = list(rbind(c(5.0, 1.5),
                                   c(1.5, 1.0)),
                             rbind(c(7.0, 3.0),
                                   c(3.0, 2.0))),
                k     = 2,
                P     = rbind(c(0.95, 0.10),
                              c(0.05, 0.90)))
# Simulate process 
simu_hmm <- simuHMM(mdl_hmm)

### ----- Simulate Markov switching Autoregressive process ----- 
# Define DGP of MS AR process
mdl_ms <- list(n     = 500, 
               mu    = c(5,10),
               sigma = c(1,1),
               phi   = c(0.75),
               k     = 2,
               P     = rbind(c(0.95, 0.10),
                             c(0.05, 0.90)))
# Simulate process 
simu_msar <- simuMSAR(mdl_ms)

### ----- Simulate Markov switching Vector Autoregressive process ----- 
# Define DGP of MS VAR process
mdl_msvar <- list(n     = 500, 
                  p     = 1,
                  q     = 2,
                  mu    = rbind(c(5, -2),
                                c(10, 2)),
                  sigma = list(rbind(c(5.0, 1.5),
                                     c(1.5, 1.0)),
                               rbind(c(7.0, 3.0),
                                     c(3.0, 2.0))),
                  phi   = rbind(c(0.50, 0.30),
                                c(0.20, 0.70)),
                  k     = 2,
                  P     = rbind(c(0.95, 0.10),
                                c(0.05, 0.90)))
# Simulate process
simu_msvar <- simuMSVAR(mdl_msvar)


# pdf(file = "simulations.pdf")
par(mfrow=c(3,2))
matplot(simu_norm$y, type = "l", ylab = "", main = "Multivariate normal process",cex.main=1)
matplot(simu_hmm$y, type = "l", ylab = "", main = "Hidden Markov process",cex.main=1)
plot(simu_ar$y, type = "l", ylab = "", main = "Autoregressive process",cex.main=1)
plot(simu_msar$y, type = "l", ylab = "", main = "Markov switching autoregressive process",cex.main=1)
matplot(simu_var$y, type = "l", ylab = "", main = "Vector autoregressive process",cex.main=1)
matplot(simu_msvar$y, type = "l", ylab = "", main = "Markov switching vector autoregressive process",cex.main=1)
# dev.off()


# =========================================================== #
## ----- Model Estimation ----- 
# =========================================================== #
set.seed(seed)
### ----- Estimate Hidden Markov model ----- 
# Set options for model estimation
control <- list(msmu   = TRUE, 
                msvar  = TRUE,
                method = "EM",
                use_diff_init = 30)
# Estimate model
mdl_est_hmm <- HMmdl(simu_hmm[["y"]], k = 2, control = control)
summary(mdl_est_hmm)


set.seed(seed)
### ----- Estimate Markov switching autoregressive model ----- 
# Set options for model estimation
control <- list(msmu   = TRUE, 
                msvar  = FALSE, 
                method = "EM",
                use_diff_init = 30)
# Estimate model
mdl_est_msar <- MSARmdl(simu_msar[["y"]], p = 1, k = 2, control = control)
summary(mdl_est_msar)


set.seed(seed)
### ----- Estimate Markov switching vector autoregressive model ----- 
# Set options for model estimation
control <- list(msmu   = TRUE, 
                msvar  = TRUE,
                method = "EM",
                use_diff_init = 30)
# Estimate model
mdl_est_msvar <- MSVARmdl(simu_msvar[["y"]], p = 1, k = 2, control = control)
summary(mdl_est_msvar)


# plot simulated process, true regime states and model estimated smoothed probabilities
# pdf(file = "MSestim_smoothedprobs.pdf")
par(mfrow=c(3,1))
matplot(simu_hmm$y, type = "l", ylab = "", main = "Hidden Markov process",cex.main=1, col = c("black", "blue"))
lines(simu_hmm$y[,2], type = "l", ylab = "", col = "blue")
par(new = TRUE)
plot(mdl_est_hmm$St[,2], type = "l", ylab = "", col = "green3", lty="dashed", axes = FALSE)
par(new = TRUE)
plot(simu_hmm$St, type = "l", ylab = "", col = "red", lty="dashed", axes = FALSE)
axis(side = 4, at = pretty(range(0,1)))

plot(simu_msar$y[,1], type = "l", ylab = "", main = "Markov switching autoregressive process",cex.main=1)
par(new = TRUE)
plot(mdl_est_msar$St[,1], type = "l", ylab = "", col = "green3", lty="dashed", axes = FALSE)
par(new = TRUE)
plot(simu_msar$St, type = "l", ylab = "", col = "red", lty="dashed", axes = FALSE)
axis(side = 4, at = pretty(range(0,1)))

matplot(simu_msvar$y, type = "l", ylab = "", main = "Markov switching vector autoregressive process",cex.main=1, col = c("black", "blue"))
par(new = TRUE)
plot(mdl_est_msvar$St[,1], type = "l", ylab = "", col = "green3", lty="dashed", axes = FALSE)
par(new = TRUE)
plot(simu_msvar$St, type = "l", ylab = "", col = "red", lty="dashed", axes = FALSE)
axis(side = 4, at = pretty(range(0,1)))
# dev.off()


# =========================================================== #
## ----- Hypothesis Testing ----- 
# =========================================================== #
### ----- Test Markov switching autoregressive process using Rodriguez-Rondon & Dufour (2025) LMC-LRT ----- 
set.seed(seed)
# Set options for testing procedure
lmc_control = list(N = 19,
                   mdl_h0_control = list(const  = TRUE, 
                                         getSE  = FALSE),
                   mdl_h1_control = list(msmu   = TRUE, 
                                         msvar  = TRUE,
                                         getSE  = FALSE,
                                         method = "EM",
                                         use_diff_init = 1))
# Perform Rodriguez-Rondon & Dufour (2025) LMC-LRT
lmclrt <- LMCLRTest(simu_msvar[["y"]], p = 1, k0 = 1 , k1 = 2, control = lmc_control)
summary(lmclrt)



### ----- Test autoregressive process using Rodriguez-Rondon & Dufour (2025) MMC-LRT ----- 
set.seed(seed)
# Set options for testing procedure
mmc_control = list(N = 19,
                   eps = 0.3,
                   threshold_stop = 0.05 + 1e-6, 
                   type = "GenSA",
                   CI_union = FALSE,
                   mdl_h0_control = list(const  = TRUE, 
                                         getSE  = FALSE),
                   mdl_h1_control = list(msmu   = TRUE, 
                                         msvar  = TRUE,
                                         getSE  = FALSE,
                                         method = "EM"),
                   maxit  = 100)
# Perform Rodriguez-Rondon & Dufour (2025) MMC-LRT
mmclrt <- MMCLRTest(simu_norm[["y"]], p = 0, k0 = 1 , k1 = 2, control = mmc_control)
summary(mmclrt)


### ----- Test Markov switching autoregressive process using Dufour & Luger (2017) LMC test ----- 
set.seed(seed)
# Set options for testing procedure
lmc_control = list(N = 99,
                   simdist_N = 10000,
                   getSE = TRUE)
# Perform Dufour & Luger (2017) LMC test 
lmcmoment <- DLMCTest(simu_msar[["y"]], p = 1, control = lmc_control)
colnames(lmcmoment$S0) <- c("M(eps)","V(eps)","S(eps)","K(eps)")
colnames(lmcmoment$F0_min) <- "F(eps)"
colnames(lmcmoment$F0_prod) <- "F(eps)"
summary(lmcmoment)


### ----- Test Markov switching autoregressive process using Dufour & Luger (2017) MMC test ----- 
set.seed(seed)
# Set options for testing procedure
mmc_control <- list(N = 99,
                    getSE = TRUE,
                    eps = 1e-9, 
                    CI_union = TRUE,
                    optim_type = "GenSA",
                    threshold_stop = 0.05 + 1e-6, 
                    maxit = 100)
# Perform Dufour & Luger (2017) MMC test 
mmcmoment <- DLMMCTest(simu_msar[["y"]], p = 1, control = mmc_control)
colnames(mmcmoment$S0_min) <- c("M(eps)","V(eps)","S(eps)","K(eps)")
colnames(mmcmoment$S0_prod) <- c("M(eps)","V(eps)","S(eps)","K(eps)")
colnames(mmcmoment$F0_min) <- "F(eps)"
colnames(mmcmoment$F0_prod) <- "F(eps)"
summary(mmcmoment)


### ----- Test autoregressive process using Carrasco, Hu, & Ploberger (2014) test ----- 
set.seed(seed)
# Set options for testing procedure
chp_control = list(N = 1000, 
                   rho_b = 0.7, 
                   msvar = FALSE)
# Perform Carrasco, Hu, & Ploberger (2014) test
pstabilitytest <- CHPTest(simu_ar[["y"]], p = 1, control = chp_control)
summary(pstabilitytest)


### ----- Test Markov switching autoregressive process using Hansen (1992) LRT ----- 
set.seed(seed)
# Set options for testing procedure
hlrt_control  <- list(msvar          = FALSE,
                      gridsize       = 20,
                      mugrid_from    = 0,
                      mugrid_by      = 1,
                      theta_null_low = c(0,-0.99,0.01),
                      theta_null_upp = c(20,0.99,20))
# Perform Hansen (1992) likelihood ratio test
hlrt <- HLRTest(simu_msar[["y"]], p = 1, control = hlrt_control)
summary(hlrt)







# --------------- END 

end <- proc.time() - start

print(round(end[3]/60,2))
print(end[3]/60 < 10)