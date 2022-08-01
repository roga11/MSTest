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

plot(y_ms_simu$y, type = 'l')
plot(y_ms_simu$St[,1], type = 'l', col = 'red')