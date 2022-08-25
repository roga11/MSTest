set.seed(1234)
# Define DGP of VAR process
mdl_3var2 <- list(n     = 1000, 
                  p     = 2,
                  q     = 3,
                  mu    = c(5, -2, 1),
                  sigma = rbind(c(5.0, 1.5, 2.5),
                                c(1.5, 1.0, 1.5),
                                c(2.5, 1.5, 4.2)),
                  phi   = rbind(c(0.70, 0.30, 0.35,  -0.50, -0.20,   0.25),
                                c(0.20, 0.40, 0.35,  -0.30,  0.30,   0.25),
                                c(0.20, 0.30, 0.50,  -0.30, -0.20,  -0.40)))

# Simulate process using simuVAR() function
y3var2_simu <- simuVAR(mdl_3var2)


plot(y3var2_simu$y[,1], type = 'l')
plot(y3var2_simu$y[,2], type = 'l')
plot(y3var2_simu$y[,3], type = 'l')