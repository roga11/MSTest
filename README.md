# MSTest

This package implements hypothesis testing procedures that can be used to identify the number of regimes in a Markov-Switching model. It includes the Monte Carlo moment-based test of Dufour & Luger (2017), the parametric bootstrap test described in Qu & Zhuo (2021) and Kasahara & Shimotsu (2018), the Monte Carlo Likelihood ratio tests of Rodriguez Rondon & Dufour (2022), the optimal test for regime switching of Carrasco, Hu, & Ploberger (2014), and the likelihood ratio test described in Hansen (1992).

In addition to testing procedures, the package also includes functions that can be used to simulate: autoregressive, vector autoregressive, Markov switching autoregressive, and Markov switching vector autoregressive processes among others. Model estimation procedures are also available.

For a more detailed description of this package see Rodriguez Rondon & Dufour (2022).

## Installation

Currently, this package is only available through github. To install it you can use the following (requires 'devtools' package): 

```r
devtools::install_github("roga11/MSTest")
```

## Examples

This first example uses the US GNP growth data from 1951Q2-1984Q4 considered in Hamilton (1989). The data is made available as 'hamilton84GNP' through this package. In Hamilton (1989), the model is estimated with four autoregressive lags and only the mean is allowed to change between two (i.e., expansionary and recessionary) regimes and it is estimated by MLE and so we begin by estimating that model.  Estimation results can be compared with those found in Hamilton (1994) p. 698. Note however, that standard errors here were obtained using a different approximation method and hence these may differ slightly. 

```r
y_gnp_gw_84 <- hamilton84GNP$GNP_logdiff


# Set options for model estimation
control <- list(msmu   = TRUE, 
                msvar  = FALSE, 
                method = "MLE",
                use_diff_init = 10)

# Estimate model with p=4 and switch in mean only as in Hamilton (1989)
hamilton89_mdl <- MSARmdl(as.matrix(y_gnp_gw_84), p = 4, k = 2, control)

hamilton89_mdl

# plot smoothed probability of recessionary state
plot(hamilton89_mdl$St[,2], type = 'l')
```
 
This next example uses the 'simuMSAR' function to simulate a Markov switching autoregressive process and then uses 'MSARmdl' to estimate the model. Estimated coefficients may be compared with the true parameters used to generate the data. A plot also shows the fit of the smoothed probabilities. 

```r
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

# Set options for model estimation
control <- list(msmu   = TRUE, 
                msvar  = TRUE, 
                method = "EM",
                use_diff_init = 10)

# Estimate model
y_ms_mdl <- MSARmdl(y_ms_simu$y, p = 2, k = 2, control)

y_ms_mdl

plot(y_ms_mdl$St[,1], type = 'l')
lines(y_ms_simu$St, col = 'red', lty = 2)
```

This third example, the 'simuMSVAR' function to simulate a bivariate Markov switching vector autoregressive process and then uses 'MSVARmdl' to estimate the model. Estimated coefficients may be compared with the true parameters used to generate the data. A plot also shows the fit of the smoothed probabilities. 

```r
set.seed(1234)
# Define DGP of MS VAR process
mdl_msvar2 <- list(n     = 1000, 
                   p     = 1,
                   mu    = rbind(c(5,-2),
                                 c(10,2)),
                   sigma = list(rbind(c(5.0, 1.5),
                                      c(1.5, 1.0)),
                                rbind(c(7.0, 3.0),
                                      c(3.0, 2.0))),
                   phi   = rbind(c(0.50, 0.30),
                                 c(0.20, 0.70)),
                   k     = 2,
                   P     = rbind(c(0.90, 0.10),
                                 c(0.10, 0.90)))

# Simulate process using simuMSVAR() function
y_msvar_simu <- simuMSVAR(mdl_msvar2)

# Set options for model estimation
control <- list(msmu   = TRUE, 
                msvar  = TRUE,
                method = "EM",
                use_diff_init = 10)
                
# Estimate model
y_msvar_mdl <- MSVARmdl(y_msvar_simu$y, p = 1, k = 2, control)

y_msvar_mdl

plot(y_msvar_mdl$St[,2], type = 'l')
lines(y_msvar_simu$St, col = 'red', lty = 2)


```

## References

Carrasco, Marine, Liang Hu, and Werner Ploberger. 2014. “Optimal test for Markov switching parameters.” Econometrica 82 (2): 765–784.

Dufour, Jean-Marie, and Richard Luger. 2017. “Identification-robust moment-based tests for Markov switching in autoregressive models.” Econometric Reviews 36 (6-9): 713–727.

Kasahara, Hiroyuki, and Katsumi Shimotsu. 2018. “Testing the number of regimes in Markov regime switching models.” arXiv preprint arXiv:1801.06862.

Hamilton, James D. 1989. “A new approach to the economic analysis of nonstationary time series and the business cycle.” Econometrica: Journal of the Econometric Society, 357–384.

Hamilton, James Douglas. 1994. Time series analysis. Princeton university press.

Hansen, Bruce E. 1992. “The likelihood ratio test under nonstandard conditions: testing the Markov switching model of GNP.” Journal of applied Econometrics 7 (S1): S61–S82.

Rodriguez Rondon, Gabriel and Jean-Marie Dufour. 2022. “Monte Carlo Likelihood Ratio Tests for Markov Switching Models.” Unpublished manuscript.

Rodriguez Rondon, Gabriel and Jean-Marie Dufour. 2022. “MSTest: An R-package for Testing Markov-Switching Models.” Unpublished manuscript.

Qu, Zhongjun, and Fan Zhuo. 2021. “Likelihood Ratio-Based Tests for Markov Regime Switching.” The Review of Economic Studies 88 (2): 937–968
