## Submission
This is a resubmission. In this version I have:

* Addressed issues highlighted by CRAN.
* Made change to MSARmdl(), MSVARmdl(), and HMmdl(). Specifically, when msmu or msvar is FALSE, output list copies single regime value for each k regime. This is needed for simulating the null when either msmu or msvar is FALSE in LMCLRTest() and MMCLRTest().
* Added option to use different number of initial values for estimating MSMs with observed data vs. for null distribution (see documentation for LMCLRTest() and MMCLRTest()).
* Updated USGNP data set to include 2022 Q3.
* MMCLRTest() now has the option to add lower and upper bounds for autoregressive coefficients and transition probabilities to help reduce errors from polyroot() when optimizing. MMC_bounds() has been updated to reflect this. 
* DLMMCTest() now has the option to add lower and upper bounds for autoregressive coefficients to help reduce errors from polyroot() when optimizing. DLMMC_bounds() has been created for this. 
* Fixed error in bootCV() used by CHPTest(). When NaN occurs, new draw is used. 
* Added option to allow user to specify 'mle_theta_low' and 'mle_theta_upp' which determine the lower and upper bounds for optimization in HMmdl(), MSARmdl(), and MSVARmdl() when "method='MLE'" is specified.
* print.CHPTest() now used two lines to print description.
* In HLRTest() user can now define entire grid for transition probabilities. 
* added new optional optimization routine for HLRparamSearch() (nloptr::slsqp() can be used now).
* Fixed bug in HLRTest() where grid for sigma options are properly used now. 
* Now using nearest_spd() from pracma instead of nearPD() from lmf package.
* Added classes for simulation functions
* Added new methods, namely: coef, fitted, predict, summary, residuals, nobs, plot
* Made changes to print method.  
* Changed methods for models, namely: logLiklihood now uses logLik, and aic and bic now use AIC and BIC methods.
* Updated USGNP dataset to include data up to end of 2023.  
* Updated `README.md` file  with usage of new methods.
* Updated `DESCRIPTION` file for changed dependencies and new version.

## Test environments
* 

## R CMD check results

There were no ERRORs, WARNINGs or NOTEs.

## Downstream dependencies

There are currently no downstream dependencies for this package.

