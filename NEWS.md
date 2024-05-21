# MSTest 0.1.2.9000
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

# MSTest 0.1.1

* Made changes to MMC LRT related functions for obtaining null distribution of statistical test.
* Added more examples. Examples that take long to complete are commented out but serve to get familiar with usage. 
* Fixed bug related to using init_theta (setting initial values of parameters) when estimating Markov models (i.e. MSARmdl(), HMmdl(), and MSVARmdl()). 
* Fixed bug in MMC_bound() when k0>1
* Updated `DESCRIPTION` file for new version.

# MSTest 0.1.0

* Added a `NEWS.md` file to track changes to the package.
* Added a `README.md` file to describe the package.
* First public version of package.
