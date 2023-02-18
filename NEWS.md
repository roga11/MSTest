# MSTest 0.1.1.9000
* Made change to MSARmdl(), MSVARmdl(), and HMmdl(). Specifically, when msmu or msvar is FALSE, output list copies single regime value for each k regime. This is needed for simulating the null when either msmu or msvar is FALSE in LMCLRTest() and MMCLRTest().
* Added option to use different number of initial values for estimating MSMs with observed data vs. for null distribution (see documentation for LMCLRTest() and MMCLRTest()).
* Updated USGNP data set to include 2022 Q3.
* MMCLRTest() now has the option to add lower and upper bounds for autoregressive coefficients to help reduce errors from polyroot() when optimizing. MMC_bounds() has been updated to reflect this. 

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
