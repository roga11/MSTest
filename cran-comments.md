## Submission
This is a resubmission. In this version I have:
* Fixed bug in MMC when null==1 - msmu and msvar undefined in 'mdledit'. result of previous edit when null>1.  
* Updated MLE estimation following deprecation of hin>=0 (inequality constraint direction) in slsqr
* Changed OLS unbiased estimates of models with k=1 to be consistent with MLE estimates. Package is for testing more than estimation so comparison with MLE-based tetsing is prioritized. 
* Changed use of arma::is_finite(X) to std::isfinite(X) because former is now deprecated.

## Test environments
* 

## R CMD check results

There were no ERRORs, WARNINGs or NOTEs.

## Downstream dependencies

There are currently no downstream dependencies for this package.

