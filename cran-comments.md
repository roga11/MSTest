## Submission
This is a resubmission. In this version I have:
* Fixed bug in MMC when null>1 - matrix out of bound error in simulation under null due to 'mdledit' not taking into account if msmu (msvar) are false or not. 
* Fixed bug in estimation of MSVAR and MSVARX when msmu=FALSE due to intercept - correction for repeating output$mu by k had to come before defining inter

## Test environments
* 

## R CMD check results

There were no ERRORs, WARNINGs or NOTEs.

## Downstream dependencies

There are currently no downstream dependencies for this package.

