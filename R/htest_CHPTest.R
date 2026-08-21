#' @title Carrasco, Hu, and Ploberger (2014) parameter stability test
#'
#' @description This function performs the CHP (2014) parameter stability test as outline in Carrasco, M., Hu, L. and Ploberger, W. (2014).
#' Original source code can be found \href{https://www.econometricsociety.org/content/supplement-optimal-test-markov-switching-parametershere}{here}.
#'
#' @param Y A (\code{T x 1}) matrix of observations.  
#' @param p Integer determining the number of autoregressive lags.
#' @param control List with test procedure options including: 
#' \itemize{
#'   \item N: Integer determining the number of Bootstrap iterations. Default is set to \code{3000} as in paper.
#'   \item rho_b: Number determining bounds for distribution of \code{rh0} (i.e. \code{rho} ~ \code{[-rho_b,rho_b]}).
#'   \item msvar: Boolean indicator. If \code{TRUE}, there is a switch in variance. If \code{FALSE} only switch in mean is considered.
#'   \item getSE: Boolean indicator. If \code{TRUE}, standard errors for restricted model are estimated. If \code{FALSE} no standard errors are estimated. Default is \code{TRUE}.
#' }
#' 
#' @return List of class \code{CHPTest} (\code{S3} object) with model attributes including: 
#' \itemize{
#'   \item mdl_h0: List with restricted model attributes. This will be of class \code{ARmdl} (\code{S3} object). See \code{\link{ARmdl}}.
#'   \item supTS: supTS test statistic value.
#'   \item expTS: expTS test statistic value.
#'   \item supTS_N: A (\code{N x 1}) vector with simulated supTS test statistics under null hypothesis.
#'   \item expTS_N: A (\code{N x 1}) vector with simulated expTS test statistics under null hypothesis.
#'   \item pval_supTS: P-value for supTS version of parameter stability test.
#'   \item pval_expTS: P-value for expTS version of parameter stability test.
#'   \item supTS_cv: Vector with 90\%, 95\%, and 99\% bootstrap critical values for supTS version of parameter stability test.
#'   \item expTS_cv: Vector with 90\%, 95\%, and 99\% bootstrap critical values for expTS version of parameter stability test.
#'   \item control: List with test procedure options used.
#' }
#' 
#' @references Carrasco, Marine, Liang Hu, and Werner Ploberger. 2014. “Optimal test for Markov switching parameters.” \emph{Econometrica} 82 (2): 765–784.
#' @example /inst/examples/CHPTest_examples.R
#' @export
CHPTest <- function(Y, p, control = list()){
  # ----- Set control values
  con <- list(N     = 3000, 
              rho_b = 0.7,
              msvar = FALSE,
              getSE = TRUE)
  # Perform some checks for controls
  nmsC <- names(con)
  con[(namc <- names(control))] <- control
  if(length(noNms <- namc[!namc %in% nmsC])){
    warning("unknown names in control: ", paste(noNms,collapse=", ")) 
  }
  # --------------- Begin by estimating model under H0
  null_control <- list(const = TRUE, getSE = con$getSE)
  mdl_h0  <- ARmdl(Y, p, null_control)
  # --------------- calculate supTS and expTS test statistic
  cv3     <- chpStat(mdl_h0, con$rho_b, con$msvar)
  supts   <- cv3[1]
  expts   <- cv3[2]
  # --------------- Bootstrap Critival Values
  SN        <- CHPbootCV(mdl_h0, con$rho_b, con$N, con$msvar)
  supTS_N   <- as.matrix(sort(SN[,1]))
  expTS_N   <- as.matrix(sort(SN[,2]))
  supTS_cv  <- supTS_N[round(c(0.90,0.95,0.99)*nrow(supTS_N)),]
  expTS_cv  <- expTS_N[round(c(0.90,0.95,0.99)*nrow(expTS_N)),]
  names(supTS_cv) <- paste0(c("0.90","0.95","0.99"), "%")
  names(expTS_cv) <- paste0(c("0.90","0.95","0.99"), "%")
  # --------------- Bootstrap p-value 
  sup_pval  <- sum(supts<supTS_N)/con$N
  exp_pval  <- sum(expts<expTS_N)/con$N
  # --------------- Save Results
  CHPTest_output <- list(mdl_h0 = mdl_h0, supTS = supts, expTS = expts, supTS_N = supTS_N, expTS_N = expTS_N, 
                         pval_supTS = sup_pval, pval_expTS = exp_pval, supTS_cv = supTS_cv, expTS_cv = expTS_cv,
                         control = con)
  class(CHPTest_output) <- "CHPTest"
  return(CHPTest_output)
}
