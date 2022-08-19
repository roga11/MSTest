# ------------------------------------------------------------------------------------------
#' @title CHP parameter stability test
#'
#' @description This function performs the CHP parameter stability test as outline in Carrasco, M., Hu, L. and Ploberger, W. (2014).
#' Original source code can be found \href{https://www.econometricsociety.org/content/supplement-optimal-test-markov-switching-parameters}{here}.
#'
#' @param Y A (\code{T x 1}) matrix of observations.  
#' @param p Integer determining the number of autoregressive lags.
#' @param control List with test procedure options including: 
#' \itemize{
#'   \item{\code{N}: }{Integer determining the number of Bootstrap iterations. Default is set to \code{3000} as in paper.}
#'   \item{\code{rho_b}: }{Number determining bounds for distribution of \code{rh0} (i.e. \code{rho} ~ \code{[-rho_b,rho_b]}).}
#'   \item{\code{msvar}: }{Boolean indicator. If \code{TRUE}, there is a switch in variance. If \code{FALSE} only switch in mean is considered.}
#'   \item{\code{getSE}: }{Boolean indicator. If \code{TRUE}, standard errors for restricted model are estimated. If \code{FALSE} no standard errors are estimated. Default is \code{TRUE}.}
#' }
#' 
#' @return List of class \code{CHPTest} (\code{S3} object) with model attributes including: 
#' \itemize{
#'   \item{\code{mdl_h0}: }{List with restricted model attributes. This will be of class \code{ARmdl} (\code{S3} object). See \code{\link{ARmdl}}.}
#'   \item{\code{supTS}: }{supTS test statistic value.}
#'   \item{\code{expTS}: }{expTS test statistic value.}
#'   \item{\code{supTS_N}: }{A (\code{N x 1}) vector with simulated supTS test statistics under null hypothesis.}
#'   \item{\code{expTS_N}: }{A (\code{N x 1}) vector with simulated expTS test statistics under null hypothesis.}
#'   \item{\code{pval_supTS}: }{P-value for supTS version of parameter stability test.}
#'   \item{\code{pval_expTS}: }{P-value for expTS version of parameter stability test.}
#'   \item{\code{supTS_cv}: }{Vector with 90\%, 95\%, and 99\% bootstrap critical values for supTS version of parameter stability test.}
#'   \item{\code{expTS_cv}: }{Vector with 90\%, 95\%, and 99\% bootstrap critical values for expTS version of parameter stability test.}
#'   \item{\code{control}: }{List with test procedure options used.}
#' }
#' 
#' @references Carrasco, Marine, Liang Hu, and Werner Ploberger. 2014. “Optimal test for Markov switching parameters.” \emph{Econometrica} 82 (2): 765–784.
#' 
#' @example /examples/CHPTest_examples.R
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
  # --------------- Get first and second derivatives 
  ltmt    <- chpDmat(mdl_h0, con$msvar)
  # --------------- calculate supTS and expTS test statistic 
  cv3     <- chpStat(mdl_h0, con$rho_b, ltmt, con$msvar)
  supts   <- cv3[1]
  expts   <- cv3[2]
  # --------------- Bootstrap Critival Values
  SN        <- bootCV(mdl_h0, con$rho_b, con$N, con$msvar)
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
# ------------------------------------------------------------------------------------------
#' @title Derivative matrix
#'
#' @description This function organizes the first and second derivatives of the log-likelihoood. 
#'
#' @param Mdl List containing output from \code{\link{ARmdl}}.
#' @param msvar Boolean indicator. If \code{TRUE}, there is a switch in variance. If \code{FALSE} only switch in mean is considered.
#' 
#' @return List containing relevant first and second derivatives of log-likelihood function.
#' 
#' @references Carrasco, Marine, Liang Hu, and Werner Ploberger. 2014. “Optimal test for Markov switching parameters.” \emph{Econometrica} 82 (2): 765–784.
#' 
#' @export
chpDmat <-function(mdl, msvar){
  # ----------------- Load Mdl variables
  u       <- as.vector(mdl$resid)
  xtj     <- mdl$x
  v0      <- as.numeric(mdl$stdev)
  nar     <- mdl$p
  phi     <- mdl$phi
  b0      <- mdl$coef
  mu0     <- mdl$mu
  output  <- list()
  # --------------- Get 1st derivative of the log likelihood, use lt as prefix
  # ----- d_lt w.r.t. mu
  ltmu  <- u*as.numeric((1-sum(phi))/(v0^2))
  # ----- d_lt w.r.t. phi parameters
  ltphi <- matrix(nrow=length(u),ncol=0)
  for (xi in 1:nar){
    ltphi <- cbind(ltphi,(u *(xtj[,xi]-mu0))/(v0^2))
  }
  # ----- d_lt w.r.t. sig2
  ltsig2  <- (u^2)/(2*v0^4)-1/(2*v0^2)
  # ----- combine 1st derivatives 
  ltmx    <- cbind(ltmu,ltphi,ltsig2)
  # --------------- Get the 2nd derivative of the log likelihood, use mt as prefix 
  # ----- d_lt_mu w.r.t. mu
  mtmu    <- -(1-sum(phi))^2/(v0^2)
  # ---------- Save first derivatives and second derivative w.r.t. mu
  output$ltmx <- ltmx
  output$mtmu <- mtmu
  if (msvar == TRUE){
    # ---------- If Mean and Var switch
    # ----- d_lt_mu w.r.t. phi
    mtmuphi <- matrix(nrow=length(u),ncol=0)
    for (xi in 1:nar){
      mtmuphi <- cbind(mtmuphi,(-u/(v0^2)+(1-sum(phi))*(mu0-xtj[,xi])/(v0^2)))
    }
    # ----- d_lt_mu w.r.t. sig2
    mtmusig2  <- -u*(1-sum(phi))/(v0^4)
    # ----- d_lt_phi w.r.t. mu
    # same as mtmuphi
    # ----- d_lt_phi w.r.t. phi
    # done in other loop
    # ----- d_lt_phi w.r.t. sig2
    mtphisig2 <- matrix(nrow=length(u),ncol=0)
    for (xi in 1:nar){
      mtphisig2 <- cbind(mtphisig2,(u*(mu0-xtj[,xi]))/(v0^4))
    }
    # ----- d_lt_sig2 w.r.t. mu
    # same as mtmusig2
    # ----- d_lt_phi w.r.t. phi
    # same as mtphisig2
    # ----- d_lt_phi w.r.t. sig2
    mtsig2 <- (1/(2*v0^4))-(u^2/(v0^6))
    # --------------- Save output
    # ---------- Save other second derivatives if variance can also switch. 
    output$ltmx       <- ltmx
    output$mtmu       <- mtmu
    output$mtmuphi    <- mtmuphi
    output$mtmusig2   <- mtmusig2
    output$mtphisig2  <- mtphisig2
    output$mtsig2     <- mtsig2
  }
  return(output)
}
