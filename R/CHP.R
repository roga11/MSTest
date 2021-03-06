# ------------------------------------------------------------------------------------------
#' @title CHP Test
#'
#' @description This function performs the CHP test as outline in Carrasco, M., Hu, L. and Ploberger, W. (2014).
#' This function can be used to recreate results from Table III of the paper. Econometrica as the original 
#' publisher and source code can be found here: https://www.econometricsociety.org/content/supplement-optimal-test-markov-switching-parameters.
#'
#' @param Y the series to be tested
#' @param p Number of autoregressive lags AR(p)
#' @param N Number of Bootstrap iterations. Default is set to 3000 as in paper. 
#' @param rho_b bound for apriori distribution of rh0 (i.e. rho ~ \([-rho_b,rho_b]\)
#' @param var_switch Indicator for switch in Variance (if = 0, only Mean is subject to switch)
#' 
#' @return 
#' 
#' @references Carrasco, Marine, Liang Hu, and Werner Ploberger. 2014. 
#' “Optimal test for Markov switch- ing parameters.” \emph{Econometrica} \bold{82 (2)}: 765–784.
#' 
#' @export
CHPtest <- function(Y, p=1, N=3000, rho_b=0.7, var_switch = 0){
  # --------------- Begin by estimating model under H0
  Mdl   <- ARmdl(Y,p)
  # --------------- Get first and second derivatives 
  ltmt  <- chpDmat(Mdl, var_switch)
  # --------------- calculate supTS and expTS test statistic 
  cv3   <- chpStat(Mdl, rho_b, ltmt, var_switch)
  supts <- cv3[1]
  expts <- cv3[2]
  # --------------- Bootstrap Critival Values
  SN    <- bootCV(Mdl, rho_b, N, var_switch)
  supb  <- SN[,1]
  expb  <- SN[,2]
  # --------------- Bootstrap p-value 
  sup_pval  <- round(mean(supts<supb),3)
  exp_pval  <- round(mean(expts<expb),3)
  # --------------- Save Results
  cv1             <- as.matrix(supb[round(c(0.90,0.95,0.99)*length(supb))])
  cv2             <- as.matrix(expb[round(c(0.90,0.95,0.99)*length(expb))])
  row.names(cv1)  <- c('0.90%','0.95%','0.99%')
  row.names(cv2)  <- c('0.90%','0.95%','0.99%')
  phi             <- t(as.matrix(Mdl$phi))
  colnames(phi)   <- paste0('ar',1:Mdl$ar)
  output        <- list(supts,t(cv1),sup_pval,phi,expts,t(cv2),exp_pval,phi)
  names(output) <- c('Test-Stat_supTS','Crit-Values_supTS','p-value_supTS','params_supTS',
                     'Test-Stat_expTS','Crit-Values_expTS','p-value_expTS','params_expTS')
  return(output)
} 
# ------------------------------------------------------------------------------------------
#' @title CHP Monte Carlo Test
#'
#' @description This function performs a Monte Carlo test by drawing from the asymptotic distribution of 
#' supTS and expTS test-statistics. It is recommended to use a large number of iterations.
#'
#' @param Y the series to be tested
#' @param p Number of autoregressive lags AR(p)
#' @param N Number of Bootstrap iterations
#' @param rho_b bound for apriori distribution of rh0 (i.e. rho ~ \([-rho_b,rho_b]\)
#' @param var_switch Indicator for switch in Variance (if = 0, only Mean is subject to switch)
#' 
#' @references Carrasco, Marine, Liang Hu, and Werner Ploberger. 2014. 
#' “Optimal test for Markov switch- ing parameters.” \emph{Econometrica} \bold{82 (2)}: 765–784.
#' 
#' @export
CHPMCtest <- function(Y, p=1, N=3000, rho_b=0.7, var_switch = 0){
  # ** Consider adding MMC to this MC approach
  # --------------- Begin by estimating model under H0
  Mdl   <- ARmdl(Y,p)
  # --------------- Get first and second derivatives 
  ltmt  <- chpDmat(Mdl, var_switch)
  # --------------- calculate supTS and expTS test statistics 
  cv3   <- chpStat(Mdl, rho_b, ltmt, var_switch)
  supts <- cv3[1]
  expts <- cv3[2]
  # --------------- Compute Assymptotic p-value from asymptotic critical values 
  asymdist    <- apply(CHPadist(N,Mdl$n,rho_b),MARGIN=2,sort)
  sup_pval    <- p_val(supts, asymdist[,1], type = 'geq')
  exp_pval    <- p_val(expts, asymdist[,2], type = 'geq')
  # --------------- Save Results
  cv1         <- as.matrix(asymdist[round(c(0.90,0.95,0.99)*nrow(asymdist)),1])
  cv2         <- as.matrix(asymdist[round(c(0.90,0.95,0.99)*nrow(asymdist)),2])
  row.names(cv1)  <- c('0.90%','0.95%','0.99%')
  row.names(cv2)  <- c('0.90%','0.95%','0.99%')
  phi             <- t(as.matrix(Mdl$phi))
  colnames(phi)   <- paste0('ar',1:Mdl$ar)
  output          <- list(supts,t(cv1),round(sup_pval,2),phi,
                          expts,t(cv2),round(exp_pval,2),phi)
  names(output)   <- c('Test-Stat_supTS','Crit-Values_supTS','p-value_supTS','params_supTS',
                        'Test-Stat_expTS','Crit-Values_expTS','p-value_expTS','params_expTS')
  return(output)
}
# ------------------------------------------------------------------------------------------
#' @title Derivative matrix
#'
#' @description This function organizes the first and second derivatives of the log Likelihoood. 
#'
#' @param Mdl is a list containing AR model components
#' Specifically, it containg y, x, X (x plus constant), residuals, coefficients, stdev,
#' logLike
#' @param var_switch is an indicator = 1 if there is a switch in both Mean and Variance 
#' and = 0 if there is only a switch in the Mean. Less Second-Order dervatives are 
#' calculated if only the Mean is subject to regime switch. 
#' 
#' @return List containing relevant first and second derivatves of log likelihood function.
#' 
#' @export
chpDmat <-function(Mdl, var_switch = 1){
  # ----------------- Load Mdl variables
  u       <- as.vector(Mdl$residuals)
  xtj     <- Mdl$x
  v0      <- as.numeric(Mdl$stdev)
  nar     <- Mdl$ar
  phi     <- Mdl$coef[2:length(Mdl$coef)]
  b0      <- Mdl$coef
  mu0     <- Mdl$coef[1]/(1-sum(phi))
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
  if (var_switch==1){
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
