MMCenv <- new.env()
# -----------------------------------------------------------------------------
#' @title Monte-Carlo Moment-based test for MS AR model
#'
#' This function performs the Local Monte-Carlo Moment-Based test for
#' MS AR models presented in Dufour & Luger (2017) (i.e when no nuissance 
#' parameters are present). 
#'
#' @param Y Series to be tested 
#' @param p Order of autoregressive components AR(p).
#' @param x exogenous variables if any. Test in Dufour & Luger is model for AR lags
#' @param N number of samples
#' @param N2 number of simulations when approximating distribution used to combine 
#' p-values (eq. 16).
#'
#' @return List with model and test results.
#' 
#' @references Dufour, J. M., & Luger, R. (2017). Identification-robust moment-based 
#' tests for Markov switching in autoregressive models. Econometric Reviews, 36(6-9), 713-727.
#' 
#' @export
DLMCtest <- function(Y, p = NULL, x = NULL, N = 100, N2 = 10000){
  # ---------------------------------------------------------------------------
  if (is.null(p)==FALSE & is.null(x)==TRUE){
    odr         <- c(p,0)
    mdl         <- estARMA(Y,odr)
    # define z (p.720 Dufour & Luger 2017 - Econometric reviews)
    z           <- mdl$y - mdl$x %*% mdl$phi
    # See if can also mix ARMA and X. 
    # --------- Get parameters from approximated distribution ----------
    params      <- approxDist(length(Y),N2)
    # -------------------------- Get P-Values --------------------------
    LMC_ans     <- calc_mcstat((z - mean(z)),N,params)
  }else if (is.null(p)==TRUE & is.null(x)==TRUE){
    z           <- Y
    # --------- Get parameters from approximated distribution ----------
    params      <- approxDist(length(Y),N2)
    # -------------------------- Get P-Values --------------------------
    LMC_ans     <- calc_mcstat(z - mean(z),N,params)
  }else if (is.null(odr)==TRUE & is.null(x)==FALSE){
    X           <- as.matrix(cbind(1,x))
    B           <- solve(t(X)%*%X)%*%t(X)%*%Y
    z           <- Y - X%*%B
    # --------- Get parameters from approximated distribution ----------
    params      <- approxDist(length(Y),N2)
    # -------------------------- Get P-Values --------------------------
    LMC_ans     <- calc_mcstat(z - mean(z),N,params)
  }
  S0            <- LMC_ans[nrow(LMC_ans),]
  SN            <- apply(LMC_ans[1:(nrow(LMC_ans)-1),],MARGIN=2, sort)
  cv            <- SN[round(c(0.90,0.95,0.99)*nrow(SN)),]
  row.names(cv) <- c('0.90%','0.95%','0.99%')
  pval_Fmin     <- p_val(S0[1], SN[,1], type = 'geq')
  pval_Fprod    <- p_val(S0[2], SN[,2], type = 'geq')
  if (is.null(odr)==FALSE & is.null(x)==TRUE){
    outp        <- list(S0[1],cv[,1],pval_Fmin,mdl$phi,S0[2],cv[,2],
                        pval_Fprod,mdl$phi)
    names(outp) <- c('Test-Stat_min','Crit-Values_min','p-value_min',
                       'params_min',
                       'Test-Stat_prod','Crit-Values_prod','p-value_prod',
                       'params_prod')
  }else if (is.null(odr)==TRUE & is.null(x)==TRUE){
    outp        <- list(S0[1],cv[,1],pval_Fmin,S0[2],cv[,2],pval_Fprod)
    names(outp) <- c('Test-Stat_min','Crit-Values_min','p-value_min',
                       'Test-Stat_prod','Crit-Values_prod','p-value_prod')
  }else if (is.null(odr)==TRUE & is.null(x)==FALSE){
    outp        <- list(S0[1],cv[,1],pval_Fmin,S0[2],cv[,2],pval_Fprod)
    names(outp) <- c('Test-Stat_min','Crit-Values_min','p-value_min',
                       'Test-Stat_prod','Crit-Values_prod','p-value_prod')
  }
  return(outp)
}
# -----------------------------------------------------------------------------
#' @title Maximized Monte-Carlo Moment-based test for MS AR model
#'
#' This function performs the MMC version of the Moment-Based test for
#' MS AR models presented in Dufour & Luger (2017). It is useful when
#' nuissance parameters are present in the null disribution. 
#'
#' @param Y Series to be tested 
#' @param p Order of autoregressive components AR(p).
#' @param x exogenous variables if any. Test in Dufour & Luger is model for AR lags
#' @param N number of samples
#' @param N2 number of simulations when approximating distribution used to combine 
#' p-values (eq. 16).
#' @param N3 number of parameter values to try for nuissance params. Used only 
#' when caling gridSearch_paramCI or randSearch_paramCI.
#' @param searchType Type of optimization algorithm when searching nuissance 
#' parameter space. Avaiable options are: GenSA, GA, PSO, randSearch_paramCI and 
#' gridSearch_paramCI. Default is set to randSearch_paramCI to match results in paper.
#'
#' @return List with model and test results.
#' 
#' @references Dufour, J. M., & Luger, R. (2017). Identification-robust moment-based 
#' tests for Markov switching in autoregressive models. Econometric Reviews, 36(6-9), 713-727.
#' 
#' @export
DLMMCtest <- function(Y,
                      p = NULL, 
                      x = NULL, 
                      N = 100, 
                      N2 = 10000, 
                      N3 = 100000,
                      searchType = 'randSearch_paramCI',
                      optimOptions = NULL){
  # ------------------------- Estimate Model -------------------------
  if (is.null(p)==FALSE & is.null(x)==TRUE){
    odr       <- c(p,0) 
    mdl       <- estARMA(Y,odr)
    # define z (p.720 Dufour & Luger 2017 - Econometric reviews)
    z         <- mdl$y - mdl$x %*% mdl$phi
    # Since test only depends on errors. See if we can extend to general 
    # model with no ARMA but X matrix of explanatory variables instead. 
    # See if can also mix ARMA and X. 
  }else if (is.null(p)==TRUE & is.null(x)==TRUE){
    stop('No explanatory variables were given. No parameters to treat as nuisance parameters.')
  }else if (is.null(p)==TRUE & is.null(x)==FALSE){
    if (is.null(optimOptions)==1){
      warning('Consider providing limits using optimOptions. Default: [-1,1]') 
    }
    X         <- as.matrix(cbind(1,x))
    B         <- solve(t(X)%*%X)%*%t(X)%*%Y
    z         <- Y - X%*%B
    se0       <- sqrt(diag(solve(t(X)%*%X)%*%t(X)%*%z%*%t(z)%*%X%*%solve(t(X)%*%X)))
  }
  # --------- Get parameters from approximated distribution ----------
  params <- approxDist(length(Y),N2)
  # --------- Assign Variable Values to MMC Environment --------------
  
  # ***** SOME OF THESE VARS MAY NOT EXIST IF NOT ARMA MODEL
  assign('N',N,envir = MMCenv)
  assign('N2',N2,envir = MMCenv)
  assign('N3',N3,envir = MMCenv)
  assign('params',params,envir = MMCenv)
  # -------------------------- Get P-Values --------------------------
  if (is.null(odr)==FALSE & is.null(x)==TRUE){
    assign('y',mdl$y,envir = MMCenv)
    assign('x',mdl$x,envir = MMCenv)
    assign('phi',mdl$phi,envir = MMCenv)
    assign('se0',mdl$se0,envir = MMCenv)
    assign('npar',mdl$npar,envir = MMCenv)
    assign('nar',mdl$nar,envir = MMCenv)
    # Calc mmc pvalue
    output  <- calc_mmcpval(searchType,optimOptions)
    if (searchType=='GenSA' | searchType=='GA' | searchType=='PSO'){
      outputlst = output 
    }else{
      cv1 <- as.matrix(output[1,3:5])
      cv2 <- as.matrix(output[2,3:5])
      param1 <- as.matrix(output[1,6:ncol(output)])
      param2 <- as.matrix(output[2,6:ncol(output)])
      rownames(cv1) <- c('0.90%','0.95%','0.99%')
      rownames(cv2) <- c('0.90%','0.95%','0.99%')
      rownames(param1) <- paste0('ar',seq(1,mdl$nar))
      rownames(param2) <- paste0('ar',seq(1,mdl$nar))
      outputlst <- list(output[1,2],t(cv1),output[1,1],t(param1),
                        output[2,2],t(cv2),output[2,1],t(param2))
      names(outputlst) <- c('Test-Stat_min','Crit-Values_min','p-value_min',
                            'params_min',
                            'Test-Stat_prod','Crit-Values_prod','p-value_prod',
                            'params_prod')
    }
  }else if (is.null(odr)==TRUE & is.null(x)==FALSE){
    assign('y',Y,envir = MMCenv)
    assign('x',X,envir = MMCenv)
    assign('phi',B,envir = MMCenv)
    assign('se0',se0,envir = MMCenv)
    assign('npar',ncol(X),envir = MMCenv)
    assign('nar',0,envir = MMCenv)
    # Calc mmc pvalue
    output  <- calc_mmcpval(searchType,optimOptions)
    if (searchType=='GenSA' | searchType=='GA' | searchType=='PSO'){
      outputlst = output 
    }else{
      outputlst <- list(output[1,2],output[1,3:5],output[1,1],
                        output[2,2],output[2,3:5],output[2,1])
      names(outputlst) <- c('Test-Stat_min','Crit-Values_min','p-value_min',
                            'Test-Stat_prod','Crit-Values_prod','p-value_prod')
    }
  }
  return(outputlst)
}
# -----------------------------------------------------------------------------
#' @title Approximate Distribuion 
#'
#' This function obtains the parameters needed in eq. 16 which is used for 
#' combining p-values.
#'
#' @param Tsize sample size
#' @param N2 number of draws (simulations)
#' @return params the paramters gamma in eq. 16 of Dufour & Luger (2017). 
#' 
#' @references Dufour, J. M., & Luger, R. (2017). Identification-robust moment-based 
#' tests for Markov switching in autoregressive models. Econometric Reviews, 36(6-9), 713-727.
#' 
#' @export
approxDist<- function(Tsize,N2){
  S_N2 <- sim_moments(Tsize,N2)
  Fx <- approx_dist_loop(S_N2)
  x <- matrix(c(apply(matrix(S_N2[,1]),1,sort),apply(matrix(S_N2[,2]),1,sort),
                apply(matrix(S_N2[,3]),1,sort),apply(matrix(S_N2[,4]),1,sort)),
              ncol=4)
  # setting starting values
  a_start=0.01
  b_start=0.01
  # initiate matrix
  a<-matrix(nrow=1,ncol=0)
  b<-matrix(nrow=1,ncol=0)
  # estimate params of ecdf for each moment statistic
  for (i in 1:4){
    mdl<-nls( Fx[,i]~exp(alpha+beta*x[,i])/(1+exp(alpha+beta*x[,i])),
              start=list(alpha=a_start,beta=b_start))
    params<-coef(mdl)
    a<-cbind(a,params[1])
    b<-cbind(b,params[2])
  }
  return(rbind(matrix(a,nrow=1,ncol=4),matrix(b,nrow=1,ncol=4)))
}
# -----------------------------------------------------------------------------
#' @title Product of p-value maximizing function 
#' 
#' maximization function for nuissance parameters. This version uses min method 
#' of combining p-values as in  Fisher (1932) and Pearson (1933).
#'
#' @param v nuissance parameter values 
#' @return p-value for given parameter value.
#' 
#' @references Dufour, J. M., & Luger, R. (2017). Identification-robust moment-based 
#' tests for Markov switching in autoregressive models. Econometric Reviews, 36(6-9), 713-727.
#' @references Tippett, L. (1931). The Method of Statistics. London: Williams & Norgate.
#' @references Wilkinson, B. (1951). A statistical consideration in psychological research. 
#' Psychology Bulletin 48:156–158.
#' @references Pearson, K. (1933). On a method of determining whether a sample of size n
#'  supposed to have been drawn from a parent population having a known probability integral has probably been drawn at random. Biometrika 25:379–410.
#' @references Fisher, R. (1932). Statistical Methods for Research Workers. Edinburgh: 
#' Oliver and Boyd.
#' 
#' @export
max_fn_prod <- function(v){
  if (min(Mod(as.complex(polyroot(rbind(1,as.matrix(v))))))<1){
    pval  <-  0 
  }else{
    y     <- get('y',envir = MMCenv)
    x     <- get('x',envir = MMCenv)
    params<- get('params',envir = MMCenv)
    N     <- get('N',envir = MMCenv)
    N3    <- get('N3',envir = MMCenv)
    z     <- y - x %*% v
    s0    <- calc_moments(z-mean(z))
    sN    <- sim_moments(length(z),N-1)
    Fx    <- combine_stat(s0,sN,params,N,type="prod")
    pval  <- p_val(Fx[length(Fx)],Fx[1:(length(Fx)-1)],type = 'geq')
  }
  return(pval)
}
# -----------------------------------------------------------------------------
#' @title Minimum p-value maximizing function 
#' maximization function for nuissance parameters. This version uses product 
#' method of combining p-values as in Tippett (1931) and Wilkinson (1951).
#'
#' @param v nuissance parameter values 
#' @return p-value for given parameter value.
#' 
#' @references Dufour, J. M., & Luger, R. (2017). Identification-robust moment-based 
#' tests for Markov switching in autoregressive models. Econometric Reviews, 36(6-9), 713-727.
#' @references Tippett, L. (1931). The Method of Statistics. London: Williams & Norgate.
#' @references Wilkinson, B. (1951). A statistical consideration in psychological research. 
#' Psychology Bulletin 48:156–158.
#' @references Pearson, K. (1933). On a method of determining whether a sample of size n
#'  supposed to have been drawn from a parent population having a known probability integral has probably been drawn at random. Biometrika 25:379–410.
#' @references Fisher, R. (1932). Statistical Methods for Research Workers. Edinburgh: 
#' Oliver and Boyd.
#' 
#' @export
max_fn_min <- function(v){
  if (min(Mod(as.complex(polyroot(rbind(1,as.matrix(v))))))<1){
    pval  <-  0 
  }else{
    y     <- get('y',envir = MMCenv)
    x     <- get('x',envir = MMCenv)
    params<- get('params',envir = MMCenv)
    N     <- get('N',envir = MMCenv)
    N3    <- get('N3',envir = MMCenv)
    z     <- y - x %*% v
    s0    <- calc_moments(z-mean(z))
    sN    <- sim_moments(length(z),N-1)
    Fx    <- combine_stat(s0,sN,params,N,type="min")
    pval  <- p_val(Fx[length(Fx)],Fx[1:(length(Fx)-1)],type = 'geq')
  }
  return(pval)
}
# -----------------------------------------------------------------------------
#' @title Product of p-value minimizing function
#' 
#' minimization function for nuissance parameters. This version uses the product 
#' method of combining p-values as in Fisher (1932) and Pearson (1933).
#'
#' @param v nuissance parameter values 
#' @return p-value for given parameter value.
#' 
#' @references Dufour, J. M., & Luger, R. (2017). Identification-robust moment-based 
#' tests for Markov switching in autoregressive models. Econometric Reviews, 36(6-9), 713-727.
#' @references Tippett, L. (1931). The Method of Statistics. London: Williams & Norgate.
#' @references Wilkinson, B. (1951). A statistical consideration in psychological research. 
#' Psychology Bulletin 48:156–158.
#' @references Pearson, K. (1933). On a method of determining whether a sample of size n
#'  supposed to have been drawn from a parent population having a known probability integral has probably been drawn at random. Biometrika 25:379–410.
#' @references Fisher, R. (1932). Statistical Methods for Research Workers. Edinburgh: 
#' Oliver and Boyd.
#' 
#' @export
min_fn_prod <- function(v){
  n_pval <- -max_fn_prod(v)
  return(n_pval)
}
# -----------------------------------------------------------------------------
#' @title Minimum p-value minimizing min function
#' 
#' minimization function for nuissance parameters. This version uses min method 
#' of combining p-values as in Tippett (1931) and Wilkinson (1951).
#'
#' @param v nuissance parameter values 
#' @return p-value for given parameter value.
#' 
#' @references Dufour, J. M., & Luger, R. (2017). Identification-robust moment-based 
#' tests for Markov switching in autoregressive models. Econometric Reviews, 36(6-9), 713-727.
#' @references Tippett, L. (1931). The Method of Statistics. London: Williams & Norgate.
#' @references Wilkinson, B. (1951). A statistical consideration in psychological research. 
#' Psychology Bulletin 48:156–158.
#' @references Pearson, K. (1933). On a method of determining whether a sample of size n
#'  supposed to have been drawn from a parent population having a known probability integral has probably been drawn at random. Biometrika 25:379–410.
#' @references Fisher, R. (1932). Statistical Methods for Research Workers. Edinburgh: 
#' Oliver and Boyd.
#' 
#' @export
min_fn_min <- function(v){
  n_pval <- -max_fn_min(v)
  return(n_pval)
}