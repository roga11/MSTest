

# ==============================================================================
#' @title Convert transition matrix to transition matrix consistent with AR model
#'
#' @param P original transition matrix
#' @param k number of regimes
#' @param ar number of autoregressive lags
#'
#' @return transformed tansition matrix
#' 
#' @references Hamilton (1994)
#' 
#' @export
transMatAR <- function(P, k, ar){
  if (ar>0){
    Pnew_tmp <- matrix(0,k*k,k)
    for (xk in 1:k){
      ind <- rep(FALSE,k)
      ind[xk] <- TRUE
      Pnew_tmp[(1+(xk-1)*k):(1+(xk-1)*k+(k-1)),xk] <- P[,ind]
    }
    repmat <- k^(ar-1)
    Pnew <- do.call(cbind,replicate(k,kronecker(diag(repmat), Pnew_tmp),simplify=FALSE)) 
  }else{
    Pnew <- P
  }
  return(Pnew)  
}
# ==============================================================================
#' @title Mu and Sigma AR grid
#'
#' @description Creates grid of mu and sigma consistent with number of Autoregressive lags
#'
#' @param mu vector (k x 1) of mu in each regime
#' @param sig vector (k x 1) of sigma in each regime
#' @param k number of regimes
#' @param ar number of autoregressive lags
#'
#' @return 
#' 
#' @references Hamilton (1994) see pg. 691
#' 
#' @export
musigGrid <- function(mu, sig, k, ar, msmu, msvar){
  if (msmu==FALSE){
    mu <- rep(mu, k)
  }
  if (msvar==FALSE){
    sig <- rep(sig, k)
  }
  mu <- as.vector(mu)
  sig <- as.vector(sig)
  # create grid of regimes
  lx <-list()
  for (xi in 1:(ar+1)) lx[[xi]] <- 1:k
  mu_stategrid <- as.matrix(expand.grid(lx))
  sig_stategrid <- as.matrix(expand.grid(lx))
  state_indicator = mu_stategrid[,1]
  # get matrix of mu and sigma
  mu_stategrid[] <- mu[mu_stategrid]
  sig_stategrid[] <- sig[sig_stategrid]
  sig_stategrid <- sig_stategrid[,1] # sig doesnt have AR seq
  musig_out <- list()
  musig_out[["mu"]] <- as.matrix(mu_stategrid)
  musig_out[["sig"]] <- as.matrix(sig_stategrid)
  musig_out[["state_ind"]] <- as.matrix(state_indicator)
  return(musig_out)
}
# ==============================================================================
#' @title Mu and Sigma VAR grid
#'
#' @description Creates grid of mu and sigma consistent with number of Autoregressive lags
#'
#' @param mu vector (k x N) of mu in each regime
#' @param sig list of k regime specific (N x N) matrices 
#' @param k number of regimes
#' @param ar number of autoregressive lags
#'
#' @return 
#' 
#' 
#' @export
musigVARGrid <- function(mu, sigma, k, ar, msmu, msvar){
  # add if msmu and msvar ==0
  #if (msmu==FALSE){
  #  N <- length(mu)
  #  mu_tmp <- mu
  #  mu <- matrix(0,k,N)
  #  for (xk in 1:k){
  #    mu[xk,] <- mu_tmp
  #  }
    #mu <- t(matrix(rep(mu,k),N,k))
  #}
  #if (msvar==FALSE){
  #  for (xk in 2:k){
  #   sigma[[xk]] <- sigma[[1]]
  #  }
  #}
  # create grid of regimes
  lx <-list()
  for (xi in 1:(ar+1)) lx[[xi]] <- 1:k
  mu_stategrid <- as.matrix(expand.grid(lx))
  state_indicator <- mu_stategrid[,1]
  M <- k^(ar+1)
  sig_stategrid <- list()
  muN_stategrid <- list()
  for (xm in 1:M){
    sig_stategrid[[xm]] <- sigma[[state_indicator[xm]]]
    muN_stategrid[[xm]] <- t(mu[mu_stategrid[xm,],])
  }
  musig_out <- list()
  musig_out[["mu"]] <- muN_stategrid
  musig_out[["sig"]] <- sig_stategrid
  musig_out[["state_ind"]] <- as.matrix(state_indicator)
  return(musig_out)
}
# ==============================================================================
#' @title Obtain Hessian matrix 
#' 
#' @description Use numDeriv() package to obtain Hessian matrix which is then used to get standard errors of parameters
#'
#' @param EMmdl_out List with model output from EM algorithm
#' @param k int number of regimes
#' @param msmu bool indicating if mean is switching with regime
#' @param msvar bool indicating if variance is switching with regime
#'
#' @return Hessian matrix
#' 
#' @export
getHess <- function(EMmdl_out, k){
  theta_EM = EMmdl_out$theta
  hess <- numDeriv::hessian(MSloglik_fun, theta_EM, method = "Richardson", mdl = EMmdl_out, k = k)
  return(hess)
}
# ==============================================================================
#' @title Approximate Distribuion For Dufour & Luger Test
#'
#' @description This function obtains the parameters needed in eq. 16 which is used for 
#' combining p-values.
#'
#' @param Tsize sample size
#' @param simdist_N number of draws (simulations)
#' 
#' @return params the paramters gamma in eq. 16 of Dufour & Luger (2017). 
#' 
#' @references Dufour, J. M., & Luger, R. (2017). Identification-robust moment-based 
#' tests for Markov switching in autoregressive models. Econometric Reviews, 36(6-9), 713-727.
#' 
#' @export
approxDistDL<- function(Tsize, simdist_N){
  S_N2  <- sim_DLmoments(Tsize, simdist_N)
  x     <- apply(S_N2, 2, sort)
  Fx    <- approx_dist_loop(x)
  a_start <- 0.01
  b_start <- 0.01
  # initiate matrix
  a<-matrix(nrow=1,ncol=0)
  b<-matrix(nrow=1,ncol=0)
  # estimate params of ecdf for each moment statistic
  for (i in 1:4){
    mdl<-nls(Fx[,i]~exp(alpha+beta*x[,i])/(1+exp(alpha+beta*x[,i])),
              start=list(alpha=a_start,beta=b_start))
    params<-coef(mdl)
    a<-cbind(a,params[1])
    b<-cbind(b,params[2])
  }
  return(rbind(matrix(a,nrow=1,ncol=4),matrix(b,nrow=1,ncol=4)))
}
# ==============================================================================
#' @title Approximate Distribuion For Quantile-Based Test
#'
#' @description This function obtains the parameters needed in eq. 16 which is used for 
#' combining p-values.
#'
#' @param Tsize sample size
#' @param simdist_N number of draws (simulations)
#' 
#' @return params the paramters gamma in eq. 16 of Dufour & Luger (2017). 
#' 
#' @references Dufour, J. M., & Luger, R. (2017). Identification-robust moment-based 
#' tests for Markov switching in autoregressive models. Econometric Reviews, 36(6-9), 713-727.
#' 
#' @export
approxDistQ <- function(Tsize, simdist_N, mdl_h0, k1){
  S_N2  <- sim_Qmoments(Tsize, simdist_N, mdl_h0, k1)
  x     <- apply(S_N2, 2, sort)
  Fx    <- approx_dist_loop(x)
  # setting starting values
  a_start <- 0.01
  b_start <- 0.01
  # initiate matrix
  a <- matrix(nrow=1,ncol=0)
  b <- matrix(nrow=1,ncol=0)
  # estimate params of ecdf for each moment statistic
  for (xi in 1:ncol(x)){
    mdl<-nls(Fx[,xi]~exp(alpha+beta*x[,xi])/(1+exp(alpha+beta*x[,xi])),
              start=list(alpha=a_start,beta=b_start))
    params<-coef(mdl)
    a <- cbind(a,params[1])
    b <- cbind(b,params[2])
  }
  return(rbind(matrix(a,nrow=1,ncol=ncol(x)),matrix(b,nrow=1,ncol=ncol(x))))
}
# ==============================================================================
#' @title Estimate Markov-Switching Model using Constrained Optimization
#' 
#' @description 
#'
#'
#' @return 
#' 
#' @export
MSARmdl_optim <- function(Y, ar, k, msmu =TRUE , msvar = TRUE, maxit = 10000, thtol = 1e-6){
  mdl <- ARmdl(Y, ar, TRUE)
  mdl$msmu <- msmu
  mdl$msvar <- msvar
  mu = mdl$mu
  stdev = mdl$stdev
  mu_0 = mu
  sig_0 = stdev^2
  if (msmu==TRUE){
    mu_0 = c(mu_0, mu_0+rnorm(k-1))
  }
  if (msvar==TRUE){
    sig_0 = c(sig_0,sig_0+runif(k-1))
  }
  theta_0 = c(mu_0, sig_0)
  if (ar>0){
    theta_0 = c(theta_0, mdl$phi)
  }
  P = randTransMat(k, mdl$n)
  theta_0 = c(theta_0, c(P))
  # ----- Optimization bounds
  theta_lowr = c(rep(mu-(stdev*10), 1+msmu*(k-1)), rep(0.000000001,1+msvar*(k-1)), rep(-1, ar), rep(0.000000001, k*k))
  theta_uppr = c(rep(mu+(stdev*10), 1+msmu*(k-1)), rep(stdev^2+(stdev^2*1000), 1+msvar*(k-1)), rep(1, ar), rep(0.999999999, k*k))
  # ----- Optimization options
  opts <- list( "algorithm"= "NLOPT_GN_ISRES",
                "xtol_rel"= thtol,
                "maxeval"= maxit)
  # Optimization
  res <- nloptr::nloptr( x0 = theta_0,
                         eval_f = MSloglik_fun_min,
                         lb = theta_lowr,
                         ub = theta_uppr,
                         eval_g_ineq = MSloglik_const_ineq,
                         eval_g_eq = MSloglik_const_eq,
                         opts = opts , 
                         mdl = mdl,
                         k = k)
  return(res)
}

