#' @title Companion Matrix
#'
#' @description This function converts the (q x 1)  vector of constants and (q x qp) matrix of autoregressive coefficients into (qp x qp) matrix belonging to the companion form
#'
#' @param phi matrix of dimension (q x qp) containing autoregressive coefficients
#' @param inter vector of dimension (q x 1) containing constants
#' @param p integer for number of autoregressive lags
#' @param q integer for number of series
#'
#' @return matrix of dimension (qp x qp) of companion form 
#' 
#' @export
companionMat <- function(phi, inter, p, q){
  interzero   <- matrix(0,q*(p-1), 1)
  comp_inter  <- as.matrix(c(c(inter),interzero))
  F_tmp       <- phi
  diagmat     <- diag(q*(p-1))
  diagzero    <- matrix(0, q*(p-1), q)
  Mn          <- cbind(diagmat,diagzero)
  compMat     <- rbind(F_tmp,Mn)
  return(compMat)
}



#' @title Autoregressive transition matrix
#'
#' @description This function converts a transition matrix to the transition matrix consistent with a Markov-switching autoregressive model.
#'
#' @param P original transition matrix
#' @param k integer determining the number of regimes
#' @param ar number of autoregressive lags
#'
#' @return transformed transition matrix
#' 
#' @export
arP <- function(P, k, ar){
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

#' @title Autoregressive moment grid
#'
#' @description This function creates a grid of mean and variance consistent with a Markov-switching autoregressive model.
#'
#' @param mu vector (k x 1) of mu in each regime
#' @param sig vector (k x 1) of sigma in each regime
#' @param k integer determining the number of regimes
#' @param ar number of autoregressive lags
#' @param msmu boolean indicator. If 'TRUE' mean is subject to change. If 'FALSE' mean is constant across regimes
#' @param msvar boolean indicator. If 'TRUE' variance is subject to change. If 'FALSE' variance is constant across regimes
#'
#' @return List with (M x ar+1) matrix of means for each regime M (where M = k^(ar+1)) and each time t,... t-ar, vector with variance for each regime M, and vector indicating the corresponded 1,..., k regime. 
#' 
#' @export
arGrid <- function(mu, sig, k, ar, msmu, msvar){
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


#' @title Vector autoregressive moment grid
#'
#' @description Creates grid of means and covariance matrices consistent  with a Markov-switching vector autoregressive model.
#'
#' @param mu (k x q) matrix of means in each regime (for k regimes and q time series)
#' @param sig list with k regime specific (q x q) covariance matrices 
#' @param k integer determining the number of regimes
#' @param ar number of autoregressive lags
#' @param msmu boolean indicator. If 'TRUE' mean is subject to change. If 'FALSE' mean is constant across regimes
#' @param msvar boolean indicator. If 'TRUE' variance is subject to change. If 'FALSE' variance is constant across regimes
#'
#' @return List with M regime specific (q x k) matrices of means, List with M regime specific covariance matrices, and vector indicating the corresponded 1,..., k regime. 
#' 
#' @export
varGrid <- function(mu, sigma, k, ar, msmu, msvar){
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

#' @title Log likelihood  
#' 
#' @description This function is used to compute the log-likelihood for a given model
#'
#' @param mdl List with model properties
#'
#' @return Log-likelihood
#' 
#' @export
logLikelihood <- function(mdl){
  UseMethod("logLikelihood", mdl)
}

#' @title Log likelihood for autoregressive model  
#' 
#' @description This function is used to compute the log-likelihood for an autoregressive model
#'
#' @param mdl List with model properties
#'
#' @return Log-likelihood
#' 
#' @export
logLikelihood.ARmdl <- function(mdl){
  logLike <- logLike_AR(mdl$theta, mdl)
  return(logLike)
}

#' @title Log likelihood for vector autoregressive model  
#' 
#' @description This function is used to compute the log-likelihood for a vector autoregressive model
#'
#' @param mdl List with model properties
#'
#' @return Log-likelihood
#' 
#' @export
logLikelihood.VARmdl <- function(mdl){
  logLike <- logLike_VAR(mdl$theta, mdl)
  return(logLike)
}


#' @title Hessian matrix 
#' 
#' @description This function is used to obtain a numerical approximation of a Hessian matrix
#'
#' @param mdl List with model properties
#'
#' @return Hessian matrix
#' 
#' @export
getHessian <- function(mdl){
  UseMethod("getHessian", mdl)
}

#' @title Hessian matrix of autoregressive model
#' 
#' @description This function is used to obtain a numerical approximation of a Hessian matrix for an autoregressive model
#'
#' @param mdl List with model properties
#'
#' @return Hessian matrix
#' 
#' @export
getHessian.ARmdl <- function(mdl){
  hess <- numDeriv::hessian(logLike_AR, mdl$theta, method = "Richardson", mdl = mdl) 
  return(hess)
}

#' @title Hessian matrix of vector autoregressive model 
#' 
#' @description This function is used to obtain a numerical approximation of a Hessian matrix for a vector autoregressive model
#'
#' @param mdl List with model properties
#'
#' @return Hessian matrix
#' 
#' @export
getHessian.VARmdl <- function(mdl){
  hess <- numDeriv::hessian(logLike_VAR, mdl$theta, method = "Richardson", mdl = mdl) 
  return(hess)
}

#' @title Hessian matrix of Markov-switching autoregressive model
#' 
#' @description This function is used to obtain a numerical approximation of a Hessian matrix for a Markov-switching autoregressive model
#'
#' @param mdl List with model properties
#'
#' @return Hessian matrix
#' 
#' @export
getHessian.MSARmdl <- function(mdl){
  hess <- numDeriv::hessian(MSloglik_fun, mdl$theta, method = "Richardson", mdl = mdl, k = mdl$k) 
  return(hess)
}

#' @title Hessian matrix of Markov-switching vector autoregressive
#' 
#' @description This function is used to obtain a numerical approximation of a Hessian matrix for a Markov-switching vector autoregressive model
#'
#' @param mdl List with model properties
#'
#' @return Hessian matrix
#' 
#' @export
getHessian.MSVARmdl <- function(mdl){
  hess <- numDeriv::hessian(MSVARloglik_fun, mdl$theta, method = "Richardson", mdl = mdl, k = mdl$k) 
  return(hess)
}



