#' @title Companion Matrix
#'
#' @description This function converts the (\code{q x 1})  vector of constants and (\code{q x qp}) matrix of autoregressive coefficients into (\code{qp x qp}) matrix belonging to the companion form
#'
#' @param phi matrix of dimension (\code{q x qp}) containing autoregressive coefficients.
#' @param p integer for number of autoregressive lags.
#' @param q integer for number of series.
#'
#' @return matrix of dimension (\code{qp x qp}) of companion form.
#' 
#' @keywords internal
#' 
#' @export
companionMat <- function(phi, p, q){
  interzero   <- matrix(0,q*(p-1), 1)
  F_tmp       <- phi
  diagmat     <- diag(q*(p-1))
  diagzero    <- matrix(0, q*(p-1), q)
  Mn          <- cbind(diagmat,diagzero)
  compMat     <- rbind(F_tmp,Mn)
  return(compMat)
}

#' @title Compute AIC
#'
#' @description This function computes the AIC when given the log-likelihood and number of parameters.
#'
#' @param logLike log-likelihood of estimated model.
#' @param param_len number of parameters in model.
#'
#' @return AIC value. 
#' 
#' @keywords internal
#' 
#' @export
aic <- function(logLike, param_len){
  aic_val <- -2*(logLike) + 2*param_len
  return(aic_val)
}

#' @title Compute BIC
#'
#' @description This function computes the BIC when given the log-likelihood, number of observations, and number of parameters.
#'
#' @param logLike log-likelihood of estimated model.
#' @param n number of time series observations.
#' @param param_len number of parameters in model.
#'
#' @return BIC value.
#' 
#' @keywords internal
#' 
#' @export
bic <- function(logLike, n, param_len){
  bic_val <- -param_len*log(n) -2*(logLike) 
  return(bic_val)
}

#' @title Autoregressive transition matrix
#'
#' @description This function converts a transition matrix to the transition matrix consistent with a Markov-switching autoregressive model.
#'
#' @param P original transition matrix.
#' @param k integer determining the number of regimes.
#' @param ar number of autoregressive lags.
#'
#' @return transformed transition matrix.
#' 
#' @keywords internal
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
#' @param mu vector (\code{k x 1}) of mu in each regime.
#' @param sig vector (\code{k x 1}) of sigma in each regime.
#' @param k integer determining the number of regimes.
#' @param ar number of autoregressive lags.
#' @param msmu Boolean indicator. If \code{TRUE} mean is subject to change. If \code{FALSE} mean is constant across regimes.
#' @param msvar Boolean indicator. If \code{TRUE} variance is subject to change. If \code{FALSE} variance is constant across regimes.
#'
#' @return List with (\code{M x ar+1}) matrix of means for each regime \code{M} (where \code{M = k^(ar+1)}) and each time \code{t,... t-ar}, vector with variance for each regime \code{M}, and vector indicating the corresponded \code{1,..., k} regime. 
#' 
#' @keywords internal
#'
#' @export
argrid_MSARmdl <- function(mu, sig, k, ar, msmu, msvar){
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
#' @description Creates grid of means and covariance matrices consistent with a Markov-switching vector autoregressive model.
#'
#' @param mu a (\code{k x q}) matrix of means in each regime (for \code{k} regimes and \code{q} time series).
#' @param sigma list with \code{k} regime specific (\code{q x q}) covariance matrices.
#' @param k integer determining the number of regimes.
#' @param ar number of autoregressive lags.
#' @param msmu Boolean indicator. If \code{TRUE} mean is subject to change. If \code{FALSE} mean is constant across regimes.
#' @param msvar Boolean indicator. If \code{TRUE} variance is subject to change. If \code{FALSE} variance is constant across regimes.
#'
#' @return List with M regime specific (\code{q x k}) matrices of means, List with \code{M} regime specific covariance matrices, and vector indicating the corresponded \code{1,..., k} regime. 
#' 
#' @keywords internal
#' 
#' @export
argrid_MSVARmdl <- function(mu, sigma, k, ar, msmu, msvar){
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
#' @description This function is used to compute the log-likelihood for a given model.
#'
#' @param mdl List with model properties.
#'
#' @return Log-likelihood value.
#' 
#' @keywords internal
#' 
#' @export
logLikelihood <- function(mdl){
  UseMethod("logLikelihood", mdl)
}

#' @title Log likelihood for Normal model  
#' 
#' @description This function is used to compute the log-likelihood for a normally distributed model.
#'
#' @param mdl List with model properties.
#'
#' @return Log-likelihood value.
#' 
#' @keywords internal
#' 
#' @export
logLikelihood.Nmdl <- function(mdl){
  logLike <- logLike_Nmdl(mdl$theta, mdl)
  return(logLike)
}

#' @title Log likelihood for autoregressive model  
#' 
#' @description This function is used to compute the log-likelihood for an autoregressive model.
#'
#' @param mdl List with model properties.
#'
#' @return Log-likelihood value.
#' 
#' @keywords internal
#' 
#' @export
logLikelihood.ARmdl <- function(mdl){
  logLike <- logLike_ARmdl(mdl$theta, mdl)
  return(logLike)
}

#' @title Log likelihood for vector autoregressive model  
#' 
#' @description This function is used to compute the log-likelihood for a vector autoregressive model.
#'
#' @param mdl List with model properties.
#'
#' @return Log-likelihood value.
#' 
#' @keywords internal
#' 
#' @export
logLikelihood.VARmdl <- function(mdl){
  logLike <- logLike_VARmdl(mdl$theta, mdl)
  return(logLike)
}

#' @title Log likelihood for Hidden Markov model  
#' 
#' @description This function is used to compute the log-likelihood for a Hidden Markov model.
#'
#' @param mdl List with model properties.
#'
#' @return Log-likelihood value. 
#' 
#' @keywords internal
#' 
#' @export
logLikelihood.HMmdl <- function(mdl){
  logLike <- logLike_HMmdl(mdl$theta, mdl, mdl$k)
  return(logLike)
}

#' @title Log likelihood for Markov-switching autoregressive model  
#' 
#' @description This function is used to compute the log-likelihood for a Markov-switching autoregressive model.
#'
#' @param mdl List with model properties.
#'
#' @return Log-likelihood value.
#' 
#' @keywords internal
#' 
#' @export
logLikelihood.MSARmdl <- function(mdl){
  logLike <- logLike_MSARmdl(mdl$theta, mdl, mdl$k)
  return(logLike)
}

#' @title Log likelihood for Markov-switching vector autoregressive model  
#' 
#' @description This function is used to compute the log-likelihood for a Markov-switching vector autoregressive model.
#'
#' @param mdl List with model properties.
#'
#' @return Log-likelihood value.
#' 
#' @keywords internal
#' 
#' @export
logLikelihood.MSVARmdl <- function(mdl){
  logLike <- logLike_MSVARmdl(mdl$theta, mdl, mdl$k)
  return(logLike)
}

#' @title Hessian matrix 
#' 
#' @description This function is used to obtain a numerical approximation of a Hessian matrix.
#'
#' @param mdl List with model properties.
#'
#' @return Hessian matrix.
#' 
#' @keywords internal
#' 
#' @export
getHessian <- function(mdl){
  UseMethod("getHessian", mdl)
}

#' @title Hessian matrix of normal model
#' 
#' @description This function is used to obtain a numerical approximation of a Hessian matrix for a normally distributed model.
#'
#' @param mdl List with model properties.
#'
#' @return Hessian matrix.
#' 
#' @keywords internal
#' 
#' @export
getHessian.Nmdl <- function(mdl){
  hess <- numDeriv::hessian(logLike_Nmdl, mdl$theta, method = "Richardson", mdl = mdl) 
  return(hess)
}

#' @title Hessian matrix of autoregressive model
#' 
#' @description This function is used to obtain a numerical approximation of a Hessian matrix for an autoregressive model.
#'
#' @param mdl List with model properties.
#'
#' @return Hessian matrix.
#' 
#' @keywords internal
#' 
#' @export
getHessian.ARmdl <- function(mdl){
  hess <- numDeriv::hessian(logLike_ARmdl, mdl$theta, method = "Richardson", mdl = mdl) 
  return(hess)
}

#' @title Hessian matrix of vector autoregressive model 
#' 
#' @description This function is used to obtain a numerical approximation of a Hessian matrix for a vector autoregressive model.
#'
#' @param mdl List with model properties.
#'
#' @return Hessian matrix.
#' 
#' @keywords internal
#' 
#' @export
getHessian.VARmdl <- function(mdl){
  hess <- numDeriv::hessian(logLike_VARmdl, mdl$theta, method = "Richardson", mdl = mdl) 
  return(hess)
}

#' @title Hessian matrix of Hidden Markov model
#' 
#' @description This function is used to obtain a numerical approximation of a Hessian matrix for a Hidden Markov model.
#'
#' @param mdl List with model properties.
#'
#' @return Hessian matrix.
#' 
#' @keywords internal
#' 
#' @export
getHessian.HMmdl <- function(mdl){
  hess <- numDeriv::hessian(logLike_HMmdl, mdl$theta, method = "Richardson", mdl = mdl, k = mdl$k) 
  return(hess)
}

#' @title Hessian matrix of Markov-switching autoregressive model
#' 
#' @description This function is used to obtain a numerical approximation of a Hessian matrix for a Markov-switching autoregressive model.
#'
#' @param mdl List with model properties.
#'
#' @return Hessian matrix.
#' 
#' @keywords internal
#' 
#' @export
getHessian.MSARmdl <- function(mdl){
  hess <- numDeriv::hessian(logLike_MSARmdl, mdl$theta, method = "Richardson", mdl = mdl, k = mdl$k) 
  return(hess)
}

#' @title Hessian matrix of Markov-switching vector autoregressive model
#' 
#' @description This function is used to obtain a numerical approximation of a Hessian matrix for a Markov-switching vector autoregressive model.
#'
#' @param mdl List with model properties.
#'
#' @return Hessian matrix.
#' 
#' @keywords internal
#' 
#' @export
getHessian.MSVARmdl <- function(mdl){
  hess <- numDeriv::hessian(logLike_MSVARmdl, mdl$theta, method = "Richardson", mdl = mdl, k = mdl$k) 
  return(hess)
}


#' @title Theta standard errors
#' 
#' @description This function computes the standard errors of the parameters in vector theta. This is done using an approximation of the Hessian matrix (using \code{\link[numDeriv]{hessian}} and \code{nearPD} if \code{info_mat} is not PD).
#'
#' @param mdl List with model properties
#'
#' @return List provided as input with additional attributes \code{HESS},\code{theta_se}, \code{info_mat}, and \code{nearPD_used}.
#' 
#' @keywords internal
#' 
#' @export
thetaSE <- function(mdl){
  Hess <- getHessian(mdl)
  info_mat <- solve(-Hess)
  nearPD_used <- FALSE
  if ((all(is.na(Hess)==FALSE)) & (any(diag(info_mat)<0))){
    info_mat <- lmf::nearPD(info_mat)
    nearPD_used <- TRUE
  }
  mdl$Hess <- Hess
  mdl$theta_se <- sqrt(diag(info_mat))
  mdl$info_mat <- info_mat
  mdl$nearPD_used <- nearPD_used 
  return(mdl)
}

#' @title Hidden Markov model maximum likelihood estimation
#' 
#' @description This function computes estimate a Hidden Markov model using MLE.
#' 
#' @param theta_0 vector containing initial values to use in optimization
#' @param mdl_in List with model properties (can be obtained from estimating linear model i.e., using \code{\link{ARmdl}})
#' @param k integer determining the number of regimes
#' @param optim_options List containing 
#' \itemize{
#'  \item{\code{maxit}: }{maximum number of iterations.}
#'  \item{\code{thtol}: }{convergence criterion.}
#'}
#'
#' @return List with model attributes
#' 
#' @keywords internal
#' 
#' @export
HMmdl_mle <- function(theta_0, mdl_in, k, optim_options){
  # ----- Define function environment
  HMmdl_mle_env <- rlang::env()
  # set environment variables
  HMmdl_mle_env$k <- k
  HMmdl_mle_env$q <- mdl_in$q
  HMmdl_mle_env$msmu <- mdl_in$msmu
  HMmdl_mle_env$msvar <- mdl_in$msvar
  HMmdl_mle_env$var_k1 <- mdl_in$sigma
  HMmdl_mle_env$mle_variance_constraint <- mdl_in$mle_variance_constraint
  # ----------- HMmdl_mle equality constraint functions
  loglik_const_eq_HMmdl <- function(theta){
    # ----- Load equality constraint parameters
    k <- get("k", envir = HMmdl_mle_env)
    # ----- constraint
    P = matrix(theta[(length(theta)-k*k+1):(length(theta))],k,k)
    constraint = colSums(P)-1
    return(constraint)
  }
  # ---------- HMmdl_mle inequality constraint functions
  loglik_const_ineq_HMmdl <- function(theta){
    # ----- Load inequality constraint parameters
    k <- get("k", envir = HMmdl_mle_env)
    q <- get("q", envir = HMmdl_mle_env)
    msmu <- get("msmu", envir = HMmdl_mle_env)
    msvar <- get("msvar", envir = HMmdl_mle_env)
    var_k1 <- get("var_k1", envir = HMmdl_mle_env)
    mle_variance_constraint <- get("mle_variance_constraint", envir = HMmdl_mle_env)
    # ----- constraint
    # transition probabilities
    Pvec <- theta[(length(theta)-k*k+1):(length(theta))]
    ineq_constraint <- c(Pvec, 1-Pvec)
    if (mle_variance_constraint>=0){
      # variances (lower bound is 1% of single regime variance)
      #vars <- theta[c(rep(0, 1 + (k-1)*msmu), rep(1, 1 + (k-1)*msvar), rep(0, k*k))==1] - mle_variance_constraint*var_k1
      #ineq_constraint <- c(vars, ineq_constraint)
      Nsig <- (q*(q+1))/2
      eigen_vals <- c()
      sig <- theta[c(rep(0, q + q*(k-1)*msmu), rep(1, Nsig + Nsig*(k-1)*msvar))==1]
      if (msvar==TRUE){
        for (xk in  1:k){
          sig_tmp <- sig[(Nsig*(xk-1)+1):(Nsig*(xk-1)+Nsig)]
          sigma <- covar_unvech(sig_tmp, q)
          eigen_vals <- c(eigen_vals,eigen(sigma)[[1]])
        } 
      }else{
        sigma <- covar_unvech(sig, q)
        eigen_vals <- c(eigen_vals,eigen(sigma)[[1]])
      }
      ineq_constraint = c(eigen_vals-mle_variance_constraint, ineq_constraint)
    }
    return(ineq_constraint)
  }
  # use nloptr optimization to minimize (maximize) likelihood
  res <- nloptr::slsqp(x0 = theta_0,
                       fn = logLike_HMmdl_min,
                       gr = NULL,
                       lower = NULL,
                       upper = NULL,
                       hin = loglik_const_ineq_HMmdl,
                       heq = loglik_const_eq_HMmdl,
                       nl.info = FALSE,
                       control = list(maxeval = optim_options$maxit, xtol_rel = optim_options$thtol),
                       mdl = mdl_in,
                       k = k) 
  output <- ExpectationM_HMmdl(res$par, mdl_in, k)
  output$iterations <- res$iter
  output$St <- output$xi_t_T
  return(output)
}

#' @title Markov-switching autoregressive maximum likelihood estimation
#' 
#' @description This function computes estimate a Markov-switching autoregressive model using MLE.
#' 
#' @param theta_0 vector containing initial values to use in optimization
#' @param mdl_in List with model properties (can be obtained from estimating linear model i.e., using \code{\link{ARmdl}})
#' @param k integer determining the number of regimes
#' @param optim_options List containing 
#' \itemize{
#'  \item{\code{maxit}: }{maximum number of iterations.}
#'  \item{\code{thtol}: }{convergence criterion.}
#'}
#'
#' @return List with model attributes
#' 
#' @keywords internal
#' 
#' @export
MSARmdl_mle <- function(theta_0, mdl_in, k, optim_options){
  # ----- Define function environment
  MSARmdl_mle_env <- rlang::env()
  # set environment variables
  MSARmdl_mle_env$p <- mdl_in$p
  MSARmdl_mle_env$k <- k
  MSARmdl_mle_env$msmu <- mdl_in$msmu
  MSARmdl_mle_env$msvar <- mdl_in$msvar
  MSARmdl_mle_env$var_k1 <- mdl_in$sigma
  MSARmdl_mle_env$mle_stationary_constraint <- mdl_in$mle_stationary_constraint
  MSARmdl_mle_env$mle_variance_constraint <- mdl_in$mle_variance_constraint
  # ----------- MSARmdl_mle equality constraint functions
  loglik_const_eq_MSARmdl <- function(theta){
    # ----- Load equality constraint parameters
    k <- get("k", envir = MSARmdl_mle_env)
    # ----- constraint
    P = matrix(theta[(length(theta)-k*k+1):(length(theta))],k,k)
    constraint = colSums(P)-1
    return(constraint)
  }
  # ---------- MSARmdl_mle inequality constraint functions
  loglik_const_ineq_MSARmdl <- function(theta){
    # ----- Load inequality constraint parameters
    p <- get("p", envir = MSARmdl_mle_env)
    k <- get("k", envir = MSARmdl_mle_env)
    msmu <- get("msmu", envir = MSARmdl_mle_env)
    msvar <- get("msvar", envir = MSARmdl_mle_env)
    var_k1 <- get("var_k1", envir = MSARmdl_mle_env)
    mle_stationary_constraint <- get("mle_stationary_constraint", envir = MSARmdl_mle_env)
    mle_variance_constraint <- get("mle_variance_constraint", envir = MSARmdl_mle_env)
    # ----- constraint
    # transition probabilities
    Pvec <- theta[(length(theta)-k*k+1):(length(theta))]
    ineq_constraint <- c(Pvec, 1-Pvec)
    if (mle_stationary_constraint==TRUE){
      # roots of characteristic function
      n_p <- 3 + msmu*(k-1)+msvar*(k-1)
      poly_fun <- c(1,-theta[n_p:(n_p+p-1)])
      roots <- Mod(polyroot(poly_fun))
      ineq_constraint = c((roots-1), ineq_constraint)
    }
    if (mle_variance_constraint>=0){
      # variances (lower bound is 1% of single regime variance)
      vars <- theta[c(rep(0, 1 + (k-1)*msmu), rep(1, 1 + (k-1)*msvar), rep(0, p + k*k))==1] - mle_variance_constraint*var_k1
      ineq_constraint <- c(vars, ineq_constraint)
    }
    return(ineq_constraint)
  }
  # use nloptr optimization to minimize (maximize) likelihood
  res <- nloptr::slsqp(x0 = theta_0,
                       fn = logLike_MSARmdl_min,
                       gr = NULL,
                       lower = NULL,
                       upper = NULL,
                       hin = loglik_const_ineq_MSARmdl,
                       heq = loglik_const_eq_MSARmdl,
                       nl.info = FALSE,
                       control = list(maxeval = optim_options$maxit, xtol_rel = optim_options$thtol),
                       mdl = mdl_in,
                       k = k) 
  output <- ExpectationM_MSARmdl(res$par, mdl_in, k)
  output$iterations <- res$iter
  output$St <- output$xi_t_T
  return(output)
}


#' @title Markov-switching vector autoregressive maximum likelihood estimation
#' 
#' @description This function computes estimate a Markov-switching vector autoregressive model using MLE.
#' 
#' @param theta_0 vector containing initial values to use in optimization
#' @param mdl_in List with model properties (can be obtained from estimating linear model i.e., using \code{\link{VARmdl}})
#' @param k integer determining the number of regimes
#' @param optim_options List containing 
#' \itemize{
#'  \item{\code{maxit}: }{maximum number of iterations.}
#'  \item{\code{thtol}: }{convergence criterion.}
#'}
#'
#' @return List with model attributes
#' 
#' @keywords internal
#' 
#' @export
MSVARmdl_mle <- function(theta_0, mdl_in, k, optim_options){
  # ---------- Define function environment
  MSVARmdl_mle_env <- rlang::env()
  # set environment variables
  MSVARmdl_mle_env$p <- mdl_in$p
  MSVARmdl_mle_env$k <- k
  MSVARmdl_mle_env$q <- mdl_in$q
  MSVARmdl_mle_env$msmu <- mdl_in$msmu
  MSVARmdl_mle_env$msvar <- mdl_in$msvar
  MSVARmdl_mle_env$mle_stationary_constraint <- mdl_in$mle_stationary_constraint
  MSVARmdl_mle_env$mle_variance_constraint <- mdl_in$mle_variance_constraint
  # ---------- MSVARmdl_mle equality constraint functions
  loglik_const_eq_MSVARmdl <- function(theta){
    # ----- Load equality constraint parameters
    k <- get("k", envir = MSVARmdl_mle_env)
    # ----- constraint
    P <- matrix(theta[(length(theta)-k*k+1):(length(theta))],k,k)
    constraint <- colSums(P)-1
    return(constraint)
  }
  # ---------- MSVARmdl_mle inequality constraint functions
  loglik_const_ineq_MSVARmdl <- function(theta){
    # ----- Load inequality constraint parameters
    p <- get("p", envir = MSVARmdl_mle_env)
    k <- get("k", envir = MSVARmdl_mle_env)
    q <- get("q", envir = MSVARmdl_mle_env)
    msmu <- get("msmu", envir = MSVARmdl_mle_env)
    msvar <- get("msvar", envir = MSVARmdl_mle_env)
    mle_stationary_constraint <- get("mle_stationary_constraint", envir = MSVARmdl_mle_env)
    mle_variance_constraint <- get("mle_variance_constraint", envir = MSVARmdl_mle_env)
    Nsig <- (q*(q+1))/2
    phi_len <- q*p*q
    # ----- constraint 
    Pvec = theta[(length(theta)-k*k+1):(length(theta))]
    ineq_constraint = c(Pvec, 1-Pvec)
    if (mle_stationary_constraint==TRUE){
      # eigen values of companion matrix
      phi <- t(matrix(theta[c(rep(0, q + q*(k-1)*msmu + Nsig + Nsig*(k-1)*msvar), rep(1, phi_len), rep(0, k*k))==1], q*p, q))
      Fmat <- companionMat(phi,p,q)
      stationary  <- abs(Mod(eigen(Fmat)[[1]]) - 1)
      ineq_constraint = c(stationary, 1-stationary, ineq_constraint)
    }
    if (mle_variance_constraint>=0){
      # eigen values of covariance matrices
      eigen_vals <- c()
      sig <- theta[c(rep(0, q + q*(k-1)*msmu), rep(1, Nsig + Nsig*(k-1)*msvar), rep(0, phi_len + k*k))==1]
      if (msvar==TRUE){
        for (xk in  1:k){
          sig_tmp <- sig[(Nsig*(xk-1)+1):(Nsig*(xk-1)+Nsig)]
          sigma <- covar_unvech(sig_tmp, q)
          eigen_vals <- c(eigen_vals,eigen(sigma)[[1]])
        } 
      }else{
        sigma <- covar_unvech(sig, q)
        eigen_vals <- c(eigen_vals,eigen(sigma)[[1]])
      }
      ineq_constraint = c(eigen_vals-mle_variance_constraint, ineq_constraint)
    }
    return(ineq_constraint)
  }
  # use nloptr optimization to minimize (maximize) likelihood
  res <- nloptr::slsqp(x0 = theta_0,
                       fn = logLike_MSVARmdl_min,
                       gr = NULL,
                       lower = NULL,
                       upper = NULL,
                       hin = loglik_const_ineq_MSVARmdl,
                       heq = loglik_const_eq_MSVARmdl,
                       nl.info = FALSE,
                       control = list(maxeval = optim_options$maxit, xtol_rel = optim_options$thtol),
                       mdl = mdl_in,
                       k = k) 
  output <- ExpectationM_MSVARmdl(res$par, mdl_in, k)
  output$iterations <- res$iter
  output$St <- output$xi_t_T
  return(output)
}


#' @title Print summary of a \code{Nmdl} object
#'
#' @description This is a method for the function \code{print()} for objects of the class \code{Nmdl}.
#' 
#' @inheritParams base::print
#'
#' @return The \code{Nmdl} object is returned invisibly.
#' 
#' @keywords internal
#' 
#' @export
print.Nmdl <- function(x, digits = getOption("digits"), ...){
  cat("\nNormally Distributed Model\n")
  frame_tmp <- data.frame(coef = x$theta)
  if (x$control$getSE==TRUE){
    frame_tmp["s.e."] <- x$theta_se
  }
  rownames(frame_tmp) <- names(x$theta)
  print(format(signif(frame_tmp, max(1L, digits - 2L))))
  cat(paste("\nlog-likelihood = "),x$logLike)
  cat(paste("\nAIC = "),x$AIC)
  cat(paste("\nBIC = "),x$BIC)
  invisible(x)
}


#' @title Print summary of an \code{ARmdl} object
#'
#' @description This is a method for the function \code{print()} for objects of the class \code{ARmdl}.
#' 
#' @inheritParams base::print
#'
#' @return The \code{ARmdl} object is returned invisibly.
#' 
#' @keywords internal
#' 
#' @export
print.ARmdl <- function(x, digits = getOption("digits"), ...){
  cat("\nAutoregressive Model\n")
  frame_tmp <- data.frame(coef = x$theta)
  if (x$control$getSE==TRUE){
    frame_tmp["s.e."] <- x$theta_se
  }
  rownames(frame_tmp) <- names(x$theta)
  print(format(signif(frame_tmp, max(1L, digits - 2L))))
  cat(paste("\nlog-likelihood = "),x$logLike)
  cat(paste("\nAIC = "),x$AIC)
  cat(paste("\nBIC = "),x$BIC)
  invisible(x)
}

#' @title Print summary of an \code{VARmdl} object
#'
#' @description This is a method for the function \code{print()} for objects of the class \code{VARmdl}.
#' 
#' @inheritParams base::print
#'
#' @return The \code{VARmdl} object is returned invisibly.
#' 
#' @keywords internal
#' 
#' @export
print.VARmdl <- function(x, digits = getOption("digits"), ...){
  cat("\nVector Autoregressive Model\n")
  frame_tmp <- data.frame(coef = x$theta)
  if (x$control$getSE==TRUE){
    frame_tmp["s.e."] <- x$theta_se
  }
  rownames(frame_tmp) <- names(x$theta)
  print(format(signif(frame_tmp, max(1L, digits - 2L))))
  cat(paste("\nlog-likelihood = "),x$logLike)
  cat(paste("\nAIC = "),x$AIC)
  cat(paste("\nBIC = "),x$BIC)
  invisible(x)
}

#' @title Print summary of a \code{HMmdl} object
#'
#' @description This is a method for the function \code{print()} for objects of the class \code{HMmdl}.
#' 
#' @inheritParams base::print
#'
#' @return The \code{HMmdl} object is returned invisibly.
#' 
#' @keywords internal
#' 
#' @export
print.HMmdl <- function(x, digits = getOption("digits"), ...){
  cat("\nHidden Markov Model\n")
  frame_tmp <- data.frame(coef = x$theta)
  if (x$control$getSE==TRUE){
    frame_tmp["s.e."] <- x$theta_se
  }
  rownames(frame_tmp) <- names(x$theta)
  print(format(signif(frame_tmp, max(1L, digits - 2L))))
  cat(paste("\nlog-likelihood = "),x$logLike)
  cat(paste("\nAIC = "),x$AIC)
  cat(paste("\nBIC = "),x$BIC)
  invisible(x)
}



#' @title Print summary of a \code{MSARmdl} object
#'
#' @description This is a method for the function \code{print()} for objects of the class \code{MSARmdl}.
#' 
#' @inheritParams base::print
#'
#' @return The \code{MSARmdl} object is returned invisibly.
#' 
#' @keywords internal
#' 
#' @export
print.MSARmdl <- function(x, digits = getOption("digits"), ...){
  cat("\nMarkov Switching Autoregressive Model\n")
  frame_tmp <- data.frame(coef = x$theta)
  if (x$control$getSE==TRUE){
    frame_tmp["s.e."] <- x$theta_se
  }
  rownames(frame_tmp) <- names(x$theta)
  print(format(signif(frame_tmp, max(1L, digits - 2L))))
  cat(paste("\nlog-likelihood = "),x$logLike)
  cat(paste("\nAIC = "),x$AIC)
  cat(paste("\nBIC = "),x$BIC)
  invisible(x)
}


#' @title Print summary of a \code{MSVARmdl} object
#'
#' @description This is a method for the function \code{print()} for objects of the class \code{MSVARmdl}.
#' 
#' @inheritParams base::print
#'
#' @return The \code{MSVARmdl} object is returned invisibly.
#' 
#' @keywords internal
#' 
#' @export
print.MSVARmdl <- function(x, digits = getOption("digits"), ...){
  cat("\nMarkov Switching Vector Autoregressive Model\n")
  frame_tmp <- data.frame(coef = x$theta)
  if (x$control$getSE==TRUE){
    frame_tmp["s.e."] <- x$theta_se
  }
  rownames(frame_tmp) <- names(x$theta)
  print(format(signif(frame_tmp, max(1L, digits - 2L))))
  cat(paste("\nlog-likelihood = "),x$logLike)
  cat(paste("\nAIC = "),x$AIC)
  cat(paste("\nBIC = "),x$BIC)
  invisible(x)
}

#' @title Print summary of a \code{CHPTest} object
#'
#' @description This is a method for the function \code{print()} for objects of the class \code{CHPTest}.
#' 
#' @inheritParams base::print
#'
#' @return The \code{CHPTest} object is returned invisibly.
#' 
#' @keywords internal
#' 
#' @export
print.HLRTest <- function(x, digits = getOption("digits"), ...){
  cat("\nRestricted Model\n")
  frame_tmp <- data.frame(coef = x$mdl_h0$theta)
  if (x$mdl_h0$control$getSE==TRUE){
    frame_tmp["s.e."] <- x$mdl_h0$theta_se
  }
  rownames(frame_tmp) <- names(x$mdl_h0$theta)
  print(format(signif(frame_tmp, max(1L, digits - 2L))))
  cat(paste("\nlog-likelihood = "),x$mdl_h0$logLike)
  cat(paste("\nAIC = "),x$mdl_h0$AIC)
  cat(paste("\nBIC = "),x$mdl_h0$BIC)
  cat("\n")
  if (x$control$msvar){
    cat("\nHansen (1992) Likelihood Ratio Bound Test -  Switch in Mean and Variance\n")
  }else{
    cat("\nHansen (1992) Likelihood Ratio Bound Test -  Switch in Mean only\n") 
  }
  out <- data.frame(cbind(x$LR0, x$LR_cv, x$pval))
  colnames(out) <- c("test-stat", colnames(x$LR_cv), "p-value")
  print(format(signif(out, max(1L, digits - 2L))))
  invisible(x)
}

#' @title Print summary of a \code{CHPTest} object
#'
#' @description This is a method for the function \code{print()} for objects of the class \code{CHPTest}.
#' 
#' @inheritParams base::print
#'
#' @return The \code{CHPTest} object is returned invisibly.
#' 
#' @keywords internal
#' 
#' @export
print.CHPTest <- function(x, digits = getOption("digits"), ...){
  cat("\nRestricted Model\n")
  frame_tmp <- data.frame(coef = x$mdl_h0$theta)
  if (x$mdl_h0$control$getSE==TRUE){
    frame_tmp["s.e."] <- x$mdl_h0$theta_se
  }
  rownames(frame_tmp) <- names(x$mdl_h0$theta)
  print(format(signif(frame_tmp, max(1L, digits - 2L))))
  cat(paste("\nlog-likelihood = "),x$mdl_h0$logLike)
  cat(paste("\nAIC = "),x$mdl_h0$AIC)
  cat(paste("\nBIC = "),x$mdl_h0$BIC)
  cat("\n")
  if (x$control$msvar){
    cat("\nCarrasco, Hu, & Ploberger (2014) Parameter Stability Test -  Switch in Mean and Variance\n")
  }else{
    cat("\nCarrasco, Hu, & Ploberger (2014) Parameter Stability Test -  Switch in Mean only\n") 
  }
  out <- data.frame(rbind(c(x$supTS, x$supTS_cv, x$pval_supTS),
                          c(x$expTS, x$expTS_cv, x$pval_expTS)))
  colnames(out) <- c("test-stat", names(x$supTS_cv), "p-value")
  rownames(out) <- c("supTS", "expTS")
  print(format(signif(out, max(1L, digits - 2L))))
  invisible(x)
}

#' @title Print summary of a \code{DLMCTest} object
#'
#' @description This is a method for the function \code{print()} for objects of the class \code{DLMCTest}.
#' 
#' @inheritParams base::print
#'
#' @return The \code{DLMCTest} object is returned invisibly.
#' 
#' @keywords internal
#' 
#' @export
print.DLMCTest <- function(x, digits = getOption("digits"), ...){
  cat("\nRestricted Model\n") 
  frame_tmp <- data.frame(coef = x$mdl_h0$theta)
  if (x$mdl_h0$control$getSE==TRUE){
    frame_tmp["s.e."] <- x$mdl_h0$theta_se
  }
  rownames(frame_tmp) <- names(x$mdl_h0$theta)
  print(format(signif(frame_tmp, max(1L, digits - 2L))))
  cat(paste("\nlog-likelihood = "),x$mdl_h0$logLike)
  cat(paste("\nAIC = "),x$mdl_h0$AIC)
  cat(paste("\nBIC = "),x$mdl_h0$BIC)
  cat("\n")
  cat("\nDufour & Luger (2017) Moment-Based Local Monte Carlo Test\n")
  out <- data.frame(rbind(c(t(x$theta),x$S0, x$F0_min, x$FN_min_cv, x$pval_min),
                          c(t(x$theta),x$S0, x$F0_prod, x$FN_prod_cv, x$pval_prod)))
  colnames(out) <- c(rownames(x$theta), colnames(x$S0), colnames(x$F0_min), names(x$FN_min_cv), "p-value")
  rownames(out) <- c("LMC_min", "LMC_prod")
  print(format(signif(out, max(1L, digits - 2L))))
  invisible(x)
}

#' @title Print summary of a \code{DLMMCTest} object
#'
#' @description This is a method for the function \code{print()} for objects of the class \code{DLMMCTest}.
#' 
#' @inheritParams base::print
#'
#' @return The \code{DLMMCTest} object is returned invisibly.
#' 
#' @keywords internal
#' 
#' @export
print.DLMMCTest <- function(x, digits = getOption("digits"), ...){
  cat("\nRestricted Model\n") 
  frame_tmp <- data.frame(coef = x$mdl_h0$theta)
  if (x$mdl_h0$control$getSE==TRUE){
    frame_tmp["s.e."] <- x$mdl_h0$theta_se
  }
  rownames(frame_tmp) <- names(x$mdl_h0$theta)
  print(format(signif(frame_tmp, max(1L, digits - 2L))))
  cat(paste("\nlog-likelihood = "),x$mdl_h0$logLike)
  cat(paste("\nAIC = "),x$mdl_h0$AIC)
  cat(paste("\nBIC = "),x$mdl_h0$BIC)
  cat("\n")
  cat("\nDufour & Luger (2017) Moment-Based Maximized Monte Carlo Test\n")
  out <- data.frame(rbind(c(t(x$theta_max_min),x$S0_min, x$F0_min, x$FN_min_cv, x$pval_min),
                          c(t(x$theta_max_prod),x$S0_prod, x$F0_prod, x$FN_prod_cv, x$pval_prod)))
  colnames(out) <- c(rownames(x$theta_max_min), colnames(x$S0_min), colnames(x$F0_min), names(x$FN_min_cv), "p-value")
  rownames(out) <- c("MMC_min","MMC_prod")
  print(format(signif(out, max(1L, digits - 2L))))
  invisible(x)
}

#' @title Print summary of a \code{LMCLRTest} object
#'
#' @description This is a method for the function \code{print()} for objects of the class \code{LMCLRTest}.
#' 
#' @inheritParams base::print
#'
#' @return The \code{LMCLRTest} object is returned invisibly.
#' 
#' @keywords internal
#' 
#' @export
print.LMCLRTest <- function(x, digits = getOption("digits"), ...){
  cat("\nRestricted Model\n") 
  frame_h0_tmp <- data.frame(coef = x$mdl_h0$theta)
  if (x$mdl_h0$control$getSE==TRUE){
    frame_h0_tmp["s.e."] <- x$mdl_h0$theta_se
  }
  rownames(frame_h0_tmp) <- names(x$mdl_h0$theta)
  print(format(signif(frame_h0_tmp, max(1L, digits - 2L))))
  cat(paste("\nlog-likelihood = "),x$mdl_h0$logLike)
  cat(paste("\nAIC = "),x$mdl_h0$AIC)
  cat(paste("\nBIC = "),x$mdl_h0$BIC)
  cat("\n")
  cat("\nUnrestricted Model\n") 
  frame_h1_tmp <- data.frame(coef = x$mdl_h1$theta)
  if (x$mdl_h1$control$getSE==TRUE){
    frame_h1_tmp["s.e."] <- x$mdl_h1$theta_se
  }
  rownames(frame_h1_tmp) <- names(x$mdl_h1$theta)
  print(format(signif(frame_h1_tmp, max(1L, digits - 2L))))
  cat(paste("\nlog-likelihood = "),x$mdl_h1$logLike)
  cat(paste("\nAIC = "),x$mdl_h1$AIC)
  cat(paste("\nBIC = "),x$mdl_h1$BIC)
  cat("\n")
  cat("\nRodriguez Rondon & Dufour (2022) Local Monte Carlo Likelihood Ratio Test\n")
  out <- data.frame(t(as.matrix(c(x$LRT_0, x$LRN_cv, x$pval))))
  colnames(out) <- c(names(x$LRT_0), names(x$LRN_cv), "p-value")
  rownames(out) <- "LMC_LRT"
  print(format(signif(out, max(1L, digits - 2L))))
  invisible(x)
}



#' @title Print summary of a \code{MMCLRTest} object
#'
#' @description This is a method for the function \code{print()} for objects of the class \code{MMCLRTest}.
#' 
#' @inheritParams base::print
#'
#' @return The \code{MMCLRTest} object is returned invisibly.
#' 
#' @keywords internal
#' 
#' @export
print.MMCLRTest <- function(x, digits = getOption("digits"), ...){
  cat("\nRestricted Model\n") 
  frame_h0_tmp <- data.frame(coef = x$mdl_h0$theta)
  if (x$mdl_h0$control$getSE==TRUE){
    frame_h0_tmp["s.e."] <- x$mdl_h0$theta_se
  }
  rownames(frame_h0_tmp) <- names(x$mdl_h0$theta)
  print(format(signif(frame_h0_tmp, max(1L, digits - 2L))))
  cat(paste("\nlog-likelihood = "),x$mdl_h0$logLike)
  cat(paste("\nAIC = "),x$mdl_h0$AIC)
  cat(paste("\nBIC = "),x$mdl_h0$BIC)
  cat("\n")
  cat("\nUnrestricted Model\n") 
  frame_h1_tmp <- data.frame(coef = x$mdl_h1$theta)
  if (x$mdl_h1$control$getSE==TRUE){
    frame_h1_tmp["s.e."] <- x$mdl_h1$theta_se
  }
  rownames(frame_h1_tmp) <- names(x$mdl_h1$theta)
  print(format(signif(frame_h1_tmp, max(1L, digits - 2L))))
  cat(paste("\nlog-likelihood = "),x$mdl_h1$logLike)
  cat(paste("\nAIC = "),x$mdl_h1$AIC)
  cat(paste("\nBIC = "),x$mdl_h1$BIC)
  cat("\n")
  cat("\nRodriguez Rondon & Dufour (2022) Maximized Monte Carlo Likelihood Ratio Test\n")
  out <- data.frame(t(as.matrix(c(x$LRT_0, x$pval))))
  colnames(out) <- c(names(x$LRT_0), "p-value")
  rownames(out) <- "MMC_LRT"
  print(format(signif(out, max(1L, digits - 2L))))
  invisible(x)
}


#' @title Print summary of a \code{BootLRTest} object
#'
#' @description This is a method for the function \code{print()} for objects of the class \code{BootLRTest}.
#' 
#' @inheritParams base::print
#'
#' @return The \code{BootLRTest} object is returned invisibly.
#' 
#' @keywords internal
#' 
#' @export
print.BootLRTest <- function(x, digits = getOption("digits"), ...){
  cat("\nRestricted Model\n") 
  frame_h0_tmp <- data.frame(coef = x$mdl_h0$theta)
  if (x$mdl_h0$control$getSE==TRUE){
    frame_h0_tmp["s.e."] <- x$mdl_h0$theta_se
  }
  rownames(frame_h0_tmp) <- names(x$mdl_h0$theta)
  print(format(signif(frame_h0_tmp, max(1L, digits - 2L))))
  cat(paste("\nlog-likelihood = "),x$mdl_h0$logLike)
  cat(paste("\nAIC = "),x$mdl_h0$AIC)
  cat(paste("\nBIC = "),x$mdl_h0$BIC)
  cat("\n")
  cat("\nUnrestricted Model\n") 
  frame_h1_tmp <- data.frame(coef = x$mdl_h1$theta)
  if (x$mdl_h1$control$getSE==TRUE){
    frame_h1_tmp["s.e."] <- x$mdl_h1$theta_se
  }
  rownames(frame_h1_tmp) <- names(x$mdl_h1$theta)
  print(format(signif(frame_h1_tmp, max(1L, digits - 2L))))
  cat(paste("\nlog-likelihood = "),x$mdl_h1$logLike)
  cat(paste("\nAIC = "),x$mdl_h1$AIC)
  cat(paste("\nBIC = "),x$mdl_h1$BIC)
  cat("\n")
  cat("\nBootstrap Likelihood Ratio Test\n")
  out <- data.frame(t(as.matrix(c(x$LRT_0, x$LRN_cv, x$pval))))
  colnames(out) <- c(names(x$LRT_0), names(x$LRN_cv), "p-value")
  rownames(out) <- "Boot_LRT"
  print(format(signif(out, max(1L, digits - 2L))))
  invisible(x)
}
