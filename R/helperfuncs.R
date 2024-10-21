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
    info_mat <- pracma::nearest_spd(info_mat)
    nearPD_used <- TRUE
  }
  mdl$Hess <- Hess
  mdl$theta_se <- sqrt(diag(info_mat))
  mdl$info_mat <- info_mat
  mdl$nearPD_used <- nearPD_used 
  return(mdl)
}

#' @title Intercept from mu for MSARmdl
#' 
#' @description This function computes the intercept for each regime k for an Markov switching AR model
#'
#' @param mdl List with model properties
#'
#' @return a \code{(k x 1)} vector of intercepts
#' 
#' @keywords internal
#' 
#' @export
interMSARmdl <- function(mdl){
  inter <- matrix(0,mdl$k,1)
  for(xk in 1:mdl$k){
    inter[xk,] <- t(as.matrix(mdl$mu[xk,]))*(1-sum(mdl$phi))
  }
  return(inter)
}

#' @title Intercept from mu for MSVARmdl
#' 
#' @description This function computes the intercept for each regime k for an Markov switching VAR model
#'
#' @param mdl List with model properties
#'
#' @return a \code{(k x q)} vector of intercepts
#' 
#' @keywords internal
#' 
#' @export
interMSVARmdl <- function(mdl){
  inter <- matrix(0,mdl$k,mdl$q)
  for(xk in 1:mdl$k){
    nu_tmp <- (diag(mdl$q*mdl$p)-mdl$Fmat)%*%as.matrix(c(mdl$mu[xk,],rep(0,mdl$q*(mdl$p-1))))
    inter[xk,] <- nu_tmp[(1:(mdl$q))]
  }
  return(inter)
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
#'  \item maxit: maximum number of iterations.
#'  \item thtol: convergence criterion.
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
                       lower = optim_options$mle_theta_low,
                       upper = optim_options$mle_theta_upp,
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
#'  \item maxit: maximum number of iterations.
#'  \item thtol: convergence criterion.
#'}
#'
#' @return List with model attributes
#' 
#' @keywords internal
#' 
#' @export
MSARmdl_mle <- function(theta_0, mdl_in, k, optim_options){
  # ----- Define function environment
  MSARmdl_mle_env         <- rlang::env()
  # set environment variables
  MSARmdl_mle_env$p       <- mdl_in$p
  MSARmdl_mle_env$k       <- k
  MSARmdl_mle_env$msmu    <- mdl_in$msmu
  MSARmdl_mle_env$msvar   <- mdl_in$msvar
  MSARmdl_mle_env$var_k1  <- mdl_in$sigma
  MSARmdl_mle_env$betaZ   <- mdl_in$betaZ
  MSARmdl_mle_env$mle_stationary_constraint <- mdl_in$mle_stationary_constraint
  MSARmdl_mle_env$mle_variance_constraint   <- mdl_in$mle_variance_constraint
  # ----------- MSARmdl_mle equality constraint functions
  loglik_const_eq_MSARmdl <- function(theta){
    # ----- Load equality constraint parameters
    k <- get("k", envir = MSARmdl_mle_env)
    # ----- constraint
    P <- matrix(theta[(length(theta)-k*k+1):(length(theta))],k,k)
    constraint <- colSums(P)-1
    return(constraint)
  }
  # ---------- MSARmdl_mle inequality constraint functions
  loglik_const_ineq_MSARmdl <- function(theta){
    # ----- Load inequality constraint parameters
    p       <- get("p", envir = MSARmdl_mle_env)
    k       <- get("k", envir = MSARmdl_mle_env)
    msmu    <- get("msmu", envir = MSARmdl_mle_env)
    msvar   <- get("msvar", envir = MSARmdl_mle_env)
    var_k1  <- get("var_k1", envir = MSARmdl_mle_env)
    betaZ   <- get("betaZ", envir = MSARmdl_mle_env)
    mle_stationary_constraint <- get("mle_stationary_constraint", envir = MSARmdl_mle_env)
    mle_variance_constraint   <- get("mle_variance_constraint", envir = MSARmdl_mle_env)
    # ----- constraint
    # transition probabilities
    Pvec <- theta[(length(theta)-k*k+1):(length(theta))]
    ineq_constraint <- c(Pvec, 1-Pvec)
    if (mle_stationary_constraint==TRUE){
      # roots of characteristic function
      n_p <- 2 + msmu*(k-1)
      poly_fun <- c(1,-theta[n_p:(n_p+p-1)])
      roots <- Mod(polyroot(poly_fun))
      ineq_constraint = c((roots-1), ineq_constraint)
    }
    if (mle_variance_constraint>=0){
      # variances (lower bound is 1% of single regime variance)
      vars <- theta[c(rep(0, 1+(k-1)*msmu+p+length(betaZ)), rep(1, 1 + (k-1)*msvar), rep(0, k*k))==1] - c(mle_variance_constraint*var_k1)
      ineq_constraint <- c(vars, ineq_constraint)
    }
    return(ineq_constraint)
  }
  # use nloptr optimization to minimize (maximize) likelihood
  if (length(mdl_in$betaZ)>0){
    res <- nloptr::slsqp(x0 = theta_0,
                         fn = logLike_MSARXmdl_min,
                         gr = NULL,
                         lower = optim_options$mle_theta_low,
                         upper = optim_options$mle_theta_upp,
                         hin = loglik_const_ineq_MSARmdl,
                         heq = loglik_const_eq_MSARmdl,
                         nl.info = FALSE,
                         control = list(maxeval = optim_options$maxit, xtol_rel = optim_options$thtol),
                         mdl = mdl_in,
                         k = k) 
    output <- ExpectationM_MSARXmdl(res$par, mdl_in, k)
  }else{
    res <- nloptr::slsqp(x0 = theta_0,
                         fn = logLike_MSARmdl_min,
                         gr = NULL,
                         lower = optim_options$mle_theta_low,
                         upper = optim_options$mle_theta_upp,
                         hin = loglik_const_ineq_MSARmdl,
                         heq = loglik_const_eq_MSARmdl,
                         nl.info = FALSE,
                         control = list(maxeval = optim_options$maxit, xtol_rel = optim_options$thtol),
                         mdl = mdl_in,
                         k = k) 
    output <- ExpectationM_MSARmdl(res$par, mdl_in, k)
  }
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
#'  \item maxit: maximum number of iterations.
#'  \item thtol: convergence criterion.
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
                       lower = optim_options$mle_theta_low,
                       upper = optim_options$mle_theta_upp,
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

