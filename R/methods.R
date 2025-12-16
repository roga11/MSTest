#' @title coef of a \code{Nmdl} object
#'
#' @description This is a method for the function \code{coef()} for objects of the class \code{Nmdl}.
#' 
#' @inheritParams stats::coef
#'
#' @return vector of coefficients. 
#' 
#' @keywords internal
#' 
#' @export
coef.Nmdl <- function(object, ...){
  return(object$theta)
}

#' @title coef of a \code{ARmdl} object
#'
#' @description This is a method for the function \code{coef()} for objects of the class \code{ARmdl}.
#' 
#' @inheritParams stats::coef
#'
#' @return vector of coefficients. 
#' 
#' @keywords internal
#' 
#' @export
coef.ARmdl <- function(object, ...){
  return(object$theta)
}

#' @title coef of a \code{VARmdl} object
#'
#' @description This is a method for the function \code{coef()} for objects of the class \code{VARmdl}.
#' 
#' @inheritParams stats::coef
#'
#' @return vector of coefficients. 
#' 
#' @keywords internal
#' 
#' @export
coef.VARmdl <- function(object, ...){
  return(object$theta)
}

#' @title coef of a \code{HMmdl} object
#'
#' @description This is a method for the function \code{coef()} for objects of the class \code{HMmdl}.
#' 
#' @inheritParams stats::coef
#'
#' @return vector of coefficients. 
#' 
#' @keywords internal
#' 
#' @export
coef.HMmdl <- function(object, ...){
  return(object$theta)
}

#' @title coef of a \code{MSARmdl} object
#'
#' @description This is a method for the function \code{coef()} for objects of the class \code{MSARmdl}.
#' 
#' @inheritParams stats::coef
#'
#' @return vector of coefficients. 
#' 
#' @keywords internal
#' 
#' @export
coef.MSARmdl <- function(object, ...){
  return(object$theta)
}

#' @title coef of a \code{MSVARmdl} object
#'
#' @description This is a method for the function \code{coef()} for objects of the class \code{MSVARmdl}.
#' 
#' @inheritParams stats::coef
#'
#' @return vector of coefficients. 
#' 
#' @keywords internal
#' 
#' @export
coef.MSVARmdl <- function(object, ...){
  return(object$theta)
}



#' @title fitted values of a \code{Nmdl} object
#' 
#' @description This is a method for the function \code{fitted()} for objects of the class \code{Nmdl}.
#' 
#' @inheritParams stats::fitted
#'
#' @return matrix with fitted values
#' 
#' @keywords internal
#' 
#' @export
fitted.Nmdl <- function(object, ...){
  return(object$fitted)
}

#' @title fitted values of a \code{ARmdl} object
#' 
#' @description This is a method for the function \code{fitted()} for objects of the class \code{ARmdl}.
#' 
#' @inheritParams stats::fitted
#'
#' @return matrix with fitted values
#' 
#' @keywords internal
#' 
#' @export
fitted.ARmdl <- function(object, ...){
  return(object$fitted)
}

#' @title fitted values of a \code{VARmdl} object
#' 
#' @description This is a method for the function \code{fitted()} for objects of the class \code{VARmdl}.
#' 
#' @inheritParams stats::fitted
#'
#' @return matrix with fitted values
#' 
#' @keywords internal
#' 
#' @export
fitted.VARmdl <- function(object, ...){
  return(object$fitted)
}


#' @title fitted values of a \code{HMmdl} object
#' 
#' @description This is a method for the function \code{fitted()} for objects of the class \code{HMmdl}.
#' 
#' @inheritParams stats::fitted
#'
#' @return matrix with fitted values
#' 
#' @keywords internal
#' 
#' @export
fitted.HMmdl <- function(object, ...){
  fittedvals <- matrix(0,object$n,object$q)
  for( xk in 1:object$k){
    fittedvals <- fittedvals + object$St[,xk] * matrix(1, object$n, 1)%*%t(as.matrix(object$mu[xk,]))  
  }
  return(fittedvals)
}


#' @title fitted values of a \code{MSARmdl} object
#' 
#' @description This is a method for the function \code{fitted()} for objects of the class \code{MSARmdl}.
#' 
#' @inheritParams stats::fitted
#'
#' @return matrix with fitted values
#' 
#' @keywords internal
#' 
#' @export
fitted.MSARmdl <- function(object, ...){
  fittedvals <- matrix(0,object$n,object$q)
  if (is.null(object$betaZ)){
    for(xk in 1:object$k){
      intercept <- t(as.matrix(object$mu[xk,]))*(1-sum(object$phi))
      fittedvals <- fittedvals + 
        object$St[,xk] * (matrix(1, object$n, 1)%*%intercept + 
                            object$x%*%object$phi)
    }  
  }else{
    for(xk in 1:object$k){
      intercept <- object$intercept[xk,]
      fittedvals <- fittedvals + 
        object$St[,xk] * (matrix(1, object$n, 1)%*%intercept + 
                            object$x%*%object$phi + object$Z%*%object$betaZ)
    }
  }
  return(fittedvals)
}



#' @title fitted values of a \code{MSVARmdl} object
#' 
#' @description This is a method for the function \code{fitted()} for objects of the class \code{MSVARmdl}.
#' 
#' @inheritParams stats::fitted
#'
#' @return matrix with fitted values
#' 
#' @keywords internal
#' 
#' @export
fitted.MSVARmdl <- function(object, ...){
  fittedvals <- matrix(0,object$n,object$q)
  for( xk in 1:object$k){
    fittedvals <- fittedvals + 
      object$St[,xk] * (matrix(1, object$n, 1)%*%t(as.matrix(object$inter[xk,])) + 
      object$x%*%t(object$phi))
  }
  return(fittedvals)
}

#' @title residuals of a \code{Nmdl} object
#'
#' @description This is a method for the function \code{residuals()} for objects of the class \code{Nmdl}.
#' 
#' @inheritParams stats::residuals
#'
#' @return vector of residuals. 
#' 
#' @keywords internal
#' 
#' @export
residuals.Nmdl <- function(object, ...){
  return(object$resid)
}

#' @title residuals of a \code{ARmdl} object
#'
#' @description This is a method for the function \code{residuals()} for objects of the class \code{ARmdl}.
#' 
#' @inheritParams stats::residuals
#'
#' @return vector of residuals. 
#' 
#' @keywords internal
#' 
#' @export
residuals.ARmdl <- function(object, ...){
  return(object$resid)
}


#' @title residuals of a \code{VARmdl} object
#'
#' @description This is a method for the function \code{residuals()} for objects of the class \code{VARmdl}.
#' 
#' @inheritParams stats::residuals
#'
#' @return vector of residuals. 
#' 
#' @keywords internal
#' 
#' @export
residuals.VARmdl <- function(object, ...){
  return(object$resid)
}


#' @title residuals of a \code{HMmdl} object
#'
#' @description This is a method for the function \code{residuals()} for objects of the class \code{HMmdl}.
#' 
#' @inheritParams stats::residuals
#'
#' @return vector of residuals. 
#' 
#' @keywords internal
#' 
#' @export
residuals.HMmdl <- function(object, ...){
  return(object$resid)
}

#' @title residuals of a \code{MSARmdl} object
#'
#' @description This is a method for the function \code{residuals()} for objects of the class \code{MSARmdl}.
#' 
#' @inheritParams stats::residuals
#'
#' @return vector of residuals. 
#' 
#' @keywords internal
#' 
#' @export
residuals.MSARmdl <- function(object, ...){
  return(object$resid)
}


#' @title residuals of a \code{MSVARmdl} object
#'
#' @description This is a method for the function \code{residuals()} for objects of the class \code{MSVARmdl}.
#' 
#' @inheritParams stats::residuals
#'
#' @return vector of residuals. 
#' 
#' @keywords internal
#' 
#' @export
residuals.MSVARmdl <- function(object, ...){
  return(object$resid)
}


#' @title Nobs of a \code{Nmdl} object
#'
#' @description This is a method for the function \code{nobs()} for objects of the class \code{Nmdl}.
#' 
#' @import stats
#' 
#' @inheritParams stats::nobs
#' 
#' @return Number of time series observations. 
#' 
#' @keywords internal
#' 
#' @export 
nobs.Nmdl <- function(object, ...){
  return(object$n)
}

#' @title Nobs of a \code{ARmdl} object
#'
#' @description This is a method for the function \code{nobs()} for objects of the class \code{ARmdl}.
#' 
#' @import stats
#' 
#' @inheritParams stats::nobs
#'
#' @return Number of time series observations. 
#' 
#' @keywords internal
#' 
#' @export 
nobs.ARmdl <- function(object, ...){
  return(object$n)
}

#' @title Nobs of a \code{VARmdl} object
#'
#' @description This is a method for the function \code{nobs()} for objects of the class \code{VARmdl}.
#' 
#' @import stats
#' 
#' @inheritParams stats::nobs
#'
#' @return Number of time series observations. 
#' 
#' @keywords internal
#' 
#' @export 
nobs.VARmdl <- function(object, ...){
  return(object$n)
}

#' @title Nobs of a \code{HMmdl} object
#'
#' @description This is a method for the function \code{nobs()} for objects of the class \code{HMmdl}.
#' 
#' @import stats
#' 
#' @inheritParams stats::nobs
#'
#' @return Number of time series observations. 
#' 
#' @keywords internal
#' 
#' @export 
nobs.HMmdl <- function(object, ...){
  return(object$n)
}


#' @title Nobs of a \code{MSARmdl} object
#'
#' @description This is a method for the function \code{nobs()} for objects of the class \code{MSARmdl}.
#'
#' @import stats
#' 
#' @inheritParams stats::nobs
#'
#' @return Number of time series observations. 
#' 
#' @keywords internal
#' 
#' @export 
nobs.MSARmdl <- function(object, ...){
  return(object$n)
}

#' @title Nobs of a \code{MSVARmdl} object
#'
#' @description This is a method for the function \code{nobs()} for objects of the class \code{MSVARmdl}.
#' 
#' @import stats
#' 
#' @inheritParams stats::nobs
#'
#' @return Number of time series observations. 
#' 
#' @keywords internal
#' 
#' @export
nobs.MSVARmdl <- function(object, ...){
  return(object$n)
}



#' @title Log likelihood for Normal model  
#' 
#' @description This function is used to compute the log-likelihood for a normally distributed model.
#'
#' @inheritParams stats::logLik
#'
#' @return Log-likelihood value.
#' 
#' @keywords internal
#' 
#' @export
logLik.Nmdl <- function(object, ...){
  logLike <- logLike_Nmdl(stats::coef(object), object)
  return(logLike)
}

#' @title Log likelihood for autoregressive model  
#' 
#' @description This function is used to compute the log-likelihood for an autoregressive model.
#'
#' @inheritParams stats::logLik
#'
#' @return Log-likelihood value.
#' 
#' @keywords internal
#' 
#' @export
logLik.ARmdl <- function(object, ...){
  if (is.null(object$control$Z)){
    logLike <- logLike_ARmdl(stats::coef(object), object)  
  }else{
    logLike <- logLike_ARXmdl(stats::coef(object), object)  
  }
  return(logLike)
}

#' @title Log likelihood for vector autoregressive model  
#' 
#' @description This function is used to compute the log-likelihood for a vector autoregressive model.
#'
#' @inheritParams stats::logLik
#'
#' @return Log-likelihood value.
#' 
#' @keywords internal
#' 
#' @export
logLik.VARmdl <- function(object, ...){
  if (is.null(object$control$Z)){
    logLike <- logLike_VARmdl(stats::coef(object), object)
  }else{
    logLike <- logLike_VARXmdl(stats::coef(object), object)
  }
  return(logLike)
}

#' @title Log likelihood for Hidden Markov model  
#' 
#' @description This function is used to compute the log-likelihood for a Hidden Markov model.
#'
#' @inheritParams stats::logLik
#'
#' @return Log-likelihood value. 
#' 
#' @keywords internal
#' 
#' @export
logLik.HMmdl <- function(object, ...){
  logLike <- logLike_HMmdl(stats::coef(object), object, object$k)
  return(logLike)
}

#' @title Log likelihood for Markov-switching autoregressive model  
#' 
#' @description This function is used to compute the log-likelihood for a Markov-switching autoregressive model.
#'
#' @inheritParams stats::logLik
#'
#' @return Log-likelihood value.
#' 
#' @keywords internal
#' 
#' @export
logLik.MSARmdl <- function(object, ...){
  if (is.null(object$control$Z)){
    logLike <- logLike_MSARmdl(stats::coef(object), object, object$k)
  }else{
    logLike <- logLike_MSARXmdl(stats::coef(object), object, object$k)
  }
  return(logLike)
}

#' @title Log likelihood for Markov-switching vector autoregressive model  
#' 
#' @description This function is used to compute the log-likelihood for a Markov-switching vector autoregressive model.
#'
#' @inheritParams stats::logLik
#'
#' @return Log-likelihood value.
#' 
#' @keywords internal
#' 
#' @export
logLik.MSVARmdl <- function(object, ...){
  if (is.null(object$control$Z)){
    logLike <- logLike_MSVARmdl(stats::coef(object), object, object$k)
  }else{
    logLike <- logLike_MSVARXmdl(stats::coef(object), object, object$k)  
  }
  return(logLike)
}



#' @title AIC of a \code{Nmdl} object
#'
#' @description This is a method for the function \code{AIC()} for objects of the class \code{Nmdl}.
#' 
#' @inheritParams stats::AIC
#'
#' @return AIC value. 
#' 
#' @keywords internal
#' 
#' @export
AIC.Nmdl <- function(object, ..., k = 2){
  aic_val <-  k*length(stats::coef(object)) - 2*(stats::logLik(object))
  return(aic_val)
}

#' @title AIC of a \code{ARmdl} object
#'
#' @description This is a method for the function \code{AIC()} for objects of the class \code{ARmdl}.
#' 
#' @inheritParams stats::AIC
#'
#' @return AIC value. 
#' 
#' @keywords internal
#' 
#' @export
AIC.ARmdl <- function(object, ..., k = 2){
  aic_val <-  k*length(stats::coef(object)) - 2*(stats::logLik(object))
  return(aic_val)
}

#' @title AIC of a \code{VARmdl} object
#'
#' @description This is a method for the function \code{AIC()} for objects of the class \code{VARmdl}.
#' 
#' @inheritParams stats::AIC
#'
#' @return AIC value. 
#' 
#' @keywords internal
#' 
#' @export
AIC.VARmdl <- function(object, ..., k = 2){
  aic_val <-  k*length(stats::coef(object)) - 2*(stats::logLik(object))
  return(aic_val)
}

#' @title AIC of a \code{HMmdl} object
#'
#' @description This is a method for the function \code{AIC()} for objects of the class \code{HMmdl}.
#' 
#' @inheritParams stats::AIC
#'
#' @return AIC value. 
#' 
#' @keywords internal
#' 
#' @export
AIC.HMmdl <- function(object, ..., k = 2){
  aic_val <-  k*length(stats::coef(object)) - 2*(stats::logLik(object))
  return(aic_val)
}


#' @title AIC of a \code{MSARmdl} object
#'
#' @description This is a method for the function \code{AIC()} for objects of the class \code{MSARmdl}.
#' 
#' @inheritParams stats::AIC
#'
#' @return AIC value. 
#' 
#' @keywords internal
#' 
#' @export
AIC.MSARmdl <- function(object, ..., k = 2){
  aic_val <-  k*length(stats::coef(object)) - 2*(stats::logLik(object))
  return(aic_val)
}


#' @title AIC of a \code{MSVARmdl} object
#'
#' @description This is a method for the function \code{AIC()} for objects of the class \code{MSVARmdl}.
#' 
#' @inheritParams stats::AIC
#'
#' @return AIC value. 
#' 
#' @keywords internal
#' 
#' @export
AIC.MSVARmdl <- function(object, ..., k = 2){
  aic_val <-  k*length(stats::coef(object)) - 2*(stats::logLik(object))
  return(aic_val)
}



#' @title BIC of a \code{Nmdl} object
#'
#' @description This is a method for the function \code{BIC()} for objects of the class \code{Nmdl}.
#' 
#' @inheritParams stats::BIC
#'
#' @return BIC value. 
#'
#' @keywords internal
#' 
#' @export
BIC.Nmdl <- function(object, ...){
  bic_val <- log(stats::nobs(object))*length(stats::coef(object)) -2*(stats::logLik(object)) 
  return(bic_val)
}

#' @title BIC of a \code{ARmdl} object
#'
#' @description This is a method for the function \code{BIC()} for objects of the class \code{ARmdl}.
#' 
#' @inheritParams stats::BIC
#'
#' @return BIC value. 
#'
#' @keywords internal
#' 
#' @export
BIC.ARmdl <- function(object, ...){
  bic_val <- log(stats::nobs(object))*length(stats::coef(object)) -2*(stats::logLik(object)) 
  return(bic_val)
}

#' @title BIC of a \code{VARmdl} object
#'
#' @description This is a method for the function \code{BIC()} for objects of the class \code{VARmdl}.
#' 
#' @inheritParams stats::BIC
#'
#' @return BIC value. 
#'
#' @keywords internal
#' 
#' @export
BIC.VARmdl <- function(object, ...){
  bic_val <- log(stats::nobs(object))*length(stats::coef(object)) -2*(stats::logLik(object)) 
  return(bic_val)
}


#' @title BIC of a \code{HMmdl} object
#'
#' @description This is a method for the function \code{BIC()} for objects of the class \code{HMmdl}.
#' 
#' @inheritParams stats::BIC
#'
#' @return BIC value. 
#'
#' @keywords internal
#' 
#' @export
BIC.HMmdl <- function(object, ...){
  bic_val <- log(stats::nobs(object))*length(stats::coef(object)) -2*(stats::logLik(object)) 
  return(bic_val)
}

#' @title BIC of a \code{MSARmdl} object
#'
#' @description This is a method for the function \code{BIC()} for objects of the class \code{MSARmdl}.
#' 
#' @inheritParams stats::BIC
#'
#' @return BIC value. 
#'
#' @keywords internal
#' 
#' @export
BIC.MSARmdl <- function(object, ...){
  bic_val <- log(stats::nobs(object))*length(stats::coef(object)) -2*(stats::logLik(object)) 
  return(bic_val)
}


#' @title BIC of a \code{MSVARmdl} object
#'
#' @description This is a method for the function \code{BIC()} for objects of the class \code{MSVARmdl}.
#' 
#' @inheritParams stats::BIC
#'
#' @return BIC value. 
#'
#' @keywords internal
#' 
#' @export
BIC.MSVARmdl <- function(object, ...){
  bic_val <- log(stats::nobs(object))*length(stats::coef(object)) -2*(stats::logLik(object)) 
  return(bic_val)
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
  if (is.null(mdl$control$Z)){
    hess <- numDeriv::hessian(logLike_ARmdl, mdl$theta, method = "Richardson", mdl = mdl)   
  }else{
    hess <- numDeriv::hessian(logLike_ARXmdl, mdl$theta, method = "Richardson", mdl = mdl)   
  }
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
  if (is.null(mdl$control$Z)){
    hess <- numDeriv::hessian(logLike_VARmdl, mdl$theta, method = "Richardson", mdl = mdl) 
  }else{
    hess <- numDeriv::hessian(logLike_VARXmdl, mdl$theta, method = "Richardson", mdl = mdl) 
  }
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
  if (is.null(mdl$control$Z)){
    hess <- numDeriv::hessian(logLike_MSARmdl, mdl$theta, method = "Richardson", mdl = mdl, k = mdl$k) 
  }else{
    hess <- numDeriv::hessian(logLike_MSARXmdl, mdl$theta, method = "Richardson", mdl = mdl, k = mdl$k) 
  }
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
  if (is.null(mdl$control$Z)){
    hess <- numDeriv::hessian(logLike_MSVARmdl, mdl$theta, method = "Richardson", mdl = mdl, k = mdl$k) 
  }else{
    hess <- numDeriv::hessian(logLike_MSVARXmdl, mdl$theta, method = "Richardson", mdl = mdl, k = mdl$k) 
  }
  return(hess)
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
  cat("\nCarrasco, Hu, & Ploberger (2014) Parameter Stability Test \n")
  if (x$control$msvar){
    cat("\n- Switch in Mean and Variance\n")
  }else{
    cat("\n- Switch in Mean only\n") 
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
  cat("\nDufour & Luger (2017) Moment-Based Maximized Monte Carlo Test\n")
  out <- data.frame(rbind(c(x$S0_min, x$F0_min, x$pval_min),
                          c(x$S0_prod, x$F0_prod, x$pval_prod)))
  colnames(out) <- c(colnames(x$S0_min), colnames(x$F0_min), "p-value")
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
  cat("\nRodriguez-Rondon & Dufour (2025) Local Monte Carlo Likelihood Ratio Test\n")
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
  cat("\nRodriguez-Rondon & Dufour (2025) Maximized Monte Carlo Likelihood Ratio Test\n")
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
  cat("\nBootstrap Likelihood Ratio Test\n")
  out <- data.frame(t(as.matrix(c(x$LRT_0, x$LRN_cv, x$pval))))
  colnames(out) <- c(names(x$LRT_0), names(x$LRN_cv), "p-value")
  rownames(out) <- "Boot_LRT"
  print(format(signif(out, max(1L, digits - 2L))))
  invisible(x)
}



#' @title Summary of a \code{Nmdl} object
#'
#' @description This is a method for the function \code{summary()} for objects of the class \code{Nmdl}.
#' 
#' @inheritParams base::summary
#'
#' @return The \code{Nmdl} object is returned invisibly.
#' 
#' @keywords internal
#' 
#' @export
summary.Nmdl <- function(object, digits = getOption("digits"), ...){
  cat("\nNormally Distributed Model\n")
  frame_tmp <- data.frame(coef = object$theta)
  if (object$control$getSE==TRUE){
    frame_tmp["s.e."] <- object$theta_se
  }
  rownames(frame_tmp) <- names(object$theta)
  print(format(signif(frame_tmp, max(1L, digits - 2L))))
  cat(paste("\nlog-likelihood = "),object$logLike)
  cat(paste("\nAIC = "),object$AIC)
  cat(paste("\nBIC = "),object$BIC)
  
  resdist <- apply(object$resid,2,quantile)
  resdist <- data.frame(t(resdist))
  colnames(resdist) <- c("Min", "1Q", "Median", "3Q", "Max")
  rownames(resdist) <- paste0("Y",1:object$q)
  cat(paste("\n"))
  cat(paste("\nResiduals:\n"))
  print(format(signif(resdist, max(1L, digits - 2L))))
  
  invisible(object)
}


#' @title Summary of an \code{ARmdl} object
#'
#' @description This is a method for the function \code{summary()} for objects of the class \code{ARmdl}.
#' 
#' @inheritParams base::summary
#'
#' @return The \code{ARmdl} object is returned invisibly.
#' 
#' @keywords internal
#' 
#' @export
summary.ARmdl <- function(object, digits = getOption("digits"), ...){
  cat("\nAutoregressive Model\n")
  frame_tmp <- data.frame(coef = object$theta)
  if (object$control$getSE==TRUE){
    frame_tmp["s.e."] <- object$theta_se
  }
  rownames(frame_tmp) <- names(object$theta)
  print(format(signif(frame_tmp, max(1L, digits - 2L))))
  cat(paste("\nlog-likelihood = "),object$logLike)
  cat(paste("\nAIC = "),object$AIC)
  cat(paste("\nBIC = "),object$BIC)
  
  resdist <- quantile(object$resid)
  resdist <- data.frame(t(resdist))
  colnames(resdist) <- c("Min", "1Q", "Median", "3Q", "Max")
  rownames(resdist) <- paste0("Y",1:object$q)
  cat(paste("\n"))
  cat(paste("\nResiduals:\n"))
  print(format(signif(resdist, max(1L, digits - 2L))))
  
  invisible(object)
}


#' @title Summary of an \code{VARmdl} object
#'
#' @description This is a method for the function \code{summary()} for objects of the class \code{VARmdl}.
#' 
#' @inheritParams base::summary
#'
#' @return The \code{VARmdl} object is returned invisibly.
#' 
#' @keywords internal
#' 
#' @export
summary.VARmdl <- function(object, digits = getOption("digits"), ...){
  cat("\nVector Autoregressive Model\n")
  frame_tmp <- data.frame(coef = object$theta)
  if (object$control$getSE==TRUE){
    frame_tmp["s.e."] <- object$theta_se
  }
  rownames(frame_tmp) <- names(object$theta)
  print(format(signif(frame_tmp, max(1L, digits - 2L))))
  cat(paste("\nlog-likelihood = "),object$logLike)
  cat(paste("\nAIC = "),object$AIC)
  cat(paste("\nBIC = "),object$BIC)
  
  resdist <- apply(object$resid,2,quantile)
  resdist <- data.frame(t(resdist))
  colnames(resdist) <- c("Min", "1Q", "Median", "3Q", "Max")
  rownames(resdist) <- paste0("Y",1:object$q)
  cat(paste("\n"))
  cat(paste("\nResiduals:\n"))
  print(format(signif(resdist, max(1L, digits - 2L))))
  
  invisible(object)
}


#' @title Summary of a \code{HMmdl} object
#' 
#' @description This is a method for the function \code{summary()} for objects of the class \code{HMmdl}.
#' 
#' @inheritParams base::summary
#'
#' @return The \code{HMmdl} object is returned invisibly.
#' 
#' @keywords internal
#' 
#' @export
summary.HMmdl <- function(object, digits = getOption("digits"), ...){
  cat("\nHidden Markov Model\n")
  frame_tmp <- data.frame(coef = object$theta)
  if (object$control$getSE==TRUE){
    frame_tmp["s.e."] <- object$theta_se
  }
  rownames(frame_tmp) <- names(object$theta)
  print(format(signif(frame_tmp, max(1L, digits - 2L))))
  cat(paste("\nlog-likelihood = "),object$logLike)
  cat(paste("\nAIC = "),object$AIC)
  cat(paste("\nBIC = "),object$BIC)
  
  resdist <- apply(object$resid,2,quantile)
  resdist <- data.frame(t(resdist))
  colnames(resdist) <- c("Min", "1Q", "Median", "3Q", "Max")
  rownames(resdist) <- paste0("Y",1:object$q)
  cat(paste("\n"))
  cat(paste("\nResiduals:\n"))
  print(format(signif(resdist, max(1L, digits - 2L))))
  
  invisible(object)
}



#' @title Summary of a \code{MSARmdl} object
#'
#' @description This is a method for the function \code{summary()} for objects of the class \code{MSARmdl}.
#' 
#' @inheritParams base::summary
#'
#' @return The \code{MSARmdl} object is returned invisibly.
#' 
#' @keywords internal
#' 
#' @export
summary.MSARmdl <- function(object, digits = getOption("digits"), ...){
  cat("\nMarkov Switching Autoregressive Model\n")
  frame_tmp <- data.frame(coef = object$theta)
  if (object$control$getSE==TRUE){
    frame_tmp["s.e."] <- object$theta_se
  }
  rownames(frame_tmp) <- names(object$theta)
  print(format(signif(frame_tmp, max(1L, digits - 2L))))
  cat(paste("\nlog-likelihood = "),object$logLike)
  cat(paste("\nAIC = "),object$AIC)
  cat(paste("\nBIC = "),object$BIC)
  
  resdist <- quantile(object$resid)
  resdist <- data.frame(t(resdist))
  colnames(resdist) <- c("Min", "1Q", "Median", "3Q", "Max")
  rownames(resdist) <- paste0("Y",1:object$q)
  cat(paste("\n"))
  cat(paste("\nResiduals:\n"))
  print(format(signif(resdist, max(1L, digits - 2L))))
  
  invisible(object)
}


#' @title Summary of a \code{MSVARmdl} object
#'
#' @description This is a method for the function \code{summary()} for objects of the class \code{MSVARmdl}.
#' 
#' @inheritParams base::summary
#'
#' @return The \code{MSVARmdl} object is returned invisibly.
#' 
#' @keywords internal
#' 
#' @export
summary.MSVARmdl <- function(object, digits = getOption("digits"), ...){
  cat("\nMarkov Switching Vector Autoregressive Model\n")
  frame_tmp <- data.frame(coef = object$theta)
  if (object$control$getSE==TRUE){
    frame_tmp["s.e."] <- object$theta_se
  }
  rownames(frame_tmp) <- names(object$theta)
  print(format(signif(frame_tmp, max(1L, digits - 2L))))
  cat(paste("\nlog-likelihood = "),object$logLike)
  cat(paste("\nAIC = "),object$AIC)
  cat(paste("\nBIC = "),object$BIC)
  
  resdist <- t(apply(object$resid,2,quantile))
  resdist <- data.frame(resdist)
  colnames(resdist) <- c("Min", "1Q", "Median", "3Q", "Max")
  rownames(resdist) <- paste0("Y",1:object$q)
  cat(paste("\n"))
  cat(paste("\nResiduals:\n"))
  print(format(signif(resdist, max(1L, digits - 2L))))
  
  invisible(object)
  
}

#' @title Summary of a \code{CHPTest} object
#'
#' @description This is a method for the function \code{summary()} for objects of the class \code{CHPTest}.
#' 
#' @inheritParams base::summary
#'
#' @return The \code{CHPTest} object is returned invisibly.
#' 
#' @keywords internal
#' 
#' @export
summary.HLRTest <- function(object, digits = getOption("digits"), ...){
  cat("\nRestricted Model\n")
  frame_tmp <- data.frame(coef = object$mdl_h0$theta)
  if (object$mdl_h0$control$getSE==TRUE){
    frame_tmp["s.e."] <- object$mdl_h0$theta_se
  }
  rownames(frame_tmp) <- names(object$mdl_h0$theta)
  print(format(signif(frame_tmp, max(1L, digits - 2L))))
  cat(paste("\nlog-likelihood = "),object$mdl_h0$logLike)
  cat(paste("\nAIC = "),object$mdl_h0$AIC)
  cat(paste("\nBIC = "),object$mdl_h0$BIC)
  cat("\n")
  if (object$control$msvar){
    cat("\nHansen (1992) Likelihood Ratio Bound Test -  Switch in Mean and Variance\n")
  }else{
    cat("\nHansen (1992) Likelihood Ratio Bound Test -  Switch in Mean only\n") 
  }
  out <- data.frame(cbind(object$LR0, object$LR_cv, object$pval))
  colnames(out) <- c("test-stat", colnames(object$LR_cv), "p-value")
  print(format(signif(out, max(1L, digits - 2L))))
  invisible(object)
}

#' @title Summary of a \code{CHPTest} object
#'
#' @description This is a method for the function \code{summary()} for objects of the class \code{CHPTest}.
#' 
#' @inheritParams base::summary
#'
#' @return The \code{CHPTest} object is returned invisibly.
#' 
#' @keywords internal
#' 
#' @export
summary.CHPTest <- function(object, digits = getOption("digits"), ...){
  cat("\nRestricted Model\n")
  frame_tmp <- data.frame(coef = object$mdl_h0$theta)
  if (object$mdl_h0$control$getSE==TRUE){
    frame_tmp["s.e."] <- object$mdl_h0$theta_se
  }
  rownames(frame_tmp) <- names(object$mdl_h0$theta)
  print(format(signif(frame_tmp, max(1L, digits - 2L))))
  cat(paste("\nlog-likelihood = "),object$mdl_h0$logLike)
  cat(paste("\nAIC = "),object$mdl_h0$AIC)
  cat(paste("\nBIC = "),object$mdl_h0$BIC)
  cat("\n")
  cat("\nCarrasco, Hu, & Ploberger (2014) Parameter Stability Test \n")
  if (object$control$msvar){
    cat("\n- Switch in Mean and Variance\n")
  }else{
    cat("\n- Switch in Mean only\n") 
  }
  out <- data.frame(rbind(c(object$supTS, object$supTS_cv, object$pval_supTS),
                          c(object$expTS, object$expTS_cv, object$pval_expTS)))
  colnames(out) <- c("test-stat", names(object$supTS_cv), "p-value")
  rownames(out) <- c("supTS", "expTS")
  print(format(signif(out, max(1L, digits - 2L))))
  invisible(object)
}

#' @title summaryummary of a \code{DLMCTest} object
#'
#' @description This is a method for the function \code{summary()} for objects of the class \code{DLMCTest}.
#' 
#' @inheritParams base::summary
#'
#' @return The \code{DLMCTest} object is returned invisibly.
#' 
#' @keywords internal
#' 
#' @export
summary.DLMCTest <- function(object, digits = getOption("digits"), ...){
  cat("\nRestricted Model\n") 
  frame_tmp <- data.frame(coef = object$mdl_h0$theta)
  if (object$mdl_h0$control$getSE==TRUE){
    frame_tmp["s.e."] <- object$mdl_h0$theta_se
  }
  rownames(frame_tmp) <- names(object$mdl_h0$theta)
  print(format(signif(frame_tmp, max(1L, digits - 2L))))
  cat(paste("\nlog-likelihood = "),object$mdl_h0$logLike)
  cat(paste("\nAIC = "),object$mdl_h0$AIC)
  cat(paste("\nBIC = "),object$mdl_h0$BIC)
  cat("\n")
  cat("\nDufour & Luger (2017) Moment-Based Local Monte Carlo Test\n")
  out <- data.frame(rbind(c(t(object$theta),object$S0, object$F0_min, object$FN_min_cv, object$pval_min),
                          c(t(object$theta),object$S0, object$F0_prod, object$FN_prod_cv, object$pval_prod)))
  colnames(out) <- c(rownames(object$theta), colnames(object$S0), colnames(object$F0_min), names(object$FN_min_cv), "p-value")
  rownames(out) <- c("LMC_min", "LMC_prod")
  print(format(signif(out, max(1L, digits - 2L))))
  invisible(object)
}

#' @title Summary of a \code{DLMMCTest} object
#'
#' @description This is a method for the function \code{summary()} for objects of the class \code{DLMMCTest}.
#' 
#' @inheritParams base::summary
#'
#' @return The \code{DLMMCTest} object is returned invisibly.
#' 
#' @keywords internal
#' 
#' @export
summary.DLMMCTest <- function(object, digits = getOption("digits"), ...){
  cat("\nRestricted Model\n") 
  frame_tmp <- data.frame(coef = object$mdl_h0$theta)
  if (object$mdl_h0$control$getSE==TRUE){
    frame_tmp["s.e."] <- object$mdl_h0$theta_se
  }
  rownames(frame_tmp) <- names(object$mdl_h0$theta)
  print(format(signif(frame_tmp, max(1L, digits - 2L))))
  cat(paste("\nlog-likelihood = "),object$mdl_h0$logLike)
  cat(paste("\nAIC = "),object$mdl_h0$AIC)
  cat(paste("\nBIC = "),object$mdl_h0$BIC)
  cat("\n")
  cat("\nDufour & Luger (2017) Moment-Based Maximized Monte Carlo Test\n")
  out <- data.frame(rbind(c(object$S0_min, object$F0_min, object$pval_min),
                          c(object$S0_prod, object$F0_prod, object$pval_prod)))
  colnames(out) <- c(colnames(object$S0_min), colnames(object$F0_min), "p-value")
  rownames(out) <- c("MMC_min","MMC_prod")
  print(format(signif(out, max(1L, digits - 2L))))
  invisible(object)
}

#' @title Summary of a \code{LMCLRTest} object
#'
#' @description This is a method for the function \code{summary()} for objects of the class \code{LMCLRTest}.
#' 
#' @inheritParams base::summary
#'
#' @return The \code{LMCLRTest} object is returned invisibly.
#' 
#' @keywords internal
#' 
#' @export
summary.LMCLRTest <- function(object, digits = getOption("digits"), ...){
  cat("\nRestricted Model\n") 
  frame_h0_tmp <- data.frame(coef = object$mdl_h0$theta)
  if (object$mdl_h0$control$getSE==TRUE){
    frame_h0_tmp["s.e."] <- object$mdl_h0$theta_se
  }
  rownames(frame_h0_tmp) <- names(object$mdl_h0$theta)
  print(format(signif(frame_h0_tmp, max(1L, digits - 2L))))
  cat(paste("\nlog-likelihood = "),object$mdl_h0$logLike)
  cat(paste("\nAIC = "),object$mdl_h0$AIC)
  cat(paste("\nBIC = "),object$mdl_h0$BIC)
  cat("\n")
  cat("\nUnrestricted Model\n") 
  frame_h1_tmp <- data.frame(coef = object$mdl_h1$theta)
  if (object$mdl_h1$control$getSE==TRUE){
    frame_h1_tmp["s.e."] <- object$mdl_h1$theta_se
  }
  rownames(frame_h1_tmp) <- names(object$mdl_h1$theta)
  print(format(signif(frame_h1_tmp, max(1L, digits - 2L))))
  cat(paste("\nlog-likelihood = "),object$mdl_h1$logLike)
  cat(paste("\nAIC = "),object$mdl_h1$AIC)
  cat(paste("\nBIC = "),object$mdl_h1$BIC)
  cat("\n")
  cat("\nRodriguez-Rondon & Dufour (2025) Local Monte Carlo Likelihood Ratio Test\n")
  out <- data.frame(t(as.matrix(c(object$LRT_0, object$LRN_cv, object$pval))))
  colnames(out) <- c(names(object$LRT_0), names(object$LRN_cv), "p-value")
  rownames(out) <- "LMC_LRT"
  print(format(signif(out, max(1L, digits - 2L))))
  invisible(object)
}



#' @title Summary of a \code{MMCLRTest} object
#'
#' @description This is a method for the function \code{summary()} for objects of the class \code{MMCLRTest}.
#' 
#' @inheritParams base::summary
#'
#' @return The \code{MMCLRTest} object is returned invisibly.
#' 
#' @keywords internal
#' 
#' @export
summary.MMCLRTest <- function(object, digits = getOption("digits"), ...){
  cat("\nRestricted Model\n") 
  frame_h0_tmp <- data.frame(coef = object$mdl_h0$theta)
  if (object$mdl_h0$control$getSE==TRUE){
    frame_h0_tmp["s.e."] <- object$mdl_h0$theta_se
  }
  rownames(frame_h0_tmp) <- names(object$mdl_h0$theta)
  print(format(signif(frame_h0_tmp, max(1L, digits - 2L))))
  cat(paste("\nlog-likelihood = "),object$mdl_h0$logLike)
  cat(paste("\nAIC = "),object$mdl_h0$AIC)
  cat(paste("\nBIC = "),object$mdl_h0$BIC)
  cat("\n")
  cat("\nUnrestricted Model\n") 
  frame_h1_tmp <- data.frame(coef = object$mdl_h1$theta)
  if (object$mdl_h1$control$getSE==TRUE){
    frame_h1_tmp["s.e."] <- object$mdl_h1$theta_se
  }
  rownames(frame_h1_tmp) <- names(object$mdl_h1$theta)
  print(format(signif(frame_h1_tmp, max(1L, digits - 2L))))
  cat(paste("\nlog-likelihood = "),object$mdl_h1$logLike)
  cat(paste("\nAIC = "),object$mdl_h1$AIC)
  cat(paste("\nBIC = "),object$mdl_h1$BIC)
  cat("\n")
  cat("\nRodriguez-Rondon & Dufour (2025) Maximized Monte Carlo Likelihood Ratio Test\n")
  out <- data.frame(t(as.matrix(c(object$LRT_0, object$pval))))
  colnames(out) <- c(names(object$LRT_0), "p-value")
  rownames(out) <- "MMC_LRT"
  print(format(signif(out, max(1L, digits - 2L))))
  invisible(object)
}


#' @title Summary of a \code{BootLRTest} object
#'
#' @description This is a method for the function \code{summary()} for objects of the class \code{BootLRTest}.
#' 
#' @inheritParams base::summary
#'
#' @return The \code{BootLRTest} object is returned invisibly.
#' 
#' @keywords internal
#' 
#' @export
summary.BootLRTest <- function(object, digits = getOption("digits"), ...){
  cat("\nRestricted Model\n") 
  frame_h0_tmp <- data.frame(coef = object$mdl_h0$theta)
  if (object$mdl_h0$control$getSE==TRUE){
    frame_h0_tmp["s.e."] <- object$mdl_h0$theta_se
  }
  rownames(frame_h0_tmp) <- names(object$mdl_h0$theta)
  print(format(signif(frame_h0_tmp, max(1L, digits - 2L))))
  cat(paste("\nlog-likelihood = "),object$mdl_h0$logLike)
  cat(paste("\nAIC = "),object$mdl_h0$AIC)
  cat(paste("\nBIC = "),object$mdl_h0$BIC)
  cat("\n")
  cat("\nUnrestricted Model\n") 
  frame_h1_tmp <- data.frame(coef = object$mdl_h1$theta)
  if (object$mdl_h1$control$getSE==TRUE){
    frame_h1_tmp["s.e."] <- object$mdl_h1$theta_se
  }
  rownames(frame_h1_tmp) <- names(object$mdl_h1$theta)
  print(format(signif(frame_h1_tmp, max(1L, digits - 2L))))
  cat(paste("\nlog-likelihood = "),object$mdl_h1$logLike)
  cat(paste("\nAIC = "),object$mdl_h1$AIC)
  cat(paste("\nBIC = "),object$mdl_h1$BIC)
  cat("\n")
  cat("\nBootstrap Likelihood Ratio Test\n")
  out <- data.frame(t(as.matrix(c(object$LRT_0, object$LRN_cv, object$pval))))
  colnames(out) <- c(names(object$LRT_0), names(object$LRN_cv), "p-value")
  rownames(out) <- "Boot_LRT"
  print(format(signif(out, max(1L, digits - 2L))))
  invisible(object)
}




#' @title Predict for a \code{Nmdl} object
#'
#' @description This is a method for the function \code{predict()} for objects of the class \code{Nmdl}.
#' 
#' @inheritParams stats::predict
#' @param h max number of prediction periods
#'
#' @return a \code{(h x q)} matrix with predicted value values.
#' 
#' @keywords internal
#' 
#' @export
predict.Nmdl <- function(object, ..., h = 10){
  # E[Y(t+h)|It]
  predictvals <- matrix(1,h,1)%*%t(as.matrix(object$mu))
  return(predictvals)
}

#' @title Predict for a \code{ARmdl} object
#'
#' @description This is a method for the function \code{predict()} for objects of the class \code{ARmdl}.
#' 
#' @inheritParams stats::predict
#' @param h max number of prediction periods
#'
#' @return a \code{(h x q)} matrix with predicted value values.
#' 
#' @keywords internal
#' 
#' @export
predict.ARmdl <- function(object, ..., h = 10){
  Fmat <- companionMat(c(object$phi), object$p, object$q) # Companion matrix
  if (object$p>1){ # Elimination vector
    Mn <- cbind(diag(object$q),matrix(0,object$q,object$q*(object$p-1)))
  }else if (object$p==1){
    Mn <- diag(object$q)
  }
  # E[Y(t)-mu|It]
  Yt <- t(object$x[object$n,,drop=F]) - c(object$mu)
  # E[Y(t+h)|It]
  predictvals <- matrix(0,h,object$q)
  for (xh in 1:h){
    predictvals[xh,] <- object$mu + Mn%*%((Fmat^xh)%*%Yt)
  }
  return(predictvals)
}


#' @title Predict for a \code{VARmdl} object
#'
#' @description This is a method for the function \code{predict()} for objects of the class \code{VARmdl}.
#' 
#' @inheritParams stats::predict
#' @param h max number of prediction periods
#'
#' @return a \code{(h x q)} matrix with predicted value values.
#' 
#' @keywords internal
#' 
#' @export
predict.VARmdl <- function(object, ..., h = 10){
  if (object$p>1){ # Elimination vector
    Mn <- cbind(diag(object$q),matrix(0,object$q,object$q*(object$p-1)))
  }else if (object$p==1){
    Mn <- diag(object$q)
  }
  # E[Y(t)-mu|It]
  Yt <- t(object$x[object$n,,drop=F]) - rep(object$mu,object$p)
  # E[Y(t+h)|It]
  predictvals <- matrix(0,h,object$q)
  for (xh in 1:h){
    predictvals[xh,] <- object$mu + Mn%*%((object$Fmat^xh)%*%Yt)
  }
  return(predictvals)
}


#' @title Predict for a \code{HMmdl} object
#'
#' @description This is a method for the function \code{predict()} for objects of the class \code{HMmdl}.
#' 
#' @inheritParams stats::predict
#' @param h max number of prediction periods
#'
#' @return a \code{(h x q)} matrix with predicted value values.
#' 
#' @keywords internal
#' 
#' @export
predict.HMmdl <- function(object, ..., h = 10){
  muy <- crossprod(object$mu,object$pinf) # long-run mean of process
  # E[Y(t+h)|It]
  predictvals <- matrix(0,h,object$q)
  predictxis <- matrix(0,h,object$k)
  for (xh in 1:h){
    predictxis_tmp <- object$pinf + (object$P^xh)%*%(t(object$St[object$n,,drop=F])-object$pinf)
    predictxis[xh,] <- predictxis_tmp/sum(predictxis_tmp) # normalized to sum to 1 (needed due to precision issue). 
    predictvals[xh,] <- muy + 
      crossprod(object$mu,(object$P^xh)%*%(t(object$St[object$n,,drop=F])-object$pinf))
  }
  return(list(predict = predictvals, predictSt = predictxis))
}

#' @title Predict for a \code{MSARmdl} object
#'
#' @description This is a method for the function \code{predict()} for objects of the class \code{MSARmdl}.
#' 
#' @inheritParams stats::predict
#' @param h max number of prediction periods
#'
#' @return a \code{(h x q)} matrix with predicted value values.
#' 
#' @keywords internal
#' 
#' @export
predict.MSARmdl <- function(object, ..., h = 10){
  Fmat <- companionMat(c(object$phi), object$p, object$q) # Companion matrix
  if (object$p>1){ # Elimination vector
    Mn <- cbind(diag(object$q),matrix(0,object$q,object$q*(object$p-1)))
  }else if (object$p==1){
    Mn <- diag(object$q)
  }
  muy <- crossprod(object$mu,object$pinf) # long-run mean of process
  # E[Y(t)-mu|It]
  Yt <- t(object$x[object$n,,drop=F]) - t(t(object$mu)%*%t(object$St[(object$n-object$p+1):object$n,,drop=F]))
  # E[Y(t+h)|It]
  predictvals <- matrix(0,h,object$q)
  predictxis <- matrix(0,h,object$k)
  for (xh in 1:h){
    predictxis_tmp <- object$pinf + (object$P^xh)%*%(t(object$St[object$n,,drop=F])-object$pinf)
    predictxis[xh,] <- predictxis_tmp/sum(predictxis_tmp) # normalized to sum to 1 (needed due to precision issue). 
    predictvals[xh,] <- muy + 
      crossprod(object$mu,(object$P^xh)%*%(t(object$St[object$n,,drop=F])-object$pinf)) + 
      Mn%*%((Fmat^xh)%*%Yt)
  }
  return(list(predict = predictvals, predictSt = predictxis))
}


#' @title Predict for a \code{MSVARmdl} object
#'
#' @description This is a method for the function \code{predict()} for objects of the class \code{MSVARmdl}.
#' 
#' @inheritParams stats::predict
#' @param h max number of prediction periods
#'
#' @return a \code{(h x q)} matrix with predicted value values.
#' 
#' @keywords internal
#' 
#' @export
predict.MSVARmdl <- function(object, ..., h = 10){
  if (object$p>1){ # Elimination vector
    Mn <- cbind(diag(object$q),matrix(0,object$q,object$q*(object$p-1)))
  }else if (object$p==1){
    Mn <- diag(object$q)
  }
  muy <- crossprod(object$mu,object$pinf) # long-run mean of process
  # E[Y(t)-mu|It]
  Yt <- t(object$x[object$n,,drop=F] - c(t(object$mu)%*%t(object$St[(object$n-object$p+1):object$n,,drop=F])))
  # E[Y(t+h)|It]
  predictvals <- matrix(0,h,object$q)
  predictxis <- matrix(0,h,object$k)
  for (xh in 1:h){
    predictxis_tmp <- object$pinf + (object$P^xh)%*%(t(object$St[object$n,,drop=F])-object$pinf)
    predictxis[xh,] <- predictxis_tmp/sum(predictxis_tmp) # normalized to sum to 1 (needed due to precision issue). 
    predictvals[xh,] <- muy + 
      crossprod(object$mu,(object$P^xh)%*%(t(object$St[object$n,,drop=F])-object$pinf)) + 
      Mn%*%((object$Fmat^xh)%*%Yt)
  }
  return(list(predict = predictvals, predictSt = predictxis))
}


#' @title Plot of a \code{simuNorm} object
#'
#' @description This is a method for the function \code{plot()} for objects of the class \code{simuNorm}.
#' 
#' @inheritParams base::plot
#'
#' @return The \code{simuNorm} object is returned invisibly.
#' 
#' @keywords internal
#' 
#' @export
plot.simuNorm <- function(x, ...){
  if (x$q==1){
    graphics::matplot(1:x$n, x$y, type = "l", ylab = "Simulated process", 
                      xlab ="Time", main = "Time series of simulated process")   
  }else if (x$q>1){
    graphics::matplot(1:x$n, x$y, type = "l", ylab = "Simulated processes", 
                      xlab ="Time", main = "Time series of simulated processes")   
  }
}


#' @title Plot of a \code{simuAR} object
#'
#' @description This is a method for the function \code{plot()} for objects of the class \code{simuAR}.
#' 
#' @inheritParams base::plot
#'
#' @return The \code{simuAR} object is returned invisibly.
#' 
#' @keywords internal
#' 
#' @export
plot.simuAR <- function(x, ...){
  graphics::matplot(1:x$n, x$y, type = "l", ylab = "Simulated process", 
                    xlab ="Time", main = "Time series of simulated process")   
}

#' @title Plot of a \code{simuARX} object
#'
#' @description This is a method for the function \code{plot()} for objects of the class \code{simuAR}.
#' 
#' @inheritParams base::plot
#'
#' @return The \code{simuARX} object is returned invisibly.
#' 
#' @keywords internal
#' 
#' @export
plot.simuARX <- function(x, ...){
  graphics::matplot(1:x$n, x$y, type = "l", ylab = "Simulated process", 
                    xlab ="Time", main = "Time series of simulated process")   
}


#' @title Plot of a \code{simuVAR} object
#'
#' @description This is a method for the function \code{plot()} for objects of the class \code{simuVAR}.
#' 
#' @inheritParams base::plot
#'
#' @return The \code{simuVAR} object is returned invisibly.
#' 
#' @keywords internal
#' 
#' @export
plot.simuVAR <- function(x, ...){
  graphics::matplot(1:x$n, x$y, type = "l", ylab = "Simulated processes", 
                    xlab ="Time", main = "Time series of simulated processes")   
}

#' @title Plot of a \code{simuVARX} object
#'
#' @description This is a method for the function \code{plot()} for objects of the class \code{simuVAR}.
#' 
#' @inheritParams base::plot
#'
#' @return The \code{simuVARX} object is returned invisibly.
#' 
#' @keywords internal
#' 
#' @export
plot.simuVARX <- function(x, ...){
  graphics::matplot(1:x$n, x$y, type = "l", ylab = "Simulated processes", 
                    xlab ="Time", main = "Time series of simulated processes")   
}


#' @title Plot of a \code{simuHMM} object
#'
#' @description This is a method for the function \code{plot()} for objects of the class \code{simuHMM}.
#' 
#' @inheritParams base::plot
#'
#' @return The \code{simuHMM} object is returned invisibly.
#' 
#' @keywords internal
#' 
#' @export
plot.simuHMM <- function(x, ...){
  oldpar <- graphics::par(no.readonly = TRUE)
  on.exit(graphics::par(oldpar)) 
  graphics::par(mfrow=c(2,1))
  if (x$q==1){
    graphics::matplot(1:x$n, x$y, type = "l", ylab = "Simulated process", 
                      xlab ="Time", main = "Time series of simulated process")   
    graphics::matplot(1:x$n, x$St+1, type = "l", ylab = "State (St)", xlab ="Time",)
  }else if (x$q>1){
    graphics::matplot(1:x$n, x$y, type = "l", ylab = "Simulated processes", 
                      xlab ="Time", main = "Time series of simulated processes")   
    graphics::matplot(1:x$n, x$St+1, type = "l", ylab = "State (St)", xlab ="Time",)
  }
}



#' @title Plot of a \code{simuMSAR} object
#'
#' @description This is a method for the function \code{plot()} for objects of the class \code{simuMSAR}.
#' 
#' @inheritParams base::plot
#'
#' @return The \code{simuMSAR} object is returned invisibly.
#' 
#' @keywords internal
#' 
#' @export
plot.simuMSAR <- function(x, ...){
  oldpar <- graphics::par(no.readonly = TRUE)
  on.exit(graphics::par(oldpar)) 
  graphics::par(mfrow=c(2,1))
  graphics::matplot(1:x$n, x$y, type = "l", ylab = "Simulated process", 
                    xlab ="Time", main = "Time series of simulated process")   
  graphics::matplot(1:x$n, x$St+1, type = "l", ylab = "State (St)", xlab ="Time",)
}


#' @title Plot of a \code{simuMSARX} object
#'
#' @description This is a method for the function \code{plot()} for objects of the class \code{simuMSAR}.
#' 
#' @inheritParams base::plot
#'
#' @return The \code{simuMSARX} object is returned invisibly.
#' 
#' @keywords internal
#' 
#' @export
plot.simuMSARX <- function(x, ...){
  oldpar <- graphics::par(no.readonly = TRUE)
  on.exit(graphics::par(oldpar)) 
  graphics::par(mfrow=c(2,1))
  graphics::matplot(1:x$n, x$y, type = "l", ylab = "Simulated process", 
                    xlab ="Time", main = "Time series of simulated process")   
  graphics::matplot(1:x$n, x$St+1, type = "l", ylab = "State (St)", xlab ="Time",)
}

#' @title Plot of a \code{simuMSVAR} object
#'
#' @description This is a method for the function \code{plot()} for objects of the class \code{simuMSVAR}.
#' 
#' @inheritParams base::plot
#'
#' @return The \code{simuMSVAR} object is returned invisibly.
#' 
#' @keywords internal
#' 
#' @export
plot.simuMSVAR <- function(x, ...){
  oldpar <- graphics::par(no.readonly = TRUE)
  on.exit(graphics::par(oldpar)) 
  graphics::par(mfrow=c(2,1))
  graphics::matplot(1:x$n, x$y, type = "l", ylab = "Simulated processes", 
                    xlab ="Time", main = "Time series of simulated processes")   
  graphics::matplot(1:x$n, x$St+1, type = "l", ylab = "State (St)", xlab ="Time",)
}

#' @title Plot of a \code{simuMSVARX} object
#'
#' @description This is a method for the function \code{plot()} for objects of the class \code{simuMSVAR}.
#' 
#' @inheritParams base::plot
#'
#' @return The \code{simuMSVARX} object is returned invisibly.
#' 
#' @keywords internal
#' 
#' @export
plot.simuMSVARX <- function(x, ...){
  oldpar <- graphics::par(no.readonly = TRUE)
  on.exit(graphics::par(oldpar)) 
  graphics::par(mfrow=c(2,1))
  graphics::matplot(1:x$n, x$y, type = "l", ylab = "Simulated processes", 
                    xlab ="Time", main = "Time series of simulated processes")   
  graphics::matplot(1:x$n, x$St+1, type = "l", ylab = "State (St)", xlab ="Time",)
}



#' @title Plot of a \code{Nmdl} object
#'
#' @description This is a method for the function \code{plot()} for objects of the class \code{Nmdl}.
#' 
#' @inheritParams base::plot
#'
#' @return The \code{Nmdl} object is returned invisibly.
#' 
#' @keywords internal
#' 
#' @export
plot.Nmdl <- function(x, ...){
  oldpar <- graphics::par(no.readonly = TRUE)
  on.exit(graphics::par(oldpar)) 
  graphics::par(mfrow=c(2,1))
  if (x$q==1){
    graphics::matplot(1:x$n, x$y, type = "l", ylab = "Observed process", 
                      xlab ="Time", main = "Time series of process & fitted values")   
    graphics::matplot(1:x$n, x$fitted, type = "l", ylab = "Fitted values", xlab ="Time")   
  }else if (x$q>1){
    graphics::matplot(1:x$n, x$y, type = "l", ylab = "Observed processes", 
                      xlab ="Time", main = "Time series of processes & fitted values")   
    graphics::matplot(1:x$n, x$fitted, type = "l", ylab = "Fitted values", xlab ="Time")   
  }
}


#' @title Plot of a \code{ARmdl} object
#'
#' @description This is a method for the function \code{plot()} for objects of the class \code{ARmdl}.
#' 
#' @inheritParams base::plot
#'
#' @return The \code{ARmdl} object is returned invisibly.
#' 
#' @keywords internal
#' 
#' @export
plot.ARmdl <- function(x, ...){
  oldpar <- graphics::par(no.readonly = TRUE)
  on.exit(graphics::par(oldpar)) 
  graphics::par(mfrow=c(2,1))
  graphics::matplot(1:x$n, x$y, type = "l", ylab = "Observed process", 
                    xlab ="Time", main = "Time series of process & fitted values")   
  graphics::matplot(1:x$n, x$fitted, type = "l", ylab = "Fitted values", xlab ="Time")   
}


#' @title Plot of a \code{VARmdl} object
#'
#' @description This is a method for the function \code{plot()} for objects of the class \code{VARmdl}.
#' 
#' @inheritParams base::plot
#'
#' @return The \code{VARmdl} object is returned invisibly.
#' 
#' @keywords internal
#' 
#' @export
plot.VARmdl <- function(x, ...){
  oldpar <- graphics::par(no.readonly = TRUE)
  on.exit(graphics::par(oldpar)) 
  graphics::par(mfrow=c(2,1))
  graphics::matplot(1:x$n, x$y, type = "l", ylab = "Observed processes", 
                    xlab ="Time", main = "Time series of processes & fitted values")   
  graphics::matplot(1:x$n, x$fitted, type = "l", ylab = "Fitted values", xlab ="Time")   
}




#' @title Plot of a \code{HMmdl} object
#'
#' @description This is a method for the function \code{plot()} for objects of the class \code{HMmdl}.
#' 
#' @inheritParams base::plot
#'
#' @return The \code{HMmdl} object is returned invisibly.
#' 
#' @keywords internal
#' 
#' @export
plot.HMmdl <- function(x, ...){
  oldpar <- graphics::par(no.readonly = TRUE)
  on.exit(graphics::par(oldpar)) 
  graphics::par(mfrow=c(3,1))
  if (x$q==1){
    graphics::matplot(1:x$n, x$y, type = "l", ylab = "Observed process", 
                      xlab ="Time", main = "Time series of process & fitted values")   
    graphics::matplot(1:x$n, x$fitted, type = "l", ylab = "Fitted values", xlab ="Time")   
    graphics::matplot(1:x$n, x$St, type = "l", ylab = "State (St)", xlab ="Time")   
  }else if (x$q>1){
    graphics::matplot(1:x$n, x$y, type = "l", ylab = "Observed processes", 
                      xlab ="Time", main = "Time series of processes & fitted values")   
    graphics::matplot(1:x$n, x$fitted, type = "l", ylab = "Fitted values", xlab ="Time")   
    graphics::matplot(1:x$n, x$St, type = "l", ylab = "State (St)", xlab ="Time")   
  }
}

#' @title Plot of a \code{MSARmdl} object
#'
#' @description This is a method for the function \code{plot()} for objects of the class \code{MSARmdl}.
#' 
#' @inheritParams base::plot
#'
#' @return The \code{MSARmdl} object is returned invisibly.
#' 
#' @keywords internal
#' 
#' @export
plot.MSARmdl <- function(x, ...){
  oldpar <- graphics::par(no.readonly = TRUE)
  on.exit(graphics::par(oldpar)) 
  graphics::par(mfrow=c(3,1))
  graphics::matplot(1:x$n, x$y, type = "l", ylab = "Observed process", 
                    xlab ="Time", main = "Time series of process & fitted values")   
  graphics::matplot(1:x$n, x$fitted, type = "l", ylab = "Fitted values", xlab ="Time")  
  graphics::matplot(1:x$n, x$St, type = "l", ylab = "State (St)", xlab ="Time")   
}



#' @title Plot of a \code{MSVARmdl} object
#'
#' @description This is a method for the function \code{plot()} for objects of the class \code{MSVARmdl}.
#' 
#' @inheritParams base::plot
#'
#' @return The \code{MSVARmdl} object is returned invisibly.
#' 
#' @keywords internal
#' 
#' @export
plot.MSVARmdl <- function(x, ...){
  oldpar <- graphics::par(no.readonly = TRUE)
  on.exit(graphics::par(oldpar)) 
  graphics::par(mfrow=c(3,1))
  graphics::matplot(1:x$n, x$y, type = "l", ylab = "Observed processes", 
                    xlab ="Time", main = "Time series of processes & fitted values")   
  graphics::matplot(1:x$n, x$fitted, type = "l", ylab = "Fitted values", xlab ="Time")  
  graphics::matplot(1:x$n, x$St, type = "l", ylab = "State (St)", xlab ="Time")   
}












