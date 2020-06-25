# -----------------------------------------------------------------------------
#' @title AR model
#'
#' @param infile Path to the input file
#' @return A matrix of the infile
#' 
#' @export
ARmdl <- function(Y, ar=NULL){
  # -------------- AR(p) model (p=0 is allowed)
  if (is.null(ar)==T || (ar==0)){
    ar  <- 0
    y   <- as.matrix(Y)
    n   <- length(y)
    x    <- as.matrix(rep(1,n))
    X   <- x
    b0  <- solve(t(X)%*%X)%*%t(X)%*%y
    phi <- NULL
    mu0 <- b0[1]
    u0  <- y - X%*%b0
    v0  <- sqrt((t(u0)%*%u0)/(length(u0)))
    se  <- as.matrix(sqrt(diag(as.numeric(v0^2)*solve(t(X)%*%X))))
  }else{
    y   <- Y[(ar+1):length(Y)]
    x   <- matrix(0,length(y),0)
    for (xi in 1:ar){
      x <- cbind(x, Y[(ar+1-xi):(length(Y)-xi)])
    }
    X   <- cbind(1,x) # include intercept in data for when needed. 
    n   <- length(y) #length of y vector   
    b0  <- solve(t(X)%*%X)%*%t(X)%*%y
    phi <- b0[2:length(b0)]
    mu0 <- b0[1]/(1-sum(phi))
    u0  <- y - X%*%b0
    v0  <- sqrt((t(u0)%*%u0)/(length(u0)))
    se  <- as.matrix(sqrt(diag(as.numeric(v0^2)*solve(t(X)%*%X))))
  }
  logLike <- sum(log(dnorm(u0,mean(u0),v0)))
  ARmdl_out <- list()
  ARmdl_out[['y']] <- y
  ARmdl_out[['X']] <- X
  ARmdl_out[['x']] <- x
  ARmdl_out[['n']] <- n
  ARmdl_out[['ar']] <- ar
  ARmdl_out[['coef']] <- b0
  ARmdl_out[['phi']] <- phi
  ARmdl_out[['mu']] <- mu0
  ARmdl_out[['residuals']] <- u0
  ARmdl_out[['stdev']] <- v0
  ARmdl_out[['se']] <- se
  ARmdl_out[['logLike']] <- logLike
  return(ARmdl_out)
}
# -----------------------------------------------------------------------------
#' @title AR MS model
#'
#' @param infile Path to the input file
#' @return A matrix of the infile
#' @export
ARMSmdl <- function(Y, ar = NULL, k = 2, method="eval", mxit=500){
  if (is.null(ar)==T || (ar==0)){
    ar  <- 0
    y   <- as.matrix(Y)
    n   <- length(y)
    x    <- as.matrix(rep(1,n))
    X   <- x
    b0  <- solve(t(X)%*%X)%*%t(X)%*%y
    phi <- 0
  }else{
    y   <- Y[(ar+1):length(Y)]
    x   <- matrix(0,length(y),0)
    for (xi in 1:ar){
      x <- cbind(x, Y[(ar+1-xi):(length(Y)-xi)])
    }
    X   <- cbind(1,x) # include intercept in data for when needed. 
    n   <- length(y) #length of y vector   
    b0  <- solve(t(X)%*%X)%*%t(X)%*%y
    phi <- b0[2:length(b0)]
  }
  # transform model
  z <- y - x%*%phi
  # create list of relevant variables
  mdlz <- list()
  mdlz$y <- z
  mdlz$n <- length(z)
  mdlz$ar <- 0
  mdlz$x <- matrix(1,mdlz$n,1)
  mdlz$X <- mdlz$x
  mdlz$b0 <- mean(z)
  mdlz$v0 <- sd(z)
  mdlz$u0 <- z - mdlz$X%*%mdlz$b0
  # estimate MS mean and variance
  ARMSmdl_tmp <- msmdl_cpp(mdlz,k,method,mxit)
  # Create list with model output
  ARMSmdl_out <- list()
  ARMSmdl_out$y <- y
  ARMSmdl_out$X <- X
  ARMSmdl_out$x <- x
  ARMSmdl_out$n <- n
  ARMSmdl_out$ar <- ar
  if (all((phi>0))){
    ARMSmdl_out$coef <- cbind(ARMSmdl_tmp$MScoef,t(matrix(rep(phi,k),nrow = length(phi),ncol=k)))
  }else{
    ARMSmdl_out$coef <- ARMSmdl_tmp$MScoef
  }
  ARMSmdl_out$stdev <- ARMSmdl_tmp$thn[5:6]
  ARMSmdl_out$transMat <- rbind(c(ARMSmdl_tmp$thn[3],1-ARMSmdl_tmp$thn[3]),
                                c(1-ARMSmdl_tmp$thn[4],ARMSmdl_tmp$thn[4]))
  ARMSmdl_out$residuals <- ARMSmdl_tmp$resid
  # *** GET parameter standard errors
  ARMSmdl_out$logLike <- ARMSmdl_tmp$f0
  ARMSmdl_out$initVal <- ARMSmdl_tmp$initVal
  ARMSmdl_out$deltath <- ARMSmdl_tmp$deltath
  ARMSmdl_out$iterations <- ARMSmdl_tmp$it
  ARMSmdl_out$tol <- ARMSmdl_tmp$tol
  return(ARMSmdl_out)
}
# -----------------------------------------------------------------------------
#' @title ARMA model
#'
#' @param infile Path to the input file
#' @return A matrix of the infile
#' @export
estARMA <- function(Y,odr){
  odr <- cbind(odr[1],0,odr[2])
  mdl <- arima(Y,order=odr,include.mean = TRUE,method="CSS")
  phi <- mdl$coef[1:length(mdl$coef)-1]
  u   <- mdl$residuals
  b0  <- mdl$coef
  if (min(diag(mdl$var.coef))<0){
    stop(" Variance-covariance matrix is not p.s.d." )
  }else{
    se0 <- sqrt(diag(mdl$var.coef))[1:length(phi)] 
  }
  Tsize <- length(Y)
  maxlag <- max(odr)
  y <- Y[(1+maxlag):Tsize]
  xtj <- matrix(nrow=length(Y)-maxlag,ncol=0)
  if (odr[1]>0){
    for ( xi in 1:odr[1]){
      xtj <- cbind(xtj, Y[(maxlag-xi+1):(Tsize-xi)])
    } 
  }
  
  utj <- matrix(nrow=length(Y)-maxlag,ncol=0) 
  if (odr[3]>0){
    for ( xi in 1:odr[3]){
      utj <- cbind(utj, u[(maxlag-xi+1):(Tsize-xi)])
    }
  }
  x <- cbind(xtj,utj)
  # number of params & number of autoregressive components
  npar <- length(phi)
  nar  <- length(grep('ar',names(phi),value=TRUE))
  # roots of process
  roots <- polyroot( c(1,-phi[1:nar]))
  croots <- as.complex(roots)
  #Mod(croots)
  # check stationarity
  if (min(Mod(croots)) < 1){
    warning("roots of series are not outside unit circle")
  }
  output <- list(y,x,u,b0,phi,se0,npar,nar)
  names(output)<-c('y','x','u','b0','phi','se0','npar','nar')
  return(output)
}