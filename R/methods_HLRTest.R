

#' @title Hansen (1992) likelihood ratio test
#'
#' @description This function performs Hansen's likelihood ratio test as described in Hansen (1992).
#' Original source code can be found \href{https://www.ssc.wisc.edu/~bhansen/progs/jae_92.html}{here}.
#'
#' @param Y A (\code{T x 1}) matrix of observations.  
#' @param p Integer determining the number of autoregressive lags.
#' @param control List with test procedure options including: 
#' \itemize{
#'   \item{\code{ix}: }{List of Markov Switching parameters. 1 = just mean c(1,2) = mean and first param, (default: 1). • iv: Idicates if Variance is also Markov Switching. 1 for true (default: 0).}
#'   \item{\code{msvar}: }{Boolean indicator. If \code{TRUE}, there is a switch in variance. If \code{FALSE} only switch in mean is considered.}
#'   \item{\code{qbound}: }{Indicator that bounds q by 1-p (default: 0).}
#'   \item{\code{gridsize}: }{}
#'   \item{\code{p_gridsize}: }{}
#'   \item{\code{p_stepsize}: }{}
#'   \item{\code{mugrid_from}: }{}
#'   \item{\code{mugrid_by}: }{}
#'   \item{\code{siggrid_from}: }{}
#'   \item{\code{siggrid_by}: }{}
#'   \item{\code{N}: }{}
#'   \item{\code{sig_min}: }{}
#' }
#' 
#' @return List of class \code{HLRTest} (\code{S3} object) with model attributes including: 
#' \itemize{
#'   \item{\code{mdl_h0}: }{List with restricted model attributes. This will be of class \code{ARmdl} (\code{S3} object). See \code{\link{ARmdl}}.}
#'   \item{\code{supTS}: }{supTS test statistic value.}
#'   \item{\code{supTS_N}: }{A (\code{N x 1}) vector with simulated supTS test statistics under null hypothesis.}
#'   \item{\code{pval_supTS}: }{P-value for supTS version of parameter stability test.}
#'   \item{\code{supTS_cv}: }{Vector with 90\%, 95\%, and 99\% bootstrap critical values for supTS version of parameter stability test.}
#'   \item{\code{control}: }{List with test procedure options used.}
#' }
#' 
#' @references Hansen, Bruce E. 1992. “The likelihood ratio test under nonstandard conditions: testing the Markov switching model of GNP.” \emph{Journal of applied Econometrics} 7 (S1): S61–S82.
#'
#' @example /inst/examples/HLRTest_examples.R
#' @export
HLRTest <- function(Y, p, control = list()){
  # ----- Set control values
  con <- list(ix = 1, 
              msvar = FALSE, 
              qbound = FALSE,
              gridsize = 20,
              p_gridsize = 12,
              p_stepsize = 0.075,
              mugrid_from = 0.1,
              mugrid_by  = 0.1,
              siggrid_from = 0.1,
              siggrid_by  = 0.1,
              N = 1000,
              nwband  = 4,
              sig_min = 0.01)
  # Perform some checks for controls
  nmsC <- names(con)
  con[(namc <- names(control))] <- control
  if(length(noNms <- namc[!namc %in% nmsC])){
    warning("unknown names in control: ", paste(noNms,collapse=", ")) 
  }
  # Estimate Model under Null Hypothesis
  if(p>0){
    mdl_h0 <- ARmdl(Y, p, list(const = TRUE, getSE = TRUE))  
  }else if(p == 0){
    mdl_h0 <- Nmdl(Y, list(const = TRUE, getSE = TRUE))
  }else{
    stop("Value for 'p' must be an integer >=0.")
  }
  # Optimization values
  HLR_opt_ls <- con
  HLR_opt_ls$y    <- mdl_h0$y
  HLR_opt_ls$x    <- mdl_h0$X
  HLR_opt_ls$k    <- ncol(mdl_h0$X) 
  HLR_opt_ls$kx   <- length(con$ix)
  HLR_opt_ls$k1   <- HLR_opt_ls$kx + con$msvar
  HLR_opt_ls$b1   <- matrix(0, HLR_opt_ls$k, 1)
  HLR_opt_ls$p    <- 0 # under null hypothesis
  HLR_opt_ls$q    <- 0 # under null hypothesis
  # Create grids
  gp    <- seq(from = 0.1, by = con$p_stepsize, to = (con$p_gridsize-1)*con$p_stepsize+0.1)
  gq    <- gp
  gmu   <- as.matrix(seq(from = con$mugrid_from, by = con$mugrid_by, length.out = con$gridsize))
  gar   <- matrix(rep(seq(from = -1, to = 1, length.out = con$gridsize), HLR_opt_ls$kx-1), nrow = con$gridsize, ncol = HLR_opt_ls$kx-1)
  gsig  <- matrix(rep(seq(from = con$mugrid_from, by = con$mugrid_by, length.out = con$gridsize), con$msvar), nrow = con$gridsize, ncol = con$msvar)
  gx    <- do.call(expand.grid, as.data.frame(cbind(gmu, gar, gsig)))
  # ---------- Calculation Under null # Regression 
  # LR from Null 
  b <- c(mdl_h0$coef,mdl_h0$stdev)
  null <- clike(b, HLR_opt_ls)
  # ---------- Perform Test
  out <- HLRparamSearch(gx, gp, gq, b, null, HLR_opt_ls)
  # ---------- Prepare output list
  LR_stat <- max(out$cs)
  coef <- as.matrix(out$coefficients[which.max(out$cs),])
  draws <- apply(t(out$draws[[which.max(out$cs)]]),2,sort)
  if (con$msvar == FALSE){
    signame_tmp_0 <-  "sig" 
    signame_tmp_1 <-  NULL
  }else{
    signame_tmp_0 <- "sig_0"  
    signame_tmp_1 <- "sig_1"  
  }
  if (length(con$ix)>1){
    rname_tmp <- c("const_0", paste0("ar", seq(1, (HLR_opt_ls$k-1)),"_0"), signame_tmp_0,
                   "const_1", paste0("ar", con$ix[2:length(con$ix)]-1,"_1"), signame_tmp_1)
  }else{
    rname_tmp <- c("const_0", paste0("ar", seq(1, (HLR_opt_ls$k-1)),"_0"), signame_tmp_0,
                   "const_1", signame_tmp_1)
  }
  rname_tmp <- c(rname_tmp, "p")
  if (con$qbound == FALSE){
    rname_tmp <- c(rname_tmp, "q")
  }
  rownames(coef) <- rname_tmp
  pvalue <- as.matrix(colMeans(draws>LR_stat))
  cr <- c(0.9, 0.95, 0.99)
  cu <- t(as.matrix(draws[round(cr*con$N),]))
  colnames(cu) <- c('0.90 %', '0.95 %', '0.99 %')
  rownames(cu) <- paste0("M = ",0:(con$nwband))
  rownames(pvalue) <- paste0("M = ",0:(con$nwband))
  colnames(draws) <- paste0("M = ",0:(con$nwband))
  HLRTest_output <- list(mdl_h0 = mdl_h0, LR0 = LR_stat, LRN = draws, LR_cv = cu, 
                         pval = pvalue, coef = coef, control = con)
  class(HLRTest_output) <- "HLRTest"
  return(HLRTest_output)
}

#' @title HLR param search
#'
#' @keywords internal
#' 
#' @export
HLRparamSearch <- function(gx, gp, gq, b, null, HLR_opt_ls){
  k         <- HLR_opt_ls$k 
  y         <- HLR_opt_ls$y
  k1        <- HLR_opt_ls$k1
  qbound    <- HLR_opt_ls$qbound
  reps      <- HLR_opt_ls$N
  nwband    <- HLR_opt_ls$nwband
  mnull     <- sum(null)
  gnn <- nrow(gx)
  cnum <- length(gp)*gnn
  if (qbound == FALSE){
    cnum <- cnum*length(gq)
  }
  ny <- length(y)
  eps <- matrix(stats::rnorm((ny+nwband)*reps), nrow=(ny+nwband), ncol=reps)
  draws <- matrix(1,nrow=1+nwband,ncol=reps)*(-1000)
  drawsLs <- list()
  c <- matrix(0, nrow = cnum, ncol = 1)
  cs <- c
  beta <- matrix(0, nrow = k + 3 + k1 - qbound, ncol = cnum) # All parameters 
  j <- 0
  HLR_opt_ls_tmp <- HLR_opt_ls
  lowb <- rep(-Inf, length(b))
  lowb[length(b)] <- HLR_opt_ls$sig_min
  uppb <- rep(Inf, length(b))
  for (i1 in 1:length(gp)){
    HLR_opt_ls_tmp$p <- gp[i1]
    i2 <- 1
    while (i2 <=length(gq)){
      if (qbound == TRUE){
        HLR_opt_ls_tmp$q <- 1 - HLR_opt_ls_tmp$p
        i2 <- length(gq)+1  
      }else{
        HLR_opt_ls_tmp$q  <- gq[i2]
        i2 <- i2 + 1
      }
      bs <- b # starting values for optim
      for (xi in 1:nrow(gx)){
        j <- j+1  
        HLR_opt_ls_tmp$b1 <- as.matrix(gx[xi,]) #value for constant in regime 2
        # Optimization (optimal theta which contains constant in regime 1, ar coefs, and variance)
        #optlst <- optim(bs, mclike, dmclike, HLR_opt_ls = HLR_opt_ls_tmp, method = 'BFGS')
        optlst <- stats::optim(bs, mclike, dmclike, HLR_opt_ls = HLR_opt_ls_tmp, method = 'L-BFGS-B', lower = lowb, upper = uppb)
        bnew <- optlst$par
        f <- optlst$value
        beta[1:(k+1),j] <- bnew
        beta[(k+2):(k+1+k1),j] <- HLR_opt_ls_tmp$b1
        beta[(k+2+k1),j] <- HLR_opt_ls_tmp$p
        if (qbound == FALSE){
          beta[(k+3+k1),j] <- HLR_opt_ls_tmp$q
        }
        diff <- null - clike(bnew, HLR_opt_ls_tmp)
        diff <- diff - mean(diff)
        se <- as.numeric(sqrt(t(diff)%*%diff))
        diff <- (diff/se)
        c[j] <- mnull - f
        cs[j] <- ((mnull-f)/se)
        diffe <- matrix(0, nrow = 1, ncol = reps)
        nw <- 0
        while (nw<=nwband){
          diffe <- diffe+t(diff)%*%eps[(1+nw):(ny+nw),]
          draws[(1+nw),] <- apply(rbind(draws[(1+nw),],(diffe/sqrt(1+nw))),2,max)
          nw <- nw+1
        }
        drawsLs[[j]] <- draws
        bs <- bnew
      }
    }
  }
  output <- list(cs = cs, draws = drawsLs, coefficients = t(beta))
  return(output)
}

#' @title clike
#'
#' @keywords internal
#' 
#' @export
clike <- function(b, HLR_opt_ls){
  th2 <- matrix(0, (HLR_opt_ls$k+1), 1)
  th2[HLR_opt_ls$ix] <- as.numeric(HLR_opt_ls$b1[1:HLR_opt_ls$kx])
  if (HLR_opt_ls$msvar == TRUE){
    th2[HLR_opt_ls$k+1] <- as.numeric(HLR_opt_ls$b1[HLR_opt_ls$k1])
  }
  if (HLR_opt_ls$qbound == TRUE){
    qs <- 1 - HLR_opt_ls$p
  }else{
    qs = HLR_opt_ls$q
  }
  ths <- c(b, th2, HLR_opt_ls$p, qs)
  lik <- marklike(ths, HLR_opt_ls)
  return(lik)
}

#' @title dmclike
#'
#' @description  Gradient optimization function i.e markov likelihood function.
#'
#' @keywords internal
#' 
#' @export
dmclike <-function(th, HLR_opt_ls){
  k   <- HLR_opt_ls$k 
  p   <- HLR_opt_ls$p
  q   <- HLR_opt_ls$q
  n   <- HLR_opt_ls$n
  x   <- HLR_opt_ls$x
  y   <- HLR_opt_ls$y
  kx  <- HLR_opt_ls$kx
  k1  <- HLR_opt_ls$k1
  ix  <- HLR_opt_ls$ix
  b1  <- HLR_opt_ls$b1
  iv  <- HLR_opt_ls$msvar
  
  
  c <- sqrt(2*pi)
  
  beta1 <- th[1:k]
  sig0 <- abs(th[k+1])
  sig1 <- sig0
  if (iv == TRUE){
    sig1 <- sig1 + abs(b1[k1])
  }
  beta2 <- matrix(0,k,1)
  beta2[ix] <- b1[1:kx]
  sig02 <- sig0^2
  sig12 <- sig1^2
  qs <- 1 - q
  
  pp <- qs/(1-p+qs)
  
  z <- y - x%*%beta1
  z2 <- z^2
  qq0 <- exp(-z2/(2*sig02))/(sig0*c)
  dqq0 <- cbind((z%*%t(as.matrix(rep(1,ncol(x)))))*x, (z2/sig0 - sig0))
  dqq0 <- dqq0*matrix(rep(qq0/sig02,length(dqq0[1,])),nrow=length(qq0),ncol=length(dqq0[1,]))
  
  z <- z - x%*%beta2
  z2 <- z^2
  qq1 <- exp(-z2/(2*sig12))/(sig1*c)
  dqq1 <- cbind((z%*%t(as.matrix(rep(1,ncol(x)))))*x, (z2/sig1 - sig1))
  dqq1 <- dqq1* matrix(rep(qq1/sig12,length(dqq0[1,])),nrow=length(qq1),ncol=length(dqq0[1,]))
  
  n <- length(y)
  fit <- matrix(0,1,1+k)
  dpp <- matrix(0,1,1+k)
  
  it <- 1
  while (it<=n){
    p1 <- qs + (p - qs)*pp
    f1 <- qq1[it]*p1
    ff <- f1 + qq0[it]*(1-p1)
    pp <- f1/ff
    
    dp1 <- (p-qs)*dpp
    df1 <- dqq1[it,]*p1 + qq1[it]*dp1
    dff <- (df1 + dqq0[it,]*(1-p1) - qq0[it]*dp1)/ff
    dpp <- (df1 - f1*dff)/ff
    fit <- fit + dff
    
    it <- it+1
  }
  nfit <- -fit
  return(nfit)
}

#' @title marklike
#'
#' @description Markov likelihood
#' @keywords internal
#' 
#' @export
marklike <- function(ths, HLR_opt_ls){
  # computes likelihood at each time t
  k <- HLR_opt_ls$k
  y <- HLR_opt_ls$y
  x <- HLR_opt_ls$x
  c <- sqrt(2*pi)
  sig0 <- abs(ths[k+1])
  sig1 <- sig0 + abs(ths[2*k+2])
  p <- ths[2*k+3]
  q1 <- 1 - ths[2*k+4]
  pp <- q1/(1-p+q1)
  
  z <- y - x%*%ths[1:k]
  qq0 <- exp(-(z^2)/(2*sig0*sig0))/(sig0*c)
  
  z <- z - x%*%ths[(k+2):(2*k+1)]
  qq1 <- exp(-(z^2)/(2*sig1*sig1))/(sig1*c)
  rm(z)
  
  n <- length(y)
  fit <- matrix(0,n,1)
  it <- 1
  while (it<=n){
    p1 <- q1 + (p - q1)*pp
    f1 <- qq1[it]*p1
    ff <- f1 + qq0[it]*(1-p1)
    pp <- f1/ff
    fit[it] <- ff
    it <- it+1
  }
  
  # some values of b0 in optimization produce qq1 and/or qq0 =0 and result
  # in pp =NaN and everything after also NaN. This line eliminates these.
  fit[is.nan(fit)]<- -9999
  
  if (min(fit)>0){
    f <- log(fit)  
  }else{
    f <- -1000
  }
  nf <- -f
  return(nf)  
}


#' @title mclike
#'
#' @keywords internal
#' 
#' @export
mclike <- function(th, HLR_opt_ls){
  # optimization function
  # gets likelihood value & value of gradient of likelihood.
  lik   <- sum(clike(th, HLR_opt_ls))
  return(lik)
}

