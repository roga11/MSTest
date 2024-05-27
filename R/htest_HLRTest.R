

#' @title Hansen (1992) likelihood ratio test
#'
#' @description This function performs Hansen's likelihood ratio test as described in Hansen (1992).
#' Original source code can be found \href{https://www.ssc.wisc.edu/~bhansen/progs/jae_92.html}{here}.
#'
#' @param Y A (\code{T x 1}) matrix of observations.  
#' @param p Integer determining the number of autoregressive lags.
#' @param control List with test procedure options including: 
#' \itemize{
#'   \item{ix: }{List of Markov Switching parameters. 1 = just mean c(1,2) = mean and first param, (default: 1).}
#'   \item{msvar: }{Boolean indicator. If \code{TRUE}, there is a switch in variance. If \code{FALSE} only switch in mean is considered. Default is \code{FALSE}.}
#'   \item{qbound: }{Indicator that bounds q by 1-p (default: \code{FALSE}).}
#'   \item{gridsize: }{Integer determining the number of grid points for markov switching parameters. Default is \code{20}.}
#'   \item{pgrid_from: }{Double determining the initial grid point for transition probabilities. Default is \code{0.1}.}
#'   \item{pgrid_by: }{Double determining the step size for grid points of transition probabilities. This, along with \code{p_gridsize} will determine the bounds of search space. Default is \code{0.075}.}
#'   \item{pgrid_to: }{Double determining the end grid point for transition probabilities. Default is \code{0.925}.}
#'   \item{mugrid_from: }{Double determining the minimum value of mean in second regime. Default is \code{0.1}.}
#'   \item{mugrid_by: }{Double determining the step size for grid points of mean in second regime. This, along with \code{gridsize} will determine the max value of mean in second regime. Default is \code{0.1}.}
#'   \item{siggrid_from: }{Double determining the minimum value of sigma in second regime (if \code{msvar = TRUE}). Default is \code{0.1}.}
#'   \item{siggrid_by: }{Double determining the step size for grid points of sigma in second regime. This, along with \code{gridsize} will determine the max value of sigma in second regime. Default is \code{0.1}.}
#'   \item{N: }{Integer determining the number of replications. Default is \code{1000}.}
#'   \item{nwband: }{Integer determining maximum bandwidth in Bartlett kernel. Critical values and p-values are returned for each bandwidth from \code{0:nwband} as suggested in Hansen (1996). Default is \code{4}.}
#'   \item{theta_null_low: }{Vector determining lower bound on parameters under the null hypothesis. Length of vector should be number of model coefficients + 1 for variance. Default is to only bound variance at \code{0.01}.}
#'   \item{theta_null_upp: }{Vector determining upper bound on parameters under the null hypothesis. Length of vector should be number of model coefficients + 1 for variance. Default is to no bounds (i.e. \code{Inf}).}
#'   \item{optim_method: }{String determining the type of optimization procedure used. Allowed options are "gp-optim" for general purpose optimization using \code{\link{optim}} from \code{\link{stats}} or "nl-optim" using \code{\link{slsqp}} from \code{\link{nloptr}}. Default is "gp-optim".}
#' }
#' 
#' @return List of class \code{HLRTest} (\code{S3} object) with model attributes including: 
#' \itemize{
#'   \item{mdl_h0: }{List with restricted model attributes. This will be of class \code{ARmdl} (\code{S3} object). See \code{\link{ARmdl}}.}
#'   \item{LR0: }{Likelihood ratio test statistic value.}
#'   \item{LRN: }{A (\code{N x 1}) vector with simulated LRT statistics under null hypothesis.}
#'   \item{pval: }{P-value.}
#'   \item{LR_cv: }{A (\code{nwband x 3}) matrix with 90\%, 95\%, and 99\% critical values in each column respectively.}
#'   \item{coef: }{Vector of coefficients from restricted model and grid search that maximized standardized LRT. }
#'   \item{control: }{List with test procedure options used.}
#' }
#' 
#' @references Hansen, Bruce E. 1992. “The likelihood ratio test under nonstandard conditions: testing the Markov switching model of GNP.” \emph{Journal of applied Econometrics} 7 (S1): S61–S82.
#' @references Hansen, Bruce E. 1996. “Erratum: The likelihood ratio test under nonstandard conditions: testing the Markov switching model of GNP.” \emph{Journal of applied Econometrics} 7 (S1): S61–S82.
#' @example /inst/examples/HLRTest_examples.R
#' @export
HLRTest <- function(Y, p, control = list()){
  # ----- Set control values
  con <- list(ix = 1, 
              msvar = FALSE, 
              qbound = FALSE,
              gridsize = 20,
              pgrid_from = 0.1,
              pgrid_by = 0.075,
              pgrid_to = 0.925,
              mugrid_from = 0.1,
              mugrid_by  = 0.1,
              siggrid_from = 0.1,
              siggrid_by  = 0.1,
              N = 1000,
              nwband  = 4,
              optim_method = "gp-optim",
              theta_null_low = NULL,
              theta_null_upp = NULL)
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
  if(is.null(con$theta_null_low)){
    con$theta_null_low <- c(rep(-Inf, length(mdl_h0$coef)),0.01)
  }
  if(is.null(con$theta_null_upp)){
    con$theta_null_upp <- c(rep(Inf, length(mdl_h0$coef)+1))
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
  gp    <- seq(from = con$pgrid_from, by = con$pgrid_by, to = con$pgrid_to)
  gq    <- gp
  gmu   <- as.matrix(seq(from = con$mugrid_from, by = con$mugrid_by, length.out = con$gridsize))
  gar   <- matrix(rep(seq(from = -1, to = 1, length.out = con$gridsize), HLR_opt_ls$kx-1), nrow = con$gridsize, ncol = HLR_opt_ls$kx-1)
  gsig  <- matrix(rep(seq(from = con$siggrid_from, by = con$siggrid_by, length.out = con$gridsize), con$msvar), nrow = con$gridsize, ncol = con$msvar)
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
#' @description This function performs the parameter grid search needed for 
#' the likelihood ratio test described in Hansen (1992).
#'
#' @param gx matrix/grid containing values for switching parameters in second regime.
#' @param gp matrix/grid containing values for probability of going from regime 1 at (\code{t}) to regime 1 at (\code{t+1}).
#' @param gq matrix/grid containing values for probability of going from regime 2 at (\code{t}) to regime 2 at (\code{t+1}) (if not bounded to be \code{1-p} i.e., \code{qbound=FALSE}).
#' @param b vector of initial parameter values of restricted model. 
#' @param null vector with likelihood under the null hypothesis. 
#' @param HLR_opt_ls List with model properties and test controls defined in \code{HLRTest()}.
#' 
#' @return List which contains:
#' \itemize{
#'   \item{cs: }{Vector with standardized LRT statistic for each grid point.}
#'   \item{draws: }{List with a (\code{nwband+1 x N} matrix for each grid point. Each row of these matrices is a vector of simulated test statistics under the null hypothesis for a value of bandwidth .}
#'   \item{coefficients: }{A  matrix with coefficients for each grid point.}
#' }
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
  gnn       <- nrow(gx)
  cnum      <- length(gp)*gnn
  if (qbound == FALSE){
    cnum <- cnum*length(gq)
  }
  ny      <- length(y)
  eps     <- matrix(stats::rnorm((ny+nwband)*reps), nrow=(ny+nwband), ncol=reps)
  draws   <- matrix(1,nrow=1+nwband,ncol=reps)*(-1000)
  drawsLs <- list()
  c       <- matrix(0, nrow = cnum, ncol = 1)
  cs      <- c
  beta    <- matrix(0, nrow = k + 3 + k1 - qbound, ncol = cnum) # All parameters 
  j       <- 0
  HLR_opt_ls_tmp <- HLR_opt_ls
  lowb <- HLR_opt_ls$theta_null_low
  uppb <- HLR_opt_ls$theta_null_upp
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
        #optlst <- stats::optim(bs, mclike, dmclike, HLR_opt_ls = HLR_opt_ls_tmp, method = 'BFGS')
        if (HLR_opt_ls$optim_method == "gp-optim"){
          optlst <- NULL
          try(
            optlst <- stats::optim(bs, mclike, dmclike, HLR_opt_ls = HLR_opt_ls_tmp, method = 'L-BFGS-B', lower = lowb, upper = uppb)   
          )
          if (is.null(optlst)){
            optlst  <- nloptr::slsqp(bs, fn = mclike, gr = dmclike, lower = lowb, upper = uppb, HLR_opt_ls = HLR_opt_ls_tmp)  
          }
        }else if(HLR_opt_ls$optim_method == "nl-optim"){
          optlst  <- nloptr::slsqp(bs, fn = mclike, gr = dmclike, lower = lowb, upper = uppb, HLR_opt_ls = HLR_opt_ls_tmp)  
        }
        bnew    <- optlst$par
        f       <- optlst$value
        beta[1:(k+1),j] <- bnew
        beta[(k+2):(k+1+k1),j] <- HLR_opt_ls_tmp$b1
        beta[(k+2+k1),j] <- HLR_opt_ls_tmp$p
        if (qbound == FALSE){
          beta[(k+3+k1),j] <- HLR_opt_ls_tmp$q
        }
        diff <- null - clike(bnew, HLR_opt_ls_tmp)
        diff <- diff - mean(diff)
        se <- as.numeric(sqrt(crossprod(diff)))
        diff <- (diff/se)
        c[j] <- mnull - f
        cs[j] <- ((mnull-f)/se)
        diffe <- matrix(0, nrow = 1, ncol = reps)
        nw <- 0
        while (nw<=nwband){
          diffe <- diffe+crossprod(diff,eps[(1+nw):(ny+nw),])
          draws[(1+nw),] <- apply(rbind(draws[(1+nw),],(diffe/sqrt(1+nw))),2,max)
          nw <- nw+1
        }
        drawsLs[[j]] <- draws
        #bs <- bnew
      }
    }
  } 
  output <- list(cs = cs, draws = drawsLs, coefficients = t(beta))
  return(output)
}

#' @title Parameter vector & likelihood function used by \code{HLRTest()}
#'
#' @description  This function combines parameters of restricted model 
#' with parameters of unrestricted model and then computes the likelihood using 
#' \code{marklike()}.
#' 
#' @param b vector of parameters from restricted model.
#' @param HLR_opt_ls List with model properties and test controls defined in \code{HLRTest()}.
#'
#' @return Value of likelihood at each period \code{t}.
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

#' @title Gradient of likelihood function.
#'
#' @description this function computes the score vector.
#' 
#' @param th vector of parameter values.
#' @param HLR_opt_ls List with model properties and test controls defined in \code{HLRTest()}.
#' 
#' @return vector with gradient of likelihood function for each parameter. Used in \code{HLRpramSearch()}.
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

#' @title Likelihood function used by \code{HLRTest()}
#'
#' @description this function computes the sum Markov likelihood
#' 
#' @param ths vector of parameter values.
#' @param HLR_opt_ls List with model properties and test controls defined in \code{HLRTest()}.
#' 
#' @return Vector of likelihood values.
#' 
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


#' @title Sum of likelihood used by \code{HLRTest()}
#'
#' @description This function computes the sum of the likelihood. 
#' 
#' @param th vector of parameter values.
#' @param HLR_opt_ls List with model properties and test controls defined in \code{HLRTest()}.
#' 
#' @return Sum of likelihood.
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

