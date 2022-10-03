

#' @title HLR Test
#'
#' @param p: Number of AR(p) (default is 1).
#' @param ix: List of Markov Switching parameters. 1 = just mean c(1,2) = mean and first param, (default: 1). • iv: Idicates if Variance is also Markov Switching. 1 for true (default: 0).
#' @param iq: Indicator that bounds q by 1-p (default: 0).
#' @param gn: grid size for switching parameters (default: 20)
#' @param pn: grid size for p (and q) (default is: 12)
#' @param pn_step_size grid step size 0.075
#' @param nwband default is 4
#' @param rep: (default: 1000)
#'
#' @example /inst/examples/HLRTest_examples.R
#' @export
HLRTest <- function(Y, p = 1, control = list()){
  # ----- Set control values
  con <- list(ix = 1, 
              iv = 0, 
              iq = 1,
              gn = 20, 
              pn = 12,
              pn_step_size = 0.075,
              nwband = 0, 
              reps = 1000)
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
  HLR_opt_ls$k1   <- HLR_opt_ls$kx + con$iv
  HLR_opt_ls$b1   <- matrix(0, HLR_opt_ls$k, 1)
  HLR_opt_ls$p    <- 0
  HLR_opt_ls$q    <- 0
  #qn      <- pn
  # Create grids
  gp      <- seq(from = 0.1, by = con$pn_step_size, to = (con$pn-1)*con$pn_step_size+0.1)
  gq      <- gp
  
  gmu <- as.matrix(seq(from = 0.1, by = 0.1, length.out = con$gn))
  gar <- matrix(rep(seq(from = -1, to = 1, length.out = con$gn), HLR_opt_ls$kx-1), nrow = con$gn, ncol = HLR_opt_ls$kx-1)
  gsg <- matrix(rep(seq(from = 0.1, by = 0.1, length.out = con$gn), con$iv), nrow = con$gn, ncol = con$iv)
  gx  <- do.call(expand.grid, as.data.frame(cbind(gmu, gar, gsg)))
  # ---------- Calculation Under null # Regression 
  # LR from Null 
  b <- c(mdl_h0$coef,mdl_h0$stdev)
  null <- clike(b, HLR_opt_ls)
  # ---------- Perform Test
  out <- HLRparamSearch(gx, gp, gq, b, null, HLR_opt_ls)
  # ---------- Prepare output list
  m <- max(out$cs)
  resbeta <- as.matrix(out$ar_coef[which.max(out$cs),])
  rownames(resbeta) <- paste0('ar', seq(1, (HLR_opt_ls$k-1)))
  draws <- t(out$draws)
  pvalue <- mean(apply(draws>m,2,mean))  
  xi <- 0
  while (xi<=con$nwband){
    draws[,(xi+1)] <- sort(draws[,(xi+1)])
    xi <- xi+1
  }
  cr <- c(0.9, 0.95, 0.99)
  cr2 <- c('0.90 %', '0.95 %', '0.99 %')
  cu <- as.matrix(rowMeans(draws[round(cr*con$reps),]))
  rownames(cu) <- cr2
  output <- list(m,t(cu),pvalue,t(resbeta))
  names(output)<-c('Test-Stat','Crit-Values','p-value','params')
  return(output)
}

#' @title HLR param search
#'
#' @keywords internal
#' 
#' @export
HLRparamSearch <- function(gx, gp, gq, b, null, HLR_opt_ls){
  # load HLR environment variables
  k         <- HLR_opt_ls$k 
  y         <- HLR_opt_ls$y
  k1        <- HLR_opt_ls$k1
  iq        <- HLR_opt_ls$iq
  nwband    <- HLR_opt_ls$nwband
  reps      <- HLR_opt_ls$reps
  mnull     <- sum(null)
  
  gnn <- nrow(gx)
  cnum <- length(gp)*gnn
  if (iq == 0){
    cnum <- cnum*length(gq)
  }
  ny <- length(y)
  eps <- matrix(rnorm((ny+nwband)*reps), nrow = (ny+nwband), ncol = reps)
  draws <- matrix(1, nrow = 1+nwband, ncol = reps)*(-1000)
  c <- matrix(0, nrow = cnum, ncol = 1)
  cs <- c
  paramval <- matrix(0, nrow = cnum, ncol = k-1)
  beta <- matrix(0, nrow = k+3+k1-iq, ncol = cnum)
  j <- 0
  HLR_opt_ls_tmp <- HLR_opt_ls
  for (i1 in 1:length(gp)){
    print(i1)
    HLR_opt_ls_tmp$p <- gp[i1]
    i2 <- 1
    while (i2 <=length(gq)){
      if (iq == 1){
        HLR_opt_ls_tmp$q <- 1 - HLR_opt_ls_tmp$p
        i2 <- length(gq)+1  
      }else{
        HLR_opt_ls_tmp$q  <- gq[i2]
        i2 <- i2 + 1
      }
      bs <- b
      for (xi in 1:nrow(gx)){
        j <- j+1  
        HLR_opt_ls_tmp$b1 <- as.matrix(gx[xi,])
        # Optimization 
        optlst <- optim(bs, mclike, dmclike, HLR_opt_ls = HLR_opt_ls_tmp, method = 'BFGS')
        bnew <- optlst$par
        paramval[j,] <- bnew[2:(length(bnew)-1)]
        f <- optlst$value
        beta[1:(k+1),j] <- bnew
        beta[(k+2):(k+1+k1),j] <- HLR_opt_ls_tmp$b1
        beta[(k+2+k1),j] <- HLR_opt_ls_tmp$p
        if (iq==0){
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
        while (nw <= nwband){
          diffe <- diffe+t(diff)%*%eps[(1+nw):(ny+nw),]
          draws[(1+nw),] <- apply(rbind(draws[(1+nw),],(diffe/sqrt(1+nw))),2,max)
          nw <- nw+1
        }
        bs <- bnew
      }
    }
  }
  output <- list(cs, draws, paramval, beta)
  names(output) <- c('cs','draws','ar_coef', "beta")
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
  if (HLR_opt_ls$iv == 1){
    th2[k+1] <- as.numeric(HLR_opt_ls$b1[HLR_opt_ls$k1])
  }
  if (HLR_opt_ls$iq == 1){
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
#' @keywords internal
#' 
#' @export
# Gradient optimization function i.e likelihood function.
dmclike <-function(th, HLR_opt_ls){
  # get global varibales 
  k   <- HLR_opt_ls$k 
  p   <- HLR_opt_ls$p
  q   <- HLR_opt_ls$q
  n   <- HLR_opt_ls$n
  x   <- HLR_opt_ls$x
  y   <- HLR_opt_ls$y
  #nar <- HLR_opt_ls$nar
  kx  <- HLR_opt_ls$kx
  k1  <- HLR_opt_ls$k1
  ix  <- HLR_opt_ls$ix
  b1  <- HLR_opt_ls$b1
  iv  <- HLR_opt_ls$iv
  iq  <- HLR_opt_ls$iq
  
  
  c <- sqrt(2*pi)
  
  beta1 <- th[1:k]
  sig0 <- abs(th[k+1])
  sig1 <- sig0
  if (iv==1){
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
  dqq0 <- cbind(matrix(rep(z,length(x[1,])),nrow=length(z),ncol=length(x[1,]))*x,
                (z2/sig0 - sig0))
  dqq0 <- dqq0*matrix(rep(qq0/sig02,length(dqq0[1,])),nrow=length(qq0),
                      ncol=length(dqq0[1,]))
  
  
  z <- z - x%*%beta2
  z2 <- z^2
  qq1 <- exp(-z2/(2*sig12))/(sig1*c)
  dqq1 <- cbind(matrix(rep(z,length(x[1,])),nrow=length(z),ncol=length(x[1,]))*x,
                (z2/sig1 - sig1))
  
  dqq1 <- dqq1* matrix(rep(qq1/sig12,length(dqq0[1,])),nrow=length(qq1),
                       ncol=length(dqq0[1,]))
  
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

