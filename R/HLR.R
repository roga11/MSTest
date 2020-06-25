# ------------------------------------------------------------------------------------------
#' @title HLR Test
#'
#' This function ...
#'
#' @param oder is the order of Autoregressive and moving average components
#' @references Dufour and Luger (2017)
#' @export
HLRenv <- new.env()
HLRtest <- function(Y,p=1,ix=1,iv=0,iq=0,gn=20,reps=1000){
  # Estimate Model under Null Hypothesis
  odr     <- c(p,0)
  Mdl     <- estARMA(Y,odr)
  b0      <- as.matrix(c(Mdl$b0[sum(odr)+1]*(1-sum(Mdl$phi)),Mdl$phi))
  u       <- Mdl$u[(odr[1]+1):length(Mdl$u)]
  v0      <- sqrt(t(u)%*%u/length(u))
  b       <- c(b0,v0)
  x       <- cbind(1,Mdl$x)
  # Other Values
  k       <- length(x[1,])
  kx      <- length(ix)
  k1      <- kx+iv
  b1      <- matrix(0,k,1)
  p       <- 0
  q       <- 0
  pn      <- 12
  gp      <- seq(from=0.1,by=.075,to= (pn-1)*0.075+0.1)
  qn      <- pn
  gq      <- gp
  nwband  <- 4        # Set the Newey-West bandwidth 
  # Assign variables to HLR environment
  assign('nwband',nwband,env=HLRenv)
  assign('reps',reps,env=HLRenv)
  assign('ix',ix,env=HLRenv)
  assign('iv',iv,env=HLRenv)
  assign('iq',iq,env=HLRenv)
  assign('gn',gn,env=HLRenv)
  assign('y',Mdl$y,env=HLRenv)
  assign('x',x,env=HLRenv)
  assign('k',k,env=HLRenv)
  assign('kx',kx,env=HLRenv)
  assign('k1',k1,env=HLRenv)
  assign('b1',b1,env=HLRenv)
  assign('p',p,env=HLRenv)
  assign('q',q,env=HLRenv)
  # create Grid
  gmu <- as.matrix(seq(from=0.1,by=0.1,length.out = gn))
  gar <- matrix(rep(seq(from=-1,to=1,length.out = gn),kx-1),nrow=gn,ncol=kx-1)
  gsg <- matrix(rep(seq(from=0.1,by=0.1,length.out = gn),iv),nrow=gn,ncol=iv)
  gx  <- do.call(expand.grid,as.data.frame(cbind(gmu,gar,gsg)))
  # ---------- Calculation Under null # Regression 
  # LR from Null 
  null <- clike(b)
  # ---------- Perform Test
  out <- HLRparamSearch(gx,gp,gq,b,null)
  # ---------- Prepare output list
  m <- max(out$cs)
  resbeta <- as.matrix(out$beta[which.max(out$cs),])
  rownames(resbeta) <- paste0('ar',seq(1,(k-1)))
  draws <- t(out$draws)
  pvalue <- mean(apply(draws>m,2,mean))  #check if mean should be what we keep
  xi <- 0
  while (xi<=nwband){
    draws[,(xi+1)] <- sort(draws[,(xi+1)])
    xi <- xi+1
  }
  cr <- c(0.9, 0.95, 0.99)
  cr2 <- c('0.90 %', '0.95 %', '0.99 %')
  cu <- as.matrix(rowMeans(draws[round(cr*reps),]))
  rownames(cu) <- cr2
  output <- list(m,t(cu),pvalue,t(resbeta))
  names(output)<-c('Test-Stat','Crit-Values','p-value','params')
  return(output)
}
# ------------------------------------------------------------------------------------------
#' @title HLR param search
#'
#' This function ...
#'
#' @param oder is the order of Autoregressive and moving average components
#' @references Dufour and Luger (2017)
#' @export
HLRparamSearch <- function(gx,gp,gq,b,null){
  # load HLR environment variables
  k         <- HLRenv$k 
  y         <- HLRenv$y
  k1        <- HLRenv$k1
  iq        <- HLRenv$iq
  nwband    <- HLRenv$nwband
  reps      <- HLRenv$reps
  mnull     <- sum(null)
  
  gnn <- nrow(gx)
  cnum <- length(gp)*gnn
  if (iq == 0){
    cnum <- cnum*length(gq)
  }
  ny <- length(y)
  e <- matrix(rnorm((ny+nwband)*reps),nrow=(ny+nwband),ncol=reps)
  draws <- matrix(1,nrow=1+nwband,ncol=reps)*(-1000)
  c <- matrix(0,nrow=cnum,ncol=1)
  cs <- c
  paramval <- matrix(0,nrow=cnum,ncol=k-1)
  beta <- matrix(0,nrow=k+3+k1-iq,ncol=cnum)
  j <- 0
  for (i1 in 1:length(gp)){
    print(i1)
    p <- gp[i1]
    assign('p',p,env=HLRenv)
    for (i2 in 1:length(gq)){
      if (iq==1){
        q <- 1-p
        assign('q',q,env=HLRenv)
        i2 <- length(gq)  
      }else{
        q <- gq[i2]
        assign('q',q,env=HLRenv)
      }
      bs <- b
      for (xi in 1:nrow(gx)){
        j <- j+1
        b1 <- as.matrix(gx[xi,])
        assign('b1',b1,env=HLRenv)
        # Optimization 
        optlst <- optim(bs,mclike,dmclike,method='BFGS')
        bnew <- optlst$par
        paramval[j,] <- bnew[2:(length(bnew)-1)]
        f <- optlst$value
        beta[1:(k+1),j] <- bnew
        beta[(k+2):(k+1+k1),j] <- b1
        beta[(k+2+k1),j] <- p
        if (iq==0){
          beta[(k+3+k1),j] <- q
        }
        diff <- null - clike(bnew)
        diff <- diff - mean(diff)
        se <- as.numeric(sqrt(t(diff)%*%diff))
        diff <- (diff/se)
        c[j] <- mnull - f
        cs[j] <- ((mnull-f)/se)
        diffe <- matrix(0,nrow=1,ncol=reps)
        nw <- 0
        while (nw<=nwband){
          diffe <- diffe+t(diff)%*%e[(1+nw):(ny+nw),]
          draws[(1+nw),] <- apply(rbind(draws[(1+nw),],(diffe/sqrt(1+nw))),2,max)
          nw <- nw+1
        }
        bs <- bnew
      }
    }
  }
  output <- list(cs,draws,paramval)
  names(output) <- c('cs','draws','beta')
  return(output)
}
# ------------------------------------------------------------------------------------------
#' @title clike
#'
#' This function ...
#'
#' @param oder is the order of Autoregressive and moving average components
#' @export
clike <- function(b){
  # gets values at which likelihood function should be evaluated
  k   <- HLRenv$k 
  ix  <- HLRenv$ix
  kx  <- HLRenv$kx
  k1  <- HLRenv$k1
  b1  <- HLRenv$b1
  iv  <- HLRenv$iv
  iq  <- HLRenv$iq
  q   <- HLRenv$q
  p   <- HLRenv$p
  
  
  th2 <- matrix(0,(k+1),1)
  th2[ix] <- as.numeric(b1[1:kx])
  if (iv == 1){
    th2[k+1] <- as.numeric(b1[k1])
  }
  if (iq == 1){
    qs <- 1-p
  }else{
    qs = q 
  }
  ths <- c(b,th2,p,qs)
  lik <- marklike(ths)
  return(lik)
}
# ------------------------------------------------------------------------------------------
#' @title glike
#'
#' This function ...
#'
#' @param oder is the order of Autoregressive and moving average components
#' @export
glike <- function(th){
  # get global varibales 
  get("k", envir = .GlobalEnv)
  get("kx", envir = .GlobalEnv)
  get("ix", envir = .GlobalEnv)
  get("iq", envir = .GlobalEnv)
  
  
  th1 <- th[1:(k+1)]
  th2 <- matrix(0,nrow=k+1,ncol=1)
  th2[ix] <- th[(k+2):(k+1+kx)]
  r <- length(th)
  if (iq==1){
    p <- th[r]
    p <- exp(p)/(1+exp(p))
    q <- 1-p
  }else{
    p <- th[r-1]
    q <- th[r]
    p <- exp(p)/(1+exp(p))
    q <- exp(q)/(1+exp(q)) 
  }
  th3 <- c(th1,th2,p,q)
  lik <- marklike(th3)
  return(lik)  
}
# ------------------------------------------------------------------------------------------
#' @title dmclike
#'
#' This function ...
#'
#' @param oder is the order of Autoregressive and moving average components
#' @export
# Gradient optimation function i.e likelihood function.
dmclike <-function(th){
  # get global varibales 
  k   <- HLRenv$k 
  p   <- HLRenv$p
  q   <- HLRenv$q
  n   <- HLRenv$n
  x   <- HLRenv$x
  y   <- HLRenv$y
  nar <- HLRenv$nar
  kx  <- HLRenv$kx
  k1  <- HLRenv$k1
  ix  <- HLRenv$ix
  b1  <- HLRenv$b1
  iv  <- HLRenv$iv
  iq  <- HLRenv$iq
  
  
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
# ------------------------------------------------------------------------------------------
#' @title marklike
#'
#' This function ...
#'
#' @param oder is the order of Autoregressive and moving average components
#' @export
marklike <- function(ths){
  # computes likelihood at each time t
  k <- HLRenv$k
  y <- HLRenv$y
  x <- HLRenv$x
  
  
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
# ------------------------------------------------------------------------------------------
#' @title mglike
#'
#' This function ...
#'
#' @param oder is the order of Autoregressive and moving average components
#' @export
mglike <- function(th){
  lik <- sum(glike(th))  
  return(lik)
}
# ------------------------------------------------------------------------------------------
#' @title mclike
#'
#' This function ...
#'
#' @param oder is the order of Autoregressive and moving average components
#' @export
mclike <- function(th){
  # optimization function
  # gets likelihood value & value of gradient of likelihood.
  lik   <- sum(clike(th))
  return(lik)
}

