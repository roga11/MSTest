# -----------------------------------------------------------------------------
#' @title Calculate MC P-Value
#' 
#' @description This function calculates a Monte-Carlo p-value
#'
#' @param test_stat test statistic under the alternative (e.g. S_0)
#' @param null_vec series with test statistic under the null (i.e. vector S)
#' @param type like of test. options are: "geq" for right-tail test, "leq" for 
#' left-tail test, "abs" for absolute vallue test and "two-tail" for two-tail test.
#' 
#' @return MC p-value of test
#' 
#' @references Dufour, J. M. (2006). Monte Carlo tests with nuisance parameters: 
#' A general approach to finite-sample inference and nonstandard asymptotics. 
#' Journal of Econometrics, 133(2), 443-477.
#' @references Dufour, J. M., & Luger, R. (2017). Identification-robust moment-based 
#' tests for Markov switching in autoregressive models. Econometric Reviews, 36(6-9), 713-727.
#' 
#' @export
p_val <- function(test_stat, null_vec, type = c("geq", "leq", "abs", "two-tail")) {
  N         <- length(null_vec)
  u         <- runif(1 + sum(null_vec == test_stat))
  test_rank      <- sum(test_stat > null_vec) + sum(u <= u[1])
  # negative pvalue for min function
  survival_pval <- ((N+1-test_rank)/N)
  # Compute the p-value
  if (type == "absolute" || type == "geq") {
    pval <- (N * survival_pval + 1)/(N + 1)
    
  } else if (type == "leq") {
    pval <- (N * (1 - survival_pval) + 1)/(N + 1)
    
  } else if (type == "two-tailed") {
    pval <- 2 * min((N * (1 - survival_pval) + 1)/(N + 1), 
                    (N * survival_pval + 1)/(N + 1))
  }
  return(pval) 
}
# -----------------------------------------------------------------------------
#' @title Calculate MMC p-value
#'
#' @description This function calculates a Maximized Monte-Carlo p-value
#'
#' @param search_type Type of optimization algorithm when searching nuissance 
#' parameter space. Avaiable options are: GenSA, GA, PSO, randSearch_paramCI and 
#' gridSearch_paramCI. Default is set to randSearch_paramCI to match results in paper.
#' @param optimOptions List containing upper "upp" and lower "low" bounds of 
#' nuissance parameter space.
#' 
#' @return MMC p-value of test
#' 
#' @references Dufour, J. M. (2006). Monte Carlo tests with nuisance parameters: 
#' A general approach to finite-sample inference and nonstandard asymptotics. 
#' Journal of Econometrics, 133(2), 443-477.
#' @references Dufour, J. M., & Luger, R. (2017). Identification-robust moment-based 
#' tests for Markov switching in autoregressive models. Econometric Reviews, 36(6-9), 713-727.
#' 
#' @export
calc_mmcpval <- function(search_type = 'randSearch_paramCI', 
                         optimOptions = NULL){
  if (is.null(optimOptions)){
    low <- -0.99
    upp <- 0.99
  }else{
    low <- optimOptions$lower
    upp <- optimOptions$upper
  } 
  y     <- get('y', envir = MMCenv)
  x     <- get('x', envir = MMCenv)
  phi   <- get('phi', envir = MMCenv)
  se0   <- get('se0', envir = MMCenv)
  npar  <- get('npar', envir = MMCenv)
  nar   <- get('nar', envir = MMCenv)
  N     <- get('N', envir = MMCenv)
  N3    <- get('N3', envir = MMCenv)
  params<- get('params', envir = MMCenv)
  mcgrid <-matrix(nrow = 0, ncol=npar)
  # ---------------------- RANDOM SEARCH OVER PARAM CI ----------------------
  if (search_type=='randSearch_paramCI'){
    # random search over fixed interval from parameter confidence interval 
    while (nrow(mcgrid)<N3){
      ttp <- N3 - nrow(mcgrid)
      phitemp <- sweep(sweep(matrix(runif(ttp*npar),nrow=ttp,ncol=npar),
                             MARGIN=2,4*se0,'*'),MARGIN=2, phi-2*se0,'+')
      if ( nar>0){
        station <- which(apply(cbind(rep(1,nrow(phitemp)),
                                     matrix(-phitemp[,1: nar],
                                            nrow=nrow(phitemp),ncol=nar)),
                               MARGIN=1,function(x) 
                                 min(Mod(as.complex(polyroot(x)))))>=1)
        mcgrid  <- rbind(mcgrid,matrix(phitemp[station,],
                                       nrow=length(station),ncol=npar))
      }else{
        mcgrid  <- rbind(mcgrid,matrix(phitemp,nrow=nrow(phitemp),ncol=npar)) 
      }
      rm(phitemp)
    }
    output <- grid_eval(y,x,mcgrid,N,params)
    # ----------------------- GRID SEARCH OVER PARAM CI -----------------------
  }else if (search_type=='gridSearch_paramCI'){
    # grid search over fixed interval from parameter confidence interval 
    ttp <- N3 - nrow(mcgrid)
    lst <- vector("list", npar)
    for (xi in 1:npar){
      lst[[xi]] <- seq(from=phi[xi]-2*se0[xi],to=phi[xi]+2*se0[xi],
                       length.out = N3^(1/npar))
    }
    phitemp <- as.matrix(expand.grid(lst))
    if (nar>0){
      station <- which(apply(cbind(rep(1,nrow(phitemp)),
                                   matrix(-phitemp[,1: nar],
                                          nrow=nrow(phitemp),ncol=nar)),
                             MARGIN=1,function(x) 
                               min(Mod(as.complex(polyroot(x)))))>=1)
      mcgrid  <- rbind(mcgrid,matrix(phitemp[station,],
                                     nrow=length(station),ncol=npar))
    }else{
      mcgrid  <- rbind(mcgrid,matrix(phitemp,nrow=nrow(phitemp),ncol=npar)) 
    }
    rm(phitemp)
    output <- grid_eval(y,x,mcgrid,N,params)
    # ------------------ Simulated Annealing Opt. Algo. -----------------------   
  }else if (search_type=='GenSA'){
    Minres <- GenSA::GenSA(par=phi,fn=min_fn_min,lower = rep(low,length(phi)),
                           upper = rep(upp,length(phi)))
    Prodres <- GenSA::GenSA(par=phi,fn=min_fn_prod,lower = rep(low,length(phi)),
                            upper = rep(upp,length(phi)))
    output <- list(-Minres$value,Minres$par,-Prodres$value,Prodres$par)
    names(output) <- c('p-value_min','params_min','p-value_prod','params_prod')
    # ------------------------ Genetic Algorithm ----------------------------- 
  }else if (search_type=='GA'){
    Minres <- GA::ga(type='real-valued',max_fn_min,lower = rep(low,length(phi)),
                     upper = rep(upp,length(phi)))
    Prodres <- GA::ga(type='real-valued',
                      max_fn_prod,lower = rep(low,length(phi)),
                      upper = rep(upp,length(phi)))
    output <- list(Minres@fitnessValue,Minres@solution,
                   Prodres@fitnessValue,Prodres@solution)
    names(output) <- c('p-value_min','params_min','p-value_prod','params_prod')
    # ------------------------ Particle Swarm ------------------------------- 
  }else if (search_type=='PSO'){
    Minres <- pso::psoptim(par = phi, fn = min_fn_min, gr = NULL, 
                           lower = rep(low,length(phi)),
                           upper = rep(upp,length(phi)))
    Prodres <- pso::psoptim(par = phi , fn = min_fn_prod, gr = NULL,
                            lower = rep(low,length(phi)),
                            upper = rep(upp,length(phi)))
    output <- list(-Minres$value, Minres$par, -Prodres$value, Prodres$par)
    names(output) <- c('p-value_min','params_min','p-value_prod','params_prod')
  }
  return(output)
}
# -----------------------------------------------------------------------------
#' @title Simulate Auto-Regressive process
#' 
#' @description This function can be used to simulate a Auto-Regressive process
#'
#' @param n sample size of process
#' @param mu mean of process. This is the E[Y] and not the intercept of a model.
#' @param std standard deviation of process (i.e. \sigma, not \(\sigma^2\))
#' @param phi vector containing auto-regressive coefficients.
#' @param burnin This parameter can be used to set the number of burnin observations. 
#' Process with mean different from 0 can be very different at begining of simulation 
#' and so dropping some of the initial observations avoids these types of issues.
#' 
#' @return AR-MS series
#' 
#' @export
simuAR = function(n, mu = 0, std = 1, phi, burnin = 200){
  # -------------- Simulate series with autoregressive data generating process. 
  nphi = length(phi)
  series = matrix(0,n+burnin,1)
  series[1:nphi] = mu*(1-sum(phi)) + rnorm(nphi,0,std)
  for (xi in (nphi+1):(n+burnin)){
    series[xi] = mu*(1-sum(phi)) + t(as.matrix(phi))%*%as.matrix(series[(xi-1):(xi-nphi)]) + rnorm(1,0,std) 
  }
  # get rid of burnin values
  series = series[(burnin + 1):(burnin + n)]
  return(series)
}
# -----------------------------------------------------------------------------
#' @title Simulate Auto-Regressive Markov-Switching process
#' 
#' @description This function can be used to simulate a Auto-Regressive Markov-Switching process
#'
#' @param n sample size of process
#' @param th vector containing parameters of process. The first two values should 
#' be the two means, the next 2 are the limiting probabilities of being in each state 
#' (i.e. p and q), and the last two are the standard deviation in each state.
#' @param burnin This parameter can be used to set the number of burnin observations. 
#' Process with mean different from 0 can be very different at begining of simulation 
#' and so dropping some of the initial observations avoids these types of issues.
#' 
#' @return AR-MS series
#' 
#' @export
simuARMS = function(n, th, burnin = 200){
  mu0 = th[1]
  mu1 = th[2]
  p = th[3]      
  q = th[4]
  v0 = th[5]
  v1 = th[6]
  rho = (1-q)/(2-p-q)
  w = c(rho,1 - rho)
  phi = th[7:length(th)]
  #U = runif(n)
  randsamp = matrix(0, n + burnin, 1)
  randsamp[1:length(phi),1] = mu0*(1-sum(phi)) + rnorm(length(phi),0,v0)
  for(xi in (1+length(phi)):(n+burnin)){
    if(runif(1)<w[1]){
      randsamp[xi] = mu0*(1-sum(phi)) + t(as.matrix(randsamp[(xi-length(phi)):(xi-1)]))%*%as.matrix(phi) + rnorm(1,0,v0)
    }else{
      randsamp[xi] = mu1*(1-sum(phi)) + t(as.matrix(randsamp[(xi-length(phi)):(xi-1)]))%*%as.matrix(phi) + rnorm(1,0,v1)
    }
  }
  randsamp = randsamp[(burnin+1):(burnin+n)]
  return(randsamp)
}