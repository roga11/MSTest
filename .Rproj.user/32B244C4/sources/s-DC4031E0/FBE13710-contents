LRMMCenv <- new.env()
# -----------------------------------------------------------------------------
#' @title LRMC Simulation
#' This function ...
#'
#' @param infile Path to the input file
#' @return A matrix of the infile
#' @export
LRMC_simu = function(mdl,k,N){
  LRT_N = matrix(0,N,1)
  if (mdl$ar>0){
    phi = mdl$coef[2:length(mdl$coef)]
  }else{
    phi = 0
  }
  for (xxi in 1:N){
    y0 = as.matrix(simu_AR_dgp(n = mdl$n, mu = mean(mdl$y), std = mdl$stdev,phi))
    mdl_h0_tmp = ARmdl(y0, mdl$ar)
    mdl_h1_tmp = ARMSmdl(y0, mdl$ar, k)
    # test stat
    LRT_N[xxi] = -2*(mdl_h0_tmp$logLike - mdl_h1_tmp$logLike)
  }
  return(sort(LRT_N))
}
# -----------------------------------------------------------------------------
#' @title LRMC simultion based on coefficients
#' 
#' This function ...
#'
#' @param infile Path to the input file
#' @return A matrix of the infile
#' @export
LRMC_fn_simu = function(w,mdl,k,N){
  # Uses "iter" (i.e. iteration method) with max iterations set to 100 for speed
  LRT_N = matrix(0,N,1)
  for (xxi in 1:N){
    y0 = simu_AR_dgp(n = mdl$n, mu = mean(mdl$y),std = mdl$stdev,w)
    mdl_h0_tmp = ARmdl(y0, mdl$ar)
    # ***NOTE: Estimation of alternative model, when performing MMC test only, is
    #          performed using 100 iterations of EM algo only, due to speed 
    #          constraint. Edit below. 
    mdl_h1_tmp = ARMSmdl(y0, mdl$ar, k, method = "iter", mxit = 100)
    # test stat
    LRT_N[xxi] = -2*(mdl_h0_tmp$logLike - mdl_h1_tmp$logLike)
  }
  return(sort(LRT_N))
}
# -----------------------------------------------------------------------------
#' @title LRMC ARMS maximization function 
#' This function ...
#'
#' @param infile Path to the input file
#' @return A matrix of the infile
#' @export
LRMMC_ARMS_fn_max = function(w){
  if (min(Mod(as.complex(polyroot(rbind(1,as.matrix(w))))))<1){
    pval =  0
  }else{
    mdl_tmp = get('mdl_h0',envir = LRMMCenv)
    LRT_0 = get('LRT_0',envir = LRMMCenv)
    k = get('k',envir = LRMMCenv)
    N = get('N',envir = LRMMCenv)
    LRT_N = LRMC_fn_simu(w,mdl_tmp,k,N-1) 
    pval = p_val(LRT_0, LRT_N, type = 'geq')
  }
  return(pval)
}
# -----------------------------------------------------------------------------
#' @title LRMC ARMMC minimization function
#' This function ...
#'
#' @param infile Path to the input file
#' @return A matrix of the infile
#' @export
LRMMC_ARMS_fn_min = function(w){
  n_pval = -LRMMC_ARMS_fn_max(w)
  return(n_pval)
}
# -----------------------------------------------------------------------------
#' @title LRMC ARMS
#' This function ...
#'
#' @param infile Path to the input file
#' @return A matrix of the infile
#' @export
LRMC_ARMS = function(Y, ar = NULL, k = 2, N = 100){
  # --------------- Likelihood-Ratio Monte-Carlo test for Autoregressive Markov-Switching 
  mdl_h0 = ARmdl(Y,ar)
  mdl_h1 = ARMSmdl(Y,ar,k)
  # null model (linear model)
  l_0 = mdl_h0$logLike
  # alternative model (MS model)
  l_1 = mdl_h1$logLike
  # calc test stat
  LRT_0 = -2*(l_0-l_1)
  # empirical Monte-Carlo simulations using null model 
  LRT_N = LRMC_simu(mdl_h0,k,N-1)
  # calculate p-value
  pval = p_val(LRT_0, LRT_N, type = 'geq')
  # format function output 
  output = list()
  output$pval = pval
  output$mdl_h0 = mdl_h0
  output$msmdl_h1 = mdl_h1
  output$LRT_0 = LRT_0
  output$LRT_N = LRT_N
  return(output)
}
# -----------------------------------------------------------------------------
#' @title LR MMC ARMS
#' This function ...
#'
#' @param infile Path to the input file
#' @return A matrix of the infile
#' @export
LRMMC_ARMS = function(Y,ar = NULL, k = 2,N = 100, method = "GenSA", controls=NULL){
  # ---------------- Likelihood-Ratio MMC test for Autoregressive Markov-Switching 
  mdl_h0 = ARmdl(Y,ar)
  mdl_h1 = ARMSmdl(Y,ar,k)
  ## ----- Calc Lik of H0 ------ ##
  l_0 = mdl_h0$logLike
  ## ----- Calc Lik of H1 ------ ##
  l_1 = mdl_h1$logLike
  # Calc Test Stat
  LRT_0 = -2*(l_0-l_1)
  # Some initial values 
  phi = mdl_h0$coef[2:length(mdl_h0$coef)]
  phiSE = mdl_h0$se[2:length(mdl_h0$coef)]
  # Assign variables to be used in optimization to global environment 
  assign('mdl_h0',mdl_h0,envir = LRMMCenv)
  assign('LRT_0',LRT_0,envir = LRMMCenv)
  assign('k',k,envir = LRMMCenv)
  assign('N',N,envir = LRMMCenv)
  # Use MMC Optimization function
  output = LRMMC_opt(phi,phiSE,method,opt_func=LRMMC_ARMS_fn_min,controls)
  # format remaining functions output variables
  output$mdl_h0 = mdl_h0
  output$msmdl_h1 = mdl_h1
  output$LR_stat = LRT_0
  return(output)
}
# -----------------------------------------------------------------------------
#' @title LRMMC optimization
#' This function ...
#'
#' @param infile Path to the input file
#' @return A matrix of the infile
#' @export
LRMMC_opt = function(w,w_se=NULL,method="GenSA",opt_func=NULL,controls=NULL){
  # Add GA, PSO and
  # ---------------- MMC Optimization procedure 
  # Optimization Controls 
  if (is.null(controls)==T){
    #maxcont = list(max.time=30)
    maxcont = list(maxit=5)
    if (is.null(w_se)){
      # entire space
      low = rep(-0.99,length(w))
      upp = rep(0.99,length(w))
    }else{
      # Consistent Estimator space
      low = round(w - 1.96*w_se,2)
      upp = round(w + 1.96*w_se,2)
    }
  }else{
    maxcont = list(max.time=controls$max.time) 
    low = controls$lower
    upp = controls$upper
  }
  # Optimization
  output = list()
  if (method=="GenSA"){
    resGSA  = GenSA::GenSA(par = w,fn = opt_func,lower = low,upper = upp,control = maxcont)
    # save output
    output[['pval']] = -resGSA$value
    output[['trace']]  = resGSA$trace
    output[['res']] = resGSA
  }else if(method=="GA"){
    
    
  }else if(method=="pso"){
    
    
  }else if(method=="randSearch"){
    N2 = 100
    if (is.null(w_se)==T){
      randPhi = rnorm(N2,w,0.5)
    }else{
      randPhi = rnorm(N2,w,w_se) 
    }
    mmcSave = matrix(0,N2,1)
    for (xi in 1:N2){
      mmcSave[xi] = LRMMC_ARMS_fn_max(randPhi[xi])
    }
    output[['pval']] = max(mmcSave)
    output[['trace']] = cbind(randPhi,mmcSave)
  }
  return(output)
}

rowsum.finite = function(mat){
  vec = matrix(0,nrow(mat),1)
  for (xi in 1:nrow(mat)){
    vec[xi,1] = sum(mat[xi,is.na(mat[xi,])==0])
  }
  return(vec)
}


calcPower_mmc = function(mu,sig, p11 = 0.9, p22 = 0.9, phi = 0.9 , sampsize = 100, N = 99, srep = 100, parcontrols = NULL,Dir=NULL){
  nar=1
  k=2
  th = expand.grid(p22,p11,sampsize,phi,mu,sig)
  th = th[rowSums(th[,(5:6)])>0,]
  th_par = cbind(0 , 0 + th[,5], th[,2], th[,1], 1, 1 + th[,6], th[,4])
  th_n = th[,3]
  alpha05 = 0.05
  alpha10 = 0.10
  allout1 = list()
  allpval1 = matrix(NA,nrow=nrow(th),ncol=srep)
  for (xrep in 1:srep){
    simuLst  = list()
    for (xth in 1:nrow(th)){
      simuLst[[xth]] = as.matrix(simu_ARMS_dgp(th_n[xth],th_par[xth,]))
    }
    st = Sys.time()
    cl <- parallel::makeCluster(parcontrols$cores)
    res = list()
    doParallel::registerDoParallel(cl)
    res = foreach(xi = 1:nrow(th),.export=c("simu_AR_dgp","simu_ARMS_dgp","ARmdl","ARMSmdl", "LRMC_fn_simu", "LRMMC_ARMS_fn_max", "LRMMC_ARMS_fn_min", "LRMMC_ARMS", "LRMMC_opt"),.packages = "Rcpp",.noexport = c("transMat_cpp", "randInitVal_cpp","ms_logLikel_cpp","EM_hamilton_cpp","msmdl_cpp")) %dopar% {
      sourceCpp(paste0(Dir,"LRMMC_ARMSmdl/code/LRMMC_cppfuncs.cpp"))
      mdl_h0 = ARmdl(simuLst[[xi]],nar)
      mdl_h1 = ARMSmdl(simuLst[[xi]],nar,k)
      out = LRMMC_ARMS(mdl_h0,mdl_h1,k,N,method="randSearch")
    }
    parallel::stopCluster(cl)
    rtime = Sys.time() - st
    print(paste0("rep: ",xrep," - Time elpased is: ",rtime))
    allout1[[xrep]] = res
    allpval1[,xrep] = as.matrix(unlist(lapply(res, `[[`, 1)))
  }
  allPower_ms_05 = as.matrix((rowsum.finite(allpval1<=alpha05)/srep)*100)
  allPower_ms_10 = as.matrix((rowsum.finite(allpval1<=alpha10)/srep)*100)
  tab_ms_05 = matrix(NA,nrow=length(p11)*length(p22),ncol=length(phi)*length(sampsize))
  tab_ms_10 = matrix(NA,nrow=length(p11)*length(p22),ncol=length(phi)*length(sampsize))
  xcol = 1
  for (xphi in 1:length(phi)){
    loc = seq(1:nrow(th))
    loc_tmp = loc[th[,4]==phi[xphi]]
    for (xsampsize in 1:length(sampsize)){
      loc_tmp2 =  loc_tmp[th[loc_tmp,3]==sampsize[xsampsize]]
      tab_ms_05[,xcol] = allPower_ms_05[loc_tmp2,]
      tab_ms_10[,xcol] = allPower_ms_10[loc_tmp2,]
      xcol=xcol+1
    }
  }
  output = list()
  output[['power_level05_ms']] = tab_ms_05
  output[['power_level10_ms']] = tab_ms_10
  output[['trace_ms']] = allout1
  return(output)
}




calcSize_mmc = function(sampsize=100, N=99, srep=100, phi=0.9, parcontrols=NULL, Dir=NULL){
  nar=1
  k=2
  th = expand.grid(sampsize,phi)
  alpha05 = 0.05
  alpha10 = 0.10
  allout1 = list()
  allpval1 = matrix(NA,nrow=nrow(th),ncol=srep)
  for (xrep in 1:srep){
    simuLst  = list()
    for (xth in 1:nrow(th)){
      simuLst[[xth]] = as.matrix(simu_AR_dgp(th[xth,1],mu=0,std=1,th[xth,2]))
    }
    st = Sys.time()
    cl <- parallel::makeCluster(parcontrols$cores)
    doParallel::registerDoParallel(cl)
    res = foreach(xi = 1:nrow(th),.export=c("simu_AR_dgp","simu_ARMS_dgp","ARmdl","ARMSmdl", "LRMC_fn_simu", "LRMMC_ARMS_fn_max", "LRMMC_ARMS_fn_min", "LRMMC_ARMS", "LRMMC_opt"),.packages = "Rcpp",.noexport = c("transMat_cpp", "randInitVal_cpp","ms_logLikel_cpp","EM_hamilton_cpp","msmdl_cpp")) %dopar% {
      sourceCpp(paste0(Dir,"LRMMC_ARMSmdl/code/LRMMC_cppfuncs.cpp"))
      mdl_h0 = ARmdl(simuLst[[xi]],nar)
      mdl_h1 = ARMSmdl(simuLst[[xi]],nar,k)
      out = LRMMC_ARMS(mdl_h0,mdl_h1,k,N,method="randSearch")
    }
    
    parallel::stopCluster(cl)
    rtime = Sys.time() - st
    print(paste0("rep: ",xrep," - Time elpased is: ",rtime))
    allout1[[xrep]] = res
    allpval1[,xrep] = as.matrix(unlist(lapply(res, `[[`, 1)))
  }
  allPower_ms_05 = as.matrix((rowsum.finite(allpval1<=alpha05)/srep)*100)
  allPower_ms_10 = as.matrix((rowsum.finite(allpval1<=alpha10)/srep)*100)
  tab_ms_05 = matrix(NA,nrow=1,ncol=length(phi)*length(sampsize))
  tab_ms_10 = matrix(NA,nrow=1,ncol=length(phi)*length(sampsize))
  xcol = 1
  for (xphi in 1:length(phi)){
    loc = seq(1:nrow(th))
    loc_tmp = loc[th[,2]==phi[xphi]]
    for (xsampsize in 1:length(sampsize)){
      loc_tmp2 =  loc_tmp[th[loc_tmp,1]==sampsize[xsampsize]]
      tab_ms_05[,xcol] = allPower_ms_05[loc_tmp2,]
      tab_ms_10[,xcol] = allPower_ms_10[loc_tmp2,]
      xcol=xcol+1
    }
  }
  output = list()
  output[['power_level05_ms']] = tab_ms_05
  output[['power_level10_ms']] = tab_ms_10
  output[['trace_ms']] = allout1
  return(output)
}






calcPower_mc = function(mu,sig, p11 = 0.9, p22 = 0.9, phi = 0.9 ,sampsize = 100, N = 99, srep = 100, parcontrols = NULL,Dir=NULL){
  nar=1
  k=2
  th = expand.grid(p22,p11,sampsize,phi,mu,sig)
  th = th[rowSums(th[,(5:6)])>0,]
  th_par = cbind(0 , 0 + th[,5], th[,2], th[,1], 1, 1 + th[,6], th[,4])
  th_n = th[,3]
  alpha05 = 0.05
  alpha10 = 0.10
  allout1 = list()
  allpval1 = matrix(NA,nrow=nrow(th),ncol=srep)
  for (xrep in 1:srep){
    simuLst  = list()
    for (xth in 1:nrow(th)){
      simuLst[[xth]] = as.matrix(simu_ARMS_dgp(th_n[xth],th_par[xth,]))
    }
    st = Sys.time()
    cl <- parallel::makeCluster(parcontrols$cores)
    res = list()
    doParallel::registerDoParallel(cl)
    res = foreach(xi = 1:nrow(th),.export=c("simu_AR_dgp","simu_ARMS_dgp","ARmdl","ARMSmdl", "LRMC_simu", "LRMC_ARMS"),.packages = "Rcpp",.noexport = c("transMat_cpp", "randInitVal_cpp","ms_logLikel_cpp","EM_hamilton_cpp","msmdl_cpp")) %dopar% {
      sourceCpp(paste0(Dir,"LRMMC_ARMSmdl/code/LRMMC_cppfuncs.cpp"))
      mdl_h0 = ARmdl(simuLst[[xi]],nar)
      mdl_h1 = ARMSmdl(simuLst[[xi]],nar,k)
      out = LRMC_ARMS(mdl_h0,mdl_h1,k,N)
    }
    rtime = Sys.time() - st
    parallel::stopCluster(cl)
    print(paste0("rep: ",xrep," - Time elpased is: ",rtime))
    allout1[[xrep]] = res
    allpval1[,xrep] = as.matrix(unlist(lapply(res, `[[`, 1)))
  }
  allPower_ms_05 = as.matrix((rowsum.finite(allpval1<=alpha05)/srep)*100)
  allPower_ms_10 = as.matrix((rowsum.finite(allpval1<=alpha10)/srep)*100)
  tab_ms_05 = matrix(NA,nrow=length(p11)*length(p22),ncol=length(phi)*length(sampsize))
  tab_ms_10 = matrix(NA,nrow=length(p11)*length(p22),ncol=length(phi)*length(sampsize))
  xcol = 1
  for (xphi in 1:length(phi)){
    loc = seq(1:nrow(th))
    loc_tmp = loc[th[,4]==phi[xphi]]
    for (xsampsize in 1:length(sampsize)){
      loc_tmp2 =  loc_tmp[th[loc_tmp,3]==sampsize[xsampsize]]
      tab_ms_05[,xcol] = allPower_ms_05[loc_tmp2,]
      tab_ms_10[,xcol] = allPower_ms_10[loc_tmp2,]
      xcol=xcol+1
    }
  }
  output = list()
  output[['power_level05_ms']] = tab_ms_05
  output[['power_level10_ms']] = tab_ms_10
  output[['trace_ms']] = allout1
  return(output)
}




calcSize_mc = function(sampsize=100, N=99, srep=100, phi=0.9, parcontrols=NULL, Dir=NULL){
  nar=1
  k=2
  th = expand.grid(sampsize,phi)
  alpha05 = 0.05
  alpha10 = 0.10
  allout1 = list()
  allpval1 = matrix(NA,nrow=nrow(th),ncol=srep)
  for (xrep in 1:srep){
    simuLst  = list()
    for (xth in 1:nrow(th)){
      simuLst[[xth]] = as.matrix(simu_AR_dgp(th[xth,1],mu=0,std=1,th[xth,2]))
    }
    st = Sys.time()
    cl <- parallel::makeCluster(parcontrols$cores)
    doParallel::registerDoParallel(cl)
    res = foreach(xi = 1:nrow(th),.export=c("simu_AR_dgp","simu_ARMS_dgp","ARmdl","ARMSmdl", "LRMC_simu", "LRMC_ARMS"),.packages = "Rcpp",.noexport = c("transMat_cpp", "randInitVal_cpp","ms_logLikel_cpp","EM_hamilton_cpp","msmdl_cpp")) %dopar% {
      sourceCpp(paste0(Dir,"LRMMC_ARMSmdl/code/LRMMC_cppfuncs.cpp"))
      mdl_h0 = ARmdl(simuLst[[xi]],nar)
      mdl_h1 = ARMSmdl(simuLst[[xi]],nar,k)
      out = LRMC_ARMS(mdl_h0,mdl_h1,k,N)
    }
    parallel::stopCluster(cl)
    rtime = Sys.time() - st
    print(paste0("rep: ",xrep," - Time elpased is: ",rtime))
    allout1[[xrep]] = res
    allpval1[,xrep] = as.matrix(unlist(lapply(res, `[[`, 1)))
  }
  allPower_ms_05 = as.matrix((rowsum.finite(allpval1<=alpha05)/srep)*100)
  allPower_ms_10 = as.matrix((rowsum.finite(allpval1<=alpha10)/srep)*100)
  tab_ms_05 = matrix(NA,nrow=1,ncol=length(phi)*length(sampsize))
  tab_ms_10 = matrix(NA,nrow=1,ncol=length(phi)*length(sampsize))
  xcol = 1
  for (xphi in 1:length(phi)){
    loc = seq(1:nrow(th))
    loc_tmp = loc[th[,2]==phi[xphi]]
    for (xsampsize in 1:length(sampsize)){
      loc_tmp2 =  loc_tmp[th[loc_tmp,1]==sampsize[xsampsize]]
      tab_ms_05[,xcol] = allPower_ms_05[loc_tmp2,]
      tab_ms_10[,xcol] = allPower_ms_10[loc_tmp2,]
      xcol=xcol+1
    }
  }
  output = list()
  output[['power_level05_ms']] = tab_ms_05
  output[['power_level10_ms']] = tab_ms_10
  output[['trace_ms']] = allout1
  return(output)
}
