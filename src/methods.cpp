#include <RcppArmadillo.h>
#include <RcppArmadilloExtensions/sample.h>
#include "models.h"
//[[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;

// ==============================================================================
//' @title Covariance to correlation matrix
//' 
//' @description This function takes an (n x n) covariance matrix and returns the associated (n x n) correlation matrix
//' 
//' @param (n x n) covariance matrix
//' 
//' @return (n x n) correlation matrix
//' 
//' @export
// [[Rcpp::export]]
arma::mat cov2corr(arma::mat cov_mat){
  arma::mat corr_mat = inv(diagmat(sqrt(cov_mat)))*cov_mat*inv(diagmat(sqrt(cov_mat)));
  return(corr_mat);
}

// ==============================================================================
//' @title Covariance vech function
//' 
//' @description This function returns the half-vectorization of an input matrix as a column vector.
//' 
//' @param mat: (n x n) covariance matrix
//' 
//' @return (n+1)*n/2 column vector 
//' 
//' @export
// [[Rcpp::export]]
arma::vec covar_vech(arma::mat mat){
  int nr = mat.n_rows;
  int nc = mat.n_cols;
  if (nc==nr){
    arma::vec sigma_vec = trans(mat.row(0));
    for (int xq = 1; xq<nr; xq++){
      sigma_vec = join_vert(sigma_vec, trans(mat.submat(xq,xq,xq,nr-1)));
    }
    return(sigma_vec);  
  }else{
    stop("Input must be a square matrix");
  }
}

// ==============================================================================
//' @title Covariance vech to matrix 
//' 
//' @description This function undoes the half-vectorization of a covariance matrix.
//' 
//' @param sig: (n+1)*n/2 vector 
//' n: integer determining shape of the orginal matrix
//' 
//' @return (n x n) covariance matrix
//' 
//' @export
// [[Rcpp::export]]
arma::mat covar_unvech(arma::vec sig, int n){
  arma::mat sigma_mat(n, n, arma::fill::zeros);
  int count = 0;
  for (int xq = 0; xq<n; xq++){
    for (int xqq = xq; xqq<n; xqq++){
      sigma_mat(xq,xqq) = sig(count);
      sigma_mat(xqq,xq) = sig(count);
      count += 1;
    }
  }
  return(sigma_mat);
}


// ==============================================================================
//' @title Random Transition Matrix
//' 
//' @description This function creates a (k x k) random transition matrix
//' 
//' @param k number of regimes. Must be greater than or equal to 2. 
//' @param n number of random sample to use. By default it is 100 but this can be set to length of TS for example
//'  
//' 
//' @return transition matrix with randomly generated entries.
//' 
//' @export
// [[Rcpp::export]]
arma::mat randTransMat(int k){
  arma::mat P = reshape(arma::randu(k*k), k, k); 
  arma::vec PcolSums = trans(arma::sum(P,0));
  for (int xk = 0; xk<k; xk++){
    P.col(xk) = P.col(xk)/PcolSums(xk);
  } 
  return(P);
}

// ==============================================================================
//' @title Ergodic (limiting) Probabilities of States
//' 
//' @description Takes a transition matrix and returns the limiting probabilities
//' 
//' @param P matrix with transition probabilities
//' 
//' @return Vector of Limiting probabilities of a transition matrix
//' 
//' @export
// [[Rcpp::export]]
arma::vec limP(arma::mat P, int k){
  arma::mat onevec(1, k, arma::fill::ones);
  arma::mat ep(1, k+1, arma::fill::zeros);
  ep(0,k) = 1;
  arma::mat Atmp = join_cols(arma::eye(k,k)-P,onevec);
  //arma::vec pinf = inv(trans(Atmp)*Atmp)*trans(Atmp)*trans(ep);
  arma::vec pinf = solve(trans(Atmp)*Atmp,trans(Atmp))*trans(ep);
  return (pinf);
}

// ==============================================================================
//' @title Lagged Time Series Data
//' 
//' @description This function takes a (Tx1) vector Y and returns the (T-px1) vector y and matrix of lagged observations.
//' 
//' @param Y vector with time series observations. Required argument. 
//' @param ar integer for the number of lags to use in estimation. Must be greater than or equal to 1. Default is 1.
//' 
//' @return List with vector y (vector of lagged Y) and matrix X of lagged observations.
//' 
//' @export
// [[Rcpp::export]]
List ts_lagged(arma::mat Y, int ar){
  int Tsize = Y.n_rows;
  int N = Y.n_cols;
  arma::mat y = Y.submat(ar,0,Tsize-1,N-1);
  int n = Tsize - ar;
  arma::mat X(n, ar*N, arma::fill::zeros);
  for (int xp = 0; xp<ar; xp++){
    X.submat(0,N*xp,n-1,N*xp+N-1) = Y.submat(ar-(xp+1),0,Tsize-(xp+1)-1,N-1);
  }
  List lagged_output;
  lagged_output["y"] = y;
  lagged_output["X"] = X;
  return(lagged_output);
}

// ==============================================================================
//' @title Parameter List
//' 
//' @description This function takes the vector of parameters of interest and converts it to a list with the parameters seperated
//' 
//' @param theta vector of parameters.
//' @param mdl List of model properties
//' @param k number of regimes. Must be greater than or equal to 2. 
//' @param msmu bool indicating is the mean switches with regime 
//' @param msvar bool indicating is the variance switches with regime 
//' 
//' @return List with the mean, variance, transition matrix and limiting probabilities
//' 
//' @export
// [[Rcpp::export]]
List paramListMS(arma::vec theta, int ar, int k, bool msmu, bool msvar){
  Rcpp::Environment mstest("package:MSTest");
  Rcpp::Function musigGrid = mstest["musigGrid"];
  Rcpp::Function transMatAR = mstest["transMatAR"];
  // ----- Mean for each regime 
  arma::vec mu = theta.subvec(0, msmu*(k-1));
  // ----- Variance for each regime 
  arma::vec sig = theta.subvec(1+msmu*(k-1),1+msmu*(k-1)+msvar*(k-1));
  // ----- Phi vector
  arma::vec phi(ar, arma::fill::zeros);
  phi =  theta.subvec(2+msmu*(k-1)+msvar*(k-1), 2+msmu*(k-1)+msvar*(k-1)+ar-1);
  // ----- Transition probabilities 
  arma::mat P = reshape(theta.subvec(2+msmu*(k-1)+msvar*(k-1) + ar, 2+msmu*(k-1)+msvar*(k-1) + ar + k*k - 1),k, k);
  // Regime limiting probabilities
  arma::vec pinf = limP(P, k);
  // ----- Obtain AR consistent grid of Mu, Sigma and State indicators
  List musig_out = musigGrid(mu, sig, k, ar, msmu, msvar);
  arma::mat muAR = as<arma::mat>(musig_out["mu"]);
  arma::mat sigAR = as<arma::mat>(musig_out["sig"]);
  arma::vec state_ind = as<arma::vec>(musig_out["state_ind"]);
  // ----- Obtain AR consistent P and pinf
  int M = pow(k, ar+1);
  arma::mat P_AR = as<arma::mat>(transMatAR(P, k, ar));
  arma::mat pinf_AR = limP(P_AR, M);
  // ----- Organize output
  List param_out;
  param_out["mu"] = mu;
  param_out["sig"] = sig;
  param_out["phi"] = phi;
  param_out["P"] = P;
  param_out["pinf"] = pinf;
  param_out["muAR"] = muAR;
  param_out["sigAR"] = sigAR;
  param_out["state_ind"] = state_ind;
  param_out["P_AR"] = P_AR;
  param_out["pinf_AR"] = pinf_AR;
  return(param_out);
}
// ==============================================================================
//' @title Parameter List
//' 
//' @description This function takes the vector of parameters of interest and converts it to a list with the parameters seperated
//' 
//' @param theta vector of parameters.
//' @param mdl List of model properties
//' @param k number of regimes. Must be greater than or equal to 2. 
//' @param msmu bool indicating is the mean switches with regime 
//' @param msvar bool indicating is the variance switches with regime 
//' 
//' @return List with the mean, variance, transition matrix and limiting probabilities
//' 
//' @export
// [[Rcpp::export]]
List paramListMSVAR(arma::vec theta, int q, int ar, int k, bool msmu, bool msvar){
  Rcpp::Environment mstest("package:MSTest");
  Rcpp::Function musigVARGrid = mstest["musigVARGrid"];
  Rcpp::Function transMatAR = mstest["transMatAR"];
  // ----- Mean for each regime 
  arma::mat mu_k(k, q, arma::fill::zeros);
  arma::vec mu = theta.subvec(0, q+q*msmu*(k-1)-1);
  if (msmu==TRUE){
    for (int xk = 0; xk<k; xk++){
      mu_k.row(xk) = trans(mu.subvec(xk*q,xk*q+q-1));
    }
  }else{
    for (int xk = 0; xk<k; xk++){
      mu_k.row(xk) = trans(mu);
    }
  }
  // ----- Variance for each regime 
  int sigN = (q*(q+1))/2;
  arma::vec sig = theta.subvec(q+q*msmu*(k-1),q+q*msmu*(k-1)+sigN+sigN*msvar*(k-1)-1);
  List sigma(k);
  if (msvar==TRUE){
    for (int xk = 0; xk<k; xk++){
      arma::vec sig_tmp = sig.subvec(sigN*xk,sigN*xk+sigN-1);
      sigma[xk] = covar_unvech(sig_tmp, q);
    } 
  }else{
    for (int xk = 0; xk<k; xk++){
      sigma[xk] = covar_unvech(sig, q);
    }
  }
  // ----- Phi vector
  int phiN = q+q*msmu*(k-1)+sigN+sigN*msvar*(k-1);
  arma::vec phi_tmp =  theta.subvec(phiN, phiN+ q*q*ar-1);
  arma::mat phi = reshape(phi_tmp, q*ar, q);
  // ----- Transition probabilities 
  int PN = q+q*msmu*(k-1)+sigN+sigN*msvar*(k-1)+q*q*ar;
  arma::mat P = reshape(theta.subvec(PN, PN + k*k - 1), k, k);
  // Regime limiting probabilities
  arma::vec pinf = limP(P, k);
  // ----- Obtain AR consistent grid of Mu, Sigma and State indicators
  List musig_out = musigVARGrid(mu_k, sigma, k, ar, msmu, msvar);
  List muAR = musig_out["mu"];
  List sigAR = musig_out["sig"];
  arma::vec state_ind = as<arma::vec>(musig_out["state_ind"]);
  // ----- Obtain AR consistent P and pinf
  int M = pow(k, ar+1);
  arma::mat P_AR = as<arma::mat>(transMatAR(P, k, ar));
  arma::mat pinf_AR = limP(P_AR, M);
  // ----- Organize output
  List param_out;
  param_out["mu"] = mu_k;
  param_out["sigma"] = sigma;
  param_out["phi"] = phi;
  param_out["P"] = P;
  param_out["pinf"] = pinf;
  param_out["muAR"] = muAR;
  param_out["sigAR"] = sigAR;
  param_out["state_ind"] = state_ind;
  param_out["P_AR"] = P_AR;
  param_out["pinf_AR"] = pinf_AR;
  return(param_out);
}
// ==============================================================================
//' @title Calculate Residuals for Markov-switching model
//' 
//' @description This function computes residuals when mean or conditional mean 
//' (i.e. when dealing with AR model) switches with regime.
//' 
//' @param mdl List containing relevant parameters.
//' @param mu vector with mean in each regime.
//' @param k number of regimes. Must be greater than or equal to 2. 
//' @param ar bool indicator for whether model is AR model.
//'  
//' 
//' @return (Txk) matrix of residuals in each regime.
//' 
//' @export
// [[Rcpp::export]]
arma::mat calcMSResid(List mdl, arma::mat mu, int k){
  // ---------- Initialize parameters
  arma::vec y = mdl["y"];
  int ar = mdl["ar"];
  int Tsize = y.n_elem;
  int M = pow(k, ar+1); // number of regimes consistent with autoregressive structure (M=k if ar=0 and M=k^(ar+1) if ar>0)
  arma::mat repvec(1, M, arma::fill::ones);
  arma::mat repmu(Tsize, 1, arma::fill::ones);
  arma::mat ms_y = y*repvec;
  arma::mat resid(Tsize, M, arma::fill::zeros);
  // ---------- Compute Residuals 
  arma::mat z = ms_y - repmu*trans(mu.col(0)); // [y(t) - mu_s(t))]
  arma::mat x = mdl["x"];
  arma::vec phi = mdl["phi"];
  arma::mat xz(Tsize, M, arma::fill::zeros);
  for (int xkp = 0; xkp<M; xkp++){
    arma::mat zx_tmp = x - repmu*mu.submat(xkp,1,xkp,ar); // [y(t-i) - mu_s(t-i))]
    xz.col(xkp) = zx_tmp*phi;
  }
  resid = z - xz;
  return(resid);
}
// ==============================================================================
//' @title Calculate Residuals for Markov-switching VAR model
//' 
//' @description This function computes residuals when mean or conditional mean 
//' (i.e. when dealing with AR model) switches with regime.
//' 
//' @param mdl List containing relevant parameters.
//' @param mu vector with mean in each regime.
//' @param k number of regimes. Must be greater than or equal to 2. 
//' @param ar bool indicator for whether model is AR model.
//'  
//' 
//' @return (Txk) matrix of residuals in each regime.
//' 
//' @export
// [[Rcpp::export]]
List calcMSVARResid(List mdl, List mu, int k){
  // ---------- Initialize parameters
  arma::mat y = mdl["y"];
  int ar = mdl["ar"];
  arma::mat x = mdl["x"];
  arma::mat phi = mdl["phi"];
  int Tsize = y.n_rows;
  int N = y.n_cols;
  int M = pow(k, ar+1); // number of regimes consistent with autoregressive structure (M=k if ar=0 and M=k^(ar+1) if ar>0)
  arma::mat repmu(Tsize, 1, arma::fill::ones);
  List z_mat(M); // [y(t) - mu_s(t))]
  List xz_mat(M); // [y(t-p) - mu_s(t-p))]
  List resid(M); 
  for (int xm = 0; xm<M; xm++){
    arma::mat mu_tmp = mu[xm]; 
    arma::mat y_tmp = y - repmu*trans(mu_tmp.col(0));
    arma::mat xz_tmp(Tsize, N*ar, arma::fill::zeros); 
    for (int xp = 0; xp<ar; xp++){
      xz_tmp.submat(0,N*xp,Tsize-1,N*xp+N-1) = x.submat(0,N*xp,Tsize-1,N*xp+N-1) - repmu*trans(mu_tmp.col(xp+1));
    }
    resid[xm] = y_tmp - xz_tmp*phi;
    z_mat[xm] = y_tmp;
    xz_mat[xm] = xz_tmp;
  }
  return(resid);
}

// ==============================================================================
//' @title generate initial values for MS model
//' 
//' 
//' @export
// [[Rcpp::export]]
arma::vec initValsMS(List mdl, int k){
  arma::vec phi = mdl["phi"];
  double mu = mdl["mu"];
  double stdev = mdl["stdev"];
  bool msmu = mdl["msmu"];
  bool msvar = mdl["msvar"];
  // pre-allocate mu and sigma vectors
  arma::vec mu_0(1+msmu*(k-1), arma::fill::zeros);
  arma::vec sig_0(1+msvar*(k-1), arma::fill::zeros);
  // Set initial values using linear model if no switch
  mu_0(0) = mu;
  sig_0(0) = pow(stdev,2);
  // initial values for mu around linear model mu when switch
  if (msmu==TRUE){
    mu_0 = mu + (3*stdev)*arma::randn<arma::vec>(k);
  }
  // initial values for stdev around linear model stdev when switch
  if (msvar==TRUE){
    sig_0 = (stdev*0.1) + ((2*stdev)-(stdev*0.1))*arma::randu<arma::vec>(k);
  }
  // create vector for initial values
  arma::vec theta_0 = join_vert(mu_0, sig_0);
  theta_0 = join_vert(theta_0, phi);
  arma::mat P_0 = randTransMat(k);
  theta_0 = join_vert(theta_0, vectorise(P_0));
  return(theta_0);
}

// ==============================================================================
//' @title generate initial values for MS-VAR model
//' 
//' 
//' @export
// [[Rcpp::export]]
arma::vec initValsMSVAR(List mdl, int k){
  arma::vec phi = mdl["phi"];
  arma::vec mu = mdl["mu"];
  arma::mat sigma = mdl["sigma"];
  bool msmu = mdl["msmu"];
  bool msvar = mdl["msvar"];
  int q = mu.n_elem;
  int sigN = (q*(q+1))/2;
  // pre-allocate mu and sigma matrices
  arma::mat mu_0(1+msmu*(k-1), q, arma::fill::zeros);
  arma::mat sigma_0(1+msvar*(k-1), sigN, arma::fill::zeros);
  // Set initial values using linear model if no switch
  mu_0.row(0) = trans(mu);
  sigma_0.row(0) = trans(covar_vech(sigma));
  // initial values for mu around linear model mu when switch
  if (msmu==TRUE){
    arma::mat repvec(k,1,arma::fill::ones);
    mu_0 =   repvec*trans(mu) + (repvec*trans(3*sqrt(sigma.diag())))%arma::mat(k,q,arma::fill::randn);
  }
  arma::vec mu_out = vectorise(trans(mu_0));
  // initial values for stdev around linear model stdev when switch
  if (msvar==TRUE){
    arma::vec sigma_vec_tmp = covar_vech(sigma);
    arma::mat sig_mat_tmp = covar_vech((sigma_vec_tmp*0.1) + ((2*sigma_vec_tmp)-(sigma_vec_tmp*0.1))%arma::randu<arma::vec>(sigN));
    sig_mat_tmp = sig_mat_tmp*trans(sig_mat_tmp);
    sig_mat_tmp = sig_mat_tmp + q*arma::speye(q,q);
    sigma_0.row(0) = trans(covar_vech(sig_mat_tmp));
    for (int xk = 1; xk<k; xk++){
      arma::mat sig_mat_tmp = covar_vech((sigma_vec_tmp*0.1) + ((2*sigma_vec_tmp)-(sigma_vec_tmp*0.1))%arma::randu<arma::vec>(sigN));
      sig_mat_tmp = sig_mat_tmp*trans(sig_mat_tmp);
      sig_mat_tmp = sig_mat_tmp + q*arma::speye(q,q);
      sigma_0.row(xk) = trans(covar_vech(sig_mat_tmp));
    }
  }
  arma::vec sigma_out = vectorise(trans(sigma_0));
  // create vector for initial values
  arma::vec theta_0 = join_vert(mu_out, sigma_out);
  theta_0 = join_vert(theta_0, vectorise(phi));
  arma::mat P_0 = randTransMat(k);
  theta_0 = join_vert(theta_0, vectorise(P_0));
  return(theta_0);
}


// ==============================================================================
//' @title Calculate MC P-Value
//' 
//' @description This function calculates a Monte-Carlo p-value
//' 
//' @param test_stat test statistic under the alternative (e.g. S_0)
//' @param null_vec series with test statistic under the null (i.e. vector S)
//' @param type like of test. options are: "geq" for right-tail test, "leq" for 
//' left-tail test, "abs" for absolute vallue test and "two-tail" for two-tail test.
//' 
//' @return MC p-value of test
//' 
//' @references Dufour, J. M. (2006). Monte Carlo tests with nuisance parameters: 
//' A general approach to finite-sample inference and nonstandard asymptotics. 
//' Journal of Econometrics, 133(2), 443-477.
//' @references Dufour, J. M., & Luger, R. (2017). Identification-robust moment-based 
//' tests for Markov switching in autoregressive models. Econometric Reviews, 36(6-9), 713-727.
//' 
//' @export
// [[Rcpp::export]]
double MCpval(double test_stat, arma::vec null_vec, Rcpp::String type = "geq"){
  int N = null_vec.n_elem;
  arma::vec u   = arma::randu(1 + sum(null_vec == test_stat)); 
  double test_rank  = sum(test_stat > null_vec) + sum(u <= u(0));
  // negative p-value can occur two-tailed test
  double survival_pval = ((N + 1 - test_rank)/ N);
  // Compute the p-value
  double pval = 999; // This will be returned if allowed type isnt specified. 
  if (type == "absolute" || type == "geq") {
    pval = (N * survival_pval + 1)/(N + 1);
  }else if (type == "leq") {
    pval = (N * (1 - survival_pval) + 1)/(N + 1);
  }else if (type == "two-tailed") {
    pval = 2 * std::min((N * (1 - survival_pval) + 1)/(N + 1),(N * survival_pval + 1)/(N + 1));
  } else{
    Rcerr << "type must be one of the following: geq, leq, two-tailed or absolute\n";
  }
  return(pval); 
}

// ==============================================================================
//' @title Simulate Autoregressive Series
//' 
//' @description This function simulates an autoregresive series
//' 
//' @param mdl_h0 List containing series properties such as n: length, ar: number of autoregressive lags, stdev: standard deviation, mu: mean of process, phi: vector of autoregressive coefficients
//' @param burnin number of simulated observations to remove from begining. If using mu not equal to 0 it is recommended to use burnin to avoid dependence on initial value. Default is 200.
//' 
//' @return List with autoregressive series and its properties
//' 
//' @export
// [[Rcpp::export]]
List simuAR(List mdl_h0, int burnin = 200){
  int n = mdl_h0["n"];
  arma::vec phi = mdl_h0["phi"];
  double mu = mdl_h0["mu"];
  double std = mdl_h0["stdev"];
  int ar = mdl_h0["ar"];
  double intercept = mu*(1-sum(phi));
  arma::vec series(n+burnin, arma::fill::zeros);
  series.subvec(0, ar-1) = intercept + arma::randn(ar)*std;
  for (int xt = ar; xt<n+burnin; xt++){
    arma::vec ytmp = flipud(series.subvec((xt-ar),(xt-1)));
    series(xt) = as_scalar(intercept + trans(ytmp)*phi + arma::randn()*std);
  }
  arma::vec series_out = series.subvec(burnin, n+burnin-1);
  List simu_output;
  simu_output["y"] = series_out;
  simu_output["n"] = n;
  simu_output["phi"] = phi;
  simu_output["mu"] = mu;
  simu_output["stdev"] = std;
  simu_output["mdl"] = mdl_h0;
  return(simu_output);
}

// ==============================================================================
//' @title Simulate Markov-switching Autoregressive Series
//' 
//' @description This function simulates a Markov-switching autoregressive series
//' 
//' @param mdl_h0 List containing series properties such as n: length, ar: number of autoregressive lags, 
//' stdev: standard deviation, mu: mean of process, phi: vector of autoregressive coefficients, P: transition matrix (if type="markov")
//' or vector of component weights (if type="mixture"), k: number of regimes.
//' @param type determines type of St is a Markov process or a random mixture. Default is "markov" 
//' @param burnin number of simulated observations to remove from beginning. By assumption, series begins in state 1 
//' so it is recommended to use burnin to avoid dependence on this assumption. Default is 200.
//' 
//' @return List with Markov-switching autoregressive series and its properties
//' 
//' @export
// [[Rcpp::export]]
List simuMS(List mdl_h0, int burnin = 200){
  // ----- Obtain parameters
  int n = mdl_h0["n"];
  arma::vec phi = mdl_h0["phi"];
  arma::vec mu = mdl_h0["mu"];
  arma::vec std = mdl_h0["stdev"];
  int ar = mdl_h0["ar"];
  int k = mdl_h0["k"];
  arma::mat P = mdl_h0["P"];
  arma::vec pinf = limP(P, k);
  // vector for mean and standard dev at each time t
  arma::vec mu_t(n+burnin, arma::fill::zeros);
  arma::vec std_t(n+burnin, arma::fill::zeros);
  // Simulate data
  arma::vec series(n+burnin,arma::fill::zeros);
  arma::vec state_series(n+burnin,arma::fill::zeros);
  // initialize assuming series begins in state 1 (use burnin to reduce dependence on this assumption)
  int state = 0;
  arma::vec repar(ar, arma::fill::ones);
  series.subvec(0,ar-1) = mu(state) + arma::randn(ar)*std(state);
  mu_t.subvec(0,ar-1) = mu(state)*repar;
  std_t.subvec(0,ar-1) = std(state)*repar;
  // ----- Simulate series
  arma::vec repvec(k,arma::fill::ones);
  arma::vec state_ind = cumsum(repvec)-1;
  for (int xt = ar; xt<n+burnin; xt++){
    // Get new state
    arma::vec w_temp = P.col(state);
    arma::vec state_mat = cumsum(w_temp);
    state = as_scalar(state_ind(find(arma::randu() < state_mat, 1, "first")));
    state_series(xt) = state;
    // generate new obs
    arma::vec ytmp = flipud(series.subvec((xt-ar),(xt-1)));
    arma::vec mu_lag = flipud(mu_t.subvec((xt-ar),(xt-1))); 
    series(xt) = as_scalar(mu(state) + (trans(ytmp-mu_lag))*phi + arma::randn()*std(state));
    mu_t(xt) = mu(state);
    std_t(xt) = std(state);
  }
  // ----- Organize output
  arma::vec series_out = series.subvec(burnin,n+burnin-1);
  arma::vec state_series_out = state_series.subvec(burnin,n+burnin-1);
  arma::vec mu_t_out = mu_t.subvec(burnin,n+burnin-1);
  arma::vec std_t_out = std_t.subvec(burnin,n+burnin-1);
  List simu_output;
  simu_output["y"] = series_out;
  simu_output["St"] = state_series_out;
  simu_output["mu_t"] = mu_t_out;
  //simu_output["stdev_t"] = std_t_out;
  simu_output["P"] = P;
  simu_output["pinf"] = pinf;
  simu_output["mdl"] = mdl_h0;
  return(simu_output);
}

// ==============================================================================
//' @title Simulate VAR Model
//' 
//' 
//' @export
// [[Rcpp::export]]
List simuVAR(List mdl_h0, int burnin = 200){
  arma::vec mu = mdl_h0["mu"];
  arma::mat cov_mat = mdl_h0["sigma"];
  int Tsize = mdl_h0["n"];
  arma::mat phimat = mdl_h0["phi"];
  int ar = mdl_h0["ar"];
  // Number of time series variables
  int N = mu.n_elem;
  // random uniform for box-muller method
  double pi = arma::datum::pi;
  arma::mat U1(Tsize+burnin, N, arma::fill::randu);
  arma::mat U2(Tsize+burnin, N, arma::fill::randu);
  arma::mat eps = trans(arma::diagmat(sqrt(cov_mat))*trans(sqrt(-2*log(U1))%cos(2*pi*U2)));
  // add correlations
  arma::mat corr_mat = cov2corr(cov_mat);
  arma::mat C = chol(corr_mat, "lower");
  arma::mat eps_corr = trans(C*trans(eps));
  // Get companion form matrix
  arma::mat diagmat = arma::eye(N*(ar-1),N*(ar-1));
  arma::mat diagzero(N*(ar-1),N,arma::fill::zeros);
  arma::mat Mn = join_rows(diagmat,diagzero);
  arma::mat F = join_cols(phimat,Mn);
  // get constant vec 
  arma::vec repmu(ar,arma::fill::ones);
  arma::vec mu_tmp = vectorise(trans(repmu*trans(mu)));
  arma::vec nu_tmp = (arma::eye(N*ar,N*ar) - F)*mu_tmp;
  arma::vec nu = nu_tmp.subvec(0,N-1);
  // Simulate process 
  arma::mat Y(Tsize+burnin, N, arma::fill::zeros);
  Y.rows(0,ar-1) = repmu*trans(nu) + eps_corr.rows(0,ar-1);
  for (int xt = ar; xt<(Tsize+burnin); xt++){
    arma::mat Ytmp = flipud(Y.rows((xt-ar),(xt-1)));
    Y.row(xt) = trans(nu) + trans(vectorise(trans(Ytmp)))*trans(phimat) + eps_corr.row(xt);
  }
  arma::mat Y_out = Y.submat(burnin,0,Tsize+burnin-1,N-1);
  arma::mat eps_out = eps.submat(burnin,0,Tsize+burnin-1,N-1);
  arma::mat eps_corr_out = eps_corr.submat(burnin,0,Tsize+burnin-1,N-1);
  // Output
  List simuVAR_out;
  simuVAR_out["y"] = Y_out;
  simuVAR_out["corr_mat"] = corr_mat;
  simuVAR_out["cholesktmat"] = C;
  simuVAR_out["F_comp"] = F;
  simuVAR_out["mdl"] = mdl_h0;
  return(simuVAR_out);
}


// ==============================================================================
//' @title Simulate Markov-Switching VAR Model
//' 
//' 
//' @export
// [[Rcpp::export]]
List simuMSVAR(List mdl_h0, int burnin = 200){
  arma::mat mu = mdl_h0["mu"];
  List cov_matLs = mdl_h0["sigma"];
  int Tsize = mdl_h0["n"];
  arma::mat phimat = mdl_h0["phi"];
  int ar = mdl_h0["ar"];
  int N = mu.n_cols;
  int k = mdl_h0["k"];
  arma::mat P = mdl_h0["P"];
  arma::vec  pinf = limP(P, k);
  // vector for mean and standard dev at each time t
  arma::mat mu_t(Tsize+burnin, N, arma::fill::zeros);
  List sigma_t(Tsize+burnin);
  // Get companion form matrix
  arma::mat diagmat = arma::eye(N*(ar-1), N*(ar-1));
  arma::mat diagzero(N*(ar-1),N,arma::fill::zeros);
  arma::mat Mn = join_rows(diagmat,diagzero);
  arma::mat F = join_cols(phimat,Mn);
  // Get correlated normal errors
  List epsLs(k);
  List corr_mat(k);
  arma::mat nu(k, N, arma::fill::zeros);
  double pi = arma::datum::pi;
  arma::vec repmu(ar,arma::fill::ones);
  for (int xk = 0; xk<k; xk++){
    arma::mat U1(Tsize+burnin, N, arma::fill::randu);
    arma::mat U2(Tsize+burnin, N, arma::fill::randu);
    arma::mat cov_mat_k = cov_matLs[xk];
    arma::mat eps_k = trans(arma::diagmat(sqrt(cov_mat_k))*trans(sqrt(-2*log(U1))%cos(2*pi*U2)));
    // add correlations
    arma::mat corr_mat_k = cov2corr(cov_mat_k);
    corr_mat[xk] = corr_mat_k;
    arma::mat C = chol(corr_mat_k, "lower");
    arma::mat eps_corr = trans(C*trans(eps_k));
    epsLs[xk] = eps_corr;
  }
  // initialize assuming series begins in state 1 (use burnin to reduce dependence on this assumption)
  arma::mat Y(Tsize+burnin, N, arma::fill::zeros);
  arma::vec state_series(Tsize+burnin,arma::fill::zeros);
  int state = 0;
  mu_t.rows(0,ar-1) = repmu*mu.row(state);
  arma::mat eps_corr_k = epsLs[state];
  Y.rows(0,ar-1) = mu_t.rows(0,ar-1) + eps_corr_k.rows(0,ar-1);
  for (int xp = 0; xp<ar; xp++){
    sigma_t[xp] = cov_matLs[state];
  } 
  // ----- Simulate series
  arma::vec repvec(k,arma::fill::ones);
  arma::vec state_ind = cumsum(repvec)-1;
  for (int xt = ar; xt<(Tsize+burnin); xt++){
    // Get new state
    arma::vec w_temp = P.col(state);
    arma::vec state_mat = cumsum(w_temp);
    state = as_scalar(state_ind(find(arma::randu() < state_mat, 1, "first")));
    state_series(xt) = state;
    arma::mat eps_corr_k = epsLs[state];
    arma::mat Ytmp = flipud(Y.rows((xt-ar),(xt-1)));
    arma::mat mu_lag = flipud(mu_t.rows((xt-ar),(xt-1))); 
    Y.row(xt) = mu.row(state) + trans(vectorise(trans(Ytmp-mu_lag)))*trans(phimat) + eps_corr_k.row(xt);
    mu_t.row(xt) = mu.row(state);
    sigma_t[xt] = cov_matLs[state];
  }
  arma::mat Y_out = Y.submat(burnin,0,Tsize+burnin-1,N-1);
  arma::vec state_series_out = state_series.subvec(burnin,Tsize+burnin-1);
  arma::mat mu_t_out = mu_t.rows(burnin,Tsize+burnin-1);
  // Output
  List simuVAR_out;
  simuVAR_out["y"] = Y_out;
  simuVAR_out["F_comp"] = F;
  simuVAR_out["St"] = state_series_out;
  simuVAR_out["mu_t"] = mu_t_out;
  simuVAR_out["P"] = P;
  simuVAR_out["pinf"] = pinf;
  simuVAR_out["mdl"] = mdl_h0;
  return(simuVAR_out);
}

