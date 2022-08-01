#include <RcppArmadillo.h>
#include <RcppArmadilloExtensions/sample.h>
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
//' @description This function creates a random transition matrix
//' 
//' @param k number of regimes. Must be greater than or equal to 2. 
//' 
//' @return transition matrix with randomly generated entries.
//' 
//' @export
// [[Rcpp::export]]
arma::mat randP(int k){
  arma::mat P = reshape(arma::randu(k*k), k, k); 
  arma::vec PcolSums = trans(arma::sum(P,0));
  for (int xk = 0; xk<k; xk++){
    P.col(xk) = P.col(xk)/PcolSums(xk);
  } 
  return(P);
}

// ==============================================================================
//' @title Ergodic (limiting) probabilities of states
//' 
//' @description Takes a transition matrix and returns the limiting probabilities
//' 
//' @param P: matrix with transition probabilities
//' 
//' @return (k x 1) vector of limiting probabilities
//' 
//' @export
// [[Rcpp::export]]
arma::vec limP(arma::mat P){
  int nr = P.n_rows;
  int nc = P.n_cols;
  if (nc==nr){
    arma::mat onevec(1, nr, arma::fill::ones);
    arma::mat ep(1, nr+1, arma::fill::zeros);
    ep(0,nr) = 1;
    arma::mat Atmp = join_cols(arma::eye(nr,nr)-P,onevec);
    //arma::vec pinf = inv(trans(Atmp)*Atmp)*trans(Atmp)*trans(ep);
    arma::vec pinf = solve(trans(Atmp)*Atmp,trans(Atmp))*trans(ep);
    return (pinf);
  }else{
    stop("Input must be a square matrix");
  }
}



// ==============================================================================
//' @title Lagged Time Series Data
//' 
//' @description This function takes a (T x 1) vector Y and returns the (T-p x 1) vector y and the (T-p x p) matrix of lagged observations.
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
//' @title Parameter list for Markov-switching autoregressive model
//' 
//' @description This function takes the parameter vector of interest and converts it to a list with specific parameter vectors needed for univariate Markov-switching functions
//' 
//' @param theta: vector of parameters.
//' @param ar: umber of autoregressive lags
//' @param k: number of regimes
//' @param msmu: bool indicating is the mean switches with regime 
//' @param msvar: bool indicating is the variance switches with regime 
//' 
//' @return List with the mean, variance, transition matrix, limiting probabilities, and a vector of state indicators
//' 
//' @export
// [[Rcpp::export]]
List paramList_MSARmdl(arma::vec theta, int p, int k, bool msmu, bool msvar){
  Rcpp::Environment mstest("package:MSTest");
  Rcpp::Function argrid_MSARmdl = mstest["argrid_MSARmdl"];
  Rcpp::Function arP = mstest["arP"];
  // ----- Mean for each regime 
  arma::vec mu = theta.subvec(0, msmu*(k-1));
  // ----- Variance for each regime 
  arma::vec sig = theta.subvec(1+msmu*(k-1),1+msmu*(k-1)+msvar*(k-1));
  // ----- Phi vector
  arma::vec phi(p, arma::fill::zeros);
  phi =  theta.subvec(2+msmu*(k-1)+msvar*(k-1), 2+msmu*(k-1)+msvar*(k-1)+p-1);
  // ----- Transition probabilities 
  arma::mat P = reshape(theta.subvec(2+msmu*(k-1)+msvar*(k-1) + p, 2+msmu*(k-1)+msvar*(k-1) + p + k*k - 1),k, k);
  // Regime limiting probabilities
  arma::vec pinf = limP(P);
  // ----- Obtain AR consistent grid of Mu, Sigma and State indicators
  List musig_out = argrid_MSARmdl(mu, sig, k, p, msmu, msvar);
  arma::mat muAR = as<arma::mat>(musig_out["mu"]);
  arma::mat sigAR = as<arma::mat>(musig_out["sig"]);
  arma::vec state_ind = as<arma::vec>(musig_out["state_ind"]);
  // ----- Obtain AR consistent P and pinf
  arma::mat P_AR = as<arma::mat>(arP(P, k, p));
  arma::mat pinf_AR = limP(P_AR);
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
//' @title Parameter list for Markov-switching vector autoregressive model
//' 
//' @description This function takes the parameter vector of interest and converts it to a list with specific parameter vectors needed for multivariate Markov-switching functions
//' 
//' @param theta: vector of parameters.
//' @param q: number of time series 
//' @param ar: umber of autoregressive lags
//' @param k: number of regimes
//' @param msmu: bool indicating is the mean switches with regime 
//' @param msvar: bool indicating is the variance switches with regime 
//' 
//' @return List with the mean, variance, transition matrix, limiting probabilities, and a vector of state indicators
//' 
//' @export
// [[Rcpp::export]]
List paramList_MSVARmdl(arma::vec theta, int q, int ar, int k, bool msmu, bool msvar){
  Rcpp::Environment mstest("package:MSTest");
  Rcpp::Function argrid_MSVARmdl = mstest["argrid_MSVARmdl"];
  Rcpp::Function arP = mstest["arP"];
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
  arma::mat phi = trans(reshape(phi_tmp, q*ar, q));
  // ----- Transition probabilities 
  int PN = q+q*msmu*(k-1)+sigN+sigN*msvar*(k-1)+q*q*ar;
  arma::mat P = reshape(theta.subvec(PN, PN + k*k - 1), k, k);
  // Regime limiting probabilities
  arma::vec pinf = limP(P);
  // ----- Obtain AR consistent grid of Mu, Sigma and State indicators
  List musig_out = argrid_MSVARmdl(mu_k, sigma, k, ar, msmu, msvar);
  List muAR = musig_out["mu"];
  List sigAR = musig_out["sig"];
  arma::vec state_ind = as<arma::vec>(musig_out["state_ind"]);
  // ----- Obtain AR consistent P and pinf
  //int M = pow(k, ar+1);
  arma::mat P_AR = as<arma::mat>(arP(P, k, ar));
  arma::mat pinf_AR = limP(P_AR);
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
//' @title Markov-switching autoregressive model residuals
//' 
//' @description This function computes residuals of a Markov-switching autoregressive model
//' 
//' @param mdl List containing relevant parameters.
//' @param mu vector with mean in each regime.
//' @param k number of regimes. Must be greater than or equal to 2. 
//' 
//' @return (TxM) matrix of residuals in each regime M where M=k^(ar+1)
//' 
//' @export
// [[Rcpp::export]]
arma::mat calcResid_MSARmdl(List mdl, arma::mat mu, int k){
  // ---------- Initialize parameters
  arma::vec y = mdl["y"];
  int ar = mdl["p"];
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
//' @title Markov-switching vector autoregressive model residuals
//' 
//' @description This function computes residuals of a Markov-switching vector autoregressive model. Note that 
//' 
//' @param mdl List containing relevant parameters.
//' @param mu vector with mean in each regime.
//' @param k number of regimes. Must be greater than or equal to 2. 
//' 
//' @return List with M (Txq) matrix of residuals in each regime M where M=k^(ar+1)
//' 
//' @export
// [[Rcpp::export]]
List calcResid_MSVARmdl(List mdl, List mu, int k){
  // ---------- Initialize parameters
  arma::mat y = mdl["y"];
  int p = mdl["p"];
  arma::mat x = mdl["x"];
  arma::mat phi = mdl["phi"];
  int Tsize = y.n_rows;
  int q = y.n_cols;
  int M = pow(k, p+1); // number of regimes consistent with autoregressive structure (M=k if ar=0 and M=k^(ar+1) if ar>0)
  arma::mat repmu(Tsize, 1, arma::fill::ones);
  List z_mat(M); // [y(t) - mu_s(t))]
  List xz_mat(M); // [y(t-p) - mu_s(t-p))]
  List resid(M); 
  for (int xm = 0; xm<M; xm++){
    arma::mat mu_tmp = mu[xm]; 
    arma::mat y_tmp = y - repmu*trans(mu_tmp.col(0));
    arma::mat xz_tmp(Tsize, q*p, arma::fill::zeros); 
    for (int xp = 0; xp<p; xp++){
      xz_tmp.submat(0,q*xp,Tsize-1,q*xp+q-1) = x.submat(0,q*xp,Tsize-1,q*xp+q-1) - repmu*trans(mu_tmp.col(xp+1));
    }
    resid[xm] = y_tmp - xz_tmp*trans(phi);
    z_mat[xm] = y_tmp;
    xz_mat[xm] = xz_tmp;
  }
  return(resid);
}


// ==============================================================================
//' @title Initial values for Hidden Markov model
//' 
//' @description This function generates a random parameter vector to be used as initial values for a Hidden Markov model
//' 
//' @param mdl: List with parameter values of simple (one-regime) model. This includes
//'   - mu: vector of means
//'   - sigma: covariance matrix
//'   - msmu: boolean indicator. If TRUE, mean is function of markov process. If FALSE, mean is constant across regimes
//'   - msvar: boolean indicator. If TRUE, standard deviation is function of markov process. If FALSE, standard deviation is constant across regimes
//' @param k: number of regimes
//' 
//' @return vector of initial parameter values
//' 
//' @export
// [[Rcpp::export]]
arma::vec initVals_HMmdl(List mdl, int k){
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
    arma::mat sig_mat_tmp = covar_unvech((sigma_vec_tmp*0.1) + ((2*sigma_vec_tmp)-(sigma_vec_tmp*0.1))%arma::randu<arma::vec>(sigN), q);
    sig_mat_tmp = sig_mat_tmp*trans(sig_mat_tmp);
    sig_mat_tmp = sig_mat_tmp + q*arma::speye(q,q);
    sigma_0.row(0) = trans(covar_vech(sig_mat_tmp));
    for (int xk = 1; xk<k; xk++){
      arma::mat sig_mat_tmp = covar_unvech((sigma_vec_tmp*0.1) + ((2*sigma_vec_tmp)-(sigma_vec_tmp*0.1))%arma::randu<arma::vec>(sigN), q);
      sig_mat_tmp = sig_mat_tmp*trans(sig_mat_tmp);
      sig_mat_tmp = sig_mat_tmp + q*arma::speye(q,q);
      sigma_0.row(xk) = trans(covar_vech(sig_mat_tmp));
    }
  }
  arma::vec sigma_out = vectorise(trans(sigma_0));
  // create vector for initial values
  arma::vec theta_0 = join_vert(mu_out, sigma_out);
  arma::mat P_0 = randP(k);
  theta_0 = join_vert(theta_0, vectorise(P_0));
  return(theta_0);
}


// ==============================================================================
//' @title Initial values for Markov-switching autoregressive model
//' 
//' @description This function generates a random parameter vector to be used as initial values for a Markov-switching autoregressive model
//' 
//' @param mdl: List with parameter values of simple (one-regime) autoregressive model. This includes
//'   - phi: vector autoregressive coefficients
//'   - mu: mean of process 
//'   - stdev: standard deviation
//'   - msmu: boolean indicator. If TRUE, mean is function of markov process. If FALSE, mean is constant across regimes
//'   - msvar: boolean indicator. If TRUE, standard deviation is function of markov process. If FALSE, standard deviation is constant across regimes
//' @param k: number of regimes
//' 
//' @return vector of initial parameter values
//' 
//' @export
// [[Rcpp::export]]
arma::vec initVals_MSARmdl(List mdl, int k){
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
  arma::mat P_0 = randP(k);
  theta_0 = join_vert(theta_0, vectorise(P_0));
  return(theta_0);
}

// ==============================================================================
//' @title Initial values for Markov-switching vector autoregressive model
//' 
//' @description This function generates a random parameter vector to be used as initial values for a Markov-switching vector autoregressive model
//' 
//' @param mdl: List with parameter values of simple (one-regime) vector autoregressive model. This includes
//'   - phi: matrix autoregressive coefficients
//'   - mu: vector of means
//'   - sigma: covariance matrix
//'   - msmu: boolean indicator. If TRUE, mean is function of markov process. If FALSE, mean is constant across regimes
//'   - msvar: boolean indicator. If TRUE, standard deviation is function of markov process. If FALSE, standard deviation is constant across regimes
//' @param k: number of regimes
//' 
//' @return vector of initial parameter values
//' 
//' @export
// [[Rcpp::export]]
arma::vec initVals_MSVARmdl(List mdl, int k){
  arma::mat phi = mdl["phi"];
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
    arma::mat sig_mat_tmp = covar_unvech((sigma_vec_tmp*0.1) + ((2*sigma_vec_tmp)-(sigma_vec_tmp*0.1))%arma::randu<arma::vec>(sigN), q);
    sig_mat_tmp = sig_mat_tmp*trans(sig_mat_tmp);
    sig_mat_tmp = sig_mat_tmp + q*arma::speye(q,q);
    sigma_0.row(0) = trans(covar_vech(sig_mat_tmp));
    for (int xk = 1; xk<k; xk++){
      arma::mat sig_mat_tmp = covar_unvech((sigma_vec_tmp*0.1) + ((2*sigma_vec_tmp)-(sigma_vec_tmp*0.1))%arma::randu<arma::vec>(sigN), q);
      sig_mat_tmp = sig_mat_tmp*trans(sig_mat_tmp);
      sig_mat_tmp = sig_mat_tmp + q*arma::speye(q,q);
      sigma_0.row(xk) = trans(covar_vech(sig_mat_tmp));
    }
  }
  arma::vec sigma_out = vectorise(trans(sigma_0));
  // create vector for initial values
  arma::vec theta_0 = join_vert(mu_out, sigma_out);
  theta_0 = join_vert(theta_0, vectorise(trans(phi)));
  arma::mat P_0 = randP(k);
  theta_0 = join_vert(theta_0, vectorise(P_0));
  return(theta_0);
}


// ==============================================================================
//' @title Monte Carlo P-value
//' 
//' @description This function computes the Monte Carlo P-value
//' 
//' @param test_stat: test statistic under the alternative (e.g. S_0)
//' @param null_vec: (N x 1) vector with test statistic under the null hypothesis
//' @param type: type of test of test. options are: "geq" for right-tail test, "leq" for 
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
//' @title Simulate autoregressive process
//' 
//' @description This function simulates an autoregresive process
//' 
//' @param mdl_h0: List containing the following DGP parameters
//'   - n: length of series
//'   - mu: mean of process 
//'   - sigma: standard deviation of process
//'   - phi: vector of autoregressive coefficients
//'   
//' @param burnin: number of simulated observations to remove from beginning. Default is 100
//' 
//' @return List with simulated autoregressive series and its DGP parameters
//' 
//' @example /examples/simuAR_examples.R
//' @export
// [[Rcpp::export]]
List simuAR(List mdl_h0, int burnin = 100){
  // ----- DGP parameter
  arma::vec phi = mdl_h0["phi"];
  int Tsize = mdl_h0["n"];
  double mu = mdl_h0["mu"];
  double sigma = mdl_h0["sigma"];
  double stdev = sqrt(sigma);
  int p = phi.n_elem; 
  double intercept = mu*(1-sum(phi));
  // ----- Start simulation
  // pre-define variables
  arma::vec Y(Tsize+burnin, arma::fill::zeros);
  // simulate errors using box-Muller method
  double pi = arma::datum::pi;
  arma::mat U1(Tsize+burnin, 1, arma::fill::randu);
  arma::mat U2(Tsize+burnin, 1, arma::fill::randu);
  arma::vec eps = trans(trans(sqrt(-2*log(U1))%cos(2*pi*U2)));
  arma::vec eps_corr = eps*stdev;
  // simulate process
  Y.subvec(0, p-1) = mu + eps_corr.subvec(0, p-1);
  for (int xt = p; xt<(Tsize+burnin); xt++){
    arma::vec ytmp = flipud(Y.subvec((xt-p),(xt-1)));
    Y(xt) = as_scalar(intercept + trans(ytmp)*phi + eps_corr(xt));
  }
  // ----- Output
  arma::vec Y_out = Y.subvec(burnin, Tsize+burnin-1);
  arma::vec eps_corr_out = eps_corr.subvec(burnin, Tsize+burnin-1);
  List simu_output = clone(mdl_h0);
  simu_output["y"] = Y_out;
  simu_output["p"] = p;
  simu_output["resid"] = eps_corr_out;
  return(simu_output);
}

// ==============================================================================
//' @title Simulate Markov-switching autoregressive process
//' 
//' @description This function simulates a Markov-switching autoregressive process
//' 
//' @param mdl_h0: List containing the following DGP parameters
//'   - n: length of series
//'   - k: number of regimes
//'   - mu: (k x 1) vector with mean of process in each regime
//'   - sigma: (k x 1) vector with standard deviation of process in each regime
//'   - phi: vector of autoregressive coefficients
//'   - P: (k x k) transition matrix (columns must sum to one)
//'   
//' @param burnin: number of simulated observations to remove from beginning. Default is 100
//' 
//' @return List with simulated Markov-switching autoregressive process and its DGP properties
//' 
//' @example /examples/simuMSAR_examples.R
//' @export
// [[Rcpp::export]]
List simuMSAR(List mdl_h0, int burnin = 100){
  // ----- DGP parameter
  arma::vec phi = mdl_h0["phi"];
  arma::vec mu = mdl_h0["mu"];
  arma::vec sigma = mdl_h0["sigma"];
  arma::vec stdev = sqrt(sigma);
  arma::mat P = mdl_h0["P"];
  arma::vec pinf = limP(P);
  int Tsize = mdl_h0["n"];
  int k = mdl_h0["k"];
  int p = phi.n_elem; 
  // ----- Perform checks on DGP
  arma::vec check_Pcolsum = trans(arma::sum(P,0));
  if (max(abs(check_Pcolsum-1))>1e-8){
    stop("Columns of transition matrix 'P' must sum to 1.");
  }
  // ----- Start Simulation
  // pre-define variables
  arma::vec mu_t(Tsize+burnin, arma::fill::zeros);
  arma::vec stdev_t(Tsize+burnin, arma::fill::zeros);
  arma::vec Y(Tsize+burnin,arma::fill::zeros);
  arma::vec resid(Tsize+burnin,arma::fill::zeros);
  arma::vec state_series(Tsize+burnin,arma::fill::zeros);
  arma::vec repar(p,arma::fill::ones);
  // simulate errors using box-Muller method
  List epsLs(k);
  double pi = arma::datum::pi;
  for (int xk = 0; xk<k; xk++){
    double stdev_k = stdev(xk);
    arma::mat U1(Tsize+burnin, 1, arma::fill::randu);
    arma::mat U2(Tsize+burnin, 1, arma::fill::randu);
    arma::vec eps_k = trans(trans(sqrt(-2*log(U1))%cos(2*pi*U2)));
    arma::vec eps_corr = eps_k*stdev_k;
    epsLs[xk] = eps_corr;
  }
  // initialize assuming series begins in state 1 (use burnin to reduce dependence on this assumption)
  int state = 0;
  mu_t.rows(0,p-1) = repar*mu.row(state);
  stdev_t.rows(0,p-1) = repar*stdev.row(state);
  arma::vec eps_corr_k = epsLs[state];
  Y.rows(0,p-1) = mu_t.rows(0,p-1) + eps_corr_k.rows(0,p-1);
  // simulate process
  arma::vec repk(k,arma::fill::ones);
  arma::vec state_ind = cumsum(repk)-1;
  for (int xt = p; xt<(Tsize+burnin); xt++){
    // Get new state
    arma::vec w_temp = P.col(state);
    arma::vec state_mat = cumsum(w_temp);
    state = as_scalar(state_ind(find(arma::randu() < state_mat, 1, "first")));
    state_series(xt) = state;
    // generate new obs
    arma::vec ytmp = flipud(Y.subvec((xt-p),(xt-1)));
    arma::vec mu_lag = flipud(mu_t.subvec((xt-p),(xt-1))); 
    arma::vec eps_corr_k = epsLs[state];
    resid(xt) = eps_corr_k(xt);
    Y(xt) = as_scalar(mu(state) + (trans(ytmp-mu_lag))*phi + resid(xt));
    mu_t(xt) = mu(state);
    stdev_t(xt) = stdev(state);
  }
  // ----- Output
  arma::vec Y_out = Y.subvec(burnin,Tsize+burnin-1);
  arma::vec resid_out = resid.subvec(burnin,Tsize+burnin-1);
  arma::vec state_series_out = state_series.subvec(burnin,Tsize+burnin-1);
  arma::vec mu_t_out = mu_t.subvec(burnin,Tsize+burnin-1);
  arma::vec stdev_t_out = stdev_t.subvec(burnin,Tsize+burnin-1);
  List simu_output = clone(mdl_h0);
  simu_output["y"] = Y_out;
  simu_output["p"] = p;
  simu_output["resid"] = resid_out;
  simu_output["St"] = state_series_out;
  simu_output["mu_t"] = mu_t_out;
  simu_output["stdev_t"] = stdev_t_out;
  simu_output["pinf"] = pinf;
  return(simu_output);
}

// ==============================================================================
//' @title Simulate VAR process
//' 
//' @description This function simulates a vector autoregresive process
//' 
//' @param mdl_h0: List containing the following DGP parameters
//'   - n: length of series
//'   - mu: (q x 1) vector of means 
//'   - sigma: (q x q) covariance matrix 
//'   - phi: (q x qp) matrix of autoregressive coefficients
//'   - ar: number of autoregressive lags (i.e, p)
//'   
//' @param burnin: number of simulated observations to remove from beginning. Default is 100
//' 
//' @return List with simulated vector autoregressive series and its DGP parameters
//' 
//' @example /examples/simuVAR_examples.R
//' @export
// [[Rcpp::export]]
List simuVAR(List mdl_h0, int burnin = 100){
  // ----- DGP parameters
  arma::vec mu = mdl_h0["mu"];
  arma::mat cov_mat = mdl_h0["sigma"];
  arma::mat phimat = mdl_h0["phi"];
  int Tsize = mdl_h0["n"];
  int p = mdl_h0["p"];
  int N = mu.n_elem;
  // companion matrix
  arma::mat diagmat = arma::eye(N*(p-1),N*(p-1));
  arma::mat diagzero(N*(p-1),N,arma::fill::zeros);
  arma::mat Mn = join_rows(diagmat,diagzero);
  arma::mat F = join_cols(phimat,Mn);
  // constant vec
  arma::vec repmu(p,arma::fill::ones);
  arma::vec mu_tmp = vectorise(trans(repmu*trans(mu)));
  arma::vec nu_tmp = (arma::eye(N*p,N*p) - F)*mu_tmp;
  arma::vec nu = nu_tmp.subvec(0,N-1);
  // ----- Simulate VAR process
  // simulate errors using box-Muller method
  double pi = arma::datum::pi;
  arma::mat U1(Tsize+burnin, N, arma::fill::randu);
  arma::mat U2(Tsize+burnin, N, arma::fill::randu);
  arma::mat eps = trans(trans(sqrt(-2*log(U1))%cos(2*pi*U2)));
  // add correlations
  arma::mat C = chol(cov_mat, "lower");
  arma::mat eps_corr = trans(C*trans(eps));
  // simulate process
  arma::mat Y(Tsize+burnin, N, arma::fill::zeros);
  Y.rows(0,p-1) = repmu*trans(nu) + eps_corr.rows(0,p-1);
  for (int xt = p; xt<(Tsize+burnin); xt++){
    arma::mat Ytmp = flipud(Y.rows((xt-p),(xt-1)));
    Y.row(xt) = trans(nu) + trans(vectorise(trans(Ytmp)))*trans(phimat) + eps_corr.row(xt);
  }
  // ----- Output
  arma::mat Y_out = Y.submat(burnin,0,Tsize+burnin-1,N-1);
  arma::mat eps_corr_out = eps_corr.submat(burnin,0,Tsize+burnin-1,N-1);
  List simuVAR_out = clone(mdl_h0);
  simuVAR_out["y"] = Y_out;
  simuVAR_out["resid"] = eps_corr_out;
  simuVAR_out["F_comp"] = F;
  return(simuVAR_out);
}


// ==============================================================================
//' @title Simulate Markov-switching vector autoregressive process
//' 
//' @description This function simulates a Markov-switching vector autoregressive process
//' 
//' @param mdl_h0: List containing the following DGP parameters
//'   - n: length of series
//'   - k: number of regimes
//'   - mu: (k x q) matrix of means 
//'   - sigma: List with k (q x q) covariance matrices
//'   - phi: (q x qp) matrix of autoregressive coefficients
//'   - ar: number of autoregressive lags (i.e, p)
//'   - P: (k x k) transition matrix (columns must sum to one)
//'   
//' @param burnin: number of simulated observations to remove from beginning. Default is 100
//' 
//' @return List with simulated vector autoregressive series and its DGP parameters
//' 
//' @example /examples/simuMSVAR_examples.R
//' @export
// [[Rcpp::export]]
List simuMSVAR(List mdl_h0, int burnin = 100){
  // ----- DGP parameters
  arma::mat mu = mdl_h0["mu"];
  List cov_matLs = mdl_h0["sigma"];
  arma::mat phimat = mdl_h0["phi"];
  arma::mat P = mdl_h0["P"];
  arma::vec  pinf = limP(P);
  int Tsize = mdl_h0["n"];
  int p = mdl_h0["p"];
  int k = mdl_h0["k"];
  int N = mu.n_cols;
  // companion form matrix
  arma::mat diagmat = arma::eye(N*(p-1), N*(p-1));
  arma::mat diagzero(N*(p-1),N,arma::fill::zeros);
  arma::mat Mn = join_rows(diagmat,diagzero);
  arma::mat F = join_cols(phimat,Mn);
  // ----- Perform checks on DGP
  arma::vec check_Pcolsum = trans(arma::sum(P,0));
  if (max(abs(check_Pcolsum-1))>1e-8){
    stop("Columns of transition matrix 'P' must sum to 1.");
  }
  // ----- Start Simulation
  // pre-define variables
  arma::mat mu_t(Tsize+burnin, N, arma::fill::zeros);
  List sigma_t(Tsize+burnin);
  arma::mat Y(Tsize+burnin, N, arma::fill::zeros);
  arma::mat resid(Tsize+burnin, N, arma::fill::zeros);
  arma::vec state_series(Tsize+burnin,arma::fill::zeros);
  arma::vec repar(p,arma::fill::ones);
  // simulate errors using box-Muller method
  List epsLs(k);
  double pi = arma::datum::pi;
  for (int xk = 0; xk<k; xk++){
    arma::mat U1(Tsize+burnin, N, arma::fill::randu);
    arma::mat U2(Tsize+burnin, N, arma::fill::randu);
    arma::mat eps_k = trans(trans(sqrt(-2*log(U1))%cos(2*pi*U2)));
    // add correlations
    arma::mat cov_mat_k = cov_matLs[xk];
    arma::mat C = chol(cov_mat_k, "lower");
    arma::mat eps_corr = trans(C*trans(eps_k));
    epsLs[xk] = eps_corr;
  }
  // initialize assuming series begins in state 1 (use burnin to reduce dependence on this assumption)
  int state = 0;
  mu_t.rows(0,p-1) = repar*mu.row(state);
  arma::mat eps_corr_k = epsLs[state];
  Y.rows(0,p-1) = mu_t.rows(0,p-1) + eps_corr_k.rows(0,p-1);
  for (int xp = 0; xp<p; xp++){
    sigma_t[xp] = cov_matLs[state];
  } 
  // simulate process
  arma::vec repk(k,arma::fill::ones);
  arma::vec state_ind = cumsum(repk)-1;
  for (int xt = p; xt<(Tsize+burnin); xt++){
    // Get new state
    arma::vec w_temp = P.col(state);
    arma::vec state_mat = cumsum(w_temp);
    state = as_scalar(state_ind(find(arma::randu() < state_mat, 1, "first")));
    state_series(xt) = state;
    arma::mat eps_corr_k = epsLs[state];
    arma::mat Ytmp = flipud(Y.rows((xt-p),(xt-1)));
    arma::mat mu_lag = flipud(mu_t.rows((xt-p),(xt-1))); 
    resid.row(xt) = eps_corr_k.row(xt);
    Y.row(xt) = mu.row(state) + trans(vectorise(trans(Ytmp-mu_lag)))*trans(phimat) + resid.row(xt);
    mu_t.row(xt) = mu.row(state);
    sigma_t[xt] = cov_matLs[state];
  }
  // ----- Output
  arma::mat Y_out = Y.submat(burnin,0,Tsize+burnin-1,N-1);
  arma::mat resid_out = resid.submat(burnin,0,Tsize+burnin-1,N-1);
  arma::vec state_series_out = state_series.subvec(burnin,Tsize+burnin-1);
  arma::mat mu_t_out = mu_t.rows(burnin,Tsize+burnin-1);
  List sigma_t_out(Tsize);
  for (int xt = burnin; xt<(Tsize+burnin); xt++){
    sigma_t_out[xt-burnin] = sigma_t[xt];
  } 
  List simuVAR_out = clone(mdl_h0);
  simuVAR_out["y"] = Y_out;
  simuVAR_out["resid"] = resid_out;
  simuVAR_out["F_comp"] = F;
  simuVAR_out["St"] = state_series_out;
  simuVAR_out["mu_t"] = mu_t_out;
  simuVAR_out["sigma_t"] = sigma_t_out;
  simuVAR_out["pinf"] = pinf;
  return(simuVAR_out);
}


// ==============================================================================
//' @title Simulate normally distributed process
//' 
//' @description This function simulates a normally distributed process
//' 
//' @param mdl_h0: List containing the following DGP parameters
//'   - n: length of series
//'   - mu: (q x 1) vector of means 
//'   - sigma: (q x q) covariance matrix. If simulating a univariate process (i.e., length of mu is 1), this should be the standard devation and not the variance.
//'   
//' @return List with simulated series and its DGP parameters
//' 
//' @example /examples/simuNorm_examples.R
//' @export
// [[Rcpp::export]]
List simuNorm(List mdl_h0){
  // ----- DGP parameter
  arma::vec mu = mdl_h0["mu"];
  int Tsize = mdl_h0["n"];
  int N = mu.n_elem;
  // ----- Start simulation
  // pre-define variables
  arma::mat Y; 
  arma::mat eps_corr;
  arma::vec repT(Tsize, arma::fill::ones);
  // simulate errors using box-Muller method
  double pi = arma::datum::pi;
  arma::mat U1(Tsize, N, arma::fill::randu);
  arma::mat U2(Tsize, N, arma::fill::randu);
  arma::mat eps = trans(trans(sqrt(-2*log(U1))%cos(2*pi*U2)));
  if (N>1){
    arma::mat cov_mat = mdl_h0["sigma"];
    // add correlations
    arma::mat C = chol(cov_mat, "lower");
    eps_corr = trans(C*trans(eps));
    // simulate process
    Y  =  repT*trans(mu) + eps_corr;  
  }else if (N==1){
    double stdev = mdl_h0["sigma"];
    // add variance
    eps_corr = eps*stdev;
    // simulate process
    Y = repT*mu + eps_corr;
  }
  // ----- Output
  List simuNorm_out = clone(mdl_h0);
  simuNorm_out["y"] = Y;
  simuNorm_out["resid"] = eps_corr;
  return(simuNorm_out);
}


// ==============================================================================
//' @title Simulate Hidden Markov model with normally distributed errors
//' 
//' @description This function simulates a Hidden Markov Model process
//' 
//' @param mdl_h0: List containing the following DGP parameters
//'   - n: length of series
//'   - k: number of regimes
//'   - mu: (k x q) vector of means 
//'   - sigma: (q x q) covariance matrix. If simulating a univariate process (i.e., length of mu is 1), this should be the standard devation and not the variance.
//'   - P: (k x k) transition matrix (columns must sum to one)
//'   
//' @param burnin: number of simulated observations to remove from beginning. Default is 100
//' 
//' @return List with simulated series and its DGP parameters
//' 
//' @example /examples/simuHMM_examples.R
//' @export
// [[Rcpp::export]]
List simuHMM(List mdl_h0, int burnin = 100){
  // ----- DGP parameter
  arma::mat mu = mdl_h0["mu"];
  List cov_matLs = mdl_h0["sigma"];
  arma::mat P = mdl_h0["P"];
  arma::vec  pinf = limP(P);
  int Tsize = mdl_h0["n"];
  int k = mdl_h0["k"];
  int N = mu.n_cols;
  // ----- Perform checks on DGP
  arma::vec check_Pcolsum = trans(arma::sum(P,0));
  if (max(abs(check_Pcolsum-1))>1e-8){
    stop("Columns of transition matrix 'P' must sum to 1.");
  }
  // ----- Start Simulation
  // pre-define variables
  arma::mat mu_t(Tsize+burnin, N, arma::fill::zeros);
  List sigma_t(Tsize+burnin);
  arma::mat Y(Tsize+burnin, N, arma::fill::zeros);
  arma::mat resid(Tsize+burnin, N, arma::fill::zeros);
  arma::vec state_series(Tsize+burnin, arma::fill::zeros);
  // simulate errors using box-Muller method
  List epsLs(k);
  double pi = arma::datum::pi;
  for (int xk = 0; xk<k; xk++){
    arma::mat U1(Tsize+burnin, N, arma::fill::randu);
    arma::mat U2(Tsize+burnin, N, arma::fill::randu);
    arma::mat eps_k = trans(trans(sqrt(-2*log(U1))%cos(2*pi*U2)));
    // add correlations
    arma::mat cov_mat_k = cov_matLs[xk];
    arma::mat C = chol(cov_mat_k, "lower");
    arma::mat eps_corr = trans(C*trans(eps_k));
    epsLs[xk] = eps_corr;
  }
  // initialize assuming series begins in state 1 (use burnin to reduce dependence on this assumption)
  int state = 0;
  // simulate process
  arma::vec repk(k,arma::fill::ones);
  arma::vec state_ind = cumsum(repk)-1;
  for (int xt = 0; xt<(Tsize+burnin); xt++){
    // Get new state
    arma::vec w_temp = P.col(state);
    arma::vec state_mat = cumsum(w_temp);
    state = as_scalar(state_ind(find(arma::randu() < state_mat, 1, "first")));
    state_series(xt) = state;
    arma::mat eps_corr_k = epsLs[state];
    resid.row(xt) = eps_corr_k.row(xt);
    Y.row(xt) = mu.row(state) + resid.row(xt);
    mu_t.row(xt) = mu.row(state);
    sigma_t[xt] = cov_matLs[state];
  }
  arma::mat Y_out = Y.submat(burnin,0,Tsize+burnin-1,N-1);
  arma::mat resid_out = resid.submat(burnin,0,Tsize+burnin-1,N-1);
  arma::vec state_series_out = state_series.subvec(burnin,Tsize+burnin-1);
  arma::mat mu_t_out = mu_t.rows(burnin,Tsize+burnin-1);
  List sigma_t_out(Tsize);
  for (int xt = burnin; xt<(Tsize+burnin); xt++){
    sigma_t_out[xt-burnin] = sigma_t[xt];
  } 
  // ----- Output
  List simuHMM_out = clone(mdl_h0);
  simuHMM_out["y"] = Y_out;
  simuHMM_out["resid"] = resid_out;
  simuHMM_out["St"] = state_series_out;
  simuHMM_out["mu_t"] = mu_t_out;
  simuHMM_out["sigma_t"] = sigma_t_out;
  simuHMM_out["pinf"] = pinf;
  return(simuHMM_out);
}

// ==============================================================================
//' @title Normal log-likelihood objective function 
//' 
//' @description This function computes the log-likelihood for a normally distributed model
//' 
//' @param theta vector of model parameters
//' @param mdl List with model attributes
//' 
//' @return log-likelihood given data
//' 
//' @export
// [[Rcpp::export]]
double logLike_Nmdl(arma::vec theta, List mdl){
  // ---------- model parameters
  arma::mat y = mdl["y"];
  int Tsize = y.n_rows;
  int q = y.n_cols;
  // ---------- pre-define variables
  int sigN = (q*(q+1))/2;
  arma::vec mu = theta.subvec(0, q-1);
  arma::vec sig = theta.subvec(q,q+sigN-1);
  arma::mat sigma = covar_unvech(sig, q);
  arma::mat repmu(Tsize, 1, arma::fill::ones);
  arma::mat resid = y-repmu*trans(mu);
  // ---------- Compute log-likehood
  double pi = arma::datum::pi;
  arma::vec f_t(Tsize, arma::fill::zeros);
  for (int xt = 0; xt<Tsize; xt++){
    f_t(xt) = as_scalar((1/sqrt(det(sigma)*pow(2*pi,q)))*
      exp(-0.5*(resid.row(xt)*inv(sigma)*trans(resid.row(xt)))));
  }
  double logLike = sum(log(f_t));
  return(logLike);  
}

// ==============================================================================
//' @title Autoregressive log-likelihood objective function
//' 
//' @description This function computes the log-likelihood for an autoregressive model
//' 
//' @param theta vector of model parameters
//' @param mdl List with model attributes
//' 
//' @return log-likelihood given data
//' 
//' @export
// [[Rcpp::export]]
double logLike_ARmdl(arma::vec theta, List mdl){
  // ---------- Initialize parameters
  arma::vec y = mdl["y"];
  arma::mat x = mdl["x"];
  int p = mdl["p"];
  int Tsize = y.n_elem;
  // ---------- Compute log-likehood
  double logLike;
  double pi = arma::datum::pi;
  arma::mat repmu(Tsize, p,arma::fill::ones);
  logLike = sum(log((1/sqrt(2*pi*theta(1)))*
    exp(-pow((y - theta(0)) - (x-(theta(0)*repmu))*theta.subvec(2, 2+p-1),2)/(2*theta(1)))));
  return(logLike);
}

// ==============================================================================
//' @title Vector autoregressive log-likelihood objective function 
//' 
//' @description This function computes the log-likelihood for a vector autoregressive model
//' 
//' @param theta vector of model parameters
//' @param mdl List with model attributes
//' 
//' @return log-likelihood given data
//' 
//' @export
// [[Rcpp::export]]
double logLike_VARmdl(arma::vec theta, List mdl){
  // ---------- model parameters
  arma::mat y = mdl["y"];
  arma::mat x = mdl["x"];
  int p = mdl["p"];
  int Tsize = y.n_rows;
  int q = y.n_cols;
  // ---------- pre-define variables
  int sigN = (q*(q+1))/2;
  arma::vec mu = theta.subvec(0, q-1);
  arma::vec sig = theta.subvec(q,q+sigN-1);
  arma::mat sigma = covar_unvech(sig, q);
  int phiN = q+sigN;
  arma::vec phi_tmp = theta.subvec(phiN, phiN+q*q*p-1);
  arma::mat phi = reshape(phi_tmp, q*p, q);
  arma::mat repmu(Tsize, 1, arma::fill::ones);
  arma::mat z = y-repmu*trans(mu);
  arma::mat xz_tmp(Tsize, q*p, arma::fill::zeros); 
  // ---------- compute residual
  for (int xp = 0; xp<p; xp++){
    xz_tmp.submat(0,q*xp,Tsize-1,q*xp+q-1) = x.submat(0,q*xp,Tsize-1,q*xp+q-1) - repmu*trans(mu);
  }
  arma::mat resid = z - xz_tmp*phi;
  // ---------- Compute log-likehood
  double pi = arma::datum::pi;
  arma::vec f_t(Tsize, arma::fill::zeros);
  for (int xt = 0; xt<Tsize; xt++){
    f_t(xt) = as_scalar((1/sqrt(det(sigma)*pow(2*pi,q)))*
      exp(-0.5*(resid.row(xt)*inv(sigma)*trans(resid.row(xt)))));
  }
  double logLike = sum(log(f_t));
  return(logLike);  
}

// ==============================================================================
//' @title Hidden Markov model log-likelihood function 
//' 
//' @description This function computes the log-likelihood for a markov-switching autoregressive model
//' 
//' @param theta vector of model parameters
//' @param mdl List with model attributes
//' @param k integer determining the number of regimes
//'  
//' @return List with log-likelihood (log-likelihood), and smoothed probabilities of each regime (xi_t_T).
//' 
//' @export
// [[Rcpp::export]]
double logLike_HMmdl(arma::vec theta, List mdl, int k){
  // ---------- Initialize parameters
  arma::mat y = mdl["y"];
  int q = y.n_cols;
  int Tsize = y.n_rows;
  bool msmu = mdl["msmu"];
  bool msvar = mdl["msvar"];
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
  // ----- Transition probabilities
  int th_len = theta.n_elem;
  arma::mat P = reshape(theta.subvec(th_len - k*k, th_len - 1), k, k);
  // Regime limiting probabilities
  arma::vec pinf = limP(P);
  // ----- Compute Residuals
  List eps(k); 
  arma::mat repmu(Tsize, 1, arma::fill::ones);
  for (int xk = 0; xk<k; xk++){
    arma::mat eps_tmp = y - repmu*mu_k.row(xk);
    eps[xk] =  eps_tmp;
  } 
  // ----- Begin Calculating log-likelihood & filtered probabilities
  arma::mat eta(Tsize, k, arma::fill::zeros);       // [eq. 22.4.2]
  arma::mat f_t(Tsize, 1, arma::fill::zeros);       // [eq. 22.4.8]
  arma::mat xi_t_t(Tsize, k, arma::fill::zeros);    // [eq. 22.4.5]
  arma::mat xi_tp1_t(Tsize, k, arma::fill::zeros);  // [eq. 22.4.6]
  double pi = arma::datum::pi;
  arma::vec xi_t_tm1 = pinf;
  for (int xt = 0; xt<Tsize; xt++){
    for (int xk = 0; xk<k; xk++){
      arma::mat eps_k = eps[xk];
      arma::mat sigma_k = sigma[xk];
      //eta(xt,xk) = as_scalar((1/(sqrt(pow(2*pi,q)*det(sigma_k))))*exp(-0.5*(eps_k.row(xt)*inv(sigma_k)*trans(eps_k.row(xt)))));
      eta(xt,xk) = as_scalar((1/(sqrt(pow(2*pi,q)*det(sigma_k))))*exp(-0.5*(eps_k.row(xt)*solve(sigma_k,trans(eps_k.row(xt))))));
    }
    arma::vec xi_eta = xi_t_tm1%trans(eta.row(xt));
    f_t.row(xt) = sum(xi_eta);
    double f_tmp = sum(xi_eta);
    xi_t_t.row(xt) = trans(xi_eta/f_tmp);
    xi_t_tm1 = P*trans(xi_t_t.row(xt));
    xi_tp1_t.row(xt) = trans(xi_t_tm1);
  }
  // ----- log-Likelihood
  arma::vec logf = log(f_t);
  double logL = sum(logf);     //[eq. 22.4.7]
  return(logL);
}

// ==============================================================================
//' @title Hidden Markov model log-likelihood function  (minimization version)
//' 
//' @description This function computes the (negative) log-likelihood for a markov-switching autoregressive model
//' 
//' @param theta vector of model parameters
//' @param mdl List with model attributes
//' @param k integer determining the number of regimes
//' 
//' @return negative log-likelihood given data
//' 
//' @export
// [[Rcpp::export]]
double logLike_HMmdl_min(arma::vec theta, List mdl, int k){
  double logL_negative = -logLike_HMmdl(theta, mdl, k);
  return(logL_negative);
}

// ==============================================================================
//' @title Markov-switching autoregressive log-likelihood objective function
//' 
//' @description This function computes the log-likelihood for a markov-switching autoregressive model
//' 
//' @param theta vector of model parameters
//' @param mdl List with model attributes
//' @param k integer determining the number of regimes
//' 
//' @return log-likelihood given data
//' 
//' @export
// [[Rcpp::export]]
double logLike_MSARmdl(arma::vec theta, List mdl, int k){
  Rcpp::Environment mstest("package:MSTest");
  Rcpp::Function arP = mstest["arP"];
  Rcpp::Function argrid_MSARmdl = mstest["argrid_MSARmdl"];
  // ---------- Initialize parameters
  arma::vec y = mdl["y"];
  arma::mat x = mdl["x"];
  int ar = mdl["p"];
  bool msmu = mdl["msmu"];
  bool msvar = mdl["msvar"];
  int Tsize = y.n_elem;
  int M = pow(k, ar+1);
  // ----- Mean for each regime 
  arma::vec mu = theta.subvec(0, msmu*(k-1));
  // ----- Variance for each regime 
  arma::vec sig = theta.subvec(1+msmu*(k-1),1+msmu*(k-1)+msvar*(k-1));
  // ----- Phi vector
  arma::vec phi =  theta.subvec(2+msmu*(k-1)+msvar*(k-1), 2+msmu*(k-1)+msvar*(k-1)+ar-1);
  // ----- Transition probabilities 
  arma::mat P = reshape(theta.subvec(2+msmu*(k-1)+msvar*(k-1) + ar, 2+msmu*(k-1)+msvar*(k-1) + ar + k*k - 1),k, k);
  // Regime limiting probabilities
  arma::vec pinf = limP(P);
  // ----- Obtain AR consistent grid of Mu, Sigma and State indicators
  List musig_out = argrid_MSARmdl(mu, sig, k, ar, msmu, msvar);
  arma::mat muAR = as<arma::mat>(musig_out["mu"]);
  arma::mat sigAR = as<arma::mat>(musig_out["sig"]);
  arma::vec state_ind = as<arma::vec>(musig_out["state_ind"]);
  // ----- Obtain AR consistent P and pinf
  arma::mat P_AR = as<arma::mat>(arP(P, k, ar));
  arma::mat pinf_AR = limP(P_AR);
  // ----- Compute residuals in each regime
  arma::mat repvec(1, M, arma::fill::ones);
  arma::mat repmu(Tsize, 1, arma::fill::ones);
  arma::mat ms_y = y*repvec;
  arma::mat eps(Tsize, M, arma::fill::zeros);
  // ---------- Compute Residuals 
  arma::mat z = ms_y - repmu*trans(muAR.col(0)); // [y(t) - mu_s(t))]
  arma::mat xz(Tsize, M, arma::fill::zeros);
  for (int xkp = 0; xkp<M; xkp++){
    arma::mat zx_tmp = x - repmu*muAR.submat(xkp,1,xkp,ar); // [y(t-i) - mu_s(t-i))]
    xz.col(xkp) = zx_tmp*phi;
  }
  eps = z - xz;
  // ----- Begin Calculating likelihood and log-likelihood
  arma::mat eta(Tsize, M, arma::fill::zeros);       // [eq. 22.4.2]
  arma::mat f_t(Tsize, 1, arma::fill::zeros);       // [eq. 22.4.8]
  arma::mat xi_t_t(Tsize, M, arma::fill::zeros);    // [eq. 22.4.5]
  arma::mat xi_tp1_t(Tsize, M, arma::fill::zeros);  // [eq. 22.4.6]
  double pi = arma::datum::pi;
  arma::vec xi_t_tm1 = pinf_AR;
  for (int xt = 0; xt<Tsize; xt++){
    eta.row(xt) = trans(1/sqrt(2*pi*(sigAR)))%exp((pow(eps.row(xt),2)%trans(-1/(2*sigAR))));
    arma::vec xi_eta = xi_t_tm1%trans(eta.row(xt));
    f_t.row(xt) = sum(xi_eta);
    double f_tmp = sum(xi_eta);
    xi_t_t.row(xt) = trans(xi_eta/f_tmp);
    xi_t_tm1 = P_AR * trans(xi_t_t.row(xt));
    xi_tp1_t.row(xt) = trans(xi_t_tm1);
  }
  // log-Likelihood
  arma::vec logf = log(f_t);
  double logL = sum(logf);     //[eq. 22.4.7]
  return(logL);
  
}

// ==============================================================================
//' @title Markov-switching autoregressive log-likelihood objective function (minimization version)
//' 
//' @description This function computes the (negative) log-likelihood for a markov-switching autoregressive model
//' 
//' @param theta vector of model parameters
//' @param mdl List with model attributes
//' @param k integer determining the number of regimes
//' 
//' @return negative log-likelihood given data
//' 
//' @export
// [[Rcpp::export]]
double logLike_MSARmdl_min(arma::vec theta, List mdl, int k){
  double logL_negative = -logLike_MSARmdl(theta, mdl, k);
  return(logL_negative);
}


// ==============================================================================
//' @title Markov-switching vector autoregressive log-likelihood objective function
//' 
//' @description This function computes the log-likelihood for a markov-switching vector autoregressive model
//' 
//' @param theta vector of model parameters
//' @param mdl List with model attributes
//' @param k integer determining the number of regimes
//' 
//' @return log-likelihood given data
//' 
//' @export
// [[Rcpp::export]]
double logLike_MSVARmdl(arma::vec theta, List mdl, int k){
  Rcpp::Environment mstest("package:MSTest");
  Rcpp::Function arP = mstest["arP"];
  Rcpp::Function argrid_MSVARmdl = mstest["argrid_MSVARmdl"];
  // ---------- Initialize parameters
  arma::mat y = mdl["y"];
  arma::mat x = mdl["x"];
  int q = y.n_cols;
  int Tsize = y.n_rows;
  int ar = mdl["p"];
  int M = pow(k, ar+1);
  bool msmu = mdl["msmu"];
  bool msvar = mdl["msvar"];
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
  arma::vec pinf = limP(P);
  // ----- Obtain AR consistent grid of Mu, Sigma and State indicators
  List musig_out = argrid_MSVARmdl(mu_k, sigma, k, ar, msmu, msvar);
  List muAR = musig_out["mu"];
  List sigAR = musig_out["sig"];
  arma::vec state_ind = as<arma::vec>(musig_out["state_ind"]);
  // ----- Obtain AR consistent P and pinf
  arma::mat PAR = as<arma::mat>(arP(P, k, ar));
  arma::mat pinfAR = limP(PAR);
  // ----- Compute residuals in each regime
  arma::mat repmu(Tsize, 1, arma::fill::ones);
  List z_mat(M); // [y(t) - mu_s(t))]
  List xz_mat(M); // [y(t-p) - mu_s(t-p))]
  List eps(M); 
  for (int xm = 0; xm<M; xm++){
    arma::mat mu_tmp = muAR[xm]; 
    arma::mat y_tmp = y - repmu*trans(mu_tmp.col(0));
    arma::mat xz_tmp(Tsize, q*ar, arma::fill::zeros); 
    for (int xp = 0; xp<ar; xp++){
      xz_tmp.submat(0,q*xp,Tsize-1,q*xp+q-1) = x.submat(0,q*xp,Tsize-1,q*xp+q-1) - repmu*trans(mu_tmp.col(xp+1));
    }
    eps[xm] = y_tmp - xz_tmp*phi;
    z_mat[xm] = y_tmp;
    xz_mat[xm] = xz_tmp;
  }
  // ----- Begin Calculating log-likelihood & filtered probabilities
  arma::mat eta(Tsize, M, arma::fill::zeros);       // [eq. 22.4.2]
  arma::mat f_t(Tsize, 1, arma::fill::zeros);       // [eq. 22.4.8]
  arma::mat xi_t_t(Tsize, k, arma::fill::zeros);    // [eq. 22.4.5]
  arma::mat xi_t_t_tmp(Tsize, M, arma::fill::zeros);    
  arma::mat xi_tp1_t(Tsize, k, arma::fill::zeros);  // [eq. 22.4.6]
  arma::mat xi_tp1_t_tmp(Tsize, M, arma::fill::zeros);  
  double pi = arma::datum::pi;
  arma::vec xi_t_tm1_AR = pinfAR;
  arma::vec xi_t_tm1 = pinf;
  for (int xt = 0; xt<Tsize; xt++){
    for (int xm = 0; xm<M; xm++){
      arma::mat eps_m = eps[xm];
      arma::mat sigma_m = sigAR[xm];
      //eta(xt,xm) = as_scalar((1/(sqrt(pow(2*pi,q)*det(sigma_m))))*exp(-0.5*(eps_m.row(xt)*inv(sigma_m)*trans(eps_m.row(xt)))));
      eta(xt,xm) = as_scalar((1/(sqrt(pow(2*pi,q)*det(sigma_m))))*exp(-0.5*(eps_m.row(xt)*solve(sigma_m,trans(eps_m.row(xt))))));
    }
    arma::vec xi_eta = xi_t_tm1_AR%trans(eta.row(xt));
    f_t.row(xt) = sum(xi_eta);
    double f_tmp = sum(xi_eta);
    xi_t_t_tmp.row(xt) = trans(xi_eta/f_tmp);
    for (int xk = 1; xk<=k; xk++){
      arma::vec xi_t_t_row = trans(xi_t_t_tmp.row(xt));
      xi_t_t.submat(xt,xk-1,xt,xk-1) = sum(xi_t_t_row.rows(find(state_ind==xk)));
    }
    xi_t_tm1_AR = PAR*trans(xi_t_t_tmp.row(xt));
    xi_t_tm1 = P*trans(xi_t_t.row(xt));
    xi_tp1_t_tmp.row(xt) = trans(xi_t_tm1_AR);
    xi_tp1_t.row(xt) = trans(xi_t_tm1);
  }
  // ----- log-Likelihood
  arma::vec logf = log(f_t);
  double logL = sum(logf);     //[eq. 22.4.7]
  return(logL);
}

// ==============================================================================
//' @title Markov-switching vector autoregressive log-likelihood objective function (minimization version)
//' 
//' @description This function computes the (negative) log-likelihood for a markov-switching vector autoregressive model
//' 
//' @param theta vector of model parameters
//' @param mdl List with model attributes
//' @param k integer determining the number of regimes
//' 
//' @return negative log-likelihood given data
//' 
//' @export
// [[Rcpp::export]]
double logLike_MSVARmdl_min(arma::vec theta, List mdl, int k){
  double logL_negative = -logLike_MSVARmdl(theta, mdl, k);
  return(logL_negative);
}


// ==============================================================================
//' @title Hidden Markov model log-likelihood function 
//' 
//' @description This function computes the log-likelihood for a Hidden MArkov model and uses the Hamilton smoother to obtain smoothed probabilities of each state. This is also the expectation step in the Expectation Maximization algorithm for a Markov-switching autoregressive model.
//' 
//' @param theta vector of model parameters
//' @param mdl List with model attributes
//' @param k integer determining the number of regimes
//' 
//' @return log-likelihood given data
//'  
//' @return List with log-likelihood (log-likelihood), and smoothed probabilities of each regime (xi_t_T).
//' 
//' @export
// [[Rcpp::export]]
List ExpectationM_HMmdl(arma::vec theta, List mdl, int k){
  // ---------- Initialize parameters
  arma::mat y = mdl["y"];
  int q = y.n_cols;
  int Tsize = y.n_rows;
  bool msmu = mdl["msmu"];
  bool msvar = mdl["msvar"];
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
  // ----- Transition probabilities
  int th_len = theta.n_elem;
  arma::mat P = reshape(theta.subvec(th_len - k*k, th_len - 1), k, k);
  // Regime limiting probabilities
  arma::vec pinf = limP(P);
  // ----- Compute Residuals
  List eps(k); 
  arma::mat repmu(Tsize, 1, arma::fill::ones);
  for (int xk = 0; xk<k; xk++){
    arma::mat eps_tmp = y - repmu*mu_k.row(xk);
    eps[xk] =  eps_tmp;
  } 
  // ----- Begin Calculating log-likelihood & filtered probabilities
  arma::mat eta(Tsize, k, arma::fill::zeros);       // [eq. 22.4.2]
  arma::mat f_t(Tsize, 1, arma::fill::zeros);       // [eq. 22.4.8]
  arma::mat xi_t_t(Tsize, k, arma::fill::zeros);    // [eq. 22.4.5]
  arma::mat xi_tp1_t(Tsize, k, arma::fill::zeros);  // [eq. 22.4.6]
  double pi = arma::datum::pi;
  arma::vec xi_t_tm1 = pinf;
  for (int xt = 0; xt<Tsize; xt++){
    for (int xk = 0; xk<k; xk++){
      arma::mat eps_k = eps[xk];
      arma::mat sigma_k = sigma[xk];
      //eta(xt,xk) = as_scalar((1/(sqrt(pow(2*pi,q)*det(sigma_k))))*exp(-0.5*(eps_k.row(xt)*inv(sigma_k)*trans(eps_k.row(xt)))));
      eta(xt,xk) = as_scalar((1/(sqrt(pow(2*pi,q)*det(sigma_k))))*exp(-0.5*(eps_k.row(xt)*solve(sigma_k,trans(eps_k.row(xt))))));
    }
    arma::vec xi_eta = xi_t_tm1%trans(eta.row(xt));
    f_t.row(xt) = sum(xi_eta);
    double f_tmp = sum(xi_eta);
    xi_t_t.row(xt) = trans(xi_eta/f_tmp);
    xi_t_tm1 = P*trans(xi_t_t.row(xt));
    xi_tp1_t.row(xt) = trans(xi_t_tm1);
  }
  // ----- log-Likelihood
  arma::vec logf = log(f_t);
  double logL = sum(logf);     //[eq. 22.4.7]
  double L = logL/Tsize; 
  // ---------- Hamilton smoother
  arma::mat xi_t_T(Tsize, k, arma::fill::zeros); // [eq. 22.4.14]
  xi_t_T.row(Tsize-1)  = xi_t_t.row(Tsize-1);
  for (int xT = Tsize-2; xT>=0; xT--){
    xi_t_T.row(xT) = xi_t_t.row(xT)%trans(trans(P)*trans(xi_t_T.row(xT+1)/xi_tp1_t.row(xT)));
  }
  // ----- Compute model residuals 
  arma::mat residuals(Tsize, q, arma::fill::zeros); 
  arma::mat repq(1, q, arma::fill::ones);
  for (int xk = 0; xk<k;xk++){
    arma::mat eps_tmp = eps[xk];
    residuals = residuals + eps_tmp%(xi_t_T.col(xk)*repq);
  }
  // ----- Organize output
  List MSloglik_output;
  MSloglik_output["logLike"] = logL;
  MSloglik_output["L"] = L;
  MSloglik_output["eta"] = eta;
  MSloglik_output["f_t"] = f_t;
  MSloglik_output["xi_t_t"] = xi_t_t;
  MSloglik_output["xi_tp1_t"] = xi_tp1_t;
  MSloglik_output["xi_t_T"] = xi_t_T;
  MSloglik_output["P"] = P;
  MSloglik_output["pinf"] = pinf;
  MSloglik_output["resid"] = eps;
  MSloglik_output["mu"] = mu;
  MSloglik_output["sigma"] = sig;
  MSloglik_output["theta"] = theta;
  MSloglik_output["residuals"] = eps;
  MSloglik_output["resid"] = residuals;
  return(MSloglik_output);
}


// ==============================================================================
//' @title Markov-switching autoregressive log-likelihood function 
//' 
//' @description This function computes the log-likelihood for a markov-switching autoregressive model and uses the Hamilton smoother to obtain smoothed probabilities of each state. This is also the expectation step in the Expectation Maximization algorithm for a Markov-switching autoregressive model.
//' 
//' @param theta vector of model parameters
//' @param mdl List with model attributes
//' @param k integer determining the number of regimes
//' 
//' @return log-likelihood given data
//'  
//' @return List with log-likelihood (log-likelihood), and smoothed probabilities of each regime (xi_t_T).
//' 
//' @export
// [[Rcpp::export]]
List ExpectationM_MSARmdl(arma::vec theta, List mdl, int k){
  // ---------- Initialize parameters
  arma::vec y = mdl["y"];
  int Tsize = y.n_elem;
  int p = mdl["p"];
  bool msmu = mdl["msmu"];
  bool msvar = mdl["msvar"];
  int M = pow(k, p+1);
  List pList = paramList_MSARmdl(theta, p, k, msmu, msvar);
  arma::mat mu = pList["mu"];
  arma::mat sig = pList["sig"];
  arma::mat P = pList["P"];
  arma::vec pinf = pList["pinf"];
  arma::mat muAR = pList["muAR"];
  arma::mat sigAR = pList["sigAR"];
  arma::vec phi = pList["phi"];
  arma::mat PAR = pList["P_AR"];
  arma::vec pinfAR = pList["pinf_AR"];
  arma::vec state_ind = pList["state_ind"];
  // ----- Compute residuals in each regime M
  List mdl_tmp  = clone(mdl);
  mdl_tmp["phi"] = phi;
  arma::mat eps = calcResid_MSARmdl(mdl_tmp, muAR, k);
  // ----- Begin Calculating log-likelihood & filtered probabilities
  arma::mat eta(Tsize, M, arma::fill::zeros);       // [eq. 22.4.2]
  arma::mat f_t(Tsize, 1, arma::fill::zeros);       // [eq. 22.4.8]
  arma::mat xi_t_t(Tsize, k, arma::fill::zeros);    // [eq. 22.4.5]
  arma::mat xi_t_t_tmp(Tsize, M, arma::fill::zeros);    
  arma::mat xi_tp1_t(Tsize, k, arma::fill::zeros);  // [eq. 22.4.6]
  arma::mat xi_tp1_t_tmp(Tsize, M, arma::fill::zeros);  
  double pi = arma::datum::pi;
  arma::vec xi_t_tm1_AR = pinfAR;
  arma::vec xi_t_tm1 = pinf;
  for (int xt = 0; xt<Tsize; xt++){
    eta.row(xt) = trans(1/sqrt(2*pi*(sigAR)))%exp((pow(eps.row(xt),2)%trans(-1/(2*sigAR))));
    arma::vec xi_eta = xi_t_tm1_AR%trans(eta.row(xt));
    f_t.row(xt) = sum(xi_eta);
    double f_tmp = sum(xi_eta);
    xi_t_t_tmp.row(xt) = trans(xi_eta/f_tmp);
    for (int xk = 1; xk<=k; xk++){
      arma::vec xi_t_t_row = trans(xi_t_t_tmp.row(xt));
      xi_t_t.submat(xt,xk-1,xt,xk-1) = sum(xi_t_t_row.rows(find(state_ind==xk)));
    }
    xi_t_tm1_AR = PAR*trans(xi_t_t_tmp.row(xt));
    xi_t_tm1 = P*trans(xi_t_t.row(xt));
    xi_tp1_t_tmp.row(xt) = trans(xi_t_tm1_AR);
    xi_tp1_t.row(xt) = trans(xi_t_tm1);
  }
  // ----- log-Likelihood
  arma::vec logf = log(f_t);
  double logL = sum(logf);     //[eq. 22.4.7]
  double L = logL/Tsize; 
  // ---------- Hamilton smoother
  arma::mat xi_t_T(Tsize, k, arma::fill::zeros); // [eq. 22.4.14]
  arma::mat xi_t_T_tmp(Tsize, M, arma::fill::zeros); 
  xi_t_T_tmp.row(Tsize-1)  = xi_t_t_tmp.row(Tsize-1);
  xi_t_T.row(Tsize-1)  = xi_t_t.row(Tsize-1);
  for (int xT = Tsize-2; xT>=0; xT--){
    xi_t_T_tmp.row(xT) = xi_t_t_tmp.row(xT)%trans(trans(PAR)*trans(xi_t_T_tmp.row(xT+1)/xi_tp1_t_tmp.row(xT)));
    xi_t_T.row(xT) = xi_t_t.row(xT)%trans(trans(P)*trans(xi_t_T.row(xT+1)/xi_tp1_t.row(xT)));
  }
  // ----- Compute model residuals 
  arma::mat residuals_tmp = eps%xi_t_T_tmp; 
  arma::mat residuals(Tsize, k);
  for (int xk = 1; xk<=k;xk++){
    residuals.col(xk-1) = arma::sum(residuals_tmp.cols(find(state_ind==xk)), 1);
  }
  residuals = arma::sum(residuals,1);
  // ----- Organize output
  List MSloglik_output;
  MSloglik_output["logLike"] = logL;
  MSloglik_output["L"] = L;
  MSloglik_output["eta"] = eta;
  MSloglik_output["f_t"] = f_t;
  MSloglik_output["xi_t_t_AR"] = xi_t_t_tmp;
  MSloglik_output["xi_tp1_t_AR"] = xi_tp1_t_tmp;
  MSloglik_output["xi_t_T_AR"] = xi_t_T_tmp;
  MSloglik_output["xi_t_t"] = xi_t_t;
  MSloglik_output["xi_tp1_t"] = xi_tp1_t;
  MSloglik_output["xi_t_T"] = xi_t_T;
  MSloglik_output["PAR"] = PAR;
  MSloglik_output["P"] = P;
  MSloglik_output["pinfAR"] = pinfAR;
  MSloglik_output["pinf"] = pinf;
  MSloglik_output["muAR"] = muAR;
  MSloglik_output["sigmaAR"] = sigAR;
  MSloglik_output["mu"] = mu;
  MSloglik_output["sigma"] = sig;
  MSloglik_output["phi"] = phi;
  MSloglik_output["theta"] = theta;
  MSloglik_output["state_ind"] = state_ind;
  MSloglik_output["residuals"] = eps;
  MSloglik_output["resid"] = residuals;
  return(MSloglik_output);
}

// ==============================================================================
//' @title Markov-switching vector autoregressive log-likelihood function
//' 
//' @description This function computes the log-likelihood for a markov-switching vector autoregressive model and uses the Hamilton smoother to obtain smoothed probabilities of each state. This is also the expectation step in the Expectation Maximization algorithm for a Markov-switching autoregressive model.
//' 
//' @param theta vector of model parameters
//' @param mdl List with model attributes
//' @param k integer determining the number of regimes
//' 
//' @return log-likelihood given data
//'  
//' @return List with log-likelihood (log-likelihood), and smoothed probabilities of each regime (xi_t_T).
//' 
//' @export
// [[Rcpp::export]]
List ExpectationM_MSVARmdl(arma::vec theta, List mdl, int k){
  // ---------- Initialize parameters
  arma::mat y = mdl["y"];
  int q = y.n_cols;
  int Tsize = y.n_rows;
  int ar = mdl["p"];
  bool msmu = mdl["msmu"];
  bool msvar = mdl["msvar"];
  int M = pow(k, ar+1);
  List pList = paramList_MSVARmdl(theta, q, ar, k, msmu, msvar);
  arma::mat mu = pList["mu"];
  List sig = pList["sigma"];
  arma::mat P = pList["P"];
  arma::vec pinf = pList["pinf"];
  List muAR = pList["muAR"];
  List sigAR = pList["sigAR"];
  arma::mat phi = pList["phi"];
  arma::mat PAR = pList["P_AR"];
  arma::vec pinfAR = pList["pinf_AR"];
  arma::vec state_ind = pList["state_ind"];
  // ----- Compute residuals in each regime
  List mdl_tmp  = clone(mdl);
  mdl_tmp["phi"] = phi;
  List eps = calcResid_MSVARmdl(mdl_tmp, muAR, k);
  // ----- Begin Calculating log-likelihood & filtered probabilities
  arma::mat eta(Tsize, M, arma::fill::zeros);       // [eq. 22.4.2]
  arma::mat f_t(Tsize, 1, arma::fill::zeros);       // [eq. 22.4.8]
  arma::mat xi_t_t(Tsize, k, arma::fill::zeros);    // [eq. 22.4.5]
  arma::mat xi_t_t_tmp(Tsize, M, arma::fill::zeros);    
  arma::mat xi_tp1_t(Tsize, k, arma::fill::zeros);  // [eq. 22.4.6]
  arma::mat xi_tp1_t_tmp(Tsize, M, arma::fill::zeros);  
  double pi = arma::datum::pi;
  arma::vec xi_t_tm1_AR = pinfAR;
  arma::vec xi_t_tm1 = pinf;
  for (int xt = 0; xt<Tsize; xt++){
    for (int xm = 0; xm<M; xm++){
      arma::mat eps_m = eps[xm];
      arma::mat sigma_m = sigAR[xm];
      //eta(xt,xm) = as_scalar((1/(sqrt(pow(2*pi,q)*det(sigma_m))))*exp(-0.5*(eps_m.row(xt)*inv(sigma_m)*trans(eps_m.row(xt)))));
      eta(xt,xm) = as_scalar((1/(sqrt(pow(2*pi,q)*det(sigma_m))))*exp(-0.5*(eps_m.row(xt)*solve(sigma_m,trans(eps_m.row(xt))))));
    }
    arma::vec xi_eta = xi_t_tm1_AR%trans(eta.row(xt));
    f_t.row(xt) = sum(xi_eta);
    double f_tmp = sum(xi_eta);
    xi_t_t_tmp.row(xt) = trans(xi_eta/f_tmp);
    for (int xk = 1; xk<=k; xk++){
      arma::vec xi_t_t_row = trans(xi_t_t_tmp.row(xt));
      xi_t_t.submat(xt,xk-1,xt,xk-1) = sum(xi_t_t_row.rows(find(state_ind==xk)));
    }
    xi_t_tm1_AR = PAR*trans(xi_t_t_tmp.row(xt));
    xi_t_tm1 = P*trans(xi_t_t.row(xt));
    xi_tp1_t_tmp.row(xt) = trans(xi_t_tm1_AR);
    xi_tp1_t.row(xt) = trans(xi_t_tm1);
  }
  // ----- log-Likelihood
  arma::vec logf = log(f_t);
  double logL = sum(logf);     //[eq. 22.4.7]
  double L = logL/Tsize; 
  // ---------- Hamilton smoother
  arma::mat xi_t_T(Tsize, k, arma::fill::zeros); // [eq. 22.4.14]
  arma::mat xi_t_T_tmp(Tsize, M, arma::fill::zeros); 
  xi_t_T_tmp.row(Tsize-1)  = xi_t_t_tmp.row(Tsize-1);
  xi_t_T.row(Tsize-1)  = xi_t_t.row(Tsize-1);
  for (int xT = Tsize-2; xT>=0; xT--){
    xi_t_T_tmp.row(xT) = xi_t_t_tmp.row(xT)%trans(trans(PAR)*trans(xi_t_T_tmp.row(xT+1)/xi_tp1_t_tmp.row(xT)));
    xi_t_T.row(xT) = xi_t_t.row(xT)%trans(trans(P)*trans(xi_t_T.row(xT+1)/xi_tp1_t.row(xT)));
  }
  // ----- Compute model residuals 
  arma::mat residuals(Tsize, q, arma::fill::zeros); 
  arma::mat repq(1, q, arma::fill::ones);
  for (int xk = 0; xk<M;xk++){
    arma::mat eps_tmp = eps[xk];
    residuals = residuals + eps_tmp%(xi_t_T_tmp.col(xk)*repq);
  }
  // ----- Organize output
  List MSloglik_output;
  MSloglik_output["logLike"] = logL;
  MSloglik_output["L"] = L;
  MSloglik_output["eta"] = eta;
  MSloglik_output["f_t"] = f_t;
  MSloglik_output["xi_t_t_AR"] = xi_t_t_tmp;
  MSloglik_output["xi_tp1_t_AR"] = xi_tp1_t_tmp;
  MSloglik_output["xi_t_T_AR"] = xi_t_T_tmp;
  MSloglik_output["xi_t_t"] = xi_t_t;
  MSloglik_output["xi_tp1_t"] = xi_tp1_t;
  MSloglik_output["xi_t_T"] = xi_t_T;
  MSloglik_output["PAR"] = PAR;
  MSloglik_output["P"] = P;
  MSloglik_output["pinfAR"] = pinfAR;
  MSloglik_output["pinf"] = pinf;
  MSloglik_output["muAR"] = muAR;
  MSloglik_output["sigmaAR"] = sigAR;
  MSloglik_output["mu"] = mu;
  MSloglik_output["sigma"] = sig;
  MSloglik_output["phi"] = phi;
  MSloglik_output["theta"] = theta;
  MSloglik_output["state_ind"] = state_ind;
  MSloglik_output["residuals"] = eps;
  MSloglik_output["resid"] = residuals;
  return(MSloglik_output);
}


// ==============================================================================
//' @title Maximization step of EM algorithm for Hidden Markov model
//' 
//' @description This function performs the maximization step of the Expectation Maximization algorithm for Hidden Markov models.
//' 
//' @param theta vector of model parameters
//' @param mdl List with model attributes
//' @param MSloglik_output List with output from Expectation step of EM algorithm for Markov-switching vector autoregressive model
//' @param k integer determining the number of regimes
//'  
//' @return List with new maximized parameters
//' 
//' @export
// [[Rcpp::export]]
List EMaximization_HMmdl(arma::vec theta, List mdl, List MSloglik_output, int k){
  Rcpp::Environment mstest("package:MSTest");
  Rcpp::Function argrid_MSVARmdl = mstest["argrid_MSVARmdl"];
  // ---------- Initialize parameters
  arma::mat y = mdl["y"];
  int q = mdl["q"];
  bool msmu = mdl["msmu"];
  bool msvar = mdl["msvar"];
  // ----- Output from log-Likelihood calculation
  arma::mat P = MSloglik_output["P"];
  arma::mat xi_t_T = MSloglik_output["xi_t_T"];
  arma::mat xi_t_t = MSloglik_output["xi_t_t"];
  arma::mat xi_tp1_t = MSloglik_output["xi_tp1_t"];
  List sig = MSloglik_output["sigma"];
  // Length of Time-Series (T) and number of variables
  int Tsize = y.n_rows;
  arma::mat repq(1,q, arma::fill::ones);
  // --------------- Update transition matrix (P) and limiting probabilities (pinf)
  // ----- Estimate Transition probs for k-state Markov-chain [eq. 22.4.16]
  arma::mat xi_t_T_tmp = xi_t_T.rows(1,Tsize-1);
  arma::mat xi_t_t_tmp = xi_t_t.rows(0,Tsize-2);
  arma::mat xi_tp1_t_tmp = xi_tp1_t.rows(0,Tsize-2);
  arma::mat p_ij(Tsize-1, k*k, arma::fill::zeros);
  for (int xk = 0; xk<k; xk++){
    p_ij.submat(0,xk*k,Tsize-2,xk*k+k-1) = (xi_t_T_tmp%trans(P.col(xk)*trans(xi_t_t_tmp.col(xk))))/xi_tp1_t_tmp;
  }
  arma::mat P_new(k, k, arma::fill::zeros);
  arma::mat p_ij_sums = arma::sum(p_ij,0);
  for (int xk = 0 ; xk<k; xk++){
    arma::vec regime_prob = xi_t_T.submat(0,xk,(Tsize-2),xk);
    double regimesum = sum(regime_prob);
    P_new.row(xk) = p_ij_sums.submat(0,(xk*k),0,(xk*k+(k-1)))/regimesum;
  }
  P_new = trans(P_new);
  // Get limiting probabilities implied by transition matrix P_0
  arma::vec pinf = limP(P_new);
  // --------------- Update estimates for mu
  // ----- Obtain state indicators
  arma::mat mu(1+msmu*(k-1), q, arma::fill::zeros);
  arma::mat mu_k(k, q, arma::fill::zeros);
  // ----- Compute updated mu
  if (msmu == TRUE){
    for (int xk = 0 ; xk<k; xk++){
      mu.row(xk) = arma::sum(y%(xi_t_T.col(xk)*repq),0)/sum(xi_t_T.col(xk));
    }
    mu_k = mu;
  }else{
    for (int xk = 0 ; xk<k; xk++){
      mu = mu + arma::sum(y%(xi_t_T.col(xk)*repq),0);
    }
    mu = mu/Tsize;
    for (int xk = 0 ; xk<k; xk++){
      mu_k.row(xk) = mu;
    }
  }
  // --------------- Update estimates for variance
  // ----- Compute Residuals
  List eps(k); 
  arma::mat repmu(Tsize, 1, arma::fill::ones);
  for (int xk = 0; xk<k; xk++){
    arma::mat eps_tmp = y - repmu*mu_k.row(xk);
    eps[xk] =  eps_tmp;
  } 
  // // ----- Compute updated sigma
  List sigma(k);
  arma::vec sigma_out;
  if (msvar == TRUE){
    for (int xk = 0 ; xk<k; xk++){
      arma::mat U = eps[xk];
      arma::mat sigma_m_tmp = (trans(U)*diagmat(xi_t_T.col(xk))*U)/sum(xi_t_T.col(xk));
      sigma[xk] = sigma_m_tmp; 
      sigma_out = join_vert(sigma_out, covar_vech(sigma_m_tmp));
    }
  }else{
    arma::mat sigma_m_tmp(q,q,arma::fill::zeros);
    for (int xk = 0 ; xk<k; xk++){
      arma::mat U = eps[xk];
      sigma_m_tmp  = sigma_m_tmp + trans(U)*diagmat(xi_t_T.col(xk))*U;
    }
    sigma_m_tmp = sigma_m_tmp/Tsize;
    for (int xk = 0; xk<k;xk++){
      sigma[xk] = sigma_m_tmp;
    }
    sigma_out = covar_vech(sigma_m_tmp);
  }
  // ----- Compute model residuals 
  arma::mat residuals(Tsize, q, arma::fill::zeros); 
  for (int xk = 0; xk<k;xk++){
    arma::mat eps_tmp = eps[xk];
    residuals = residuals + eps_tmp%(xi_t_T.col(xk)*repq);
  }
  // --------------- Organize output
  // ----- Produce new theta vector
  arma::vec theta_new = join_vert(vectorise(trans(mu)), sigma_out);
  theta_new = join_vert(theta_new, vectorise(P_new));
  // ----- Fill output List
  List Maxim_output;
  Maxim_output["theta"] = theta_new;
  Maxim_output["P"] = P_new;
  Maxim_output["pinf"] = pinf;
  Maxim_output["mu"] = mu;
  Maxim_output["sigma"] = sigma;
  Maxim_output["residuals"] = eps;
  Maxim_output["resid"] = residuals;
  return(Maxim_output);
}


// ==============================================================================
//' @title Maximization step of EM algorithm for Markov-switching autoregressive model
//' 
//' @description This function performs the maximization step of the Expectation Maximization algorithm for Markov-switching autoregressive model.
//' 
//' @param theta vector of model parameters
//' @param mdl List with model attributes
//' @param MSloglik_output List with output from Expectation step of EM algorithm for Markov-switching autoregressive model
//' @param k integer determining the number of regimes
//' 
//' @return List with new maximized parameters
//' 
//' @export
// [[Rcpp::export]]
List EMaximization_MSARmdl(arma::vec theta, List mdl, List MSloglik_output, int k){
  Rcpp::Environment mstest("package:MSTest");
  Rcpp::Function argrid_MSARmdl = mstest["argrid_MSARmdl"];
  // ---------- Initialize parameters
  arma::vec y = mdl["y"];
  int ar = mdl["p"];
  bool msmu = mdl["msmu"];
  bool msvar = mdl["msvar"];
  // ----- Output from log-Likelihood calculation
  arma::mat P = MSloglik_output["P"];
  arma::mat xi_t_T = MSloglik_output["xi_t_T"];
  arma::mat xi_t_t = MSloglik_output["xi_t_t"];
  arma::mat xi_tp1_t = MSloglik_output["xi_tp1_t"];
  arma::mat xi_t_T_AR = MSloglik_output["xi_t_T_AR"];
  arma::mat xi_t_t_AR = MSloglik_output["xi_t_t_AR"];
  arma::mat xi_tp1_t_AR = MSloglik_output["xi_tp1_t_AR"];
  arma::vec state_ind = MSloglik_output["state_ind"];
  arma::vec phi = MSloglik_output["phi"];
  arma::vec sig = MSloglik_output["sigma"];
  // Length of Time-Series (T)
  int Tsize = y.n_elem;
  // Number of regimes consistent with AR lags
  int M = pow(k, ar+1);
  // --------------- Update transition matrix (P) and limiting probabilities (pinf)
  // ----- Estimate Transition probs for k-state Markov-chain [eq. 22.4.16]
  arma::mat xi_t_T_tmp = xi_t_T.rows(1,Tsize-1);
  arma::mat xi_t_t_tmp = xi_t_t.rows(0,Tsize-2);
  arma::mat xi_tp1_t_tmp = xi_tp1_t.rows(0,Tsize-2);
  arma::mat p_ij(Tsize-1, k*k, arma::fill::zeros);
  for (int xk = 0; xk<k; xk++){
    p_ij.submat(0,xk*k,Tsize-2,xk*k+k-1) = (xi_t_T_tmp%trans(P.col(xk)*trans(xi_t_t_tmp.col(xk))))/xi_tp1_t_tmp;
  }
  arma::mat P_new(k, k, arma::fill::zeros);
  arma::mat p_ij_sums = arma::sum(p_ij,0);
  for (int xk = 0 ; xk<k; xk++){
    arma::vec regime_prob = xi_t_T.submat(0,xk,(Tsize-2),xk);
    double regimesum = sum(regime_prob);
    P_new.row(xk) = p_ij_sums.submat(0,(xk*k),0,(xk*k+(k-1)))/regimesum;
  }
  P_new = trans(P_new);
  // Get limiting probabilities implied by transition matrix P_0
  arma::vec pinf = limP(P_new);
  // --------------- Update estimates for mu
  // ----- Obtain state indicators
  arma::vec mu(1+msmu*(k-1), arma::fill::zeros);
  arma::vec mu_tmp(M, arma::fill::zeros);
  arma::vec sum_tmp(M, arma::fill::zeros);
  // ----- Compute updated mu
  if (msmu == TRUE){
    for (int xk = 0 ; xk<M; xk++){
      mu_tmp(xk) = sum(y%xi_t_T_AR.col(xk));
      sum_tmp(xk) = sum(xi_t_T_AR.col(xk));
    }
    for (int xk = 1; xk<=k;xk++){
      mu(xk-1) = as_scalar(sum(mu_tmp.rows(find(state_ind==xk)))/sum(sum_tmp.rows(find(state_ind==xk))));
    }
  }else{
    for (int xk = 0 ; xk<M; xk++){
      mu = mu + sum(y%xi_t_T_AR.col(xk));
    }
    mu = mu/Tsize;
  }
  // ----- Update estimates for phi
  arma::mat muAR(M, ar+1,arma::fill::zeros);
  arma::vec sigAR(M, 1,arma::fill::zeros);
  List mugrid = argrid_MSARmdl(mu, sig, k, ar, msmu, msvar);
  muAR = as<arma::mat>(mugrid["mu"]);
  sigAR = as<arma::vec>(mugrid["sig"]);
  arma::mat x = mdl["x"];
  arma::vec repmu(Tsize,arma::fill::ones);
  arma::mat denom(ar,ar,arma::fill::zeros);
  arma::mat num(ar,1,arma::fill::zeros);
  for (int xm = 0; xm<M; xm++){
    arma::vec y_tmp = y - repmu*muAR(xm,0);
    arma::mat x_tmp = x - repmu*muAR.submat(xm,1,xm,ar);
    denom = denom + (trans(x_tmp)*diagmat(xi_t_T_AR.col(xm)/sigAR(xm))*x_tmp);
    num = num + (trans(x_tmp)*diagmat(xi_t_T_AR.col(xm)/sigAR(xm))*y_tmp);
  }
  arma::vec phi_new = inv(denom)*num;
  // ----- Update estimates for variance
  // ----- Compute new residuals
  List mdl_tmp = clone(mdl);
  mdl_tmp["phi"] = phi_new;
  arma::mat eps = calcResid_MSARmdl(mdl_tmp, muAR, k);
  // ----- Compute updated sigma
  arma::vec sigma(1+msvar*(k-1), arma::fill::zeros);
  arma::vec sigma_tmp(M, arma::fill::zeros);
  if (msvar == TRUE){
    for (int xk = 0 ; xk<M; xk++){
      sigma_tmp(xk) = sum(eps.col(xk)%eps.col(xk)%xi_t_T_AR.col(xk)); 
      sum_tmp(xk) = sum(xi_t_T_AR.col(xk)); 
    }
    for (int xk = 1; xk<=k;xk++){
      sigma(xk-1) = as_scalar(sum(sigma_tmp.rows(find(state_ind==xk)))/sum(sum_tmp.rows(find(state_ind==xk))));
    }
  }else{
    sigma = arma::sum(arma::sum(eps%eps%xi_t_T_AR,0),1)/Tsize;
  }
  // ----- Compute model residuals 
  arma::mat residuals(Tsize, k);
  arma::mat residuals_tmp = eps%xi_t_T_AR; 
  for (int xk = 1; xk<=k;xk++){
    residuals.col(xk-1) = arma::sum(residuals_tmp.cols(find(state_ind==xk)), 1);
  }
  residuals = arma::sum(residuals,1);
  // --------------- Organize output
  // ----- Produce new theta vector
  arma::vec theta_new = join_vert(mu, sigma);
  theta_new = join_vert(theta_new, phi_new);
  theta_new = join_vert(theta_new, vectorise(P_new));
  // ----- Fill output List
  List Maxim_output;
  Maxim_output["theta"] = theta_new;
  Maxim_output["P"] = P_new;
  Maxim_output["pinf"] = pinf;
  Maxim_output["mu"] = mu;
  Maxim_output["sigma"] = sigma;
  Maxim_output["phi"] = phi_new;
  Maxim_output["residuals"] = eps;
  Maxim_output["resid"] = residuals;
  return(Maxim_output);
}

// ==============================================================================
//' @title Maximization step of EM algorithm for Markov-switching vector autoregressive model
//' 
//' @description This function performs the maximization step of the Expectation Maximization algorithm for Markov-switching vector autoregressive model.
//' 
//' @param theta vector of model parameters
//' @param mdl List with model attributes
//' @param MSloglik_output List with output from Expectation step of EM algorithm for Markov-switching vector autoregressive model
//' @param k integer determining the number of regimes
//'  
//' @return List with new maximized parameters
//' 
//' @export
// [[Rcpp::export]]
List EMaximization_MSVARmdl(arma::vec theta, List mdl, List MSloglik_output, int k){
  Rcpp::Environment mstest("package:MSTest");
  Rcpp::Function argrid_MSVARmdl = mstest["argrid_MSVARmdl"];
  // ---------- Initialize parameters
  arma::mat y = mdl["y"];
  int ar = mdl["p"];
  int q = mdl["q"];
  bool msmu = mdl["msmu"];
  bool msvar = mdl["msvar"];
  // ----- Output from log-Likelihood calculation
  arma::mat P = MSloglik_output["P"];
  arma::mat xi_t_T = MSloglik_output["xi_t_T"];
  arma::mat xi_t_t = MSloglik_output["xi_t_t"];
  arma::mat xi_tp1_t = MSloglik_output["xi_tp1_t"];
  arma::mat xi_t_T_AR = MSloglik_output["xi_t_T_AR"];
  arma::mat xi_t_t_AR = MSloglik_output["xi_t_t_AR"];
  arma::mat xi_tp1_t_AR = MSloglik_output["xi_tp1_t_AR"];
  arma::vec state_ind = MSloglik_output["state_ind"];
  arma::mat phi = MSloglik_output["phi"];
  List sig = MSloglik_output["sigma"];
  // Length of Time-Series (T) and number of variables
  int Tsize = y.n_rows;
  arma::mat repq(1,q, arma::fill::ones);
  // Number of regimes consistent with AR lags
  int M = pow(k, ar+1);
  // --------------- Update transition matrix (P) and limiting probabilities (pinf)
  // ----- Estimate Transition probs for k-state Markov-chain [eq. 22.4.16]
  arma::mat xi_t_T_tmp = xi_t_T.rows(1,Tsize-1);
  arma::mat xi_t_t_tmp = xi_t_t.rows(0,Tsize-2);
  arma::mat xi_tp1_t_tmp = xi_tp1_t.rows(0,Tsize-2);
  arma::mat p_ij(Tsize-1, k*k, arma::fill::zeros);
  for (int xk = 0; xk<k; xk++){
    p_ij.submat(0,xk*k,Tsize-2,xk*k+k-1) = (xi_t_T_tmp%trans(P.col(xk)*trans(xi_t_t_tmp.col(xk))))/xi_tp1_t_tmp;
  }
  arma::mat P_new(k, k, arma::fill::zeros);
  arma::mat p_ij_sums = arma::sum(p_ij,0);
  for (int xk = 0 ; xk<k; xk++){
    arma::vec regime_prob = xi_t_T.submat(0,xk,(Tsize-2),xk);
    double regimesum = sum(regime_prob);
    P_new.row(xk) = p_ij_sums.submat(0,(xk*k),0,(xk*k+(k-1)))/regimesum;
  }
  P_new = trans(P_new);
  // Get limiting probabilities implied by transition matrix P_0
  arma::vec pinf = limP(P_new);
  // --------------- Update estimates for mu
  // ----- Obtain state indicators
  arma::mat mu(1+msmu*(k-1), q, arma::fill::zeros);
  arma::mat mu_k(k, q, arma::fill::zeros);
  arma::mat mu_tmp(M, q, arma::fill::zeros);
  arma::vec sum_tmp(M, arma::fill::zeros);
  // ----- Compute updated mu
  if (msmu == TRUE){
    for (int xk = 0 ; xk<M; xk++){
      mu_tmp.row(xk) = arma::sum(y%(xi_t_T_AR.col(xk)*repq),0);
      sum_tmp(xk) = sum(xi_t_T_AR.col(xk));
    }
    for (int xk = 1; xk<=k;xk++){
      mu.row(xk-1) = arma::sum(mu_tmp.rows(find(state_ind==xk)),0)/(sum(sum_tmp.rows(find(state_ind==xk)))*repq);
    }
    mu_k = mu;
  }else{
    for (int xk = 0 ; xk<M; xk++){
      mu = mu + arma::sum(y%(xi_t_T_AR.col(xk)*repq),0);
    }
    mu = mu/Tsize;
    for (int xk = 0 ; xk<k; xk++){
      mu_k.row(xk) = mu;
    }
  }
  // --------------- Update estimates for phi
  List mugrid = argrid_MSVARmdl(mu_k, sig, k, ar, msmu, msvar);
  List muAR = mugrid["mu"];
  List sigAR = mugrid["sig"];
  arma::mat x = mdl["x"];
  arma::vec repmu(Tsize, arma::fill::ones);
  arma::mat denom(q*q*ar, q*q*ar, arma::fill::zeros);
  arma::mat num(q*q*ar, 1, arma::fill::zeros);
  for (int xm = 0; xm<M; xm++){
    arma::mat muAR_tmp = muAR[xm]; 
    arma::mat y_tmp = y - repmu*trans(muAR_tmp.col(0));
    arma::mat xz_tmp(Tsize, q*ar, arma::fill::zeros); 
    for (int xp = 0; xp<ar; xp++){
      xz_tmp.submat(0,q*xp,Tsize-1,q*xp+q-1) = x.submat(0,q*xp,Tsize-1,q*xp+q-1) - repmu*trans(muAR_tmp.col(xp+1));
    }
    arma::mat sigma_m = sigAR[xm];
    denom = denom + kron(inv(sigma_m),trans(xz_tmp)*diagmat(xi_t_T_AR.col(xm))*xz_tmp);
    num = num + kron(inv(sigma_m),trans(xz_tmp)*diagmat(xi_t_T_AR.col(xm)))*vectorise(y_tmp);
  }
  arma::mat phi_new = reshape(inv(denom)*num, q*ar, q); 
  
  // Companion Mat
  arma::mat F_tmp = trans(phi_new);
  arma::mat diag_mat = arma::eye(q*(ar-1),q*(ar-1));
  arma::mat diag_zero(q*(ar-1), q, arma::fill::zeros);
  arma::mat Mn = join_rows(diag_mat,diag_zero);
  arma::mat F = join_cols(F_tmp,Mn); 
  // --------------- Update estimates for variance
  // ----- Compute new residuals
  List mdl_tmp  = clone(mdl);
  mdl_tmp["phi"] = trans(phi_new);
  List eps = calcResid_MSVARmdl(mdl_tmp, muAR, k);
  // // ----- Compute updated sigma
  List sigma(k);
  for (int xk = 0; xk<k;xk++){
    arma::mat sigma_m_tmp(q,q,arma::fill::zeros);
    sigma[xk] = sigma_m_tmp;
  }
  int sigN = (q*(q+1))/2;
  arma::vec sigma_out(sigN+sigN*msvar*(k-1), arma::fill::zeros);
  arma::vec Tsize_k(k, arma::fill::zeros);
  if (msvar == TRUE){
    for (int xk = 0 ; xk<M; xk++){
      arma::mat U = eps[xk];
      arma::mat sigma_m_tmp  = sigma[state_ind(xk)-1];
      sigma_m_tmp = sigma_m_tmp + trans(U)*diagmat(xi_t_T_AR.col(xk))*U;
      sigma[state_ind(xk)-1] = sigma_m_tmp; 
      Tsize_k(state_ind(xk)-1) = Tsize_k(state_ind(xk)-1) + sum(xi_t_T_AR.col(xk));
    }
    for (int xk = 0; xk<k; xk++){
      arma::mat sigma_m_tmp  = sigma[xk];
      sigma_m_tmp = sigma_m_tmp/Tsize_k(xk);
      sigma[xk] = sigma_m_tmp;
      arma::vec sigma_out_tmp = trans(sigma_m_tmp.row(0));
      for (int xn = 1; xn<q; xn++){
        sigma_out_tmp = join_vert(sigma_out_tmp,trans(sigma_m_tmp.submat(xn,xn,xn,q-1)));
      }
      sigma_out.subvec(xk*sigN,xk*sigN+sigN-1) = sigma_out_tmp;
    }
  }else{
    arma::mat sigma_m_tmp(q,q,arma::fill::zeros);
    for (int xk = 0 ; xk<M; xk++){
      arma::mat U = eps[xk];
      sigma_m_tmp  = sigma_m_tmp + trans(U)*diagmat(xi_t_T_AR.col(xk))*U;
    }
    sigma_m_tmp = sigma_m_tmp/Tsize;
    for (int xk = 0; xk<k;xk++){
      sigma[xk] = sigma_m_tmp;
    }
    arma::vec sigma_out_tmp = trans(sigma_m_tmp.row(0));
    for (int xn = 1; xn<q; xn++){
      sigma_out_tmp = join_vert(sigma_out_tmp, trans(sigma_m_tmp.submat(xn,xn,xn,q-1)));
    }
    sigma_out = sigma_out_tmp;
  }
  // ----- Compute model residuals 
  arma::mat residuals(Tsize, q, arma::fill::zeros); 
  for (int xk = 0; xk<M;xk++){
    arma::mat eps_tmp = eps[xk];
    residuals = residuals + eps_tmp%(xi_t_T_AR.col(xk)*repq);
  }
  // --------------- Organize output
  // ----- Produce new theta vector
  arma::vec theta_new = join_vert(vectorise(trans(mu)), sigma_out);
  theta_new = join_vert(theta_new,vectorise(phi_new));
  theta_new = join_vert(theta_new, vectorise(P_new));
  // ----- Fill output List
  List Maxim_output;
  Maxim_output["theta"] = theta_new;
  Maxim_output["P"] = P_new;
  Maxim_output["pinf"] = pinf;
  Maxim_output["mu"] = mu;
  Maxim_output["sigma"] = sigma;
  Maxim_output["denom"] = denom;
  Maxim_output["phi"] = trans(phi_new);
  Maxim_output["companion_mat"] = F; 
  Maxim_output["residuals"] = eps;
  Maxim_output["resid"] = residuals;
  return(Maxim_output);
}

// ==============================================================================
//' @title EM algorithm iteration for Hidden Markov model
//' 
//' @description This function performs the one iteration (E-step and M-step) of the Expectation Maximization algorithm for a Hidden Markov model.
//' 
//' @param mdl List with model attributes
//' @param EMest_output List with attributes from previous iteration
//' @param k integer determining the number of regimes
//' 
//' @return List with attributes from new iteration
//' @export
// [[Rcpp::export]]
List EMiter_HMmdl(List mdl, List EMest_output, int k){
  
  arma::vec theta = EMest_output["theta"];
  // ---------- Expectation
  List HMMloglik_output = ExpectationM_HMmdl(theta, mdl, k);
  // ---------- Maximization
  List Maxim_output = EMaximization_HMmdl(theta, mdl, HMMloglik_output, k);
  // ---------- Organize output
  arma::vec thl(3, arma::fill::zeros);
  double loglik_new = HMMloglik_output["logLike"];
  double loglik_old = EMest_output["logLike"];
  thl(0) = loglik_new;
  thl(1) = loglik_new - loglik_old;
  arma::mat P_n = Maxim_output["P"];
  arma::vec theta_n = Maxim_output["theta"];
  arma::vec delta = theta_n - theta;
  double deltath = max(abs(delta));
  thl(2) = deltath;
  List EM_output;
  EM_output["theta"] = theta_n; 
  EM_output["P"] = P_n;
  EM_output["pinf"] = Maxim_output["pinf"];
  EM_output["mu"] = Maxim_output["mu"];
  EM_output["sigma"] = Maxim_output["sigma"]; 
  EM_output["St"] = HMMloglik_output["xi_t_T"];
  EM_output["eta"] = HMMloglik_output["eta"];
  EM_output["thl"] = thl;
  EM_output["deltath"] = deltath;
  EM_output["logLike"] = loglik_new;
  EM_output["residuals"] = Maxim_output["residuals"];
  EM_output["resid"] = Maxim_output["resid"];
  return(EM_output);
}

// ==============================================================================
//' @title EM algorithm iteration for Markov-switching autoregressive model
//' 
//' @description This function performs the one iteration (E-step and M-step) of the Expectation Maximization algorithm for Markov-switching autoregressive model.
//' 
//' @param mdl List with model attributes
//' @param EMest_output List with attributes from previous iteration
//' @param k integer determining the number of regimes
//' 
//' @return List with attributes from new iteration
//' 
//' @export
// [[Rcpp::export]]
List EMiter_MSARmdl(List mdl, List EMest_output, int k){
  
  arma::vec theta = EMest_output["theta"];
  // ---------- Expectation
  List MSloglik_output = ExpectationM_MSARmdl(theta, mdl, k);
  // ---------- Maximization
  List Maxim_output = EMaximization_MSARmdl(theta, mdl, MSloglik_output, k);
  // ---------- Organize output
  arma::vec thl(3, arma::fill::zeros);
  double loglik_new = MSloglik_output["logLike"];
  double loglik_old = EMest_output["logLike"];
  thl(0) = loglik_new;
  thl(1) = loglik_new - loglik_old;
  arma::mat P_n = Maxim_output["P"];
  arma::vec theta_n = Maxim_output["theta"];
  arma::vec delta = theta_n - theta;
  double deltath = max(abs(delta));
  thl(2) = deltath;
  List EM_output;
  EM_output["theta"] = theta_n; 
  EM_output["P"] = P_n;
  EM_output["pinf"] = Maxim_output["pinf"];
  EM_output["mu"] = Maxim_output["mu"];
  EM_output["sigma"] = Maxim_output["sigma"]; 
  EM_output["phi"] = Maxim_output["phi"]; 
  EM_output["St"] = MSloglik_output["xi_t_T"];
  EM_output["eta"] = MSloglik_output["eta"];
  EM_output["thl"] = thl;
  EM_output["deltath"] = deltath;
  EM_output["logLike"] = loglik_new;
  EM_output["residuals"] = Maxim_output["residuals"];
  EM_output["resid"] = Maxim_output["resid"];
  return(EM_output);
}


// ==============================================================================
//' @title EM algorithm iteration for Markov-switching vector autoregressive model
//' 
//' @description This function performs the one iteration (E-step and M-step) of the Expectation Maximization algorithm for Markov-switching vector autoregressive model.
//' 
//' @param mdl List with model attributes
//' @param EMest_output List with attributes from previous iteration
//' @param k integer determining the number of regimes
//' 
//' @return List with attributes from new iteration
//' @export
// [[Rcpp::export]]
List EMiter_MSVARmdl(List mdl, List EMest_output, int k){
  
  arma::vec theta = EMest_output["theta"];
  // ---------- Expectation
  List MSVARloglik_output = ExpectationM_MSVARmdl(theta, mdl, k);
  // ---------- Maximization
  List Maxim_output = EMaximization_MSVARmdl(theta, mdl, MSVARloglik_output, k);
  // ---------- Organize output
  arma::vec thl(3, arma::fill::zeros);
  double loglik_new = MSVARloglik_output["logLike"];
  double loglik_old = EMest_output["logLike"];
  thl(0) = loglik_new;
  thl(1) = loglik_new - loglik_old;
  arma::mat P_n = Maxim_output["P"];
  arma::vec theta_n = Maxim_output["theta"];
  arma::vec delta = theta_n - theta;
  double deltath = max(abs(delta));
  thl(2) = deltath;
  List EM_output;
  EM_output["theta"] = theta_n; 
  EM_output["P"] = P_n;
  EM_output["pinf"] = Maxim_output["pinf"];
  EM_output["mu"] = Maxim_output["mu"];
  EM_output["sigma"] = Maxim_output["sigma"]; 
  EM_output["phi"] = Maxim_output["phi"]; 
  EM_output["companion_mat"] = Maxim_output["companion_mat"]; 
  EM_output["St"] = MSVARloglik_output["xi_t_T"];
  EM_output["eta"] = MSVARloglik_output["eta"];
  EM_output["thl"] = thl;
  EM_output["deltath"] = deltath;
  EM_output["logLike"] = loglik_new;
  EM_output["residuals"] = Maxim_output["residuals"];
  EM_output["resid"] = Maxim_output["resid"];
  return(EM_output);
}

// ==============================================================================
//' @title Estimation of Hidden Markov model by EM Algorithm 
//' 
//' @description Estimate Hidden Markov model by EM algorithm. This function is used by \code{\link{HMmdl}} which organizes the output and takes raw data as input.
//' 
//' @param theta_0 vector with initital values for parameters
//' @param mdl List with model attributes
//' @param k integer determining the number of regimes
//' @param optim_options List with optimiztion options
//' 
//' @return List with model results
//' 
//' @export
// [[Rcpp::export]]
List HMmdl_em(arma::vec theta_0, List mdl, int k, List optim_options){
  // ---------- Get optimization options
  int maxit = optim_options["maxit"];
  double thtol = optim_options["thtol"];
  // ---------- Begin EM Algo
  // ----- Initial values
  List EMest_output;
  EMest_output["theta"] = theta_0;  // initial values for parameters
  EMest_output["logLike"] = 0;     // initial log-likelihood
  EMest_output = EMiter_HMmdl(mdl, EMest_output, k);
  int it = 1;
  double deltath = EMest_output["deltath"];
  // iterate to unitl convergence criterion is met or max number of iterations (maxit) 
  while ((it<maxit) & (deltath>thtol)){
    EMest_output = EMiter_HMmdl(mdl, EMest_output, k);
    deltath = EMest_output["deltath"];
    it = it + 1;
  }
  EMest_output["iterations"] = it; 
  return(EMest_output);
}



// ==============================================================================
//' @title Estimation of Markov-switching autoregressive model by EM Algorithm 
//' 
//' @description Estimate Markov-switching autoregressive model by EM algorithm. This function is used by \code{\link{MSARmdl}} which organizes the output and takes raw data as input.
//' 
//' @param theta_0 vector with initital values for parameters
//' @param mdl List with model attributes
//' @param k integer determining the number of regimes
//' @param optim_options List with optimiztion options
//' 
//' @return List with model results
//' 
//' @export
// [[Rcpp::export]]
List MSARmdl_em(arma::vec theta_0, List mdl, int k, List optim_options){
  // ---------- Get optimization options
  int maxit = optim_options["maxit"];
  double thtol = optim_options["thtol"];
  // ---------- Begin EM Algo
  // ----- Initial values
  List EMest_output;
  EMest_output["theta"] = theta_0;  // initial values for parameters
  EMest_output["logLike"] = 0;     // initial log-likelihood
  EMest_output = EMiter_MSARmdl(mdl, EMest_output, k);
  int it = 1;
  double deltath = EMest_output["deltath"];
  // iterate to unitl convergence criterion is met or max number of iterations (maxit) 
  while ((it<maxit) & (deltath>thtol)){
    EMest_output = EMiter_MSARmdl(mdl, EMest_output, k);
    deltath = EMest_output["deltath"];
    it = it + 1;
  }
  EMest_output["iterations"] = it; 
  return(EMest_output);
}


// ==============================================================================
//' @title Estimation of Markov-switching vector autoregressive model by EM Algorithm 
//' 
//' @description Estimate Markov-switching vector autoregressive model by EM algorithm. This function is used by \code{\link{MSVARmdl}} which organizes the output and takes raw data as input.
//' 
//' @param theta_0 vector with initital values for parameters
//' @param mdl List with model attributes
//' @param k integer determining the number of regimes
//' @param optim_options List with optimiztion options
//' 
//' @return List with model results
//' 
//' @export
// [[Rcpp::export]]
List MSVARmdl_em(arma::vec theta_0, List mdl, int k, List optim_options){
  // ---------- Get optimization options
  int maxit = optim_options["maxit"];
  double thtol = optim_options["thtol"];
  // ---------- Begin EM Algo
  // ----- Initial values
  List EMest_output;
  EMest_output["theta"] = theta_0;  // initial values for parameters
  EMest_output["logLike"] = 0;     // initial log-likelihood
  EMest_output = EMiter_MSVARmdl(mdl, EMest_output, k);
  int it = 1;
  double deltath = EMest_output["deltath"];
  // iterate to unitl convergence criterion is met or max number of iterations (maxit) 
  while ((it<maxit) & (deltath>thtol)){
    EMest_output = EMiter_MSVARmdl(mdl, EMest_output, k);
    deltath = EMest_output["deltath"];
    it = it + 1;
  }
  EMest_output["iterations"] = it; 
  return(EMest_output);
}

