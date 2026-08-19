#include <RcppArmadillo.h>
#include <RcppArmadilloExtensions/sample.h>
//[[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;


// ==============================================================================
//' @title Covariance to correlation matrix
//' 
//' @description This function takes an (\code{n x n}) covariance matrix and returns the associated (\code{n x n}) correlation matrix.
//' 
//' @param cov_mat A (\code{n x n}) covariance matrix.
//' 
//' @return A (\code{n x n}) correlation matrix.
//' 
//' @keywords internal
//' 
//' @export
// [[Rcpp::export]]
arma::mat cov2corr(arma::mat cov_mat){
  arma::vec diag_inv = 1.0 / sqrt(cov_mat.diag());
  arma::mat corr_mat = diagmat(diag_inv) * cov_mat * diagmat(diag_inv);
  return(corr_mat);
}

// ==============================================================================
//' @title Covariance vech function
//' 
//' @description This function returns the half-vectorization of an input matrix as a column vector.
//' 
//' @param mat A (\code{n x n}) covariance matrix.
//' 
//' @return A \code{(n+1)*n/2} column vector.
//' 
//' @keywords internal
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
//' @param sig A (n+1)*n/2 vector.
//' @param n Integer determining shape of the original matrix.
//' 
//' @return A (\code{n x n}) covariance matrix.
//' 
//' @keywords internal
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
//' @description This function creates a random transition matrix.
//' 
//' @param k Number of regimes. Must be greater than or equal to \code{2}. 
//' 
//' @return Transition matrix with randomly generated entries.
//' 
//' @keywords internal
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
//' @description Takes a transition matrix and returns the limiting probabilities.
//' 
//' @param P Matrix with transition probabilities.
//' 
//' @return A (\code{k x 1}) vector of limiting probabilities.
//' 
//' @keywords internal
//' 
//' @export
// [[Rcpp::export]]
arma::vec limP(arma::mat P){
  int nr = P.n_rows;
  int nc = P.n_cols;
  if (nc==nr){
    // Guard: if P is not finite, return uniform ergodic distribution
    if (!P.is_finite()) {
      arma::vec pinf_unif(nr, arma::fill::ones);
      return (pinf_unif / nr);
    }
    arma::mat onevec(1, nr, arma::fill::ones);
    arma::mat ep(1, nr+1, arma::fill::zeros);
    ep(0,nr) = 1.0;
    arma::mat Atmp = join_cols(arma::eye(nr,nr)-P,onevec);
    arma::vec pinf = solve(trans(Atmp)*Atmp,trans(Atmp), arma::solve_opts::allow_ugly)*trans(ep);
    // If solve returned non-finite (degenerate P), fall back to uniform
    if (!pinf.is_finite()) {
      pinf.fill(1.0 / nr);
    }
    // A stationary distribution is nonnegative by definition; the least-squares
    // solve can return tiny negative entries at boundary/rank-deficient P, which
    // corrupt the filter/smoother downstream (support monotonicity fails).
    // Conditional: a no-op (bit-identical) whenever the solution is nonnegative.
    if (pinf.min() < 0.0) {
      pinf = arma::clamp(pinf, 0.0, arma::datum::inf);
      double s = arma::accu(pinf);
      if (s > 0.0 && std::isfinite(s)) pinf /= s; else pinf.fill(1.0/nr);
    }
    return (pinf);
  }else{
    stop("Input must be a square matrix");
  }
}

// ==============================================================================
// O(M) forward prediction: computes P_AR * xi without forming the M×M matrix.
//
// Extended state n = i_0 + i_1*k + ... + i_p*k^p where i_j in {0,...,k-1}.
// P_AR(n', n) = P(j_0, j_1) * I[j_1=i_0] * I[j_2=i_1] * ... * I[j_p=i_{p-1}].
//
// Identity: [P_AR * xi](n') = P(j_0, j_1) * marginal(base)
//   where base = j_1 + j_2*k + ... + j_p*k^(p-1)  (k^p possible values)
//         n' = j_0 + k*base
//         j_1 = base % k
//         marginal(base) = sum_{i_p=0}^{k-1} xi(base + i_p * k^p)
//
// Cost: O(k^p * k) + O(k^(p+1)) = O(M).
arma::vec predict_fast(const arma::vec& xi, const arma::mat& P, int k, int ar) {
  int M = xi.n_elem;       // k^(ar+1)
  int kp = M / k;          // k^ar = number of "base" values
  // Step 1: marginalize over i_p (the highest-order component)
  arma::vec marginal(kp, arma::fill::zeros);
  for (int base = 0; base < kp; base++) {
    for (int ip = 0; ip < k; ip++) {
      marginal(base) += xi(base + ip * kp);
    }
  }
  // Step 2: multiply by transition probability
  arma::vec result(M);
  for (int base = 0; base < kp; base++) {
    int j1 = base % k;  // regime at lag 1 in the destination state
    for (int j0 = 0; j0 < k; j0++) {
      result(j0 + k * base) = P(j0, j1) * marginal(base);
    }
  }
  return result;
}

// ==============================================================================
// O(M) smoother step: computes P_AR' * w without forming the M×M matrix.
//
// Identity: [P_AR' * w](n) = sum_{j_0=0}^{k-1} P(j_0, i_0) * w(j_0 + i_0*k + i_1*k^2 + ... + i_{p-1}*k^p)
// The result is independent of i_p, so we compute per "tail base" and broadcast.
//
// For state n = i_0 + i_1*k + ... + i_p*k^p:
//   tail_base = i_0 + i_1*k + ... + i_{p-1}*k^(p-1)  (= n % k^p)
//   w_index = j_0 + k * (i_0 + i_1*k + ... + i_{p-1}*k^(p-1)) = j_0 + k*tail_base
//
// Cost: O(k^p * k) + O(M) = O(M).
arma::vec smooth_fast(const arma::vec& w, const arma::mat& P, int k, int ar) {
  int M = w.n_elem;        // k^(ar+1)
  int kp = M / k;          // k^ar
  // Step 1: compute sum for each tail_base
  arma::vec tail_sum(kp, arma::fill::zeros);
  for (int tb = 0; tb < kp; tb++) {
    int i0 = tb % k;  // regime at lag 0 (i_0)
    for (int j0 = 0; j0 < k; j0++) {
      tail_sum(tb) += P(j0, i0) * w(j0 + k * tb);
    }
  }
  // Step 2: broadcast — result is independent of i_p
  arma::vec result(M);
  for (int tb = 0; tb < kp; tb++) {
    for (int ip = 0; ip < k; ip++) {
      result(tb + ip * kp) = tail_sum(tb);
    }
  }
  return result;
}

// ==============================================================================
//' @title Lagged Time Series Data
//' 
//' @description This function takes a (\code{T x 1}) vector \code{Y} and returns the (\code{T-p x 1}) vector \code{y} and the (\code{T-p x p}) matrix of lagged observations.
//' 
//' @param Y Vector with time series observations.
//' @param p integer for the number of lags to use in estimation. Must be greater than or equal to \code{1}.
//' 
//' @return List with vector \code{y} (vector of lagged \code{Y}) and matrix \code{X} of lagged observations.
//' 
//' @keywords internal
//' 
//' @export
// [[Rcpp::export]]
List ts_lagged(arma::mat Y, int p){
  int Tsize = Y.n_rows;
  int N = Y.n_cols;
  arma::mat y = Y.submat(p,0,Tsize-1,N-1);
  int n = Tsize - p;
  arma::mat X(n, p*N, arma::fill::zeros);
  for (int xp = 0; xp<p; xp++){
    X.submat(0,N*xp,n-1,N*xp+N-1) = Y.submat(p-(xp+1),0,Tsize-(xp+1)-1,N-1);
  }
  List lagged_output;
  lagged_output["y"] = y;
  lagged_output["X"] = X;
  return(lagged_output);
}

// ==============================================================================
// Autoregressive moment grid and extended transition matrix (C++ ports) —
// definitions in argrid_ports.cpp; declarations here for the call sites.
List argrid_MSARmdl_cpp(arma::vec mu, arma::vec sig, int k, int ar);
List argrid_MSVARmdl_cpp(arma::mat mu, List sigma, int k, int ar);
arma::mat arP_cpp(arma::mat P, int k, int ar);


// ==============================================================================
//' @title Parameter list for Markov-switching autoregressive model
//' 
//' @description This function takes the parameter vector of interest and converts it to a list with specific parameter vectors needed for univariate Markov-switching functions.
//' 
//' @param theta Vector of parameters.
//' @param p Number of autoregressive lags.
//' @param k Number of regimes.
//' @param msmu Boolean indicating if the mean switches with regime.
//' @param msvar Boolean indicating if the variance switches with regime. 
//' 
//' @return List with the mean, variance, transition matrix, limiting probabilities, and a vector of state indicators.
//' 
//' @keywords internal
//' 
//' @export
// [[Rcpp::export]]
List paramList_MSARmdl(arma::vec theta, int p, int k, bool msmu, bool msvar){
  // ----- Mean for each regime 
  arma::vec mu        = theta.subvec(0, msmu*(k-1));
  // ----- Phi vector
  arma::vec phi       =  theta.subvec(1+msmu*(k-1), 1+msmu*(k-1)+p-1);
  // ----- Variance for each regime 
  arma::vec sig       = theta.subvec(1+msmu*(k-1)+p,1+msmu*(k-1)+p+msvar*(k-1));
  // ----- Transition probabilities 
  arma::mat P         = reshape(theta.subvec(2+msmu*(k-1)+p+msvar*(k-1), 2+msmu*(k-1)+p+msvar*(k-1) + k*k - 1),k, k);
  // Regime limiting probabilities
  arma::vec pinf      = limP(P);
  // ----- Obtain AR consistent grid of Mu, Sigma and State indicators
  List musig_out      = argrid_MSARmdl_cpp(mu, sig, k, p);
  arma::mat muAR      = as<arma::mat>(musig_out["mu"]);
  arma::mat sigAR     = as<arma::mat>(musig_out["sig"]);
  arma::vec state_ind = as<arma::vec>(musig_out["state_ind"]);
  // ----- Obtain AR consistent P and pinf
  arma::mat P_AR      = arP_cpp(P, k, p);
  arma::mat pinf_AR   = limP(P_AR);
  // ----- Organize output
  List param_out;
  param_out["mu"]         = mu;
  param_out["phi"]        = phi;
  param_out["sig"]        = sig;
  param_out["P"]          = P;
  param_out["pinf"]       = pinf;
  param_out["muAR"]       = muAR;
  param_out["sigAR"]      = sigAR;
  param_out["state_ind"]  = state_ind;
  param_out["P_AR"]       = P_AR;
  param_out["pinf_AR"]    = pinf_AR;
  return(param_out);
}


// ==============================================================================
//' @title Parameter list for Markov-switching ARX model
//' 
//' @description This function takes the parameter vector of interest and converts it to a list with specific parameter vectors needed for univariate Markov-switching functions.
//' 
//' @param theta Vector of parameters.
//' @param p Number of autoregressive lags.
//' @param k Number of regimes.
//' @param qz Number of exogenous variables.
//' @param msmu Boolean indicating if the mean switches with regime.
//' @param msvar Boolean indicating if the variance switches with regime. 
//' 
//' @return List with the mean, variance, transition matrix, limiting probabilities, and a vector of state indicators.
//' 
//' @keywords internal
//' 
//' @export
// [[Rcpp::export]]
List paramList_MSARXmdl(arma::vec theta, int p, int k, int qz, bool msmu, bool msvar){
  // ----- Mean for each regime 
  arma::vec mu          = theta.subvec(0, msmu*(k-1));
  // ----- Phi vector
  arma::vec phi         = theta.subvec(1+msmu*(k-1), 1+msmu*(k-1)+p-1);
  // ----- betaZ vector
  arma::vec betaZ       = theta.subvec(1+msmu*(k-1)+p, 1+msmu*(k-1)+p+qz-1);
  // ----- Variance for each regime 
  arma::vec sig         = theta.subvec(1+msmu*(k-1)+p+qz,1+msmu*(k-1)+p+qz+msvar*(k-1));
  // ----- Transition probabilities 
  arma::mat P           = reshape(theta.subvec(2+msmu*(k-1)+p+qz+msvar*(k-1), 2+msmu*(k-1)+p+qz+msvar*(k-1) + k*k - 1),k, k);
  // Regime limiting probabilities
  arma::vec pinf        = limP(P);
  // ----- Obtain AR consistent grid of Mu, Sigma and State indicators
  List musig_out        = argrid_MSARmdl_cpp(mu, sig, k, p);
  arma::mat muAR        = as<arma::mat>(musig_out["mu"]);
  arma::mat sigAR       = as<arma::mat>(musig_out["sig"]);
  arma::vec state_ind   = as<arma::vec>(musig_out["state_ind"]);
  // ----- Obtain AR consistent P and pinf
  arma::mat P_AR        = arP_cpp(P, k, p);
  arma::mat pinf_AR     = limP(P_AR);
  // ----- Organize output
  List param_out;
  param_out["mu"]         = mu;
  param_out["phi"]        = phi;
  param_out["betaZ"]      = betaZ;
  param_out["sig"]        = sig;
  param_out["P"]          = P;
  param_out["pinf"]       = pinf;
  param_out["muAR"]       = muAR;
  param_out["sigAR"]      = sigAR;
  param_out["state_ind"]  = state_ind;
  param_out["P_AR"]       = P_AR;
  param_out["pinf_AR"]    = pinf_AR;
  return(param_out);
}
// ==============================================================================
//' @title Parameter list for Markov-switching vector autoregressive model
//' 
//' @description This function takes the parameter vector of interest and converts it to a list with specific parameter vectors needed for multivariate Markov-switching functions.
//' 
//' @param theta Vector of parameters.
//' @param q Number of time series.
//' @param p Number of autoregressive lags.
//' @param k Number of regimes.
//' @param msmu Boolean indicating if the mean switches with regime.
//' @param msvar Boolean indicating if the variance switches with regime. 
//' 
//' @return List with the mean, variance, transition matrix, limiting probabilities, and a vector of state indicators.
//' 
//' @keywords internal
//' 
//' @export
// [[Rcpp::export]]
List paramList_MSVARmdl(arma::vec theta, int q, int p, int k, bool msmu, bool msvar){
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
  // ----- Phi vector
  arma::vec phi_tmp =  theta.subvec(q+q*msmu*(k-1), q+q*msmu*(k-1) + q*q*p - 1);
  arma::mat phi = trans(reshape(phi_tmp, q*p, q));
  // ----- Variance for each regime 
  int sigN = (q*(q+1))/2;
  arma::vec sig = theta.subvec(q+q*msmu*(k-1) + (q*q*p) ,q+q*msmu*(k-1) + q*q*p + (sigN + sigN*msvar*(k-1)) - 1);
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
  int PN = theta.n_elem - (k*k);
  arma::mat P = reshape(theta.subvec(PN, PN + k*k - 1), k, k);
  // Regime limiting probabilities
  arma::vec pinf = limP(P);
  // ----- Obtain AR consistent grid of Mu, Sigma and State indicators
  List musig_out = argrid_MSVARmdl_cpp(mu_k, sigma, k, p);
  List muAR = musig_out["mu"];
  List sigAR = musig_out["sig"];
  arma::vec state_ind = as<arma::vec>(musig_out["state_ind"]);
  // ----- Obtain AR consistent P and pinf
  //int M = pow(k, ar+1);
  arma::mat P_AR = arP_cpp(P, k, p);
  arma::mat pinf_AR = limP(P_AR);
  // ----- Organize output
  List param_out;
  param_out["mu"]         = mu_k;
  param_out["phi"]        = phi;
  param_out["sigma"]      = sigma;
  param_out["P"]          = P;
  param_out["pinf"]       = pinf;
  param_out["muAR"]       = muAR;
  param_out["sigAR"]      = sigAR;
  param_out["state_ind"]  = state_ind;
  param_out["P_AR"]       = P_AR;
  param_out["pinf_AR"]    = pinf_AR;
  return(param_out);
}
// ==============================================================================
//' @title Parameter list for Markov-switching VARX model
//' 
//' @description This function takes the parameter vector of interest and converts it to a list with specific parameter vectors needed for multivariate Markov-switching functions.
//' 
//' @param theta Vector of parameters.
//' @param q Number of time series.
//' @param p Number of autoregressive lags.
//' @param k Number of regimes.
//' @param qz Number of exogenous variables.
//' @param msmu Boolean indicating if the mean switches with regime.
//' @param msvar Boolean indicating if the variance switches with regime. 
//' 
//' @return List with the mean, variance, transition matrix, limiting probabilities, and a vector of state indicators.
//' 
//' @keywords internal
//' 
//' @export
// [[Rcpp::export]]
List paramList_MSVARXmdl(arma::vec theta, int q, int p, int k, int qz, bool msmu, bool msvar){
  // ----- Mean for each regime 
  arma::mat mu_k(k, q, arma::fill::zeros);
  arma::vec mu = theta.subvec(0, q + q*msmu*(k-1)-1);
  if (msmu==TRUE){
    for (int xk = 0; xk<k; xk++){
      mu_k.row(xk) = trans(mu.subvec(xk*q,xk*q+q-1));
    }
  }else{
    for (int xk = 0; xk<k; xk++){
      mu_k.row(xk) = trans(mu);
    }
  }
  // ----- Phi vector
  arma::vec phi_tmp =  theta.subvec(q + q*msmu*(k-1), q + q*msmu*(k-1) + (q*q*p) - 1);
  arma::mat phi = trans(reshape(phi_tmp, q*p, q));
  // ----- betaZ vector
  arma::vec betaZ_tmp  = theta.subvec(q + q*msmu*(k-1) + (q*q*p), q + q*msmu*(k-1) + (q*q*p) + qz*q - 1);
  arma::mat betaZ = reshape(betaZ_tmp,qz,q);
  // ----- Variance for each regime 
  int sigN = (q*(q+1))/2;
  arma::vec sig = theta.subvec(q + q*msmu*(k-1) + (q*q*p) + qz*q,q + q*msmu*(k-1) + (q*q*p) + qz*q + (sigN + sigN*msvar*(k-1)) - 1);
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
  int PN = theta.n_elem - (k*k);
  arma::mat P = reshape(theta.subvec(PN, PN + k*k - 1), k, k);
  // Regime limiting probabilities
  arma::vec pinf = limP(P);
  // ----- Obtain AR consistent grid of Mu, Sigma and State indicators
  List musig_out = argrid_MSVARmdl_cpp(mu_k, sigma, k, p);
  List muAR = musig_out["mu"];
  List sigAR = musig_out["sig"];
  arma::vec state_ind = as<arma::vec>(musig_out["state_ind"]);
  // ----- Obtain AR consistent P and pinf
  //int M = pow(k, ar+1);
  arma::mat P_AR = arP_cpp(P, k, p);
  arma::mat pinf_AR = limP(P_AR);
  // ----- Organize output
  List param_out;
  param_out["mu"]         = mu_k;
  param_out["phi"]        = phi;
  param_out["betaZ"]      = betaZ;
  param_out["sigma"]      = sigma;
  param_out["P"]          = P;
  param_out["pinf"]       = pinf;
  param_out["muAR"]       = muAR;
  param_out["sigAR"]      = sigAR;
  param_out["state_ind"]  = state_ind;
  param_out["P_AR"]       = P_AR;
  param_out["pinf_AR"]    = pinf_AR;
  return(param_out);
}
// ==============================================================================
//' @title Markov-switching autoregressive model residuals
//' 
//' @description This function computes residuals of a Markov-switching autoregressive model.
//' 
//' @param mdl List containing relevant parameters.
//' @param mu Vector with mean in each regime.
//' @param k number of regimes. Must be greater than or equal to \code{2}. 
//' 
//' @return A (\code{TxM}) matrix of residuals in each regime \code{M} where \code{M=k^(ar+1)}.
//' 
//' @keywords internal
//' 
//' @export
// [[Rcpp::export]]
arma::mat calcResid_MSARmdl(List mdl, arma::mat mu, int k){
  // ---------- Initialize parameters
  arma::vec y   = mdl["y"];
  arma::mat x   = mdl["x"];
  arma::vec phi = mdl["phi"];
  int ar        = mdl["p"];
  int Tsize     = y.n_elem;
  int M         = pow(k, ar+1); // number of regimes consistent with autoregressive structure (M=k if ar=0 and M=k^(ar+1) if ar>0)
  arma::mat repvec(1, M, arma::fill::ones);
  arma::mat repmu(Tsize, 1, arma::fill::ones);
  arma::mat ms_y = y*repvec;
  arma::mat resid(Tsize, M, arma::fill::zeros);
  // ---------- Compute Residuals 
  arma::mat z   = ms_y - repmu*trans(mu.col(0)); // [y(t) - mu_s(t))]
  arma::mat xz(Tsize, M, arma::fill::zeros);
  for (int xkp = 0; xkp<M; xkp++){
    arma::mat zx_tmp = x - repmu*mu.submat(xkp,1,xkp,ar); // [y(t-i) - mu_s(t-i))]
    xz.col(xkp) = zx_tmp*phi;
  }
  resid = z - xz;
  return(resid);
}
// ==============================================================================
//' @title Markov-switching autoregressive model residuals
//' 
//' @description This function computes residuals of a Markov-switching autoregressive model.
//' 
//' @param mdl List containing relevant parameters.
//' @param mu Vector with mean in each regime.
//' @param k number of regimes. Must be greater than or equal to \code{2}. 
//' 
//' @return A (\code{TxM}) matrix of residuals in each regime \code{M} where \code{M=k^(ar+1)}.
//' 
//' @keywords internal
//' 
//' @export
// [[Rcpp::export]]
arma::mat calcResid_MSARXmdl(List mdl, arma::mat mu, int k){
   // ---------- Initialize parameters
   arma::vec y        = mdl["y"];
   int ar             = mdl["p"];
   arma::mat x        = mdl["x"];
   List mdl_con       = mdl["control"];
   arma::mat Z        = mdl_con["Z"];
   arma::vec phi      = mdl["phi"];
   arma::vec betaZ    = mdl["betaZ"];
   int Tsize          = y.n_elem;
   int M              = pow(k, ar+1); // number of regimes consistent with autoregressive structure (M=k if ar=0 and M=k^(ar+1) if ar>0)
   Z = Z.rows(ar,Tsize+ar-1);
   arma::rowvec zbar  = arma::mean(Z,0);
   arma::mat repvec(1, M, arma::fill::ones);
   arma::mat repmu(Tsize, 1, arma::fill::ones);
   arma::mat ms_y     = y*repvec;
   arma::mat ms_zdm   = ((Z-repmu*zbar)*betaZ)*repvec;
   arma::mat resid(Tsize, M, arma::fill::zeros);
   // ---------- Compute Residuals 
   arma::mat ms_ydm = ms_y - repmu*trans(mu.col(0)); // [y(t) - mu_s(t))]
   arma::mat ms_xdm(Tsize, M, arma::fill::zeros);
   for (int xkp = 0; xkp<M; xkp++){
     arma::mat zx_tmp = x - repmu*mu.submat(xkp,1,xkp,ar); // [y(t-i) - mu_s(t-i))]
     ms_xdm.col(xkp) = zx_tmp*phi;
   }
   resid = ms_ydm - ms_xdm - ms_zdm;
   return(resid);
}
// ==============================================================================
//' @title Markov-switching vector autoregressive model residuals
//' 
//' @description This function computes residuals of a Markov-switching vector autoregressive model. 
//' 
//' @param mdl List containing relevant parameters.
//' @param mu Vector with mean in each regime.
//' @param k Number of regimes. Must be greater than or equal to \code{2}. 
//' 
//' @return List with \code{M} (\code{Txq}) matrices of residuals in each regime \code{M} where \code{M=k^(ar+1)}.
//' 
//' @keywords internal
//' 
//' @export
// [[Rcpp::export]]
List calcResid_MSVARmdl(List mdl, List mu, int k){
  // ---------- Initialize parameters
  arma::mat y   = mdl["y"];
  int p         = mdl["p"];
  arma::mat x   = mdl["x"];
  arma::mat phi = mdl["phi"];
  int Tsize     = y.n_rows;
  int q         = y.n_cols;
  int M         = pow(k, p+1); // number of regimes consistent with autoregressive structure (M=k if ar=0 and M=k^(ar+1) if ar>0)
  arma::mat repmu(Tsize, 1, arma::fill::ones);
  List resid(M); 
  for (int xm = 0; xm<M; xm++){
    arma::mat mu_tmp  = mu[xm]; 
    arma::mat y_tmp   = y - repmu*trans(mu_tmp.col(0));
    arma::mat xz_tmp(Tsize, q*p, arma::fill::zeros); 
    for (int xp = 0; xp<p; xp++){
      xz_tmp.submat(0,q*xp,Tsize-1,q*xp+q-1) = x.submat(0,q*xp,Tsize-1,q*xp+q-1) - repmu*trans(mu_tmp.col(xp+1));
    }
    resid[xm]   = y_tmp - xz_tmp*trans(phi);
  }
  return(resid);
}
// ==============================================================================
//' @title Markov-switching VARX model residuals
//' 
//' @description This function computes residuals of a Markov-switching VARX model. 
//' 
//' @param mdl List containing relevant parameters.
//' @param mu Vector with mean in each regime.
//' @param k Number of regimes. Must be greater than or equal to \code{2}. 
//' 
//' @return List with \code{M} (\code{Txq}) matrices of residuals in each regime \code{M} where \code{M=k^(ar+1)}.
//' 
//' @keywords internal
//' 
//' @export
// [[Rcpp::export]]
List calcResid_MSVARXmdl(List mdl, List mu, int k){
  // ---------- Initialize parameters
  arma::mat y       = mdl["y"];
  int p             = mdl["p"];
  arma::mat x       = mdl["x"];
  List mdl_con      = mdl["control"];
  arma::mat Z       = mdl_con["Z"];
  arma::mat phi     = mdl["phi"];
  arma::mat betaZ   = mdl["betaZ"];
  int Tsize         = y.n_rows;
  int q             = y.n_cols;
  int M             = pow(k, p+1); // number of regimes consistent with autoregressive structure (M=k if ar=0 and M=k^(ar+1) if ar>0)
  Z = Z.rows(p,Tsize+p-1);
  arma::rowvec zbar  = arma::mean(Z,0);
  arma::mat repmu(Tsize, 1, arma::fill::ones);
  arma::mat ms_zdm   = ((Z-repmu*zbar)*betaZ);
  List resid(M); 
   for (int xm = 0; xm<M; xm++){
     arma::mat mu_tmp  = mu[xm]; 
     arma::mat y_tmp   = y - repmu*trans(mu_tmp.col(0));
     arma::mat xz_tmp(Tsize, q*p, arma::fill::zeros); 
     for (int xp = 0; xp<p; xp++){
       xz_tmp.submat(0,q*xp,Tsize-1,q*xp+q-1) = x.submat(0,q*xp,Tsize-1,q*xp+q-1) - repmu*trans(mu_tmp.col(xp+1));
     }
     resid[xm]   = y_tmp - xz_tmp*trans(phi) - ms_zdm;
   }
   return(resid);
}
// ==============================================================================
//' @title Initial values for Hidden Markov model
//' 
//' @description This function generates a random parameter vector to be used as initial values for a Hidden Markov model.
//' 
//' @param mdl List with parameter values of simple (one-regime) model. This includes:
//' \itemize{
//'   \item mu: Vector of means.
//'   \item sigma: covariance matrix.
//'   \item msmu: Boolean indicator. If \code{TRUE}, mean is function of markov process. If \code{FALSE}, mean is constant across regimes.
//'   \item msvar: Boolean indicator. If \code{TRUE}, standard deviation is function of markov process. If \code{FALSE}, standard deviation is constant across regimes.
//' }
//' @param k Number of regimes.
//' 
//' @return Vector of initial parameter values.
//' 
//' @keywords internal
//' 
//' @export
// [[Rcpp::export]]
arma::vec initVals_HMmdl(List mdl, int k){
  arma::vec mu = mdl["mu"];
  arma::mat sigma = mdl["sigma"];
  bool msmu = mdl["msmu"];
  bool msvar = mdl["msvar"];
  bool exog = mdl["exog"];
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
  //arma::vec mu_out = vectorise(trans(mu_0));
  // initial values for stdev around linear model stdev when switch
  if (msvar==TRUE){
    arma::vec sigma_vec_tmp = covar_vech(sigma);
    arma::mat sig_mat_tmp = covar_unvech((sigma_vec_tmp*0.1) + ((2.0*sigma_vec_tmp)-(sigma_vec_tmp*0.1))%arma::randu<arma::vec>(sigN), q);
    sig_mat_tmp = sig_mat_tmp*trans(sig_mat_tmp);
    sig_mat_tmp = sig_mat_tmp + q*arma::speye(q,q);
    sigma_0.row(0) = trans(covar_vech(sig_mat_tmp));
    for (int xk = 1; xk<k; xk++){
      arma::mat sig_mat_tmp = covar_unvech((sigma_vec_tmp*0.1) + ((2.0*sigma_vec_tmp)-(sigma_vec_tmp*0.1))%arma::randu<arma::vec>(sigN), q);
      sig_mat_tmp = sig_mat_tmp*trans(sig_mat_tmp);
      sig_mat_tmp = sig_mat_tmp + q*arma::speye(q,q);
      sigma_0.row(xk) = trans(covar_vech(sig_mat_tmp));
    }
  }
  arma::vec sigma_out = vectorise(trans(sigma_0));
  // use estimates of exog regrssors as initial values (if exog is true)
  arma::vec theta_0 = vectorise(trans(mu_0));
  if (exog==TRUE){
    arma::vec betaZ = vectorise(as<arma::mat>(mdl["betaZ"]));
    theta_0 = join_vert(theta_0,betaZ);
  }
  // create vector for initial values
  theta_0 = join_vert(theta_0, sigma_out);
  arma::mat P_0 = randP(k);
  theta_0 = join_vert(theta_0, vectorise(P_0));
  return(theta_0);
}


// ==============================================================================
//' @title Initial values for Markov-switching autoregressive model
//' 
//' @description This function generates a random parameter vector to be used as initial values for a Markov-switching autoregressive model.
//' 
//' @param mdl List with parameter values of simple (one-regime) autoregressive model. This includes:
//' \itemize{
//'   \item phi: Vector autoregressive coefficients.
//'   \item mu: Mean of process.
//'   \item stdev: Standard deviation.
//'   \item msmu: Boolean indicator. If \code{TRUE}, mean is function of markov process. If \code{FALSE}, mean is constant across regimes.
//'   \item msvar: Boolean indicator. If \code{TRUE}, standard deviation is function of markov process. If \code{FALSE}, standard deviation is constant across regimes.
//' }
//' @param k Number of regimes.
//' 
//' @return Vector of initial parameter values.
//' 
//' @keywords internal
//' 
//' @export
// [[Rcpp::export]]
arma::vec initVals_MSARmdl(List mdl, int k){
  arma::vec phi = mdl["phi"];
  double mu     = mdl["mu"];
  double stdev  = mdl["stdev"];
  bool msmu     = mdl["msmu"];
  bool msvar    = mdl["msvar"];
  // pre-allocate mu and sigma vectors
  arma::vec mu_0(1+msmu*(k-1), arma::fill::zeros);
  arma::vec sig_0(1+msvar*(k-1), arma::fill::zeros);
  // Set initial values using linear model if no switch
  mu_0(0) = mu;
  sig_0(0) = pow(stdev,2.0);
  // initial values for mu around linear model mu when switch
  if (msmu==TRUE){
    mu_0 = mu + (3*stdev)*arma::randn<arma::vec>(k);
  }
  // initial values for stdev around linear model stdev when switch
  if (msvar==TRUE){
    sig_0 = (stdev*0.1) + ((2*stdev)-(stdev*0.1))*arma::randu<arma::vec>(k);
  }
  // create vector for initial values
  arma::vec theta_0 = join_vert(mu_0, phi);
  theta_0 = join_vert(theta_0, sig_0);
  arma::mat P_0 = randP(k);
  theta_0 = join_vert(theta_0, vectorise(P_0));
  return(theta_0);
}


// ==============================================================================
//' @title Initial values for Markov-switching ARX model
//' 
//' @description This function generates a random parameter vector to be used as initial values for a Markov-switching ARX model.
//' 
//' @param mdl List with parameter values of simple (one-regime) autoregressive model. This includes:
//' \itemize{
//'   \item phi: Vector autoregressive coefficients.
//'   \item mu: Mean of process.
//'   \item betaZ: vector of coefficients for exogenous regressors
//'   \item stdev: Standard deviation.
//'   \item msmu: Boolean indicator. If \code{TRUE}, mean is function of markov process. If \code{FALSE}, mean is constant across regimes.
//'   \item msvar: Boolean indicator. If \code{TRUE}, standard deviation is function of markov process. If \code{FALSE}, standard deviation is constant across regimes.
//' }
//' @param k Number of regimes.
//' 
//' @return Vector of initial parameter values.
//' 
//' @keywords internal
//' 
//' @export
// [[Rcpp::export]]
arma::vec initVals_MSARXmdl(List mdl, int k){
  arma::vec phi    = mdl["phi"];
  double mu        = mdl["mu"];
  double stdev     = mdl["stdev"];
  bool msmu        = mdl["msmu"];
  bool msvar       = mdl["msvar"];
  arma::vec betaz  = mdl["betaZ"];
  // pre-allocate mu and sigma vectors
  arma::vec mu_0(1+msmu*(k-1), arma::fill::zeros);
  arma::vec sig_0(1+msvar*(k-1), arma::fill::zeros);
  // Set initial values using linear model if no switch
  mu_0(0) = mu;
  sig_0(0) = pow(stdev,2.0);
  // initial values for mu around linear model mu when switch
  if (msmu==TRUE){
    mu_0 = mu + (3*stdev)*arma::randn<arma::vec>(k);
  }
  // initial values for stdev around linear model stdev when switch
  if (msvar==TRUE){
    sig_0 = (stdev*0.1) + ((2*stdev)-(stdev*0.1))*arma::randu<arma::vec>(k);
  }
  // create vector for initial values
  arma::vec theta_0 = join_vert(mu_0, phi);
  theta_0 = join_vert(theta_0, betaz);
  theta_0 = join_vert(theta_0, sig_0);
  arma::mat P_0 = randP(k);
  theta_0 = join_vert(theta_0, vectorise(P_0));
  return(theta_0);
}

// ==============================================================================
//' @title Initial values for Markov-switching vector autoregressive model
//' 
//' @description This function generates a random parameter vector to be used as initial values for a Markov-switching vector autoregressive model.
//' 
//' @param mdl List with parameter values of simple (one-regime) vector autoregressive model. This includes:
//'   \itemize{
//'    \item phi: Matrix autoregressive coefficients.
//'    \item mu: Vector of means.
//'    \item sigma: Covariance matrix.
//'    \item msmu: Boolean indicator. If \code{TRUE}, mean is function of markov process. If \code{FALSE}, mean is constant across regimes.
//'    \item msvar: Boolean indicator. If \code{TRUE}, standard deviation is function of markov process. If \code{FALSE}, standard deviation is constant across regimes.
//' }
//' @param k Number of regimes.
//' 
//' @return Vector of initial parameter values.
//' 
//' @keywords internal
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
    arma::mat sig_mat_tmp = covar_unvech((sigma_vec_tmp*0.1) + ((2.0*sigma_vec_tmp)-(sigma_vec_tmp*0.1))%arma::randu<arma::vec>(sigN), q);
    sig_mat_tmp = sig_mat_tmp*trans(sig_mat_tmp);
    sig_mat_tmp = sig_mat_tmp + q*arma::speye(q,q);
    sigma_0.row(0) = trans(covar_vech(sig_mat_tmp));
    for (int xk = 1; xk<k; xk++){
      arma::mat sig_mat_tmp = covar_unvech((sigma_vec_tmp*0.1) + ((2.0*sigma_vec_tmp)-(sigma_vec_tmp*0.1))%arma::randu<arma::vec>(sigN), q);
      sig_mat_tmp = sig_mat_tmp*trans(sig_mat_tmp);
      sig_mat_tmp = sig_mat_tmp + q*arma::speye(q,q);
      sigma_0.row(xk) = trans(covar_vech(sig_mat_tmp));
    }
  }
  arma::vec sigma_out = vectorise(trans(sigma_0));
  // create vector for initial values
  arma::vec theta_0 = join_vert(mu_out, vectorise(trans(phi)));
  theta_0 = join_vert(theta_0, sigma_out);
  arma::mat P_0 = randP(k);
  theta_0 = join_vert(theta_0, vectorise(P_0));
  return(theta_0);
}

// ==============================================================================
//' @title Initial values for Markov-switching VARX model
//' 
//' @description This function generates a random parameter vector to be used as initial values for a Markov-switching VARX model.
//' 
//' @param mdl List with parameter values of simple (one-regime) VARX model. This includes:
//'   \itemize{
//'    \item phi: Matrix autoregressive coefficients.
//'    \item mu: Vector of means.
//'    \item betaZ: vector of coefficients for exogenous regressors
//'    \item sigma: Covariance matrix.
//'    \item msmu: Boolean indicator. If \code{TRUE}, mean is function of markov process. If \code{FALSE}, mean is constant across regimes.
//'    \item msvar: Boolean indicator. If \code{TRUE}, standard deviation is function of markov process. If \code{FALSE}, standard deviation is constant across regimes.
//' }
//' @param k Number of regimes.
//' 
//' @return Vector of initial parameter values.
//' 
//' @keywords internal
//' 
//' @export
// [[Rcpp::export]]
arma::vec initVals_MSVARXmdl(List mdl, int k){
  arma::mat phi    = mdl["phi"];
  arma::vec mu     = mdl["mu"];
  arma::mat sigma  = mdl["sigma"];
  bool msmu        = mdl["msmu"];
  bool msvar       = mdl["msvar"];
  arma::mat betaz  = mdl["betaZ"];
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
    arma::mat sig_mat_tmp = covar_unvech((sigma_vec_tmp*0.1) + ((2.0*sigma_vec_tmp)-(sigma_vec_tmp*0.1))%arma::randu<arma::vec>(sigN), q);
    sig_mat_tmp = sig_mat_tmp*trans(sig_mat_tmp);
    sig_mat_tmp = sig_mat_tmp + q*arma::speye(q,q);
    sigma_0.row(0) = trans(covar_vech(sig_mat_tmp));
    for (int xk = 1; xk<k; xk++){
      arma::mat sig_mat_tmp = covar_unvech((sigma_vec_tmp*0.1) + ((2.0*sigma_vec_tmp)-(sigma_vec_tmp*0.1))%arma::randu<arma::vec>(sigN), q);
      sig_mat_tmp = sig_mat_tmp*trans(sig_mat_tmp);
      sig_mat_tmp = sig_mat_tmp + q*arma::speye(q,q);
      sigma_0.row(xk) = trans(covar_vech(sig_mat_tmp));
    }
  }
  arma::vec sigma_out = vectorise(trans(sigma_0));
  // create vector for initial values
  arma::vec theta_0 = join_vert(mu_out, vectorise(trans(phi)));
  theta_0 = join_vert(theta_0, vectorise(betaz));
  theta_0 = join_vert(theta_0, sigma_out);
  arma::mat P_0 = randP(k);
  theta_0 = join_vert(theta_0, vectorise(P_0));
  return(theta_0);
}


// ==============================================================================
//' @title Monte Carlo P-value
//' 
//' @description This function computes the Monte Carlo P-value.
//' 
//' @param test_stat Test statistic under the alternative (e.g. \code{S_0}).
//' @param null_vec A (\code{N x 1}) vector with test statistic under the null hypothesis.
//' @param type String determining the type of test. Options are \code{"geq"} (right/upper-tail test; the default), \code{"leq"} (left/lower-tail test), \code{"two-tailed"} (or \code{"two-tail"}; two-sided test), and \code{"absolute"} (or \code{"abs"}; the upper-tail rule, intended for non-negative statistics). An unrecognized value raises an error.
//' 
//' @return MC p-value of test
//' 
//' @references Dufour, Jean-Marie 2006. "Monte Carlo tests with nuisance parameters: A general approach to finite-sample inference and nonstandard asymptotics". \emph{Journal of Econometrics}, 133(2), 443-477.
//' @references Dufour, Jean-Marie, and Richard Luger. 2017. "Identification-robust moment-based tests for Markov switching in autoregressive models". \emph{Econometric Reviews}, 36(6-9), 713-727.
//' 
//' @export
// [[Rcpp::export]]
double MCpval(double test_stat, arma::vec null_vec, Rcpp::String type = "geq"){
  null_vec = null_vec.elem(arma::find_finite(null_vec)); // drop non-finite draws (NaN/Inf) before ranking
  int N = null_vec.n_elem;
  if (N == 0) {
    Rcpp::warning("MCpval: no finite null draws available; returning NA.");
    return R_NaN;
  }
  arma::vec u   = arma::randu(1 + sum(null_vec == test_stat)); 
  double test_rank  = sum(test_stat > null_vec) + sum(u <= u(0));
  // negative p-value can occur two-tailed test
  double survival_pval = ((N + 1 - test_rank)/ N);
  // Compute the p-value. Accepted 'type' values (long and short spellings):
  //   "geq"                         upper-tail (reject for large statistic; the default)
  //   "absolute" / "abs"            same upper-tail rule (for non-negative statistics)
  //   "leq"                         lower-tail
  //   "two-tailed" / "two-tail"     two-sided
  double pval;
  if (type == "geq" || type == "absolute" || type == "abs") {
    pval = (N * survival_pval + 1)/(N + 1);
  }else if (type == "leq") {
    pval = (N * (1 - survival_pval) + 1)/(N + 1);
  }else if (type == "two-tailed" || type == "two-tail") {
    pval = 2 * std::min((N * (1 - survival_pval) + 1)/(N + 1),(N * survival_pval + 1)/(N + 1));
  } else{
    Rcpp::stop("MCpval: 'type' must be one of 'geq', 'leq', 'two-tailed'/'two-tail', or 'absolute'/'abs'.");
  }
  return(pval);
}

// ==============================================================================
//' @title Standard normal errors using box Muller 
//' 
//' @description This function generates uncorrelated standard normal processes using box Muller method.
//' 
//' @param n Integer determining the length of the process to be simulated
//' @param q  Integer determining the number of processes to be simulated
//' 
//' @return A (\code{T x q}) matrix of standard normal distributed errors
//' 
//' @keywords internal
//' 
//' @export
// [[Rcpp::export]]
arma::mat randSN(int n, int q){
  // Use arma::randn: respects R's set.seed(), avoids log(0) risk from U1=0,
  // and is ~3-5x faster than Box-Muller.
  return arma::randn(n, q);
}


// ==============================================================================
//' @title Simulate autoregressive process
//' 
//' @description This function simulates an autoregresive process.
//' 
//' @param mdl_h0 List containing the following DGP parameters
//' \itemize{
//'   \item n: Length of series.
//'   \item mu: Mean of process.
//'   \item sigma: variance of process.
//'   \item phi: Vector of autoregressive coefficients.
//'   \item eps: An optional (\code{T+burnin x q}) matrix with standard normal errors to be used. Errors will be generated if not provided.
//' }
//' @param burnin Number of simulated observations to remove from beginning. Default is \code{100}.
//' 
//' @return List with simulated autoregressive series and its DGP parameters.
//' 
//' @keywords internal
//' 
//' @export
// [[Rcpp::export]]
List simuAR_cpp(List mdl_h0, int burnin = 100){
  // ----- DGP parameter
  arma::vec phi = mdl_h0["phi"];
  int Tsize     = mdl_h0["n"];
  double mu     = mdl_h0["mu"];
  double sigma  = mdl_h0["sigma"];
  int p         = phi.n_elem; 
  arma::mat eps;
  if(mdl_h0.containsElementNamed("eps") ){
    eps = as<arma::mat>(mdl_h0["eps"]);
  } else {
    eps = randSN(Tsize+burnin, 1); 
  }
  double stdev      = sqrt(sigma);
  double intercept  = mu*(1-sum(phi));
  // ----- Start simulation
  // pre-define variables
  arma::vec Y(Tsize+burnin, arma::fill::zeros);

  arma::vec eps_corr = eps*stdev;
  // simulate process
  Y.subvec(0, p-1) = mu + eps_corr.subvec(0, p-1);
  for (int xt = p; xt<(Tsize+burnin); xt++){
    arma::vec ytmp = flipud(Y.subvec((xt-p),(xt-1)));
    Y(xt) = arma::as_scalar(intercept + trans(ytmp)*phi + eps_corr(xt));
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
//' @title Simulate autoregressive process with exogenous regressors
//' 
//' @description This function simulates an ARX process.
//' 
//' @param mdl_h0 List containing the following DGP parameters
//' \itemize{
//'   \item n: Length of series.
//'   \item mu: Mean of process.
//'   \item sigma: variance of process.
//'   \item phi: Vector of autoregressive coefficients.
//'   \item eps: An optional (\code{T+burnin x q}) matrix with standard normal errors to be used. Errors will be generated if not provided.
//'   \item Z: A (\code{T x qz}) matrix with exogenous regressors (Optional) and where qz is the number of exogenous variables.
//'   \item betaZ: A (\code{qz x 1}) matrix  true coefficients on exogenous regressors (Optional) and where qz is the number of exogenous variables.
//' }
//' @param burnin Number of simulated observations to remove from beginning. Default is \code{100}.
//' 
//' @return List with simulated autoregressive series and its DGP parameters.
//' 
//' @keywords internal
//' 
//' @export
// [[Rcpp::export]]
List simuARX_cpp(List mdl_h0, int burnin = 100){
  // ----- DGP parameter
  arma::vec phi   = mdl_h0["phi"];
  int Tsize       = mdl_h0["n"];
  double mu       = mdl_h0["mu"];
  double sigma    = mdl_h0["sigma"];
  int p           = phi.n_elem; 
  arma::mat X     = mdl_h0["Z"];
  arma::rowvec Xbar  = arma::mean(X,0);
  arma::mat repv(Tsize,1,arma::fill::ones); 
  arma::mat Xdm   = X - repv*Xbar;
  arma::mat beta0 = mdl_h0["betaZ"];
  arma::mat eps;
  if(mdl_h0.containsElementNamed("eps")){
    eps = as<arma::mat>(mdl_h0["eps"]);
  } else {
    eps = randSN(Tsize+burnin, 1); 
  }
  double stdev      = sqrt(sigma);
  double intercept  = mu*(1-sum(phi));
  // ----- Start simulation
  // pre-define variables
  arma::vec Y(Tsize+burnin, arma::fill::zeros);
  arma::vec eps_corr = eps*stdev;
  // simulate process
  Y.subvec(0, p-1) = mu + eps_corr.subvec(0, p-1);
  for (int xt = p; xt<(Tsize+burnin); xt++){
    arma::vec ytmp = flipud(Y.subvec((xt-p),(xt-1)));
    Y(xt) = arma::as_scalar(intercept + trans(ytmp)*phi + eps_corr(xt));
    if (xt>=burnin){
      Y(xt) = Y(xt) + arma::as_scalar(Xdm.row(xt-burnin)*beta0);
    }
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
//' @description This function simulates a Markov-switching autoregressive process.
//' 
//' @param mdl_h0 List containing the following DGP parameters
//' \itemize{
//'   \item n: Length of series.
//'   \item k: Number of regimes.
//'   \item mu: A (\code{k x 1}) vector with mean of process in each regime.
//'   \item sigma: A (\code{k x 1}) vector with variance of process in each regime.
//'   \item phi: Vector of autoregressive coefficients.
//'   \item P: A (\code{k x k}) transition matrix (columns must sum to one).
//'   \item eps: An optional (\code{T+burnin x q}) matrix with standard normal errors to be used. Errors will be generated if not provided.
//' }
//' @param burnin Number of simulated observations to remove from beginning. Default is \code{100}.
//' 
//' @return List with simulated Markov-switching autoregressive process and its DGP properties.
//' 
//' @keywords internal
//' 
//' @export
// [[Rcpp::export]]
List simuMSAR_cpp(List mdl_h0, int burnin = 100){
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
  arma::mat eps;
  if(mdl_h0.containsElementNamed("eps") ){
    eps = as<arma::mat>(mdl_h0["eps"]);
  } else {
    eps = randSN(Tsize+burnin, 1);
  }
  // Pre-drawn state transition uniforms (for fixed-error MMC, Dufour 2006)
  arma::vec state_rand_vec;
  bool use_predrawn_states = mdl_h0.containsElementNamed("state_rand");
  if (use_predrawn_states) {
    state_rand_vec = as<arma::vec>(mdl_h0["state_rand"]);
  }
  // ----- Perform checks on DGP
  arma::vec check_Pcolsum = trans(arma::sum(P,0));
  if ((!P.is_finite()) || (max(abs(check_Pcolsum-1))>1e-8)){
    stop("Transition matrix 'P' must be finite with columns summing to 1.");
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
  for (int xk = 0; xk<k; xk++){
    double stdev_k = stdev(xk);
    arma::vec eps_corr = eps*stdev_k;
    epsLs[xk] = eps_corr;
  }
  // initialize assuming series begins in state 1 (use burnin to reduce dependence on this assumption)
  int state = 0;
  mu_t.rows(0,p-1) = repar*mu.row(state);
  stdev_t.rows(0,p-1) = repar*stdev.row(state);
  arma::vec eps_corr_k = epsLs[state];
  Y.rows(0,p-1) = mu_t.rows(0,p-1) + eps_corr_k.rows(0,p-1);
  // simulate process
  for (int xt = p; xt<(Tsize+burnin); xt++){
    // Get new state
    // Inverse CDF state selection (faster than find())
    double u_state = use_predrawn_states ? state_rand_vec(xt) : arma::randu<double>();
    int new_state = 0;
    double cum_p = P(0, state);
    while (u_state > cum_p && new_state < k - 1) { new_state++; cum_p += P(new_state, state); }
    state = new_state;
    state_series(xt) = state;
    // generate new obs
    arma::vec ytmp = flipud(Y.subvec((xt-p),(xt-1)));
    arma::vec mu_lag = flipud(mu_t.subvec((xt-p),(xt-1))); 
    arma::vec eps_corr_k = epsLs[state];
    resid(xt) = eps_corr_k(xt);
    Y(xt) = arma::as_scalar(mu(state) + (trans(ytmp-mu_lag))*phi + resid(xt));
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
//' @title Simulate Markov-switching ARX process
//' 
//' @description This function simulates a Markov-switching ARX process.
//' 
//' @param mdl_h0 List containing the following DGP parameters
//' \itemize{
//'   \item n: Length of series.
//'   \item k: Number of regimes.
//'   \item mu: A (\code{k x 1}) vector with mean of process in each regime.
//'   \item sigma: A (\code{k x 1}) vector with variance of process in each regime.
//'   \item phi: Vector of autoregressive coefficients.
//'   \item P: A (\code{k x k}) transition matrix (columns must sum to one).
//'   \item eps: An optional (\code{T+burnin x q}) matrix with standard normal errors to be used. Errors will be generated if not provided.
//'   \item Z: A (\code{T x qz}) matrix with exogenous regressors (Optional) and where qz is the number of exogenous variables.
//'   \item betaZ: A (\code{qz x 1}) matrix  true coefficients on exogenous regressors (Optional) and where qz is the number of exogenous variables.
//' }
//' @param burnin Number of simulated observations to remove from beginning. Default is \code{100}.
//' 
//' @return List with simulated Markov-switching autoregressive process and its DGP properties.
//' 
//' @keywords internal
//' 
//' @export
// [[Rcpp::export]]
List simuMSARX_cpp(List mdl_h0, int burnin = 100){
  // ----- DGP parameter
  arma::vec phi    = mdl_h0["phi"];
  arma::vec mu     = mdl_h0["mu"];
  arma::vec sigma  = mdl_h0["sigma"];
  arma::mat P      = mdl_h0["P"];
  int Tsize        = mdl_h0["n"];
  int k            = mdl_h0["k"];
  int p            = phi.n_elem; 
  arma::mat X      = mdl_h0["Z"];
  arma::rowvec Xbar  = arma::mean(X,0);
  arma::mat repmu(Tsize, 1, arma::fill::ones);
  arma::mat Xdm    = X - (repmu*Xbar);
  arma::mat beta0  = mdl_h0["betaZ"];
  arma::vec stdev  = sqrt(sigma);
  arma::vec pinf   = limP(P);
  arma::mat eps;
  if(mdl_h0.containsElementNamed("eps") ){
    eps = as<arma::mat>(mdl_h0["eps"]);
  }else{
    eps = randSN(Tsize+burnin, 1);
  }
  // Pre-drawn state transition uniforms (for fixed-error MMC, Dufour 2006)
  arma::vec state_rand_vec;
  bool use_predrawn_states = mdl_h0.containsElementNamed("state_rand");
  if (use_predrawn_states) {
    state_rand_vec = as<arma::vec>(mdl_h0["state_rand"]);
  }
  // ----- Perform checks on DGP
  arma::vec check_Pcolsum = trans(arma::sum(P,0));
  if ((!P.is_finite()) || (max(abs(check_Pcolsum-1))>1e-8)){
    stop("Transition matrix 'P' must be finite with columns summing to 1.");
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
  for (int xk = 0; xk<k; xk++){
    double stdev_k     = stdev(xk);
    arma::vec eps_corr = eps*stdev_k;
    epsLs[xk]          = eps_corr;
  }
  // initialize assuming series begins in state 1 (use burnin to reduce dependence on this assumption)
  int state            = 0;
  mu_t.rows(0,p-1)     = repar*mu.row(state);
  stdev_t.rows(0,p-1)  = repar*stdev.row(state);
  arma::vec eps_corr_k = epsLs[state];
  Y.rows(0,p-1)        = mu_t.rows(0,p-1) + eps_corr_k.rows(0,p-1);
  // simulate process
  for (int xt = p; xt<(Tsize+burnin); xt++){
    // Get new state
    // Inverse CDF state selection (faster than find())
    double u_state = use_predrawn_states ? state_rand_vec(xt) : arma::randu<double>();
    int new_state = 0;
    double cum_p = P(0, state);
    while (u_state > cum_p && new_state < k - 1) { new_state++; cum_p += P(new_state, state); }
    state = new_state;
    state_series(xt)     = state;
    // generate new obs
    arma::vec ytmp       = flipud(Y.subvec((xt-p),(xt-1)));
    arma::vec mu_lag     = flipud(mu_t.subvec((xt-p),(xt-1))); 
    arma::vec eps_corr_k = epsLs[state];
    resid(xt)            = eps_corr_k(xt);
    Y(xt)    = arma::as_scalar(mu(state) + (trans(ytmp-mu_lag))*phi + resid(xt));
    if (xt>=burnin){
      Y(xt)  = Y(xt) + arma::as_scalar(Xdm.row(xt-burnin)*beta0);
    }
    mu_t(xt)     = mu(state);
    stdev_t(xt)  = stdev(state);
  }
  // ----- Output
  arma::vec Y_out        = Y.subvec(burnin,Tsize+burnin-1);
  arma::vec resid_out    = resid.subvec(burnin,Tsize+burnin-1);
  arma::vec state_series_out = state_series.subvec(burnin,Tsize+burnin-1);
  arma::vec mu_t_out     = mu_t.subvec(burnin,Tsize+burnin-1);
  arma::vec stdev_t_out  = stdev_t.subvec(burnin,Tsize+burnin-1);
  List simu_output       = clone(mdl_h0);
  simu_output["y"]       = Y_out;
  simu_output["p"]       = p;
  simu_output["resid"]   = resid_out;
  simu_output["St"]      = state_series_out;
  simu_output["mu_t"]    = mu_t_out;
  simu_output["stdev_t"] = stdev_t_out;
  simu_output["pinf"]    = pinf;
  return(simu_output);
}
// ==============================================================================
//' @title Simulate VAR process
//' 
//' @description This function simulates a vector autoregresive process.
//' 
//' @param mdl_h0 List containing the following DGP parameters
//' \itemize{
//'   \item n: Length of series.
//'   \item mu: A (\code{q x 1}) vector of means.
//'   \item sigma: A (\code{q x q}) covariance matrix.
//'   \item phi:  A (\code{q x qp}) matrix of autoregressive coefficients.
//'   \item p: Number of autoregressive lags.
//'   \item q: Number of series.
//'   \item eps: An optional (\code{T+burnin x q}) matrix with standard normal errors to be used. Errors will be generated if not provided.
//' }
//' @param burnin Number of simulated observations to remove from beginning. Default is \code{100}.
//' 
//' @return List with simulated vector autoregressive series and its DGP parameters.
//' 
//' @keywords internal
//' 
//' @export
// [[Rcpp::export]]
List simuVAR_cpp(List mdl_h0, int burnin = 100){
  // ----- DGP parameters
  arma::vec mu = mdl_h0["mu"];
  arma::mat cov_mat = mdl_h0["sigma"];
  arma::mat phimat = mdl_h0["phi"];
  int Tsize = mdl_h0["n"];
  int p = mdl_h0["p"];
  int q = mdl_h0["q"];
  arma::mat eps;
  if(mdl_h0.containsElementNamed("eps") ){
    eps = as<arma::mat>(mdl_h0["eps"]);
  } else {
    eps = randSN(Tsize+burnin, q); 
  }
  // companion matrix
  arma::mat diagmat = arma::eye(q*(p-1),q*(p-1));
  arma::mat diagzero(q*(p-1),q,arma::fill::zeros);
  arma::mat Mn = join_rows(diagmat,diagzero);
  arma::mat F = join_cols(phimat,Mn);
  // constant vec
  arma::vec repmu(p,arma::fill::ones);
  arma::vec mu_tmp = vectorise(trans(repmu*trans(mu)));
  arma::vec nu_tmp = (arma::eye(q*p,q*p) - F)*mu_tmp;
  arma::vec nu = nu_tmp.subvec(0,q-1);
  // ----- Simulate VAR process
  // add standard devs & correlations
  arma::mat C = chol(cov_mat, "lower");
  arma::mat eps_corr = trans(C*trans(eps));
  // simulate process
  arma::mat Y(Tsize+burnin, q, arma::fill::zeros);
  Y.rows(0,p-1) = repmu*trans(nu) + eps_corr.rows(0,p-1);
  for (int xt = p; xt<(Tsize+burnin); xt++){
    arma::mat Ytmp = flipud(Y.rows((xt-p),(xt-1)));
    Y.row(xt) = trans(nu) + trans(vectorise(trans(Ytmp)))*trans(phimat) + eps_corr.row(xt);
  }
  // ----- Output
  arma::mat Y_out = Y.submat(burnin,0,Tsize+burnin-1,q-1);
  arma::mat eps_corr_out = eps_corr.submat(burnin,0,Tsize+burnin-1,q-1);
  List simuVAR_out = clone(mdl_h0);
  simuVAR_out["y"] = Y_out;
  simuVAR_out["resid"] = eps_corr_out;
  simuVAR_out["F_comp"] = F;
  return(simuVAR_out);
}


// ==============================================================================
//' @title Simulate VARX process
//' 
//' @description This function simulates a VARX process.
//' 
//' @param mdl_h0 List containing the following DGP parameters
//' \itemize{
//'   \item n: Length of series.
//'   \item mu: A (\code{q x 1}) vector of means.
//'   \item sigma: A (\code{q x q}) covariance matrix.
//'   \item phi:  A (\code{q x qp}) matrix of autoregressive coefficients.
//'   \item p: Number of autoregressive lags.
//'   \item q: Number of series.
//'   \item eps: An optional (\code{T+burnin x q}) matrix with standard normal errors to be used. Errors will be generated if not provided.
//'   \item Z: A (\code{T x qz}) matrix with exogenous regressors (Optional) and where qz is the number of exogenous variables.
//'   \item betaZ: A (\code{qz x q}) matrix  true coefficients on exogenous regressors (Optional) and where qz is the number of exogenous variables.
//' }
//' @param burnin Number of simulated observations to remove from beginning. Default is \code{100}.
//' 
//' @return List with simulated vector autoregressive series and its DGP parameters.
//' 
//' @keywords internal
//' 
//' @export
// [[Rcpp::export]]
List simuVARX_cpp(List mdl_h0, int burnin = 100){
  // ----- DGP parameters
  arma::vec mu      = mdl_h0["mu"];
  arma::mat cov_mat = mdl_h0["sigma"];
  arma::mat phimat  = mdl_h0["phi"];
  int Tsize         = mdl_h0["n"];
  int p             = mdl_h0["p"];
  int q             = mdl_h0["q"];
  arma::mat X       = mdl_h0["Z"];
  arma::rowvec Xbar = arma::mean(X,0);
  arma::mat repmu(Tsize, 1, arma::fill::ones);
  arma::mat Xdm     = X - (repmu*Xbar);
  arma::mat beta0   = mdl_h0["betaZ"];
  arma::mat eps;
  if(mdl_h0.containsElementNamed("eps") ){
    eps = as<arma::mat>(mdl_h0["eps"]);
  } else {
    eps = randSN(Tsize+burnin, q); 
  }
  // companion matrix
  arma::mat diagmat = arma::eye(q*(p-1),q*(p-1));
  arma::mat diagzero(q*(p-1),q,arma::fill::zeros);
  arma::mat Mn = join_rows(diagmat,diagzero);
  arma::mat F = join_cols(phimat,Mn);
  // constant vec (need p-element ones vector, not the T-element repmu used for X demeaning)
  arma::vec repmu_p(p, arma::fill::ones);
  arma::vec mu_tmp = vectorise(trans(repmu_p*trans(mu)));
  arma::vec nu_tmp = (arma::eye(q*p,q*p) - F)*mu_tmp;
  arma::vec nu = nu_tmp.subvec(0,q-1);
  // ----- Simulate VAR process
  // add standard devs & correlations
  arma::mat C = chol(cov_mat, "lower");
  arma::mat eps_corr = trans(C*trans(eps));
  // simulate process
  arma::mat Y(Tsize+burnin, q, arma::fill::zeros);
  Y.rows(0,p-1) = repmu_p*trans(nu) + eps_corr.rows(0,p-1);
  for (int xt = p; xt<(Tsize+burnin); xt++){
    arma::mat Ytmp = flipud(Y.rows((xt-p),(xt-1)));
    Y.row(xt) = trans(nu) + trans(vectorise(trans(Ytmp)))*trans(phimat) + eps_corr.row(xt);
    if (xt>=burnin){
      Y.row(xt) = Y.row(xt) + Xdm.row(xt-burnin)*beta0;
    }
  }
  // ----- Output
  arma::mat Y_out = Y.submat(burnin,0,Tsize+burnin-1,q-1);
  arma::mat eps_corr_out = eps_corr.submat(burnin,0,Tsize+burnin-1,q-1);
  List simuVAR_out = clone(mdl_h0);
  simuVAR_out["y"] = Y_out;
  simuVAR_out["resid"] = eps_corr_out;
  simuVAR_out["F_comp"] = F;
  return(simuVAR_out);
}


// ==============================================================================
//' @title Simulate Markov-switching vector autoregressive process
//' 
//' @description This function simulates a Markov-switching vector autoregressive process.
//' 
//' @param mdl_h0 List containing the following DGP parameters
//' \itemize{
//'   \item n: Length of series.
//'   \item k: Number of regimes.
//'   \item mu: A (\code{k x q}) matrix of means.
//'   \item sigma: List with \code{k} (\code{q x q}) covariance matrices.
//'   \item phi: A (\code{q x qp}) matrix of autoregressive coefficients.
//'   \item p: Number of autoregressive lags.
//'   \item q: Number of series.
//'   \item P: A (\code{k x k}) transition matrix (columns must sum to one).
//'   \item eps: An optional (\code{T+burnin x q}) matrix with standard normal errors to be used. Errors will be generated if not provided.
//' }
//' @param burnin Number of simulated observations to remove from beginning. Default is \code{100}.
//' 
//' @return List with simulated vector autoregressive series and its DGP parameters.
//' 
//' @keywords internal
//' 
//' @export
// [[Rcpp::export]]
List simuMSVAR_cpp(List mdl_h0, int burnin = 100){
  // ----- DGP parameters
  arma::mat mu = mdl_h0["mu"];
  List cov_matLs = mdl_h0["sigma"];
  arma::mat phimat = mdl_h0["phi"];
  arma::mat P = mdl_h0["P"];
  arma::vec  pinf = limP(P);
  int Tsize = mdl_h0["n"];
  int p = mdl_h0["p"];
  int k = mdl_h0["k"];
  int q = mdl_h0["q"];
  arma::mat eps;
  if(mdl_h0.containsElementNamed("eps") ){
    eps = as<arma::mat>(mdl_h0["eps"]);
  } else {
    eps = randSN(Tsize+burnin, q);
  }
  // Pre-drawn state transition uniforms (for fixed-error MMC, Dufour 2006)
  arma::vec state_rand_vec;
  bool use_predrawn_states = mdl_h0.containsElementNamed("state_rand");
  if (use_predrawn_states) {
    state_rand_vec = as<arma::vec>(mdl_h0["state_rand"]);
  }
  // companion form matrix
  arma::mat diagmat = arma::eye(q*(p-1), q*(p-1));
  arma::mat diagzero(q*(p-1), q, arma::fill::zeros);
  arma::mat Mn = join_rows(diagmat,diagzero);
  arma::mat F = join_cols(phimat,Mn);
  // ----- Perform checks on DGP
  arma::vec check_Pcolsum = trans(arma::sum(P,0));
  if ((!P.is_finite()) || (max(abs(check_Pcolsum-1))>1e-8)){
    stop("Transition matrix 'P' must be finite with columns summing to 1.");
  }
  // ----- Start Simulation
  // pre-define variables
  arma::mat mu_t(Tsize+burnin, q, arma::fill::zeros);
  List sigma_t(Tsize+burnin);
  arma::mat Y(Tsize+burnin, q, arma::fill::zeros);
  arma::mat resid(Tsize+burnin, q, arma::fill::zeros);
  arma::vec state_series(Tsize+burnin,arma::fill::zeros);
  arma::vec repar(p,arma::fill::ones);
  // simulate errors using box-Muller method
  List epsLs(k);
  for (int xk = 0; xk<k; xk++){
    // add correlations
    arma::mat cov_mat_k = cov_matLs[xk];
    arma::mat C = chol(cov_mat_k, "lower");
    arma::mat eps_corr = trans(C*trans(eps));
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
  for (int xt = p; xt<(Tsize+burnin); xt++){
    // Get new state
    // Inverse CDF state selection (faster than find())
    double u_state = use_predrawn_states ? state_rand_vec(xt) : arma::randu<double>();
    int new_state = 0;
    double cum_p = P(0, state);
    while (u_state > cum_p && new_state < k - 1) { new_state++; cum_p += P(new_state, state); }
    state = new_state;
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
  arma::mat Y_out = Y.submat(burnin,0,Tsize+burnin-1,q-1);
  arma::mat resid_out = resid.submat(burnin,0,Tsize+burnin-1,q-1);
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
//' @title Simulate Markov-switching VARX process
//' 
//' @description This function simulates a Markov-switching VARX process.
//' 
//' @param mdl_h0 List containing the following DGP parameters
//' \itemize{
//'   \item n: Length of series.
//'   \item k: Number of regimes.
//'   \item mu: A (\code{k x q}) matrix of means.
//'   \item sigma: List with \code{k} (\code{q x q}) covariance matrices.
//'   \item phi: A (\code{q x qp}) matrix of autoregressive coefficients.
//'   \item p: Number of autoregressive lags.
//'   \item q: Number of series.
//'   \item P: A (\code{k x k}) transition matrix (columns must sum to one).
//'   \item eps: An optional (\code{T+burnin x q}) matrix with standard normal errors to be used. Errors will be generated if not provided.
//'   \item Z: A (\code{T x qz}) matrix with exogenous regressors (Optional) and where qz is the number of exogenous variables.
//'   \item betaZ: A (\code{qz x q}) matrix  true coefficients on exogenous regressors (Optional) and where qz is the number of exogenous variables.
//' }
//' @param burnin Number of simulated observations to remove from beginning. Default is \code{100}.
//' 
//' @return List with simulated vector autoregressive series and its DGP parameters.
//' 
//' @keywords internal
//' 
//' @export
// [[Rcpp::export]]
List simuMSVARX_cpp(List mdl_h0, int burnin = 100){
   // ----- DGP parameters
   arma::mat mu       = mdl_h0["mu"];
   List cov_matLs     = mdl_h0["sigma"];
   arma::mat phimat   = mdl_h0["phi"];
   arma::mat P        = mdl_h0["P"];
   arma::vec  pinf    = limP(P);
   int Tsize          = mdl_h0["n"];
   int p              = mdl_h0["p"];
   int k              = mdl_h0["k"];
   int q              = mdl_h0["q"];
   arma::mat X        = mdl_h0["Z"];
   arma::rowvec Xbar  = arma::mean(X,0);
   arma::mat repmu(Tsize, 1, arma::fill::ones);
   arma::mat Xdm      = X - (repmu*Xbar);
   arma::mat beta0    = mdl_h0["betaZ"];
   arma::mat eps;
   if(mdl_h0.containsElementNamed("eps") ){
     eps = as<arma::mat>(mdl_h0["eps"]);
   } else {
     eps = randSN(Tsize+burnin, q);
   }
   // Pre-drawn state transition uniforms (for fixed-error MMC, Dufour 2006)
   arma::vec state_rand_vec;
   bool use_predrawn_states = mdl_h0.containsElementNamed("state_rand");
   if (use_predrawn_states) {
     state_rand_vec = as<arma::vec>(mdl_h0["state_rand"]);
   }
   // companion form matrix
   arma::mat diagmat = arma::eye(q*(p-1), q*(p-1));
   arma::mat diagzero(q*(p-1), q, arma::fill::zeros);
   arma::mat Mn = join_rows(diagmat,diagzero);
   arma::mat F = join_cols(phimat,Mn);
   // ----- Perform checks on DGP
   arma::vec check_Pcolsum = trans(arma::sum(P,0));
   if ((!P.is_finite()) || (max(abs(check_Pcolsum-1))>1e-8)){
     stop("Transition matrix 'P' must be finite with columns summing to 1.");
   }
   // ----- Start Simulation
   // pre-define variables
   arma::mat mu_t(Tsize+burnin, q, arma::fill::zeros);
   List sigma_t(Tsize+burnin);
   arma::mat Y(Tsize+burnin, q, arma::fill::zeros);
   arma::mat resid(Tsize+burnin, q, arma::fill::zeros);
   arma::vec state_series(Tsize+burnin,arma::fill::zeros);
   arma::vec repar(p,arma::fill::ones);
   // simulate errors using box-Muller method
   List epsLs(k);
   for (int xk = 0; xk<k; xk++){
     // add correlations
     arma::mat cov_mat_k = cov_matLs[xk];
     arma::mat C = chol(cov_mat_k, "lower");
     arma::mat eps_corr = trans(C*trans(eps));
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
   for (int xt = p; xt<(Tsize+burnin); xt++){
     // Get new state
     // Inverse CDF state selection (faster than find())
     double u_state = use_predrawn_states ? state_rand_vec(xt) : arma::randu<double>();
     int new_state = 0;
     double cum_p = P(0, state);
     while (u_state > cum_p && new_state < k - 1) { new_state++; cum_p += P(new_state, state); }
     state = new_state;
     state_series(xt) = state;
     arma::mat eps_corr_k = epsLs[state];
     arma::mat Ytmp = flipud(Y.rows((xt-p),(xt-1)));
     arma::mat mu_lag = flipud(mu_t.rows((xt-p),(xt-1)));
     resid.row(xt) = eps_corr_k.row(xt);
     Y.row(xt) = mu.row(state) + trans(vectorise(trans(Ytmp-mu_lag)))*trans(phimat) + resid.row(xt);
     if (xt>=burnin){
       Y.row(xt) = Y.row(xt) + Xdm.row(xt-burnin)*beta0;
     }
     mu_t.row(xt) = mu.row(state);
     sigma_t[xt] = cov_matLs[state];
   }
   // ----- Output
   arma::mat Y_out = Y.submat(burnin,0,Tsize+burnin-1,q-1);
   arma::mat resid_out = resid.submat(burnin,0,Tsize+burnin-1,q-1);
   arma::vec state_series_out = state_series.subvec(burnin,Tsize+burnin-1);
   arma::mat mu_t_out = mu_t.rows(burnin,Tsize+burnin-1);
   List sigma_t_out(Tsize);
   for (int xt = burnin; xt<(Tsize+burnin); xt++){
     sigma_t_out[xt-burnin] = sigma_t[xt];
   }
   List simuVAR_out       = clone(mdl_h0);
   simuVAR_out["y"]       = Y_out;
   simuVAR_out["p"]       = p;
   simuVAR_out["resid"]   = resid_out;
   simuVAR_out["F_comp"]  = F;
   simuVAR_out["St"]      = state_series_out;
   simuVAR_out["mu_t"]    = mu_t_out;
   simuVAR_out["sigma_t"] = sigma_t_out;
   simuVAR_out["pinf"]    = pinf;
   return(simuVAR_out);
}

// ==============================================================================
//' @title Simulate normally distributed process
//' 
//' @description This function simulates a normally distributed process.
//' 
//' @param mdl_h0 List containing the following DGP parameters
//' \itemize{
//'   \item n: Length of series.
//'   \item mu: A (\code{q x 1}) vector of means.
//'   \item sigma: A (\code{q x q}) covariance matrix.
//'   \item q: Number of series.
//'   \item eps: An optional (\code{T+burnin x q}) matrix with standard normal errors to be used. Errors will be generated if not provided.
//'   \item Z: A (\code{T x qz}) matrix with exogenous regressors (Optional) and where qz is the number of exogenous variables.
//'   \item betaZ: A (\code{qz x q}) matrix  true coefficients on exogenous regressors (Optional) and where qz is the number of exogenous variables.
//' }
//' @param burnin Number of simulated observations to remove from beginning. Default is \code{100}.
//' @param exog bool determining if there are exogenous variables (\code{true}) or not (\code{false}). Default is \code{false}.
//' 
//' @return List with simulated series and its DGP parameters.
//' 
//' @keywords internal
//' 
//' @export
// [[Rcpp::export]]
List simuNorm_cpp(List mdl_h0, int burnin = 0, bool exog = false){
  // ----- DGP parameter
  arma::vec mu  = mdl_h0["mu"];
  int Tsize     = mdl_h0["n"];
  int q         = mdl_h0["q"];
  arma::mat eps;
  if(mdl_h0.containsElementNamed("eps") ){
    eps = as<arma::mat>(mdl_h0["eps"]);
  } else {
    eps = randSN(Tsize+burnin, q); 
  }
  // ----- Start simulation
  // pre-define variables
  arma::vec repT(Tsize+burnin, arma::fill::ones);
  // add standard dev & correlations
  arma::mat cov_mat   = mdl_h0["sigma"];
  arma::mat C         = chol(cov_mat, "lower");
  arma::mat eps_corr  = trans(C*trans(eps));
  // simulate process
  arma::mat Y  =  repT*trans(mu) + eps_corr;  
  arma::mat Y_out = Y.submat(burnin,0,Tsize+burnin-1,q-1);
  if (exog==TRUE){
    arma::mat ZZ      = mdl_h0["Z"];
    arma::rowvec zbar  = arma::mean(ZZ,0);
    arma::mat repmu(Tsize, 1, arma::fill::ones);
    arma::mat Zdm = ZZ - (repmu*zbar);
    arma::mat beta0   = mdl_h0["betaZ"];
    Y_out = Y_out + Zdm*beta0; 
  }
  arma::mat resid_out = eps_corr.submat(burnin,0,Tsize+burnin-1,q-1);
  // ----- Output
  List simuNorm_out = clone(mdl_h0);
  simuNorm_out["y"] = Y_out;
  simuNorm_out["resid"] = resid_out;
  return(simuNorm_out);
}


// ==============================================================================
//' @title Simulate Hidden Markov model with normally distributed errors
//' 
//' @description This function simulates a Hidden Markov Model process.
//' 
//' @param mdl_h0 List containing the following DGP parameters
//' \itemize{
//'   \item n: Length of series.
//'   \item k: Number of regimes.
//'   \item mu: A (\code{k x q}) vector of means.
//'   \item sigma: A (\code{q x q}) covariance matrix.
//'   \item q: Number of series.
//'   \item P: A (\code{k x k}) transition matrix (columns must sum to one).
//'   \item eps: An optional (\code{T+burnin x q}) matrix with standard normal errors to be used. Errors will be generated if not provided.
//'   \item Z: A (\code{T x qz}) matrix with exogenous regressors (Optional) and where qz is the number of exogenous variables.
//'   \item betaZ: A (\code{qz x q}) matrix true coefficients on exogenous regressors (Optional) and where qz is the number of exogenous variables.
//' }
//' @param burnin Number of simulated observations to remove from beginning. Default is \code{100}.
//' @param exog bool determining if there are exogenous variables (\code{true}) or not (\code{false}). Default is \code{false}.
//' 
//' @return List with simulated series and its DGP parameters.
//' 
//' @keywords internal
//' 
//' @export
// [[Rcpp::export]]
List simuHMM_cpp(List mdl_h0, int burnin = 100, bool exog = false){
  // ----- DGP parameter
  arma::mat mu    = mdl_h0["mu"];
  List cov_matLs  = mdl_h0["sigma"];
  arma::mat P     = mdl_h0["P"];
  arma::vec  pinf = limP(P);
  int Tsize       = mdl_h0["n"];
  int k           = mdl_h0["k"];
  int q           = mdl_h0["q"];
  arma::mat Zdm(Tsize,1,arma::fill::zeros);
  arma::mat beta0(1,q,arma::fill::zeros);
  if (exog==TRUE){
    arma::mat ZZ  = as<arma::mat>(mdl_h0["Z"]);
    beta0         = as<arma::mat>(mdl_h0["betaZ"]);
    arma::rowvec zbar  = arma::mean(ZZ,0);
    arma::mat repmu(Tsize, 1, arma::fill::ones);
    Zdm = ZZ - (repmu*zbar);
  }
  arma::mat eps;
  if(mdl_h0.containsElementNamed("eps") ){
    eps = as<arma::mat>(mdl_h0["eps"]);
  } else {
    eps = randSN(Tsize+burnin, q);
  }
  // Pre-drawn state transition uniforms (for fixed-error MMC, Dufour 2006)
  arma::vec state_rand_vec;
  bool use_predrawn_states = mdl_h0.containsElementNamed("state_rand");
  if (use_predrawn_states) {
    state_rand_vec = as<arma::vec>(mdl_h0["state_rand"]);
  }
  // ----- Perform checks on DGP
  arma::vec check_Pcolsum = trans(arma::sum(P,0));
  if ((!P.is_finite()) || (max(abs(check_Pcolsum-1))>1e-8)){
    stop("Transition matrix 'P' must be finite with columns summing to 1.");
  }
  // ----- Start Simulation
  // pre-define variables
  arma::mat mu_t(Tsize+burnin, q, arma::fill::zeros);
  List sigma_t(Tsize+burnin);
  arma::mat Y(Tsize+burnin, q, arma::fill::zeros);
  arma::mat resid(Tsize+burnin, q, arma::fill::zeros);
  arma::vec state_series(Tsize+burnin, arma::fill::zeros);
  // simulate errors using box-Muller method
  List epsLs(k);
  for (int xk = 0; xk<k; xk++){
    // add correlations
    arma::mat cov_mat_k = cov_matLs[xk];
    arma::mat C = chol(cov_mat_k, "lower");
    arma::mat eps_corr = trans(C*trans(eps));
    epsLs[xk] = eps_corr;
  }
  // initialize assuming series begins in state 1 (use burnin to reduce dependence on this assumption)
  int state = 0;
  // simulate process
  for (int xt = 0; xt<(Tsize+burnin); xt++){
    // Get new state
    // Inverse CDF state selection (faster than find())
    double u_state = use_predrawn_states ? state_rand_vec(xt) : arma::randu<double>();
    int new_state = 0;
    double cum_p = P(0, state);
    while (u_state > cum_p && new_state < k - 1) { new_state++; cum_p += P(new_state, state); }
    state = new_state;
    state_series(xt) = state;
    arma::mat eps_corr_k = epsLs[state];
    resid.row(xt) = eps_corr_k.row(xt);
    Y.row(xt) = mu.row(state) + resid.row(xt);
    if (xt>=burnin){
      Y.row(xt) = Y.row(xt) + Zdm.row(xt-burnin)*beta0;
    }
    mu_t.row(xt) = mu.row(state);
    sigma_t[xt] = cov_matLs[state];
  }
  arma::mat Y_out = Y.submat(burnin,0,Tsize+burnin-1,q-1);
  arma::mat resid_out = resid.submat(burnin,0,Tsize+burnin-1,q-1);
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
//' @description This function computes the log-likelihood for a normally distributed model.
//' 
//' @param theta Vector of model parameters.
//' @param mdl List with model attributes.
//' 
//' @return Log-likelihood value.
//' 
//' @keywords internal
//' 
//' @export
// [[Rcpp::export]]
double logLike_Nmdl(arma::vec theta, List mdl){
  // ---------- model parameters
  arma::mat y = mdl["y"];
  arma::mat xz = mdl["X"];
  arma::uvec beta_ind = arma::find(as<arma::uvec>(mdl["theta_beta_ind"])==1);
  arma::uvec sig_ind = arma::find(as<arma::uvec>(mdl["theta_sig_ind"])==1);
  int Tsize = y.n_rows;
  int q = y.n_cols;
  int npar = xz.n_cols;
  // ---------- pre-define variables
  arma::vec betavec = theta.elem(beta_ind);
  arma::mat beta = trans(reshape(betavec, q, npar));
  arma::vec sig = theta.elem(sig_ind);
  arma::mat sigma = covar_unvech(sig, q);
  // ---------- Compute log-likehood
  arma::mat resid = y-xz*beta;
  double pi = arma::datum::pi;
  arma::vec f_t(Tsize, arma::fill::zeros);
  for (int xt = 0; xt<Tsize; xt++){
    f_t(xt) = arma::as_scalar((1.0/sqrt(det(sigma)*pow(2.0*pi,q)))*
      exp(-0.5*(resid.row(xt)*inv(sigma)*trans(resid.row(xt)))));
  }
  double logLike = sum(log(f_t));
  return(logLike);  
}


// ==============================================================================
//' @title Autoregressive log-likelihood objective function
//' 
//' @description This function computes the log-likelihood for an autoregressive model.
//' 
//' @param theta Vector of model parameters.
//' @param mdl List with model attributes.
//' 
//' @return Log-likelihood value.
//' 
//' @keywords internal
//' 
//' @export
// [[Rcpp::export]]
double logLike_ARmdl(arma::vec theta, List mdl){
  // ---------- Initialize parameters
  arma::vec y = mdl["y"];
  arma::mat x = mdl["x"];
  arma::uvec mu_ind   = arma::find(as<arma::uvec>(mdl["theta_mu_ind"])==1);
  arma::uvec phi_ind  = arma::find(as<arma::uvec>(mdl["theta_phi_ind"])==1);
  arma::uvec sig_ind  = arma::find(as<arma::uvec>(mdl["theta_sig_ind"])==1);
  double mu     = arma::as_scalar(theta.elem(mu_ind));
  arma::vec phi = theta.elem(phi_ind);
  double sigma  = arma::as_scalar(theta.elem(sig_ind));
  int Tsize     = y.n_elem;
  // ---------- Compute log-likehood
  double logLike;
  double pi = arma::datum::pi;
  arma::mat repmu(Tsize, phi.n_elem,arma::fill::ones);
  logLike = sum(log((1.0/sqrt(2.0*pi*sigma))*
    exp(-pow((y - mu) - (x-(mu*repmu))*phi,2.0)/(2.0*sigma))));
  return(logLike);
}




// ==============================================================================
//' @title ARX log-likelihood objective function
//' 
//' @description This function computes the log-likelihood for an autoregressive model with exogenous regressors.
//' 
//' @param theta Vector of model parameters.
//' @param mdl List with model attributes.
//' 
//' @return Log-likelihood value.
//' 
//' @keywords internal
//' 
//' @export
// [[Rcpp::export]]
double logLike_ARXmdl(arma::vec theta, List mdl){
  // ---------- Initialize parameters
  arma::vec y = mdl["y"];
  arma::mat x = mdl["x"];
  List con    = mdl["control"];
  arma::mat z = con["Z"];
  arma::uvec mu_ind     = arma::find(as<arma::uvec>(mdl["theta_mu_ind"])==1);
  arma::uvec phi_ind    = arma::find(as<arma::uvec>(mdl["theta_phi_ind"])==1);
  arma::uvec betaZ_ind  = arma::find(as<arma::uvec>(mdl["theta_x_ind"])==1);
  arma::uvec sig_ind    = arma::find(as<arma::uvec>(mdl["theta_sig_ind"])==1);
  double mu         = arma::as_scalar(theta.elem(mu_ind));
  arma::vec phi     = theta.elem(phi_ind);
  arma::vec betaZ   = theta.elem(betaZ_ind);
  double sigma      = arma::as_scalar(theta.elem(sig_ind));
  int Tsize         = y.n_elem;
  z = z.rows(phi.n_elem,Tsize+phi.n_elem-1);
  arma::rowvec zbar = arma::mean(z,0);
  // ---------- Compute log-likehood
  double logLike;
  double pi = arma::datum::pi;
  arma::mat repmu(Tsize, phi.n_elem,arma::fill::ones);
  arma::mat repzb(Tsize, 1,arma::fill::ones);
  arma::vec resid = (y - mu) - (x-(mu*repmu))*phi - (z-(repzb*zbar))*betaZ;
  logLike = arma::as_scalar(sum(log((1.0/sqrt(2.0*pi*sigma))*exp(-pow(resid,2.0)/(2.0*sigma)))));
  return(logLike);
}
// ==============================================================================
//' @title Vector autoregressive log-likelihood objective function 
//' 
//' @description This function computes the log-likelihood for a vector autoregressive model.
//' 
//' @param theta Vector of model parameters.
//' @param mdl List with model attributes.
//' 
//' @return Log-likelihood value.
//' 
//' @keywords internal
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
  arma::uvec mu_ind   = arma::find(as<arma::uvec>(mdl["theta_mu_ind"])==1);
  arma::uvec phi_ind  = arma::find(as<arma::uvec>(mdl["theta_phi_ind"])==1);
  arma::uvec sig_ind  = arma::find(as<arma::uvec>(mdl["theta_sig_ind"])==1);
  arma::vec mu        = theta.elem(mu_ind);
  arma::vec phi_tmp   = theta.elem(phi_ind);
  arma::mat phi       = reshape(phi_tmp, q*p, q);
  arma::vec sigma_tmp = theta.elem(sig_ind);
  arma::mat sigma     = covar_unvech(sigma_tmp, q);
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
    f_t(xt) = arma::as_scalar((1.0/sqrt(det(sigma)*pow(2.0*pi,q)))*
      exp(-0.5*(resid.row(xt)*inv(sigma)*trans(resid.row(xt)))));
  }
  double logLike = sum(log(f_t));
  return(logLike);  
}
// ==============================================================================
//' @title VARX log-likelihood objective function 
//' 
//' @description This function computes the log-likelihood for a VARX model.
//' 
//' @param theta Vector of model parameters.
//' @param mdl List with model attributes.
//' 
//' @return Log-likelihood value.
//' 
//' @keywords internal
//' 
//' @export
// [[Rcpp::export]]
double logLike_VARXmdl(arma::vec theta, List mdl){
   // ---------- model parameters
   arma::mat y = mdl["y"];
   arma::mat x = mdl["x"];
   List con    = mdl["control"];
   arma::mat z = con["Z"];
   int p = mdl["p"];
   int Tsize = y.n_rows;
   int q = y.n_cols;
   arma::uvec mu_ind    = arma::find(as<arma::uvec>(mdl["theta_mu_ind"])==1);
   arma::uvec phi_ind   = arma::find(as<arma::uvec>(mdl["theta_phi_ind"])==1);
   arma::uvec betaZ_ind = arma::find(as<arma::uvec>(mdl["theta_x_ind"])==1);
   arma::uvec sig_ind   = arma::find(as<arma::uvec>(mdl["theta_sig_ind"])==1);
   arma::vec mu         = theta.elem(mu_ind);
   arma::vec phi_tmp    = theta.elem(phi_ind);
   arma::mat phi        = reshape(phi_tmp, q*p, q);
   arma::vec betaZ_tmp  = theta.elem(betaZ_ind);
   arma::mat betaZ      = reshape(betaZ_tmp,z.n_cols,q);
   arma::vec sigma_tmp  = theta.elem(sig_ind);
   arma::mat sigma      = covar_unvech(sigma_tmp, q);
   arma::mat repmu(Tsize, 1, arma::fill::ones);
   arma::mat ytild      = y - repmu*trans(mu);
   z = z.rows(p,Tsize+p-1);
   arma::rowvec zbar    = arma::mean(z,0);
   arma::mat repzb(Tsize, 1,arma::fill::ones);
   arma::mat xz_tmp(Tsize, q*p, arma::fill::zeros); 
   // ---------- compute residual
   for (int xp = 0; xp<p; xp++){
     xz_tmp.submat(0,q*xp,Tsize-1,q*xp+q-1) = x.submat(0,q*xp,Tsize-1,q*xp+q-1) - repmu*trans(mu);
   }
   arma::mat resid = ytild - xz_tmp*phi - (z-(repzb*zbar))*betaZ;
   // ---------- Compute log-likehood
   double pi = arma::datum::pi;
   arma::vec f_t(Tsize, arma::fill::zeros);
   for (int xt = 0; xt<Tsize; xt++){
     f_t(xt) = arma::as_scalar((1.0/sqrt(det(sigma)*pow(2.0*pi,q)))*
       exp(-0.5*(resid.row(xt)*inv(sigma)*trans(resid.row(xt)))));
   }
   double logLike = sum(log(f_t));
   return(logLike);  
 }
// ==============================================================================
//' @title Hidden Markov model log-likelihood function 
//' 
//' @description This function computes the log-likelihood for a markov-switching autoregressive model.
//' 
//' @param theta Vector of model parameters.
//' @param mdl List with model attributes.
//' @param k Integer determining the number of regimes.
//'  
//' @return Log-likelihood value.
//' 
//' @keywords internal
//' 
//' @export
// [[Rcpp::export]]
double logLike_HMmdl(arma::vec theta, List mdl, int k){
  // ---------- Initialize parameters
  arma::mat y = mdl["y"];
  int q       = y.n_cols;
  int Tsize   = y.n_rows;
  bool msmu   = mdl["msmu"];
  bool msvar  = mdl["msvar"];
  bool exog   = mdl["exog"];
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
  // ----- parmas of exogenous regressors (if present)
  int qz = 0;
  arma::mat Zdm(Tsize,1,arma::fill::zeros);
  arma::mat betaZ(1,q,arma::fill::zeros);
  if (exog==TRUE){
    arma::mat Z         = as<arma::mat>(mdl["Z"]);
    arma::rowvec zbar   = arma::mean(Z,0);
    arma::mat repmu(Tsize, 1, arma::fill::ones);
    Zdm = Z - (repmu*zbar);
    qz                  = Z.n_cols;
    arma::vec betaZtmp  = theta.subvec(q + q*msmu*(k-1), q + q*msmu*(k-1) + (qz*q) - 1);
    betaZ               = reshape(betaZtmp,qz,q); 
  }
  // ----- Variance for each regime 
  int sigN = (q*(q+1))/2;
  arma::vec sig = theta.subvec(q + q*msmu*(k-1) + (qz*q),q + q*msmu*(k-1) + (qz*q) + sigN + (sigN*msvar*(k-1)) - 1);
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
    arma::mat eps_tmp = y - repmu*mu_k.row(xk) - Zdm*betaZ;
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
      //eta(xt,xk) = arma::as_scalar((1.0/(sqrt(pow(2.0*pi,q)*det(sigma_k))))*exp(-0.5*(eps_k.row(xt)*inv(sigma_k)*trans(eps_k.row(xt)))));
      eta(xt,xk) = arma::as_scalar((1.0/(sqrt(pow(2.0*pi,q)*det(sigma_k))))*exp(-0.5*(eps_k.row(xt)*solve(sigma_k,trans(eps_k.row(xt)), arma::solve_opts::allow_ugly))));
    }
    arma::vec xi_eta = xi_t_tm1%trans(eta.row(xt));
    f_t.row(xt) = sum(xi_eta);
    double f_tmp = sum(xi_eta);
    if (f_tmp > 0 && std::isfinite(f_tmp)) {
      xi_t_t.row(xt) = trans(xi_eta/f_tmp);
      xi_t_tm1 = P*trans(xi_t_t.row(xt));
    } else {
      // Underflow: all regime densities collapsed to zero.
      // Reset to ergodic distribution to break NaN cascade.
      // f_t stores true 0; log(0)=-Inf makes logL=-Inf (theta ruled out).
      xi_t_t.row(xt) = trans(pinf);
      xi_t_tm1 = pinf;
    }
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
//' @description This function computes the (negative) log-likelihood for a markov-switching autoregressive model.
//' 
//' @param theta Vector of model parameters.
//' @param mdl List with model attributes.
//' @param k integer determining the number of regimes.
//' 
//' @return Negative log-likelihood value.
//' 
//' @keywords internal
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
//' @description This function computes the log-likelihood for a markov-switching autoregressive model.
//' 
//' @param theta Vector of model parameters.
//' @param mdl List with model attributes.
//' @param k Integer determining the number of regimes.
//' 
//' @return Log-likelihood value.
//' 
//' @keywords internal
//' 
//' @export
// [[Rcpp::export]]
double logLike_MSARmdl(arma::vec theta, List mdl, int k){
  // ---------- Initialize parameters
  arma::vec y     = mdl["y"];
  arma::mat x     = mdl["x"];
  int ar          = mdl["p"];
  bool msmu       = mdl["msmu"];
  bool msvar      = mdl["msvar"];
  int Tsize       = y.n_elem;
  int M           = pow(k, ar+1);
  // ----- Mean for each regime 
  arma::vec mu    = theta.subvec(0, msmu*(k-1));
  // ----- Phi vector
  arma::vec phi   =  theta.subvec(1+msmu*(k-1), 1+msmu*(k-1)+ar-1);
  // ----- Variance for each regime 
  arma::vec sig   = theta.subvec(1+msmu*(k-1)+ar,1+msmu*(k-1)+ar+msvar*(k-1));
  // ----- Transition probabilities 
  arma::mat P     = reshape(theta.subvec(2+msmu*(k-1)+msvar*(k-1) + ar, 2+msmu*(k-1)+msvar*(k-1) + ar + k*k - 1),k, k);
  // Regime limiting probabilities
  arma::vec pinf  = limP(P);
  // ----- Obtain AR consistent grid of Mu, Sigma and State indicators
  List musig_out  = argrid_MSARmdl_cpp(mu, sig, k, ar);
  arma::mat muAR  = as<arma::mat>(musig_out["mu"]);
  arma::mat sigAR = as<arma::mat>(musig_out["sig"]);
  arma::vec state_ind = as<arma::vec>(musig_out["state_ind"]);
  // ----- Obtain AR consistent P and pinf
  arma::mat P_AR = arP_cpp(P, k, ar);
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
    eta.row(xt) = trans(1.0/sqrt(2.0*pi*(sigAR)))%exp((pow(eps.row(xt),2.0)%trans(-1.0/(2.0*sigAR))));
    arma::vec xi_eta = xi_t_tm1%trans(eta.row(xt));
    f_t.row(xt) = sum(xi_eta);
    double f_tmp = sum(xi_eta);
    if (f_tmp > 0 && std::isfinite(f_tmp)) {
      xi_t_t.row(xt) = trans(xi_eta/f_tmp);
      xi_t_tm1 = predict_fast(trans(xi_t_t.row(xt)), P, k, ar);
    } else {
      // Underflow: all regime densities collapsed to zero.
      // Reset to ergodic distribution to break NaN cascade.
      // f_t stores true 0; log(0)=-Inf makes logL=-Inf (theta ruled out).
      xi_t_t.row(xt) = trans(pinf_AR);
      xi_t_tm1 = pinf_AR;
    }
    xi_tp1_t.row(xt) = trans(xi_t_tm1);
  }
  // log-Likelihood
  arma::vec logf = log(f_t);
  double logL = sum(logf);     //[eq. 22.4.7]
  return(logL);
}

// ==============================================================================
//' @title Markov-switching ARX log-likelihood objective function
//' 
//' @description This function computes the log-likelihood for a markov-switching ARX model.
//' 
//' @param theta Vector of model parameters.
//' @param mdl List with model attributes.
//' @param k Integer determining the number of regimes.
//' 
//' @return Log-likelihood value.
//' 
//' @keywords internal
//' 
//' @export
// [[Rcpp::export]]
double logLike_MSARXmdl(arma::vec theta, List mdl, int k){
  // ---------- Initialize parameters
  arma::vec y     = mdl["y"];
  arma::mat x     = mdl["x"];
  int ar          = mdl["p"];
  bool msmu       = mdl["msmu"];
  bool msvar      = mdl["msvar"];
  List mdl_con    = mdl["control"];
  arma::mat Z     = mdl_con["Z"];
  int Tsize       = y.n_elem;
  int M           = pow(k, ar+1);
  int qz          = Z.n_cols;
  Z               = Z.rows(ar,Tsize+ar-1);
  // ----- Mean for each regime 
  arma::vec mu          = theta.subvec(0, msmu*(k-1));
  // ----- Phi vector
  arma::vec phi         = theta.subvec(1+msmu*(k-1), 1+msmu*(k-1)+ar-1);
  // ----- betaZ vector
  arma::vec betaZ       = theta.subvec(1+msmu*(k-1)+ar, 1+msmu*(k-1)+ar+qz-1);
  // ----- Variance for each regime 
  arma::vec sig         = theta.subvec(1+msmu*(k-1)+ar+qz,1+msmu*(k-1)+ar+qz+msvar*(k-1));
  // ----- Transition probabilities 
  arma::mat P           = reshape(theta.subvec(2+msmu*(k-1)+ar+qz+msvar*(k-1), 2+msmu*(k-1)+ar+qz+msvar*(k-1) + k*k - 1),k, k);
  // Regime limiting probabilities
  arma::vec pinf        = limP(P);
  // ----- Obtain AR consistent grid of Mu, Sigma and State indicators
  List musig_out        = argrid_MSARmdl_cpp(mu, sig, k, ar);
  arma::mat muAR        = as<arma::mat>(musig_out["mu"]);
  arma::mat sigAR       = as<arma::mat>(musig_out["sig"]);
  arma::vec state_ind   = as<arma::vec>(musig_out["state_ind"]);
  // ----- Obtain AR consistent P and pinf
  arma::mat P_AR        = arP_cpp(P, k, ar);
  arma::mat pinf_AR     = limP(P_AR);
  // ----- Compute residuals in each regime
  arma::rowvec zbar  = arma::mean(Z,0);
  arma::mat repvec(1, M, arma::fill::ones);
  arma::mat repmu(Tsize, 1, arma::fill::ones);
  arma::mat ms_y     = y*repvec;
  arma::mat ms_zdm   = ((Z-repmu*zbar)*betaZ)*repvec;
  arma::mat eps(Tsize, M, arma::fill::zeros);
  arma::mat ms_ydm = ms_y - repmu*trans(muAR.col(0)); // [y(t) - mu_s(t))]
  arma::mat ms_xdm(Tsize, M, arma::fill::zeros);
  for (int xkp = 0; xkp<M; xkp++){
    arma::mat zx_tmp = x - repmu*muAR.submat(xkp,1,xkp,ar); // [y(t-i) - mu_s(t-i))]
    ms_xdm.col(xkp) = zx_tmp*phi;
  }
  eps = ms_ydm - ms_xdm - ms_zdm;
  // ----- Begin Calculating likelihood and log-likelihood
  arma::mat eta(Tsize, M, arma::fill::zeros);       // [eq. 22.4.2]
  arma::mat f_t(Tsize, 1, arma::fill::zeros);       // [eq. 22.4.8]
  arma::mat xi_t_t(Tsize, M, arma::fill::zeros);    // [eq. 22.4.5]
  arma::mat xi_tp1_t(Tsize, M, arma::fill::zeros);  // [eq. 22.4.6]
  double pi = arma::datum::pi;
  arma::vec xi_t_tm1 = pinf_AR;
  for (int xt = 0; xt<Tsize; xt++){
    eta.row(xt) = trans(1.0/sqrt(2.0*pi*(sigAR)))%exp((pow(eps.row(xt),2.0)%trans(-1.0/(2.0*sigAR))));
    arma::vec xi_eta = xi_t_tm1%trans(eta.row(xt));
    f_t.row(xt) = sum(xi_eta);
    double f_tmp = sum(xi_eta);
    if (f_tmp > 0 && std::isfinite(f_tmp)) {
      xi_t_t.row(xt) = trans(xi_eta/f_tmp);
      xi_t_tm1 = predict_fast(trans(xi_t_t.row(xt)), P, k, ar);
    } else {
      // Underflow: all regime densities collapsed to zero.
      // Reset to ergodic distribution to break NaN cascade.
      // f_t stores true 0; log(0)=-Inf makes logL=-Inf (theta ruled out).
      xi_t_t.row(xt) = trans(pinf_AR);
      xi_t_tm1 = pinf_AR;
    }
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
//' @description This function computes the (negative) log-likelihood for a markov-switching autoregressive model.
//' 
//' @param theta Vector of model parameters.
//' @param mdl List with model attributes.
//' @param k integer determining the number of regimes.
//' 
//' @return Negative log-likelihood value.
//' 
//' @keywords internal
//' 
//' @export
// [[Rcpp::export]]
double logLike_MSARmdl_min(arma::vec theta, List mdl, int k){
  double logL_negative = -logLike_MSARmdl(theta, mdl, k);
  return(logL_negative);
}

// ==============================================================================
//' @title Markov-switching ARX log-likelihood objective function (minimization version)
//' 
//' @description This function computes the (negative) log-likelihood for a markov-switching ARX model.
//' 
//' @param theta Vector of model parameters.
//' @param mdl List with model attributes.
//' @param k integer determining the number of regimes.
//' 
//' @return Negative log-likelihood value.
//' 
//' @keywords internal
//' 
//' @export
// [[Rcpp::export]]
double logLike_MSARXmdl_min(arma::vec theta, List mdl, int k){
   double logL_negative = -logLike_MSARXmdl(theta, mdl, k);
   return(logL_negative);
}


// ==============================================================================
//' @title Markov-switching vector autoregressive log-likelihood objective function
//' 
//' @description This function computes the log-likelihood for a markov-switching vector autoregressive model.
//' 
//' @param theta Vector of model parameters.
//' @param mdl List with model attributes.
//' @param k Integer determining the number of regimes.
//' 
//' @return Log-likelihood value.
//' 
//' @keywords internal
//' 
//' @export
// [[Rcpp::export]]
double logLike_MSVARmdl(arma::vec theta, List mdl, int k){
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
  // ----- Phi vector
  arma::vec phi_tmp =  theta.subvec(q+q*msmu*(k-1), q+q*msmu*(k-1) + q*q*ar - 1);
  arma::mat phi = trans(reshape(phi_tmp, q*ar, q));
  // ----- Variance for each regime 
  int sigN = (q*(q+1))/2;
  arma::vec sig = theta.subvec(q+q*msmu*(k-1) + (q*q*ar) ,q+q*msmu*(k-1) + q*q*ar + (sigN + sigN*msvar*(k-1)) - 1);
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
  int PN = theta.n_elem - (k*k);
  arma::mat P = reshape(theta.subvec(PN, PN + k*k - 1), k, k);
  // Regime limiting probabilities
  arma::vec pinf = limP(P);
  // ----- Obtain AR consistent grid of Mu, Sigma and State indicators
  List musig_out = argrid_MSVARmdl_cpp(mu_k, sigma, k, ar);
  List muAR = musig_out["mu"];
  List sigAR = musig_out["sig"];
  arma::vec state_ind = as<arma::vec>(musig_out["state_ind"]);
  // ----- Obtain AR consistent P and pinf
  arma::mat PAR = arP_cpp(P, k, ar);
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
    eps[xm] = y_tmp - xz_tmp*trans(phi);
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
  // Pre-compute per-regime density constants (sigAR constant during T-loop)
  std::vector<arma::mat> inv_sigma_precomp(M);
  arma::vec norm_const_precomp(M);
  for (int xm_pre = 0; xm_pre < M; xm_pre++) {
    arma::mat sm = sigAR[xm_pre];
    double det_val = arma::det(sm);
    norm_const_precomp(xm_pre) = (det_val > 0) ? 1.0 / sqrt(pow(2.0*pi, q) * det_val) : 0.0;
    arma::mat inv_sm;
    bool ok = arma::inv(inv_sm, sm);
    inv_sigma_precomp[xm_pre] = ok ? inv_sm : arma::zeros<arma::mat>(q, q);
  }
  arma::vec xi_t_tm1_AR = pinfAR;
  arma::vec xi_t_tm1 = pinf;
  for (int xt = 0; xt<Tsize; xt++){
    for (int xm = 0; xm<M; xm++){
      arma::mat eps_m = eps[xm];
      arma::vec e_vec = trans(eps_m.row(xt));
      eta(xt,xm) = norm_const_precomp(xm) * exp(-0.5 * arma::dot(e_vec, inv_sigma_precomp[xm] * e_vec));
    }
    arma::vec xi_eta = xi_t_tm1_AR%trans(eta.row(xt));
    f_t.row(xt) = sum(xi_eta);
    double f_tmp = sum(xi_eta);
    if (f_tmp > 0 && std::isfinite(f_tmp)) {
      xi_t_t_tmp.row(xt) = trans(xi_eta/f_tmp);
    } else {
      // Underflow: all regime densities collapsed to zero.
      // Reset to ergodic distribution to break NaN cascade.
      // f_t stores true 0; log(0)=-Inf makes logL=-Inf (theta ruled out).
      xi_t_t_tmp.row(xt) = trans(pinfAR);
    }
    for (int xk = 1; xk<=k; xk++){
      arma::vec xi_t_t_row = trans(xi_t_t_tmp.row(xt));
      xi_t_t.submat(xt,xk-1,xt,xk-1) = sum(xi_t_t_row.rows(find(state_ind==xk)));
    }
    xi_t_tm1_AR = predict_fast(trans(xi_t_t_tmp.row(xt)), P, k, ar);
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
//' @title Markov-switching VARX log-likelihood objective function
//' 
//' @description This function computes the log-likelihood for a markov-switching VARX model.
//' 
//' @param theta Vector of model parameters.
//' @param mdl List with model attributes.
//' @param k Integer determining the number of regimes.
//' 
//' @return Log-likelihood value.
//' 
//' @keywords internal
//' 
//' @export
// [[Rcpp::export]]
double logLike_MSVARXmdl(arma::vec theta, List mdl, int k){
  // ---------- Initialize parameters
  arma::mat y   = mdl["y"];
  arma::mat x   = mdl["x"];
  int ar        = mdl["p"];
  bool msmu     = mdl["msmu"];
  bool msvar    = mdl["msvar"];
  List mdl_con  = mdl["control"];
  arma::mat Z   = mdl_con["Z"];
  int q         = y.n_cols;
  int Tsize     = y.n_rows;
  int M         = pow(k, ar+1);
  int qz        = Z.n_cols;
  Z             = Z.rows(ar,Tsize+ar-1);
  // ----- Mean for each regime 
  arma::mat mu_k(k, q, arma::fill::zeros);
  arma::vec mu = theta.subvec(0, q + q*msmu*(k-1)-1);
  if (msmu==TRUE){
    for (int xk = 0; xk<k; xk++){
      mu_k.row(xk) = trans(mu.subvec(xk*q,xk*q+q-1));
    }
  }else{
    for (int xk = 0; xk<k; xk++){
      mu_k.row(xk) = trans(mu);
    }
  } 
  // ----- Phi vector
  arma::vec phi_tmp =  theta.subvec(q + q*msmu*(k-1), q + q*msmu*(k-1) + (q*q*ar) - 1);
  arma::mat phi = trans(reshape(phi_tmp, q*ar, q));
  // ----- betaZ vector
  arma::vec betaZ_tmp  = theta.subvec(q + q*msmu*(k-1) + (q*q*ar), q + q*msmu*(k-1) + (q*q*ar) + qz*q - 1);
  arma::mat betaZ = reshape(betaZ_tmp,qz,q);
  // ----- Variance for each regime 
  int sigN = (q*(q+1))/2;
  arma::vec sig = theta.subvec(q + q*msmu*(k-1) + (q*q*ar) + qz*q,q + q*msmu*(k-1) + (q*q*ar) + qz*q + (sigN + sigN*msvar*(k-1)) - 1);
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
  int PN = theta.n_elem - (k*k);
  arma::mat P = reshape(theta.subvec(PN, PN + k*k - 1), k, k);
  // Regime limiting probabilities
  arma::vec pinf = limP(P);
  // ----- Obtain AR consistent grid of Mu, Sigma and State indicators
  List musig_out = argrid_MSVARmdl_cpp(mu_k, sigma, k, ar);
  List muAR = musig_out["mu"];
  List sigAR = musig_out["sig"];
  arma::vec state_ind = as<arma::vec>(musig_out["state_ind"]);
  // ----- Obtain AR consistent P and pinf
  arma::mat PAR = arP_cpp(P, k, ar);
  arma::mat pinfAR = limP(PAR);
  // ----- Compute residuals in each regime
  arma::rowvec zbar  = arma::mean(Z,0);
  arma::mat repmu(Tsize, 1, arma::fill::ones);
  arma::mat ms_zdm   = ((Z-repmu*zbar)*betaZ);
  List eps(M); 
  for (int xm = 0; xm<M; xm++){
    arma::mat mu_tmp  = muAR[xm]; 
    arma::mat y_tmp   = y - repmu*trans(mu_tmp.col(0));
    arma::mat xz_tmp(Tsize, q*ar, arma::fill::zeros); 
    for (int xp = 0; xp<ar; xp++){
      xz_tmp.submat(0,q*xp,Tsize-1,q*xp+q-1) = x.submat(0,q*xp,Tsize-1,q*xp+q-1) - repmu*trans(mu_tmp.col(xp+1));
    }
    eps[xm] = y_tmp - xz_tmp*trans(phi) - ms_zdm;
  }
  // ----- Begin Calculating log-likelihood & filtered probabilities
  arma::mat eta(Tsize, M, arma::fill::zeros);       // [eq. 22.4.2]
  arma::mat f_t(Tsize, 1, arma::fill::zeros);       // [eq. 22.4.8]
  arma::mat xi_t_t(Tsize, k, arma::fill::zeros);    // [eq. 22.4.5]
  arma::mat xi_t_t_tmp(Tsize, M, arma::fill::zeros);    
  arma::mat xi_tp1_t(Tsize, k, arma::fill::zeros);  // [eq. 22.4.6]
  arma::mat xi_tp1_t_tmp(Tsize, M, arma::fill::zeros);  
  double pi = arma::datum::pi;
  // Pre-compute per-regime density constants (sigAR constant during T-loop)
  std::vector<arma::mat> inv_sigma_precomp(M);
  arma::vec norm_const_precomp(M);
  for (int xm_pre = 0; xm_pre < M; xm_pre++) {
    arma::mat sm = sigAR[xm_pre];
    double det_val = arma::det(sm);
    norm_const_precomp(xm_pre) = (det_val > 0) ? 1.0 / sqrt(pow(2.0*pi, q) * det_val) : 0.0;
    arma::mat inv_sm;
    bool ok = arma::inv(inv_sm, sm);
    inv_sigma_precomp[xm_pre] = ok ? inv_sm : arma::zeros<arma::mat>(q, q);
  }
  arma::vec xi_t_tm1_AR = pinfAR;
  arma::vec xi_t_tm1 = pinf;
  for (int xt = 0; xt<Tsize; xt++){
    for (int xm = 0; xm<M; xm++){
      arma::mat eps_m = eps[xm];
      arma::vec e_vec = trans(eps_m.row(xt));
      eta(xt,xm) = norm_const_precomp(xm) * exp(-0.5 * arma::dot(e_vec, inv_sigma_precomp[xm] * e_vec));
    } 
    arma::vec xi_eta = xi_t_tm1_AR%trans(eta.row(xt));
    f_t.row(xt) = sum(xi_eta);
    double f_tmp = sum(xi_eta);
    if (f_tmp > 0 && std::isfinite(f_tmp)) {
      xi_t_t_tmp.row(xt) = trans(xi_eta/f_tmp);
    } else {
      // Underflow: all regime densities collapsed to zero.
      // Reset to ergodic distribution to break NaN cascade.
      // f_t stores true 0; log(0)=-Inf makes logL=-Inf (theta ruled out).
      xi_t_t_tmp.row(xt) = trans(pinfAR);
    }
    for (int xk = 1; xk<=k; xk++){
      arma::vec xi_t_t_row = trans(xi_t_t_tmp.row(xt));
      xi_t_t.submat(xt,xk-1,xt,xk-1) = sum(xi_t_t_row.rows(find(state_ind==xk)));
    } 
    xi_t_tm1_AR = predict_fast(trans(xi_t_t_tmp.row(xt)), P, k, ar);
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
//' @param theta Vector of model parameters.
//' @param mdl List with model attributes.
//' @param k Integer determining the number of regimes.
//' 
//' @return Negative log-likelihood value.
//' 
//' @keywords internal
//' 
//' @export
// [[Rcpp::export]]
double logLike_MSVARmdl_min(arma::vec theta, List mdl, int k){
  double logL_negative = -logLike_MSVARmdl(theta, mdl, k);
  return(logL_negative);
}

// ==============================================================================
//' @title Markov-switching VARX log-likelihood objective function (minimization version)
//' 
//' @description This function computes the (negative) log-likelihood for a markov-switching VARX model
//' 
//' @param theta Vector of model parameters.
//' @param mdl List with model attributes.
//' @param k Integer determining the number of regimes.
//' 
//' @return Negative log-likelihood value.
//' 
//' @keywords internal
//' 
//' @export
// [[Rcpp::export]]
double logLike_MSVARXmdl_min(arma::vec theta, List mdl, int k){
  double logL_negative = -logLike_MSVARXmdl(theta, mdl, k);
  return(logL_negative);
}

// ==============================================================================
//' @title Hidden Markov model log-likelihood function 
//' 
//' @description This function computes the log-likelihood for a Hidden Markov model and uses the Hamilton smoother to obtain smoothed probabilities of each state. This is also the expectation step in the Expectation Maximization algorithm for a Markov-switching autoregressive model.
//' 
//' @param theta Vector of model parameters.
//' @param mdl List with model attributes.
//' @param k Integer determining the number of regimes.
//'  
//' @return List which includes log-likelihood value and smoothed probabilities of each regime.
//' 
//' @keywords internal
//' 
//' @export
// [[Rcpp::export]]
List ExpectationM_HMmdl(arma::vec theta, List mdl, int k){
  // ---------- Initialize parameters
  arma::mat y = mdl["y"];
  int q       = y.n_cols;
  int Tsize   = y.n_rows;
  bool msmu   = mdl["msmu"];
  bool msvar  = mdl["msvar"];
  bool exog   = mdl["exog"];
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
  // ----- parmas of exogenous regressors (if present)
  int qz = 0;
  arma::mat Zdm(Tsize,1,arma::fill::zeros);
  arma::mat betaZ(1,q,arma::fill::zeros);
  if (exog==TRUE){
    arma::mat Z         = as<arma::mat>(mdl["Z"]);
    arma::rowvec zbar  = arma::mean(Z,0);
    arma::mat repmu(Tsize, 1, arma::fill::ones);
    Zdm = Z - (repmu*zbar);
    qz                  = Z.n_cols;
    arma::vec betaZtmp  = theta.subvec(q + q*msmu*(k-1), q + q*msmu*(k-1) + (qz*q) - 1);
    betaZ               = reshape(betaZtmp,qz,q); 
  }
  // ----- Variance for each regime 
  int sigN = (q*(q+1))/2;
  arma::vec sig = theta.subvec(q + q*msmu*(k-1) + (qz*q),q + q*msmu*(k-1) + (qz*q) + sigN + (sigN*msvar*(k-1)) - 1);
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
    arma::mat eps_tmp = y - repmu*mu_k.row(xk) - Zdm*betaZ;
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
      //eta(xt,xk) = arma::as_scalar((1.0/(sqrt(pow(2.0*pi,q)*det(sigma_k))))*exp(-0.5*(eps_k.row(xt)*inv(sigma_k)*trans(eps_k.row(xt)))));
      eta(xt,xk) = arma::as_scalar((1.0/(sqrt(pow(2.0*pi,q)*det(sigma_k))))*exp(-0.5*(eps_k.row(xt)*solve(sigma_k,trans(eps_k.row(xt)), arma::solve_opts::allow_ugly))));
    }
    arma::vec xi_eta = xi_t_tm1%trans(eta.row(xt));
    f_t.row(xt) = sum(xi_eta);
    double f_tmp = sum(xi_eta);
    if (f_tmp > 0 && std::isfinite(f_tmp)) {
      xi_t_t.row(xt) = trans(xi_eta/f_tmp);
      xi_t_tm1 = P*trans(xi_t_t.row(xt));
    } else {
      // Underflow: all regime densities collapsed to zero.
      // Reset to ergodic distribution to break NaN cascade.
      // f_t stores true 0; log(0)=-Inf makes logL=-Inf (theta ruled out).
      xi_t_t.row(xt) = trans(pinf);
      xi_t_tm1 = pinf;
    }
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
    arma::rowvec xi_tp1_safe = xi_tp1_t.row(xT);
    xi_tp1_safe.clamp(1e-300, 1.0);  // guard against exact-zero denominator
    xi_t_T.row(xT) = xi_t_t.row(xT)%trans(trans(P)*trans(xi_t_T.row(xT+1)/xi_tp1_safe));
    // Non-finite guard: if smoother fails (should never happen for k-state), use filtered.
    if (!xi_t_T.row(xT).is_finite()) {
      xi_t_T.row(xT) = xi_t_t.row(xT);
    }
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
  MSloglik_output["mu"] = mu_k;
  MSloglik_output["betaZ"] = betaZ;
  MSloglik_output["sigma"] = sigma;
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
//' @param theta Vector of model parameters.
//' @param mdl List with model attributes.
//' @param k Integer determining the number of regimes.
//'  
//' @return List which includes log-likelihood and smoothed probabilities of each regime.
//' 
//' @keywords internal
//' 
//' @export
// [[Rcpp::export]]
List ExpectationM_MSARmdl(arma::vec theta, List mdl, int k){
  // ---------- Initialize parameters
  arma::vec y = mdl["y"];
  int Tsize   = y.n_elem;
  int p       = mdl["p"];
  bool msmu   = mdl["msmu"];
  bool msvar  = mdl["msvar"];
  int M       = pow(k, p+1);
  List pList  = paramList_MSARmdl(theta, p, k, msmu, msvar);
  arma::mat mu        = pList["mu"];
  arma::mat sig       = pList["sig"];
  arma::mat P         = pList["P"];
  arma::vec pinf      = pList["pinf"];
  arma::mat muAR      = pList["muAR"];
  arma::mat sigAR     = pList["sigAR"];
  arma::vec phi       = pList["phi"];
  arma::mat PAR       = pList["P_AR"];
  arma::vec pinfAR    = pList["pinf_AR"];
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
    eta.row(xt) = trans(1.0/sqrt(2.0*pi*(sigAR)))%exp((pow(eps.row(xt),2.0)%trans(-1.0/(2.0*sigAR))));
    arma::vec xi_eta = xi_t_tm1_AR%trans(eta.row(xt));
    f_t.row(xt) = sum(xi_eta);
    double f_tmp = sum(xi_eta);
    if (f_tmp > 0 && std::isfinite(f_tmp)) {
      xi_t_t_tmp.row(xt) = trans(xi_eta/f_tmp);
    } else {
      // Underflow: all regime densities collapsed to zero.
      // Reset to ergodic distribution to break NaN cascade.
      // f_t stores true 0; log(0)=-Inf makes logL=-Inf (theta ruled out).
      xi_t_t_tmp.row(xt) = trans(pinfAR);
    }
    for (int xk = 1; xk<=k; xk++){
      arma::vec xi_t_t_row = trans(xi_t_t_tmp.row(xt));
      xi_t_t.submat(xt,xk-1,xt,xk-1) = sum(xi_t_t_row.rows(find(state_ind==xk)));
    }
    xi_t_tm1_AR = predict_fast(trans(xi_t_t_tmp.row(xt)), P, k, p);
    xi_t_tm1 = P*trans(xi_t_t.row(xt));
    xi_tp1_t_tmp.row(xt) = trans(xi_t_tm1_AR);
    xi_tp1_t.row(xt) = trans(xi_t_tm1);
  }
  // ----- log-Likelihood
  arma::vec logf = log(f_t);
  double logL = sum(logf);     //[eq. 22.4.7]
  double L = logL/Tsize; 
  // ---------- Hamilton smoother
  arma::mat xi_t_T(Tsize, k, arma::fill::zeros);
  arma::mat xi_t_T_tmp(Tsize, M, arma::fill::zeros);
  xi_t_T_tmp.row(Tsize-1)  = xi_t_t_tmp.row(Tsize-1);
  for (int xT = Tsize-2; xT>=0; xT--){
    arma::rowvec xi_tp1_tmp_safe = xi_tp1_t_tmp.row(xT);
    xi_tp1_tmp_safe.clamp(1e-300, 1.0);  // guard against exact-zero denominator
    xi_t_T_tmp.row(xT) = xi_t_t_tmp.row(xT)%trans(smooth_fast(trans(xi_t_T_tmp.row(xT+1)/xi_tp1_tmp_safe), P, k, p));
    double row_sum_M = accu(xi_t_T_tmp.row(xT));
    if (row_sum_M > 0 && std::isfinite(row_sum_M)) {
      xi_t_T_tmp.row(xT) = xi_t_T_tmp.row(xT) / row_sum_M;
    } else {
      xi_t_T_tmp.row(xT) = xi_t_t_tmp.row(xT);
    }
  }
  // Marginal (k-state) smoothed probabilities: the conditional density of an MS-AR depends
  // on the date-t and lagged regimes, so smoothed inference lives on the extended state
  // (Hamilton 1994, sec. 22.4); the marginal is its aggregate over the date-t regime.
  for (int xk = 1; xk<=k; xk++){
    xi_t_T.col(xk-1) = arma::sum(xi_t_T_tmp.cols(find(state_ind==xk)), 1);
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
  MSloglik_output["logLike"]    = logL;
  MSloglik_output["L"]          = L;
  MSloglik_output["eta"]        = eta;
  MSloglik_output["f_t"]        = f_t;
  MSloglik_output["xi_t_t_AR"]  = xi_t_t_tmp;
  MSloglik_output["xi_tp1_t_AR"]= xi_tp1_t_tmp;
  MSloglik_output["xi_t_T_AR"]  = xi_t_T_tmp;
  MSloglik_output["xi_t_t"]     = xi_t_t;
  MSloglik_output["xi_tp1_t"]   = xi_tp1_t;
  MSloglik_output["xi_t_T"]     = xi_t_T;
  MSloglik_output["PAR"]        = PAR;
  MSloglik_output["P"]          = P;
  MSloglik_output["pinfAR"]     = pinfAR;
  MSloglik_output["pinf"]       = pinf;
  MSloglik_output["muAR"]       = muAR;
  MSloglik_output["sigmaAR"]    = sigAR;
  MSloglik_output["mu"]         = mu;
  MSloglik_output["sigma"]      = sig;
  MSloglik_output["phi"]        = phi;
  MSloglik_output["beta"]       = phi;
  MSloglik_output["theta"]      = theta;
  MSloglik_output["state_ind"]  = state_ind;
  MSloglik_output["residuals"]  = eps;
  MSloglik_output["resid"]      = residuals;
  return(MSloglik_output);
}

// ==============================================================================
//' @title Markov-switching ARX log-likelihood function 
//' 
//' @description This function computes the log-likelihood for a markov-switching autoregressive model and uses the Hamilton smoother to obtain smoothed probabilities of each state. This is also the expectation step in the Expectation Maximization algorithm for a Markov-switching ARX model.
//' 
//' @param theta Vector of model parameters.
//' @param mdl List with model attributes.
//' @param k Integer determining the number of regimes.
//'  
//' @return List which includes log-likelihood and smoothed probabilities of each regime.
//' 
//' @keywords internal
//' 
//' @export
// [[Rcpp::export]]
List ExpectationM_MSARXmdl(arma::vec theta, List mdl, int k){
  // ---------- Initialize parameters
  arma::vec y   = mdl["y"];
  int Tsize     = y.n_elem;
  int p         = mdl["p"];
  bool msmu     = mdl["msmu"];
  bool msvar    = mdl["msvar"];
  int M         = pow(k, p+1);
  List mdl_con  = mdl["control"];
  arma::mat Z   = mdl_con["Z"];
  int qz        = Z.n_cols;
  List pList    = paramList_MSARXmdl(theta, p, k, qz, msmu, msvar);
  arma::mat mu        = pList["mu"];
  arma::vec phi       = pList["phi"];
  arma::vec betaZ     = pList["betaZ"];
  arma::mat sig       = pList["sig"];
  arma::mat P         = pList["P"];
  arma::vec pinf      = pList["pinf"];
  arma::mat muAR      = pList["muAR"];
  arma::mat sigAR     = pList["sigAR"];
  arma::mat PAR       = pList["P_AR"];
  arma::vec pinfAR    = pList["pinf_AR"];
  arma::vec state_ind = pList["state_ind"];
  // ----- Compute residuals in each regime M
  List mdl_tmp        = clone(mdl);
  mdl_tmp["phi"]      = phi;
  mdl_tmp["betaZ"]    = betaZ;
  arma::mat eps       = calcResid_MSARXmdl(mdl_tmp, muAR, k);
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
    eta.row(xt) = trans(1.0/sqrt(2.0*pi*(sigAR)))%exp((pow(eps.row(xt),2.0)%trans(-1.0/(2.0*sigAR))));
    arma::vec xi_eta = xi_t_tm1_AR%trans(eta.row(xt));
    f_t.row(xt) = sum(xi_eta);
    double f_tmp = sum(xi_eta);
    if (f_tmp > 0 && std::isfinite(f_tmp)) {
      xi_t_t_tmp.row(xt) = trans(xi_eta/f_tmp);
    } else {
      // Underflow: all regime densities collapsed to zero.
      // Reset to ergodic distribution to break NaN cascade.
      // f_t stores true 0; log(0)=-Inf makes logL=-Inf (theta ruled out).
      xi_t_t_tmp.row(xt) = trans(pinfAR);
    }
    for (int xk = 1; xk<=k; xk++){
      arma::vec xi_t_t_row = trans(xi_t_t_tmp.row(xt));
      xi_t_t.submat(xt,xk-1,xt,xk-1) = sum(xi_t_t_row.rows(find(state_ind==xk)));
    }
    xi_t_tm1_AR = predict_fast(trans(xi_t_t_tmp.row(xt)), P, k, p);
    xi_t_tm1 = P*trans(xi_t_t.row(xt));
    xi_tp1_t_tmp.row(xt) = trans(xi_t_tm1_AR);
    xi_tp1_t.row(xt) = trans(xi_t_tm1);
  }
  // ----- log-Likelihood
  arma::vec logf = log(f_t);
  double logL = sum(logf);     //[eq. 22.4.7]
  double L = logL/Tsize; 
  // ---------- Hamilton smoother
  arma::mat xi_t_T(Tsize, k, arma::fill::zeros);
  arma::mat xi_t_T_tmp(Tsize, M, arma::fill::zeros);
  xi_t_T_tmp.row(Tsize-1)  = xi_t_t_tmp.row(Tsize-1);
  for (int xT = Tsize-2; xT>=0; xT--){
    arma::rowvec xi_tp1_tmp_safe = xi_tp1_t_tmp.row(xT);
    xi_tp1_tmp_safe.clamp(1e-300, 1.0);  // guard against exact-zero denominator
    xi_t_T_tmp.row(xT) = xi_t_t_tmp.row(xT)%trans(smooth_fast(trans(xi_t_T_tmp.row(xT+1)/xi_tp1_tmp_safe), P, k, p));
    double row_sum_M = accu(xi_t_T_tmp.row(xT));
    if (row_sum_M > 0 && std::isfinite(row_sum_M)) {
      xi_t_T_tmp.row(xT) = xi_t_T_tmp.row(xT) / row_sum_M;
    } else {
      xi_t_T_tmp.row(xT) = xi_t_t_tmp.row(xT);
    }
  }
  // Marginal (k-state) smoothed probabilities: the conditional density of an MS-AR depends
  // on the date-t and lagged regimes, so smoothed inference lives on the extended state
  // (Hamilton 1994, sec. 22.4); the marginal is its aggregate over the date-t regime.
  for (int xk = 1; xk<=k; xk++){
    xi_t_T.col(xk-1) = arma::sum(xi_t_T_tmp.cols(find(state_ind==xk)), 1);
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
  MSloglik_output["logLike"]    = logL;
  MSloglik_output["L"]          = L;
  MSloglik_output["eta"]        = eta;
  MSloglik_output["f_t"]        = f_t;
  MSloglik_output["xi_t_t_AR"]  = xi_t_t_tmp;
  MSloglik_output["xi_tp1_t_AR"]= xi_tp1_t_tmp;
  MSloglik_output["xi_t_T_AR"]  = xi_t_T_tmp;
  MSloglik_output["xi_t_t"]     = xi_t_t;
  MSloglik_output["xi_tp1_t"]   = xi_tp1_t;
  MSloglik_output["xi_t_T"]     = xi_t_T;
  MSloglik_output["PAR"]        = PAR;
  MSloglik_output["P"]          = P;
  MSloglik_output["pinfAR"]     = pinfAR;
  MSloglik_output["pinf"]       = pinf;
  MSloglik_output["muAR"]       = muAR;
  MSloglik_output["sigmaAR"]    = sigAR;
  MSloglik_output["mu"]         = mu;
  MSloglik_output["phi"]        = phi;
  MSloglik_output["betaZ"]      = betaZ;
  MSloglik_output["beta"]       = join_vert(phi,betaZ);
  MSloglik_output["sigma"]      = sig;
  MSloglik_output["theta"]      = theta;
  MSloglik_output["state_ind"]  = state_ind;
  MSloglik_output["residuals"]  = eps;
  MSloglik_output["resid"]      = residuals;
  return(MSloglik_output);
}

// ==============================================================================
//' @title Markov-switching vector autoregressive log-likelihood function
//' 
//' @description This function computes the log-likelihood for a markov-switching vector autoregressive model and uses the Hamilton smoother to obtain smoothed probabilities of each state. This is also the expectation step in the Expectation Maximization algorithm for a Markov-switching autoregressive model.
//' 
//' @param theta Vector of model parameters.
//' @param mdl List with model attributes.
//' @param k Integer determining the number of regimes.
//'  
//' @return List which includes log-likelihood and smoothed probabilities of each regime.
//' 
//' @keywords internal
//' 
//' @export
// [[Rcpp::export]]
List ExpectationM_MSVARmdl(arma::vec theta, List mdl, int k){
  // ---------- Initialize parameters
  arma::mat y = mdl["y"];
  int q       = y.n_cols;
  int Tsize   = y.n_rows;
  int ar      = mdl["p"];
  bool msmu   = mdl["msmu"];
  bool msvar  = mdl["msvar"];
  int M       = pow(k, ar+1);
  List pList  = paramList_MSVARmdl(theta, q, ar, k, msmu, msvar);
  arma::mat mu        = pList["mu"];
  List sig            = pList["sigma"];
  arma::mat P         = pList["P"];
  arma::vec pinf      = pList["pinf"];
  List muAR           = pList["muAR"];
  List sigAR          = pList["sigAR"];
  arma::mat phi       = pList["phi"];
  arma::mat PAR       = pList["P_AR"];
  arma::vec pinfAR    = pList["pinf_AR"];
  arma::vec state_ind = pList["state_ind"];
  // ----- Compute residuals in each regime
  List mdl_tmp    = clone(mdl);
  mdl_tmp["phi"]  = phi;
  List eps        = calcResid_MSVARmdl(mdl_tmp, muAR, k);
  // ----- Begin Calculating log-likelihood & filtered probabilities
  arma::mat eta(Tsize, M, arma::fill::zeros);       // [eq. 22.4.2]
  arma::mat f_t(Tsize, 1, arma::fill::zeros);       // [eq. 22.4.8]
  arma::mat xi_t_t(Tsize, k, arma::fill::zeros);    // [eq. 22.4.5]
  arma::mat xi_t_t_tmp(Tsize, M, arma::fill::zeros);    
  arma::mat xi_tp1_t(Tsize, k, arma::fill::zeros);  // [eq. 22.4.6]
  arma::mat xi_tp1_t_tmp(Tsize, M, arma::fill::zeros);  
  double pi = arma::datum::pi;
  // Pre-compute per-regime density constants (sigAR constant during T-loop)
  std::vector<arma::mat> inv_sigma_precomp(M);
  arma::vec norm_const_precomp(M);
  for (int xm_pre = 0; xm_pre < M; xm_pre++) {
    arma::mat sm = sigAR[xm_pre];
    double det_val = arma::det(sm);
    norm_const_precomp(xm_pre) = (det_val > 0) ? 1.0 / sqrt(pow(2.0*pi, q) * det_val) : 0.0;
    arma::mat inv_sm;
    bool ok = arma::inv(inv_sm, sm);
    inv_sigma_precomp[xm_pre] = ok ? inv_sm : arma::zeros<arma::mat>(q, q);
  }
  arma::vec xi_t_tm1_AR = pinfAR;
  arma::vec xi_t_tm1 = pinf;
  for (int xt = 0; xt<Tsize; xt++){
    for (int xm = 0; xm<M; xm++){
      arma::mat eps_m = eps[xm];
      arma::vec e_vec = trans(eps_m.row(xt));
      eta(xt,xm) = norm_const_precomp(xm) * exp(-0.5 * arma::dot(e_vec, inv_sigma_precomp[xm] * e_vec));
    }
    arma::vec xi_eta = xi_t_tm1_AR%trans(eta.row(xt));
    f_t.row(xt) = sum(xi_eta);
    double f_tmp = sum(xi_eta);
    if (f_tmp > 0 && std::isfinite(f_tmp)) {
      xi_t_t_tmp.row(xt) = trans(xi_eta/f_tmp);
    } else {
      // Underflow: all regime densities collapsed to zero.
      // Reset to ergodic distribution to break NaN cascade.
      // f_t stores true 0; log(0)=-Inf makes logL=-Inf (theta ruled out).
      xi_t_t_tmp.row(xt) = trans(pinfAR);
    }
    for (int xk = 1; xk<=k; xk++){
      arma::vec xi_t_t_row = trans(xi_t_t_tmp.row(xt));
      xi_t_t.submat(xt,xk-1,xt,xk-1) = sum(xi_t_t_row.rows(find(state_ind==xk)));
    }
    xi_t_tm1_AR = predict_fast(trans(xi_t_t_tmp.row(xt)), P, k, ar);
    xi_t_tm1 = P*trans(xi_t_t.row(xt));
    xi_tp1_t_tmp.row(xt) = trans(xi_t_tm1_AR);
    xi_tp1_t.row(xt) = trans(xi_t_tm1);
  }
  // ----- log-Likelihood
  arma::vec logf = log(f_t);
  double logL = sum(logf);     //[eq. 22.4.7]
  double L = logL/Tsize; 
  // ---------- Hamilton smoother
  arma::mat xi_t_T(Tsize, k, arma::fill::zeros);
  arma::mat xi_t_T_tmp(Tsize, M, arma::fill::zeros);
  xi_t_T_tmp.row(Tsize-1)  = xi_t_t_tmp.row(Tsize-1);
  for (int xT = Tsize-2; xT>=0; xT--){
    arma::rowvec xi_tp1_tmp_safe = xi_tp1_t_tmp.row(xT);
    xi_tp1_tmp_safe.clamp(1e-300, 1.0);  // guard against exact-zero denominator
    xi_t_T_tmp.row(xT) = xi_t_t_tmp.row(xT)%trans(smooth_fast(trans(xi_t_T_tmp.row(xT+1)/xi_tp1_tmp_safe), P, k, ar));
    double row_sum_M = accu(xi_t_T_tmp.row(xT));
    if (row_sum_M > 0 && std::isfinite(row_sum_M)) {
      xi_t_T_tmp.row(xT) = xi_t_T_tmp.row(xT) / row_sum_M;
    } else {
      xi_t_T_tmp.row(xT) = xi_t_t_tmp.row(xT);
    }
  }
  // Marginal (k-state) smoothed probabilities: the conditional density of an MS-VAR depends
  // on the date-t and lagged regimes, so smoothed inference lives on the extended state
  // (Hamilton 1994, sec. 22.4); the marginal is its aggregate over the date-t regime.
  for (int xk = 1; xk<=k; xk++){
    xi_t_T.col(xk-1) = arma::sum(xi_t_T_tmp.cols(find(state_ind==xk)), 1);
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
//' @title Markov-switching VARX log-likelihood function
//' 
//' @description This function computes the log-likelihood for a markov-switching VARX model and uses the Hamilton smoother to obtain smoothed probabilities of each state. This is also the expectation step in the Expectation Maximization algorithm for a Markov-switching autoregressive model.
//' 
//' @param theta Vector of model parameters.
//' @param mdl List with model attributes.
//' @param k Integer determining the number of regimes.
//'  
//' @return List which includes log-likelihood and smoothed probabilities of each regime.
//' 
//' @keywords internal
//' 
//' @export
// [[Rcpp::export]]
List ExpectationM_MSVARXmdl(arma::vec theta, List mdl, int k){
  // ---------- Initialize parameters
  arma::mat y   = mdl["y"];
  int q         = y.n_cols;
  int Tsize     = y.n_rows;
  int ar        = mdl["p"];
  bool msmu     = mdl["msmu"];
  bool msvar    = mdl["msvar"];
  int M         = pow(k, ar+1);
  List mdl_con  = mdl["control"];
  arma::mat Z   = mdl_con["Z"];
  int qz        = Z.n_cols;
  List pList  = paramList_MSVARXmdl(theta, q, ar, k, qz, msmu, msvar);
  arma::mat mu        = pList["mu"];
  arma::mat phi       = pList["phi"];
  arma::mat betaZ     = pList["betaZ"];
  List sig            = pList["sigma"];
  arma::mat P         = pList["P"];
  arma::vec pinf      = pList["pinf"];
  List muAR           = pList["muAR"];
  List sigAR          = pList["sigAR"];
  arma::mat PAR       = pList["P_AR"];
  arma::vec pinfAR    = pList["pinf_AR"];
  arma::vec state_ind = pList["state_ind"];
  // ----- Compute residuals in each regime
  List mdl_tmp        = clone(mdl);
  mdl_tmp["phi"]      = phi;
  mdl_tmp["betaZ"]    = betaZ;
  List eps            = calcResid_MSVARXmdl(mdl_tmp, muAR, k);
  // ----- Begin Calculating log-likelihood & filtered probabilities
  arma::mat eta(Tsize, M, arma::fill::zeros);       // [eq. 22.4.2]
  arma::mat f_t(Tsize, 1, arma::fill::zeros);       // [eq. 22.4.8]
  arma::mat xi_t_t(Tsize, k, arma::fill::zeros);    // [eq. 22.4.5]
  arma::mat xi_t_t_tmp(Tsize, M, arma::fill::zeros);    
  arma::mat xi_tp1_t(Tsize, k, arma::fill::zeros);  // [eq. 22.4.6]
  arma::mat xi_tp1_t_tmp(Tsize, M, arma::fill::zeros);  
  double pi = arma::datum::pi;
  // Pre-compute per-regime density constants (sigAR constant during T-loop)
  std::vector<arma::mat> inv_sigma_precomp(M);
  arma::vec norm_const_precomp(M);
  for (int xm_pre = 0; xm_pre < M; xm_pre++) {
    arma::mat sm = sigAR[xm_pre];
    double det_val = arma::det(sm);
    norm_const_precomp(xm_pre) = (det_val > 0) ? 1.0 / sqrt(pow(2.0*pi, q) * det_val) : 0.0;
    arma::mat inv_sm;
    bool ok = arma::inv(inv_sm, sm);
    inv_sigma_precomp[xm_pre] = ok ? inv_sm : arma::zeros<arma::mat>(q, q);
  }
  arma::vec xi_t_tm1_AR = pinfAR;
  arma::vec xi_t_tm1 = pinf;
  for (int xt = 0; xt<Tsize; xt++){
    for (int xm = 0; xm<M; xm++){
      arma::mat eps_m = eps[xm];
      arma::vec e_vec = trans(eps_m.row(xt));
      eta(xt,xm) = norm_const_precomp(xm) * exp(-0.5 * arma::dot(e_vec, inv_sigma_precomp[xm] * e_vec));
    }
    arma::vec xi_eta = xi_t_tm1_AR%trans(eta.row(xt));
    f_t.row(xt) = sum(xi_eta);
    double f_tmp = sum(xi_eta);
    if (f_tmp > 0 && std::isfinite(f_tmp)) {
      xi_t_t_tmp.row(xt) = trans(xi_eta/f_tmp);
    } else {
      // Underflow: all regime densities collapsed to zero.
      // Reset to ergodic distribution to break NaN cascade.
      // f_t stores true 0; log(0)=-Inf makes logL=-Inf (theta ruled out).
      xi_t_t_tmp.row(xt) = trans(pinfAR);
    }
    for (int xk = 1; xk<=k; xk++){
      arma::vec xi_t_t_row = trans(xi_t_t_tmp.row(xt));
      xi_t_t.submat(xt,xk-1,xt,xk-1) = sum(xi_t_t_row.rows(find(state_ind==xk)));
    }
    xi_t_tm1_AR = predict_fast(trans(xi_t_t_tmp.row(xt)), P, k, ar);
    xi_t_tm1 = P*trans(xi_t_t.row(xt));
    xi_tp1_t_tmp.row(xt) = trans(xi_t_tm1_AR);
    xi_tp1_t.row(xt) = trans(xi_t_tm1);
  }
  // ----- log-Likelihood
  arma::vec logf = log(f_t);
  double logL = sum(logf);     //[eq. 22.4.7]
  double L = logL/Tsize; 
  // ---------- Hamilton smoother
  arma::mat xi_t_T(Tsize, k, arma::fill::zeros);
  arma::mat xi_t_T_tmp(Tsize, M, arma::fill::zeros);
  xi_t_T_tmp.row(Tsize-1)  = xi_t_t_tmp.row(Tsize-1);
  for (int xT = Tsize-2; xT>=0; xT--){
    arma::rowvec xi_tp1_tmp_safe = xi_tp1_t_tmp.row(xT);
    xi_tp1_tmp_safe.clamp(1e-300, 1.0);  // guard against exact-zero denominator
    xi_t_T_tmp.row(xT) = xi_t_t_tmp.row(xT)%trans(smooth_fast(trans(xi_t_T_tmp.row(xT+1)/xi_tp1_tmp_safe), P, k, ar));
    double row_sum_M = accu(xi_t_T_tmp.row(xT));
    if (row_sum_M > 0 && std::isfinite(row_sum_M)) {
      xi_t_T_tmp.row(xT) = xi_t_T_tmp.row(xT) / row_sum_M;
    } else {
      xi_t_T_tmp.row(xT) = xi_t_t_tmp.row(xT);
    }
  }
  // Marginal (k-state) smoothed probabilities: the conditional density of an MS-VAR depends
  // on the date-t and lagged regimes, so smoothed inference lives on the extended state
  // (Hamilton 1994, sec. 22.4); the marginal is its aggregate over the date-t regime.
  for (int xk = 1; xk<=k; xk++){
    xi_t_T.col(xk-1) = arma::sum(xi_t_T_tmp.cols(find(state_ind==xk)), 1);
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
  MSloglik_output["logLike"]    = logL;
  MSloglik_output["L"]          = L;
  MSloglik_output["eta"]        = eta;
  MSloglik_output["f_t"]        = f_t;
  MSloglik_output["xi_t_t_AR"]  = xi_t_t_tmp;
  MSloglik_output["xi_tp1_t_AR"]= xi_tp1_t_tmp;
  MSloglik_output["xi_t_T_AR"]  = xi_t_T_tmp;
  MSloglik_output["xi_t_t"]     = xi_t_t;
  MSloglik_output["xi_tp1_t"]   = xi_tp1_t;
  MSloglik_output["xi_t_T"]     = xi_t_T;
  MSloglik_output["PAR"]        = PAR;
  MSloglik_output["P"]          = P;
  MSloglik_output["pinfAR"]     = pinfAR;
  MSloglik_output["pinf"]       = pinf;
  MSloglik_output["muAR"]       = muAR;
  MSloglik_output["sigmaAR"]    = sigAR;
  MSloglik_output["mu"]         = mu;
  MSloglik_output["sigma"]      = sig;
  MSloglik_output["phi"]        = phi;
  MSloglik_output["betaZ"]      = betaZ;
  MSloglik_output["theta"]      = theta;
  MSloglik_output["state_ind"]  = state_ind;
  MSloglik_output["residuals"]  = eps;
  MSloglik_output["resid"]      = residuals;
  return(MSloglik_output);
}

// ==============================================================================

// ============================================================================== #
// EM regime-variance floor (constrained ML). The floor sigma^2_k >= c*sigma^2_lin
// is set per dataset in the R constructors (control 'em_variance_constraint') and
// passed in the model list as 'sigma_floor' (0 or absent = disabled, which
// reproduces the unconstrained legacy behaviour for direct *_em callers).
// Rationale: the Gaussian mixture / Markov-switching likelihood is unbounded --
// a regime variance driven to zero on a few observations sends the likelihood to
// infinity (Kiefer and Wolfowitz 1956; Day 1969; Hamilton 1994, p. 689) -- so the
// unconstrained MLE does not exist. Constraining the regime variances restores a
// well-defined maximizer (Hathaway 1985, 1986, who bounds variance ratios; the
// data-dependent floor used here follows Kasahara and Shimotsu 2018).
// ============================================================================== #

// Read the per-dataset variance floor from the model list (0.0 = disabled).
double get_sigma_floor(List mdl){
  if (mdl.containsElementNamed("sigma_floor")){
    Rcpp::RObject sf = mdl["sigma_floor"];
    if (!sf.isNULL()){
      double val = Rcpp::as<double>(sf);
      if (std::isfinite(val) && val > 0.0) return val;
    }
  }
  return 0.0;
}

// Floor a covariance matrix at sig_floor*I by eigenvalue clipping. This is the
// EXACT maximizer of the M-step objective over {Sigma : Sigma >= sig_floor*I}
// (the objective is orthogonally invariant and separates over the eigenvalues),
// not a ridge-type approximation.
arma::mat floor_cov(arma::mat S, double sig_floor){
  if (sig_floor <= 0.0) return S;
  arma::vec eigval;
  arma::mat eigvec;
  arma::mat S_sym = arma::symmatu(S);
  if (!arma::eig_sym(eigval, eigvec, S_sym)) return S;   // keep as-is on failure
  if (eigval.min() >= sig_floor) return S;               // constraint slack
  eigval = arma::clamp(eigval, sig_floor, arma::datum::inf);
  return arma::symmatu(eigvec * arma::diagmat(eigval) * trans(eigvec));
}

// Apply the floor to the multivariate M-step output, keeping the per-regime List
// 'sigma' and the vech-stacked vector 'sigma_out' consistent.
void apply_sigma_floor_mv(List& sigma, arma::vec& sigma_out, double sig_floor,
                          int q, int k, bool msvar){
  if (sig_floor <= 0.0) return;
  int sigN = (q*(q+1))/2;
  int nblk = (msvar ? k : 1);
  for (int xb = 0; xb < nblk; xb++){
    arma::mat S_f = floor_cov(as<arma::mat>(sigma[xb]), sig_floor);
    sigma_out.subvec(xb*sigN, (xb+1)*sigN - 1) = covar_vech(S_f);
    if (msvar){
      sigma[xb] = S_f;
    }else{
      for (int xk = 0; xk < k; xk++) sigma[xk] = S_f;
    }
  }
}

// Project the sigma block of a STARTING value onto the variance floor, so the EM
// iteration begins from a feasible point (monotone ascent of the constrained
// likelihood then holds from the first iteration). Layout-agnostic: in every
// model class the sigma block sits immediately before the k*k transition
// probabilities at the end of theta.
arma::vec project_theta_sigma_floor(arma::vec theta_0, List mdl, int k){
  double sig_floor = get_sigma_floor(mdl);
  if (sig_floor <= 0.0) return theta_0;
  int q = mdl.containsElementNamed("q") ? Rcpp::as<int>(mdl["q"]) : 1;
  bool msvar = mdl.containsElementNamed("msvar") ? Rcpp::as<bool>(mdl["msvar"]) : true;
  int sigN = (q*(q+1))/2;
  int nblk = 1 + (msvar ? (k-1) : 0);
  int len  = sigN * nblk;
  int n    = (int) theta_0.n_elem;
  int end  = n - k*k - 1;          // last sigma entry (P block is last)
  int start = end - len + 1;
  if (start < 0 || end >= n || end < start) return theta_0;   // unexpected layout: leave as-is
  if (q == 1){
    theta_0.subvec(start, end) = arma::clamp(theta_0.subvec(start, end), sig_floor, arma::datum::inf);
  }else{
    for (int xb = 0; xb < nblk; xb++){
      int s0 = start + xb*sigN;
      arma::mat S_f = floor_cov(covar_unvech(theta_0.subvec(s0, s0 + sigN - 1), q), sig_floor);
      theta_0.subvec(s0, s0 + sigN - 1) = covar_vech(S_f);
    }
  }
  return theta_0;
}

// ------------------------------------------------------------------------------
// Transition-probability bound (constrained EM for P; control
// em_transition_constraint). The bound eps is passed in the model list as
// 'P_bound' (0 or absent = disabled).
double get_P_bound(List mdl){
  if (mdl.containsElementNamed("P_bound")){
    Rcpp::RObject pb = mdl["P_bound"];
    if (!pb.isNULL()){
      double eps = Rcpp::as<double>(pb);
      if (std::isfinite(eps) && (eps > 0.0)) return eps;
    }
  }
  return 0.0;
}

// Exact constrained M-step for one transition-matrix column: maximizes
// sum_j n_j log p_j subject to sum_j p_j = 1 and p_j >= eps (the upper bound
// 1-eps is implied). KKT water-filling solution: entries with small counts pin
// at eps, free entries stay proportional to their counts; the free set only
// shrinks, so at most k passes are needed. Feasibility requires eps < 1/k.
// Returns fallback when the counts are unusable (non-finite, non-positive sum)
// or eps is out of range.
// [[Rcpp::export]]
arma::vec constrain_P_col(arma::vec counts, double eps, arma::vec fallback){
  int k = counts.n_elem;
  double s = arma::accu(counts);
  if ((!counts.is_finite()) || (!std::isfinite(s)) || (!(s > 0.0)) || (!(eps > 0.0)) ||
      (!(eps < 1.0/k)) || (arma::any(counts < 0.0))){
    return fallback;
  }
  arma::vec q = counts / s;
  if (arma::all(q >= eps)) return q;   // unconstrained optimum feasible => it is the solution
  std::vector<bool> pinned(k, false);
  int nb = 0;
  for (int pass = 0; pass < k; pass++){
    double nf = 0.0;
    for (int j = 0; j < k; j++){ if (!pinned[j]) nf += counts(j); }
    double mass = 1.0 - nb * eps;
    bool changed = false;
    for (int j = 0; j < k; j++){
      // pin iff n_j / lambda < eps with lambda = nf / mass
      if ((!pinned[j]) && (counts(j) * mass < eps * nf)){ pinned[j] = true; nb++; changed = true; }
    }
    if (!changed){
      arma::vec pcol(k);
      for (int j = 0; j < k; j++){ pcol(j) = pinned[j] ? eps : counts(j) * mass / nf; }
      return pcol;
    }
  }
  return fallback;  // unreachable for eps < 1/k (the free set cannot empty)
}

// Project the transition-matrix block of a parameter vector onto the P bound
// (feasible starting values; the last k*k entries of theta are vec(P),
// column-major). Non-finite columns fall back to the uniform column.
arma::vec project_theta_P_bound(arma::vec theta_0, List mdl, int k){
  double eps = get_P_bound(mdl);
  if (!(eps > 0.0)) return theta_0;
  int npar = theta_0.n_elem;
  if (npar < k*k) return theta_0;   // unexpected layout: leave as-is
  arma::vec unif(k);
  unif.fill(1.0/k);
  for (int xk = 0; xk < k; xk++){
    int i0 = npar - k*k + xk*k;
    arma::vec col = theta_0.subvec(i0, i0 + k - 1);
    // feasible columns stay bitwise untouched (renormalizing a feasible start
    // would perturb it by ~1 ulp and change the whole EM trajectory); feasible
    // means a valid probability column within the bound, so malformed user
    // starts (entries > 1, column sum away from 1) are repaired here
    if (col.is_finite() && arma::all(col >= eps) && arma::all(col <= 1.0) &&
        (std::abs(arma::accu(col) - 1.0) <= 1e-8)) continue;
    theta_0.subvec(i0, i0 + k - 1) = constrain_P_col(col, eps, unif);
  }
  return theta_0;
}

//' @title Maximization step of EM algorithm for Hidden Markov model
//' 
//' @description This function performs the maximization step of the Expectation Maximization algorithm for Hidden Markov models.
//' 
//' @param theta Vector of model parameters.
//' @param mdl List with model attributes.
//' @param MSloglik_output List with output from \code{\link{ExpectationM_HMmdl}}.
//' @param k Integer determining the number of regimes.
//'  
//' @return List with new maximized parameters.
//' 
//' @keywords internal
//' 
//' @export
// [[Rcpp::export]]
List EMaximization_HMmdl(arma::vec theta, List mdl, List MSloglik_output, int k){
  // ---------- Initialize parameters
  arma::mat y = mdl["y"];
  int q       = mdl["q"];
  bool msmu   = mdl["msmu"];
  bool msvar  = mdl["msvar"];
  bool exog   = mdl["exog"];
  // ----- Output from log-Likelihood calculation
  arma::mat P         = MSloglik_output["P"];
  arma::mat xi_t_T    = MSloglik_output["xi_t_T"];
  arma::mat xi_t_t    = MSloglik_output["xi_t_t"];
  arma::mat xi_tp1_t  = MSloglik_output["xi_tp1_t"];
  List sig            = MSloglik_output["sigma"];
  // Length of Time-Series (T) and number of variables
  int Tsize   = y.n_rows;
  arma::mat repq(1,q, arma::fill::ones);
  // --------------- Update transition matrix (P) and limiting probabilities (pinf)
  // ----- Estimate Transition probs for k-state Markov-chain [eq. 22.4.16]
  arma::mat xi_t_T_tmp    = xi_t_T.rows(1,Tsize-1);
  arma::mat xi_t_t_tmp    = xi_t_t.rows(0,Tsize-2);
  arma::mat xi_tp1_t_safe = arma::clamp(xi_tp1_t.rows(0,Tsize-2), 1e-300, 1.0);
  arma::mat p_ij(Tsize-1, k*k, arma::fill::zeros);
  for (int xk = 0; xk<k; xk++){
    p_ij.submat(0,xk*k,Tsize-2,xk*k+k-1) = (xi_t_T_tmp%trans(P.col(xk)*trans(xi_t_t_tmp.col(xk))))/xi_tp1_t_safe;
  }
  arma::mat P_new = trans(P);  // Initialize from previous P; degenerate regimes keep their column
  arma::mat p_ij_sums     = arma::sum(p_ij,0);
  double P_bound = get_P_bound(mdl);
  for (int xk = 0 ; xk<k; xk++){
    arma::vec regime_prob = xi_t_T.submat(0,xk,(Tsize-2),xk);
    double regimesum      = sum(regime_prob);
    if (regimesum < 1e-6 * Tsize) continue;  // degenerate: keep previous P column
    arma::rowvec Pcol_uc  = p_ij_sums.submat(0,(xk*k),0,(xk*k+(k-1)))/regimesum;
    if ((P_bound > 0.0) && Pcol_uc.is_finite() && arma::any(vectorise(Pcol_uc) < P_bound)){
      // exact constrained M-step for this column (see constrain_P_col)
      P_new.row(xk) = trans(constrain_P_col(trans(Pcol_uc), P_bound, trans(P_new.row(xk))));
    } else {
      P_new.row(xk) = Pcol_uc;
    }
  }
  P_new = trans(P_new);
  // Guard against NaN in P_new (can occur when p_ij numerator is also 0)
  for (int xk = 0; xk < k; xk++) {
    if (!P_new.col(xk).is_finite()) P_new.col(xk) = P.col(xk);  // keep previous
  }
  // Enforce valid transition probabilities: clamp to [0,1] and normalize columns to sum to 1.
  // In exact arithmetic, P_new is already in [0,1] with unit column sums (Krolzig 1997, eq. 6.14).
  // This guard defends against floating-point drift only.
  P_new = arma::clamp(P_new, 0.0, 1.0);
  for (int xk = 0; xk < k; xk++) {
    double colsum = arma::accu(P_new.col(xk));
    if (colsum > 0) {
      P_new.col(xk) /= colsum;
    } else {
      P_new.col(xk) = P.col(xk);  // fallback: keep previous
    }
  }
  // Get limiting probabilities implied by transition matrix P_0
  arma::vec pinf = limP(P_new);
  // --------------- Update estimates for mu
  // Initialize from previous theta (used as fallback for degenerate regimes)
  int mu_len = q * (1 + msmu*(k-1));
  arma::mat mu = arma::reshape(theta.head(mu_len), q, 1 + msmu*(k-1)).t();
  arma::mat mu_k(k, q, arma::fill::zeros);
  // ----- Compute updated mu
  if (msmu == TRUE){
    for (int xk = 0 ; xk<k; xk++){
      double mu_denom = sum(xi_t_T.col(xk));
      if (mu_denom < 1e-6 * Tsize) continue;  // degenerate: keep previous mu
      mu.row(xk) = arma::sum(y%(xi_t_T.col(xk)*repq),0)/mu_denom;
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
  // --------------- Update estimates params of exogenous regressors (if present)
  int qz = 0;
  arma::mat repmu(Tsize, 1, arma::fill::ones);
  arma::mat Zdm(Tsize,1,arma::fill::zeros);
  arma::mat betaZ_new(1,q,arma::fill::zeros);
  if (exog==TRUE){
    arma::mat Z = as<arma::mat>(mdl["Z"]);
    qz = Z.n_cols;
    // Initialize betaZ_new from previous theta (fallback for solve failure)
    betaZ_new = reshape(theta.subvec(mu_len, mu_len + qz*q - 1), qz, q);
    arma::rowvec zbar  = arma::mean(Z,0);
    Zdm = Z - (repmu*zbar);
    arma::mat denom(q*qz, q*qz, arma::fill::zeros);
    arma::mat num(q*qz, 1, arma::fill::zeros);
    for (int xk = 0; xk<k; xk++){
      arma::mat mu_tmp = mu_k.row(xk);
      arma::mat y_tmp = y - repmu*mu_tmp;
      arma::mat sigma_m = sig[xk];
      arma::mat sigma_m_inv = arma::solve(sigma_m, arma::eye(q,q), arma::solve_opts::allow_ugly);
      denom = denom + kron(sigma_m_inv,trans(Zdm)*diagmat(xi_t_T.col(xk))*Zdm);
      num = num + kron(sigma_m_inv,trans(Zdm)*diagmat(xi_t_T.col(xk)))*vectorise(y_tmp);
    }
    try {
      arma::vec bz_vec = arma::solve(denom, num, arma::solve_opts::allow_ugly);
      if (bz_vec.is_finite()) {
        betaZ_new = reshape(bz_vec, qz, q);
      }
      // else: keep betaZ_new from previous theta
    } catch (...) {
      // keep betaZ_new from previous theta
    }
  }
  // --------------- Update estimates for variance
  // ----- Compute Residuals
  List eps(k);
  for (int xk = 0; xk<k; xk++){
    arma::mat eps_tmp = y - repmu*mu_k.row(xk) - Zdm*betaZ_new;
    eps[xk] =  eps_tmp;
  } 
  // // ----- Compute updated sigma
  List sigma(k);
  int sigN = (q * (q + 1)) / 2;
  arma::vec sigma_out;
  if (msvar == TRUE){
    sigma_out.set_size(sigN * k);
    for (int xk = 0 ; xk<k; xk++){
      arma::mat U = eps[xk];
      double sig_denom = sum(xi_t_T.col(xk));
      arma::mat sigma_m_tmp;
      if (sig_denom >= 1e-6 * Tsize) {
        sigma_m_tmp = (trans(U)*diagmat(xi_t_T.col(xk))*U)/sig_denom;
      } else {
        sigma_m_tmp = as<arma::mat>(sig[xk]);  // degenerate: keep previous sigma
      }
      sigma[xk] = sigma_m_tmp;
      sigma_out.subvec(xk * sigN, (xk+1)*sigN - 1) = covar_vech(sigma_m_tmp);
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
  // ----- Apply the EM regime-variance floor (constrained ML; see get_sigma_floor)
  apply_sigma_floor_mv(sigma, sigma_out, get_sigma_floor(mdl), q, k, msvar);
  // ----- Compute model residuals 
  arma::mat residuals(Tsize, q, arma::fill::zeros); 
  for (int xk = 0; xk<k;xk++){
    arma::mat eps_tmp = eps[xk];
    residuals = residuals + eps_tmp%(xi_t_T.col(xk)*repq);
  }
  // --------------- Organize output
  // ----- Produce new theta vector
  arma::vec theta_new = vectorise(trans(mu));
  if (exog==TRUE){
    theta_new = join_vert(theta_new, vectorise(betaZ_new));  
  }
  theta_new = join_vert(theta_new, sigma_out);
  theta_new = join_vert(theta_new, vectorise(P_new));
  // ----- Fill output List
  List Maxim_output;
  Maxim_output["theta"] = theta_new;
  Maxim_output["P"] = P_new;
  Maxim_output["pinf"] = pinf;
  Maxim_output["mu"] = mu;
  Maxim_output["betaZ"] = betaZ_new;
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
//' @param theta Vector of model parameters.
//' @param mdl List with model attributes.
//' @param MSloglik_output List with output from \code{\link{ExpectationM_MSARmdl}}.
//' @param k Integer determining the number of regimes.
//' 
//' @return List with new maximized parameters.
//' 
//' @keywords internal
//' 
//' @export
// [[Rcpp::export]]
List EMaximization_MSARmdl(arma::vec theta, List mdl, List MSloglik_output, int k){
  // ---------- Initialize parameters
  arma::vec y = mdl["y"];
  arma::mat x = mdl["x"];
  int ar      = mdl["p"];
  bool msmu   = mdl["msmu"];
  bool msvar  = mdl["msvar"];
  // ----- Output from log-Likelihood calculation
  arma::mat P           = MSloglik_output["P"];
  arma::mat xi_t_T      = MSloglik_output["xi_t_T"];
  arma::mat xi_t_t      = MSloglik_output["xi_t_t"];
  arma::mat xi_tp1_t    = MSloglik_output["xi_tp1_t"];
  arma::mat xi_t_T_AR   = MSloglik_output["xi_t_T_AR"];
  arma::mat xi_t_t_AR   = MSloglik_output["xi_t_t_AR"];
  arma::mat xi_tp1_t_AR = MSloglik_output["xi_tp1_t_AR"];
  arma::vec state_ind   = MSloglik_output["state_ind"];
  arma::vec phi         = MSloglik_output["phi"];
  arma::vec sig         = MSloglik_output["sigma"];
  int Tsize = y.n_elem; // Length of Time-Series (T)
  int M     = pow(k, ar+1); // Number of regimes consistent with AR lags
  // --------------- Update transition matrix (P) and limiting probabilities (pinf)
  // ----- Smoothed transition counts [Hamilton 1994, eq. 22.4.16; Krolzig 1997, eq. 6.14].
  // The extended state encodes (s_t, s_{t-1}, ...), so Pr(s_t = j, s_{t-1} = i | Y) is the
  // sum of xi_t_T_AR over extended states with date-t regime j and lag-1 regime i.
  arma::mat counts(k, k, arma::fill::zeros);
  for (int xn = 0; xn < M; xn++){
    int reg_i0 = xn % k;        // regime at date t
    int reg_i1 = (xn / k) % k;  // regime at date t-1
    counts(reg_i0, reg_i1) += arma::accu(xi_t_T_AR.submat(1, xn, Tsize-1, xn));
  }
  arma::mat P_new = P;  // degenerate regimes keep their previous column
  double P_bound = get_P_bound(mdl);
  for (int xk = 0; xk < k; xk++){
    double countsum = arma::accu(counts.col(xk));
    if (std::isfinite(countsum) && countsum >= 1e-6 * Tsize && counts.col(xk).is_finite()){
      arma::vec Pcol_uc = counts.col(xk) / countsum;
      if ((P_bound > 0.0) && arma::any(Pcol_uc < P_bound)){
        // exact constrained M-step for this column (see constrain_P_col)
        P_new.col(xk) = constrain_P_col(counts.col(xk), P_bound, P.col(xk));
      } else {
        P_new.col(xk) = Pcol_uc;
      }
    }
  }
  // Floating-point guard: keep P_new a valid column-stochastic matrix.
  P_new = arma::clamp(P_new, 0.0, 1.0);
  for (int xk = 0; xk < k; xk++) {
    double colsum = arma::accu(P_new.col(xk));
    if (colsum > 0) {
      P_new.col(xk) /= colsum;
    } else {
      P_new.col(xk) = P.col(xk);  // fallback: keep previous
    }
  }
  // Get limiting probabilities implied by transition matrix P_0
  arma::vec pinf = limP(P_new);
  // --------------- Update estimates for mu
  // ----- Obtain state indicators
  arma::vec mu = theta.head(1 + msmu*(k-1));  // Initialize from previous theta (fallback)
  arma::vec sum_tmp(M, arma::fill::zeros);
  // ----- Compute updated mu (Krolzig 1997, Table 9.19 — MSM GLS estimator)
  {
    // ỹ_t = y_t - sum_j phi(j)*x(t,j)  [regime-independent, T×1]
    arma::vec ytilde = y;
    for (int j = 0; j < ar; j++) ytilde -= phi(j) * x.col(j);
    // Expand sigma to length k (handles msvar=T and msvar=F)
    arma::vec sig_k(k);
    for (int xk = 0; xk < k; xk++) sig_k(xk) = (sig.n_elem > 1) ? sig(xk) : sig(0);
    if (msmu == TRUE){
      // Build c_{n,k} [M×k]: c = +I[i0=k] - sum_j phi(j)*I[ij=k]  (Krolzig 1997, Table 9.19)
      arma::mat d_mat(M, k, arma::fill::zeros);
      for (int xn = 0; xn < M; xn++){
        int kpow = 1;
        for (int j = 0; j <= ar; j++){
          int reg_ij = (xn / kpow) % k;
          d_mat(xn, reg_ij) += (j == 0) ? +1.0 : -phi(j-1);
          kpow *= k;
        }
      }
      // Build A (k×k) and b (k×1)
      arma::mat A(k, k, arma::fill::zeros);
      arma::vec b_vec(k, arma::fill::zeros);
      for (int xn = 0; xn < M; xn++){
        double sig_xn = sig_k(xn % k);
        if (!std::isfinite(sig_xn) || sig_xn < 1e-14) continue;
        double wt_sum    = arma::sum(xi_t_T_AR.col(xn)) / sig_xn;
        double ytilde_wt = arma::dot(xi_t_T_AR.col(xn) / sig_xn, ytilde);
        for (int xk = 0; xk < k; xk++){
          b_vec(xk) += d_mat(xn, xk) * ytilde_wt;
          for (int xl = 0; xl < k; xl++)
            A(xk, xl) += d_mat(xn, xk) * d_mat(xn, xl) * wt_sum;
        }
      }
      // Solve A·mu = b; fall back to previous mu if A is near-singular
      if (A.is_finite() && arma::rcond(A) > 1e-12){
        try {
          arma::vec mu_new = arma::solve(A, b_vec, arma::solve_opts::allow_ugly);
          if (mu_new.is_finite()) mu = mu_new;
        } catch (...) { /* keep previous mu */ }
      }
    }else{
      // Single common mean: mu = (sum_tn xi*ytilde/sig) / (c * sum_tn xi/sig)
      // where c = 1 - sum_j phi(j)
      double c_val = 1.0;
      for (int j = 0; j < ar; j++) c_val -= phi(j);
      double numer = 0.0, denom_val = 0.0;
      for (int xn = 0; xn < M; xn++){
        double sig_xn = sig_k(xn % k);
        if (!std::isfinite(sig_xn) || sig_xn < 1e-14) continue;
        numer     += arma::dot(xi_t_T_AR.col(xn) / sig_xn, ytilde);
        denom_val += arma::sum(xi_t_T_AR.col(xn)) / sig_xn;
      }
      if (std::isfinite(numer) && std::abs(c_val * denom_val) > 1e-14)
        mu(0) = numer / (c_val * denom_val);
      // else: keep previous mu
    }
  }
  // ----- Update estimates for phi
  arma::mat muAR(M, ar+1,arma::fill::zeros);
  arma::vec sigAR(M, 1,arma::fill::zeros);
  List mugrid = argrid_MSARmdl_cpp(mu, sig, k, ar); // new mu & old sigma are used
  muAR        = as<arma::mat>(mugrid["mu"]);
  sigAR       = as<arma::vec>(mugrid["sig"]);
  arma::vec repmu(Tsize,arma::fill::ones);
  arma::mat denom(ar,ar,arma::fill::zeros);
  arma::mat num(ar,1,arma::fill::zeros);
  for (int xm = 0; xm<M; xm++){
    double sig_xm = sigAR(xm);
    if (!std::isfinite(sig_xm) || sig_xm < 1e-14) continue;  // degenerate sigma: skip
    arma::vec y_tmp = y - repmu*muAR(xm,0);
    arma::mat x_tmp = x - repmu*muAR.submat(xm,1,xm,ar);
    denom = denom + (trans(x_tmp)*diagmat(xi_t_T_AR.col(xm)/sig_xm)*x_tmp);
    num = num + (trans(x_tmp)*diagmat(xi_t_T_AR.col(xm)/sig_xm)*y_tmp);
  }
  arma::vec phi_new;
  try {
    phi_new = arma::solve(denom, num, arma::solve_opts::allow_ugly);
    if (!phi_new.is_finite()) phi_new = phi;
  } catch (...) {
    phi_new = phi;  // Singular denom: keep previous phi
  }
  // Stationarity check via companion matrix eigenvalues
  if (ar > 0) {
    arma::mat F_comp(ar, ar, arma::fill::zeros);
    F_comp.row(0) = phi_new.t();
    if (ar > 1) F_comp.submat(1, 0, ar-1, ar-2) = arma::eye(ar-1, ar-1);
    arma::cx_vec eigvals = arma::eig_gen(F_comp);
    if (arma::max(arma::abs(eigvals)) >= 1.0) {
      phi_new = phi;  // revert to previous (stationary) phi
    }
  }
  // ----- Update estimates for variance
  // ----- Compute new residuals
  List mdl_tmp    = clone(mdl);
  mdl_tmp["phi"]  = phi_new;
  arma::mat eps   = calcResid_MSARmdl(mdl_tmp, muAR, k);
  // ----- Compute updated sigma
  arma::vec sigma = sig.head(1 + msvar*(k-1));  // Initialize from previous sigma
  arma::vec sigma_tmp(M, arma::fill::zeros);
  if (msvar == TRUE){
    for (int xk = 0 ; xk<M; xk++){
      sigma_tmp(xk) = sum(eps.col(xk)%eps.col(xk)%xi_t_T_AR.col(xk));
      sum_tmp(xk) = sum(xi_t_T_AR.col(xk));
    }
    for (int xk = 1; xk<=k;xk++){
      double sig_denom = arma::as_scalar(sum(sum_tmp.rows(find(state_ind==xk))));
      if (sig_denom < 1e-6 * Tsize) continue;  // degenerate: keep previous sigma
      sigma(xk-1) = arma::as_scalar(sum(sigma_tmp.rows(find(state_ind==xk)))/sig_denom);
    }
  }else{
    sigma = arma::sum(arma::sum(eps%eps%xi_t_T_AR,0),1)/Tsize;
  }
  // ----- Compute model residuals 
  arma::mat residuals(Tsize, k);
  arma::mat residuals_tmp = eps%xi_t_T_AR; 
  for (int xk = 1; xk<=k; xk++){
    residuals.col(xk-1) = arma::sum(residuals_tmp.cols(find(state_ind==xk)), 1);
  }
  residuals = arma::sum(residuals,1);
  // --------------- Organize output
  // ----- Apply the EM regime-variance floor (constrained ML; see get_sigma_floor).
  // Univariate: clipping the variance update at the floor is the exact
  // constrained M-step maximizer (the M-step objective is unimodal in sigma^2).
  {
    double sig_floor_uv = get_sigma_floor(mdl);
    if (sig_floor_uv > 0.0) sigma = arma::clamp(sigma, sig_floor_uv, arma::datum::inf);
  }
  // ----- Produce new theta vector
  arma::vec theta_new = join_vert(mu, phi_new);
  theta_new = join_vert(theta_new, sigma);
  theta_new = join_vert(theta_new, vectorise(P_new));
  // ----- Fill output List
  List Maxim_output;
  Maxim_output["theta"]     = theta_new;
  Maxim_output["mu"]        = mu;
  Maxim_output["phi"]       = phi_new;
  Maxim_output["beta"]      = phi_new;
  Maxim_output["sigma"]     = sigma;
  Maxim_output["P"]         = P_new;
  Maxim_output["pinf"]      = pinf;
  Maxim_output["residuals"] = eps;
  Maxim_output["resid"]     = residuals;
  return(Maxim_output);
}

// ==============================================================================
//' @title Maximization step of EM algorithm for Markov-switching ARX model
//' 
//' @description This function performs the maximization step of the Expectation Maximization algorithm for Markov-switching ARX model.
//' 
//' @param theta Vector of model parameters.
//' @param mdl List with model attributes.
//' @param MSloglik_output List with output from \code{\link{ExpectationM_MSARmdl}}.
//' @param k Integer determining the number of regimes.
//' 
//' @return List with new maximized parameters.
//' 
//' @keywords internal
//' 
//' @export
// [[Rcpp::export]]
List EMaximization_MSARXmdl(arma::vec theta, List mdl, List MSloglik_output, int k){
  // ---------- Initialize parameters
  arma::vec y   = mdl["y"];
  arma::mat x   = mdl["x"];
  List mdl_con  = mdl["control"];
  arma::mat Z   = mdl_con["Z"];
  int ar        = mdl["p"];
  bool msmu     = mdl["msmu"];
  bool msvar    = mdl["msvar"];
  // ----- Output from log-Likelihood calculation
  arma::mat P           = MSloglik_output["P"];
  arma::mat xi_t_T      = MSloglik_output["xi_t_T"];
  arma::mat xi_t_t      = MSloglik_output["xi_t_t"];
  arma::mat xi_tp1_t    = MSloglik_output["xi_tp1_t"];
  arma::mat xi_t_T_AR   = MSloglik_output["xi_t_T_AR"];
  arma::mat xi_t_t_AR   = MSloglik_output["xi_t_t_AR"];
  arma::mat xi_tp1_t_AR = MSloglik_output["xi_tp1_t_AR"];
  arma::vec state_ind   = MSloglik_output["state_ind"];
  arma::vec phi         = MSloglik_output["phi"];
  arma::vec betaZ       = MSloglik_output["betaZ"];
  arma::vec sig         = MSloglik_output["sigma"];
  int Tsize             = y.n_elem; // Length of Time-Series (T)
  int M                 = pow(k, ar+1); // Number of regimes consistent with AR lags
  int qz                = Z.n_cols;
  Z                     = Z.rows(ar,Tsize+ar-1);
  arma::rowvec zbar     = arma::mean(Z,0);
  // --------------- Update transition matrix (P) and limiting probabilities (pinf)
  // ----- Smoothed transition counts [Hamilton 1994, eq. 22.4.16; Krolzig 1997, eq. 6.14].
  // The extended state encodes (s_t, s_{t-1}, ...), so Pr(s_t = j, s_{t-1} = i | Y) is the
  // sum of xi_t_T_AR over extended states with date-t regime j and lag-1 regime i.
  arma::mat counts(k, k, arma::fill::zeros);
  for (int xn = 0; xn < M; xn++){
    int reg_i0 = xn % k;        // regime at date t
    int reg_i1 = (xn / k) % k;  // regime at date t-1
    counts(reg_i0, reg_i1) += arma::accu(xi_t_T_AR.submat(1, xn, Tsize-1, xn));
  }
  arma::mat P_new = P;  // degenerate regimes keep their previous column
  double P_bound = get_P_bound(mdl);
  for (int xk = 0; xk < k; xk++){
    double countsum = arma::accu(counts.col(xk));
    if (std::isfinite(countsum) && countsum >= 1e-6 * Tsize && counts.col(xk).is_finite()){
      arma::vec Pcol_uc = counts.col(xk) / countsum;
      if ((P_bound > 0.0) && arma::any(Pcol_uc < P_bound)){
        // exact constrained M-step for this column (see constrain_P_col)
        P_new.col(xk) = constrain_P_col(counts.col(xk), P_bound, P.col(xk));
      } else {
        P_new.col(xk) = Pcol_uc;
      }
    }
  }
  // Floating-point guard: keep P_new a valid column-stochastic matrix.
  P_new = arma::clamp(P_new, 0.0, 1.0);
  for (int xk = 0; xk < k; xk++) {
    double colsum = arma::accu(P_new.col(xk));
    if (colsum > 0) {
      P_new.col(xk) /= colsum;
    } else {
      P_new.col(xk) = P.col(xk);  // fallback: keep previous
    }
  }
  // Get limiting probabilities implied by transition matrix P_0
  arma::vec pinf = limP(P_new);
  // --------------- Update estimates for mu
  // ----- Obtain state indicators
  arma::vec mu = theta.head(1 + msmu*(k-1));  // Initialize from previous theta (fallback)
  arma::vec sum_tmp(M, arma::fill::zeros);
  // Pre-compute zdm here so mu update can adjust ỹ_t for exogenous contribution
  arma::vec repmu(Tsize, arma::fill::ones);
  arma::mat zdm = (Z - repmu*zbar);
  // ----- Compute updated mu (Krolzig 1997, Table 9.19 — MSM GLS estimator)
  {
    // ỹ_t = y_t - betaZ'*zdm_t - sum_j phi(j)*x(t,j)  [regime-independent, T×1]
    arma::vec ytilde = y - zdm * betaZ;
    for (int j = 0; j < ar; j++) ytilde -= phi(j) * x.col(j);
    // Expand sigma to length k (handles msvar=T and msvar=F)
    arma::vec sig_k(k);
    for (int xk = 0; xk < k; xk++) sig_k(xk) = (sig.n_elem > 1) ? sig(xk) : sig(0);
    if (msmu == TRUE){
      // Build c_{n,k} [M×k]: c = +I[i0=k] - sum_j phi(j)*I[ij=k]  (Krolzig 1997, Table 9.19)
      arma::mat d_mat(M, k, arma::fill::zeros);
      for (int xn = 0; xn < M; xn++){
        int kpow = 1;
        for (int j = 0; j <= ar; j++){
          int reg_ij = (xn / kpow) % k;
          d_mat(xn, reg_ij) += (j == 0) ? +1.0 : -phi(j-1);
          kpow *= k;
        }
      }
      // Build A (k×k) and b (k×1)
      arma::mat A(k, k, arma::fill::zeros);
      arma::vec b_vec(k, arma::fill::zeros);
      for (int xn = 0; xn < M; xn++){
        double sig_xn = sig_k(xn % k);
        if (!std::isfinite(sig_xn) || sig_xn < 1e-14) continue;
        double wt_sum    = arma::sum(xi_t_T_AR.col(xn)) / sig_xn;
        double ytilde_wt = arma::dot(xi_t_T_AR.col(xn) / sig_xn, ytilde);
        for (int xk = 0; xk < k; xk++){
          b_vec(xk) += d_mat(xn, xk) * ytilde_wt;
          for (int xl = 0; xl < k; xl++)
            A(xk, xl) += d_mat(xn, xk) * d_mat(xn, xl) * wt_sum;
        }
      }
      // Solve A·mu = b; fall back to previous mu if A is near-singular
      if (A.is_finite() && arma::rcond(A) > 1e-12){
        try {
          arma::vec mu_new = arma::solve(A, b_vec, arma::solve_opts::allow_ugly);
          if (mu_new.is_finite()) mu = mu_new;
        } catch (...) { /* keep previous mu */ }
      }
    }else{
      // Single common mean: mu = (sum_tn xi*ytilde/sig) / (c * sum_tn xi/sig)
      double c_val = 1.0;
      for (int j = 0; j < ar; j++) c_val -= phi(j);
      double numer = 0.0, denom_val = 0.0;
      for (int xn = 0; xn < M; xn++){
        double sig_xn = sig_k(xn % k);
        if (!std::isfinite(sig_xn) || sig_xn < 1e-14) continue;
        numer     += arma::dot(xi_t_T_AR.col(xn) / sig_xn, ytilde);
        denom_val += arma::sum(xi_t_T_AR.col(xn)) / sig_xn;
      }
      if (std::isfinite(numer) && std::abs(c_val * denom_val) > 1e-14)
        mu(0) = numer / (c_val * denom_val);
    }
  }
  // ----- Update estimates for beta (phi & betaZ)
  arma::mat muAR(M, ar+1,arma::fill::zeros);
  arma::vec sigAR(M, 1,arma::fill::zeros);
  List mugrid   = argrid_MSARmdl_cpp(mu, sig, k, ar); // new mu & old sigma are used
  muAR          = as<arma::mat>(mugrid["mu"]);
  sigAR         = as<arma::vec>(mugrid["sig"]);
  // (repmu and zdm already computed above) 
  arma::mat denom(ar+qz,ar+qz,arma::fill::zeros);
  arma::mat num(ar+qz,1,arma::fill::zeros);
  for (int xm = 0; xm<M; xm++){
     double sig_xm = sigAR(xm);
     if (!std::isfinite(sig_xm) || sig_xm < 1e-14) continue;  // degenerate sigma: skip
     arma::vec y_tmp  = y - repmu*muAR(xm,0);
     arma::mat x_tmp  = x - repmu*muAR.submat(xm,1,xm,ar);
     arma::mat xz_tmp = join_rows(x_tmp,zdm);
     denom            = denom + (trans(xz_tmp)*diagmat(xi_t_T_AR.col(xm)/sig_xm)*xz_tmp);
     num              = num + (trans(xz_tmp)*diagmat(xi_t_T_AR.col(xm)/sig_xm)*y_tmp);
  }
  arma::vec prev_beta = join_vert(phi, betaZ);
  arma::vec beta_new;
  try {
    beta_new = arma::solve(denom, num, arma::solve_opts::allow_ugly);
    if (!beta_new.is_finite()) beta_new = prev_beta;
  } catch (...) {
    beta_new = prev_beta;  // Singular denom: keep previous beta
  }
  arma::vec phi_new   = beta_new.rows(0,ar-1);
  arma::vec betaZ_new = beta_new.rows(ar,ar+qz-1);
  // Stationarity check via companion matrix eigenvalues
  if (ar > 0) {
    arma::mat F_comp(ar, ar, arma::fill::zeros);
    F_comp.row(0) = phi_new.t();
    if (ar > 1) F_comp.submat(1, 0, ar-1, ar-2) = arma::eye(ar-1, ar-1);
    arma::cx_vec eigvals = arma::eig_gen(F_comp);
    if (arma::max(arma::abs(eigvals)) >= 1.0) {
      phi_new = phi;      // revert to previous (stationary) phi
      betaZ_new = betaZ;  // revert betaZ too (joint solve)
    }
  }
  // ----- Update estimates for variance
  // ----- Compute new residuals
  List mdl_tmp        = clone(mdl);
  mdl_tmp["phi"]      = phi_new;
  mdl_tmp["betaZ"]    = betaZ_new;
  arma::mat eps       = calcResid_MSARXmdl(mdl_tmp, muAR, k);
  // ----- Compute updated sigma
  arma::vec sigma = sig.head(1 + msvar*(k-1));  // Initialize from previous sigma
  arma::vec sigma_tmp(M, arma::fill::zeros);
  if (msvar == TRUE){
    for (int xk = 0 ; xk<M; xk++){
      sigma_tmp(xk) = sum(eps.col(xk)%eps.col(xk)%xi_t_T_AR.col(xk));
      sum_tmp(xk)   = sum(xi_t_T_AR.col(xk));
    }
    for (int xk = 1; xk<=k;xk++){
      double sig_denom = arma::as_scalar(sum(sum_tmp.rows(find(state_ind==xk))));
      if (sig_denom < 1e-6 * Tsize) continue;  // degenerate: keep previous sigma
      sigma(xk-1)   = arma::as_scalar(sum(sigma_tmp.rows(find(state_ind==xk)))/sig_denom);
    }
  }else{
    sigma = arma::sum(arma::sum(eps%eps%xi_t_T_AR,0),1)/Tsize;
  }
  // ----- Compute model residuals
  arma::mat residuals(Tsize, k);
  arma::mat residuals_tmp = eps%xi_t_T_AR;
  for (int xk = 1; xk<=k; xk++){
    residuals.col(xk-1) = arma::sum(residuals_tmp.cols(find(state_ind==xk)), 1);
  }
  residuals = arma::sum(residuals,1);
  // --------------- Organize output
  // ----- Apply the EM regime-variance floor (constrained ML; see get_sigma_floor).
  // Univariate: clipping the variance update at the floor is the exact
  // constrained M-step maximizer (the M-step objective is unimodal in sigma^2).
  {
    double sig_floor_uv = get_sigma_floor(mdl);
    if (sig_floor_uv > 0.0) sigma = arma::clamp(sigma, sig_floor_uv, arma::datum::inf);
  }
  // ----- Produce new theta vector
  arma::vec theta_new = join_vert(mu, phi_new);
  theta_new = join_vert(theta_new, betaZ_new);
  theta_new = join_vert(theta_new, sigma);
  theta_new = join_vert(theta_new, vectorise(P_new));
  // ----- Fill output List
  List Maxim_output;
  Maxim_output["theta"]     = theta_new;
  Maxim_output["mu"]        = mu;
  Maxim_output["phi"]       = phi_new;
  Maxim_output["betaZ"]     = betaZ_new;
  Maxim_output["beta"]      = beta_new;
  Maxim_output["sigma"]     = sigma;
  Maxim_output["P"]         = P_new;
  Maxim_output["pinf"]      = pinf;
  Maxim_output["residuals"] = eps;
  Maxim_output["resid"]     = residuals;
  return(Maxim_output);
}

// ==============================================================================
//' @title Maximization step of EM algorithm for Markov-switching vector autoregressive model
//' 
//' @description This function performs the maximization step of the Expectation Maximization algorithm for Markov-switching vector autoregressive model.
//' 
//' @param theta Vector of model parameters.
//' @param mdl List with model attributes.
//' @param MSloglik_output List with output from \code{\link{ExpectationM_MSVARmdl}}.
//' @param k Integer determining the number of regimes.
//'  
//' @return List with new maximized parameters.
//' 
//' @keywords internal
//' 
//' @export
// [[Rcpp::export]]
List EMaximization_MSVARmdl(arma::vec theta, List mdl, List MSloglik_output, int k){
  // ---------- Initialize parameters
  arma::mat y = mdl["y"];
  arma::mat x = mdl["x"];
  int ar      = mdl["p"];
  int q       = mdl["q"];
  bool msmu   = mdl["msmu"];
  bool msvar  = mdl["msvar"];
  // ----- Output from log-Likelihood calculation
  arma::mat P           = MSloglik_output["P"];
  arma::mat xi_t_T      = MSloglik_output["xi_t_T"];
  arma::mat xi_t_t      = MSloglik_output["xi_t_t"];
  arma::mat xi_tp1_t    = MSloglik_output["xi_tp1_t"];
  arma::mat xi_t_T_AR   = MSloglik_output["xi_t_T_AR"];
  arma::mat xi_t_t_AR   = MSloglik_output["xi_t_t_AR"];
  arma::mat xi_tp1_t_AR = MSloglik_output["xi_tp1_t_AR"];
  arma::vec state_ind   = MSloglik_output["state_ind"];
  arma::mat phi         = MSloglik_output["phi"];
  List sig              = MSloglik_output["sigma"];
  // Length of Time-Series (T) and number of variables
  int Tsize = y.n_rows;
  arma::mat repq(1,q, arma::fill::ones);
  // Number of regimes consistent with AR lags
  int M = pow(k, ar+1);
  // --------------- Update transition matrix (P) and limiting probabilities (pinf)
  // ----- Smoothed transition counts [Hamilton 1994, eq. 22.4.16; Krolzig 1997, eq. 6.14].
  // The extended state encodes (s_t, s_{t-1}, ...), so Pr(s_t = j, s_{t-1} = i | Y) is the
  // sum of xi_t_T_AR over extended states with date-t regime j and lag-1 regime i.
  arma::mat counts(k, k, arma::fill::zeros);
  for (int xn = 0; xn < M; xn++){
    int reg_i0 = xn % k;        // regime at date t
    int reg_i1 = (xn / k) % k;  // regime at date t-1
    counts(reg_i0, reg_i1) += arma::accu(xi_t_T_AR.submat(1, xn, Tsize-1, xn));
  }
  arma::mat P_new = P;  // degenerate regimes keep their previous column
  double P_bound = get_P_bound(mdl);
  for (int xk = 0; xk < k; xk++){
    double countsum = arma::accu(counts.col(xk));
    if (std::isfinite(countsum) && countsum >= 1e-6 * Tsize && counts.col(xk).is_finite()){
      arma::vec Pcol_uc = counts.col(xk) / countsum;
      if ((P_bound > 0.0) && arma::any(Pcol_uc < P_bound)){
        // exact constrained M-step for this column (see constrain_P_col)
        P_new.col(xk) = constrain_P_col(counts.col(xk), P_bound, P.col(xk));
      } else {
        P_new.col(xk) = Pcol_uc;
      }
    }
  }
  // Floating-point guard: keep P_new a valid column-stochastic matrix.
  P_new = arma::clamp(P_new, 0.0, 1.0);
  for (int xk = 0; xk < k; xk++) {
    double colsum = arma::accu(P_new.col(xk));
    if (colsum > 0) {
      P_new.col(xk) /= colsum;
    } else {
      P_new.col(xk) = P.col(xk);  // fallback: keep previous
    }
  }
  // Get limiting probabilities implied by transition matrix P_0
  arma::vec pinf = limP(P_new);
  // --------------- Update estimates for mu
  // Initialize from previous theta (used as fallback if linear system is ill-conditioned)
  int mu_len = q * (1 + msmu*(k-1));
  arma::mat mu = arma::reshape(theta.head(mu_len), q, 1 + msmu*(k-1)).t();
  arma::mat mu_k(k, q, arma::fill::zeros);
  // ----- Compute updated mu (Krolzig 1997, Table 9.19 — MSM GLS estimator)
  {
    // ỹ_t = y_t - sum_j A_j * y_{t-j}  [regime-independent, T×q]
    arma::mat ytilde = y;
    for (int j = 0; j < ar; j++){
      arma::mat Aj = phi.submat(0, q*j, q-1, q*(j+1)-1);  // q×q
      ytilde -= x.submat(0, q*j, Tsize-1, q*(j+1)-1) * Aj.t();
    }
    if (msmu == TRUE){
      // Build A_lsys (qk×qk) and b_vec (qk×1): normal equations for stacked mu vector
      // mu_flat = [μ₀(0),...,μ₀(q-1), μ₁(0),...,μ_{k-1}(q-1)]  (regime 0-indexed)
      arma::mat A_lsys(k*q, k*q, arma::fill::zeros);
      arma::vec b_vec(k*q, arma::fill::zeros);
      for (int xn = 0; xn < M; xn++){
        arma::mat sigma_inv = arma::solve(as<arma::mat>(sig[xn % k]), arma::eye(q,q), arma::solve_opts::allow_ugly);
        double wt_sum = arma::sum(xi_t_T_AR.col(xn));
        arma::vec ytilde_wt = ytilde.t() * xi_t_T_AR.col(xn);  // q×1
        for (int xm1 = 0; xm1 < k; xm1++){
          // Build X^μ_{n,xm1} [q×q]: I_q*I[i₀=xm1] - sum_j A_j*I[i_j=xm1]
          arma::mat Xmu1(q, q, arma::fill::zeros);
          if ((xn % k) == xm1) Xmu1 += arma::eye<arma::mat>(q, q);
          int kpow = k;
          for (int j = 1; j <= ar; j++){
            if (((xn / kpow) % k) == xm1)
              Xmu1 -= phi.submat(0, q*(j-1), q-1, q*j-1);
            kpow *= k;
          }
          b_vec.rows(xm1*q, (xm1+1)*q-1) += Xmu1.t() * sigma_inv * ytilde_wt;
          for (int xm2 = 0; xm2 < k; xm2++){
            arma::mat Xmu2(q, q, arma::fill::zeros);
            if ((xn % k) == xm2) Xmu2 += arma::eye<arma::mat>(q, q);
            int kpow2 = k;
            for (int j = 1; j <= ar; j++){
              if (((xn / kpow2) % k) == xm2)
                Xmu2 -= phi.submat(0, q*(j-1), q-1, q*j-1);
              kpow2 *= k;
            }
            A_lsys.submat(xm1*q, xm2*q, (xm1+1)*q-1, (xm2+1)*q-1) +=
              wt_sum * Xmu1.t() * sigma_inv * Xmu2;
          }
        }
      }
      // Solve A_lsys * mu_flat = b_vec; recover (k×q) matrix
      if (A_lsys.is_finite() && arma::rcond(A_lsys) > 1e-12){
        try {
          arma::vec mu_flat = arma::solve(A_lsys, b_vec, arma::solve_opts::allow_ugly);
          if (mu_flat.is_finite()){
            for (int xm = 0; xm < k; xm++)
              mu.row(xm) = mu_flat.rows(xm*q, (xm+1)*q-1).t();
          }
        } catch (...) { /* keep previous mu */ }
      }
      mu_k = mu;
    }else{
      // Single common mean C = I_q - sum_j A_j
      arma::mat C = arma::eye<arma::mat>(q, q);
      for (int j = 0; j < ar; j++) C -= phi.submat(0, q*j, q-1, q*(j+1)-1);
      // Pre-compute sigma-weighted sums over extended states
      arma::mat A_scalar(q, q, arma::fill::zeros);
      arma::mat b_scalar(q, 1, arma::fill::zeros);
      for (int xn = 0; xn < M; xn++){
        arma::mat sigma_inv = arma::solve(as<arma::mat>(sig[xn % k]), arma::eye(q,q), arma::solve_opts::allow_ugly);
        double wt_sum = arma::sum(xi_t_T_AR.col(xn));
        arma::mat ytilde_wt = ytilde.t() * xi_t_T_AR.col(xn);  // q×1
        A_scalar += wt_sum * C.t() * sigma_inv * C;
        b_scalar += C.t() * sigma_inv * ytilde_wt;
      }
      if (A_scalar.is_finite() && arma::rcond(A_scalar) > 1e-12){
        try {
          arma::mat mu_new = arma::solve(A_scalar, b_scalar, arma::solve_opts::allow_ugly);
          if (mu_new.is_finite()) mu.row(0) = mu_new.t();
        } catch (...) { /* keep previous mu */ }
      }
      for (int xm = 0; xm < k; xm++) mu_k.row(xm) = mu.row(0);
    }
  }
  // --------------- Update estimates for phi
  List mugrid = argrid_MSVARmdl_cpp(mu_k, sig, k, ar);
  List muAR   = mugrid["mu"];
  List sigAR  = mugrid["sig"];
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
    arma::mat sigma_m_inv = arma::solve(sigma_m, arma::eye(q,q), arma::solve_opts::allow_ugly);
    denom = denom + kron(sigma_m_inv,trans(xz_tmp)*diagmat(xi_t_T_AR.col(xm))*xz_tmp);
    num = num + kron(sigma_m_inv,trans(xz_tmp)*diagmat(xi_t_T_AR.col(xm)))*vectorise(y_tmp);
  }
  arma::mat phi_new;
  try {
    arma::vec phi_vec = arma::solve(denom, num, arma::solve_opts::allow_ugly);
    if (!phi_vec.is_finite()) {
      phi_new = trans(phi);  // keep previous (phi is q×(q*ar), phi_new is (q*ar)×q)
    } else {
      phi_new = reshape(phi_vec, q*ar, q);
    }
  } catch (...) {
    phi_new = trans(phi);  // Singular denom: keep previous phi
  }

  // Companion Mat
  arma::mat F_tmp = trans(phi_new);
  arma::mat diag_mat = arma::eye(q*(ar-1),q*(ar-1));
  arma::mat diag_zero(q*(ar-1), q, arma::fill::zeros);
  arma::mat Mn = join_rows(diag_mat,diag_zero);
  arma::mat F = join_cols(F_tmp,Mn);
  // Stationarity check via companion matrix eigenvalues
  arma::cx_vec eigvals = arma::eig_gen(F);
  if (arma::max(arma::abs(eigvals)) >= 1.0) {
    phi_new = trans(phi);  // revert to previous (stationary) phi
    F_tmp = trans(phi_new);
    F = join_cols(F_tmp, Mn);
  }
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
      arma::mat sigma_m_tmp;
      if (Tsize_k(xk) < 1e-6 * Tsize) {
        sigma_m_tmp = as<arma::mat>(sig[xk]);  // degenerate: keep previous sigma
      } else {
        sigma_m_tmp = as<arma::mat>(sigma[xk]);
        sigma_m_tmp = sigma_m_tmp / Tsize_k(xk);
        // PD regularization: eigenvalue floor proportional to trace
        arma::mat U_eig; arma::vec eigval;
        arma::eig_sym(eigval, U_eig, sigma_m_tmp);
        double eps_rel = 1e-6 * arma::trace(sigma_m_tmp) / q;
        if (eigval.min() < eps_rel) {
          eigval.elem(arma::find(eigval < eps_rel)).fill(eps_rel);
          sigma_m_tmp = U_eig * arma::diagmat(eigval) * U_eig.t();
        }
      }
      sigma[xk] = sigma_m_tmp;
      sigma_out.subvec(xk*sigN, xk*sigN+sigN-1) = covar_vech(sigma_m_tmp);
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
    sigma_out = covar_vech(sigma_m_tmp);
  }
  // ----- Apply the EM regime-variance floor (constrained ML; see get_sigma_floor)
  apply_sigma_floor_mv(sigma, sigma_out, get_sigma_floor(mdl), q, k, msvar);
  // ----- Compute model residuals 
  arma::mat residuals(Tsize, q, arma::fill::zeros); 
  for (int xk = 0; xk<M;xk++){
    arma::mat eps_tmp = eps[xk];
    residuals = residuals + eps_tmp%(xi_t_T_AR.col(xk)*repq);
  }
  // --------------- Organize output
  // ----- Produce new theta vector
  arma::vec theta_new = join_vert(vectorise(trans(mu)), vectorise(phi_new));
  theta_new = join_vert(theta_new, sigma_out);
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
//' @title Maximization step of EM algorithm for Markov-switching VARX model
//' 
//' @description This function performs the maximization step of the Expectation Maximization algorithm for Markov-switching VARX model.
//' 
//' @param theta Vector of model parameters.
//' @param mdl List with model attributes.
//' @param MSloglik_output List with output from \code{\link{ExpectationM_MSVARmdl}}.
//' @param k Integer determining the number of regimes.
//'  
//' @return List with new maximized parameters.
//' 
//' @keywords internal
//' 
//' @export
// [[Rcpp::export]]
List EMaximization_MSVARXmdl(arma::vec theta, List mdl, List MSloglik_output, int k){
  // ---------- Initialize parameters
  arma::mat y   = mdl["y"];
  arma::mat x   = mdl["x"];
  List mdl_con  = mdl["control"];
  arma::mat Z   = mdl_con["Z"];
  int ar        = mdl["p"];
  int q         = mdl["q"];
  bool msmu     = mdl["msmu"];
  bool msvar    = mdl["msvar"];
  // ----- Output from log-Likelihood calculation
  arma::mat P           = MSloglik_output["P"];
  arma::mat xi_t_T      = MSloglik_output["xi_t_T"];
  arma::mat xi_t_t      = MSloglik_output["xi_t_t"];
  arma::mat xi_tp1_t    = MSloglik_output["xi_tp1_t"];
  arma::mat xi_t_T_AR   = MSloglik_output["xi_t_T_AR"];
  arma::mat xi_t_t_AR   = MSloglik_output["xi_t_t_AR"];
  arma::mat xi_tp1_t_AR = MSloglik_output["xi_tp1_t_AR"];
  arma::vec state_ind   = MSloglik_output["state_ind"];
  arma::mat phi         = MSloglik_output["phi"];
  arma::mat betaZ       = MSloglik_output["betaZ"];
  List sig              = MSloglik_output["sigma"];
  // Length of Time-Series (T) and number of variables
  int Tsize = y.n_rows;
  arma::mat repq(1,q, arma::fill::ones);
  // Number of regimes consistent with AR lags
  int M                 = pow(k, ar+1);
  int qz                = Z.n_cols;
  Z                     = Z.rows(ar,Tsize+ar-1);
  arma::rowvec zbar     = arma::mean(Z,0);
  // --------------- Update transition matrix (P) and limiting probabilities (pinf)
  // ----- Smoothed transition counts [Hamilton 1994, eq. 22.4.16; Krolzig 1997, eq. 6.14].
  // The extended state encodes (s_t, s_{t-1}, ...), so Pr(s_t = j, s_{t-1} = i | Y) is the
  // sum of xi_t_T_AR over extended states with date-t regime j and lag-1 regime i.
  arma::mat counts(k, k, arma::fill::zeros);
  for (int xn = 0; xn < M; xn++){
    int reg_i0 = xn % k;        // regime at date t
    int reg_i1 = (xn / k) % k;  // regime at date t-1
    counts(reg_i0, reg_i1) += arma::accu(xi_t_T_AR.submat(1, xn, Tsize-1, xn));
  }
  arma::mat P_new = P;  // degenerate regimes keep their previous column
  double P_bound = get_P_bound(mdl);
  for (int xk = 0; xk < k; xk++){
    double countsum = arma::accu(counts.col(xk));
    if (std::isfinite(countsum) && countsum >= 1e-6 * Tsize && counts.col(xk).is_finite()){
      arma::vec Pcol_uc = counts.col(xk) / countsum;
      if ((P_bound > 0.0) && arma::any(Pcol_uc < P_bound)){
        // exact constrained M-step for this column (see constrain_P_col)
        P_new.col(xk) = constrain_P_col(counts.col(xk), P_bound, P.col(xk));
      } else {
        P_new.col(xk) = Pcol_uc;
      }
    }
  }
  // Floating-point guard: keep P_new a valid column-stochastic matrix.
  P_new = arma::clamp(P_new, 0.0, 1.0);
  for (int xk = 0; xk < k; xk++) {
    double colsum = arma::accu(P_new.col(xk));
    if (colsum > 0) {
      P_new.col(xk) /= colsum;
    } else {
      P_new.col(xk) = P.col(xk);  // fallback: keep previous
    }
  }
  // Get limiting probabilities implied by transition matrix P_0
  arma::vec pinf = limP(P_new);
  // --------------- Update estimates for mu
  // Initialize from previous theta (used as fallback if linear system is ill-conditioned)
  int mu_len = q * (1 + msmu*(k-1));
  arma::mat mu = arma::reshape(theta.head(mu_len), q, 1 + msmu*(k-1)).t();
  arma::mat mu_k(k, q, arma::fill::zeros);
  // Pre-compute zdm here so mu update can adjust ỹ_t for exogenous contribution
  arma::vec repmu(Tsize, arma::fill::ones);
  arma::mat zdm = (Z - repmu*zbar);
  // ----- Compute updated mu (Krolzig 1997, Table 9.19 — MSM GLS estimator)
  {
    // ỹ_t = y_t - zdm*betaZ - sum_j A_j * y_{t-j}  [regime-independent, T×q]
    arma::mat ytilde = y - zdm * betaZ;
    for (int j = 0; j < ar; j++){
      arma::mat Aj = phi.submat(0, q*j, q-1, q*(j+1)-1);
      ytilde -= x.submat(0, q*j, Tsize-1, q*(j+1)-1) * Aj.t();
    }
    if (msmu == TRUE){
      arma::mat A_lsys(k*q, k*q, arma::fill::zeros);
      arma::vec b_vec(k*q, arma::fill::zeros);
      for (int xn = 0; xn < M; xn++){
        arma::mat sigma_inv = arma::solve(as<arma::mat>(sig[xn % k]), arma::eye(q,q), arma::solve_opts::allow_ugly);
        double wt_sum = arma::sum(xi_t_T_AR.col(xn));
        arma::vec ytilde_wt = ytilde.t() * xi_t_T_AR.col(xn);  // q×1
        for (int xm1 = 0; xm1 < k; xm1++){
          arma::mat Xmu1(q, q, arma::fill::zeros);
          if ((xn % k) == xm1) Xmu1 += arma::eye<arma::mat>(q, q);
          int kpow = k;
          for (int j = 1; j <= ar; j++){
            if (((xn / kpow) % k) == xm1)
              Xmu1 -= phi.submat(0, q*(j-1), q-1, q*j-1);
            kpow *= k;
          }
          b_vec.rows(xm1*q, (xm1+1)*q-1) += Xmu1.t() * sigma_inv * ytilde_wt;
          for (int xm2 = 0; xm2 < k; xm2++){
            arma::mat Xmu2(q, q, arma::fill::zeros);
            if ((xn % k) == xm2) Xmu2 += arma::eye<arma::mat>(q, q);
            int kpow2 = k;
            for (int j = 1; j <= ar; j++){
              if (((xn / kpow2) % k) == xm2)
                Xmu2 -= phi.submat(0, q*(j-1), q-1, q*j-1);
              kpow2 *= k;
            }
            A_lsys.submat(xm1*q, xm2*q, (xm1+1)*q-1, (xm2+1)*q-1) +=
              wt_sum * Xmu1.t() * sigma_inv * Xmu2;
          }
        }
      }
      if (A_lsys.is_finite() && arma::rcond(A_lsys) > 1e-12){
        try {
          arma::vec mu_flat = arma::solve(A_lsys, b_vec, arma::solve_opts::allow_ugly);
          if (mu_flat.is_finite()){
            for (int xm = 0; xm < k; xm++)
              mu.row(xm) = mu_flat.rows(xm*q, (xm+1)*q-1).t();
          }
        } catch (...) { /* keep previous mu */ }
      }
      mu_k = mu;
    }else{
      arma::mat C = arma::eye<arma::mat>(q, q);
      for (int j = 0; j < ar; j++) C -= phi.submat(0, q*j, q-1, q*(j+1)-1);
      arma::mat A_scalar(q, q, arma::fill::zeros);
      arma::vec b_scalar(q, arma::fill::zeros);
      for (int xn = 0; xn < M; xn++){
        arma::mat sigma_inv = arma::solve(as<arma::mat>(sig[xn % k]), arma::eye(q,q), arma::solve_opts::allow_ugly);
        double wt_sum = arma::sum(xi_t_T_AR.col(xn));
        arma::vec ytilde_wt = ytilde.t() * xi_t_T_AR.col(xn);
        A_scalar += wt_sum * C.t() * sigma_inv * C;
        b_scalar += C.t() * sigma_inv * ytilde_wt;
      }
      if (A_scalar.is_finite() && arma::rcond(A_scalar) > 1e-12){
        try {
          arma::vec mu_new = arma::solve(A_scalar, b_scalar, arma::solve_opts::allow_ugly);
          if (mu_new.is_finite()) mu.row(0) = mu_new.t();
        } catch (...) { /* keep previous mu */ }
      }
      for (int xm = 0; xm < k; xm++) mu_k.row(xm) = mu.row(0);
    }
  }
  // ----- Update estimates for beta (phi & betaZ)
  List mugrid = argrid_MSVARmdl_cpp(mu_k, sig, k, ar);
  List muAR   = mugrid["mu"];
  List sigAR  = mugrid["sig"];
  // (repmu and zdm already computed above) 
  arma::mat denom(q*q*ar + qz*q, q*q*ar  + qz*q, arma::fill::zeros);
  arma::mat num(q*q*ar  + qz*q, 1, arma::fill::zeros);
  for (int xm = 0; xm<M; xm++){
    arma::mat muAR_tmp = muAR[xm]; 
    arma::mat y_tmp = y - repmu*trans(muAR_tmp.col(0));
    arma::mat x_tmp(Tsize, q*ar, arma::fill::zeros); 
    for (int xp = 0; xp<ar; xp++){
       x_tmp.submat(0,q*xp,Tsize-1,q*xp+q-1) = x.submat(0,q*xp,Tsize-1,q*xp+q-1) - repmu*trans(muAR_tmp.col(xp+1));
    }
    arma::mat xz_tmp = join_rows(x_tmp,zdm);
    arma::mat sigma_m = sigAR[xm];
    arma::mat sigma_m_inv = arma::solve(sigma_m, arma::eye(q,q), arma::solve_opts::allow_ugly);
    denom = denom + kron(sigma_m_inv,trans(xz_tmp)*diagmat(xi_t_T_AR.col(xm))*xz_tmp);
    num = num + kron(sigma_m_inv,trans(xz_tmp)*diagmat(xi_t_T_AR.col(xm)))*vectorise(y_tmp);
  }
  arma::mat beta_new;
  arma::mat phi_new;
  arma::mat betaZ_new;
  try {
    arma::vec beta_vec = arma::solve(denom, num, arma::solve_opts::allow_ugly);
    if (!beta_vec.is_finite()) {
      phi_new = trans(phi);  // keep previous (phi is q×(q*ar), phi_new is (q*ar)×q)
      betaZ_new = betaZ;     // keep previous
    } else {
      beta_new  = reshape(beta_vec, q*ar+qz, q);
      phi_new   = beta_new.rows(0, q*ar-1);
      betaZ_new = beta_new.rows(q*ar, q*ar+qz-1);
    }
  } catch (...) {
    phi_new = trans(phi);  // Singular denom: keep previous
    betaZ_new = betaZ;     // keep previous
  }
  // Companion Mat
  arma::mat F_tmp = trans(phi_new);
  arma::mat diag_mat = arma::eye(q*(ar-1),q*(ar-1));
  arma::mat diag_zero(q*(ar-1), q, arma::fill::zeros);
  arma::mat Mn = join_rows(diag_mat,diag_zero);
  arma::mat F = join_cols(F_tmp,Mn);
  // Stationarity check via companion matrix eigenvalues
  arma::cx_vec eigvals = arma::eig_gen(F);
  if (arma::max(arma::abs(eigvals)) >= 1.0) {
    phi_new = trans(phi);  // revert to previous (stationary) phi
    betaZ_new = betaZ;     // revert betaZ too (joint solve)
    F_tmp = trans(phi_new);
    F = join_cols(F_tmp, Mn);
  }
  // --------------- Update estimates for variance
  // ----- Compute new residuals
  List mdl_tmp    = clone(mdl);
  mdl_tmp["phi"]    = trans(phi_new);
  mdl_tmp["betaZ"]  = betaZ_new;
  List eps = calcResid_MSVARXmdl(mdl_tmp, muAR, k);
  // ----- Compute updated sigma
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
      arma::mat sigma_m_tmp;
      if (Tsize_k(xk) < 1e-6 * Tsize) {
        sigma_m_tmp = as<arma::mat>(sig[xk]);  // degenerate: keep previous sigma
      } else {
        sigma_m_tmp = as<arma::mat>(sigma[xk]);
        sigma_m_tmp = sigma_m_tmp / Tsize_k(xk);
        // PD regularization: eigenvalue floor proportional to trace
        arma::mat U_eig; arma::vec eigval;
        arma::eig_sym(eigval, U_eig, sigma_m_tmp);
        double eps_rel = 1e-6 * arma::trace(sigma_m_tmp) / q;
        if (eigval.min() < eps_rel) {
          eigval.elem(arma::find(eigval < eps_rel)).fill(eps_rel);
          sigma_m_tmp = U_eig * arma::diagmat(eigval) * U_eig.t();
        }
      }
      sigma[xk] = sigma_m_tmp;
      sigma_out.subvec(xk*sigN, xk*sigN+sigN-1) = covar_vech(sigma_m_tmp);
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
    sigma_out = covar_vech(sigma_m_tmp);
  }
  // ----- Apply the EM regime-variance floor (constrained ML; see get_sigma_floor)
  apply_sigma_floor_mv(sigma, sigma_out, get_sigma_floor(mdl), q, k, msvar);
  // ----- Compute model residuals
  arma::mat residuals(Tsize, q, arma::fill::zeros);
  for (int xk = 0; xk<M;xk++){
    arma::mat eps_tmp = eps[xk];
    residuals = residuals + eps_tmp%(xi_t_T_AR.col(xk)*repq);
  }
  // --------------- Organize output
  // ----- Produce new theta vector
  arma::vec theta_new = join_vert(vectorise(trans(mu)), vectorise(phi_new));
  theta_new = join_vert(theta_new, vectorise(betaZ_new));
  theta_new = join_vert(theta_new, sigma_out);
  theta_new = join_vert(theta_new, vectorise(P_new));
  // ----- Fill output List
  List Maxim_output;
  Maxim_output["theta"]         = theta_new;
  Maxim_output["mu"]            = mu;
  Maxim_output["phi"]           = trans(phi_new);
  Maxim_output["betaZ"]         = betaZ_new;
  Maxim_output["sigma"]         = sigma;
  Maxim_output["P"]             = P_new;
  Maxim_output["pinf"]          = pinf;
  Maxim_output["companion_mat"] = F;
  Maxim_output["residuals"]     = eps;
  Maxim_output["resid"]         = residuals;
  return(Maxim_output);
}
// ==============================================================================
//' @title EM algorithm iteration for Hidden Markov model
//' 
//' @description This function performs the one iteration (E-step and M-step) of the Expectation Maximization algorithm for a Hidden Markov model.
//' 
//' @param mdl List with model attributes.
//' @param EMest_output List with attributes from previous iteration.
//' @param k Integer determining the number of regimes.
//' 
//' @return List with attributes from new iteration.
//' 
//' @keywords internal
//' 
//' @export
// [[Rcpp::export]]
List EMiter_HMmdl(List mdl, List EMest_output, int k){
  arma::vec theta = EMest_output["theta"];
  // ---------- Expectation
  List HMMloglik_output = ExpectationM_HMmdl(theta, mdl, k);
  // ---------- Maximization
  List Maxim_output     = EMaximization_HMmdl(theta, mdl, HMMloglik_output, k);
  // ---------- Organize output
  arma::vec thl(3, arma::fill::zeros);
  double loglik_new = HMMloglik_output["logLike"];
  double loglik_old = EMest_output["logLike"];
  thl(0) = loglik_new;
  thl(1) = loglik_new - loglik_old;
  arma::mat P_n     = Maxim_output["P"];
  arma::vec theta_n = Maxim_output["theta"];
  arma::vec delta_rel = arma::abs(theta_n - theta) / (arma::abs(theta) + 1e-10);
  double deltath    = arma::max(delta_rel);
  if (!std::isfinite(deltath)) deltath = arma::datum::inf;
  thl(2)            = deltath;
  List EM_output;
  EM_output["theta"]      = theta_n; 
  EM_output["P"]          = P_n;
  EM_output["pinf"]       = Maxim_output["pinf"];
  EM_output["mu"]         = Maxim_output["mu"];
  EM_output["betaZ"]      = Maxim_output["betaZ"];
  EM_output["sigma"]      = Maxim_output["sigma"]; 
  EM_output["St"]         = HMMloglik_output["xi_t_T"];
  EM_output["eta"]        = HMMloglik_output["eta"];
  EM_output["thl"]        = thl;
  EM_output["deltath"]    = deltath;
  EM_output["logLike"]    = loglik_new;
  EM_output["residuals"]  = Maxim_output["residuals"];
  EM_output["resid"]      = Maxim_output["resid"];
  return(EM_output);
}

// ==============================================================================
//' @title EM algorithm iteration for Markov-switching autoregressive model
//' 
//' @description This function performs the one iteration (E-step and M-step) of the Expectation Maximization algorithm for Markov-switching autoregressive model.
//' 
//' @param mdl List with model attributes.
//' @param EMest_output List with attributes from previous iteration.
//' @param k Integer determining the number of regimes.
//' 
//' @return List with attributes from new iteration.
//' 
//' @keywords internal
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
  thl(0)            = loglik_new;
  thl(1)            = loglik_new - loglik_old;
  arma::mat P_n     = Maxim_output["P"];
  arma::vec theta_n = Maxim_output["theta"];
  arma::vec delta_rel = arma::abs(theta_n - theta) / (arma::abs(theta) + 1e-10);
  double deltath    = arma::max(delta_rel);
  if (!std::isfinite(deltath)) deltath = arma::datum::inf;
  thl(2) = deltath;
  List EM_output;
  EM_output["theta"]      = theta_n; 
  EM_output["mu"]         = Maxim_output["mu"];
  EM_output["phi"]        = Maxim_output["phi"]; 
  EM_output["beta"]       = Maxim_output["beta"]; 
  EM_output["sigma"]      = Maxim_output["sigma"]; 
  EM_output["P"]          = P_n;
  EM_output["pinf"]       = Maxim_output["pinf"];
  EM_output["St"]         = MSloglik_output["xi_t_T"];
  EM_output["eta"]        = MSloglik_output["eta"];
  EM_output["thl"]        = thl;
  EM_output["deltath"]    = deltath;
  EM_output["logLike"]    = loglik_new;
  EM_output["residuals"]  = Maxim_output["residuals"];
  EM_output["resid"]      = Maxim_output["resid"];
  return(EM_output);
}

// ==============================================================================
//' @title EM algorithm iteration for Markov-switching ARX model
//' 
//' @description This function performs the one iteration (E-step and M-step) of the Expectation Maximization algorithm for Markov-switching ARX model.
//' 
//' @param mdl List with model attributes.
//' @param EMest_output List with attributes from previous iteration.
//' @param k Integer determining the number of regimes.
//' 
//' @return List with attributes from new iteration.
//' 
//' @keywords internal
//' 
//' @export
// [[Rcpp::export]]
List EMiter_MSARXmdl(List mdl, List EMest_output, int k){
  arma::vec theta = EMest_output["theta"];
  // ---------- Expectation
  List MSloglik_output = ExpectationM_MSARXmdl(theta, mdl, k);
  // ---------- Maximization
  List Maxim_output = EMaximization_MSARXmdl(theta, mdl, MSloglik_output, k);
  // ---------- Organize output
  arma::vec thl(3, arma::fill::zeros);
  double loglik_new = MSloglik_output["logLike"];
  double loglik_old = EMest_output["logLike"];
  thl(0)            = loglik_new;
  thl(1)            = loglik_new - loglik_old;
  arma::mat P_n     = Maxim_output["P"];
  arma::vec theta_n = Maxim_output["theta"];
  arma::vec delta_rel = arma::abs(theta_n - theta) / (arma::abs(theta) + 1e-10);
  double deltath    = arma::max(delta_rel);
  if (!std::isfinite(deltath)) deltath = arma::datum::inf;
  thl(2) = deltath;
  List EM_output;
  EM_output["theta"]      = theta_n; 
  EM_output["mu"]         = Maxim_output["mu"];
  EM_output["phi"]        = Maxim_output["phi"]; 
  EM_output["betaZ"]      = Maxim_output["betaZ"]; 
  EM_output["beta"]       = Maxim_output["beta"]; 
  EM_output["sigma"]      = Maxim_output["sigma"]; 
  EM_output["P"]          = P_n;
  EM_output["pinf"]       = Maxim_output["pinf"];
  EM_output["St"]         = MSloglik_output["xi_t_T"];
  EM_output["eta"]        = MSloglik_output["eta"];
  EM_output["thl"]        = thl;
  EM_output["deltath"]    = deltath;
  EM_output["logLike"]    = loglik_new;
  EM_output["residuals"]  = Maxim_output["residuals"];
  EM_output["resid"]      = Maxim_output["resid"];
  return(EM_output);
}


// ==============================================================================
//' @title EM algorithm iteration for Markov-switching vector autoregressive model
//' 
//' @description This function performs the one iteration (E-step and M-step) of the Expectation Maximization algorithm for Markov-switching vector autoregressive model.
//' 
//' @param mdl List with model attributes.
//' @param EMest_output List with attributes from previous iteration.
//' @param k Integer determining the number of regimes.
//' 
//' @return List with attributes from new iteration.
//' 
//' @keywords internal
//' 
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
  arma::vec delta_rel = arma::abs(theta_n - theta) / (arma::abs(theta) + 1e-10);
  double deltath = arma::max(delta_rel);
  if (!std::isfinite(deltath)) deltath = arma::datum::inf;
  thl(2) = deltath;
  List EM_output;
  EM_output["theta"]          = theta_n; 
  EM_output["P"]              = P_n;
  EM_output["pinf"]           = Maxim_output["pinf"];
  EM_output["mu"]             = Maxim_output["mu"];
  EM_output["phi"]            = Maxim_output["phi"]; 
  EM_output["sigma"]          = Maxim_output["sigma"]; 
  EM_output["companion_mat"]  = Maxim_output["companion_mat"]; 
  EM_output["St"]             = MSVARloglik_output["xi_t_T"];
  EM_output["eta"]            = MSVARloglik_output["eta"];
  EM_output["thl"]            = thl;
  EM_output["deltath"]        = deltath;
  EM_output["logLike"]        = loglik_new;
  EM_output["residuals"]      = Maxim_output["residuals"];
  EM_output["resid"]          = Maxim_output["resid"];
  return(EM_output);
}

// ==============================================================================
//' @title EM algorithm iteration for Markov-switching VARX model
//' 
//' @description This function performs the one iteration (E-step and M-step) of the Expectation Maximization algorithm for Markov-switching VARX model.
//' 
//' @param mdl List with model attributes.
//' @param EMest_output List with attributes from previous iteration.
//' @param k Integer determining the number of regimes.
//' 
//' @return List with attributes from new iteration.
//' 
//' @keywords internal
//' 
//' @export
// [[Rcpp::export]]
List EMiter_MSVARXmdl(List mdl, List EMest_output, int k){
  arma::vec theta = EMest_output["theta"];
  // ---------- Expectation
  List MSVARloglik_output = ExpectationM_MSVARXmdl(theta, mdl, k);
  // ---------- Maximization
  List Maxim_output = EMaximization_MSVARXmdl(theta, mdl, MSVARloglik_output, k);
  // ---------- Organize output
  arma::vec thl(3, arma::fill::zeros);
  double loglik_new = MSVARloglik_output["logLike"];
  double loglik_old = EMest_output["logLike"];
  thl(0) = loglik_new;
  thl(1) = loglik_new - loglik_old;
  arma::mat P_n = Maxim_output["P"];
  arma::vec theta_n = Maxim_output["theta"];
  arma::vec delta_rel = arma::abs(theta_n - theta) / (arma::abs(theta) + 1e-10);
  double deltath = arma::max(delta_rel);
  if (!std::isfinite(deltath)) deltath = arma::datum::inf;
  thl(2) = deltath;
  List EM_output;
  EM_output["theta"]          = theta_n; 
  EM_output["P"]              = P_n;
  EM_output["pinf"]           = Maxim_output["pinf"];
  EM_output["mu"]             = Maxim_output["mu"];
  EM_output["phi"]            = Maxim_output["phi"]; 
  EM_output["betaZ"]          = Maxim_output["betaZ"];
  EM_output["sigma"]          = Maxim_output["sigma"]; 
  EM_output["companion_mat"]  = Maxim_output["companion_mat"]; 
  EM_output["St"]             = MSVARloglik_output["xi_t_T"];
  EM_output["eta"]            = MSVARloglik_output["eta"];
  EM_output["thl"]            = thl;
  EM_output["deltath"]        = deltath;
  EM_output["logLike"]        = loglik_new;
  EM_output["residuals"]      = Maxim_output["residuals"];
  EM_output["resid"]          = Maxim_output["resid"];
  return(EM_output);
}

// ==============================================================================
// Internal helper: EM convergence test. Returns true when the chosen criterion
// `conv` is satisfied. Not exported; used by the *_em drivers below.
//   thl = [logLike_new, logLike_new - logLike_old, deltath] from EMiter_*.
//   conv: "theta"    relative parameter change (Hamilton 1994; Krolzig 1997),
//         "loglik"   relative log-likelihood change, floored (Krolzig 1997, 6.4.4),
//         "both"     theta AND loglik (Krolzig 1997 default),
//         "loglik-A" / "both-A"  Aitken-accelerated log-likelihood
//                    (Bohning et al. 1994; McLachlan & Krishnan 2008).
//   thtol = parameter tolerance; ltol = log-likelihood tolerance.
//   inc_prev/linf_prev/aitken_streak carry Aitken state across iterations.
bool em_converged(const std::string& conv, double deltath, const arma::vec& thl,
                  double thtol, double ltol, int it,
                  double& inc_prev, double& linf_prev, int& aitken_streak){
  double L_new = thl(0);          // log-likelihood at iteration it
  double inc   = thl(1);          // L_it - L_{it-1}
  double L_old = L_new - inc;     // log-likelihood at iteration it-1
  // --- relative parameter change (legacy; valid from it == 1) ---
  bool conv_theta = (deltath <= thtol);
  // --- relative log-likelihood change, floored (Krolzig 1997, 6.4.4) ---
  //     at it == 1, L_{it-1} is the 0.0 placeholder, so require it >= 2.
  bool conv_ll = false;
  if (it >= 2){
    double denom = std::fabs(L_old);
    if (denom < 1.0) denom = 1.0;
    conv_ll = ((std::fabs(inc) / denom) <= ltol);
  }
  // --- Aitken-accelerated log-likelihood (Bohning et al. 1994) ---
  //     needs it >= 3, a valid prior increment, rate a in (0,1) bounded away
  //     from 1, and two consecutive valid Aitken steps (so linf_prev is the
  //     immediately preceding extrapolated limit).
  bool conv_llA = false;
  double linf_cur = arma::datum::nan;
  bool aitken_ok = false;
  if (it >= 3 && std::isfinite(inc_prev) &&
      std::fabs(inc_prev) > 1e-12 * (std::fabs(L_old) + 1.0)){
    double a = inc / inc_prev;
    if (a > 0.0 && a < 1.0 && (1.0 - a) > 1e-8){
      linf_cur = L_old + inc / (1.0 - a);
      if (std::isfinite(linf_cur)){
        aitken_ok = true;
        if (aitken_streak >= 1 && std::isfinite(linf_prev)){
          double denomA = std::fabs(linf_prev);
          if (denomA < 1.0) denomA = 1.0;
          conv_llA = ((std::fabs(linf_cur - linf_prev) / denomA) <= ltol);
        }
      }
    }
  }
  // --- combine per mode ---
  bool converged;
  if      (conv == "theta")    converged = conv_theta;
  else if (conv == "loglik")   converged = conv_ll;
  else if (conv == "both")     converged = conv_theta && conv_ll;
  else if (conv == "loglik-A") converged = conv_llA;
  else if (conv == "both-A")   converged = conv_theta && conv_llA;
  else                         converged = conv_theta;   // unknown -> legacy
  // --- carry state to next iteration ---
  if (aitken_ok){ linf_prev = linf_cur; aitken_streak += 1; }
  else          { aitken_streak = 0; }
  if (it >= 2)  { inc_prev = inc; }   // never store the placeholder it == 1 increment
  return converged;
}

// ==============================================================================
//' @title Estimation of Hidden Markov model by EM Algorithm
//' 
//' @description Estimate Hidden Markov model by EM algorithm. This function is used by \code{\link{HMmdl}} which organizes the output and takes raw data as input.
//' 
//' @param theta_0 vector with initial values for parameters.
//' @param mdl List with model attributes.
//' @param k Integer determining the number of regimes.
//' @param optim_options List with optimization options.
//' 
//' @return List with model results.
//' 
//' @references Dempster, A. P., N. M. Laird, and D. B. Rubin. 1977. “Maximum Likelihood from Incomplete Data via the EM Algorithm.” \emph{Journal of the Royal Statistical Society}. Series B 39 (1): 1–38.
//' 
//' @keywords internal
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
  // ----- Start from a feasible point w.r.t. the EM variance floor
  theta_0 = project_theta_sigma_floor(theta_0, mdl, k);
  theta_0 = project_theta_P_bound(theta_0, mdl, k);
  EMest_output["theta"] = theta_0;  // initial values for parameters
  EMest_output["logLike"] = 0.0;     // initial log-likelihood
  EMest_output = EMiter_HMmdl(mdl, EMest_output, k);
  int it = 1;
  double deltath = EMest_output["deltath"];
  arma::vec thl = EMest_output["thl"];
  // convergence controls (ltol/conv/em_best_iterate may be absent for legacy direct callers)
  double ltol = optim_options.containsElementNamed("ltol") ? Rcpp::as<double>(optim_options["ltol"]) : 1e-7;
  std::string conv = optim_options.containsElementNamed("conv") ? Rcpp::as<std::string>(optim_options["conv"]) : "theta";
  bool best_iterate = optim_options.containsElementNamed("em_best_iterate") ? Rcpp::as<bool>(optim_options["em_best_iterate"]) : false;
  double inc_prev = arma::datum::nan, linf_prev = arma::datum::nan;
  int aitken_streak = 0;
  // Track the best likelihood visited: the ergodic filter initialization makes the
  // transition-probability update a generalized EM step, so single iterations can lose
  // likelihood; with em_best_iterate the best visited iterate is returned.
  arma::vec best_theta = theta_0;
  double best_ll = thl(0);
  double ll_seen = thl(0);
  int descents = 0;
  bool converged = em_converged(conv, deltath, thl, thtol, ltol, it, inc_prev, linf_prev, aitken_streak);
  // iterate until convergence criterion is met or max number of iterations (maxit)
  while ((it<maxit) & (!converged)){
    arma::vec theta_in = EMest_output["theta"];
    EMest_output = EMiter_HMmdl(mdl, EMest_output, k);
    deltath = EMest_output["deltath"];
    thl = Rcpp::as<arma::vec>(EMest_output["thl"]);
    if (std::isfinite(thl(0))){
      if (std::isfinite(ll_seen) && thl(0) < ll_seen - 1e-9) descents++;
      if (!std::isfinite(best_ll) || thl(0) > best_ll){ best_ll = thl(0); best_theta = theta_in; }
      ll_seen = thl(0);
    }
    it = it + 1;
    converged = em_converged(conv, deltath, thl, thtol, ltol, it, inc_prev, linf_prev, aitken_streak);
  }
  if (best_iterate){
    // One E-step at the returned iterate: logLike/St/residuals correspond to it exactly.
    arma::vec theta_last = EMest_output["theta"];
    List E_fin = ExpectationM_HMmdl(theta_last, mdl, k);
    double ll_fin = E_fin["logLike"];
    if (std::isfinite(ll_seen) && std::isfinite(ll_fin) && ll_fin < ll_seen - 1e-9) descents++;
    arma::vec theta_w = theta_last; List E_w = E_fin;
    if (std::isfinite(best_ll) && !(std::isfinite(ll_fin) && ll_fin >= best_ll)){
      theta_w = best_theta; E_w = ExpectationM_HMmdl(best_theta, mdl, k);
    }
    EMest_output["theta"]     = theta_w;
    EMest_output["mu"]        = E_w["mu"];
    EMest_output["betaZ"]     = E_w["betaZ"];
    EMest_output["sigma"]     = E_w["sigma"];
    EMest_output["P"]         = E_w["P"];
    EMest_output["pinf"]      = E_w["pinf"];
    EMest_output["St"]        = E_w["xi_t_T"];
    EMest_output["eta"]       = E_w["eta"];
    EMest_output["logLike"]   = E_w["logLike"];
    EMest_output["residuals"] = E_w["residuals"];
    EMest_output["resid"]     = E_w["resid"];
  }
  EMest_output["descents"] = descents;
  EMest_output["iterations"] = it;  // EM iterations performed, including the pre-loop iteration
  EMest_output["converged"] = converged;
  return(EMest_output);
}



// ==============================================================================
//' @title Estimation of Markov-switching autoregressive model by EM Algorithm 
//' 
//' @description Estimate Markov-switching autoregressive model by EM algorithm. This function is used by \code{\link{MSARmdl}} which organizes the output and takes raw data as input.
//' 
//' @param theta_0 vector with initial values for parameters.
//' @param mdl List with model attributes.
//' @param k Integer determining the number of regimes.
//' @param optim_options List with optimization options.
//' 
//' @return List with model results.
//' 
//' @keywords internal
//' 
//' @references Dempster, A. P., N. M. Laird, and D. B. Rubin. 1977. “Maximum Likelihood from Incomplete Data via the EM Algorithm.” \emph{Journal of the Royal Statistical Society}. Series B 39 (1): 1–38.
//' @references Hamilton, James D. 1990. “Analysis of time series subject to changes in regime.” \emph{Journal of econometrics}, 45 (1-2): 39–70.
//' 
//' @export
// [[Rcpp::export]]
List MSARmdl_em(arma::vec theta_0, List mdl, int k, List optim_options){
  // ---------- Get optimization options
  int maxit     = optim_options["maxit"];
  double thtol  = optim_options["thtol"];
  // ---------- Begin EM Algo
  // ----- Initial values
  List EMest_output;
  // ----- Start from a feasible point w.r.t. the EM variance floor
  theta_0 = project_theta_sigma_floor(theta_0, mdl, k);
  theta_0 = project_theta_P_bound(theta_0, mdl, k);
  EMest_output["theta"] = theta_0;  // initial values for parameters
  EMest_output["logLike"] = 0.0;     // initial log-likelihood
  EMest_output = EMiter_MSARmdl(mdl, EMest_output, k);
  int it = 1;
  double deltath = EMest_output["deltath"];
  arma::vec thl = EMest_output["thl"];
  // convergence controls (ltol/conv/em_best_iterate may be absent for legacy direct callers)
  double ltol = optim_options.containsElementNamed("ltol") ? Rcpp::as<double>(optim_options["ltol"]) : 1e-7;
  std::string conv = optim_options.containsElementNamed("conv") ? Rcpp::as<std::string>(optim_options["conv"]) : "theta";
  bool best_iterate = optim_options.containsElementNamed("em_best_iterate") ? Rcpp::as<bool>(optim_options["em_best_iterate"]) : false;
  double inc_prev = arma::datum::nan, linf_prev = arma::datum::nan;
  int aitken_streak = 0;
  // Track the best likelihood visited: the ergodic filter initialization makes the
  // transition-probability update a generalized EM step, so single iterations can lose
  // likelihood; with em_best_iterate the best visited iterate is returned.
  arma::vec best_theta = theta_0;
  double best_ll = thl(0);
  double ll_seen = thl(0);
  int descents = 0;
  bool converged = em_converged(conv, deltath, thl, thtol, ltol, it, inc_prev, linf_prev, aitken_streak);
  // iterate until convergence criterion is met or max number of iterations (maxit)
  while ((it<maxit) & (!converged)){
    arma::vec theta_in = EMest_output["theta"];
    EMest_output = EMiter_MSARmdl(mdl, EMest_output, k);
    deltath = EMest_output["deltath"];
    thl = Rcpp::as<arma::vec>(EMest_output["thl"]);
    if (std::isfinite(thl(0))){
      if (std::isfinite(ll_seen) && thl(0) < ll_seen - 1e-9) descents++;
      if (!std::isfinite(best_ll) || thl(0) > best_ll){ best_ll = thl(0); best_theta = theta_in; }
      ll_seen = thl(0);
    }
    it = it + 1;
    converged = em_converged(conv, deltath, thl, thtol, ltol, it, inc_prev, linf_prev, aitken_streak);
  }
  if (best_iterate){
    // One E-step at the returned iterate: logLike/St/residuals correspond to it exactly.
    arma::vec theta_last = EMest_output["theta"];
    List E_fin = ExpectationM_MSARmdl(theta_last, mdl, k);
    double ll_fin = E_fin["logLike"];
    if (std::isfinite(ll_seen) && std::isfinite(ll_fin) && ll_fin < ll_seen - 1e-9) descents++;
    arma::vec theta_w = theta_last; List E_w = E_fin;
    if (std::isfinite(best_ll) && !(std::isfinite(ll_fin) && ll_fin >= best_ll)){
      theta_w = best_theta; E_w = ExpectationM_MSARmdl(best_theta, mdl, k);
    }
    EMest_output["theta"]     = theta_w;
    EMest_output["mu"]        = E_w["mu"];
    EMest_output["phi"]       = E_w["phi"];
    EMest_output["beta"]      = E_w["beta"];
    EMest_output["sigma"]     = E_w["sigma"];
    EMest_output["P"]         = E_w["P"];
    EMest_output["pinf"]      = E_w["pinf"];
    EMest_output["St"]        = E_w["xi_t_T"];
    EMest_output["eta"]       = E_w["eta"];
    EMest_output["logLike"]   = E_w["logLike"];
    EMest_output["residuals"] = E_w["residuals"];
    EMest_output["resid"]     = E_w["resid"];
  }
  EMest_output["descents"] = descents;
  EMest_output["iterations"] = it;  // EM iterations performed, including the pre-loop iteration
  EMest_output["converged"] = converged;
  return(EMest_output);
}

// ==============================================================================
//' @title Estimation of Markov-switching ARX model by EM Algorithm 
//' 
//' @description Estimate Markov-switching ARX model by EM algorithm. This function is used by \code{\link{MSARmdl}} which organizes the output and takes raw data as input.
//' 
//' @param theta_0 vector with initial values for parameters.
//' @param mdl List with model attributes.
//' @param k Integer determining the number of regimes.
//' @param optim_options List with optimization options.
//' 
//' @return List with model results.
//' 
//' @keywords internal
//' 
//' @references Dempster, A. P., N. M. Laird, and D. B. Rubin. 1977. “Maximum Likelihood from Incomplete Data via the EM Algorithm.” \emph{Journal of the Royal Statistical Society}. Series B 39 (1): 1–38.
//' @references Hamilton, James D. 1990. “Analysis of time series subject to changes in regime.” \emph{Journal of econometrics}, 45 (1-2): 39–70.
//' 
//' @export
// [[Rcpp::export]]
List MSARXmdl_em(arma::vec theta_0, List mdl, int k, List optim_options){
  // ---------- Get optimization options
  int maxit     = optim_options["maxit"];
  double thtol  = optim_options["thtol"];
  // ---------- Begin EM Algo
  // ----- Initial values
  List EMest_output;
  // ----- Start from a feasible point w.r.t. the EM variance floor
  theta_0 = project_theta_sigma_floor(theta_0, mdl, k);
  theta_0 = project_theta_P_bound(theta_0, mdl, k);
  EMest_output["theta"] = theta_0;  // initial values for parameters
  EMest_output["logLike"] = 0.0;     // initial log-likelihood
  EMest_output = EMiter_MSARXmdl(mdl, EMest_output, k);
  int it = 1;
  double deltath = EMest_output["deltath"];
  arma::vec thl = EMest_output["thl"];
  // convergence controls (ltol/conv/em_best_iterate may be absent for legacy direct callers)
  double ltol = optim_options.containsElementNamed("ltol") ? Rcpp::as<double>(optim_options["ltol"]) : 1e-7;
  std::string conv = optim_options.containsElementNamed("conv") ? Rcpp::as<std::string>(optim_options["conv"]) : "theta";
  bool best_iterate = optim_options.containsElementNamed("em_best_iterate") ? Rcpp::as<bool>(optim_options["em_best_iterate"]) : false;
  double inc_prev = arma::datum::nan, linf_prev = arma::datum::nan;
  int aitken_streak = 0;
  // Track the best likelihood visited: the ergodic filter initialization makes the
  // transition-probability update a generalized EM step, so single iterations can lose
  // likelihood; with em_best_iterate the best visited iterate is returned.
  arma::vec best_theta = theta_0;
  double best_ll = thl(0);
  double ll_seen = thl(0);
  int descents = 0;
  bool converged = em_converged(conv, deltath, thl, thtol, ltol, it, inc_prev, linf_prev, aitken_streak);
  // iterate until convergence criterion is met or max number of iterations (maxit)
  while ((it<maxit) & (!converged)){
    arma::vec theta_in = EMest_output["theta"];
    EMest_output = EMiter_MSARXmdl(mdl, EMest_output, k);
    deltath = EMest_output["deltath"];
    thl = Rcpp::as<arma::vec>(EMest_output["thl"]);
    if (std::isfinite(thl(0))){
      if (std::isfinite(ll_seen) && thl(0) < ll_seen - 1e-9) descents++;
      if (!std::isfinite(best_ll) || thl(0) > best_ll){ best_ll = thl(0); best_theta = theta_in; }
      ll_seen = thl(0);
    }
    it = it + 1;
    converged = em_converged(conv, deltath, thl, thtol, ltol, it, inc_prev, linf_prev, aitken_streak);
  }
  if (best_iterate){
    // One E-step at the returned iterate: logLike/St/residuals correspond to it exactly.
    arma::vec theta_last = EMest_output["theta"];
    List E_fin = ExpectationM_MSARXmdl(theta_last, mdl, k);
    double ll_fin = E_fin["logLike"];
    if (std::isfinite(ll_seen) && std::isfinite(ll_fin) && ll_fin < ll_seen - 1e-9) descents++;
    arma::vec theta_w = theta_last; List E_w = E_fin;
    if (std::isfinite(best_ll) && !(std::isfinite(ll_fin) && ll_fin >= best_ll)){
      theta_w = best_theta; E_w = ExpectationM_MSARXmdl(best_theta, mdl, k);
    }
    EMest_output["theta"]     = theta_w;
    EMest_output["mu"]        = E_w["mu"];
    EMest_output["phi"]       = E_w["phi"];
    EMest_output["betaZ"]     = E_w["betaZ"];
    EMest_output["beta"]      = E_w["beta"];
    EMest_output["sigma"]     = E_w["sigma"];
    EMest_output["P"]         = E_w["P"];
    EMest_output["pinf"]      = E_w["pinf"];
    EMest_output["St"]        = E_w["xi_t_T"];
    EMest_output["eta"]       = E_w["eta"];
    EMest_output["logLike"]   = E_w["logLike"];
    EMest_output["residuals"] = E_w["residuals"];
    EMest_output["resid"]     = E_w["resid"];
  }
  EMest_output["descents"] = descents;
  EMest_output["iterations"] = it;  // EM iterations performed, including the pre-loop iteration
  EMest_output["converged"] = converged;
  return(EMest_output);
}

// ==============================================================================
//' @title Estimation of Markov-switching vector autoregressive model by EM Algorithm 
//' 
//' @description Estimate Markov-switching vector autoregressive model by EM algorithm. This function is used by \code{\link{MSVARmdl}} which organizes the output and takes raw data as input.
//' 
//' @param theta_0 vector with initial values for parameters.
//' @param mdl List with model attributes.
//' @param k Integer determining the number of regimes.
//' @param optim_options List with optimization options.
//' 
//' @return List with model results.
//' 
//' @keywords internal
//' 
//' @references Dempster, A. P., N. M. Laird, and D. B. Rubin. 1977. “Maximum Likelihood from Incomplete Data via the EM Algorithm.” \emph{Journal of the Royal Statistical Society}. Series B 39 (1): 1–38.
//' @references Krolzig, Hans-Martin. 1997. “The markov-switching vector autoregressive model.”. Springer.
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
  // ----- Start from a feasible point w.r.t. the EM variance floor
  theta_0 = project_theta_sigma_floor(theta_0, mdl, k);
  theta_0 = project_theta_P_bound(theta_0, mdl, k);
  EMest_output["theta"] = theta_0;  // initial values for parameters
  EMest_output["logLike"] = 0.0;     // initial log-likelihood
  EMest_output = EMiter_MSVARmdl(mdl, EMest_output, k);
  int it = 1;
  double deltath = EMest_output["deltath"];
  arma::vec thl = EMest_output["thl"];
  // convergence controls (ltol/conv/em_best_iterate may be absent for legacy direct callers)
  double ltol = optim_options.containsElementNamed("ltol") ? Rcpp::as<double>(optim_options["ltol"]) : 1e-7;
  std::string conv = optim_options.containsElementNamed("conv") ? Rcpp::as<std::string>(optim_options["conv"]) : "theta";
  bool best_iterate = optim_options.containsElementNamed("em_best_iterate") ? Rcpp::as<bool>(optim_options["em_best_iterate"]) : false;
  double inc_prev = arma::datum::nan, linf_prev = arma::datum::nan;
  int aitken_streak = 0;
  // Track the best likelihood visited: the ergodic filter initialization makes the
  // transition-probability update a generalized EM step, so single iterations can lose
  // likelihood; with em_best_iterate the best visited iterate is returned.
  arma::vec best_theta = theta_0;
  double best_ll = thl(0);
  double ll_seen = thl(0);
  int descents = 0;
  bool converged = em_converged(conv, deltath, thl, thtol, ltol, it, inc_prev, linf_prev, aitken_streak);
  // iterate until convergence criterion is met or max number of iterations (maxit)
  while ((it<maxit) & (!converged)){
    arma::vec theta_in = EMest_output["theta"];
    EMest_output = EMiter_MSVARmdl(mdl, EMest_output, k);
    deltath = EMest_output["deltath"];
    thl = Rcpp::as<arma::vec>(EMest_output["thl"]);
    if (std::isfinite(thl(0))){
      if (std::isfinite(ll_seen) && thl(0) < ll_seen - 1e-9) descents++;
      if (!std::isfinite(best_ll) || thl(0) > best_ll){ best_ll = thl(0); best_theta = theta_in; }
      ll_seen = thl(0);
    }
    it = it + 1;
    converged = em_converged(conv, deltath, thl, thtol, ltol, it, inc_prev, linf_prev, aitken_streak);
  }
  if (best_iterate){
    // One E-step at the returned iterate: logLike/St/residuals correspond to it exactly.
    arma::vec theta_last = EMest_output["theta"];
    List E_fin = ExpectationM_MSVARmdl(theta_last, mdl, k);
    double ll_fin = E_fin["logLike"];
    if (std::isfinite(ll_seen) && std::isfinite(ll_fin) && ll_fin < ll_seen - 1e-9) descents++;
    arma::vec theta_w = theta_last; List E_w = E_fin;
    if (std::isfinite(best_ll) && !(std::isfinite(ll_fin) && ll_fin >= best_ll)){
      theta_w = best_theta; E_w = ExpectationM_MSVARmdl(best_theta, mdl, k);
    }
    EMest_output["theta"]     = theta_w;
    EMest_output["mu"]        = E_w["mu"];
    EMest_output["phi"]       = E_w["phi"];
    EMest_output["sigma"]     = E_w["sigma"];
    EMest_output["P"]         = E_w["P"];
    EMest_output["pinf"]      = E_w["pinf"];
    EMest_output["St"]        = E_w["xi_t_T"];
    EMest_output["eta"]       = E_w["eta"];
    EMest_output["logLike"]   = E_w["logLike"];
    EMest_output["residuals"] = E_w["residuals"];
    EMest_output["resid"]     = E_w["resid"];
  }
  EMest_output["descents"] = descents;
  EMest_output["iterations"] = it;  // EM iterations performed, including the pre-loop iteration
  EMest_output["converged"] = converged;
  return(EMest_output);
}

// ==============================================================================
//' @title Estimation of Markov-switching VARX model by EM Algorithm 
//' 
//' @description Estimate Markov-switching VARX model by EM algorithm. This function is used by \code{\link{MSVARmdl}} which organizes the output and takes raw data as input.
//' 
//' @param theta_0 vector with initial values for parameters.
//' @param mdl List with model attributes.
//' @param k Integer determining the number of regimes.
//' @param optim_options List with optimization options.
//' 
//' @return List with model results.
//' 
//' @keywords internal
//' 
//' @references Dempster, A. P., N. M. Laird, and D. B. Rubin. 1977. “Maximum Likelihood from Incomplete Data via the EM Algorithm.” \emph{Journal of the Royal Statistical Society}. Series B 39 (1): 1–38.
//' @references Krolzig, Hans-Martin. 1997. “The markov-switching vector autoregressive model.”. Springer.
//' 
//' @export
// [[Rcpp::export]]
List MSVARXmdl_em(arma::vec theta_0, List mdl, int k, List optim_options){
  // ---------- Get optimization options
  int maxit = optim_options["maxit"];
  double thtol = optim_options["thtol"];
  // ---------- Begin EM Algo
  // ----- Initial values
  List EMest_output;
  // ----- Start from a feasible point w.r.t. the EM variance floor
  theta_0 = project_theta_sigma_floor(theta_0, mdl, k);
  theta_0 = project_theta_P_bound(theta_0, mdl, k);
  EMest_output["theta"] = theta_0;  // initial values for parameters
  EMest_output["logLike"] = 0.0;     // initial log-likelihood
  EMest_output = EMiter_MSVARXmdl(mdl, EMest_output, k);
  int it = 1;
  double deltath = EMest_output["deltath"];
  arma::vec thl = EMest_output["thl"];
  // convergence controls (ltol/conv/em_best_iterate may be absent for legacy direct callers)
  double ltol = optim_options.containsElementNamed("ltol") ? Rcpp::as<double>(optim_options["ltol"]) : 1e-7;
  std::string conv = optim_options.containsElementNamed("conv") ? Rcpp::as<std::string>(optim_options["conv"]) : "theta";
  bool best_iterate = optim_options.containsElementNamed("em_best_iterate") ? Rcpp::as<bool>(optim_options["em_best_iterate"]) : false;
  double inc_prev = arma::datum::nan, linf_prev = arma::datum::nan;
  int aitken_streak = 0;
  // Track the best likelihood visited: the ergodic filter initialization makes the
  // transition-probability update a generalized EM step, so single iterations can lose
  // likelihood; with em_best_iterate the best visited iterate is returned.
  arma::vec best_theta = theta_0;
  double best_ll = thl(0);
  double ll_seen = thl(0);
  int descents = 0;
  bool converged = em_converged(conv, deltath, thl, thtol, ltol, it, inc_prev, linf_prev, aitken_streak);
  // iterate until convergence criterion is met or max number of iterations (maxit)
  while ((it<maxit) & (!converged)){
    arma::vec theta_in = EMest_output["theta"];
    EMest_output = EMiter_MSVARXmdl(mdl, EMest_output, k);
    deltath = EMest_output["deltath"];
    thl = Rcpp::as<arma::vec>(EMest_output["thl"]);
    if (std::isfinite(thl(0))){
      if (std::isfinite(ll_seen) && thl(0) < ll_seen - 1e-9) descents++;
      if (!std::isfinite(best_ll) || thl(0) > best_ll){ best_ll = thl(0); best_theta = theta_in; }
      ll_seen = thl(0);
    }
    it = it + 1;
    converged = em_converged(conv, deltath, thl, thtol, ltol, it, inc_prev, linf_prev, aitken_streak);
  }
  if (best_iterate){
    // One E-step at the returned iterate: logLike/St/residuals correspond to it exactly.
    arma::vec theta_last = EMest_output["theta"];
    List E_fin = ExpectationM_MSVARXmdl(theta_last, mdl, k);
    double ll_fin = E_fin["logLike"];
    if (std::isfinite(ll_seen) && std::isfinite(ll_fin) && ll_fin < ll_seen - 1e-9) descents++;
    arma::vec theta_w = theta_last; List E_w = E_fin;
    if (std::isfinite(best_ll) && !(std::isfinite(ll_fin) && ll_fin >= best_ll)){
      theta_w = best_theta; E_w = ExpectationM_MSVARXmdl(best_theta, mdl, k);
    }
    EMest_output["theta"]     = theta_w;
    EMest_output["mu"]        = E_w["mu"];
    EMest_output["phi"]       = E_w["phi"];
    EMest_output["betaZ"]     = E_w["betaZ"];
    EMest_output["sigma"]     = E_w["sigma"];
    EMest_output["P"]         = E_w["P"];
    EMest_output["pinf"]      = E_w["pinf"];
    EMest_output["St"]        = E_w["xi_t_T"];
    EMest_output["eta"]       = E_w["eta"];
    EMest_output["logLike"]   = E_w["logLike"];
    EMest_output["residuals"] = E_w["residuals"];
    EMest_output["resid"]     = E_w["resid"];
  }
  EMest_output["descents"] = descents;
  EMest_output["iterations"] = it;  // EM iterations performed, including the pre-loop iteration
  EMest_output["converged"] = converged;
  return(EMest_output);
 }

