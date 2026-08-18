#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;

// ==============================================================================
// Autoregressive moment grid and extended transition matrix (C++ ports)
// ==============================================================================
// C++ ports of the R helpers argrid_MSARmdl, argrid_MSVARmdl, and arP, called
// directly from the estimation code (the R versions remain the exported public
// API and the reference implementations for the identity tests). The grid is a
// base-k digit expansion of the extended state: digit j of 0-indexed state xn
// is (xn / k^j) % k, with digit 0 = the date-t regime; R's expand.grid orders
// states with its first column varying fastest, the same convention.

//' @title Autoregressive moment grid (C++)
//'
//' @description C++ port of \code{\link{argrid_MSARmdl}}: means, variances, and
//' regime indicator for each extended state. A length-1 \code{mu} or \code{sig}
//' applies to all regimes (the non-switching case).
//'
//' @param mu vector of regime means (length 1 or k).
//' @param sig vector of regime variances (length 1 or k).
//' @param k number of regimes.
//' @param ar number of autoregressive lags.
//'
//' @return List with M x (ar+1) matrix of means, M x 1 matrix of variances, and
//' M x 1 regime indicator, where M = k^(ar+1).
//'
//' @keywords internal
//'
//' @export
// [[Rcpp::export]]
List argrid_MSARmdl_cpp(arma::vec mu, arma::vec sig, int k, int ar){
  double Md = std::pow((double)k, ar+1);
  if (!(k >= 1) || !(ar >= 0) || Md > 16777216.0){
    Rcpp::stop("k^(ar+1) grid size is invalid or too large.");
  }
  int M = (int) Md;
  arma::mat muAR(M, ar+1);
  arma::mat sigAR(M, 1);
  arma::vec state_ind(M);
  for (int xn = 0; xn < M; xn++){
    int rem = xn;
    for (int j = 0; j <= ar; j++){
      int dj = rem % k;
      rem /= k;
      muAR(xn, j) = (mu.n_elem == 1) ? mu(0) : mu(dj);
    }
    int d0 = xn % k;
    sigAR(xn, 0) = (sig.n_elem == 1) ? sig(0) : sig(d0);
    state_ind(xn) = d0 + 1;
  }
  List musig_out;
  musig_out["mu"] = muAR;
  musig_out["sig"] = sigAR;
  musig_out["state_ind"] = state_ind;
  return(musig_out);
}

//' @title Vector autoregressive moment grid (C++)
//'
//' @description C++ port of \code{\link{argrid_MSVARmdl}}: mean matrices,
//' covariance matrices, and regime indicator for each extended state. Callers
//' pass \code{mu} as a full (k x q) matrix and \code{sigma} as a length-k list
//' (rows/elements duplicated by the caller when a parameter does not switch).
//'
//' @param mu (k x q) matrix of regime means.
//' @param sigma List of k (q x q) regime covariance matrices.
//' @param k number of regimes.
//' @param ar number of autoregressive lags.
//'
//' @return List with M mean matrices (q x (ar+1); column j is the regime mean at
//' lag j), M covariance matrices, and M x 1 regime indicator, M = k^(ar+1).
//'
//' @keywords internal
//'
//' @export
// [[Rcpp::export]]
List argrid_MSVARmdl_cpp(arma::mat mu, List sigma, int k, int ar){
  double Md = std::pow((double)k, ar+1);
  if (!(k >= 1) || !(ar >= 0) || Md > 16777216.0){
    Rcpp::stop("k^(ar+1) grid size is invalid or too large.");
  }
  int M = (int) Md;
  int q = mu.n_cols;
  if (sigma.size() < k){
    Rcpp::stop("sigma must be a list with at least k covariance matrices.");
  }
  List muAR(M);
  List sigAR(M);
  arma::vec state_ind(M);
  for (int xn = 0; xn < M; xn++){
    arma::mat mu_n(q, ar+1);
    int rem = xn;
    for (int j = 0; j <= ar; j++){
      int dj = rem % k;
      rem /= k;
      mu_n.col(j) = trans(mu.row(dj));
    }
    muAR[xn] = mu_n;
    int d0 = xn % k;
    sigAR[xn] = sigma[d0];
    state_ind(xn) = d0 + 1;
  }
  List musig_out;
  musig_out["mu"] = muAR;
  musig_out["sig"] = sigAR;
  musig_out["state_ind"] = state_ind;
  return(musig_out);
}

//' @title Extended transition matrix (C++)
//'
//' @description C++ port of \code{\link{arP}}: the M x M transition matrix of
//' the extended state (a shift register of the last ar+1 regimes). With
//' 0-indexed row r, column c, and d_j the base-k digit j:
//' out(r, c) = P(d_0(r), d_0(c)) when d_(j+1)(r) = d_j(c) for j = 0..ar-1, and
//' 0 otherwise (the oldest lag is marginalized). \code{ar = 0} returns P.
//'
//' @param P (k x k) column-stochastic transition matrix.
//' @param k number of regimes.
//' @param ar number of autoregressive lags.
//'
//' @return (M x M) extended transition matrix, M = k^(ar+1).
//'
//' @keywords internal
//'
//' @export
// [[Rcpp::export]]
arma::mat arP_cpp(arma::mat P, int k, int ar){
  if (ar <= 0) return(P);
  double Md = std::pow((double)k, ar+1);
  if (!(k >= 1) || Md > 16777216.0){
    Rcpp::stop("k^(ar+1) grid size is invalid or too large.");
  }
  int Mar = 1;
  for (int j = 0; j < ar; j++) Mar *= k;
  int M = Mar * k;
  arma::mat Pnew(M, M, arma::fill::zeros);
  for (int c = 0; c < M; c++){
    int d0c = c % k;
    int rbase = (c % Mar) * k;
    for (int d0r = 0; d0r < k; d0r++){
      Pnew(rbase + d0r, c) = P(d0r, d0c);
    }
  }
  return(Pnew);
}
