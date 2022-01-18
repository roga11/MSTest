#include <RcppArmadillo.h>
#include "methods.h"
//[[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;
// ==============================================================================
//' @title Fitting Autoregressive Model
//' 
//' @description This function estimates an autoregresive series
//' 
//' @param Y vector of observations
//' @param ar number of lags
//' @param intercept bool determines whether or not to include constant in estimation. Default is TRUE.
//' 
//' @return List with model characteristics which include
//' - y: transformed process
//' - X: matrix of lagged observations (with or without vector of 1s depending on const=1 or const=0)
//' - x: matrix of lagged observations without vector of 1s
//' - n: number of observations (T-ar)
//' - ar: number of autoregressive parameters
//' - coef: coefficient estimates. This is the same as phi is const=0
//' - phi: autoregressive coefficient estimates
//' - mu: the mean of the process
//' - residuals: vectot of residuals
//' - stdev: standard deviations
//' - se: standard errors of parameter estimtes
//' - logLike: the log-likelihood 
//' 
//' @export
// [[Rcpp::export]]
List ARmdl(arma::vec Y, int ar, bool intercept = 1){
  // ----- transform data
  List lagged_vals = ts_lagged(Y, ar);
  arma::vec y = lagged_vals["y"];
  arma::mat x = lagged_vals["X"];
  arma::mat X = lagged_vals["X"];
  arma::vec b0 (ar + intercept, arma::fill::zeros);
  arma::vec phi (ar , arma::fill::zeros);
  int n = y.n_elem;
  int npar = ar;
  // ----- estimate model
  if (intercept==TRUE){
    arma::mat onevec(n, 1, arma::fill::ones);
    X = join_rows(onevec,X);
    npar = npar + 1;
    //b0  = inv(trans(X)*X)*trans(X)*y;
    b0  = solve(trans(X)*X,trans(X))*y;
    phi = b0.subvec(1,ar);
  }else{
    //b0  = inv(trans(X)*X)*trans(X)*y;
    b0  = solve(trans(X)*X,trans(X))*y;
    phi = b0;
  }
  // ----- Obtain variables of interest
  double phisum = sum(phi);
  double mu = b0(0)/(1-phisum);
  arma::vec u  = y - X*b0;
  double stdev  = as_scalar(sqrt((trans(u)*u)/(n-1)));
  double mu_u = mean(u);
  arma::mat repvec(n,1,arma::fill::ones);
  arma::mat setmp = pow(stdev,2)*inv(trans(X)*X);
  arma::vec se = sqrt(setmp.diag());
  // ----- log-likelihood
  double logLike = sum(log(dnorm(as<NumericVector>(wrap(u)), mu_u, stdev)));
  // ----- Organize output
  double sigma = pow(stdev,2);
  arma::vec theta(2,arma::fill::zeros);
  theta(0) = mu;
  theta(1) = sigma;
  if (ar>0){
    theta = join_vert(theta, phi);
  }
  List ARmdl_out;
  ARmdl_out["y"] = y;
  ARmdl_out["X"] = X;
  ARmdl_out["x"] = x;
  ARmdl_out["n"] = n;
  ARmdl_out["k"] = 1;
  ARmdl_out["ar"] = ar;
  ARmdl_out["coef"] = b0;
  ARmdl_out["phi"] = phi;
  ARmdl_out["mu"] = mu;
  ARmdl_out["residuals"] = u;
  ARmdl_out["stdev"] = stdev;
  ARmdl_out["sigma"] = sigma;
  ARmdl_out["se"] = se;
  ARmdl_out["theta"] = theta;
  ARmdl_out["logLike"] = logLike;
  return(ARmdl_out);
}
// ==============================================================================
//' @title Fitting Vector Autoregressive Model
//' 
//' @description This function estimates a vector autoregresive series
//' 
//' @export
// [[Rcpp::export]]
List VARmdl(arma::mat Y, int ar, bool intercept = 1){
  // ----- transform data
  int N = Y.n_cols;
  List lagged_vals = ts_lagged(Y, ar);
  arma::mat y = lagged_vals["y"];
  arma::mat x = lagged_vals["X"];
  arma::mat X = lagged_vals["X"];
  arma::mat b0 (N*ar + intercept, N, arma::fill::zeros);
  arma::mat phi (N*ar , N, arma::fill::zeros);
  arma::vec inter(N, arma::fill::zeros);
  int n = y.n_rows;
  int npar = ar;
  // ----- estimate model
  if (intercept==TRUE){
    arma::mat onevec(n, 1, arma::fill::ones);
    X = join_rows(onevec,X);
    npar = npar + 1;
    b0  = solve(trans(X)*X,trans(X))*y;
    inter = trans(b0.submat(0,0,0,N-1));
    phi = b0.submat(1,0,N*ar,N-1);
  }else{
    b0  = solve(trans(X)*X,trans(X))*y;
    phi = b0;
  }
  // ----- Obtain residuals
  arma::mat yhat = X*b0;
  arma::mat resid = y - yhat;
  // ----- Obtain variance and standard errors
  arma::mat sigma = (trans(resid)*resid)/n;
  arma::vec std_error = sqrt(arma::diagvec(kron(sigma,(inv(trans(X)*X)))));
  // ----- Obtain companion form
  //arma::mat residzero(n, N*(ar-1), arma::fill::zeros);
  //arma::mat comp_resid = join_rows(resid, residzero);
  arma::vec interzero(N*(ar-1), arma::fill::zeros);
  arma::vec comp_inter = join_vert(inter,interzero);
  //arma::mat F_tmp(N, N*ar);
  //for (int xp = 0; xp<ar; xp++){
  //  F_tmp.submat(0,xp*N,N-1,xp*N+N-1) = phi.submat(xp*N,0,xp*N+N-1,N-1);  
  //}
  arma::mat F_tmp = trans(phi);
  arma::mat diagmat = arma::eye(N*(ar-1),N*(ar-1));
  arma::mat diagzero(N*(ar-1),N,arma::fill::zeros);
  arma::mat Mn = join_rows(diagmat,diagzero);
  arma::mat F = join_cols(F_tmp,Mn);
  // // ----- Obtain variables of interest
  arma::vec eigenval = abs(arma::eig_gen(F));
  arma::vec mu_tmp = inv(arma::eye(N*ar,N*ar)-F)*comp_inter;
  arma::vec mu = mu_tmp.subvec(0,N-1);
  // arma::vec mu_u = mean(resid, 0);
  arma::vec stdev = sqrt(arma::diagvec(sigma));
  // // ----- log-likelihood
  double pi = arma::datum::pi;
  arma::vec f_t(n, arma::fill::zeros);
  for (int xt = 0; xt<n; xt++){
    f_t(xt) = as_scalar((1/sqrt(det(sigma)*pow(2*pi,N)))*exp(-0.5*(resid.row(xt)*inv(sigma)*trans(resid.row(xt)))));
  }
  double logLike = sum(log(f_t));
  // ----- Organize output
  List VARmdl_out;
  VARmdl_out["y"] = y;
  VARmdl_out["X"] = X;
  VARmdl_out["x"] = x;
  VARmdl_out["n"] = n;
  VARmdl_out["k"] = 1;
  VARmdl_out["ar"] = ar;
  VARmdl_out["coef"] = b0;
  VARmdl_out["standard_errors"] = std_error;
  VARmdl_out["phi"] = phi;
  VARmdl_out["mu"] = mu;
  VARmdl_out["residuals"] = resid;
  VARmdl_out["stdev"] = stdev;
  VARmdl_out["sigma"] = sigma;
  VARmdl_out["logLike"] = logLike;
  VARmdl_out["abs_eigval"] = eigenval;
  VARmdl_out["companion_mat"] = F;
  return(VARmdl_out);
}
// ==============================================================================
//' @title Markov-switching AR Log-likelihood objective function (used to find Hessian)
//' 
//' 
//' @export
// [[Rcpp::export]]
double MSloglik_fun(arma::vec theta, List mdl, int k){
  Rcpp::Environment mstest("package:MSTest");
  Rcpp::Function transMatAR = mstest["transMatAR"];
  Rcpp::Function musigGrid = mstest["musigGrid"];
  // ============================================================================
  // ---------- Initialize parameters
  // ============================================================================
  arma::vec y = mdl["y"];
  int ar = mdl["ar"];
  bool msmu = mdl["msmu"];
  bool msvar = mdl["msvar"];
  int Tsize = y.n_elem;
  int M = pow(k, ar+1);
  // ----- Mean for each regime 
  arma::vec mu = theta.subvec(0, msmu*(k-1));
  // ----- Variance for each regime 
  arma::vec sig = theta.subvec(1+msmu*(k-1),1+msmu*(k-1)+msvar*(k-1));
  // ----- Phi vector
  arma::vec phi(std::max(ar,1), arma::fill::zeros);
  if (ar>0){
    phi =  theta.subvec(2+msmu*(k-1)+msvar*(k-1), 2+msmu*(k-1)+msvar*(k-1)+ar-1);
  }
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
  arma::mat P_AR = as<arma::mat>(transMatAR(P, k, ar));
  arma::mat pinf_AR = limP(P_AR, M);
  // ============================================================================
  // ----- Compute residuals in each regime
  // ============================================================================
  arma::mat repvec(1, M, arma::fill::ones);
  arma::mat repmu(Tsize, 1, arma::fill::ones);
  arma::mat ms_y = y*repvec;
  arma::mat eps(Tsize, M, arma::fill::zeros);
  // ---------- Compute Residuals 
  if(ar>0){
    arma::mat z = ms_y - repmu*trans(muAR.col(0)); // [y(t) - mu_s(t))]
    arma::mat x = mdl["x"];
    arma::mat xz(Tsize, M, arma::fill::zeros);
    for (int xkp = 0; xkp<M; xkp++){
      arma::mat zx_tmp = x - repmu*muAR.submat(xkp,1,xkp,ar); // [y(t-i) - mu_s(t-i))]
      xz.col(xkp) = zx_tmp*phi;
    }
    eps = z - xz;
  }else{
    eps = ms_y - repmu*trans(muAR);
  }
  // ============================================================================
  // ----- Begin Calculating likelihood and log-likelihood
  // ============================================================================
  arma::mat eta(Tsize, M, arma::fill::zeros);           // [eq. 22.4.2]
  arma::mat f_t(Tsize, 1, arma::fill::zeros);           // [eq. 22.4.8]
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
//' @title Markov-switching AR Log-likelihood objective function Minization
//' 
//' 
//' @export
// [[Rcpp::export]]
double MSloglik_fun_min(arma::vec theta, List mdl, int k){
  double logL_negative = -MSloglik_fun(theta, mdl, k);
  return(logL_negative);
}
// ==============================================================================
//' @title Markov-switching Model equality constraint
//' 
//' 
//' @export
// [[Rcpp::export]]
arma::vec MSloglik_const_eq(arma::vec theta, List mdl, int k){
  // ============================================================================
  // ---------- Initialize parameters
  // ============================================================================
  int ar = mdl["ar"];
  bool msmu = mdl["msmu"];
  bool msvar = mdl["msvar"];
  int n_p = msmu*(k-1)+msvar*(k-1)+ar+2;
  arma::mat P = reshape(theta.subvec(n_p, n_p+k*k-1), k, k);
  arma::vec constraint = trans(arma::sum(P,0))-1;
  return(constraint);
}
// ==============================================================================
//' @title Markov-switching Model Inequality Constraint
//' 
//' 
//' @export
// [[Rcpp::export]]
double MSloglik_const_ineq(arma::vec theta, List mdl, int k){
  // ============================================================================
  // ---------- Initialize parameters
  // ============================================================================
  Rcpp::Function polyroot("polyroot");  
  Rcpp::Function Mod("Mod");  
  int ar = mdl["ar"];
  bool msmu = mdl["msmu"];
  bool msvar = mdl["msvar"];
  int n_p = msmu*(k-1)+msvar*(k-1)+2;
  arma::vec phi =  theta.subvec(n_p, n_p+ar-1);
  arma::vec poly_fun(ar+1, arma::fill::ones);
  poly_fun.subvec(1,ar) = -phi;
  arma::vec roots = as<arma::vec>(Mod(wrap(as<ComplexVector>(polyroot(wrap(poly_fun))))));
  double ineq_constraint = 1-roots.min();
  return(ineq_constraint);
}
// ==============================================================================
//' @title Markov-switching loglikelihood and Hamilton smoother
//' 
//' @description This function computes loglikelihood for Markov-switching models and computes the hamilton smoother.
//' 
//' @param mdl List containing relevant parameters.
//' @param theta vector with mean and variance in each regime.
//' @param P matrix with transition probabilities.
//' @param pinf limiting probabilities of each regime.
//' @param k number of regimes. Must be greater than or equal to 2. 
//'  
//' @return List with log-likelihood (loglik), and smoothed probabilities of each regime (xi_t_T).
//' 
//' @export
// [[Rcpp::export]]
List MSloglik(arma::vec theta, List mdl, int k){
  // ============================================================================
  // ---------- Initialize parameters
  // ============================================================================
  arma::vec y = mdl["y"];
  int Tsize = y.n_elem;
  int ar = mdl["ar"];
  bool msmu = mdl["msmu"];
  bool msvar = mdl["msvar"];
  int M = pow(k, ar+1);
  List pList = paramList(theta, ar, k, msmu, msvar);
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
  // ============================================================================
  // ----- Compute residuals in each regime
  // ============================================================================
  mdl["phi"] = phi;
  arma::mat eps = calcMSResid(mdl, muAR, k);
  // ============================================================================
  // ----- Begin Calculating log-likelihood & filtered probabilities
  // ============================================================================
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
  // ============================================================================
  // ---------- Hamilton smoother
  // ============================================================================
  arma::mat xi_t_T(Tsize, k, arma::fill::zeros); // [eq. 22.4.14]
  arma::mat xi_t_T_tmp(Tsize, M, arma::fill::zeros); 
  xi_t_T_tmp.row(Tsize-1)  = xi_t_t_tmp.row(Tsize-1);
  xi_t_T.row(Tsize-1)  = xi_t_t.row(Tsize-1);
  for (int xT = Tsize-2; xT>=0; xT--){
    xi_t_T_tmp.row(xT) = xi_t_t_tmp.row(xT)%trans(trans(PAR)*trans(xi_t_T_tmp.row(xT+1)/xi_tp1_t_tmp.row(xT)));
    xi_t_T.row(xT) = xi_t_t.row(xT)%trans(trans(P)*trans(xi_t_T.row(xT+1)/xi_tp1_t.row(xT)));
  }
  // ============================================================================
  // ----- Organize output
  // ============================================================================
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
  MSloglik_output["resid"] = eps;
  MSloglik_output["muAR"] = muAR;
  MSloglik_output["sigmaAR"] = sigAR;
  MSloglik_output["mu"] = mu;
  MSloglik_output["sigma"] = sig;
  MSloglik_output["phi"] = phi;
  MSloglik_output["theta"] = theta;
  MSloglik_output["state_ind"] = state_ind;
  MSloglik_output["residuals"] = eps;
  return(MSloglik_output);
}
// ==============================================================================
//' @title EM Algo Maximization step
//' 
//' @description This function performs the maximization step of the Expectation Maximization (EM) algorithm.
//' 
//' @export
// [[Rcpp::export]]
List EMaximization(arma::vec theta, List mdl, List MSloglik_output, int k){
  Rcpp::Environment mstest("package:MSTest");
  Rcpp::Function musigGrid = mstest["musigGrid"];
  // ============================================================================
  // ---------- Initialize parameters
  // ============================================================================
  arma::vec y = mdl["y"];
  int ar = mdl["ar"];
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
  // ============================================================================
  // ----- Update transition matrix (P) and limiting probabilities (pinf)
  // ============================================================================
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
  arma::vec pinf = limP(P_new, k);
  // ============================================================================
  // ----- Update estimates for mu
  // ============================================================================
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
    mu = theta(0); // Change this to update with xi_t_T so that a good estimate is obtained even if bad initial values are given
  }
  // ============================================================================
  // ----- Update estimates for phi
  // ============================================================================
  arma::vec phi_new(std::max(1,ar), arma::fill::zeros);
  arma::mat muAR(M, ar+1,arma::fill::zeros);
  arma::vec sigAR(M, 1,arma::fill::zeros);
  List mugrid = musigGrid(mu, sig, k, ar, msmu, msvar);
  muAR = as<arma::mat>(mugrid["mu"]);
  sigAR = as<arma::vec>(mugrid["sig"]);
  if (ar>0){
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
    phi_new = inv(denom)*num;
  }
  // ============================================================================
  // ----- Update estimates for variance
  // ============================================================================
  // ----- Compute new residuals
  mdl["phi"] = phi_new;
  arma::mat resid = calcMSResid(mdl, muAR, k);
  // ----- Compute updated sigma
  arma::vec sigma(1+msvar*(k-1), arma::fill::zeros);
  arma::vec sigma_tmp(M, arma::fill::zeros);
  if (msvar == TRUE){
    for (int xk = 0 ; xk<M; xk++){
      sigma_tmp(xk) = sum(resid.col(xk)%resid.col(xk)%xi_t_T_AR.col(xk)); 
      sum_tmp(xk) = sum(xi_t_T_AR.col(xk)); 
    }
    for (int xk = 1; xk<=k;xk++){
      sigma(xk-1) = as_scalar(sum(sigma_tmp.rows(find(state_ind==xk)))/sum(sum_tmp.rows(find(state_ind==xk))));
    }
  }else{
    sigma = arma::sum(arma::sum(resid%resid%xi_t_T_AR,0),1)/Tsize;
  }
  // ============================================================================
  // ----- Organize output
  // ============================================================================
  // ----- Produce new theta vector
  arma::vec theta_new = join_vert(mu, sigma);
  if (ar>0){
    theta_new = join_vert(theta_new, phi_new);
  }
  theta_new = join_vert(theta_new, vectorise(P_new));
  // ----- Fill output List
  List Maxim_output;
  Maxim_output["theta"] = theta_new;
  Maxim_output["P"] = P_new;
  Maxim_output["pinf"] = pinf;
  Maxim_output["mu"] = mu;
  Maxim_output["sigma"] = sigma;
  Maxim_output["phi"] = phi_new;
  Maxim_output["residuals"] = resid;
  return(Maxim_output);
}
// ==============================================================================
//' @title Markov-switching loglikelihood and Hamilton smoother
//' 
//' @description This function computes loglikelihood for Markov-switching models and computes the hamilton smoother.
//' 
//' @param mdl List containing relevant parameters.
//' @param theta vector with mean and variance in each regime.
//' @param P matrix with transition probabilities.
//' @param pinf limiting probabilities of each regime.
//' @param k number of regimes. Must be greater than or equal to 2. 
//'  
//' @return List with log-likelihood (loglik), and smoothed probabilities of each regime (xi_t_T).
//' 
//' @export
// [[Rcpp::export]]
List MSVARloglik(arma::vec theta, List mdl, int k){
  // ============================================================================
  // ---------- Initialize parameters
  // ============================================================================
  arma::mat y = mdl["y"];
  int N = y.n_cols;
  int Tsize = y.n_rows;
  int ar = mdl["ar"];
  bool msmu = mdl["msmu"];
  bool msvar = mdl["msvar"];
  int M = pow(k, ar+1);
  List pList = VARparamList(theta, N, ar, k, msmu, msvar);
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
  // ============================================================================
  // ----- Compute residuals in each regime
  // ============================================================================
  mdl["phi"] = phi;
  List eps = calcMSVARResid(mdl, muAR, k);
  // ============================================================================
  // ----- Begin Calculating log-likelihood & filtered probabilities
  // ============================================================================
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
      eta(xt,xm) = as_scalar((1/(sqrt(pow(2*pi,N)*det(sigma_m))))*exp(-0.5*(eps_m.row(xt)*inv(sigma_m)*trans(eps_m.row(xt)))));
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
  // ============================================================================
  // ---------- Hamilton smoother
  // ============================================================================
  arma::mat xi_t_T(Tsize, k, arma::fill::zeros); // [eq. 22.4.14]
  arma::mat xi_t_T_tmp(Tsize, M, arma::fill::zeros); 
  xi_t_T_tmp.row(Tsize-1)  = xi_t_t_tmp.row(Tsize-1);
  xi_t_T.row(Tsize-1)  = xi_t_t.row(Tsize-1);
  for (int xT = Tsize-2; xT>=0; xT--){
    xi_t_T_tmp.row(xT) = xi_t_t_tmp.row(xT)%trans(trans(PAR)*trans(xi_t_T_tmp.row(xT+1)/xi_tp1_t_tmp.row(xT)));
    xi_t_T.row(xT) = xi_t_t.row(xT)%trans(trans(P)*trans(xi_t_T.row(xT+1)/xi_tp1_t.row(xT)));
  }
  // ============================================================================
  // ----- Organize output
  // ============================================================================
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
  MSloglik_output["resid"] = eps;
  MSloglik_output["muAR"] = muAR;
  MSloglik_output["sigmaAR"] = sigAR;
  MSloglik_output["mu"] = mu;
  MSloglik_output["sigma"] = sig;
  MSloglik_output["phi"] = phi;
  MSloglik_output["theta"] = theta;
  MSloglik_output["state_ind"] = state_ind;
  MSloglik_output["residuals"] = eps;
  return(MSloglik_output);
}
// ==============================================================================
//' @title EM Algo Maximization step
//' 
//' @description This function performs the maximization step of the Expectation Maximization (EM) algorithm.
//' 
//' @export
// [[Rcpp::export]]
List VAREMaximization(arma::vec theta, List mdl, List MSloglik_output, int k){
  Rcpp::Environment mstest("package:MSTest");
  Rcpp::Function musigVARGrid = mstest["musigVARGrid"];
  // ============================================================================
  // ---------- Initialize parameters
  // ============================================================================
  arma::mat y = mdl["y"];
  int ar = mdl["ar"];
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
  int N = y.n_cols;
  // Number of regimes consistent with AR lags
  int M = pow(k, ar+1);
  // ============================================================================
  // ----- Update transition matrix (P) and limiting probabilities (pinf)
  // ============================================================================
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
  arma::vec pinf = limP(P_new, k);
  // ============================================================================
  // ----- Update estimates for mu
  // ============================================================================
  // ----- Obtain state indicators
  arma::mat mu(1+msmu*(k-1), N, arma::fill::zeros);
  arma::mat mu_k(k, N, arma::fill::zeros);
  arma::mat mu_tmp(M, N, arma::fill::zeros);
  arma::vec sum_tmp(M, arma::fill::zeros);
  // ----- Compute updated mu
  arma::mat repN(1,N, arma::fill::ones);
  if (msmu == TRUE){
    for (int xk = 0 ; xk<M; xk++){
      mu_tmp.row(xk) = arma::sum(y%(xi_t_T_AR.col(xk)*repN),0);
      sum_tmp(xk) = sum(xi_t_T_AR.col(xk));
    }
    for (int xk = 1; xk<=k;xk++){
      mu.row(xk-1) = arma::sum(mu_tmp.rows(find(state_ind==xk)),0)/(sum(sum_tmp.rows(find(state_ind==xk)))*repN);
    }
    mu_k = mu;
  }else{
    for (int xk = 0 ; xk<M; xk++){
      mu = mu + arma::sum(y%(xi_t_T_AR.col(xk)*repN),0);
    }
    mu = mu/Tsize;
    for (int xk = 0 ; xk<k; xk++){
      mu_k.row(xk) = mu;
    }
  }
  // ============================================================================
  // ----- Update estimates for phi
  // ============================================================================
  arma::mat phi_new(N*ar, N,  arma::fill::zeros);
  List mugrid = musigVARGrid(mu_k, sig, k, ar, msmu, msvar);
  List muAR = mugrid["mu"];
  List sigAR = mugrid["sig"];
  arma::mat x = mdl["x"];
  arma::vec repmu(Tsize, arma::fill::ones);
  arma::mat denom(N*N*ar, N*N*ar, arma::fill::zeros);
  arma::mat num(N*N*ar, 1, arma::fill::zeros);
  for (int xm = 0; xm<M; xm++){
    arma::mat muAR_tmp = muAR[xm]; 
    arma::mat y_tmp = y - repmu*trans(muAR_tmp.col(0));
    arma::mat xz_tmp(Tsize, N*ar, arma::fill::zeros); 
    for (int xp = 0; xp<ar; xp++){
      xz_tmp.submat(0,N*xp,Tsize-1,N*xp+N-1) = x.submat(0,N*xp,Tsize-1,N*xp+N-1) - repmu*trans(muAR_tmp.col(xp+1));
    }
    arma::mat sigma_m = sigAR[xm];
    //denom = denom + kron(trans(xz_tmp)*diagmat(xi_t_T_AR.col(xm))*xz_tmp,inv(sigma_m));
    denom = denom + kron(inv(sigma_m),trans(xz_tmp)*diagmat(xi_t_T_AR.col(xm))*xz_tmp);
    //num = num + kron(trans(xz_tmp)*diagmat(xi_t_T_AR.col(xm)),inv(sigma_m))*vectorise(y_tmp);
    num = num + kron(inv(sigma_m),trans(xz_tmp)*diagmat(xi_t_T_AR.col(xm)))*vectorise(y_tmp);
  }
  phi_new = reshape(inv(denom)*num, N*ar, N); 
  // ============================================================================
  // ----- Update estimates for variance
  // ============================================================================
  // ----- Compute new residuals
  mdl["phi"] = phi_new;
  List resid = calcMSVARResid(mdl, muAR, k);
  // // ----- Compute updated sigma
  List sigma(k);
  for (int xk = 0; xk<k;xk++){
    arma::mat sigma_m_tmp(N,N,arma::fill::zeros);
    sigma[xk] = sigma_m_tmp;
  }
  int sigN = (N*(N+1))/2;
  arma::vec sigma_out(sigN+sigN*msvar*(k-1), arma::fill::zeros);
  arma::vec Tsize_k(k, arma::fill::zeros);
  if (msvar == TRUE){
    for (int xk = 0 ; xk<M; xk++){
      arma::mat U = resid[xk];
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
      for (int xn = 1; xn<N; xn++){
        sigma_out_tmp = join_vert(sigma_out_tmp,trans(sigma_m_tmp.submat(xn,xn,xn,N-1)));
      }
      sigma_out.subvec(xk*sigN,xk*sigN+sigN-1) = sigma_out_tmp;
    }
  }else{
    arma::mat sigma_m_tmp(N,N,arma::fill::zeros);
    for (int xk = 0 ; xk<M; xk++){
      arma::mat U = resid[xk];
      sigma_m_tmp  = sigma_m_tmp + trans(U)*diagmat(xi_t_T_AR.col(xk))*U;
    }
    sigma_m_tmp = sigma_m_tmp/Tsize;
    for (int xk = 0; xk<k;xk++){
      sigma[xk] = sigma_m_tmp;
    }
    arma::vec sigma_out_tmp = trans(sigma_m_tmp.row(0));
    for (int xn = 1; xn<N; xn++){
      sigma_out_tmp = join_vert(sigma_out_tmp, trans(sigma_m_tmp.submat(xn,xn,xn,N-1)));
    }
    sigma_out = sigma_out_tmp;
  }
  // ============================================================================
  // ----- Organize output
  // ============================================================================
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
  Maxim_output["phi"] = phi_new;
  Maxim_output["residuals"] = resid;
  return(Maxim_output);
}
// ==============================================================================
//' @title EM Algorithm 
//' 
//' @description This function performs a single iteration of the Expectation Maximization (EM) algorithm
//' 
//' @export
// [[Rcpp::export]]
List EMiter(List mdl, List EMest_output, int k){
  // ============================================================================
  // ---------- Expectation
  // ============================================================================
  arma::vec theta = EMest_output["theta"];
  //arma::mat P = EMest_output["P"];
  List MSloglik_output = MSloglik(theta, mdl, k);  
  // ============================================================================
  // ---------- Maximization
  // ============================================================================
  List Maxim_output = EMaximization(theta, mdl, MSloglik_output, k);
  // ============================================================================
  // ---------- Organize output
  // ============================================================================
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
  return(EM_output);
}
// ==============================================================================
//' @title EM Algorithm Iteration for VAR model
//' 
//' @description This function performs a single iteration of the Expectation Maximization (EM) algorithm
//' 
//' @export
// [[Rcpp::export]]
List EMiter_VAR(List mdl, List EMest_output, int k){
  // ============================================================================
  // ---------- Expectation
  // ============================================================================
  arma::vec theta = EMest_output["theta"];
  //arma::mat P = EMest_output["P"];
  List MSloglik_output = MSVARloglik(theta, mdl, k);  
  // ============================================================================
  // ---------- Maximization
  // ============================================================================
  List Maxim_output = VAREMaximization(theta, mdl, MSloglik_output, k);
  // ============================================================================
  // ---------- Organize output
  // ============================================================================
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
  return(EM_output);
}
// ==============================================================================
//' @title Estimation by EM Algorithm 
//' 
//' @description Estimated by the EM algorithm by calling the function EM() 
//' multiple times until convergence. Convergence is determined by a max number 
//' of iterations (maxit) or a convergence criterion. Specifically, the convergence 
//' criterion is the absolute difference between parameter estimates in adjacent 
//' iterations be less than 1.e-6 in absolute value. These can be changed by 
//' defining them in the optim_options List. 
//' 
//' @param z vector with observations.
//' @param k number of regimes. Must be greater than or equal to 2. 
//' @param msmu bool indicator for switch in mean (TRUE) or constant mean (FALSE). Default is TRUE.
//' @param msvar bool indicator for switch in variance (TRUE) or constant variance (FALSE). Default is TRUE.
//' @param theta_0  vector with initial values for mean and variance parameters in each state (e.g. theta_0 = c(mu_1, mu_2,..., mu_K, var_1, var_2, ..., var_K)). If NULL then initVal() is used to generate initial values. Default is NULL.
//' @param P_0 matrix with initial values for transition matrix P. Is NULL, randTransMat() is used to generate a random transition matrix. Default is NULL.
//' @param optim_options List with optimization options. When not specified (i.e. NULL), the max number of EM iterations is 200 and the convergence criterion is 1.e-6. Default is NULL.
//' 
//' @return List with model results
//' 
//' @export
// [[Rcpp::export]]
List EMest(arma::vec theta_0, List mdl, int k, List optim_options){
  // ============================================================================
  // ---------- Get optimization options
  // ============================================================================
  int maxit = optim_options["maxit"];
  double thtol = optim_options["thtol"];
  // ============================================================================
  // ---------- Begin EM Algo
  // ============================================================================
  // ----- Initial values
  List EMest_output;
  EMest_output["theta"] = theta_0;  // initial values for parameters
  EMest_output["logLike"] = 0;     // initial log-likelihood
  EMest_output = EMiter(mdl, EMest_output, k);
  int it = 1;
  double deltath = EMest_output["deltath"];
  // iterate to unitl convergence criterion is met or max number of iterations (maxit) 
  while ((it<maxit) & (deltath>thtol)){
    EMest_output = EMiter(mdl, EMest_output, k);
    deltath = EMest_output["deltath"];
    it = it + 1;
  }
  EMest_output["iterations"] = it; 
  return(EMest_output);
}
// ==============================================================================
//' @title Estimation of MSVAR by EM Algorithm 
//' 
//' 
//' @export
// [[Rcpp::export]]
List EMest_VAR(arma::vec theta_0, List mdl, int k, List optim_options){
  // ============================================================================
  // ---------- Get optimization options
  // ============================================================================
  int maxit = optim_options["maxit"];
  double thtol = optim_options["thtol"];
  // ============================================================================
  // ---------- Begin EM Algo
  // ============================================================================
  // ----- Initial values
  List EMest_output;
  EMest_output["theta"] = theta_0;  // initial values for parameters
  EMest_output["logLike"] = 0;     // initial log-likelihood
  EMest_output = EMiter_VAR(mdl, EMest_output, k);
  int it = 1;
  double deltath = EMest_output["deltath"];
  // iterate to unitl convergence criterion is met or max number of iterations (maxit) 
  while ((it<maxit) & (deltath>thtol)){
    EMest_output = EMiter_VAR(mdl, EMest_output, k);
    deltath = EMest_output["deltath"];
    it = it + 1;
  }
  EMest_output["iterations"] = it; 
  return(EMest_output);
}
// ==============================================================================
//' @title Markov-switching Autoregressive (AR) model
//' 
//' @description This function estimates a Markov-switching autoregressive (AR) model. 
//' This function uses the EM algorithm and has the option to use an optimization method 
//' as a second step to improve estimates. The function also allows the mean and/or the 
//' variance to switch but the autoregressive parameters remain constant accross regimes.
//' 
//' @param Y (Tx1) vector with observational data. Required argument.
//' @param p integer for the number of lags to use in estimation. Must be greater than or equal to 0. Default is 0.
//' @param k integer for the number of regimes to use in estimation. Must be greater than or equal to 2. Default is 2.
//' @param msmu bool indicator for switch in mean (TRUE) or constant mean (FALSE). Default is TRUE.
//' @param msvar bool indicator for switch in variance (TRUE) or constant variance (FALSE). Default is TRUE.
//' @param theta_0  vector with initial values for mean and variance parameters in each state (e.g. theta_0 = c(mu_1, mu_2,..., mu_K, var_1, var_2, ..., var_K)). If NULL then initVal() is used to generate initial values. Default is NULL.
//' @param P_0 matrix with initial values for transition matrix P. Is NULL, randTransMat() is used to generate a random transition matrix. Default is NULL.
//' @param optim_options List with optimization options. When not specified (i.e. NULL), the max number of EM iterations is 200 and the convergence criterion is 1.e-6. Default is NULL.
//' 
//' @return List with model characteristics including:
//' - y: transformed process
//' - X: matrix of lagged observations (with or without vector of 1s depending on const=1 or const=0)
//' - x: matrix of lagged observations without vector of 1s
//' - n: number of observations (T-ar)
//' - ar: number of autoregressive parameters
//' - coef: coefficient estimates. This is the same as phi is const=0
//' - phi: autoregressive coefficient estimates
//' - mu: the mean of the process
//' - residuals: vectot of residuals
//' - stdev: standard deviations
//' - se: standard errors of parameter estimtes
//' - logLike: the log-likelihood 
//' 
//' @export
// [[Rcpp::export]]
List MSARmdl(arma::vec Y, int ar = 0, int k = 2, bool msmu = 1, bool msvar = 1, int maxit = 10000, double thtol = 1.e-8, bool getHess = 0, int max_init= 100){
  Rcpp::Environment mstest("package:MSTest");
  Rcpp::Function hessian = mstest["getHess"];
  Rcpp::Function transMatAR = mstest["transMatAR"];
  // =============================================================================
  // ---------- Initialize options
  // =============================================================================
  List optim_options;
  optim_options["maxit"] = maxit; // max iterations
  optim_options["thtol"] = thtol; // Convergence criterion/tolerance
  // =============================================================================
  // ---------- Transform model
  // =============================================================================
  List mdl_out;
  if (ar>0){
    bool intercept = TRUE;
    mdl_out = ARmdl(Y, ar, intercept);
  }else{
    int n =  Y.n_elem;
    double mu = mean(Y);
    double stdev = as_scalar(sqrt(trans(Y- mu)*(Y- mu)/(n-1)));
    arma::vec X(n, arma::fill::ones);
    arma::vec phi_tmp(1, arma::fill::zeros);
    mdl_out["ar"] = 0;
    mdl_out["y"] = Y;
    mdl_out["n"] = n;
    mdl_out["x"] = R_NilValue;
    mdl_out["X"] = X;
    mdl_out["mu"] = mu;
    mdl_out["phi"] = phi_tmp;
    mdl_out["stdev"] = stdev;
  }
  mdl_out["msmu"] = msmu;
  mdl_out["msvar"] = msvar;
  // ----- Use OLS output in initial values
  arma::vec y = mdl_out["y"];
  arma::vec phi = mdl_out["phi"];
  double mu = mdl_out["mu"];
  double stdev = mdl_out["stdev"];
  int Tsize =  mdl_out["n"];
  arma::vec theta(2, arma::fill::zeros);
  theta(0) = mu;
  theta(1) = pow(stdev, 2);
  // =============================================================================
  // ---------- Estimate model 
  // =============================================================================
  List EM_output;
  int init_used = 0;
  //bool finite_check = TRUE;
  //bool initVal_check = TRUE;
  bool converge_check= TRUE;
  while (converge_check==TRUE){
    // ----- Initial values
    arma::vec theta_0 = initVals(theta, k, msmu, msvar);
    arma::mat P_0 = randTransMat(k, Tsize);
    if (ar>0){
      theta_0 = join_vert(theta_0, phi);
    }
    theta_0 = join_vert(theta_0, vectorise(P_0));
    init_used = init_used + 1;
    // ----- Estimate using EM algorithm 
    EM_output = EMest(theta_0, mdl_out, k, optim_options);
    double logLike_tmp = EM_output["logLike"];
    arma::vec theta_tmp = EM_output["theta"];
    //finite_check =  (arma::is_finite(logLike_tmp)==FALSE) and (init_used<=max_init); 
    converge_check = ((arma::is_finite(logLike_tmp)==FALSE) or (theta_tmp.is_finite()==FALSE)) and (init_used<max_init);
    //initVal_check = init_used<max_init;
    //converge_check = finite_check and initVal_check;
  }
  EM_output["init_used"] = init_used;  
  // =============================================================================
  // ---------- organize output
  // =============================================================================
  List MSARmdl_output;
  //MSARmdl_output["check1"] = finite_check;
  //MSARmdl_output["check2"] = initVal_check;
  MSARmdl_output["theta"] = EM_output["theta"];
  MSARmdl_output["mu"] = EM_output["mu"];
  MSARmdl_output["sigma"] = EM_output["sigma"];
  MSARmdl_output["P"] = EM_output["P"];
  MSARmdl_output["phi"] = EM_output["phi"];
  MSARmdl_output["pinf"] = EM_output["pinf"];
  MSARmdl_output["St"] = EM_output["St"];
  MSARmdl_output["logLike"] = EM_output["logLike"];
  MSARmdl_output["residuals"] = EM_output["residuals"];
  MSARmdl_output["eta"] = EM_output["eta"];
  MSARmdl_output["thl"] = EM_output["thl"];
  MSARmdl_output["deltath"] = EM_output["deltath"];
  MSARmdl_output["y"] = mdl_out["y"];
  MSARmdl_output["ar"] = mdl_out["ar"];
  MSARmdl_output["n"] = mdl_out["n"];
  MSARmdl_output["k"] = k;
  MSARmdl_output["x"] = mdl_out["x"];
  MSARmdl_output["X"] = mdl_out["X"];
  MSARmdl_output["msmu"] = msmu;
  MSARmdl_output["msvar"] = msvar;
  arma::vec sigma = EM_output["sigma"];
  MSARmdl_output["stdev"] = sqrt(sigma);
  MSARmdl_output["optim_options"] = optim_options;
  MSARmdl_output["iterations"] = EM_output["iterations"];
  MSARmdl_output["initVals_used"] = EM_output["init_used"];
  if (getHess==TRUE){
    arma::mat Hess = as<arma::mat>(hessian(MSARmdl_output, k));
    MSARmdl_output["Hess"] = Hess; 
  }
  return(MSARmdl_output);
}
// ==============================================================================
//' @title Markov-switching Vector Autoregressive (VAR) model
//' 
//' @export
// [[Rcpp::export]]
List MSVARmdl(arma::mat Y, int ar = 1, int k = 2, bool msmu = 1, bool msvar = 1, int maxit = 10000, double thtol = 1.e-8, bool getHess = 0, int max_init = 100){
  Rcpp::Environment mstest("package:MSTest");
  Rcpp::Function hessian = mstest["getHess"];
  Rcpp::Function transMatAR = mstest["transMatAR"];
  // =============================================================================
  // ---------- Initialize options
  // =============================================================================
  List optim_options;
  optim_options["maxit"] = maxit; // max iterations
  optim_options["thtol"] = thtol; // Convergence criterion/tolerance
  // =============================================================================
  // ---------- Transform model
  // =============================================================================
  List mdl_out;
  bool intercept = TRUE;
  mdl_out = VARmdl(Y, ar, intercept);
  mdl_out["msmu"] = msmu;
  mdl_out["msvar"] = msvar;
  // ----- Use OLS output in initial values
  arma::mat y = mdl_out["y"];
  arma::mat phi = mdl_out["phi"];
  arma::vec mu = mdl_out["mu"];
  arma::mat sigma = mdl_out["sigma"];
  int Tsize =  y.n_rows;
  // =============================================================================
  // ---------- Estimate model 
  // =============================================================================
  List EM_output;
  int init_used = 0;
  bool finite_check = TRUE;
  while (finite_check==TRUE){
    // ----- Initial values
    arma::vec theta_0 = initValsVAR(mu, sigma, k, msmu, msvar);
    arma::mat P_0 = randTransMat(k, Tsize);
    theta_0 = join_vert(theta_0, vectorise(phi));
    theta_0 = join_vert(theta_0, vectorise(P_0));
    init_used = init_used + 1;
    // ----- Estimate using EM algorithm 
    EM_output = EMest_VAR(theta_0, mdl_out, k, optim_options);
    double logLike_tmp = EM_output["logLike"];
    finite_check =  (arma::is_finite(logLike_tmp)==FALSE) and init_used<=max_init; 
  }
  EM_output["init_used"] = init_used;
  return(EM_output);
}




