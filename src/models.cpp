#include <RcppArmadillo.h>
#include "methods.h"
//[[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;


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
double AR_loglik_fun(arma::vec theta, List mdl){
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
//' @title Autoregressive Model
//' 
//' @description This function estimates an autoregresive model
//' 
//' @param Y vector of observations
//' @param ar number of lags
//' @param intercept boolean indicator determining whether or not to include constant. Default is TRUE.
//' @param getSE boolean indicator determining whether or not to estimate standard errors
//' 
//' @return List with model attributes which include:
//' \itemize{
//'   \item{y - }{a ((T-ar) x 1) vector of observations}
//'   \item{X - }{matrix of lagged observations (with or without vector of 1s depending on const=1 or const=0)}
//'   \item{x - }{matrix of lagged observations without vector of 1s}
//'   \item{n - }{number of observations (T-ar)}
//'   \item{ar - }{number of autoregressive parameters}
//'   \item{coef - }{coefficient estimates. This is the same as phi is const=0}
//'   \item{phi - }{autoregressive coefficient estimates}
//'   \item{mu - }{mean of the process}
//'   \item{resid - }{vectot of residuals}
//'   \item{stdev - }{standard deviations}
//'   \item{se - }{standard errors of parameter estimtes}
//'   \item{logLike - }{log-likelihood}
//' }
//' 
//' @export
// [[Rcpp::export]]
List ARmdl_cpp(arma::vec Y, int ar, bool intercept = 1, bool getSE = 0){
  Rcpp::Environment mstest("package:MSTest");
  Rcpp::Function hessian = mstest["getHess"];
  Rcpp::Environment lmf("package:lmf");
  Rcpp::Function nearPD = lmf["nearPD"];
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
    b0  = solve(trans(X)*X,trans(X))*y;
    phi = b0.subvec(1,ar);
  }else{
    b0  = solve(trans(X)*X,trans(X))*y;
    phi = b0;
  }
  // ----- Obtain variables of interest
  double phisum = sum(phi);
  double mu = b0(0)/(1-phisum);
  arma::vec u  = y - X*b0;
  double stdev  = as_scalar(sqrt((trans(u)*u)/(n-1)));
  // ----- Organize output
  double sigma = pow(stdev,2);
  arma::vec theta(2,arma::fill::zeros);
  theta(0) = mu;
  theta(1) = sigma;
  theta = join_vert(theta, phi);
  arma::vec theta_mu_ind = join_vert(arma::vec(1,arma::fill::ones),arma::vec(1+phi.n_elem));
  arma::vec theta_sig_ind = join_vert(join_vert(arma::vec(1,arma::fill::zeros),arma::vec(1,arma::fill::ones)),arma::vec(phi.n_elem,arma::fill::zeros));
  arma::vec theta_phi_ind = join_vert(arma::vec(2,arma::fill::zeros),arma::vec(phi.n_elem,arma::fill::ones));
  List ARmdl_out;
  ARmdl_out["y"] = y;
  ARmdl_out["X"] = X;
  ARmdl_out["x"] = x;
  ARmdl_out["n"] = n;
  ARmdl_out["q"] = 1;
  ARmdl_out["k"] = 1;
  ARmdl_out["p"] = ar;
  ARmdl_out["coef"] = b0;
  ARmdl_out["phi"] = phi;
  ARmdl_out["mu"] = mu;
  ARmdl_out["resid"] = u;
  ARmdl_out["stdev"] = stdev;
  ARmdl_out["sigma"] = sigma;
  ARmdl_out["theta"] = theta;
  ARmdl_out["theta_mu_ind"] = theta_mu_ind;
  ARmdl_out["theta_sig_ind"] = theta_sig_ind;
  ARmdl_out["theta_phi_ind"] = theta_phi_ind;
  double logLike = AR_loglik_fun(theta,ARmdl_out);
  ARmdl_out["logLike"] = logLike;
  if (getSE==TRUE){
    arma::mat Hess = as<arma::mat>(hessian(ARmdl_out, 1));
    arma::mat info_mat = inv(-Hess);
    bool nearPD_used = FALSE;
    if (any(diagvec(info_mat)<0)){
      info_mat = as<arma::mat>(nearPD(info_mat));
      nearPD_used = TRUE;
    }
    arma::vec theta_stderr = sqrt(diagvec(inv(-Hess)));
    ARmdl_out["Hess"] = Hess; 
    ARmdl_out["theta_stderr"] = theta_stderr;
    ARmdl_out["info_mat"] = info_mat;
    ARmdl_out["nearPD_used"] = nearPD_used;
  }
  return(ARmdl_out);
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
double VAR_loglik_fun(arma::vec theta, List mdl){
  // ---------- model parameters
  arma::mat y = mdl["y"];
  arma::mat x = mdl["x"];
  int ar = mdl["ar"];
  int Tsize = y.n_rows;
  int q = y.n_cols;
  // ---------- pre-define variables
  int sigN = (q*(q+1))/2;
  arma::vec mu = theta.subvec(0, q-1);
  arma::vec sig = theta.subvec(q,q+sigN-1);
  arma::mat sigma = covar_unvech(sig, q);
  int phiN = q+sigN;
  arma::vec phi_tmp = theta.subvec(phiN, phiN+q*q*ar-1);
  arma::mat phi = reshape(phi_tmp, q*ar, q);
  arma::mat repmu(Tsize, 1, arma::fill::ones);
  arma::mat z = y-repmu*trans(mu);
  arma::mat xz_tmp(Tsize, q*ar, arma::fill::zeros); 
  // ---------- compute residual
  for (int xp = 0; xp<ar; xp++){
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
//' @title Vector autoregressive model
//' 
//' @description This function estimates a vector autoregresive model
//' 
//' @param Y is a (T x q) matrix of observations
//' @param ar number of lags
//' @param intercept boolean indicator determining whether or not to include constant. Default is TRUE.
//' @param getSE boolean indicator determining whether or not to estimate standard errors
//' 
//' @return List with model attributes
//' 
//' @export
// [[Rcpp::export]]
List VARmdl(arma::mat Y, int ar, bool intercept = 1, bool getSE = 0){
  Rcpp::Environment mstest("package:MSTest");
  Rcpp::Function hessian = mstest["getHess"];
  Rcpp::Environment lmf("package:lmf");
  Rcpp::Function nearPD = lmf["nearPD"];
  // ----- transform data
  int q = Y.n_cols;
  List lagged_vals = ts_lagged(Y, ar);
  arma::mat y = lagged_vals["y"];
  arma::mat x = lagged_vals["X"];
  arma::mat X = lagged_vals["X"];
  arma::mat b0 (q*ar + intercept, q, arma::fill::zeros);
  arma::mat phi (q*ar , q, arma::fill::zeros);
  arma::vec inter(q, arma::fill::zeros);
  int n = y.n_rows;
  int npar = ar;
  // ----- estimate model
  if (intercept==TRUE){
    arma::mat onevec(n, 1, arma::fill::ones);
    X = join_rows(onevec,X);
    npar = npar + 1;
    b0  = solve(trans(X)*X,trans(X))*y;
    inter = trans(b0.submat(0,0,0,q-1));
    phi = b0.submat(1,0,q*ar,q-1);
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
  arma::vec interzero(q*(ar-1), arma::fill::zeros);
  arma::vec comp_inter = join_vert(inter,interzero);
  arma::mat F_tmp = trans(phi);
  arma::mat diagmat = arma::eye(q*(ar-1),q*(ar-1));
  arma::mat diagzero(q*(ar-1), q, arma::fill::zeros);
  arma::mat Mn = join_rows(diagmat,diagzero);
  arma::mat F = join_cols(F_tmp,Mn);
  // // ----- Obtain variables of interest
  arma::vec eigenval = abs(arma::eig_gen(F));
  arma::vec mu_tmp = inv(arma::eye(q*ar, q*ar)-F)*comp_inter;
  arma::mat mu = trans(mu_tmp.subvec(0, q-1));
  // arma::vec mu_u = mean(resid, 0);
  arma::vec stdev = sqrt(arma::diagvec(sigma));
  // // ----- log-likelihood
  double pi = arma::datum::pi;
  arma::vec f_t(n, arma::fill::zeros);
  for (int xt = 0; xt<n; xt++){
    f_t(xt) = as_scalar((1/sqrt(det(sigma)*pow(2*pi,q)))*exp(-0.5*(resid.row(xt)*inv(sigma)*trans(resid.row(xt)))));
  }
  double logLike = sum(log(f_t));
  // ----- Organize output
  arma::vec sigma_vec = covar_vech(sigma);
  arma::vec phi_vec = vectorise(phi);
  arma::vec theta = join_vert(join_vert(trans(mu),sigma_vec),phi_vec);
  arma::vec theta_mu_ind = join_vert(join_vert(arma::vec(mu.n_elem, arma::fill::ones),arma::vec(sigma_vec.n_elem,arma::fill::zeros)),arma::vec(phi_vec.n_elem,arma::fill::zeros));
  arma::vec theta_sig_ind = join_vert(join_vert(arma::vec(mu.n_elem, arma::fill::zeros), arma::vec(sigma_vec.n_elem,arma::fill::ones)),arma::vec(phi_vec.n_elem,arma::fill::zeros));
  arma::vec theta_var_ind = join_vert(join_vert(arma::vec(mu.n_elem, arma::fill::zeros), covar_vech(arma::eye(q,q))),arma::vec(phi_vec.n_elem,arma::fill::zeros));
  arma::vec theta_phi_ind = join_vert(arma::vec(mu.n_elem + sigma_vec.n_elem,arma::fill::zeros),arma::vec(phi_vec.n_elem,arma::fill::ones));
  List VARmdl_out;
  VARmdl_out["y"] = y;
  VARmdl_out["X"] = X;
  VARmdl_out["x"] = x;
  VARmdl_out["n"] = n;
  VARmdl_out["q"] = q;
  VARmdl_out["k"] = 1;
  VARmdl_out["ar"] = ar;
  VARmdl_out["coef"] = b0;
  VARmdl_out["phi"] = phi; // should transpose this before outputting (what functions would be affected by this?)
  VARmdl_out["mu"] = mu;
  VARmdl_out["residuals"] = resid;
  VARmdl_out["stdev"] = stdev;
  VARmdl_out["sigma"] = sigma;
  VARmdl_out["theta"] = theta;
  VARmdl_out["theta_mu_ind"] = theta_mu_ind;
  VARmdl_out["theta_sig_ind"] = theta_sig_ind;
  VARmdl_out["theta_var_ind"] = theta_var_ind;
  VARmdl_out["theta_phi_ind"] = theta_phi_ind;
  VARmdl_out["abs_eigval"] = eigenval;
  VARmdl_out["companion_mat"] = F;
  double logLike2 = VAR_loglik_fun(theta, VARmdl_out);
  VARmdl_out["logLike"] = logLike;
  VARmdl_out["logLike2"] = logLike2;
  if (getSE==TRUE){
    arma::mat Hess = as<arma::mat>(hessian(VARmdl_out, 1));
    arma::mat info_mat = inv(-Hess);
    bool nearPD_used = FALSE;
    if (any(diagvec(info_mat)<0)){
      info_mat = as<arma::mat>(nearPD(info_mat));
      nearPD_used = TRUE;
    }
    arma::vec theta_stderr = sqrt(diagvec(inv(-Hess)));
    VARmdl_out["Hess"] = Hess; 
    VARmdl_out["theta_stderr"] = theta_stderr;
    VARmdl_out["info_mat"] = info_mat;
    VARmdl_out["nearPD_used"] = nearPD_used;
  }
  VARmdl_out["standard_errors"] = std_error;
  return(VARmdl_out);
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
double MSloglik_fun(arma::vec theta, List mdl, int k){
  Rcpp::Environment mstest("package:MSTest");
  Rcpp::Function arP = mstest["arP"];
  Rcpp::Function arGrid = mstest["arGrid"];
  // ---------- Initialize parameters
  arma::vec y = mdl["y"];
  arma::mat x = mdl["x"];
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
  arma::vec phi =  theta.subvec(2+msmu*(k-1)+msvar*(k-1), 2+msmu*(k-1)+msvar*(k-1)+ar-1);
  // ----- Transition probabilities 
  arma::mat P = reshape(theta.subvec(2+msmu*(k-1)+msvar*(k-1) + ar, 2+msmu*(k-1)+msvar*(k-1) + ar + k*k - 1),k, k);
  // Regime limiting probabilities
  arma::vec pinf = limP(P);
  // ----- Obtain AR consistent grid of Mu, Sigma and State indicators
  List musig_out = arGrid(mu, sig, k, ar, msmu, msvar);
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
double MSloglik_fun_min(arma::vec theta, List mdl, int k){
  double logL_negative = -MSloglik_fun(theta, mdl, k);
  return(logL_negative);
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
List MSloglik(arma::vec theta, List mdl, int k){
  // ---------- Initialize parameters
  arma::vec y = mdl["y"];
  int Tsize = y.n_elem;
  int ar = mdl["ar"];
  bool msmu = mdl["msmu"];
  bool msvar = mdl["msvar"];
  int M = pow(k, ar+1);
  List pList = paramListMS(theta, ar, k, msmu, msvar);
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
  // ----- Compute residuals in each regime
  List mdl_tmp  = clone(mdl);
  mdl_tmp["phi"] = phi;
  arma::mat eps = calcMSResid(mdl_tmp, muAR, k);
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
double MSVARloglik_fun(arma::vec theta, List mdl, int k){
  Rcpp::Environment mstest("package:MSTest");
  Rcpp::Function arP = mstest["arP"];
  Rcpp::Function varGrid = mstest["varGrid"];
  // ---------- Initialize parameters
  arma::mat y = mdl["y"];
  arma::mat x = mdl["x"];
  int q = y.n_cols;
  int Tsize = y.n_rows;
  int ar = mdl["ar"];
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
  List musig_out = varGrid(mu_k, sigma, k, ar, msmu, msvar);
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
      eta(xt,xm) = as_scalar((1/(sqrt(pow(2*pi,q)*det(sigma_m))))*
        exp(-0.5*(eps_m.row(xt)*inv(sigma_m)*trans(eps_m.row(xt)))));
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
double MSVARloglik_fun_min(arma::vec theta, List mdl, int k){
  double logL_negative = -MSVARloglik_fun(theta, mdl, k);
  return(logL_negative);
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
List MSVARloglik(arma::vec theta, List mdl, int k){
  // ---------- Initialize parameters
  arma::mat y = mdl["y"];
  int N = y.n_cols;
  int Tsize = y.n_rows;
  int ar = mdl["ar"];
  bool msmu = mdl["msmu"];
  bool msvar = mdl["msvar"];
  int M = pow(k, ar+1);
  List pList = paramListMSVAR(theta, N, ar, k, msmu, msvar);
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
  List eps = calcMSVARResid(mdl_tmp, muAR, k);
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
  // ---------- Hamilton smoother
  arma::mat xi_t_T(Tsize, k, arma::fill::zeros); // [eq. 22.4.14]
  arma::mat xi_t_T_tmp(Tsize, M, arma::fill::zeros); 
  xi_t_T_tmp.row(Tsize-1)  = xi_t_t_tmp.row(Tsize-1);
  xi_t_T.row(Tsize-1)  = xi_t_t.row(Tsize-1);
  for (int xT = Tsize-2; xT>=0; xT--){
    xi_t_T_tmp.row(xT) = xi_t_t_tmp.row(xT)%trans(trans(PAR)*trans(xi_t_T_tmp.row(xT+1)/xi_tp1_t_tmp.row(xT)));
    xi_t_T.row(xT) = xi_t_t.row(xT)%trans(trans(P)*trans(xi_t_T.row(xT+1)/xi_tp1_t.row(xT)));
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
List MS_EMaximization(arma::vec theta, List mdl, List MSloglik_output, int k){
  Rcpp::Environment mstest("package:MSTest");
  Rcpp::Function arGrid = mstest["arGrid"];
  // ---------- Initialize parameters
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
    mu = theta(0); // Change this to update with xi_t_T so that a good estimate is obtained even if bad initial values are given
  }
  // ----- Update estimates for phi
  arma::mat muAR(M, ar+1,arma::fill::zeros);
  arma::vec sigAR(M, 1,arma::fill::zeros);
  List mugrid = arGrid(mu, sig, k, ar, msmu, msvar);
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
  arma::mat resid = calcMSResid(mdl_tmp, muAR, k);
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
  Maxim_output["residuals"] = resid;
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
List MSVAR_EMaximization(arma::vec theta, List mdl, List MSloglik_output, int k){
  Rcpp::Environment mstest("package:MSTest");
  Rcpp::Function varGrid = mstest["varGrid"];
  // ---------- Initialize parameters
  arma::mat y = mdl["y"];
  int ar = mdl["ar"];
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
  int N = y.n_cols;
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
  // --------------- Update estimates for phi
  List mugrid = varGrid(mu_k, sig, k, ar, msmu, msvar);
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
  arma::mat phi_new = reshape(inv(denom)*num, N*ar, N); 
  
  // Companion Mat
  arma::mat F_tmp = trans(phi_new);
  arma::mat diag_mat = arma::eye(q*(ar-1),q*(ar-1));
  arma::mat diag_zero(q*(ar-1), q, arma::fill::zeros);
  arma::mat Mn = join_rows(diag_mat,diag_zero);
  arma::mat F = join_cols(F_tmp,Mn); 
  // --------------- Update estimates for variance
  // ----- Compute new residuals
  List mdl_tmp  = clone(mdl);
  mdl_tmp["phi"] = phi_new;
  List resid = calcMSVARResid(mdl_tmp, muAR, k);
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
  Maxim_output["phi"] = phi_new;
  Maxim_output["companion_mat"] = F; 
  Maxim_output["residuals"] = resid;
  return(Maxim_output);
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
List MS_EMiter(List mdl, List EMest_output, int k){
  // ---------- Expectation
  arma::vec theta = EMest_output["theta"];
  List MSloglik_output = MSloglik(theta, mdl, k);
  // ---------- Maximization
  List Maxim_output = MS_EMaximization(theta, mdl, MSloglik_output, k);
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
List MSVAR_EMiter(List mdl, List EMest_output, int k){
  // ---------- Expectation
  arma::vec theta = EMest_output["theta"];
  //arma::mat P = EMest_output["P"];
  List MSVARloglik_output = MSVARloglik(theta, mdl, k);
  // ---------- Maximization
  List Maxim_output = MSVAR_EMaximization(theta, mdl, MSVARloglik_output, k);
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
  return(EM_output);
}


// ==============================================================================
//' @title Estimation of Markov-switching autoregressive model by EM Algorithm 
//' 
//' @description Estimate Markov-switching autoregressive model by EM algorithm. This function is used by MSmdl_EM() which organizes the output and takes raw data as input.
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
List MS_EMest(arma::vec theta_0, List mdl, int k, List optim_options){
  // ---------- Get optimization options
  int maxit = optim_options["maxit"];
  double thtol = optim_options["thtol"];
  // ---------- Begin EM Algo
  // ----- Initial values
  List EMest_output;
  EMest_output["theta"] = theta_0;  // initial values for parameters
  EMest_output["logLike"] = 0;     // initial log-likelihood
  EMest_output = MS_EMiter(mdl, EMest_output, k);
  int it = 1;
  double deltath = EMest_output["deltath"];
  // iterate to unitl convergence criterion is met or max number of iterations (maxit) 
  while ((it<maxit) & (deltath>thtol)){
    EMest_output = MS_EMiter(mdl, EMest_output, k);
    deltath = EMest_output["deltath"];
    it = it + 1;
  }
  EMest_output["iterations"] = it; 
  return(EMest_output);
}


// ==============================================================================
//' @title Estimation of Markov-switching vector autoregressive model by EM Algorithm 
//' 
//' @description Estimate Markov-switching vector autoregressive model by EM algorithm. This function is used by MSmdl_EM() which organizes the output and takes raw data as input.
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
List MSVAR_EMest(arma::vec theta_0, List mdl, int k, List optim_options){
  // ---------- Get optimization options
  int maxit = optim_options["maxit"];
  double thtol = optim_options["thtol"];
  // ---------- Begin EM Algo
  // ----- Initial values
  List EMest_output;
  EMest_output["theta"] = theta_0;  // initial values for parameters
  EMest_output["logLike"] = 0;     // initial log-likelihood
  EMest_output = MSVAR_EMiter(mdl, EMest_output, k);
  int it = 1;
  double deltath = EMest_output["deltath"];
  // iterate to unitl convergence criterion is met or max number of iterations (maxit) 
  while ((it<maxit) & (deltath>thtol)){
    EMest_output = MSVAR_EMiter(mdl, EMest_output, k);
    deltath = EMest_output["deltath"];
    it = it + 1;
  }
  EMest_output["iterations"] = it; 
  return(EMest_output);
}

