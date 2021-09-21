#include <RcppArmadillo.h>
#include <RcppArmadilloExtensions/sample.h>
#include "models.h"
//[[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;
// ==============================================================================
//' @title Sum of only finite values of columns of a matrix (of vector)
//' 
//' @description 
//' 
//' 
//' @return 
//' 
//' @export
// [[Rcpp::export]]
arma::vec sumfinite(arma::mat x, int ncol = 1){
  arma::vec finite_sums(ncol, arma::fill::zeros);
  for (int xc = 0; xc<ncol; xc++){
    arma::vec vec_tmp = x.col(xc);
    arma::vec finite_vec = vec_tmp.rows(find_finite(vec_tmp));
    finite_sums(xc) = sum(finite_vec);
  }
  return(finite_sums);
}
// ==============================================================================
//' @title returns finite components of matrix
//' 
//' @description 
//' 
//' 
//' @return 
//' 
//' @export
// [[Rcpp::export]]
arma::mat finitemat(arma::mat x){
  int nrow = x.n_rows;
  int ncol = x.n_cols;
  arma::mat finite_matrix(nrow, ncol, arma::fill::zeros);
  for (int xc = 0; xc<ncol; xc++){
    arma::vec finite_vector(nrow, arma::fill::zeros);
    arma::vec vec_tmp = x.col(xc);
    finite_vector.rows(find_finite(vec_tmp)) = vec_tmp.rows(find_finite(vec_tmp));
    finite_matrix.col(xc) = finite_vector;
  }
  return(finite_matrix);
}
// ==============================================================================
//' @title musigGrid cpp version (TESTING FOR SPEED - SLOWER)
//' 
//' 
//' @export
// [[Rcpp::export]]
List musigGrid_cpp(arma::vec mu, arma::vec sig, int k, int ar){
  Rcpp::Function expGrid("expand.grid");
  Rcpp::Function asMatrix("as.matrix");
  // create grid of regimes
  arma::vec repk(k,arma::fill::ones);
  repk = cumsum(repk);
  List mu_lx(ar+1);
  List sig_lx(ar+1);
  List state_lx(ar+1);
  for (int xi = 0; xi<(ar+1);xi++){
    mu_lx[xi] = mu;
    sig_lx[xi] = sig;
    state_lx[xi] = repk;
  }
  arma::mat mu_stategrid = as<arma::mat>(asMatrix(expGrid(mu_lx)));
  arma::mat sig_stategrid = as<arma::mat>(asMatrix(expGrid(sig_lx)));
  arma::mat state_indicator = as<arma::mat>(asMatrix(expGrid(state_lx)));
  // keep only relevant parts
  arma::vec sig_stategrid_out = sig_stategrid.col(0);
  arma::vec state_indicator_out = state_indicator.col(0);
  List musig_out;
  musig_out["mu"] = mu_stategrid;
  musig_out["sig"] = sig_stategrid_out;
  musig_out["state_ind"] = state_indicator_out;
  return(musig_out);
}
// ==============================================================================
//' @title Limiting Probabilities of States
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
List paramList(arma::vec theta, int ar, int k, bool msmu, bool msvar){
  Rcpp::Environment mstest("package:MSTest");
  Rcpp::Function musigGrid = mstest["musigGrid"];
  Rcpp::Function transMatAR = mstest["transMatAR"];
  arma::vec repk(k, arma::fill::ones);
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
List VARparamList(arma::vec theta, int N, int ar, int k, bool msmu, bool msvar){
  Rcpp::Environment mstest("package:MSTest");
  Rcpp::Function musigVARGrid = mstest["musigVARGrid"];
  Rcpp::Function transMatAR = mstest["transMatAR"];
  arma::vec repk(k, arma::fill::ones);
  // ----- Mean for each regime 
  arma::mat mu_k(k, N, arma::fill::zeros);
  arma::vec mu = theta.subvec(0, N+N*msmu*(k-1)-1);
  if (msmu==TRUE){
    for (int xk = 0; xk<k; xk++){
      mu_k.row(xk) = trans(mu.subvec(xk*N,xk*N+N-1));
    }
  }else{
    for (int xk = 0; xk<k; xk++){
      mu_k.row(xk) = trans(mu);
    }
  }
  // ----- Variance for each regime 
  int sigN = (N*(N+1))/2;
  arma::vec sig = theta.subvec(N+N*msmu*(k-1),N+N*msmu*(k-1)+sigN+sigN*msvar*(k-1)-1);
  List sigma(k);
  if (msvar==TRUE){
    for (int xk = 0; xk<k; xk++){
      arma::mat sigma_tmp(N,N,arma::fill::zeros);
      arma::vec sig_tmp = sig.subvec(sigN*xk,sigN*xk+sigN-1);
      int count = 0;
      for (int xn = 0; xn<N; xn++){
        for (int xnn = xn; xnn<N; xnn++){
          sigma_tmp(xn,xnn) = sig_tmp(count);
          sigma_tmp(xnn,xn) = sig_tmp(count);
          count=count+1;
        }
      }
      sigma[xk] = sigma_tmp;
    } 
  }else{
    arma::mat sigma_tmp(N,N,arma::fill::zeros);
    int count = 0;
    for (int xn = 0; xn<N; xn++){
      for (int xnn = xn; xnn<N; xnn++){
        sigma_tmp(xn,xnn) = sig(count);
        sigma_tmp(xnn,xn) = sig(count);
        count=count+1;
      }
    }
    for (int xk = 0; xk<k; xk++){
      sigma[xk] = sigma_tmp;
    }
  }
  // ----- Phi vector
  int phiN = N+N*msmu*(k-1)+sigN+sigN*msvar*(k-1);
  arma::vec phi_tmp =  theta.subvec(phiN, phiN+ N*N*ar-1);
  arma::mat phi = reshape(phi_tmp, N*ar, N);
  // ----- Transition probabilities 
  int PN = N+N*msmu*(k-1)+sigN+sigN*msvar*(k-1)+N*N*ar;
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
  if(ar>0){
    arma::mat z = ms_y - repmu*trans(mu.col(0)); // [y(t) - mu_s(t))]
    arma::mat x = mdl["x"];
    arma::vec phi = mdl["phi"];
    arma::mat xz(Tsize, M, arma::fill::zeros);
    for (int xkp = 0; xkp<M; xkp++){
      arma::mat zx_tmp = x - repmu*mu.submat(xkp,1,xkp,ar); // [y(t-i) - mu_s(t-i))]
      xz.col(xkp) = zx_tmp*phi;
    }
    resid = z - xz;
  }else{
    resid = ms_y - repmu*trans(mu);
  }
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
//' @title Random Transition Matrix
//' 
//' @description This function creates a (kxk) random transition matrix
//' 
//' @param k number of regimes. Must be greater than or equal to 2. 
//' @param n number of random sample to use. By default it is 100 but this can be set to length of TS for example
//'  
//' 
//' @return transition matrix with randomly generated entries.
//' 
//' @export
// [[Rcpp::export]]
arma::mat randTransMat(int k, int n = 200){
  arma::vec seqk = arma::linspace(1,k,k);
  arma::vec ind = RcppArmadillo::sample(seqk,n,TRUE);
  arma::mat transmat(k,k,arma::fill::zeros);
  arma::vec tmp = diff(ind);
  arma::vec ind2 = ind.rows(1,n-1);
  arma::vec ind3 = ind.rows(0,n-2);
  for(int xi=1; xi<=k;xi++){
    double n2 = sum(ind3==xi);
    for(int xxi=1; xxi<=k; xxi++){
      transmat(xi-1,xxi-1) = sum(ind2==xxi && tmp==xxi-xi)/n2;
    }
  }
  return(trans(transmat));
}
// ==============================================================================
//' @title generate initial values for EM Algorithm 
//' 
//' 
//' @export
// [[Rcpp::export]]
arma::vec initVals(arma::vec theta, int k, bool msmu, bool msvar){
  double mu = theta(0);
  double sig = theta(1);
  arma::vec mu_0(1+msmu*(k-1), arma::fill::zeros);
  arma::vec sig_0(1+msvar*(k-1), arma::fill::zeros);
  mu_0(0) = mu;
  sig_0(0) = sig;
  if (msmu==TRUE){
    arma::vec mu_k = mu + arma::randn<arma::vec>(k-1);
    mu_0.subvec(1,k-1) = mu_k;
  }
  if (msvar==TRUE){
    arma::vec sig_k = sig + sig*arma::randu<arma::vec>(k-1);
    sig_0.subvec(1,k-1) = sig_k;
  }
  arma::vec theta_0 = join_vert(mu_0, sig_0);
  return(theta_0);
}
// ==============================================================================
//' @title generate initial values for EM Algorithm 
//' 
//' 
//' @export
// [[Rcpp::export]]
arma::vec initValsVAR(arma::vec mu, arma::mat sigma, int k, bool msmu, bool msvar){
  int N = mu.n_elem;
  arma::mat mu_k(1+msmu*(k-1), N, arma::fill::zeros);
  arma::mat mu_k_tmp(k-1, N, arma::fill::zeros);
  if (msmu==TRUE){
    for(int xk = 0; xk<(k-1);xk++){
      mu_k_tmp.row(xk) = trans(mu + arma::randn());
    }
    mu_k= join_cols(trans(mu),mu_k_tmp);
  }else{
    mu_k = mu;
  }
  arma::vec mu_out = vectorise(trans(mu_k));
  int sigN = (N*(N+1))/2;
  arma::mat sigma_k(1+msvar*(k-1), sigN, arma::fill::zeros);
  arma::vec sigma_k1 = trans(sigma.row(0));
  for (int xn = 1; xn<N; xn++){
    sigma_k1 = join_vert(sigma_k1, trans(sigma.submat(xn,xn,xn,N-1)));
  }
  if (msvar==TRUE){
    arma::mat sigma_k_tmp(k-1, sigN, arma::fill::zeros);
    for(int xk = 0; xk<(k-1);xk++){
      sigma_k_tmp.row(xk) = trans(sigma_k1*(2+xk));
    }
    sigma_k = join_cols(trans(sigma_k1), sigma_k_tmp);
    }else{
      sigma_k = trans(sigma_k1);
    }
  arma::vec sigma_out = vectorise(trans(sigma_k));
  arma::vec theta_0 = join_vert(mu_out, sigma_out);
  return(theta_0);
}
// ==============================================================================
//' @title generate initial values for EM Algorithm using Kmeans algorithm 
//' 
//' 
//' @export
// [[Rcpp::export]]
List initValsKM(arma::vec Y, int k, bool msmu, bool msvar){
  Rcpp::Function kmeans("kmeans");
  int Tsize = Y.n_elem;
  // ----- use k-mean to get initial values
  List kclust = kmeans(wrap(Y),k);
  // predefine mean and standard dev. variables
  arma::vec mu(k, arma::fill::zeros);
  arma::vec sig(k, arma::fill::zeros);
  arma::vec p(k,arma::fill::zeros);
  // use loop to fill these variables
  arma::vec state_ind = kclust["cluster"];
  arma::vec repone(Tsize,arma::fill::ones);
  for (int xi = 1; xi<=k; xi++){
    mu(xi-1) = as_scalar(mean(Y.rows(find(state_ind==xi))));
    sig(xi-1) = as_scalar(pow(stddev(Y.rows(find(state_ind==xi))),2));
    p(xi-1) = as_scalar(sum(repone.rows(find(state_ind==xi))))/Tsize;
  }
  if (msmu == FALSE){
    arma::vec repmu(k, arma::fill::ones);
    mu = repmu*mean(Y);
  }
  if (msvar == FALSE){
    arma::vec repsig(k, arma::fill::ones);
    sig = repsig*pow(stddev(Y),2);
  }
  arma::mat P(k, k, arma::fill::zeros);
  for (int xk = 0; xk<k; xk++){
    arma::vec ptmp(k,arma::fill::ones);
    ptmp = ptmp*((1-p(xk))/(k-1));
    ptmp(xk) = p(xk);
    P.col(xk) = ptmp;
  }
  List init;
  init["mu"] = mu;
  init["sig"] = sig;
  init["p"] = p;
  init["P"] = P;
  return(init);
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
List simuMSAR(List mdl_h0, Rcpp::String type = "markov", int burnin = 200){
  // ----- Obtain parameters
  int n = mdl_h0["n"];
  arma::vec phi = mdl_h0["phi"];
  arma::vec mu = mdl_h0["mu"];
  arma::vec std = mdl_h0["stdev"];
  int ar = mdl_h0["ar"];
  int k = mdl_h0["k"];
  arma::vec intercept = mu*(1-sum(phi));
  arma::mat P = mdl_h0["P"];
  // pre-fill limiting probabilities (if type = 'mixture' P and pinf will be the same)
  arma::vec pinf(k, arma::fill::zeros);
  // Simulate data
  arma::vec series(n+burnin,arma::fill::zeros);
  arma::vec state_series(n+burnin,arma::fill::zeros);
  // initialize assuming series begins in state 1 (use burnin to reduce dependence on this assumption)
  int state = 0;
  series.subvec(0,ar-1) = intercept(state) + arma::randn(ar)*std(state);
  // ----- Simulate series
  arma::vec repvec(k,arma::fill::ones);
  arma::vec state_ind = cumsum(repvec)-1;
  if (type == "markov"){
    for (int xt = ar; xt<n+burnin; xt++){
      // Get new state
      arma::vec w_temp = P.col(state);
      arma::vec state_mat = cumsum(w_temp);
      state = as_scalar(state_ind(find(arma::randu() < state_mat, 1, "first")));
      state_series(xt) = state;
      // generate new obs
      arma::vec ytmp = flipud(series.subvec((xt-ar),(xt-1)));
      series(xt) = as_scalar(intercept(state) + trans(ytmp)*phi + arma::randn()*std(state));
    }
    pinf = limP(P, k);
  }
  if (type == "mixture"){
    pinf = P;
    for (int xt = ar; xt<n+burnin; xt++){
      // Get new state
      arma::vec state_mat = cumsum(pinf);
      state = as_scalar(state_ind(find(arma::randu() < state_mat, 1, "first")));
      state_series(xt) = state; 
      // generate new obs
      arma::vec ytmp = flipud(series.subvec((xt-ar),(xt-1)));
      series(xt) = as_scalar(intercept(state) + trans(ytmp)*phi + arma::randn<arma::vec>(1)*std(state));
    }
  }
  // ----- Organize output
  arma::vec series_out = series.subvec(burnin,n+burnin-1);
  arma::vec state_series_out = state_series.subvec(burnin,n+burnin-1);
  List simu_output;
  simu_output["y"] = series_out;
  simu_output["St"] = state_series_out;
  simu_output["n"] = n;
  simu_output["phi"] = phi;
  simu_output["mu"] = mu;
  simu_output["stdev"] = std;
  simu_output["P"] = P;
  simu_output["pinf"] = pinf;
  return(simu_output);
}

// ==============================================================================
//' @title Simulate Markov-switching  Series
//' 
//' @description This function simulates a Markov-switching  series
//' 
//' @param mdl_h0 List containing series properties such as n: length, ar: number of autoregressive lags, 
//' stdev: standard deviation, mu: mean of process, phi: vector of autoregressive coefficients, P: transition matrix (if type="markov")
//' or vector of component weights (if type="mixture"), k: number of regimes.
//' @param type determines type of St is a Markov process or a random mixture. Default is "markov" 
//' @param burnin number of simulated observations to remove from beginning. By assumption, series begins in state 1 
//' so it is recommended to use burnin to avoid dependence on this assumption. Default is 200.
//' 
//' @return List with Markov-switching series and its properties
//' 
//' @export
// [[Rcpp::export]]
List simuMS(List mdl_h0, Rcpp::String type = "markov", int burnin = 200){
  // Obtain parameters
  int n = mdl_h0["n"];
  arma::vec phi = mdl_h0["phi"];
  arma::vec mu = mdl_h0["mu"];
  arma::vec std = mdl_h0["stdev"];
  int k = mdl_h0["k"];
  arma::mat P = mdl_h0["P"];
  // pre-fill limiting probabilities (if type = 'mixture' P and pinf will be the same)
  arma::vec pinf(k, arma::fill::zeros);
  // Simulate data
  arma::vec series(n+burnin,arma::fill::zeros);
  arma::vec state_series(n+burnin,arma::fill::zeros);
  // initialize assuming series begins in state 1 (use burnin to reduce dependence on this assumption)
  int state = 0;
  // Simulate series
  arma::vec repvec(k, arma::fill::ones);
  arma::vec state_ind = cumsum(repvec)-1;
  if (type == "markov"){
    for (int xt = 0; xt<n+burnin; xt++){
      // Get new state
      arma::vec w_temp = P.col(state);
      arma::vec state_mat = cumsum(w_temp);
      state = as_scalar(state_ind(find(arma::randu() < state_mat, 1, "first")));
      state_series(xt) = state;
      // generate new obs
      series(xt) = mu(state) + arma::randn()*std(state);
    }
    arma::vec pinf = limP(P, k);
  }
  if (type == "mixture"){
    pinf = P;
    for (int xt = 0; xt<n+burnin; xt++){
      // Get new state
      arma::vec state_mat = cumsum(pinf);
      state = as_scalar(state_ind(find(arma::randu() < state_mat, 1, "first")));
      state_series(xt) = state; 
      // generate new obs
      series(xt) = mu(state) + arma::randn()*std(state);
    }
  }
  // organize output
  arma::vec series_out = series.subvec(burnin,n+burnin-1);
  arma::vec state_series_out = state_series.subvec(burnin,n+burnin-1);
  List simu_output;
  simu_output["y"] = series_out;
  simu_output["St"] = state_series_out;
  simu_output["n"] = n;
  simu_output["phi"] = phi;
  simu_output["mu"] = mu;
  simu_output["stdev"] = std;
  simu_output["P"] = P;
  simu_output["pinf"] = pinf;
  return(simu_output);
}

// ==============================================================================
//' @title convert covariance matrix to correlation matrix
//' 
//' 
//' @export
// [[Rcpp::export]]
arma::mat cov2corr(arma::mat cov_mat){
  arma::mat corr_mat = inv(diagmat(sqrt(cov_mat)))*cov_mat*inv(diagmat(sqrt(cov_mat)));
  return(corr_mat);
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
  arma::mat Z(Tsize+burnin-ar,N*ar, arma::fill::zeros);
  arma::mat Y(Tsize+burnin-ar,N, arma::fill::zeros);
  arma::vec Zt = vectorise(trans(eps_corr.rows(0,ar-1)));
  for (int xt = 0; xt<(Tsize+burnin-ar);xt++){
    Z.row(xt) = trans(Zt);
    Y.row(xt) = trans(nu) + Z.row(xt)*trans(phimat) + eps_corr.row(xt+ar);
    if (ar==1){
      Zt = trans(Y.row(xt));
    }else{
      arma::vec Zt_tmp = Zt.subvec(0,N*(ar-1)-1);
      Zt = join_vert(trans(Y.row(xt)), Zt_tmp);
    }
  }
  arma::mat Z_out = Z.submat(burnin,0,Tsize+burnin-ar-1,N*ar-1);
  arma::mat Y_out = Y.submat(burnin,0,Tsize+burnin-ar-1,N-1);
  // Output
  List simuVAR_out;
  simuVAR_out["eps"] = eps;
  simuVAR_out["corr_mat"] = corr_mat;
  simuVAR_out["cholesktmat"] = C;
  simuVAR_out["eps_corr"] = eps_corr;
  simuVAR_out["Z"] = Z_out;
  simuVAR_out["Y"] = Y_out;
  simuVAR_out["F_comp"] = F;
  simuVAR_out["coef"] = join_rows(nu,phimat);
  return(simuVAR_out);
}


// ==============================================================================
//' @title Simulate Markov-Switching VAR Model
//' 
//' 
//' @export
// [[Rcpp::export]]
List simuMSVAR(List mdl_h0, Rcpp::String type = "markov", int burnin = 200){
  arma::mat mu = mdl_h0["mu"];
  List cov_matLs = mdl_h0["sigma"];
  int Tsize = mdl_h0["n"];
  arma::mat phimat = mdl_h0["phi"];
  int ar = mdl_h0["ar"];
  // Number of time series variables
  int N = mu.n_cols;
  // Number of regimes
  int k = mdl_h0["k"];
  //Transition matrix
  arma::mat P = mdl_h0["P"];
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
    arma::mat U1(Tsize+burnin+ar, N, arma::fill::randu);
    arma::mat U2(Tsize+burnin+ar, N, arma::fill::randu);
    arma::mat cov_mat_k = cov_matLs[xk];
    arma::mat eps_k = trans(arma::diagmat(sqrt(cov_mat_k))*trans(sqrt(-2*log(U1))%cos(2*pi*U2)));
    // add correlations
    arma::mat corr_mat_k = cov2corr(cov_mat_k);
    corr_mat[xk] = corr_mat_k;
    arma::mat C = chol(corr_mat_k, "lower");
    arma::mat eps_corr = trans(C*trans(eps_k));
    epsLs[xk] = eps_corr;
    // get constant vec 
    arma::vec mu_tmp = vectorise(trans(repmu*mu.row(xk)));
    arma::vec nu_k_tmp = (arma::eye(N*ar,N*ar) - F)*mu_tmp;
    arma::vec nu_k = nu_k_tmp.subvec(0,N-1);
    nu.row(xk) = trans(nu_k);
  }
  // Simulate process 
  int state = 0;
  arma::vec repvec(k, arma::fill::ones);
  arma::vec state_ind = cumsum(repvec)-1;
  arma::vec state_series(Tsize+burnin, arma::fill::zeros);
  arma::mat Z(Tsize+burnin, N*ar, arma::fill::zeros);
  arma::mat Y(Tsize+burnin, N, arma::fill::zeros);
  arma::mat eps(Tsize+burnin, N, arma::fill::zeros);
  arma::mat eps_corr_k = epsLs[state];
  arma::vec Zt = vectorise(trans(eps_corr_k.rows(0,ar-1)));
  if (type == "markov"){
    for (int xt = 0; xt<(Tsize+burnin);xt++){
      // Get new state
      arma::vec w_temp = P.col(state);
      arma::vec state_mat = cumsum(w_temp);
      state = as_scalar(state_ind(find(arma::randu() < state_mat, 1, "first")));
      state_series(xt) = state;
      arma::mat eps_corr_k = epsLs[state];
      arma::mat nu_k = nu.row(state);
      Z.row(xt) = trans(Zt);
      Y.row(xt) = nu_k + Z.row(xt)*trans(phimat) + eps_corr_k.row(xt+ar);
      eps.row(xt) = eps_corr_k.row(xt+ar);
      if (ar==1){
        Zt = trans(Y.row(xt));
      }else{
        arma::vec Zt_tmp = Zt.subvec(0,N*(ar-1)-1);
        Zt = join_vert(trans(Y.row(xt)), Zt_tmp);
      }
    }
  }
  arma::mat Z_out = Z.submat(burnin,0,Tsize+burnin-1,N*ar-1);
  arma::mat Y_out = Y.submat(burnin,0,Tsize+burnin-1,N-1);
  arma::mat eps_out = eps.submat(burnin,0,Tsize+burnin-1,N-1);
  arma::vec state_out = state_series.subvec(burnin,Tsize+burnin-1);
  // Output
  List simuVAR_out;
  simuVAR_out["corr_mat"] = corr_mat;
  simuVAR_out["eps"] = eps_out;
  simuVAR_out["epsList"] = epsLs;
  simuVAR_out["Z"] = Z_out;
  simuVAR_out["Y"] = Y_out;
  simuVAR_out["F_comp"] = F;
  simuVAR_out["coef"] = join_rows(trans(nu),phimat);
  simuVAR_out["St"] = state_out;
  return(simuVAR_out);
}
