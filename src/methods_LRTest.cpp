#include <RcppArmadillo.h>
#include "models.h"
#include "methods.h"
//[[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;

// ==============================================================================f
//' @title Likelihood Ratio Test Statistic Sample Distribution
//' 
//' 
//' @export
// [[Rcpp::export]]
arma::vec LR_samp_dist(List mdl_h0, int k1, bool msmu, bool msvar, int N, int maxit, double thtol, int burnin, int max_init, int dist_converge_iter){
  int k0 = mdl_h0["k"];
  arma::vec LRT_N(N,arma::fill::zeros);
  int ar = mdl_h0["ar"];
  bool getHess = FALSE;
  if (k0 == 1){
    bool inter = TRUE;
    for (int xn = 0; xn<N; xn++){
      double LRT_i;
      bool LRT_finite = FALSE;
      bool LRT_converge = FALSE;
      int converge_iter = 0;
      while ((LRT_converge==FALSE) and (converge_iter<dist_converge_iter)){
        while (LRT_finite==FALSE){
          List y0_out = simuAR(mdl_h0, burnin);
          arma::vec y0 = y0_out["y"];
          List mdl_h0_tmp = ARmdl(y0, ar, inter);
          List mdl_h1_tmp = MSARmdl(y0, ar, k1, msmu, msvar, maxit, thtol, getHess, max_init);
          // test stat
          double l_0 = mdl_h0_tmp["logLike"];
          double l_1 = mdl_h1_tmp["logLike"];
          int mdl_h1_iter = mdl_h1_tmp["iterations"];
          LRT_i = -2*(l_0 - l_1);
          // LRT_finite = arma::is_finite(LRT_i);
          LRT_finite = ((arma::is_finite(LRT_i)) and LRT_i>=0);
          LRT_converge = (mdl_h1_iter<maxit);
        }
        converge_iter +=1;
      }
      LRT_N(xn) = LRT_i;
      if (converge_iter==dist_converge_iter){
        warning("Warning: Some simulations had models that did not converge. Try using higher 'maxit'.");
      }
    } 
  }else if (k0>1){
    if (msmu == FALSE){
      arma::vec muk(k0, arma::fill::ones);
      double mu_h0 = mdl_h0["mu"];
      mdl_h0["mu"] = muk*mu_h0;
    }
    if (msvar == FALSE){
      arma::vec stdevk(k0, arma::fill::ones);
      double stdev_h0 = mdl_h0["stdev"];
      mdl_h0["stdev"] = stdevk*stdev_h0;
    }
    for (int xn = 0; xn<N; xn++){
      double LRT_i;
      bool LRT_finite = FALSE;
      bool LRT_converge = FALSE;
      int converge_iter = 0;
      while ((LRT_converge==FALSE) and (converge_iter<dist_converge_iter)){
        while (LRT_finite==FALSE){
          List y0_out = simuMSAR(mdl_h0, "markov", burnin);
          arma::vec y0 = y0_out["y"];
          List mdl_h0_tmp = MSARmdl(y0, ar, k0, msmu, msvar, maxit, thtol, getHess, max_init);
          List mdl_h1_tmp = MSARmdl(y0, ar, k1, msmu, msvar, maxit, thtol, getHess, max_init);
          // test stat
          double l_0 = mdl_h0_tmp["logLike"];
          double l_1 = mdl_h1_tmp["logLike"];
          int mdl_h0_iter = mdl_h0_tmp["iterations"];
          int mdl_h1_iter = mdl_h1_tmp["iterations"];
          LRT_i = -2*(l_0 - l_1);
          // LRT_finite = arma::is_finite(LRT_i);
          LRT_finite = ((arma::is_finite(LRT_i)) and LRT_i>=0);
          LRT_converge = ((mdl_h0_iter<maxit) and (mdl_h1_iter<maxit));
        }
        converge_iter +=1;
      }
      LRT_N(xn) = LRT_i;
      if (converge_iter==dist_converge_iter){
        warning("Warning: Some simulations had models that did not converge. Try using higher 'maxit'.");
      }
    } 
  }
  return(LRT_N);
}

// ==============================================================================
//' @title Monte Carlo Likelihood Ratio Test P-value Function 
//' 
//' 
//' @export
// [[Rcpp::export]]
double MMCLRpval_fun(arma::vec theta, List mdl_h0, List mdl_h1, bool msmu, bool msvar, int ar, int N, int maxit, double thtol, int burnin, bool stationary_ind, double lambda, int max_init, int dist_converge_iter){
  // initialize variables 
  double pval;
  double logL0;
  double logL1;
  arma::mat P_h0;
  arma::mat P_h1;
  // define required variables from inputs
  bool non_stationary_const = FALSE;
  bool P_h0_colsum_const = FALSE;
  bool P_h1_colsum_const = FALSE;
  int k0 = mdl_h0["k"];
  int k1 = mdl_h1["k"];
  arma::vec theta_h0_tmp = mdl_h0["theta"];
  arma::vec theta_h1_tmp = mdl_h1["theta"];
  int theta_h0_length = theta_h0_tmp.n_elem;
  int theta_h1_length = theta_h1_tmp.n_elem;
  arma::vec theta_h0 = theta.subvec(0,theta_h0_length-1);
  arma::vec theta_h1 = theta.subvec(theta_h0_length,theta_h0_length+theta_h1_length-1);
  // ----- Stationary constraint (i.e., only consider theta that result in stationary process) 
  if ((stationary_ind==TRUE) and (ar>0)){
    Rcpp::Function polyroot("polyroot");  
    Rcpp::Function Mod("Mod");  
    arma::vec poly_fun_h0(ar+1, arma::fill::ones);
    arma::vec poly_fun_h1(ar+1, arma::fill::ones);
    poly_fun_h0.subvec(1,ar) = -theta_h0.subvec(2+msmu*(k0-1)+msvar*(k0-1), 2+msmu*(k0-1)+msvar*(k0-1)+ar-1);
    poly_fun_h1.subvec(1,ar) = -theta_h1.subvec(2+msmu*(k1-1)+msvar*(k1-1), 2+msmu*(k1-1)+msvar*(k1-1)+ar-1);
    arma::vec roots_h0 = as<arma::vec>(Mod(wrap(as<ComplexVector>(polyroot(wrap(poly_fun_h0))))));
    arma::vec roots_h1 = as<arma::vec>(Mod(wrap(as<ComplexVector>(polyroot(wrap(poly_fun_h1))))));
    non_stationary_const = ((roots_h0.min()<=1) or (roots_h1.min()<=1));
  }
  // ----- Checking that P matrix columns sum to 1
  if (k0>1){
    P_h0 = reshape(theta_h0.subvec(theta_h0_length-k0*k0, theta_h0_length-1), k0, k0);
    P_h0_colsum_const = any(abs(arma::sum(P_h0,0)-1)>thtol); 
  }
  P_h1 = reshape(theta_h1.subvec(theta_h1_length-k1*k1, theta_h1_length-1), k1, k1); 
  P_h1_colsum_const = any(abs(arma::sum(P_h1,0)-1)>thtol);
  // ----- Compute pval 
  if ((P_h0_colsum_const==TRUE) or (P_h1_colsum_const==TRUE) or (non_stationary_const==TRUE)){
    // If either transition matrix columns do not sum to 1 OR (stationary_ind == TRUE AND non_stationary_const == TRUE), pval is a positive constant
    pval = lambda*(P_h0_colsum_const + P_h1_colsum_const + non_stationary_const);
  }else{
    // If transition matrices' columns sum to 1 AND (stationary_ind==FALSE OR non_stationary_const == FALSE), pval is computed. 
    List mdl_h0_tmp = mdl_h0;
    mdl_h0_tmp["theta"] = theta_h0;
    mdl_h0_tmp["mu"] = theta_h0.subvec(0, msmu*(k0-1));
    mdl_h0_tmp["stdev"] = sqrt(theta_h0.subvec(1+msmu*(k0-1),1+msmu*(k0-1)+msvar*(k0-1)));
    if (ar>0){
      mdl_h0_tmp["phi"] = theta_h0.subvec(2+msmu*(k0-1)+msvar*(k0-1), 2+msmu*(k0-1)+msvar*(k0-1)+ar-1);
    }
    // compute test stat
    if (k0==1){
      logL0 = loglik_fun(theta_h0, mdl_h0);  
    }else{
      logL0 = MSloglik_fun(theta_h0, mdl_h0, k0);
      mdl_h0_tmp["P"] = P_h0;
    }
    logL1 = MSloglik_fun(theta_h1, mdl_h1, k1);
    double LRT_0 = -2*(logL0-logL1);
    // simulate under null hypothesis
    arma::vec LRN_tmp = LR_samp_dist(mdl_h0_tmp, k1, msmu, msvar, N, maxit, thtol, burnin, max_init, dist_converge_iter);
    pval = -MCpval(LRT_0, LRN_tmp, "geq");
  }
  return(pval);
}


// ==============================================================================
//' @title Monte Carlo Likelihood Ratio Test P-value Function 
//' 
//' 
//' @export
// [[Rcpp::export]]
double MMCLRpval_fun_max(arma::vec theta, List mdl_h0, List mdl_h1, bool msmu, bool msvar, int ar, int N, int maxit, double thtol, int burnin,  bool stationary_ind, double lambda, int max_init, int dist_converge_iter){
  double pval = -MMCLRpval_fun(theta, mdl_h0, mdl_h1, msmu, msvar, ar, N, maxit, thtol, burnin, stationary_ind, lambda, max_init, dist_converge_iter);
  return(pval);
}


