#include <RcppArmadillo.h>
#include "methods.h"
//[[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;

// ==============================================================================f
//' @title Likelihood Ratio Test Statistic Sample Distribution
//' 
//' 
//' @export
// [[Rcpp::export]]
arma::vec LR_samp_dist(List mdl_h0, int k1, bool msmu, bool msvar, int N, int maxit, double thtol, 
                       int burnin, int max_init, int dist_converge_iter, int init_val_try_dist){
  // =============================================================================
  // ---------- Load R functions
  // =============================================================================
  Rcpp::Environment mstest("package:MSTest");
  Rcpp::Function MSmdl_EM = mstest["MSmdl_EM"];
  Rcpp::Function MSVARmdl_EM = mstest["MSVARmdl_EM"];
  Rcpp::Function ARmdl = mstest["ARmdl"];
  Rcpp::Function VARmdl = mstest["VARmdl"];
  // =============================================================================
  // ---------- Define required parameters
  // =============================================================================
  int k0 = mdl_h0["k"];
  int ar = mdl_h0["ar"];
  int q = mdl_h0["q"];
  bool getSE = FALSE;
  // =============================================================================
  // ---------- Optimization options
  // =============================================================================
  List control;
  control["msmu"] = msmu;
  control["msvar"] = msvar;
  control["maxit"] = maxit;
  control["thtol"] = thtol;
  control["getSE"] = getSE;
  control["max_init"] = max_init;
  control["max_init"] = max_init;
  control["use_diff_init"] = init_val_try_dist;
  // =============================================================================
  // ---------- Simulate test statistic under null hypothesis
  // =============================================================================
  arma::vec LRT_N(N,arma::fill::zeros);
  if ((k0 == 1) and (q==1)){
    // ----- Simulate linear AR model and estimate MS and AR model
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
          //List mdl_h0_tmp = ARmdl_cpp(y0, ar, inter, getSE);
          List mdl_h0_tmp = ARmdl(y0, ar);
          List mdl_h1_tmp = MSmdl_EM(y0, ar, k1, control);
          // test stat
          double l_0 = mdl_h0_tmp["logLike"];
          double l_1 = mdl_h1_tmp["logLike"];
          int mdl_h1_iter = mdl_h1_tmp["iterations"];
          LRT_i = -2*(l_0 - l_1);
          // LRT_finite = arma::is_finite(LRT_i);
          LRT_finite = ((arma::is_finite(LRT_i)) and (LRT_i>=0));
          LRT_converge = (mdl_h1_iter<maxit);
        }
        converge_iter +=1;
      }
      LRT_N(xn) = LRT_i;
      if (converge_iter==dist_converge_iter){
        warning("Warning: Some simulations had models that did not converge. Try using higher 'maxit'.");
      }
    } 
  }else if ((k0>1) and (q==1)){
    // ----- Simulate MS model and estimate MS models
    List mdl_h0_c = clone(mdl_h0);
    if (msmu == FALSE){
      arma::vec muk(k0, arma::fill::ones);
      double mu_h0 = mdl_h0["mu"];
      mdl_h0_c["mu"] = muk*mu_h0;
    }
    if (msvar == FALSE){
      arma::vec stdevk(k0, arma::fill::ones);
      double stdev_h0 = mdl_h0["stdev"];
      mdl_h0_c["stdev"] = stdevk*stdev_h0;
    }
    for (int xn = 0; xn<N; xn++){
      double LRT_i;
      bool LRT_finite = FALSE;
      bool LRT_converge = FALSE;
      int converge_iter = 0;
      while ((LRT_converge==FALSE) and (converge_iter<dist_converge_iter)){
        while (LRT_finite==FALSE){
          List y0_out = simuMSAR(mdl_h0_c, burnin);
          arma::vec y0 = y0_out["y"];
          List mdl_h0_tmp = MSmdl_EM(y0, ar, k0, control);
          List mdl_h1_tmp = MSmdl_EM(y0, ar, k1, control);
          // test stat
          double l_0 = mdl_h0_tmp["logLike"];
          double l_1 = mdl_h1_tmp["logLike"];
          int mdl_h0_iter = mdl_h0_tmp["iterations"];
          int mdl_h1_iter = mdl_h1_tmp["iterations"];
          LRT_i = -2*(l_0 - l_1);
          // LRT_finite = arma::is_finite(LRT_i);
          LRT_finite = ((arma::is_finite(LRT_i)) and (LRT_i>=0));
          LRT_converge = ((mdl_h0_iter<maxit) and (mdl_h1_iter<maxit));
        }
        converge_iter +=1;
      }
      LRT_N(xn) = LRT_i;
      if (converge_iter==dist_converge_iter){
        warning("Warning: Some simulations had models that did not converge. Try using higher 'maxit'.");
      }
    } 
  }else if ((k0==1) and (q>1)){
    // ----- Simulate linear VAR model and estimate MSVAR and VAR model
    bool inter = TRUE;
    for (int xn = 0; xn<N; xn++){
      double LRT_i;
      bool LRT_finite = FALSE;
      bool LRT_converge = FALSE;
      int converge_iter = 0;
      while ((LRT_converge==FALSE) and (converge_iter<dist_converge_iter)){
        while (LRT_finite==FALSE){
          List y0_out = simuVAR(mdl_h0, burnin);
          arma::mat y0 = y0_out["y"];
          //List mdl_h0_tmp = VARmdl(y0, ar, inter, getSE);
          List mdl_h0_tmp = VARmdl(y0, ar);
          List mdl_h1_tmp = MSVARmdl_EM(y0, ar, k1, control);
          // test stat
          double l_0 = mdl_h0_tmp["logLike"];
          double l_1 = mdl_h1_tmp["logLike"];
          int mdl_h1_iter = mdl_h1_tmp["iterations"];
          LRT_i = -2*(l_0 - l_1);
          // LRT_finite = arma::is_finite(LRT_i);
          LRT_finite = ((arma::is_finite(LRT_i)) and (LRT_i>=0));
          LRT_converge = (mdl_h1_iter<maxit);
        }
        converge_iter +=1;
      }
      LRT_N(xn) = LRT_i;
      if (converge_iter==dist_converge_iter){
        warning("Warning: Some simulations had models that did not converge. Try using higher 'maxit'.");
      }
    } 
  }else if ((k0>1) and (q>1)){
    // ----- Simulate  MS-VAR model and estimate MS-VAR models
    List mdl_h0_c = clone(mdl_h0);
    if (msmu == FALSE){
      arma::vec repk(k0,arma::fill::ones);
      arma::mat muk(k0, q, arma::fill::ones);
      arma::vec mu_h0 = mdl_h0["mu"];
      muk.rows(0,k0-1) = repk*trans(mu_h0);
      mdl_h0_c["mu"] = muk;
    }
    if (msvar == FALSE){
      List sigmak(k0);
      for (int xk = 0; xk<k0; xk++){
        sigmak[xk] = mdl_h0["sigma"];
      }
      mdl_h0_c["sigma"] = sigmak;
    }
    for (int xn = 0; xn<N; xn++){
      double LRT_i;
      bool LRT_finite = FALSE;
      bool LRT_converge = FALSE;
      int converge_iter = 0;
      while ((LRT_converge==FALSE) and (converge_iter<dist_converge_iter)){
        while (LRT_finite==FALSE){
          List y0_out = simuMSVAR(mdl_h0_c, burnin);
          arma::mat y0 = y0_out["y"];
          List mdl_h0_tmp = MSVARmdl_EM(y0, ar, k0, control);
          List mdl_h1_tmp = MSVARmdl_EM(y0, ar, k1, control);
          // test stat
          double l_0 = mdl_h0_tmp["logLike"];
          double l_1 = mdl_h1_tmp["logLike"];
          int mdl_h0_iter = mdl_h0_tmp["iterations"];
          int mdl_h1_iter = mdl_h1_tmp["iterations"];
          LRT_i = -2*(l_0 - l_1);
          // LRT_finite = arma::is_finite(LRT_i);
          LRT_finite = ((arma::is_finite(LRT_i)) and (LRT_i>=0));
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
double MMCLRpval_fun(arma::vec theta, List mdl_h0, List mdl_h1, bool msmu, bool msvar, int ar, int N, int maxit,
                     double thtol, int burnin, bool stationary_ind, double lambda, int max_init, int dist_converge_iter, 
                     int init_val_try_dist, int workers){
  Rcpp::Environment mstest("package:MSTest");
  Rcpp::Function LR_samp_dist_par = mstest["LR_samp_dist_par"];
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
  int q = mdl_h0["q"];
  arma::vec theta_h0_tmp = mdl_h0["theta"];
  arma::vec theta_h1_tmp = mdl_h1["theta"];
  int theta_h0_length = theta_h0_tmp.n_elem;
  int theta_h1_length = theta_h1_tmp.n_elem;
  arma::vec theta_h0 = theta.subvec(0,theta_h0_length-1);
  arma::vec theta_h1 = theta.subvec(theta_h0_length,theta_h0_length+theta_h1_length-1);
  // ----- Stationary constraint (i.e., only consider theta that result in stationary process) 
  if ((stationary_ind==TRUE) and (ar>0) and (q==1)){
    Rcpp::Function polyroot("polyroot");  
    Rcpp::Function Mod("Mod");  
    arma::vec poly_fun_h0(ar+1, arma::fill::ones);
    arma::vec poly_fun_h1(ar+1, arma::fill::ones);
    poly_fun_h0.subvec(1,ar) = -theta_h0.subvec(2+msmu*(k0-1)+msvar*(k0-1), 2+msmu*(k0-1)+msvar*(k0-1)+ar-1);
    poly_fun_h1.subvec(1,ar) = -theta_h1.subvec(2+msmu*(k1-1)+msvar*(k1-1), 2+msmu*(k1-1)+msvar*(k1-1)+ar-1);
    arma::vec roots_h0 = as<arma::vec>(Mod(wrap(as<ComplexVector>(polyroot(wrap(poly_fun_h0))))));
    arma::vec roots_h1 = as<arma::vec>(Mod(wrap(as<ComplexVector>(polyroot(wrap(poly_fun_h1))))));
    non_stationary_const = ((roots_h0.min()<=1) or (roots_h1.min()<=1));
  }else if ((stationary_ind==TRUE) and (ar>0) and (q>1)){
    Rcpp::Function Mod("Mod");  
    // phi indicators
    arma::vec theta_phi_ind_h0 = mdl_h0["theta_phi_ind"];
    arma::vec theta_phi_ind_h1 = mdl_h1["theta_phi_ind"];
    // phi vectors
    arma::vec phi_vec_h0 = theta_h0.elem(find(theta_phi_ind_h0));
    arma::vec phi_vec_h1 = theta_h1.elem(find(theta_phi_ind_h1)); 
    // phi matrices
    arma::mat phi_h0 = reshape(phi_vec_h0, q*ar, q);
    arma::mat phi_h1 = reshape(phi_vec_h1, q*ar, q);
    // companion matrices
    arma::mat F0_tmp = trans(phi_h0);
    arma::mat F1_tmp = trans(phi_h1);
    arma::mat diagmat = arma::eye(q*(ar-1),q*(ar-1));
    arma::mat diagzero(q*(ar-1), q, arma::fill::zeros);
    arma::mat Mn = join_rows(diagmat,diagzero);
    arma::mat F0 = join_cols(F0_tmp,Mn);
    arma::mat F1 = join_cols(F1_tmp,Mn);
    // eigen values
    arma::cx_vec eig_0 = eig_gen(F0);
    arma::cx_vec eig_1 = eig_gen(F1);
    bool stationary_0 = all(as<arma::vec>(Mod(eig_0))<1);
    bool stationary_1 = all(as<arma::vec>(Mod(eig_1))<1);
    non_stationary_const = ((stationary_0==FALSE) or (stationary_1==FALSE));
  }
  // ----- Checking that P matrix columns sum to 1
  if (k0>1){
    arma::vec theta_P_ind_h0 = mdl_h0["theta_P_ind"];
    arma::vec P_vec_h0 = theta_h0.elem(find(theta_P_ind_h0));
    P_h0 = reshape(P_vec_h0, k0, k0);
    P_h0_colsum_const = any(abs(arma::sum(P_h0,0)-1)>thtol); 
  }
  arma::vec theta_P_ind_h1 = mdl_h1["theta_P_ind"];
  arma::vec P_vec_h1 = theta_h1.elem(find(theta_P_ind_h1));
  P_h1 = reshape(P_vec_h1, k1, k1); 
  P_h1_colsum_const = any(abs(arma::sum(P_h1,0)-1)>thtol);
  // ----- Compute pval 
  if ((P_h0_colsum_const==TRUE) or (P_h1_colsum_const==TRUE) or (non_stationary_const==TRUE)){
    // If either transition matrix columns do not sum to 1 OR (stationary_ind == TRUE AND non_stationary_const == TRUE), pval is a positive constant
    pval = lambda*(P_h0_colsum_const + P_h1_colsum_const + non_stationary_const);
  }else{
    // If transition matrices' columns sum to 1 AND (stationary_ind==FALSE OR non_stationary_const == FALSE), pval is computed. 
    List mdl_h0_tmp = clone(mdl_h0);
    // edit theta vector
    mdl_h0_tmp["theta"] = theta_h0;
    arma::vec theta_mu_ind_h0 = mdl_h0["theta_mu_ind"];
    arma::vec theta_mu_h0 = theta_h0.elem(find(theta_mu_ind_h0));
    arma::vec theta_sig_ind_h0 = mdl_h0["theta_sig_ind"];
    arma::vec theta_sig_h0 = theta_h0.elem(find(theta_sig_ind_h0));
    if (q==1){
      mdl_h0_tmp["mu"] = theta_mu_h0;
      mdl_h0_tmp["stdev"] = sqrt(theta_sig_h0);
      if (ar>0){
        arma::vec theta_phi_ind_h0 = mdl_h0["theta_phi_ind"];
        // phi vector
        arma::vec phi_new = theta_h0.elem(find(theta_phi_ind_h0));
        mdl_h0_tmp["phi"] = phi_new;
      }
    }else if (q>1){
      arma::mat mu_new = trans(reshape(theta_mu_h0,q,k0));
      mdl_h0_tmp["mu"] = mu_new;
      if (k0==1){
        mdl_h0_tmp["sigma"] = covar_unvech(theta_sig_h0,q);
      }else if (k0>1){
        arma::mat sig_new_tmp = trans(reshape(theta_sig_h0, (q*(q+1))/2, k0));
        List sig_new(k0);
        for (int xk = 0; xk<k0; xk++){
          arma::vec sig_tmp_k = trans(sig_new_tmp.row(xk));
          sig_new[xk] = covar_unvech(sig_tmp_k, q);
        }
        mdl_h0_tmp["sigma"] = sig_new;
      } 
      if (ar>0){
        arma::vec theta_phi_ind_h0 = mdl_h0["theta_phi_ind"];
        // phi vectors
        arma::vec phi_vec_h0 = theta_h0.elem(find(theta_phi_ind_h0));
        // phi matrices
        arma::mat phi_h0 = reshape(phi_vec_h0, q*ar, q);
        mdl_h0_tmp["phi"] = phi_h0;
      } 
    }
    // compute test stat
    if ((k0==1) and (q==1)){
      // MS model & linear model under null hypothesis
      logL0 = logLike_AR(theta_h0, mdl_h0);  
      logL1 = MSloglik_fun(theta_h1, mdl_h1, k1);
    }else if ((k0>1) and (q==1)){
      // MS models
      logL0 = MSloglik_fun(theta_h0, mdl_h0, k0);
      mdl_h0_tmp["P"] = P_h0;
      logL1 = MSloglik_fun(theta_h1, mdl_h1, k1);
    }else if ((k0==1) & (q>1)){
      // MSVAR model & linear model under null hypothesis
      logL0 = logLike_VAR(theta_h0, mdl_h0);  
      logL1 = MSVARloglik_fun(theta_h1, mdl_h1, k1);
    }else if ((k0>1) & (q>1)){
      // MSVAR models
      logL0 = MSVARloglik_fun(theta_h0, mdl_h0, k0);
      mdl_h0_tmp["P"] = P_h0;
      logL1 = MSVARloglik_fun(theta_h1, mdl_h1, k1);
    }else{
      stop("Number of regimes under the Null and number of columns of Y must be >=1");
    }
    double LRT_0 = -2*(logL0-logL1);
    // simulate under null hypothesis
    arma::vec LRN_tmp;
    if(workers>0){
      LRN_tmp = as<arma::vec>(LR_samp_dist_par(mdl_h0_tmp, k1, msmu, msvar, N, maxit, thtol, burnin, max_init, dist_converge_iter, init_val_try_dist, workers));
    }else{
      LRN_tmp = LR_samp_dist(mdl_h0_tmp, k1, msmu, msvar, N, maxit, thtol, burnin, max_init, dist_converge_iter, init_val_try_dist);
    }
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
double MMCLRpval_fun_max(arma::vec theta, List mdl_h0, List mdl_h1, bool msmu, bool msvar, int ar, int N, int maxit, 
                         double thtol, int burnin,  bool stationary_ind, double lambda, int max_init, int dist_converge_iter, 
                         int init_val_try_dist, int workers){
  double pval = -MMCLRpval_fun(theta, mdl_h0, mdl_h1, msmu, msvar, ar, N, maxit, thtol, burnin, stationary_ind, lambda, max_init, dist_converge_iter, init_val_try_dist, workers);
  return(pval);
}


