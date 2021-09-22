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
arma::vec LR_samp_dist(List mdl_h0, int k1, bool msmu, bool msvar, int N, int maxit = 500, double thtol = 1e-6){
  int k0 = mdl_h0["k"];
  arma::vec LRT_N(N,arma::fill::zeros);
  int ar = mdl_h0["ar"];
  bool getHess = FALSE;
  if (k0 == 1){
    bool inter = TRUE;
    int max_init = 10;
    for (int xn = 0; xn<N; xn++){
      double LRT_i;
      bool LRT_finite = FALSE;
      while (LRT_finite==FALSE){
        List y0_out = simuAR(mdl_h0);
        arma::vec y0 = y0_out["y"];
        List mdl_h0_tmp = ARmdl(y0, ar, inter);
        List mdl_h1_tmp = MSARmdl(y0, ar, k1, msmu, msvar, maxit, thtol, getHess, max_init);
        // test stat
        double l_0 = mdl_h0_tmp["logLike"];
        double l_1 = mdl_h1_tmp["logLike"];
        LRT_i = -2*(l_0 - l_1);
        LRT_finite = arma::is_finite(LRT_i);
      }
      LRT_N(xn) = LRT_i;
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
      while (LRT_finite==FALSE){
      List y0_out = simuMSAR(mdl_h0);
        arma::vec y0 = y0_out["y"];
        List mdl_h0_tmp = MSARmdl(y0, ar, k0, msmu, msvar, maxit, thtol, getHess);
        List mdl_h1_tmp = MSARmdl(y0, ar, k1, msmu, msvar, maxit, thtol, getHess);
        // test stat
        double l_0 = mdl_h0_tmp["logLike"];
        double l_1 = mdl_h1_tmp["logLike"];
        LRT_i = -2*(l_0 - l_1);
        LRT_finite = arma::is_finite(LRT_i);
      }
      LRT_N(xn) = LRT_i;
    } 
  }
  return(LRT_N);
}
// ==============================================================================
//' @title Monte Carlo Likelihood Ratio Test
//' 
//' 
//' @export
// [[Rcpp::export]]
List LR_MCTest(arma::vec Y, int ar, int k0 = 1, int k1 = 2, bool msmu = 1, bool msvar = 1, int N = 99, int maxit = 500, double thtol = 1e-6){
  List mdl_h0;
  List mdl_h1;
  bool getHess = FALSE;
  if (k0==1){
    mdl_h0 = ARmdl(Y, ar);
    mdl_h1 = MSARmdl(Y, ar, k1, msmu, msvar, maxit, thtol, getHess);  
  }else if (k0>1){
    mdl_h0 = MSARmdl(Y, ar, k0, msmu, msvar, maxit, thtol, getHess);  
    mdl_h1 = MSARmdl(Y, ar, k1, msmu, msvar, maxit, thtol, getHess);  
  }
  double logL0 = mdl_h0["logLike"];
  double logL1 = mdl_h1["logLike"];
  double LRT_0 = -2*(logL0-logL1);
  if (arma::is_finite(LRT_0)==FALSE){
    stop("LRT is not finite. Please check series\n");
  }
  arma::vec LRN = LR_samp_dist(mdl_h0, k1, msmu, msvar, N, maxit, thtol);
  double pval = MCpval(LRT_0, LRN, "geq");
  List LRMCTest_output;
  LRMCTest_output["mdl_h0"] = mdl_h0;
  LRMCTest_output["mdl_h1"] = mdl_h1;
  LRMCTest_output["LRT_0"] = LRT_0;
  LRMCTest_output["LRN"] = LRN;
  LRMCTest_output["pval"] = pval;
  return(LRMCTest_output);
}
// ==============================================================================
//' @title Bootstrap Likelihood Ratio Test
//' 
//' 
//' @export
// [[Rcpp::export]]
List LR_BootTest(arma::vec Y, int ar, int k0 = 1, int k1 = 2, bool msmu = 1, bool msvar = 1, int N = 1000, int maxit = 500, double thtol = 1e-6){
  List mdl_h0;
  List mdl_h1;
  bool getHess = FALSE;
  if (k0==1){
    mdl_h0 = ARmdl(Y, ar);
    mdl_h1 = MSARmdl(Y, ar, k1, msmu, msvar, maxit, thtol, getHess);  
  }else if (k0>1){
    mdl_h0 = MSARmdl(Y, ar, k0, msmu, msvar, maxit, thtol, getHess);  
    mdl_h1 = MSARmdl(Y, ar, k1, msmu, msvar, maxit, thtol, getHess);  
  }
  double logL0 = mdl_h0["logLike"];
  double logL1 = mdl_h1["logLike"];
  double LRT_0 = -2*(logL0-logL1);
  if (arma::is_finite(LRT_0)==FALSE){
    stop("LRT is not finite. Please check series\n");
  }
  arma::vec LRN = LR_samp_dist(mdl_h0, k1, msmu, msvar, N, maxit, thtol);
  double B = N;
  double pval = sum(LRN>LRT_0)/B; // [eq. 4.62] (Davidson & MacKinnon, 2004)
  List LRBootTest_output;
  LRBootTest_output["mdl_h0"] = mdl_h0;
  LRBootTest_output["mdl_h1"] = mdl_h1;
  LRBootTest_output["LRT_0"] = LRT_0;
  LRBootTest_output["LRN"] = LRN;
  LRBootTest_output["pval"] = pval;
  return(LRBootTest_output);
}


