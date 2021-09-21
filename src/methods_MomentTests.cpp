#include <RcppArmadillo.h>
#include "models.h"
#include "methods.h"
//[[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;
// ==============================================================================
//' @title Calculate Dufour & Luger (2017) Moment-Based Test-Statistics 
//'
//' @description This function computes the four momment-based test-statistics (eq. 11 - 14 in paper) 
//' for a given series. The series should be the residuals from an AR model. 
//' 
//'
//' @param ehat vector containing series of residuals from model.
//' 
//' @return The four test statistics (eq. 11 - 14 in paper)
//' 
//' @references Dufour, J. M., & Luger, R. (2017). Identification-robust moment-based 
//' tests for Markov switching in autoregressive models. Econometric Reviews, 36(6-9), 713-727.
//' 
//' @export
// [[Rcpp::export]]
arma::vec calc_DLmoments(arma::vec ehat){
  arma::vec stats(4, arma::fill::zeros);
  // Mean Moment  
  arma::vec idx1 = ehat.elem(find(ehat < 0));
  arma::vec idx2 = ehat.elem(find(ehat > 0));
  double m1 = sum(idx1)/idx1.n_elem;
  double m2 = sum(idx2)/idx2.n_elem;
  double s1 = sum( (idx1 - m1) % (idx1 - m1) )/ idx1.n_elem;
  double s2 = sum( (idx2 - m2) % (idx2 - m2) )/ idx2.n_elem;
  stats(0) = std::abs(m2 - m1)/sqrt(s2 + s1);
  // Variance Moment 
  arma::vec ehat2 = ehat % ehat;
  double var = mean(ehat2);
  arma::vec idx3 = ehat2(find(ehat2 < var));
  arma::vec idx4 = ehat2(find(ehat2 > var));
  double v1 = sum(idx3)/idx3.n_elem;
  double v2 = sum(idx4)/idx4.n_elem;
  stats(1) = v2/v1;
  // Skewness Moment 
  arma::vec z = ehat/sqrt(var);
  stats(2) = std::abs( mean(z % z % z) ) ;
  // Kurtosis Moment
  stats(3) = std::abs( mean(z % z % z % z) - 3);
  // output 
  return(stats);
}
// ==============================================================================
//' @title Calculate Quantile Moment-Based Test-Statistics 
//'
//' @description This function computes the four momment-based test-statistics (eq. 11 - 14 in paper) 
//' for a given series. The series should be the residuals from an AR model. 
//' 
//' @export
// [[Rcpp::export]]
arma::vec calc_Qmoments(arma::vec eps, int k1){
  Rcpp::Function quant("quantile");
  // cluster observations
  arma::vec qseq(k1, arma::fill::ones);
  qseq = cumsum(qseq/k1);
  arma::vec lead0(1, arma::fill::zeros);
  qseq = join_vert(lead0, qseq);
  arma::vec quants = quantile(eps, qseq);
  // Get stats based on means and variance
  List indx(k1);
  arma::vec mu(k1, arma::fill::zeros);
  arma::vec std(k1, arma::fill::zeros);
  arma::vec eta(k1, arma::fill::zeros);
  for (int xi = 0; xi<k1; xi++){
    arma::vec vec_tmp = eps.rows(find((eps>=quants(xi)) and (eps<(quants(xi+1)+0.0000000001))));
    indx[xi] = vec_tmp;
    int len_tmp = vec_tmp.n_elem;
    mu(xi) = sum(vec_tmp)/len_tmp;
    std(xi) = sum((vec_tmp-mu(xi))%(vec_tmp-mu(xi)))/len_tmp;
    eta(xi) = sum(vec_tmp%vec_tmp)/len_tmp;
  }
  arma::vec S1((k1*(k1-1))/2, arma::fill::zeros);
  arma::vec S2((k1*(k1-1))/2, arma::fill::zeros);
  int count = 0;
  for (int xi = 0; xi<k1; xi++){
    for (int xj = (xi+1); xj<k1; xj++){
      S1(count) = abs(mu(xj)-mu(xi))/sqrt(std(xj)+std(xi));
      S2(count) = eta(xj)/eta(xi);
      count = count+1;
    }
  }
  arma::vec S = join_vert(S1, S2);
  arma::vec Ssk(2, arma::fill::zeros);
  // Stat based on Skewness
  int Tsize = eps.n_elem;
  //double var = sum((eps-mean(eps))%(eps-mean(eps)))/(Tsize-1);
  double var = sum((eps)%(eps))/(Tsize-1);
  arma::vec z = eps/sqrt(var);
  Ssk(0) = abs(mean(z%z%z));
  // State based on Kurtosis 
  Ssk(1) = abs(mean(z%z%z%z) - 3);
  S = join_vert(S, Ssk);
  return(S);
}
// ==============================================================================
//' @title Test-statistics from simulated data
//'
//' @description This function computes simulates residuals from a standard normal distribution
//' and calculates the four momment-based test-statistics (eq. 11 - 14 in paper) under the 
//' null hypothesis using this simulated data.
//' 
//' @param t length of sample size for simulation 
//' @param number of simulated samples
//' 
//' @return The four test statistics (eq. 11 - 14 in paper)
//' 
//' @references Dufour, J. M., & Luger, R. (2017). Identification-robust moment-based 
//' tests for Markov switching in autoregressive models. Econometric Reviews, 36(6-9), 713-727.
//' 
//' @export
// [[Rcpp::export]]
arma::mat sim_DLmoments(int Tsize,int N){
  arma::mat stats(N, 4, arma::fill::zeros);
  // loop and perform simultion each time 
  for (int isim = 0; isim <N; isim++){
    arma::vec esim = arma::randn<arma::vec>(Tsize);
    stats.row(isim) = trans(calc_DLmoments(esim - mean(esim)));
  }
  return(stats);
}
// ==============================================================================
//' @title Quantile-Based Test-Statistics From Simulated Data  
//'
//' @description 
//' 
//' @export
// [[Rcpp::export]]
arma::mat sim_Qmoments(int Tsize, int N, List mdl_h0, int k1){
  int k0 = mdl_h0["k"];
  int ar = mdl_h0["ar"];
  arma::mat SN(N,k1*(k1-1)+2,arma::fill::zeros);
  if (k0==1){
    double mu_h0 = mdl_h0["mu"];
    double sig_h0 = mdl_h0["sigma"];
    double std_h0 = sqrt(sig_h0);
    for (int ik = 0; ik<N; ik++){
      arma::vec z_null = mu_h0 + std_h0*arma::randn<arma::vec>(Tsize);
      arma::vec eps_null = z_null - mean(z_null);
      arma::vec mmtmp = calc_Qmoments(eps_null, k1);
      SN.row(ik) = trans(mmtmp);
    }
  }else{
    if (ar>0){
      for (int ik = 0; ik<N; ik++){
        List simudat  = simuMSAR(mdl_h0, "markov");
        arma::vec y_null= simudat["y"];
        bool intercept = TRUE;
        List mdl_ar = ARmdl(y_null, ar, intercept);
        arma::vec ytmp = mdl_ar["y"];
        arma::mat xtmp = mdl_ar["x"];
        arma::vec phi = mdl_ar["phi"];
        arma::vec z_null = ytmp - xtmp*phi;
        arma::vec eps_null = z_null - mean(z_null);
        arma::vec mmtmp = calc_Qmoments(eps_null, k1);
        SN.row(ik) = trans(mmtmp);
      }
    }else{
      for (int ik = 0; ik<N; ik++){
        List simudat   = simuMS(mdl_h0, "markov");
        arma::vec z_null = simudat["y"];
        arma::vec eps_null = z_null - mean(z_null);
        arma::vec mmtmp = calc_Qmoments(eps_null, k1);
        SN.row(ik)=trans(mmtmp);
      }
    }
  }
  return(SN);
}
// ==============================================================================
//' @title Combine p-values 
//'
//' @description This function is used to combine the p-values as in eq 17 and 18 of Dufour & Luger (2017).
//' The input parameter \emph{type} can be used to used to specify the method for combining 
//' the pvalues. If set to "min" the min method of combining p-values is used as in Fisher (1932) 
//' and Pearson (1933). If set to "prod" the product of p-values is used as in Tippett (1931) 
//' and Wilkinson (1951).
//' 
//' @param s0 test-statistic under the alternative
//' @param sN test-statistics under the null
//' @param param output from \emph{approxDist} which obtains the parameters needed in eq. 16 which is used for 
//' combining p-values.
//' @param N total number of test statistics (i.e. simulated + observed = N)
//' @param type the type of method used to combine p-values. If set to "min" the min method of combining p-values is used as in Fisher (1932) 
//' and Pearson (1933). If set to "prod" the product of p-values is used as in Tippett (1931) 
//' and Wilkinson (1951).
//' 
//' @return The four test statistics (eq. 11 - 14 in paper)
//' 
//' @references Dufour, J. M., & Luger, R. (2017). Identification-robust moment-based 
//' tests for Markov switching in autoregressive models. Econometric Reviews, 36(6-9), 713-727.
//' @references Tippett, L. (1931). The Method of Statistics. London: Williams & Norgate.
//' @references Wilkinson, B. (1951). A statistical consideration in psychological research. 
//' Psychology Bulletin 48:156–158.
//' @references Pearson, K. (1933). On a method of determining whether a sample of size n
//'  supposed to have been drawn from a parent population having a known probability integral has probably been drawn at random. Biometrika 25:379–410.
//' @references Fisher, R. (1932). Statistical Methods for Research Workers. Edinburgh: 
//' Oliver and Boyd.
//' 
//' @export
// [[Rcpp::export]]
arma::vec combine_stat(arma::vec s0,arma::mat sN, arma::mat params, std::string type){
  int N = sN.n_rows;
  int nc = sN.n_cols;
  arma::mat G0(1,nc,arma::fill::zeros);
  arma::mat Gx(N,nc,arma::fill::zeros);
  for (int im = 0; im<nc; im++){
    G0(0,im) = 1 - (exp(params(0,im)+params(1,im)*s0(im))/(1+exp(params(0,im)+params(1,im)*s0(im))));
    Gx.col(im) = 1 - (exp(params(0,im)+params(1,im)*sN.col(im))/(1+exp(params(0,im)+params(1,im)*sN.col(im))));
  }
  arma::vec Fx(N+1,arma::fill::zeros);
  double F0 = 0;
  if (type=="min"){
    // MC p-value using Tippett (1931) &  Wilkinson (1951) combination method. 
    F0 = 1 - min(G0.row(0));
    for (int isim = 0; isim<N; isim++){
      Fx(isim) = 1 - min(Gx.row(isim));
    }
  }
  if (type=="prod"){
    // MC p-value using Fisher (1932) &  Pearson (1933) combination method. 
    F0 = 1 - (G0(0)*G0(1)*G0(2)*G0(3));
    for (int isim = 0; isim<N; isim++){
      Fx(isim) = 1 - (Gx(isim,0)*Gx(isim,1)*Gx(isim,2)*Gx(isim,3));
    }
  }
  Fx(N) = F0;
  return(Fx);
}
// ==============================================================================
//' @title Calculate combined p-value test statistic 
//'
//' @description This function computes test-statistics for observed and simulated samples. 
//' It begins by calculating it for the obsted data, then generates test-statistics from 
//' simulated data to approximate null distribution and then calculates the p-values and 
//' combines them according to eq. 17 or 18 of Dufour & Luger (2017).
//' 
//' @param ezt observed series
//' @param N total number of test statistics (i.e. simulated + observed = N)
//' @param param output from \emph{approxDist} which obtains the parameters needed in eq. 16 which is used for 
//' combining p-values.
//' 
//' @return The four test statistics (eq. 11 - 14 in paper)
//' 
//' @references Dufour, J. M., & Luger, R. (2017). Identification-robust moment-based 
//' tests for Markov switching in autoregressive models. Econometric Reviews, 36(6-9), 713-727.
//' 
//' @export
// [[Rcpp::export]]
arma::vec calc_DLmcstat(arma::vec ezt, int N, arma::mat params, Rcpp::String type = "min"){
  // Get length of series 
  int Tsize = ezt.n_elem;
  // calulate moments of data 
  arma::vec S0 = calc_DLmoments(ezt);
  // Simulated data 
  arma::mat SN = sim_DLmoments(Tsize,N); // must leaves as N-1.
  // Get individual moment p-values 
  arma::vec Fx(N+1, arma::fill::zeros);
  if (type=="prod"){
    Fx = combine_stat(S0, SN, params, "prod");
  }else if (type=="min"){
    Fx  = combine_stat(S0, SN, params, "min");   
  }else{
    Rcerr << "type must be: 'min' or 'prod'.\n";
  }
  return(Fx);
}
// ==============================================================================
//' @title Calculate Quantile-Based Combined P-Value Test Statistic 
//'
//' 
//' @param ezt observed series
//' 
//' @export
// [[Rcpp::export]]
arma::vec calc_Qmcstat(arma::vec ezt, int N, arma::mat params, List mdl_h0, int k1, Rcpp::String type = "min"){
  int Tsize = ezt.n_elem;
  arma::vec S0 = calc_Qmoments(ezt, k1);
  // Simulated data 
  arma::mat SN = sim_Qmoments(Tsize, N, mdl_h0, k1);
  // Get individual moment p-values 
  arma::vec Fx(N+1, arma::fill::zeros);
  if (type=="prod"){
    Fx = combine_stat(S0, SN, params, "prod");
  }else if (type=="min"){
    Fx  = combine_stat(S0, SN, params, "min");   
  }else{
    Rcerr << "type must be: 'min' or 'prod'.\n";
  }
  return(Fx);
}
// ==============================================================================
//' @title Loop for \emph{approxDist}
//'
//' @description This function performs the loop in \emph{approxDist}. It is written in C++ for better 
//' performance in terms of speed.
//' 
//' @param SN2 matrix of Tx4 test-statistics
//' 
//' @return The four test statistics from smulated data in [0,1]. Used for NLS to get 
//' params needed to combine p-values
//' 
//' @references Dufour, J. M., & Luger, R. (2017). Identification-robust moment-based 
//' tests for Markov switching in autoregressive models. Econometric Reviews, 36(6-9), 713-727.
//' 
//' @export    
// [[Rcpp::export]]   
arma::mat approx_dist_loop(arma::mat SN2){
  double N2 = SN2.col(0).n_elem;
  double npar = SN2.row(0).n_elem;
  arma::mat Fx(N2,npar,arma::fill::ones);
  for (int ik = 0; ik<N2; ik++){
    for (int ikk =0; ikk<npar; ikk++){
      Fx(ik,ikk)=sum(SN2(ik,ikk)>SN2.col(ikk))/N2;
    }
  }
  return(Fx);
}
// ==============================================================================
//' @title Monte-Carlo Moment-based test for MS AR model
//'
//' This function performs the Local Monte-Carlo Moment-Based test for
//' MS AR models presented in Dufour & Luger (2017) (i.e when no nuissance 
//' parameters are present). 
//'
//' @param Y Series to be tested 
//' @param p Order of autoregressive components AR(p).
//' @param x exogenous variables if any. Test in Dufour & Luger is model for AR lags
//' @param N number of samples
//' @param N2 number of simulations when approximating distribution used to combine 
//' p-values (eq. 16).
//'
//' @return List with model and test results.
//' 
//' @references Dufour, J. M., & Luger, R. (2017). Identification-robust moment-based 
//' tests for Markov switching in autoregressive models. Econometric Reviews, 36(6-9), 713-727.
//' 
//' @export
// [[Rcpp::export]]
List DLMCtest(arma::vec Y, int ar = 0, int N = 99, int simdist_N = 10000){
  Rcpp::Environment mstest("package:MSTest");
  Rcpp::Function approxDistDL = mstest["approxDistDL"];
  int Tsize = Y.n_elem;
  arma::vec z(Tsize-ar, arma::fill::zeros);
  List mdl_out;
  if(ar>0){
    bool intercept = TRUE;
    mdl_out = ARmdl(Y, ar, intercept);
    arma::vec y = mdl_out["y"];
    arma::mat x = mdl_out["x"];
    arma::vec phi = mdl_out["phi"];
    z = y - x*phi;
  }else{
    z = Y;
  }
  // --------- Get parameters from approximated distribution ----------
  arma::mat params = as<arma::mat>(approxDistDL(Tsize-ar, simdist_N));
  // -------------------------- Get P-Values --------------------------
  arma::vec eps = z - mean(z);
  arma::vec Fmin = calc_DLmcstat(eps, N, params, "min");
  arma::vec Fprod = calc_DLmcstat(eps, N, params, "prod");
  // ----- Obtain p-value
  double Fmin0 = Fmin(N);
  double Fprod0 = Fprod(N);
  arma::vec FminSim = Fmin.subvec(0,N-1);
  arma::vec FprodSim = Fprod.subvec(0,N-1);
  double pval_min = MCpval(Fmin0, FminSim, "geq");
  double pval_prod = MCpval(Fprod0, FprodSim, "geq");
  arma::mat LMC_ans = join_rows(Fmin,Fprod);
  List DLtest_output;
  DLtest_output["eps"] = eps;
  if (ar>0){
    DLtest_output["ARmdl"] = mdl_out;
  }
  DLtest_output["params"] = params;
  DLtest_output["LMC_ans"] = LMC_ans;
  DLtest_output["Fmin"] = Fmin0;
  DLtest_output["Fprod"] = Fprod0;
  DLtest_output["p-value_min"] = pval_min;
  DLtest_output["p-value_prod"] = pval_prod;
  return(DLtest_output);
}
// ==============================================================================
//' @title Monte-Carlo Moment-based test for MS AR model
//'
//' This function performs the Local Monte-Carlo Moment-Based test for
//' MS AR models presented in Dufour & Luger (2017) (i.e when no nuissance 
//' parameters are present). 
//'
//' @param Y Series to be tested 
//' @param p Order of autoregressive components AR(p).
//' @param x exogenous variables if any. Test in Dufour & Luger is model for AR lags
//' @param N number of samples
//' @param N2 number of simulations when approximating distribution used to combine 
//' p-values (eq. 16).
//'
//' @return List with model and test results.
//' 
//' @references Dufour, J. M., & Luger, R. (2017). Identification-robust moment-based 
//' tests for Markov switching in autoregressive models. Econometric Reviews, 36(6-9), 713-727.
//' 
//' @export
// [[Rcpp::export]]
List QMCtest(arma::vec Y, int ar = 0, int k0 = 1, int k1 = 2, int N = 99, int simdist_N = 10000, bool msmu = 1, bool msvar = 1, int maxit = 500, double thtol = 1e-8){
  Rcpp::Environment mstest("package:MSTest");
  Rcpp::Function approxDistQ = mstest["approxDistQ"];
  int Tsize = Y.n_elem;
  arma::vec z(Tsize-ar, arma::fill::zeros);
  List mdl_h0;
  List mdl_h1;
  if(ar>0){
    if (k0==1){
      bool intercept = TRUE;
      mdl_h0 = ARmdl(Y, ar, intercept);
      mdl_h0["k"] = k0;
      arma::vec y = mdl_h0["y"];
      arma::mat x = mdl_h0["x"];
      arma::vec phi = mdl_h0["phi"];
      z = y - x*phi;  
    }else{
      bool getHess = FALSE;
      mdl_h0 = MSARmdl(Y, ar, k0, msmu, msvar, maxit, thtol, getHess);
      arma::vec y = mdl_h0["y"];
      arma::mat x = mdl_h0["x"];
      arma::vec phi = mdl_h0["phi"];
      z = y - x*phi;  
    }
  }else{
    z = Y;
    mdl_h0["k"] = k0;
    mdl_h0["ar"] = 0;
    mdl_h0["mu"] = mean(z);
    mdl_h0["sigma"] = sum((z-mean(z))%(z-mean(z)))/(Tsize-1);
  }
  // --------- Get parameters from approximated distribution ----------
  arma::mat params = as<arma::mat>(approxDistQ(Tsize-ar, simdist_N, mdl_h0, k1));
  // -------------------------- Get P-Values --------------------------
  arma::vec eps = z - mean(z);
  arma::vec Fmin = calc_Qmcstat(eps, N, params, mdl_h0, k1, "min");
  arma::vec Fprod = calc_Qmcstat(eps, N, params, mdl_h0, k1, "prod");
  // ----- Obtain p-value
  double Fmin0 = Fmin(N);
  double Fprod0 = Fprod(N);
  arma::vec FminSim = Fmin.subvec(0,N-1);
  arma::vec FprodSim = Fprod.subvec(0,N-1);
  double pval_min = MCpval(Fmin0, FminSim, "geq");
  double pval_prod = MCpval(Fprod0, FprodSim, "geq");
  arma::mat LMC_ans = join_rows(Fmin,Fprod);
  List Qtest_output;
  Qtest_output["eps"] = eps;
  if (ar>0){
    Qtest_output["ARmdl"] = mdl_h0;
  }
  Qtest_output["params"] = params;
  Qtest_output["LMC_ans"] = LMC_ans;
  Qtest_output["Fmin"] = Fmin0;
  Qtest_output["Fprod"] = Fprod0;
  Qtest_output["p-value_min"] = pval_min;
  Qtest_output["p-value_prod"] = pval_prod;
  return(Qtest_output);
}
// ==============================================================================
//' @title MMC pvalue Function
//'
//' 
//' @references Dufour, J. M., & Luger, R. (2017). Identification-robust moment-based 
//' tests for Markov switching in autoregressive models. Econometric Reviews, 36(6-9), 713-727.
//' 
//' @export
// [[Rcpp::export]]
double DLMMCpval_fun(arma::vec theta, arma::vec y, arma::mat x, int N = 99, Rcpp::String type = "min", int simdist_N = 10000){
  Rcpp::Function polyroot("polyroot");  
  Rcpp::Function Mod("Mod");  
  Rcpp::Environment mstest("package:MSTest");
  Rcpp::Function approxDistDL = mstest["approxDistDL"];
  // ----- Stationary Inequality constraint
  int ar = theta.n_elem;
  arma::vec poly_fun(ar+1, arma::fill::ones);
  poly_fun.subvec(1,ar) = -theta;
  arma::vec roots = as<arma::vec>(Mod(wrap(as<ComplexVector>(polyroot(wrap(poly_fun))))));
  bool ineq_constraint = roots.min()<1;
  // ----- Transform data
  arma::vec z = y - x*theta;
  int Tsize = z.n_elem;
  arma::mat params = as<arma::mat>(approxDistDL(Tsize, simdist_N));
  // ----- Compute test stats
  arma::vec eps = z - mean(z);
  arma::vec Fx = calc_DLmcstat(eps, N, params, type);
  // ----- Obtain p-value
  double F0 = Fx(N);
  arma::vec FN = Fx.subvec(0,N-1);
  double pval = -MCpval(F0, FN, "geq") + ineq_constraint*100;
  return(pval);
}

// ==============================================================================
//' @title Maximized Monte-Carlo Moment-based test for MS AR model
//'
//' 
//' @references Dufour, J. M., & Luger, R. (2017). Identification-robust moment-based 
//' tests for Markov switching in autoregressive models. Econometric Reviews, 36(6-9), 713-727.
//' 
//' @export
// [[Rcpp::export]]
List DLMMCtest(arma::vec Y, int ar = 1, int N = 99, Rcpp::String method = "GenSA", int simdist_N = 10000){
  Rcpp::Environment mstest("package:MSTest");
  Rcpp::Environment gensa("package:GenSA");
  Rcpp::Function approxDistDL = mstest["approxDistDL"];
  Rcpp::Function GenSA = gensa["GenSA"];
  int Tsize = Y.n_elem;
  arma::vec z(Tsize-ar, arma::fill::zeros);
  List mdl_out;
  if(ar>0){
    bool intercept = TRUE;
    mdl_out = ARmdl(Y, ar, intercept);
    arma::vec y = mdl_out["y"];
    arma::mat x = mdl_out["x"];
    arma::vec phi = mdl_out["phi"];
    z = y - x*phi;
  }else{
    Rcerr << "No Nuisance parameters is model is not Autoregressive. Number of lags must be greater than 0.\n";
  }
  // --------- Get parameters from approximated distribution ----------
  arma::mat params = as<arma::mat>(approxDistDL(Tsize-ar, simdist_N));
  // -------------------------- Get P-Values --------------------------
  arma::vec eps = z - mean(z);
  arma::mat LMC_ans = calc_DLmcstat(eps, N, params);
  // ----- Obtain p-value
  double Fmin0 = LMC_ans(N,0);
  double Fprod0 = LMC_ans(N,1);
  arma::vec Fmin = LMC_ans.submat(0,0,N-1,0);
  arma::vec Fprod = LMC_ans.submat(0,1,N-1,1);
  double pval_min = MCpval(Fmin0, Fmin, "geq");
  double pval_prod = MCpval(Fprod0, Fprod, "geq");
  List DLtest_output;
  DLtest_output["eps"] = eps;
  if (ar>0){
    DLtest_output["ARmdl"] = mdl_out;
  }
  DLtest_output["params"] = params;
  DLtest_output["LMC_ans"] = LMC_ans;
  DLtest_output["Fmin"] = Fmin0;
  DLtest_output["Fprod"] = Fprod0;
  DLtest_output["p-value_min"] = pval_min;
  DLtest_output["p-value_prod"] = pval_prod;
  return(DLtest_output);
}

