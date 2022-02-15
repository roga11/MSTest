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
    // if (ar>0){
    //   bool msmu = mdl_h0["msmu"];
    //   bool msvar = mdl_h0["msvar"];
    //   List mdl_options = mdl_h0["optim_options"];
    //   int maxit = mdl_options["maxit"];
    //   double thtol = mdl_options["thtol"];
    //   bool getSE = mdl_h0["getSE"];
    //   for (int ik = 0; ik<N; ik++){
    //     List simudat  = simuMSAR(mdl_h0, "markov")
    //     arma::vec y_null = simudat["y"];
    //     bool intercept = TRUE;
    //     List mdl_ar = ARmdl(y_null, ar, intercept);
    //     //List mdl_ar = MSARmdl(y_null, ar, k0, msmu, msvar, maxit, thtol, getSE);
    //     arma::vec ytmp = mdl_ar["y"];
    //     arma::mat xtmp = mdl_ar["x"];
    //     arma::vec phi = mdl_ar["phi"];
    //     arma::vec z_null = ytmp - xtmp*phi;
    //     arma::vec eps_null = z_null - mean(z_null);
    //     arma::vec mmtmp = calc_Qmoments(eps_null, k1);
    //     SN.row(ik) = trans(mmtmp);
    //   }
    // }else{
    for (int ik = 0; ik<N; ik++){
      List simudat   = simuMS(mdl_h0, "markov");
      arma::vec z_null = simudat["y"];
      arma::vec eps_null = z_null - mean(z_null);
      arma::vec mmtmp = calc_Qmoments(eps_null, k1);
      SN.row(ik)=trans(mmtmp);
    }
    //}
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
    F0 = 1 - prod(G0.row(0));
    for (int isim = 0; isim<N; isim++){
      Fx(isim) = 1 - prod(Gx.row(isim));
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
//' @title Dufour & Luger (2017) moment-based MMC test p-value function to be minimized
//'
//' 
//' @references Dufour, J. M., & Luger, R. (2017). Identification-robust moment-based 
//' tests for Markov switching in autoregressive models. Econometric Reviews, 36(6-9), 713-727.
//' 
//' @export
// [[Rcpp::export]]
double DLMMCpval_fun(arma::vec theta, arma::vec y, arma::mat x, int N, int simdist_N, Rcpp::String pval_type, bool stationary_ind, double lambda){
  Rcpp::Environment mstest("package:MSTest");
  Rcpp::Function approxDistDL = mstest["approxDistDL"];
  bool stationary_constraint = FALSE;
  double pval;
  // ----- Stationary constraint (i.e., only consider theta that result in stationary process) 
  if (stationary_ind==TRUE){
    Rcpp::Function polyroot("polyroot");  
    Rcpp::Function Mod("Mod");  
    int ar = theta.n_elem;
    arma::vec poly_fun(ar+1, arma::fill::ones);
    poly_fun.subvec(1,ar) = -theta;
    arma::vec roots = as<arma::vec>(Mod(wrap(as<ComplexVector>(polyroot(wrap(poly_fun))))));
    stationary_constraint = roots.min()<=1;
  }
  if (stationary_constraint){
    // If stationary_ind == TRUE AND ineq_constraint == TRUE (i.e. non-stationary process), then pval = constraint
    pval = lambda*stationary_constraint;
  }else{
    // If stationary_ind == FALSE OR ineq_constraint == FALSE (i.e. non-stationary process), then pval = -pval
    // ----- Transform data
    arma::vec z = y - x*theta;
    int Tsize = z.n_elem;
    arma::mat params = as<arma::mat>(approxDistDL(Tsize, simdist_N));
    // ----- Compute test stats
    arma::vec eps = z - mean(z);
    arma::vec Fx = calc_DLmcstat(eps, N, params, pval_type);
    // ----- Obtain p-value
    double F0 = Fx(N);
    arma::vec FN = Fx.subvec(0,N-1);
    pval = -MCpval(F0, FN, "geq");
  }
  return(pval);
}

// ==============================================================================
//' @title Dufour & Luger (2017) moment-based MMC test p-value function to be maximized
//'
//' 
//' @references Dufour, J. M., & Luger, R. (2017). Identification-robust moment-based 
//' tests for Markov switching in autoregressive models. Econometric Reviews, 36(6-9), 713-727.
//' 
//' @export
// [[Rcpp::export]]
double DLMMCpval_fun_max(arma::vec theta, arma::vec y, arma::mat x, int N, int simdist_N, Rcpp::String pval_type, bool stationary_ind, double lambda){
  double pval = -DLMMCpval_fun(theta, y, x, N, simdist_N, pval_type, stationary_ind, lambda);
  return(pval);
}