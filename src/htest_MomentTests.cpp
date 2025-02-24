#include <RcppArmadillo.h>
#include "methods.h"
//[[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;
// ==============================================================================
//' @title Moment-based test statistics 
//'
//' @description This function computes the four moment-based test statistics (eq. \code{11} - \code{14}) discussed in Dufour & Luger 2017.
//'
//' @param ehat A (\code{T x 1}) vector of restricted model residuals.
//' 
//' @return Vector containing the four test statistics.
//' 
//' @keywords internal
//' 
//' @references Dufour, J. M., & Luger, R. 2017. "Identification-robust moment-based 
//' tests for Markov switching in autoregressive models." \emph{Econometric Reviews}, 36(6-9), 713-727.
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
//' @title Simulated moment-based test statistics
//'
//' @description This function computes the four moment-based test statistics (eq. \code{11} - \code{14}) discussed in Dufour & Luger 2017 for \code{N} different simulated series.
//' 
//' @param Tsize Length of sample size for simulation.
//' @param N Number of simulated samples.
//' 
//' @return A (\code{N x 4}) matrix with \code{N} different simulated moment-based test statistics.
//' 
//' @keywords internal
//' 
//' @references Dufour, J. M., & Luger, R. 2017. "Identification-robust moment-based 
//' tests for Markov switching in autoregressive models." \emph{Econometric Reviews}, 36(6-9), 713-727.
//' 
//' @export
// [[Rcpp::export]]
arma::mat sim_DLmoments(int Tsize, int N){
  arma::mat stats(N, 4, arma::fill::zeros);
  arma::mat esim(Tsize, N, arma::fill::randn);
  // loop and perform simulation each time 
  for (int isim = 0; isim <N; isim++){
    stats.row(isim) = trans(calc_DLmoments(esim.col(isim) - mean(esim.col(isim))));
  }
  return(stats);
}


// ==============================================================================
//' @title Combine p-values 
//'
//' @description This function is used to combine the four moment-based p-values as in eq. \code{17} and \code{18} of Dufour & Luger 2017.
//' 
//' @param stats A (\code{l x 4}) matrix where \code{l} is the number of moment-based test statistics.
//' @param params A (\code{2 x 4}) matrix with parameters to combine test statistics. See \code{\link{approxDistDL}}.
//' @param type String determining the type of method used to combine p-values. If set to "min" the min method of combining p-values 
//' is used as in Fisher 1932 and Pearson 1933. If set to "prod" the product of p-values is used as in Tippett 1931 and Wilkinson 1951.
//' 
//' @return A (\code{N x 1}) vector with test statistics. The last element is the test statistic from observed data.
//' 
//' @keywords internal
//' 
//' @references Dufour, J. M., & Luger, R. 2017. "Identification-robust moment-based 
//' tests for Markov switching in autoregressive models." \emph{Econometric Reviews}, 36(6-9), 713-727.
//' @references Tippett, L. 1931. "The Method of Statistics". London: Williams & Norgate.
//' @references Wilkinson, B. 1951. "A statistical consideration in psychological research." \emph{Psychology Bulletin} 48:156–158.
//' @references Pearson, K. 1933. "On a method of determining whether a sample of size n
//'  supposed to have been drawn from a parent population having a known probability integral has probably been drawn at random". \emph{Biometrika} 25:379–410.
//' @references Fisher, R. 1932. "Statistical Methods for Research Workers." Edinburgh: Oliver and Boyd.
//' 
//' @export
// [[Rcpp::export]]
arma::vec combine_stat(arma::mat stats, arma::mat params, std::string type){
  int N = stats.n_rows;
  int nc = stats.n_cols;
  arma::mat Gx(N, nc, arma::fill::zeros);
  for (int im = 0; im<nc; im++){
    Gx.col(im) = 1 - (exp(params(0,im)+params(1,im)*stats.col(im))/(1+exp(params(0,im)+params(1,im)*stats.col(im))));
  }
  arma::vec Fx(N, arma::fill::zeros);
  if (type=="min"){
    // MC p-value using Tippett (1931) &  Wilkinson (1951) combination method. 
    for (int isim = 0; isim<N; isim++){
      Fx(isim) = 1 - min(Gx.row(isim));
    }
  }
  if (type=="prod"){
    // MC p-value using Fisher (1932) &  Pearson (1933) combination method. 
    for (int isim = 0; isim<N; isim++){
      Fx(isim) = 1 - prod(Gx.row(isim));
    }
  }
  return(Fx);
}

// ==============================================================================
//' @title Loop for \code{\link{approxDistDL}}
//'
//' @description This function performs the loop in required in \code{\link{approxDistDL}}. 
//' 
//' @param SN2 A (\code{T x 4}) matrix of  test-statistics.
//' 
//' @return The test statistics from simulated data. Used for NLS to get \code{params} needed to combine p-values.
//' 
//' @keywords internal
//' 
//' @references Dufour, J. M., & Luger, R. 2017. "Identification-robust moment-based 
//' tests for Markov switching in autoregressive models." \emph{Econometric Reviews}, 36(6-9), 713-727.
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
//' @title Moment-based MMC test p-value 
//'
//' @description This functions is used by numerical optimization algorithms for find maximum p-value given parameter vector \code{theta}.
//'
//' @param theta Value of nuisance parameters. Specifically, these are the consistent estimates of nuisance parameters as discussed in Dufour & Luger (2017) LMC procedure.
//' @param y series being tested.
//' @param x lagged values of series.
//' @param params A (\code{2 x 4}) matrix with parameters to combine test statistics. See \code{\link{approxDistDL}}.
//' @param sim_stats A (\code{N x 1}) vector with test statistics. The last element is the test statistic from observed data.
//' @param pval_type String determining the type of method used to combine p-values. If set to "min" the min method of combining p-values is used as in Fisher 1932 and Pearson 1933. If set to "prod" the product of p-values is used as in Tippett 1931 and Wilkinson 1951.
//' @param stationary_ind Boolean indicator determining if only stationary solutions should be considered if \code{TRUE} or any solution can be considered if \code{FALSE}. Default is \code{TRUE}.
//' @param lambda Numeric value for penalty on stationary constraint not being met. Default is \code{100}.
//' 
//' @return Maximized Monte Carlo p-value.
//' 
//' @keywords internal
//' 
//' @references Dufour, J. M., & Luger, R. 2017. "Identification-robust moment-based 
//' tests for Markov switching in autoregressive models." \emph{Econometric Reviews}, 36(6-9), 713-727.
//' 
//' @export
// [[Rcpp::export]]
double DLMMCpval_fun(arma::vec theta, arma::vec y, arma::mat x,
                     arma::mat params, arma::vec sim_stats,
                     Rcpp::String pval_type, bool stationary_ind, double lambda){
  bool stationary_constraint = FALSE;
  double F0 = 0;
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
    // ----- Compute test stats
    arma::vec eps = z - mean(z);
    arma::mat S0 = trans(calc_DLmoments(eps));
    // combine moment stats
    if (pval_type=="prod"){
      F0 = arma::as_scalar(combine_stat(S0, params, "prod"));
    }else if (pval_type=="min"){
      F0 = arma::as_scalar(combine_stat(S0, params, "min"));   
    }else{
      Rcerr << "type must be: 'min' or 'prod'.\n";
    }
    // ----- Get individual moment p-values 
    pval = MCpval(F0, sim_stats, "geq");
  }
  return(pval);
}

// ==============================================================================
//' @title Moment-based MMC test (negative) p-value 
//'
//' @description This functions is used by numerical optimization algorithms for find negative of maximum p-value given parameter vector \code{theta}.
//'
//' @param theta Value of nuisance parameters. Specifically, these are the consistent estimates of nuisance parameters as discussed in Dufour & Luger (2017) LMC procedure.
//' @param y series being tested.
//' @param x lagged values of series.
//' @param params A (\code{2 x 4}) matrix with parameters to combine test statistics. See \code{\link{approxDistDL}}.
//' @param sim_stats A (\code{N x 1}) vector with test statistics. The last element is the test statistic from observed data.
//' @param pval_type String determining the type of method used to combine p-values. If set to "min" the min method of combining p-values is used as in Fisher 1932 and Pearson 1933. If set to "prod" the product of p-values is used as in Tippett 1931 and Wilkinson 1951.
//' @param stationary_ind Boolean indicator determining if only stationary solutions should be considered if \code{TRUE} or any solution can be considered if \code{FALSE}. Default is \code{TRUE}.
//' @param lambda Numeric value for penalty on stationary constraint not being met. Default is \code{100}.
//'
//' @return Negative Maximized Monte Carlo p-value. 
//' 
//' @keywords internal
//' 
//' @references Dufour, J. M., & Luger, R. 2017. "Identification-robust moment-based 
//' tests for Markov switching in autoregressive models." \emph{Econometric Reviews}, 36(6-9), 713-727.
//' 
//' @export
// [[Rcpp::export]]
double DLMMCpval_fun_min(arma::vec theta, arma::vec y, arma::mat x, 
                         arma::mat params, arma::vec sim_stats,
                         Rcpp::String pval_type, bool stationary_ind, double lambda){
  double pval = -DLMMCpval_fun(theta, y, x, params, sim_stats, pval_type, stationary_ind, lambda);
  return(pval);
}