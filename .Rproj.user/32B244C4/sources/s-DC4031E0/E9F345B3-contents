#include <RcppArmadillo.h>
//[[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;
// ------------------------------------------------------------------------------------------
//' @title Calculate Moment-based test-statistics 
//'
//' This function computes the four momment-based test-statistics (eq. 11 - 14 in paper) 
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
arma::vec calc_moments(arma::vec ehat){
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
// ------------------------------------------------------------------------------------------
//' @title Test-statistics from simulated data
//'
//' This function computes simulates residuals from a standard normal distribution
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
arma::mat sim_moments(int t,int totsim){
  /* Pre-allocate matrix */
  arma::mat stats(totsim,4, arma::fill::zeros);
  /* loop and perform simultion each time */
  for (int isim = 0; isim <= totsim-1; isim++){
    arma::vec esim = arma::randn<arma::vec>(t);
    stats.row(isim) = arma::trans(calc_moments(esim - mean(esim)));
  }
  return(stats);
}
// ------------------------------------------------------------------------------------------
//' @title Combine p-values 
//'
//' This function is used to combine the p-values as in eq 17 and 18 of Dufour & Luger (2017).
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
arma::vec combine_stat(arma::vec s0,arma::mat sN, arma::mat params, int N, std::string type){
  arma::mat G0(1,4,arma::fill::zeros);
  arma::mat Gx(N-1,4,arma::fill::zeros);
  for (int im = 0; im<4; im++){
    G0(0,im) = 1 - (exp(params(0,im)+params(1,im)*s0(im))/(1+exp(params(0,im)+params(1,im)*s0(im))));
    Gx.col(im) = 1 - (exp(params(0,im)+params(1,im)*sN.col(im))/(1+exp(params(0,im)+params(1,im)*sN.col(im))));
  }
  arma::vec Fx(N,arma::fill::zeros);
  double F0 = 0;
  if (type=="min"){
    // MC p-value using Tippett (1931) &  Wilkinson (1951) combination method. 
    F0 = 1 - min(G0.row(0));
    for (int isim = 0; isim<N-1; isim++){
      Fx(isim) = 1 - min(Gx.row(isim));
    }
  }
  if (type=="prod"){
    // MC p-value using Tippett (1931) &  Wilkinson (1951) combination method. 
    F0 = 1 - (G0(0)*G0(1)*G0(2)*G0(3));
    for (int isim = 0; isim<N-1; isim++){
      Fx(isim) = 1 - (Gx(isim,0)*Gx(isim,1)*Gx(isim,2)*Gx(isim,3));
    }
  }
  Fx(N-1) = F0;
  return(Fx);
}
// ------------------------------------------------------------------------------------------
//' @title Calculate combined p-value test statistic 
//'
//' This function computes test-statistics for observed and simulated samples. 
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
arma::mat calc_mcstat(arma::vec ezt, int N, arma::mat params){
  // Get length of series 
  int Tsize = ezt.n_elem;
  // calulate moments of data 
  arma::vec S0 = calc_moments(ezt);
  // Simulated data 
  arma::mat SN = sim_moments(Tsize,N-1); // must leaves as N-1.
  // Get individual moment p-values 
  arma::vec Fprodx = combine_stat(S0, SN, params, N, "prod");
  arma::vec Fminx  = combine_stat(S0, SN, params, N, "min");  
  arma::mat Fvals = join_rows(Fminx,Fprodx);
  return(Fvals);
}
// ------------------------------------------------------------------------------------------
//' @title Loop for \emph{approxDist}
//'
//' This function performs the loop in \emph{approxDist}. It is written in C++ for better 
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
// ------------------------------------------------------------------------------------------
//' @title Evaluate Grid for GridSearch optimization. 
//'
//' This function performs the gridsearch loop over nuisance parameter space in 
//' \emph{calc_mmcpval}. It is written in C++ for better  performance in terms of speed.
//' It is only used when calling gridSearch_paramCI or randSearch_paramCI as optimization 
//' method.
//' 
//' @param ext observed series
//' @param N total number of test statistics (i.e. simulated + observed = N)
//' 
//' @return maximized p-value over parameter space.
//' 
//' @references Dufour, J. M., & Luger, R. (2017). Identification-robust moment-based 
//' tests for Markov switching in autoregressive models. Econometric Reviews, 36(6-9), 713-727.
//' 
//' @export
// [[Rcpp::export]]  
arma::mat grid_eval(arma::vec y, arma::mat x,arma::mat mcgrid,int N,arma::mat params){
  int totpoints=mcgrid.col(0).n_elem;
  arma::mat outputmin(totpoints,5,arma::fill::zeros);
  arma::mat outputprod(totpoints,5,arma::fill::zeros);
  for (int iv = 0; iv<totpoints; iv++){
    arma::vec z = y - x * arma::trans(mcgrid.row(iv));
    // Saved as: pval [1], test-stat[2], crit-val[3,4,5]
    //arma::mat outputtmp(N,2,arma::fill::zeros);
    arma::mat outputtmp = calc_mcstat((z-mean(z)),N,params);
    // find Ranks 
    int nlen = outputtmp.col(0).n_elem;
    arma::vec u1 = arma::randu<arma::vec>(nlen);
    arma::vec u2 = arma::randu<arma::vec>(nlen);
    arma::vec Fmin = outputtmp.col(0);
    arma::vec Fprod = outputtmp.col(1);
    double Fmin0 = Fmin(nlen - 1);
    arma::vec Fminx = arma::sort(Fmin.subvec(0,nlen - 2));
    double Fprod0 = Fprod(nlen - 1);
    arma::vec Fprodx = arma::sort(Fprod.subvec(0,nlen - 2));
    double RankFmin=1+sum(Fmin0>Fminx)+sum(u1(N-1)>=u1.elem(find(Fminx==Fmin0)));
    double RankFprod=1+sum(Fprod0>Fprodx)+sum(u2(N-1)>=u2.elem(find(Fprodx==Fprod0)));
    // Save pvals 
    outputmin(iv,0) = (N+1-RankFmin)/N;
    outputprod(iv,0) = (N+1-RankFprod)/N;
    // Save test stat 
    outputmin(iv,1) = Fmin0;
    outputprod(iv,1) = Fprod0;
    // Save crit-vals
    outputmin(iv,2) = Fminx(round(0.90*N) - 1);
    outputmin(iv,3) = Fminx(round(0.95*N) - 1);
    outputmin(iv,4) = Fminx(round(0.99*N) - 1);
    outputprod(iv,2) = Fprodx(round(0.90*N) - 1);
    outputprod(iv,3) = Fprodx(round(0.95*N) - 1);
    outputprod(iv,4) = Fprodx(round(0.99*N) - 1);
  }
  // find rows with max pvals.
  arma::mat outputminmaxtmp = outputmin.rows(find(outputmin.col(0)==max(outputmin.col(0))));
  // join with respective params from mcgrid. 
  arma::mat outputminparam = mcgrid.rows(find(outputmin.col(0)==max(outputmin.col(0))));
  arma::mat outputminmax = join_rows(outputminmaxtmp,outputminparam);
  // Find max test-stat among max pvals in case there are more than one set of 
  // params which give the same pvals.
  arma::mat outputminmaxstat = outputminmax.rows(find(outputminmax.col(1)==max(outputminmax.col(1))));
  // Repeat for prod type test- stat 
  arma::mat outputprodmaxtmp = outputprod.rows(find(outputprod.col(0)==max(outputprod.col(0))));
  arma::mat outputprodparam = mcgrid.rows(find(outputprod.col(0)==max(outputprod.col(0))));
  arma::mat outputprodmax = join_rows(outputprodmaxtmp,outputprodparam);
  arma::mat outputprodmaxstat = outputprodmax.rows(find(outputprodmax.col(1)==max(outputprodmax.col(1))));
  arma::mat output = join_cols(outputminmaxstat,outputprodmaxstat);
  return(output);
}

