#include <RcppArmadillo.h>
#include "methods.h"
//[[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;

//  ==============================================================================
//' @title Test statistic for switch in mean and variance
//'
//' @description This function computes part of the test statistic given by 
//' eq. 2.5 of CHP 2014 when the alternative has switching mean and variance. 
//' The output is used in \code{\link{chpStat}} which computes the full test
//' statistics.
//'
//' @param mdl List containing model attributes (see \code{\link{ARmdl}}).
//' @param rho Number determining value of \code{rho}.
//' @param ltmt List containing derivatives output from \code{\link{chpDmat}}.
//' @param hv Number determining value of \code{h}.
//' 
//' @return Part of test statistic given \code{rho} and \code{hv} value. 
//' 
//' @keywords internal
//' 
//' @references Carrasco, Marine, Liang Hu, and Werner Ploberger. 2014. “Optimal 
//' test for Markov switching parameters.” \emph{Econometrica} 82 (2): 765–784.
//' 
//' @export
// [[Rcpp::export]]
arma::vec calc_mu2t_mv(List mdl, double rho, List ltmt, arma::vec hv){
  // calc mu2t (eq. 2.5) for test with switch in mean and variance 
  int nn              = mdl["n"];
  int nar             = mdl["p"];
  double mu0          = mdl["mu"];
  arma::mat xtj       = mdl["x"];
  arma::mat ltmx      = ltmt["ltmx"];
  double mtmu         = ltmt["mtmu"];
  arma::vec mtmusig2  = ltmt["mtmusig2"];
  arma::mat mtmuphi   = ltmt["mtmuphi"];
  arma::mat mtphisig2 = ltmt["mtphisig2"];
  arma::vec mtsig2    = ltmt["mtsig2"];
  arma::vec mu2t(nn,arma::fill::zeros);
  arma::mat xs(nn,nar+2,arma::fill::zeros);
  xs.row(1)           = rho * ltmx.row(0);
  mu2t.row(1)         = arma::trans(hv) * arma::trans(ltmx.row(1)) * xs.row(1) * hv;
  for (int xi=0; xi<nn; xi++){
    arma::mat m(nar+2,nar+2,arma::fill::zeros);
    arma::mat mtphi = arma::trans(-(xtj.row(0)-mu0)) * (xtj.row(0)-mu0);
    
    int lst = nar + 1;
    m.submat(0,0,0,0) = mtmu;
    m.submat(0,1,0,nar) = mtmuphi.row(xi);
    m.submat(0,lst,0,lst) = mtmusig2(xi);
    
    m.submat(1,0,nar,0) = arma::trans(mtmuphi.row(xi));
    m.submat(1,1,nar,nar) = mtphi;
    m.submat(1,lst,nar,lst) = arma::trans(mtphisig2.row(xi));
    
    m.submat(lst,0,lst,0) = mtmusig2(xi);
    m.submat(lst,1,lst,nar) = mtphisig2.row(xi);
    m.submat(lst,lst,lst,lst) = mtsig2(xi);
    
    arma::mat lx = arma::trans(ltmx.row(xi)) * ltmx.row(xi);
    mu2t(xi) = arma::as_scalar(mu2t(xi) + arma::trans(hv) * (m + lx) * hv) / 2;
    
    if (xi > 2){
      xs.row(xi) = rho * (xs.row(xi- 1) + ltmx.row(xi- 1));
      mu2t(xi) = mu2t(xi) + arma::as_scalar(arma::trans(hv) * arma::trans(ltmx.row(xi)) * xs.row(xi) * hv);
    }
  }
  return(mu2t);
}

//  ==============================================================================
//' @title Test statistic for switch in mean only 
//'
//' @description This function computes part of the test statistic given by 
//' eq. 2.5 of CHP 2014 when the alternative has switching mean only. The output 
//' is used in \code{\link{chpStat}} which computes the full test statistics.
//'
//' @param mdl List containing model attributes (see \code{\link{ARmdl}}).
//' @param rho Number determining value of \code{rho}.
//' @param ltmt List containing derivatives output from \code{\link{chpDmat}}.
//' 
//' @return Part of test statistic given \code{rho} and \code{hv} value. 
//' 
//' @keywords internal
//' 
//' @references Carrasco, Marine, Liang Hu, and Werner Ploberger. 2014. “Optimal test for Markov switching parameters.” \emph{Econometrica} 82 (2): 765–784.
//' 
//' @export
// [[Rcpp::export]]
arma::vec calc_mu2t(List mdl, double rho, List ltmt){
  int nn              = mdl["n"];
  arma::mat ltmx      = ltmt["ltmx"];
  double mtmu         = ltmt["mtmu"];
  arma::vec ltmu      = ltmx.col(0);
  arma::vec mu2t(nn,arma::fill::zeros);
  arma::vec xs(nn,arma::fill::zeros);
  xs(1)               = rho * ltmu(0);
  mu2t(1)             = ltmu(1) * xs(1);
  for (int xi=2; xi<nn; xi++){
    xs(xi) = rho * (xs(xi- 1) + ltmu(xi- 1));
    mu2t(xi) = ltmu(xi) * xs(xi);
  }
  mu2t = (mtmu+pow(ltmu,2))/2+mu2t;
  return(mu2t);
}

//  ==============================================================================
//' @title Test statistic for CHP 2014 parameter stability test
//' 
//' @description This function computes the supTS and expTS test-statistics 
//' proposed in CHP 2014.
//'
//' @param mdl List containing model attributes (see \code{\link{ARmdl}}).
//' @param rho_b Number determining bounds for distribution of \code{rh0} (i.e. \code{rho} ~ \code{[-rho_b,rho_b]}).
//' @param ltmt List containing derivatives output from \code{\link{chpDmat}}.
//' @param msvar Boolean indicator. If \code{TRUE}, there is a switch in variance. If \code{FALSE} only switch in mean is considered.
//' 
//' @return A (\code{2 x 1}) vector with supTS test statistic as first element and expTS test-statistics as second element.
//' 
//' @keywords internal
//' 
//' @references Carrasco, Marine, Liang Hu, and Werner Ploberger. 2014. “Optimal 
//' test for Markov switching parameters.” \emph{Econometrica} 82 (2): 765–784.
//' 
//' @export
// [[Rcpp::export]]
arma::vec chpStat(List mdl, double rho_b, List ltmt,bool msvar){
  int nn              = mdl["n"];
  int nar             = mdl["p"];
  int tt              = nn + nar;
  arma::mat ltmx      = ltmt["ltmx"];
  double seqlen       = (rho_b-(-rho_b))*100 + 1;
  arma::vec chps(2,arma::fill::zeros);
  if (msvar==0){
    // If only Mean can Switch
    arma::vec cv(seqlen,1,arma::fill::zeros); // stores supTS critical values for each rho
    arma::vec cv2(seqlen,1,arma::fill::zeros); // stores expTS critical values for each rho
    arma::vec rhotmp = arma::linspace<arma::vec>(-rho_b,rho_b,seqlen);
    for (int ir=0; ir<seqlen; ir++){
      double rhotmp2 = arma::as_scalar(rhotmp(ir));
      arma::vec mu2t = calc_mu2t(mdl,rhotmp2,ltmt);
      double gamma_e = sum(mu2t)/sqrt(tt);
      // error from mu2t - projection of mu2t on lt1
      arma::mat tmp = arma::trans(ltmx) * ltmx;
      arma::vec epsilont = mu2t - (ltmx * arma::inv(tmp) * arma::trans(ltmx) * mu2t);
      double esqe = mean(arma::square(epsilont));
      double tspe = gamma_e/sqrt(esqe);
      double tol = 0.00001; 
      if (esqe<tol){
        cv(ir) = 0;
        cv2(ir) = 1;
      }else {
        // supTS test statistic 
        arma::vec tspetmp(2,arma::fill::zeros);
        tspetmp(0) = tspe;
        cv(ir) = pow(max(tspetmp),2)/2;
        // expTS test statistic  
        double tspetmp2 = tspe - 1;
        cv2(ir) = sqrt(2*arma::datum::pi)*exp(pow(tspetmp2,2)/2)*arma::normcdf(tspetmp2);
      }
      chps(0) = arma::as_scalar(max(cv));
      chps(1) = arma::as_scalar(mean(cv2));
    }
  }else{
    // If both Mean and Variance can Switch
    arma::vec mtmusig2  = ltmt["mtmusig2"];
    arma::mat mtmuphi   = ltmt["mtmuphi"];
    arma::mat mtphisig2 = ltmt["mtphisig2"];
    arma::vec mtsig2    = ltmt["mtsig2"];
    arma::mat cv(seqlen,100,arma::fill::zeros); // stores supTS critical values for each rho
    arma::mat cv2(seqlen,100,arma::fill::zeros); // stores expTS critical values for each rho
    arma::mat h(100,2,arma::fill::randn);
    // loop  over h 
    for (int ih=0; ih<100; ih++){
      arma::vec hv((nar + 2),arma::fill::zeros);
      arma::mat h2 = h.row(ih)%h.row(ih);
      arma::mat hu = h.row(ih) / h2;   // h-vector uniformly over the unit sphere
      hv(0) = hu(0,0);
      hv(nar + 1) = hu(0,1);
      arma::vec rhotmp = arma::linspace<arma::vec>(-rho_b,rho_b,seqlen);
      for (int ir=0; ir<seqlen; ir++){
        double rhotmp2 = arma::as_scalar(rhotmp(ir));
        arma::vec mu2t = calc_mu2t_mv(mdl,rhotmp2,ltmt,hv);
        double gamma_e = sum(mu2t)/sqrt(tt);
        // error from mu2t - projection of mu2t on lt1
        arma::mat tmp = arma::trans(ltmx) * ltmx;
        arma::vec epsilont = mu2t - (ltmx * arma::inv(tmp) * arma::trans(ltmx) * mu2t);
        double esqe = mean(arma::square(epsilont));
        double tspe = gamma_e/sqrt(esqe);
        double tol = 0.00001; 
        if (esqe<tol){
          cv(ir,ih) = 0;
          cv2(ir,ih) = 1;
        }else {
          // supTS test statistic 
          arma::vec tspetmp(2,arma::fill::zeros);
          tspetmp(0) = tspe;
          cv(ir,ih) = pow(max(tspetmp),2)/2;
          // expTS test statistic  
          double tspetmp2 = tspe - 1;
          cv2(ir,ih) = sqrt(2*arma::datum::pi)*exp(pow(tspetmp2,2)/2)*arma::normcdf(tspetmp2);
        }
      }
    }
    arma::vec cvtmp =  arma::vectorise(cv);
    arma::vec cv2tmp =  arma::vectorise(cv2);
    chps(0) = arma::as_scalar(max(cvtmp));
    chps(1) = arma::as_scalar(mean(cv2tmp));
  }
  return(chps);
}

//  ==============================================================================
//' @title Bootstrap critical values for CHP 2014 parameter stability test
//'
//' @description This bootstrap procedure is described on pg. 771 of CHP 2014.
//'
//'
//' @param mdl List containing model attributes (see \code{\link{ARmdl}}).
//' @param rho_b Number determining bounds for distribution of \code{rh0} (i.e. \code{rho} ~ \code{[-rho_b,rho_b]}).
//' @param N Number of bootstrap simulations.
//' @param msvar Boolean indicator. If \code{TRUE}, there is a switch in variance. If \code{FALSE} only switch in mean is considered.
//' 
//' @return Bootstrap critical values
//' 
//' @references Carrasco, Marine, Liang Hu, and Werner Ploberger. 2014. “Optimal 
//' test for Markov switching parameters.” \emph{Econometrica} 82 (2): 765–784.
//' 
//' @export
// [[Rcpp::export]]
arma::mat CHPbootCV(List mdl, double rho_b, int N, bool msvar){
  // calling required R functions 
  Function simuAR("simuAR");
  Function ARmdl("ARmdl");
  Function chpDmat("chpDmat");
  // define vars from Model list
  int ar =  mdl["p"];
  arma::vec supb(N,arma::fill::zeros);  // stores the bootstrapped supTS critical value
  arma::vec expb(N,arma::fill::zeros);  // stores the bootstrapped expTS critical value
  List y_out_tmp;
  arma::vec y0;
  List Mdl_tmp;
  List ltmtb;
  arma::vec cv4;
  List armdl_con;
  bool const_con = TRUE;
  bool getSE_con = FALSE;
  armdl_con["const"] = const_con;
  armdl_con["getSE"] = getSE_con;
  int itb = 0;
  while (itb<N){
    // simulate the series N times according to ML estimators 
    y_out_tmp = simuAR(mdl);
    y0        = as<arma::vec>(y_out_tmp["y"]);
    Mdl_tmp   = ARmdl(y0,ar,armdl_con);
    ltmtb     = chpDmat(Mdl_tmp, msvar);
    cv4       = chpStat(Mdl_tmp, rho_b, ltmtb, msvar);
    if (cv4.has_nan()==FALSE){
      supb(itb)   = cv4(0);
      expb(itb)   = cv4(1);
      itb = itb + 1;
    }
  }
  arma::mat boot_out = join_rows(supb,expb);
  return(boot_out);  
}










