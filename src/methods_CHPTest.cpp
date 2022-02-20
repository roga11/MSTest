#include <RcppArmadillo.h>
//[[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;
// ------------------------------------------------------------------------------------------
//' @title Calc eq. 2.5 from CHP (2014) where both mean and variance can switch
//'
//' @description When alternative has bith switching mean and variance, we take the second 
//' derivative of the log likelihood function w.r.t mu, phi and sigma.
//' 
//' Output from this function is used as input in \emph{chpStat}
//'
//' @param mdl List containing model information
//' @param rho_b bound for rho (nuisance param space)
//' @param ltmt List conatining derivatives (i.e. is the output when using \emph{chpDmat})
//' 
//' @return mu_2t from eq. 2.5 and used in test-statistic caluclation
//' 
//' @references Carrasco, Marine, Liang Hu, and Werner Ploberger. 2014. 
//' “Optimal test for Markov switch- ing parameters.” \emph{Econometrica} \bold{82 (2)}: 765–784.
//' 
//' @export
// [[Rcpp::export]]
arma::vec calc_mu2t_mv(List mdl, double rho, List ltmt, arma::vec hv){
  // calc mu2t (eq. 2.5) for test with switch in mean and variance 
  int nn              = mdl["n"];
  int nar             = mdl["ar"];
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
// ------------------------------------------------------------------------------------------
//' @title Calc eq. 2.5 from CHP (2014) where only mean can switch
//'
//' @description When alternative only has Switching mean (and not variance), we only take the second 
//' derivative of the log likelihood function w.r.t mu and not phi or sigma.
//' 
//' Output from this function is used as input in \emph{chpStat}
//'
//' @param mdl List containing model information
//' @param rho_b bound for rho (nuisance param space)
//' @param ltmt List conatining derivatives (i.e. is the output when using \emph{chpDmat})
//' 
//' @return mu_2t from eq. 2.5 and used in test-statistic caluclation
//' 
//' @references Carrasco, Marine, Liang Hu, and Werner Ploberger. 2014. 
//' “Optimal test for Markov switch- ing parameters.” \emph{Econometrica} \bold{82 (2)}: 765–784.
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
// ------------------------------------------------------------------------------------------
//' @title CHP Test Statistic
//' 
//' @description Calculate supTS and expTS test-statistics from CHP (2014).
//'
//' @param mdl List containing model information
//' @param rho_b bound for rho (nuisance param space)
//' @param ltmt List containing  relevant first and second derivatives of log likelihood function.
//' @param var_switch variance switch indicator
//' 
//' @return Test Statistic
//' 
//' @references Carrasco, Marine, Liang Hu, and Werner Ploberger. 2014. 
//' “Optimal test for Markov switch- ing parameters.” \emph{Econometrica} \bold{82 (2)}: 765–784.
//' 
//' @export
// [[Rcpp::export]]
arma::vec chpStat(List mdl, double rho_b, List ltmt,int var_switch){
  int nn              = mdl["n"];
  int nar             = mdl["ar"];
  int tt              = nn + nar;
  arma::mat ltmx      = ltmt["ltmx"];
  double seqlen       = (rho_b-(-rho_b))*100 + 1;
  arma::vec chps(2,arma::fill::zeros);
  if (var_switch==0){
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

// ------------------------------------------------------------------------------------------
//' @title Bootstrap Critival Values CHP Test
//'
//' @description This bootstrap procedure is described on page 771 of CHP (2014)
//'
//' @param mdl List containing model information
//' @param rho_b bound for rho (nuisance param space)
//' @param N number of simulations
//' @param var_switch variance switch indicator
//' 
//' @return Bootstrap critical values
//' 
//' @references Carrasco, Marine, Liang Hu, and Werner Ploberger. 2014. 
//' “Optimal test for Markov switch- ing parameters.” \emph{Econometrica} \bold{82 (2)}: 765–784.
//' 
//' @export
// [[Rcpp::export]]
arma::mat bootCV(List mdl,double rho_b, int N, int var_switch){
  // calling required R functions 
  Function simuAR("simuAR");
  Function ARmdl("ARmdl");
  Function chpDmat("chpDmat");
  // define vars from Model list
  int ar =  mdl["ar"];
  arma::vec supb(N,arma::fill::zeros);  // stores the bootstrapped supTS critical value
  arma::vec expb(N,arma::fill::zeros);  // stores the bootstrapped expTS critical value
  for (int itb=0; itb<N; itb++){
    // first simulate the series N (e.g., 3000) times according to ML estimators 
    List y_out_tmp = simuAR(mdl);
    arma::vec y0 = y_out_tmp["y"];
    List  Mdl_tmp  = ARmdl(y0,ar);
    List ltmtb  = chpDmat(Mdl_tmp,var_switch);
    arma::vec cv4  = chpStat(Mdl_tmp, rho_b, ltmtb, var_switch);
    supb(itb)   = cv4(0);
    expb(itb)   = cv4(1);
  }
  arma::mat boot_out = join_rows(sort(supb),sort(expb));
  return(boot_out);  
}
