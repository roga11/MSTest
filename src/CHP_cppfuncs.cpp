#include <RcppArmadillo.h>
//[[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;

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

// [[Rcpp::export]]
arma::vec calc_mu2t(List mdl, double rho, List ltmt){
  // calc mu2t (eq. 2.5) for test with switch in mean only
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

//' @title CHP Test Statistic
//'
//' @param mdl List containing model information
//' @param rho_b bound for rho (nuisance param space)
//' @param ltmt List containing  relevant first and second derivatives of log likelihood function.
//' @param var_switch variance switch indicator
//' 
//' @return Test Statistic
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


//' @title Bootstrap Critival Values CHP Test
//'
//' @param mdl List containing model information
//' @param rho_b bound for rho (nuisance param space)
//' @param N number of simulations
//' @param var_switch variance switch indicator
//' 
//' @return Bootstrap critical values
//' @export
// [[Rcpp::export]]
arma::mat bootCV(List mdl,double rho_b, int N, int var_switch){
  // calling required R functions 
  Function simu_AR_dgp("simu_AR_dgp");
  Function ARmdl("ARmdl");
  Function chpDmat("chpDmat");
  // define vars from Model list
  int nn = mdl["n"];
  int nar =  mdl["ar"];
  int tt = nn+nar;
  double mu = mdl["mu"];
  double v0 = mdl["stdev"];
  arma::vec phi = mdl["phi"];
  arma::vec supb(N,arma::fill::zeros);  // stores the bootstrapped supTS critical value
  arma::vec expb(N,arma::fill::zeros);  // stores the bootstrapped expTS critical value
  for (int itb=0; itb<N; itb++){
    // first simulate the series N (e.g., 3000) times according to ML estimators 
    NumericVector ys = simu_AR_dgp(tt,mu,v0,phi);
    List  Mdlb  = ARmdl(ys,nar);
    List ltmtb  = chpDmat(Mdlb,var_switch);
    arma::vec cv4  = chpStat(Mdlb, rho_b, ltmtb, var_switch);
    supb(itb)   = cv4(0);
    expb(itb)   = cv4(1);
  }
  //arma::mat boot_out = join_rows(supb,expb);
  arma::mat boot_out = join_rows(sort(supb),sort(expb));
  return(boot_out);  
}

// [[Rcpp::export]]
arma::mat CHPadist(int itn,int tr,double rho_b){
  arma::vec crt1(itn,arma::fill::zeros);  // stores the supTS statistic
  arma::vec crt2(itn,arma::fill::zeros);  // stores the expTS statistic
  arma::mat z(itn,tr + 1,arma::fill::randn);
  double seqlen = (rho_b-(-rho_b))*100 + 1;
  arma::vec rhotmp = arma::linspace<arma::vec>(-rho_b,rho_b,seqlen);
  // Bootstrap Iterations
  for ( int xxi=0; xxi<itn; xxi++){
    arma::vec cv(seqlen,arma::fill::zeros); // stores supTS asymptotic critical values for each rho
    arma::vec cv2(seqlen,arma::fill::zeros);// stores expTS asymptotic critical values for each rho
    // loop over differnt values of rho
    for (int xi=0; xi<seqlen; xi++){ 
      double rho = arma::as_scalar(rhotmp(xi));
      double sign = (rho >0)-(rho <0) ; // sign=1 if pos, =-1 if neg, =0 if null 
      double su = z(xxi,0); 
      for (int xr=0; xr<tr; xr++){
        su = su + pow(rho,xr + 1)*z(xxi,xr + 1);
      }
      double su2 = sign*su*sqrt(1-pow(rho,2));
      // calc stats 
      arma::vec tmp(2,arma::fill::zeros);
      tmp(0)  = su2;
      cv(xi)  = pow(max(tmp),2)/2; // supTS test statistic
      cv2(xi) = sqrt(2*arma::datum::pi)*exp( pow((su2 - 1),2)/2)*arma::normcdf(su2 - 1); // expTS test statistic
    }
    crt1(xxi) = arma::as_scalar(max(cv));
    crt2(xxi) = arma::as_scalar(mean(cv2));
  }
  arma::mat asymdist = join_rows(crt1,crt2);
  return(asymdist);
}