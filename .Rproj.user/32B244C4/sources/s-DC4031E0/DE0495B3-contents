// --------------------------------------------------------------------------------------------- // 
// Title:            Likelihood Ratio Maximized Monet-Carlo test for AR MS-models - C++ functions
// Author(s):        Gabriel Rodriguez Rondon
// Initial version:  12/MAR/2020
// This version:     25/MAR/2020
// Designed for:     R/R-Studio
// --------------------------------------------------------------------------------------------- // 
// Purpose: 
// C++ functions used to improve speed of ikelihood Ratio Maximized Monet-Carlo test
//
// --------------------------------------------------------------------------------------------- // 
// Notes: 
// --------------------------------------------------------------------------------------------- // 
// ------------------------------------ Requirements ------------------------------------------- //
#include <RcppArmadillo.h>
#include <RcppArmadilloExtensions/sample.h>
//[[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;
// ---------------------------------- Generic Functions ---------------------------------------- //

// ---------------- Genrate Random Initial Transition Matrix
// [[Rcpp::export]]
arma::mat transMat_cpp(arma::vec ind, int k){
  double n = ind.n_elem;
  arma::mat transmat(k,k,arma::fill::zeros);
  arma::vec tmp = diff(ind);
  arma::vec ind2 = ind.rows(1,n-1);
  arma::vec ind3 = ind.rows(0,n-2);
  for(int xi=1; xi<=k;xi++){
    double n2 = sum(ind3==xi);
    for(int xxi=1; xxi<=k; xxi++){
      transmat(xi-1,xxi-1) = sum(ind2==xxi && tmp==xxi-xi)/n2;
    }
  }
  return(transmat);
}

// ---------------- Generate Random Initial Values
// [[Rcpp::export]]
List randInitVal_cpp(List mdl,int k){
  arma::vec y = mdl["y"];
  arma::mat x = mdl["X"];
  int n = y.n_elem;
  int ar = mdl["ar"];
  int npar = ar + 1;
  arma::vec seqk = arma::linspace(1,k,k);
  arma::vec ind = RcppArmadillo::sample(seqk,n,TRUE);
  arma::mat MScoef(npar,k,arma::fill::zeros);
  arma::mat MSstd(1,k,arma::fill::zeros);
  for(int xi=0; xi<k; xi++){
    arma::mat ytmp = y.elem(find(ind==xi+1));
    arma::mat xtmp = x.rows(find(ind==xi+1));
    double ntmp = ytmp.n_elem;
    if (ar>0){
      MScoef.col(xi) = inv(arma::trans(x)*x)*arma::trans(x)*y;
      arma::vec int_tmp = inv(arma::trans(xtmp)*xtmp)*arma::trans(xtmp)*ytmp;
      MScoef(0,xi) = int_tmp(0);
    }else{
      MScoef.col(xi) = inv(arma::trans(xtmp)* xtmp )*arma::trans(xtmp)*ytmp;
    }
    arma::mat condmean =  xtmp * MScoef.col(xi);
    arma::mat residtmp = ytmp - condmean;
    MSstd.col(xi)= sqrt((arma::trans(residtmp)*residtmp)/ntmp);
  }
  arma::mat transmat = transMat_cpp(ind,k);
  arma::mat  th(3*k,1,arma::fill::zeros);
  for(int xi=0; xi<k; xi++){
    th(xi,0) = MScoef(0,xi);
    th(k+xi,0) = transmat(xi,xi);
    th(2*k+xi,0) = MSstd(0,xi);
  }
  List out;
  out["th"] = th;
  out["MScoef"] = arma::trans(MScoef);
  return(out);
}


// ------------------------- Hamilton Estimation Method Functions -------------------------- //

// ---------------- Calculate Log-likelihood for Markov-Switching Model
// [[Rcpp::export]]
List ms_logLikel_cpp(arma::vec th,arma::vec y, arma::mat X,arma::mat MScoef){
  int n = y.n_elem;
  double p = th(2);         
  double q = th(3);   
  double pi = arma::datum::pi;
  arma::vec sig(2,arma::fill::zeros);
  sig(0) = pow(th(4),2);
  sig(1) = pow(th(5),2);
  double rho = (1-q)/(2-p-q);               
  arma::vec pa(2,arma::fill::zeros);
  pa(0) = rho;
  pa(1) = 1 - rho;
  arma::vec p1(4,arma::fill::zeros);
  arma::mat pax(n,4,arma::fill::zeros);
  arma::vec pfx(n,arma::fill::zeros);
  arma::mat yxx(n,2,arma::fill::zeros);
  for (int itt=0; itt<n; itt++){
    yxx.row(itt) =  arma::trans(1/sqrt(2*pi*(sig)))%exp((pow(y(itt) - X.row(itt)*arma::trans(MScoef),2)%arma::trans(-1/(2*sig))));
  }
  yxx = join_rows(join_rows(join_rows(p*yxx.col(0),(1-p)*yxx.col(1)),(1-q)*yxx.col(0)),q*yxx.col(1));
  double f = 0;
  for (int its=0; its<n; its++){
    p1(0) = pa(0)*yxx(its,0);
    p1(1) = pa(0)*yxx(its,1);
    p1(2) = pa(1)*yxx(its,2);
    p1(3) = pa(1)*yxx(its,3);
    pfx(its) = sum(p1);
    f = f + log(pfx(its));
    p1 = p1/pfx(its);
    pax.row(its) = arma::trans(p1); 
    pa(0) = p1(0) + p1(2);
    pa(1) = p1(1) + p1(3);
  }
  long double f0 = f;
  List ret;
  ret["f0"] = f0;
  ret["yxx"] = yxx;
  ret["pax"] = pax;
  ret["pfx"] = pfx;
  ret["rho"] = rho;
  return ret;
}

// ---------------- Expected Maximization (EM) step with Hamilton Smoother
// [[Rcpp::export]]
List EM_hamilton_cpp(arma::vec th,arma::vec y, arma::mat X,arma::mat MScoef,long double vof){
  int n = y.n_elem;
  // log like
  List ret = ms_logLikel_cpp(th,y,X,MScoef);
  long double f0 = ret["f0"];
  double rho = ret["rho"];
  arma::mat yxx = ret["yxx"];
  arma::mat pax = ret["pax"];
  arma::vec pfx = ret["pfx"];
  // smoother
  arma::mat qax(n,8,arma::fill::zeros);
  qax.row(0) = join_rows(pax.row(0),pax.row(0));
  qax(0,1) = 0;
  qax(0,3) = 0;
  qax(0,4) = 0;
  qax(0,6) = 0;
  for (int its=1; its<n; its++){
    arma::mat tmp1 = ((yxx(its,0)*qax.cols(0,3)) + (yxx(its,2)*qax.cols(4,7)));
    arma::mat tmp2 = ((yxx(its,1)*qax.cols(0,3)) + (yxx(its,3)*qax.cols(4,7)));
    qax = join_rows(tmp1,tmp2)/pfx(its,0);
    qax.row(its) = join_rows(pax.row(its),pax.row(its));
    qax(its,1) = 0;
    qax(its,3) = 0;
    qax(its,4) = 0;
    qax(its,6) = 0;
  }
  qax = qax.cols(0,3) + qax.cols(4,7);
  // Calculate filter and smoother probs that st=1 (col. 5) and st-1=1 (col. 6) 
  pax = join_rows(pax,join_rows(pax.col(0)+pax.col(2),pax.col(0)+pax.col(1)));
  qax = join_rows(qax,join_rows(qax.col(0)+qax.col(2),qax.col(0)+qax.col(1)));
  //------ @ Produce a new iteration on normal equations, if desired  @
  arma::mat p = sum(qax.submat(1,0,n-1,0))/(sum(qax.submat(1,5,n-1,5)) + rho - qax(0,4));
  arma::mat q = sum(qax.submat(1,3,n-1,3))/(sum(1-qax.submat(1,5,n-1,5)) + qax(0,4) - rho);
  // mean
  arma::mat mu(1,2,arma::fill::zeros);
  mu(0,0) = sum(y%qax.col(4))/(sum(qax.col(4)));
  mu(0,1) = sum(y%(1-qax.col(4)))/(sum(1-qax.col(4)));
  // mean vectorized 
  arma::mat CondMean = X*mu;
  // Standard deviation
  arma::mat sig = join_rows(pow(trans(pow((y-CondMean.col(0)),2))*(qax.col(4)/sum(qax.col(4))),0.5),
                            pow(trans(pow((y-CondMean.col(1)),2))*((1-qax.col(4))/sum(1-qax.col(4))),0.5));
  // residuals
  arma::mat resid = join_rows(y-CondMean.col(0),y-CondMean.col(1));
  //  group
  arma::mat thn = join_rows(mu,p,q,sig);
  arma::mat thl(1,3,arma::fill::zeros);
  thl(0,0) = f0;
  thl(0,1)  = f0 - vof;
  arma::mat deltah   = trans(thn)-th; 
  arma::vec tmpvar = max(abs(deltah));
  thl(0,2) = tmpvar(0);
  //use fabs() ?
  List out;
  out["thn"] = thn;
  out["MScoef"] = trans(mu);
  out["thl"] = thl;
  out["deltath"] = f0 - vof;
  out["f0"] = f0;
  out["resid"] = resid;
  return out;
}


// ---------------- Expected Maximization (EM) algorithm with Hamilton Smoother
// [[Rcpp::export]]
List msmdl_cpp(List mdl, int k,String method="eval",int mxit=500){
  arma::vec y = mdl["y"];
  arma::vec X = mdl["X"];
  List initVal = randInitVal_cpp(mdl,k);
  arma::vec th      = initVal["th"];
  arma::mat MScoef  = initVal["MScoef"];
  // initialize algorithm
  double vof = 0;
  List out = EM_hamilton_cpp(th,y,X,MScoef,vof);
  double thtol = 1.e-8;
  int it = 1;
  if (method=="eval"){
    // iterate to find EM
    double deltath = out["deltath"];
    while(it< mxit && fabs(deltath)>thtol){
      out = EM_hamilton_cpp(out["thn"],y,X,out["MScoef"],out["f0"]);
      deltath = out["deltath"];
      it++;
    }  
  }else if (method=="iter"){
    // iterate to find EM
    for (int itx=1; itx<mxit;itx++){
      out = EM_hamilton_cpp(out["thn"],y,X,out["MScoef"],out["f0"]);
      it++;
    }
    thtol = out["deltath"];
  }
  out["it"] = it;
  out["tol"] = thtol;
  out["y"] = y;
  out["X"] = X;
  out["x"] = mdl["x"];
  out["n"] = mdl["n"];
  out["initVal"] = initVal;
  return out;
}

// ---------------------------------------- END -------------------------------------------- //