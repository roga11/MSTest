#ifndef METHODS_H
#define METHODS_H

#include <RcppArmadillo.h>
using namespace Rcpp;

arma::mat cov2corr(arma::mat cov_mat);

arma::vec covar_vech(arma::mat mat);

arma::mat covar_unvech(arma::vec sig, int q);

arma::mat randP(int k);

arma::vec limP(arma::mat P);

List ts_lagged(arma::mat Y, int ar);

List paramListMS(arma::vec theta, int ar, int k, bool msmu, bool msvar);

List paramListMSVAR(arma::vec theta, int q, int ar, int k, bool msmu, bool msvar);

arma::mat calcMSResid(List mdl, arma::mat mu, int k);

List calcMSVARResid(List mdl, List mu, int k);

arma::vec initValsMS(List mdl, int k);

arma::vec initValsMSVAR(List mdl, int k);

double MCpval(double test_stat, arma::vec null_vec, Rcpp::String type = "geq");

List simuAR(List mdl_h0, int burnin = 200);

List simuMS(List mdl_h0, int burnin = 200);

List simuVAR(List mdl_h0, int burnin = 200);

List simuMSVAR(List mdl_h0, int burnin = 200);

#endif