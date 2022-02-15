#ifndef METHODS_H
#define METHODS_H

#include <RcppArmadillo.h>
using namespace Rcpp;

arma::vec sumfinite(arma::mat x, int ncol = 1);

arma::mat finitemat(arma::mat x);

arma::vec limP(arma::mat P, int k);

List paramList(arma::vec theta, int ar, int k, bool msmu, bool msvar);

List VARparamList(arma::vec theta, int N, int ar, int k, bool msmu, bool msvar);

arma::mat calcMSResid(List mdl, arma::mat mu, int k);

List calcMSVARResid(List mdl, List mu, int k);

List ts_lagged(arma::mat Y, int ar);

arma::mat randTransMat(int k, int n = 200);

arma::vec initVals(List mdl, int k, bool msmu, bool msvar);

arma::vec initVals2(arma::vec theta, int k, bool msmu, bool msvar, arma::vec Y);

arma::vec initValsVAR(arma::vec mu, arma::mat sigma, int k, bool msmu, bool msvar);

List initValsKM(arma::vec Y, int k, bool msmu, bool msvar);

double MCpval(double test_stat, arma::vec null_vec, Rcpp::String type = "geq");

List simuAR(List mdl_h0, int burnin = 200);

List simuMSAR(List mdl_h0, Rcpp::String type = "markov", int burnin = 200);

List simuMS(List mdl_h0, Rcpp::String type = "markov", int burnin = 200);

#endif