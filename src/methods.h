#ifndef METHODS_H
#define METHODS_H

#include <RcppArmadillo.h>
using namespace Rcpp;

arma::mat cov2corr(arma::mat cov_mat);

arma::vec covar_vech(arma::mat mat);

arma::mat covar_unvech(arma::vec sig, int q);

double logLike_AR(arma::vec theta, List mdl);

double logLike_VAR(arma::vec theta, List mdl);

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

List simuAR(List mdl_h0, int burnin = 100);

List simuMSAR(List mdl_h0, int burnin = 100);

List simuVAR(List mdl_h0, int burnin = 100);

List simuMSVAR(List mdl_h0, int burnin = 100);

List simuNorm(List mdl_h0);

List simuHMM(List mdl_h0, int burnin = 100);

double MSloglik_fun(arma::vec theta, List mdl, int k);

double MSloglik_fun_min(arma::vec theta, List mdl, int k);

List MSloglik(arma::vec theta, List mdl, int k);

double MSVARloglik_fun(arma::vec theta, List mdl, int k);

double MSVARloglik_fun_min(arma::vec theta, List mdl, int k);

List MSVARloglik(arma::vec theta, List mdl, int k);

List MS_EMaximization(arma::vec theta, List mdl, List MSloglik_output, int k);

List MSVAR_EMaximization(arma::vec theta, List mdl, List MSloglik_output, int k);

List MS_EMiter(List mdl, List EMest_output, int k);

List MSVAR_EMiter(List mdl, List EMest_output, int k);

List MS_EMest(arma::vec theta_0, List mdl, int k, List optim_options);

List MSVAR_EMest(arma::vec theta_0, List mdl, int k, List optim_options);

#endif