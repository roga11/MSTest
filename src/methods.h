#ifndef METHODS_H
#define METHODS_H

#include <RcppArmadillo.h>
using namespace Rcpp;

arma::mat cov2corr(arma::mat cov_mat);

arma::vec covar_vech(arma::mat mat);

arma::mat covar_unvech(arma::vec sig, int n);

arma::mat randP(int k);

arma::vec limP(arma::mat P);

List ts_lagged(arma::mat Y, int ar);

List paramList_MSARmdl(arma::vec theta, int p, int k, bool msmu, bool msvar);

List paramList_MSVARmdl(arma::vec theta, int q, int ar, int k, bool msmu, bool msvar);

arma::mat calcResid_MSARmdl(List mdl, arma::mat mu, int k);

List calcResid_MSVARmdl(List mdl, List mu, int k);

arma::vec initVals_HMmdl(List mdl, int k);

arma::vec initVals_MSARmdl(List mdl, int k);

arma::vec initVals_MSVARmdl(List mdl, int k);

double MCpval(double test_stat, arma::vec null_vec, Rcpp::String type = "geq");

arma::mat randSN(int n, int q);

List simuAR_cpp(List mdl_h0, int burnin = 100);

List simuARX_cpp(List mdl_h0, int burnin = 100);

List simuMSAR_cpp(List mdl_h0, int burnin = 100);

List simuMSARX_cpp(List mdl_h0, int burnin = 100);

List simuVAR_cpp(List mdl_h0, int burnin = 100);

List simuVARX_cpp(List mdl_h0, int burnin = 100);

List simuMSVAR_cpp(List mdl_h0, int burnin = 100);

List simuMSVARX_cpp(List mdl_h0, int burnin = 100);

List simuNorm_cpp(List mdl_h0, int burnin = 0);

List simuHMM_cpp(List mdl_h0, int burnin = 100);

double logLike_Nmdl(arma::vec theta, List mdl);

double logLike_ARmdl(arma::vec theta, List mdl);

double logLike_ARXmdl(arma::vec theta, List mdl);

double logLike_VARmdl(arma::vec theta, List mdl);

double logLike_VARXmdl(arma::vec theta, List mdl);

double logLike_HMmdl(arma::vec theta, List mdl, int k);

double logLike_HMmdl_min(arma::vec theta, List mdl, int k);

double logLike_MSARmdl(arma::vec theta, List mdl, int k);

double logLike_MSARXmdl(arma::vec theta, List mdl, int k);

double logLike_MSARmdl_min(arma::vec theta, List mdl, int k);

double logLike_MSARXmdl_min(arma::vec theta, List mdl, int k);

double logLike_MSVARmdl(arma::vec theta, List mdl, int k);

double logLike_MSVARXmdl(arma::vec theta, List mdl, int k);

double logLike_MSVARmdl_min(arma::vec theta, List mdl, int k);

double logLike_MSVARXmdl_min(arma::vec theta, List mdl, int k);

List ExpectationM_HMmdl(arma::vec theta, List mdl, int k);

List ExpectationM_MSARmdl(arma::vec theta, List mdl, int k);

List ExpectationM_MSVARmdl(arma::vec theta, List mdl, int k);

List EMaximization_HMmdl(arma::vec theta, List mdl, List MSloglik_output, int k);

List EMaximization_MSARmdl(arma::vec theta, List mdl, List MSloglik_output, int k);

List EMaximization_MSVARmdl(arma::vec theta, List mdl, List MSloglik_output, int k);

List EMiter_HMmdl(List mdl, List EMest_output, int k);

List EMiter_MSARmdl(List mdl, List EMest_output, int k);

List EMiter_MSVARmdl(List mdl, List EMest_output, int k);

List HMmdl_em(arma::vec theta_0, List mdl, int k, List optim_options);

List MSARmdl_em(arma::vec theta_0, List mdl, int k, List optim_options);

List MSVARmdl_em(arma::vec theta_0, List mdl, int k, List optim_options);

#endif