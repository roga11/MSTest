#ifndef MODELS_H
#define MODELS_H

#include <RcppArmadillo.h>
using namespace Rcpp;

double AR_loglik_fun(arma::vec theta, List mdl);
  
List ARmdl(arma::vec Y, int ar, bool intercept = 1, bool getSE = 0);

double VAR_loglik_fun(arma::vec theta, List mdl);

List VARmdl(arma::mat Y, int ar, bool intercept = 1, bool getSE = 0);

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

List MSARmdl(arma::vec Y, int ar, int k, bool msmu = 1, bool msvar = 1, int maxit = 10000, double thtol = 1.e-6, bool getSE = 0, int max_init = 500, int use_diff_init = 1, Nullable<NumericVector> init_value = R_NilValue);

List MSVARmdl(arma::mat Y, int ar, int k, bool msmu = 1, bool msvar = 1, int maxit = 10000, double thtol = 1.e-8, bool getSE = 0, int max_init = 500, int use_diff_init = 1, Nullable<NumericVector> init_value = R_NilValue);

#endif