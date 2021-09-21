#ifndef MODELS_H
#define MODELS_H

#include <RcppArmadillo.h>
using namespace Rcpp;

List ARmdl(arma::vec Y, int ar, bool intercept = 1);

double MSloglik_fun(arma::vec theta, List mdl, int k);

List MSloglik(arma::vec theta, List mdl, int k);

List MSVARloglik(arma::vec theta, List mdl, int k);

List EMaximization(arma::vec theta, List mdl, List MSloglik_output, int k);

List VAREMaximization(arma::vec theta, List mdl, List MSloglik_output, int k);

List EMiter(List mdl, List EMest_output, int k);

List EMest(arma::vec theta_0, List mdl, int k, List optim_options);

List MSARmdl(arma::vec Y, int ar = 0, int k = 2, bool msmu = 1, bool msvar = 1, int maxit = 10000, double thtol = 1.e-8, bool getHess = 0, int max_init= 100);

#endif