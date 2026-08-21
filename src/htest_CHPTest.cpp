#include <RcppArmadillo.h>
#include "methods.h"
//[[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;

// =============================================================================
// Carrasco, Hu & Ploberger (2014) parameter stability test.
//
// The nuisance direction h = (hu(1), 0, ..., 0, hu(2)) always has zero entries
// at the AR-coefficient positions (Section 6.2 of CHP: "the other elements are
// all set to be 0" -- phi is a nuisance parameter under H0 but is not itself
// tested for switching). Sandwiching any (mu, phi, sigma2)-indexed matrix by
// such an h therefore only ever picks out its (mu, sigma2) block; the phi
// cross-terms cancel identically regardless of their value. mu2t below is
// built directly from that reduced (mu, sigma2) block rather than assembling
// and sandwiching the full matrix, which both simplifies the code and removes
// its dominant cost. The AR-coefficient score (ltphi) is still needed and
// still computed in full: mu2t is projected onto the full nuisance score space
// (mu, phi, sigma2) via ltmx, which is what actually orthogonalizes the
// statistic against phi.
// =============================================================================

// ---- lean AR(p) OLS fit (no standard errors, no S3 bookkeeping) -----------
struct ARfit {
  arma::vec phi;
  double mu;
  double stdev;
  arma::vec resid;
  arma::mat x;
};

ARfit ar_ols_fit(const arma::vec& Y, int p){
  List lagged = ts_lagged(arma::mat(Y), p);
  arma::vec y = lagged["y"];
  arma::mat x = lagged["X"];
  int n = y.n_elem;
  arma::mat zz = arma::join_rows(arma::ones(n, 1), x);
  arma::vec beta = arma::solve(zz.t() * zz, zz.t() * y, arma::solve_opts::allow_ugly);
  double inter = beta(0);
  arma::vec phi = beta.subvec(1, p);
  arma::vec resid = y - zz * beta;
  double stdev = std::sqrt(arma::as_scalar(resid.t() * resid) / n);
  ARfit out;
  out.phi = phi;
  out.mu = inter / (1.0 - arma::sum(phi));
  out.stdev = stdev;
  out.resid = resid;
  out.x = x;
  return out;
}

// ---- mu2t sandwich, per (rho, hu1, hu2), across the whole sample -----------
// Recurrence follows Carrasco, Hu & Ploberger (2014), Section 6, eq. (3.3):
//   self term (t) = ( hu'*M_t*hu + (hu.l_t)^2 ) / 2
//   cross term (t > 1, i.e. Gauss's second sample point onward) accumulates
//   an exponentially-weighted (by rho) sum of earlier scores.
// mtmu, ltmu, mtmusig2, mtsig2, ltsig2 are all length nn (or mtmu a scalar);
// with msvar = false, hu2 = 0 and the sigma2 terms drop out identically.
arma::vec mu2t_vec(double rho, double hu1, double hu2, double mtmu,
                   const arma::vec& ltmu, const arma::vec& mtmusig2,
                   const arma::vec& mtsig2, const arma::vec& ltsig2, int nn){
  arma::vec mu2t(nn, arma::fill::zeros);
  arma::vec l(nn);   // hu . (score at t), t = 0..nn-1
  for (int t = 0; t < nn; t++){
    double lt = hu1 * ltmu(t) + hu2 * (ltsig2.n_elem ? ltsig2(t) : 0.0);
    l(t) = lt;
    double mt = hu1 * hu1 * mtmu;
    if (mtmusig2.n_elem){
      mt += 2.0 * hu1 * hu2 * mtmusig2(t) + hu2 * hu2 * mtsig2(t);
    }
    mu2t(t) = (mt + lt * lt) / 2.0;
  }
  // cross term: xs(t) = rho * (xs(t-1) + l(t-1)) for t >= 2 (Gauss 1-indexed),
  // with xs(1) preset to rho*l(0) and its contribution added at t = 1 only.
  double xs = rho * l(0);
  if (nn > 1) mu2t(1) += l(1) * xs;
  for (int t = 2; t < nn; t++){
    xs = rho * (xs + l(t - 1));
    mu2t(t) += l(t) * xs;
  }
  return mu2t;
}

// sqrt(2*pi) * exp(z^2/2) * Phi(z), z = tspe - 1 -- the expTS integrand.
// Computed in log-space: for very negative z, exp(z^2/2) overflows to Inf
// while Phi(z) underflows to 0, and the naive product evaluates to Inf*0 =
// NaN even though the true limit is a small finite number (the two terms
// cancel at the same rate, per the standard Mills-ratio tail expansion of the
// normal CDF). R::pnorm's own log.p path is accurate in that tail.
double exp_ts_term(double tspe){
  double z = tspe - 1.0;
  double log_term = 0.5 * std::log(2.0 * M_PI) + z * z / 2.0 + R::pnorm(z, 0.0, 1.0, 1, 1);
  return std::exp(log_term);
}

// ---- supTS/expTS from a fitted AR(p) model -------------------------------
arma::vec chp_stat_core(const arma::vec& resid, const arma::mat& x, double mu,
                        const arma::vec& phi, double stdev, double rho_b, bool msvar){
  int nn  = resid.n_elem;
  int nar = phi.n_elem;
  int tt  = nn + nar;
  double v2 = stdev * stdev;
  double phisum = arma::sum(phi);

  arma::vec ltmu = resid * (1.0 - phisum) / v2;
  double mtmu = -std::pow(1.0 - phisum, 2.0) / v2;
  arma::mat ltphi(nn, nar);
  for (int j = 0; j < nar; j++) ltphi.col(j) = (resid % (x.col(j) - mu)) / v2;

  // ltsig2 (the sigma2 score) is needed in the projection matrix ltmx below
  // regardless of msvar -- CHP (2014), Remark 1, p. 770: epsilon-hat is the
  // residual from projecting mu2t on the WHOLE score vector l_t^(1), including
  // the coefficients not under test. mtmusig2/mtsig2 (curvature terms) are only
  // needed for the mu2t sandwich itself, which only has a sigma2 component when
  // msvar = TRUE (hu2 = 0 otherwise, so they never contribute either way).
  arma::vec ltsig2 = arma::square(resid) / (2.0 * v2 * v2) - 1.0 / (2.0 * v2);
  arma::vec mtmusig2, mtsig2;
  if (msvar){
    mtmusig2 = -resid * (1.0 - phisum) / (v2 * v2);
    mtsig2   = 1.0 / (2.0 * v2 * v2) - arma::square(resid) / (v2 * v2 * v2);
  }

  int ncol = nar + 2;
  arma::mat ltmx(nn, ncol);
  ltmx.col(0) = ltmu;
  ltmx.cols(1, nar) = ltphi;
  ltmx.col(nar + 1) = ltsig2;

  arma::mat A = ltmx.t() * ltmx;
  arma::mat Ainv;
  bool inv_ok = arma::inv_sympd(Ainv, A);
  if (!inv_ok) Ainv = arma::pinv(A);

  int seqlen = static_cast<int>(std::round((rho_b - (-rho_b)) * 100)) + 1;
  arma::vec rhogrid = arma::linspace<arma::vec>(-rho_b, rho_b, seqlen);
  double tol = 1e-5;

  double supTS;
  double expTS;

  if (!msvar){
    arma::vec cv(seqlen), cv2(seqlen);
    for (int ir = 0; ir < seqlen; ir++){
      arma::vec mu2t = mu2t_vec(rhogrid(ir), 1.0, 0.0, mtmu, ltmu, mtmusig2, mtsig2, ltsig2, nn);
      double gamma_e = arma::sum(mu2t) / std::sqrt((double) tt);
      arma::vec sol = Ainv * (ltmx.t() * mu2t);
      arma::vec eps = mu2t - ltmx * sol;
      double esqe = arma::mean(arma::square(eps));
      if (esqe < tol){
        cv(ir) = 0.0; cv2(ir) = 1.0;
      } else {
        double tspe = gamma_e / std::sqrt(esqe);
        cv(ir)  = std::pow(std::max(tspe, 0.0), 2.0) / 2.0;
        cv2(ir) = exp_ts_term(tspe);
      }
    }
    supTS = cv.max();
    expTS = arma::mean(cv2);
  } else {
    arma::vec cv(100 * seqlen), cv2(100 * seqlen);
    int idx = 0;
    for (int ih = 0; ih < 100; ih++){
      double hu1, hu2, sumsq;
      do {
        arma::vec draw = arma::randn(2);
        sumsq = arma::as_scalar(draw.t() * draw);
        hu1 = draw(0) / sumsq;
        hu2 = draw(1) / sumsq;
      } while (!std::isfinite(hu1) || !std::isfinite(hu2) || sumsq < 1e-10);
      for (int ir = 0; ir < seqlen; ir++){
        arma::vec mu2t = mu2t_vec(rhogrid(ir), hu1, hu2, mtmu, ltmu, mtmusig2, mtsig2, ltsig2, nn);
        double gamma_e = arma::sum(mu2t) / std::sqrt((double) tt);
        arma::vec sol = Ainv * (ltmx.t() * mu2t);
        arma::vec eps = mu2t - ltmx * sol;
        double esqe = arma::mean(arma::square(eps));
        if (esqe < tol){
          cv(idx) = 0.0; cv2(idx) = 1.0;
        } else {
          double tspe = gamma_e / std::sqrt(esqe);
          cv(idx)  = std::pow(std::max(tspe, 0.0), 2.0) / 2.0;
          cv2(idx) = exp_ts_term(tspe);
        }
        idx++;
      }
    }
    supTS = cv.max();
    expTS = arma::mean(cv2);
  }

  arma::vec out(2);
  out(0) = supTS;
  out(1) = expTS;
  return out;
}

//' @title Test statistic for CHP 2014 parameter stability test
//'
//' @description Computes the supTS and expTS test statistics of Carrasco, Hu
//'   & Ploberger (2014) for a fitted AR(p) model under the null of parameter
//'   constancy.
//'
//' @param mdl List containing model attributes (see \code{\link{ARmdl}}):
//'   \code{resid}, \code{x}, \code{mu}, \code{phi}, \code{stdev}.
//' @param rho_b Number determining bounds for the grid search over \code{rho}
//'   (i.e. \code{rho} in \code{[-rho_b, rho_b]}, step \code{0.01}).
//' @param msvar Boolean. If \code{TRUE}, the alternative has switching mean
//'   and variance; if \code{FALSE}, switching mean only.
//'
//' @return A (\code{2 x 1}) vector with the supTS statistic first and the
//'   expTS statistic second.
//'
//' @keywords internal
//'
//' @references Carrasco, Marine, Liang Hu, and Werner Ploberger. 2014. “Optimal
//'   test for Markov switching parameters.” \emph{Econometrica} 82 (2): 765–784.
//'
//' @export
// [[Rcpp::export]]
arma::vec chpStat(List mdl, double rho_b, bool msvar){
  arma::vec resid = mdl["resid"];
  arma::mat x      = mdl["x"];
  double mu        = mdl["mu"];
  arma::vec phi    = mdl["phi"];
  double stdev     = mdl["stdev"];
  return chp_stat_core(resid, x, mu, phi, stdev, rho_b, msvar);
}

//' @title Bootstrap critical values for CHP 2014 parameter stability test
//'
//' @description Parametric bootstrap described on p. 771 of Carrasco, Hu &
//'   Ploberger (2014): simulate under the null at the fitted parameters,
//'   re-estimate, and recompute the test statistic, entirely in C++ (no calls
//'   back into R per replication).
//'
//' @param mdl List containing model attributes (see \code{\link{ARmdl}}):
//'   \code{n}, \code{mu}, \code{sigma}, \code{phi}.
//' @param rho_b Number determining bounds for the grid search over \code{rho}.
//' @param N Number of bootstrap replications.
//' @param msvar Boolean. If \code{TRUE}, the alternative has switching mean
//'   and variance; if \code{FALSE}, switching mean only.
//' @param burnin Burn-in length used to simulate each bootstrap sample from
//'   its stationary distribution. Default is \code{100}.
//'
//' @return A (\code{N x 2}) matrix with the bootstrapped supTS statistic in
//'   the first column and expTS in the second.
//'
//' @references Carrasco, Marine, Liang Hu, and Werner Ploberger. 2014. “Optimal
//'   test for Markov switching parameters.” \emph{Econometrica} 82 (2): 765–784.
//'
//' @export
// [[Rcpp::export]]
arma::mat CHPbootCV(List mdl, double rho_b, int N, bool msvar, int burnin = 100){
  int ar = mdl["p"];
  // mdl$n (from ARmdl) is the AR-adjusted observation count, T - p, not the
  // original series length T; each bootstrap draw is simulated at the full T
  // so that re-estimating AR(p) on it recovers the same T - p usable
  // observations the original fit used.
  List sim_dgp = clone(mdl);
  int Tfull = as<int>(mdl["n"]) + ar;
  sim_dgp["n"] = Tfull;
  arma::vec supb(N), expb(N);
  int itb = 0;
  while (itb < N){
    List sim = simuAR_cpp(sim_dgp, burnin);
    arma::vec y0 = sim["y"];
    ARfit fit = ar_ols_fit(y0, ar);
    arma::vec cv4 = chp_stat_core(fit.resid, fit.x, fit.mu, fit.phi, fit.stdev, rho_b, msvar);
    if (cv4.is_finite()){
      supb(itb) = cv4(0);
      expb(itb) = cv4(1);
      itb++;
    }
  }
  return arma::join_rows(supb, expb);
}
