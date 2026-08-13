# MSTest 0.2.0

*(Development release for the Rodriguez-Rondon & Dufour (2026) simulation study; tagged
`v0.2.0` on GitHub, not submitted to CRAN.)*

## IMPORTANT: results of the Monte Carlo LR tests change relative to 0.1.9
* For a given seed, `LMCLRTest()` and `MMCLRTest()` produce different numerical results than
  version 0.1.9 (p-values, critical values, and possibly the observed statistic, which can only
  weakly increase). This is a consequence of the warm-start, LR >= 0, and common-random-numbers
  changes described below. Everything else is unchanged bit-for-bit: model estimation
  (`ARmdl()`, `MSARmdl()`, etc. with default controls), the simulation functions (`simu*`), and
  the other testing procedures (`DLMCTest()`, `DLMMCTest()`, `CHPTest()`, `HLRTest()`) reproduce
  0.1.9 exactly. The R Journal article results are reproduced exactly by CRAN version 0.1.9
  (see the note in `inst/examples/article.R`). Setting `init_method = "random"` and `crn = FALSE`
  recovers the legacy behavior of the LR tests.

## Warm-started estimation of the alternative model (new default)
* `LMCLRTest()` and `MMCLRTest()` gain an `init_method` control (default `"warmstart"`). The
  alternative (`k1`-regime) model is now also initialized at deterministic warm starts built by
  embedding the fitted null model into the `k1`-regime parameter space (duplicating regimes with
  a split transition column, so the embedded point attains exactly the null log-likelihood, then
  applying a small antisymmetric perturbation; one start per way of distributing the extra
  regimes). These compete with the existing `use_diff_init` random starts, and the identical
  procedure is applied to the observed data and to every simulated draw, preserving the
  exchangeability underlying the Monte Carlo p-value. See `warmstart_theta()` (new, exported).
  The embedding follows the null-hypothesis geometry of Kasahara & Shimotsu (2018) and the
  component-splitting initialization of the EM-test literature (Chen & Li 2009).
* Model constructors (`HMmdl()`, `MSARmdl()`, `MSARXmdl()`, `MSVARmdl()`, `MSVARXmdl()`) gain an
  additive `init_theta_extra` control: a list of extra deterministic starting vectors that
  compete with the random multi-start draws (unlike `init_theta`, which bypasses the multi-start
  search). Default behavior is unchanged bit-for-bit.

## LR statistic numerical policy: floored at 0, never an error
* For nested models the LR statistic is non-negative by construction; a negative computed value
  is a numerical optimization artifact. `LMCLRTest()` and `MMCLRTest()` no longer stop with
  "LRT_0 is negative": the observed statistic is set to 0 with a warning (when the magnitude is
  non-negligible), and simulated null statistics are treated identically inside `LR_samp_dist()`
  (previously, negative simulated draws were silently redrawn). Applying the same rule on both
  sides preserves exchangeability; the induced ties at 0 are handled by the randomized
  tie-breaking rule of Dufour (2006) already implemented in `MCpval()` and can only make the
  test conservative. This also removes a source of bias in simulation studies, where
  discarded negative-statistic replications previously inflated empirical rejection rates.
* For `k0 > 1`, the null log-likelihood is floored at the one-regime (linear) fit — a feasible
  point of the null space — identically on the observed and simulated sides.

## Common random numbers (CRN) across the MMC nuisance-parameter search
* `MMCLRTest()` gains a `crn` control (default `TRUE`): an identical estimation-RNG stream (EM
  starting values; master and parallel-worker streams) is replayed at every candidate theta
  evaluation, so that together with the pre-drawn innovations the MC p-value is a deterministic
  function of theta during the search — implementing Dufour (2006, Prop. 4.2) with the
  estimation randomness folded into the fixed disturbance vector u. This removes optimizer
  objective noise (which is size-safe but costs power and can trigger `threshold_stop`
  spuriously) and makes the search fully reproducible. The RNG state is saved and restored
  around each evaluation, so the optimizer's own search randomness and the user's session RNG
  are unaffected.
* The p-value at the initial (estimated) nuisance parameters is now evaluated with the same
  worker configuration and CRN stream as the optimizer's evaluations and returned as `pval_0`
  (additive output field); under `crn = TRUE` the reported `pval` satisfies `pval >= pval_0` by
  construction.

## Seeding scheme
* Internal seeds (optimizer trajectory, worker RNG streams, CRN stream) are now drawn from the
  `mc_seed` stream via `sample.int()` instead of being derived as `mc_seed + 1`, `mc_seed + 2`.
  The arithmetic scheme made replication r's optimizer stream identical to replication (r+1)'s
  data/estimation stream whenever an outer simulation study used consecutive `mc_seed` values —
  correlating replications. Drawn seeds remove this systematic collision (a draw-then-reseed
  construction keeps the observed-data estimation and the pre-drawn innovations bit-identical
  to 0.1.9 for a given `mc_seed`; only optimizer/worker streams change). When `mc_seed` is
  `NULL`, internal seeds are drawn from the ambient RNG state, so a script-level `set.seed()`
  makes both sequential and parallel runs fully reproducible.
* Documentation note: GenSA is deterministic given the fixed starting value used by
  `MMCLRTest()`, so its search trajectory does not depend on any seed; pso and GA consume the
  optimizer seed.

## Other
* `MMCLRTest()`: user-supplied `mdl_h0_control`/`mdl_h1_control` sublists are now merged with
  their defaults instead of replacing them, fixing an "argument is of length zero" error when a
  supplied sublist omitted `getSE`.
* `LMCLRTest()`: the p-value is computed before the null draws are sorted (sorting is now only
  used for the reported critical values); numerically equivalent, clearer semantics.

# MSTest 0.1.9

## Reproducibility note
* The internal standard-normal generator `randSN()` now uses `arma::randn` instead of a Box-Muller
  transform (faster, and avoids a `log(0)` edge case). IMPORTANT: this changes the random numbers
  drawn for a given `set.seed()`, so simulated series (the `simu*` functions) and Monte Carlo
  p-values (`LMCLRTest()`, `MMCLRTest()`, `DLMCTest()`, `DLMMCTest()`, `CHPTest()`) differ
  numerically from earlier versions. Test conclusions and validity are unaffected.

## Testing procedures
* `MMCLRTest()`: the observed likelihood ratio statistic is now held fixed across the
  nuisance-parameter search, following Dufour (2006), instead of being recomputed at each candidate
  value. This corrects the maximized Monte Carlo p-value (which could previously be slightly liberal).
* `MMCLRTest()`: pre-drawn simulation innovations are now held fixed across the optimization
  (fixed-error Monte Carlo; Dufour 2006, Prop. 4.2).
* `LMCLRTest()` and `MMCLRTest()`: added an `mc_seed` control for fully reproducible Monte Carlo
  p-values; the number of parallel workers is capped at `N` when more workers than replications are
  requested.
* `LMCLRTest()` and `MMCLRTest()`: robust handling of failed simulated draws. Non-finite draws are
  dropped before the p-value is computed. If the simulated null distribution cannot be built at all
  (every draw fails even after the re-draw safety), the functions now stop with an informative message
  (increase `use_diff_init` or inspect the fit) instead of returning an invalid value or crashing the
  optimizer. In `MMCLRTest()`, a candidate parameter value whose null cannot be simulated is penalized
  so the optimizer avoids it; the initial value (`theta_0`) is exempt so that the MMC p-value remains
  at least as large as the LMC p-value.
* `MMCLRTest()`: the GA optimizer now uses `popSize = 10` by default (GA's own default of 50 makes the
  search perform 50 x `maxit` expensive evaluations); additional GA controls can be passed through
  `optim_control`.
* `MCpval()`: the `type` argument now accepts both long and short spellings
  (`"geq"`, `"leq"`, `"two-tailed"`/`"two-tail"`, `"absolute"`/`"abs"`) and raises an informative error
  for an unrecognized value (previously it silently returned a sentinel). The documentation now matches
  the accepted values.
* `DLMMCTest()`: fixed the sign of the stationarity penalty so non-stationary candidate parameters
  are correctly avoided by the optimizer (previously, in some cases the maximized p-value could be
  returned as the raw penalty constant).
* `DLMCTest()` and `DLMMCTest()`: corrected the sample size used to simulate the null distribution
  of the moment-based statistics. The simulated moments are now computed from samples of length
  `T - p`, matching the number of AR(`p`) residuals used for the observed statistic and for the
  p-value calibration in `approxDistDL()` (previously `T`). This restores the exact exchangeability
  underlying the Monte Carlo p-value (Dufour & Luger 2017); the effect is negligible for large
  samples and grows with `p/T`.
* `print()`/`summary()` for `DLMMCTest` objects now display the nuisance parameter value that
  maximizes the Monte Carlo p-value (`theta_max_min`/`theta_max_prod`) alongside the moment
  statistics, matching the output of `DLMCTest` (display only; the returned object is unchanged).

## Estimation (EM)
* New `conv` control for the EM stopping criterion in `HMmdl()`, `MSARmdl()`, `MSARXmdl()`,
  `MSVARmdl()`, and `MSVARXmdl()`, with options `"loglik"` (relative log-likelihood change; the new
  default), `"theta"` (relative parameter change, the previous behavior), `"both"` (both, following
  Krolzig 1997), and `"loglik-A"`/`"both-A"` (Aitken-accelerated log-likelihood; Böhning et al. 1994,
  McLachlan and Krishnan 2008). A separate `ltol` control sets the log-likelihood tolerance
  (default `1e-7`); `thtol` (default `1e-6`) remains the parameter tolerance. NOTE: the default change
  from parameter-based stopping to `"loglik"` can change estimates slightly for a given `set.seed()`;
  set `conv = "theta"` to recover the previous behavior. Simulations show `"loglik"` matches the size
  and power of `"theta"` while converging substantially faster. Each fitted model now also reports a
  `converged` flag.
* Corrected the M-step mean update for Markov-switching AR and VAR models following Krolzig (1997)
  (the previous weighted-average update is exact only for hidden Markov models).
* Comprehensive numerical-stability improvements to the EM algorithm: guarded matrix inversions
  (`solve()` with fallback), Hamilton-filter underflow handling, degenerate-regime guards,
  positive-definite covariance regularization, and stationarity checks. Several data sets that
  previously triggered `solve(): solution not found` now estimate successfully.
* Faster Hamilton filter and Kim smoother (O(M) forward/backward recursions).
* Correct Kim-smoother handling for the transition-probability update (renormalization fix) and
  additional guards against transition-matrix corruption in degenerate replications.
* Fixed `MSVARXmdl()` to call `MSVARXmdl_em()` (not `MSVARmdl_em()`) in the single-initial-value EM
  branch (used when `use_diff_init = 1` with a user-supplied `init_theta`), so exogenous regressors
  are handled correctly on that path.

## Standard errors
* New `se_method = "louis"` control for the Markov-switching model constructors (`HMmdl()`,
  `MSARmdl()`, `MSARXmdl()`, `MSVARmdl()`, `MSVARXmdl()`), computing expected-complete-data
  standard errors (Louis 1982); automatic fallback from the observed-information Hessian to Louis
  when the Hessian is ill-conditioned.
* Transition-matrix standard errors now use a reduced (free-parameter) parameterization, fixing
  inflated or `NA` standard errors for `k >= 2`. Per-parameter step sizes are used in the numerical
  Hessian to keep transition probabilities within `[0, 1]`.

## Other
* Model constructors now validate `use_diff_init >= 1` (and `maxit_converge >= 1`) with a clear error
  instead of failing cryptically.
* C++ sources: made integer/double literal usage explicit throughout (no change to results).
* Added a unit-test suite (`testthat`) and an introductory vignette.

# MSTest 0.1.8
* made changes to some function examples.

# MSTest 0.1.7
* Made change in trycatch error message for HLRTest() to better describe optimization issues. 
* Made minor changes to article.R example file.


# MSTest 0.1.6
* Updated MLE estimation following deprecation of hin>=0 (inequality constraint direction) in slsqr
* Changed OLS unbiased estimates of models with k=1 to be consistent with MLE estimates. Package is for testing more than estimation so comparison with MLE-based tetsing is prioritized. 
* Changed use of arma::is_finite(X) to std::isfinite(X) because former is now deprecated.


# MSTest 0.1.5
* patch to compile following (Rcpp) Armadillo update (i.e., added proper namespace scopeas_scalar).


# MSTest 0.1.2
* Made change to MSARmdl(), MSVARmdl(), and HMmdl(). Specifically, when msmu or msvar is FALSE, output list copies single regime value for each k regime. This is needed for simulating the null when either msmu or msvar is FALSE in LMCLRTest() and MMCLRTest().
* Added option to use different number of initial values for estimating MSMs with observed data vs. for null distribution (see documentation for LMCLRTest() and MMCLRTest()).
* Updated USGNP data set to include 2022 Q3.
* MMCLRTest() now has the option to add lower and upper bounds for autoregressive coefficients and transition probabilities to help reduce errors from polyroot() when optimizing. MMC_bounds() has been updated to reflect this. 
* DLMMCTest() now has the option to add lower and upper bounds for autoregressive coefficients to help reduce errors from polyroot() when optimizing. DLMMC_bounds() has been created for this. 
* Fixed error in bootCV() used by CHPTest(). When NaN occurs, new draw is used. 
* Added option to allow user to specify 'mle_theta_low' and 'mle_theta_upp' which determine the lower and upper bounds for optimization in HMmdl(), MSARmdl(), and MSVARmdl() when "method='MLE'" is specified.
* print.CHPTest() now used two lines to print description.
* In HLRTest() user can now define entire grid for transition probabilities. 
* added new optional optimization routine for HLRparamSearch() (nloptr::slsqp() can be used now).
* Fixed bug in HLRTest() where grid for sigma options are properly used now. 
* Now using nearest_spd() from pracma instead of nearPD() from lmf package.
* Added classes for simulation functions
* Added new methods, namely: coef, fitted, predict, summary, residuals, nobs, plot
* Made changes to print method.  
* Changed methods for models, namely: logLiklihood now uses logLik, and aic and bic now use AIC and BIC methods.
* Updated USGNP dataset to include data up to end of 2023.  
* Updated `README.md` file  with usage of new methods.
* Updated `DESCRIPTION` file for changed dependencies and new version.

# MSTest 0.1.1

* Made changes to MMC LRT related functions for obtaining null distribution of statistical test.
* Added more examples. Examples that take long to complete are commented out but serve to get familiar with usage. 
* Fixed bug related to using init_theta (setting initial values of parameters) when estimating Markov models (i.e. MSARmdl(), HMmdl(), and MSVARmdl()). 
* Fixed bug in MMC_bound() when k0>1
* Updated `DESCRIPTION` file for new version.

# MSTest 0.1.0

* Added a `NEWS.md` file to track changes to the package.
* Added a `README.md` file to describe the package.
* First public version of package.
