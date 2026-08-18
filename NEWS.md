# MSTest 0.1.9.9006

## Constrained EM for the transition matrix (`em_transition_constraint`)

* New control `em_transition_constraint` (default `0` = off) for the five
  Markov-switching model classes (`HMmdl()`, `MSARmdl()`, `MSARXmdl()`,
  `MSVARmdl()`, `MSVARXmdl()`): bounds every transition probability below by the
  given value during EM estimation (`p_ij >= eps`; because columns sum to 1 this is
  equivalent to the two-sided bound `eps <= p_ij <= 1-eps` used by Kasahara and
  Shimotsu (2018) and Qu and Zhuo (2021)). The constrained update is the exact
  maximizer of the expected complete-data log-likelihood over the bounded simplex
  (a closed-form water-filling solution), so the EM monotone-ascent property is
  preserved (restricted EM in the sense of Kim and Taylor, 1995). Starting values
  are projected into the bound. The default `0` (or a negative value) disables the
  constraint and is bit-identical to the previous behavior. Must be less than
  `1/k`; note the bound implies each diagonal entry is at most
  `1 - (k-1)*em_transition_constraint`, so choose it with this cap in mind for
  `k > 2`.
* New additive outputs: `em_transition_constraint` (echoed) and
  `P_constraint_binding` (k x k logical, TRUE where a fitted entry sits at the
  bound). Entries at the bound are boundary estimates and their standard errors
  are `NA`. Both apply to `method = "EM"` only; MLE fits return
  `P_constraint_binding = NULL` (use `mle_theta_low`/`mle_theta_upp` to bound the
  MLE path).
* A user-supplied `init_theta` is now validated: its transition-matrix block (the
  last `k*k` entries, `vec(P)` column-major) must be finite, in `[0,1]`, and
  column-stochastic (columns summing to 1). Previously malformed starting values
  were accepted and could be returned as invalid "fitted" transition matrices;
  they now raise an informative error.

# MSTest 0.1.9.9005

## MLE estimation (`method = "MLE"`): corrected optimizer-result handling and constraints

* Fixed a defect that made every `method = "MLE"` fit return its starting values: the
  result guard tested `res$status`, but `nloptr::slsqp()` returns `res$convergence`
  (positive codes are successes). MLE fits now return the optimized solution. All five
  Markov-switching classes were affected (`HMmdl()`, `MSARmdl()`, `MSARXmdl()`,
  `MSVARmdl()`, `MSVARXmdl()`).
* `MSVARmdl_mle()` now passes `deprecatedBehavior = FALSE` to `nloptr::slsqp()`, matching
  the other classes; previously its inequality constraints were interpreted with the
  reversed (deprecated) sign convention.
* The MSVAR stationarity constraint now bounds the companion-matrix eigenvalue moduli by
  1 (`1 - Mod(eigenvalue) >= 0`); the previous form permitted moduli up to 2, admitting
  explosive solutions. Non-finite eigenvalues map to a finite penalty so the optimizer's
  finite-difference Jacobian stays defined.
* `HMmdl()` with exogenous regressors: the variance floor was applied to the `betaZ`
  block of the parameter vector (pinning coefficients near zero) because the constraint
  selector did not skip it; the selector now covers the full parameter layout.
* `mle_variance_constraint` is now relative in every class: the smallest covariance
  eigenvalue is bounded below by `mle_variance_constraint * tr(Sigma_lin)/q`, where
  `Sigma_lin` is the one-regime fit covariance — the same anchoring as
  `em_variance_constraint`. Previously `HMmdl()`, `MSVARmdl()`, and `MSVARXmdl()` used an
  absolute floor (default `1e-3`); their defaults are now `0.01` on the relative scale.
* Starting values are projected into the user box (`mle_theta_low`/`mle_theta_upp`)
  before optimization, with the transition-probability block projected onto the
  column-sum-to-one simplex inside the box; infeasible random starts no longer error and
  burn retry attempts.
* A failed `slsqp()` run whose final iterate is finite, feasible, and better than the
  starting value is now kept (with a warning) instead of discarded.
* Transition-probability entries returned by the optimizer are clipped into [0, 1] (with
  column renormalization) before the final smoothing pass when boundary optima produce
  entries a rounding error outside the interval.
* MLE fits now return `converged` (codes 1-4) and `convergence_code` (the `nloptr`
  status); `estimMdl()` with `converge_check` no longer errors on MLE fits.
* The exported `HMmdl_mle()`, `MSARmdl_mle()`, and `MSVARmdl_mle()` raise an informative
  error when `mdl_in$sigma` is missing and a variance floor is requested.

# MSTest 0.1.9.9004

## EM estimation: corrected transition-matrix update and smoothed probabilities for Markov-switching AR/VAR models

* The EM update of the transition matrix in `MSARmdl()`, `MSARXmdl()`, `MSVARmdl()`, and
  `MSVARXmdl()` now computes the smoothed transition counts from the extended
  (regime-history) state, as required when the conditional density depends on lagged
  regimes (Hamilton 1994, sec. 22.4, eq. 22.4.16; Krolzig 1997, eq. 6.14 and sec. 9.3.5).
  Previously the counts were computed from marginal single-regime probabilities, a formula
  that is valid for hidden Markov models (p = 0) but not for switching-mean AR models;
  this could cause the EM iteration to reduce the likelihood and converge to points that
  are not likelihood optima.
* The smoothed regime probabilities (`St`) returned by the same four model classes are now
  the aggregate of the extended-state smoother over the date-t regime. Previously they were
  computed by a separate marginal recursion subject to the same limitation.
* `limP()` now guarantees a valid probability vector: at boundary or rank-deficient
  transition matrices the least-squares solve could return tiny negative entries, which
  corrupted the filter and smoother downstream; such solutions are now clipped at zero and
  renormalized. The guard is conditional, so results are bit-identical whenever the
  solution was already nonnegative (all interior transition matrices).
* The EM drivers now execute exactly `maxit` iterations and the returned `iterations`
  reports the true count (previously one additional iteration ran and the count was one
  too low).
* Louis (1982) standard errors: the expected transition counts now sum over the same
  sample range as the EM update (an order 1/T correction to transition-probability
  standard errors).
* New control `em_best_iterate` (default `TRUE`) for the five Markov-switching model
  classes: the EM returns the iterate with the highest log-likelihood visited, with the
  smoothed probabilities, residuals, and log-likelihood evaluated at that iterate (one
  additional E-step). Set to `FALSE` for the previous behavior (return the final iterate,
  whose reported log-likelihood and smoothed probabilities lag the parameters by one
  E-step). The new output element `descents` counts likelihood decreases observed along
  the EM path.
* The default of the `maxit_converge` control (maximum number of fresh starting values
  attempted when an estimation attempt returns errors or non-finite values) is now `10`
  (was `500`). The retry loop is silent and was never observed to need more than one
  attempt in extensive testing; the old default allowed a pathological worst case of
  hundreds of wasted estimations per start.
* The `MMCLRTest()` nuisance-parameter search box now applies the `CI_union` rule per
  parameter: each parameter's bounds are the union of its `eps` band and its plus/minus two
  standard-error interval where the standard error is finite. Previously a single
  non-finite standard error (for example a variance at the EM floor, whose standard error
  is `NA` by design) discarded the standard-error information for every parameter,
  shrinking the whole search box to the `eps` band.
* If every starting value fails estimation (each exhausting its `maxit_converge`
  attempts), the Markov-switching constructors now raise an informative error naming the
  model class and the relevant controls; previously this crashed with an unrelated
  low-level error (`replacement has length zero`) part-way through the multistart loop.
* The `converged` flag of EM fits obtained through `estimMdl()` now reports convergence
  under the criterion selected by the `conv` control, as the model constructors already
  did; previously the dispatcher replaced it with the parameter-change test regardless of
  the chosen criterion. MLE fits, which do not set the flag, keep the parameter-change
  fallback.
* New regression tests: EM likelihood ascent, `limP` validity at boundary transition
  matrices, and smoothed-probability consistency with the extended-state smoother.
* `HMmdl()` (p = 0) is unaffected: the marginal formulas are exact there and unchanged.
  Fits with a common (non-switching) mean are numerically unaffected beyond floating-point
  noise. Estimates for switching-mean fits change; attained log-likelihoods typically
  increase (small decreases of order 1e-2 are possible on individual datasets because the
  filter is initialized at the ergodic distribution implied by the transition matrix, whose
  dependence on the transition probabilities the update does not account for).

# MSTest 0.1.9.9002

## Pre-simulation freeze
* This development version is the frozen state used for the Rodriguez-Rondon & Dufour (2026)
  simulation study; it is not on CRAN. It will be released to CRAN as version 0.2.0 once the
  simulation study is finalized. The GitHub tag matches the package version: `v0.1.9.9002`
  is this state (pre-simulation freeze plus the EM regime-variance floor); `v0.1.9.9001` was
  the freeze alone, before the variance floor.

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

## Regime-variance floor in the EM algorithm (new default; constrained ML)
* The Gaussian Markov-switching likelihood is unbounded: a regime variance driven to zero on a
  few observations sends the likelihood to infinity, so the unconstrained maximum likelihood
  estimate does not exist (Kiefer and Wolfowitz 1956; Day 1969; Hamilton 1994, p. 689). With
  thorough multi-start estimation these degenerate modes are found often enough to distort
  Monte Carlo LR tests: the simulated null statistics acquire spurious spikes, inflating the
  simulated critical values and destroying power. (Earlier versions were partly shielded by
  numerical fragility: the degenerate paths tended to abort before converging; the numerical
  stability work in 0.1.9 made them reachable.)
* The Markov-switching constructors (`HMmdl()`, `MSARmdl()`, `MSARXmdl()`, `MSVARmdl()`,
  `MSVARXmdl()`) now constrain the regime variances during EM estimation, via the new control
  `em_variance_constraint` (default `0.01`): `sigma_k^2 >= 0.01 * sigma_linear^2`, where
  `sigma_linear^2` is the variance of the one-regime fit of the same data (for multivariate
  models the bound applies to the smallest eigenvalue of each regime covariance). Setting the
  control to `0` disables it and reproduces the previous unconstrained behaviour. The bound
  follows Kasahara and Shimotsu (2018) (on the variance rather than the standard-deviation
  scale) and is in the spirit of the constrained maximum likelihood formulation of Hathaway
  (1985, 1986).
* Clipping the variance update at the bound in each EM iteration is the exact constrained
  maximization step, so monotone ascent of the (constrained) likelihood is retained; starting
  values are projected onto the bound so this holds from the first iteration.
* Fitted models gain the additive fields `em_variance_constraint`, `sigma_floor`, and
  `sigma_floor_binding` (per-regime logical). Standard errors of variances that end at the
  bound are reported as `NA`, since they are boundary estimates.
* The constraint applies identically to the observed data and to every simulated sample inside
  `LMCLRTest()` and `MMCLRTest()` (each sample's bound comes from its own linear fit), so the
  Monte Carlo p-values remain exact: Monte Carlo validity holds for any statistic computed by
  the same rule on both sides, and the constrained-likelihood ratio is such a statistic. NOTE:
  transition probabilities remain unconstrained, so the tests continue to allow boundary values
  of the transition matrix.

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
  makes both sequential and parallel runs fully reproducible. In the sequential legacy
  configuration (`mc_seed = NULL`, `crn = FALSE`, `workers = 0`) no internal seeding happens at
  all, so a script-seeded legacy run reproduces version 0.1.9 bit-for-bit (verified on the
  R Journal article's MMC example).
* Documentation note: GenSA is deterministic given the fixed starting value used by
  `MMCLRTest()`, so its search trajectory does not depend on any seed; pso and GA consume the
  optimizer seed.

## Other
* `MMCLRTest()`: when the p-value at the initial values already satisfies
  `pval_0 >= threshold_stop`, the nuisance-parameter search is skipped entirely and
  `pval = pval_0` is returned (`mmc_optimout` carries an early-stop marker). Previously the
  optimizer was always launched and stopped via its own threshold rule, wasting one evaluation
  (GenSA) or a full first swarm (pso). With `threshold_stop` slightly above the significance
  level, most true-null replications in a simulation study now cost no more than a Local MC run.
* `MMCLRTest()`: user-supplied `mdl_h0_control`/`mdl_h1_control` sublists are now merged with
  their defaults instead of replacing them, fixing an "argument is of length zero" error when a
  supplied sublist omitted `getSE`.
* `LMCLRTest()`: the p-value is computed before the null draws are sorted (sorting is now only
  used for the reported critical values); numerically equivalent, clearer semantics.
* Documentation fix: the `sigma` element of the DGP list for the univariate simulators
  (`simuAR()`, `simuARX()`, `simuMSAR()`, `simuMSARX()`) was documented as a standard
  deviation; it is (and always was) a **variance** (per-regime variances for the
  Markov-switching versions), consistent with the multivariate simulators and the estimation
  output convention. Documentation only; no change to any computation.

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
