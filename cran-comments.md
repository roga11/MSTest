## Submission

This is a minor release (0.1.9) with bug fixes, numerical-stability improvements, and
documentation updates. All changes are backward compatible: existing function signatures,
argument names and order, and returned object structures are unchanged. New functionality is
provided through new optional arguments whose defaults reproduce the previous behavior.

Selected user-visible changes since 0.1.8 (full list in NEWS.md):

* `MMCLRTest()`: the observed likelihood-ratio statistic is now held fixed across the
  nuisance-parameter search, following Dufour (2006); the simulated innovations are also held
  fixed across the optimization.
* `DLMCTest()` / `DLMMCTest()`: corrected the stationarity-penalty sign and the sample size
  used to simulate the moment-based null distribution.
* Estimation (EM): corrected the M-step mean update for Markov-switching AR and VAR models
  following Krolzig (1997), together with numerical-stability and speed improvements.
* New optional controls (`conv`, `ltol`, `se_method`, `mc_seed`) with backward-compatible
  defaults.
* Added a unit-test suite (testthat) and an introductory vignette.

## Maintainer email change

The maintainer email address has changed from gabriel.rodriguezrondon@mail.mcgill.ca to
gabrodriguezrondon@gmail.com. This is the same maintainer, Gabriel Rodriguez-Rondon; the
institutional (McGill) address is being retired. I can be reached at the new address and can
confirm the change from the previous address if required.

## Test environments

* local: macOS (aarch64), R 4.4.0
* win-builder: R-release (R 4.6.1) and R-devel

## R CMD check results

0 errors | 0 warnings | 1 note

On win-builder (both R-release and R-devel) the only NOTE is the maintainer email change
described above (new maintainer gabrodriguezrondon@gmail.com; old maintainer
gabriel.rodriguezrondon@mail.mcgill.ca). This change is intentional, and the maintainer is the
same person.

## Reverse dependencies

There are no reverse dependencies for this package.
