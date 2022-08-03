# MSTest

This package implements hypothesis testing procedures that can be used to identify the number of regimes in a Markov-Switching model. It includes the likelihood ratio test described in Hansen (1992), the optimal test for regime switching of Carrasco, Hu, & Ploberger (2014), the Monte Carlo moment-based test of Dufour & Luger (2017), the parametric bootstrap test described in Qu & Zhuo (2021) and Kasahara & Shimotsu (2018), and finally the Monte Carlo Likelihood ratio tests of Rodriguez Rondon & Dufour (2022). 

In addition to testing procedures, the package also includes functions that can be used to simulate: autoregressive, vector autoregressive, Markov switching autoregressive, and Markov switching vector autoregressive processes among others. Model estimation procedures are also available.

For a more detailed description of this package see Rodriguez Rondon & Dufour (2022).

## Installation

Currently, this package is only available through github. To install it you can use the following (requires 'devtools' package): 

```{r}
devtools::install_github("roga11/MSTest")
```

## Examples



## References

Carrasco, Marine, Liang Hu, and Werner Ploberger. 2014. “Optimal test for Markov switching parameters.” Econometrica 82 (2): 765–784.

Dufour, Jean-Marie, and Richard Luger. 2017. “Identification-robust moment-based tests for Markov switching in autoregressive models.” Econometric Reviews 36 (6-9): 713–727.

Kasahara, Hiroyuki, and Katsumi Shimotsu. 2018. “Testing the number of regimes in Markov regime switching models.” arXiv preprint arXiv:1801.06862.

Hansen, Bruce E. 1992. “The likelihood ratio test under nonstandard conditions: testing the Markov switching model of GNP.” Journal of applied Econometrics 7 (S1): S61–S82.

Rodriguez Rondon, Gabriel and Jean-Marie Dufour. 2022. “Monte Carlo Likelihood Ratio Tests for Markov Switching Models.” Unpublished manuscript.

Rodriguez Rondon, Gabriel and Jean-Marie Dufour. 2022. “MSTest: An R-package for Testing Markov-Switching Models.” Unpublished manuscript.

Qu, Zhongjun, and Fan Zhuo. 2021. “Likelihood Ratio-Based Tests for Markov Regime Switching.” The Review of Economic Studies 88 (2): 937–968
