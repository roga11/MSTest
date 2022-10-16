#' @title US GNP data 1947Q2 - 2022Q2 
#' 
#' @format 
#' \describe{This data is used in Rodriguez Rondon & Dufour (2022) working paper. The series ranges from 
#' 1947Q2 to 2022Q2.
#'   \item{DATE}{Vector of dates}
#'   \item{GNP_logdiff}{log difference of US GDP series}
#'   \item{GNP}{US GNP seies}
#' }
#' @source \url{https://fred.stlouisfed.org/graph/?g=UKDQ}
#' @references Rodriguez Rondon, Gabriel and Jean-Marie Dufour. 2022. “Monte Carlo Likelihood Ratio Tests for Markov Switching Models.” \emph{Unpublished manuscript}.
"USGNP"
# -----------------------------------------------------------------------------
#' @title Carrasco, Hu, & Ploberger 2010 GNP data
#' 
#' @format 
#' \describe{This data is the extension of the GNP series used in CHP (2014), Econometrica. This series ranges from 
#' 1951Q2 to 2010Q4. 
#'   \item{DATE}{Vector of dates}
#'   \item{GNP_logdiff}{log difference of US GNP series}
#'   \item{GNP}{US GNP series}
#' }
#' @source \url{https://www.econometricsociety.org/content/supplement-optimal-test-markov-switching-parameters}
#' 
#' @references Carrasco, Marine, Liang Hu, and Werner Ploberger. 2014. “Optimal test for Markov switching parameters.” \emph{Econometrica} 82 (2): 765–784.
#' 
"chp10GNP"
# -----------------------------------------------------------------------------
#' @title Hamilton 1984 & Hansen 1992 GNP data
#' 
#' @format 
#' \describe{This data set is used in Hansen (1992) to test the US GNP model proposed by 
#' Hamilton (1989). This series ranges from 1951Q2 to 1984Q4. 
#'   \item{DATE}{Vector of dates}
#'   \item{GNP_logdiff}{US GNP log difference}
#'   \item{GNP}{US GNP series}
#' }
#' @source \url{https://www.ssc.wisc.edu/~bhansen/progs/jae_92.html}
#' 
#' @references Hansen, Bruce E. 1992. “The likelihood ratio test under nonstandard conditions: testing the Markov switching model of GNP.” \emph{Journal of applied Econometrics} 7 (S1): S61–S82.
#' @references Hamilton, James D. 1989. “A new approach to the economic analysis of nonstationary time series and the business cycle.” \emph{Econometrica} 57 (2): 357–384.
"hamilton84GNP"
