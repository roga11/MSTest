#' @title US Real GDP data 1947Q2 - 2024Q2
#' 
#' @format 
#' \describe{This data is used in Rodriguez-Rondon & Dufour (2024). The series ranges from 
#' 1947Q2 to 2024Q2.
#'   \item{Date}{Vector of dates}
#'   \item{RGDP}{US Real GDP series}
#'   \item{RGDP_gr}{log difference of US Real GDP series}
#' }
#' @source \url{https://fred.stlouisfed.org/series/GDPC1}
#' @references U.S. Bureau of Economic Analysis, Real Gross Domestic Product [GDPC1], retrieved from FRED, Federal Reserve Bank of St. Louis; \emph{\url{https://fred.stlouisfed.org/series/GDPC1}}, September, 2024.
#' @references Rodriguez-Rondon, Gabriel and Jean-Marie Dufour. 2024. “Monte Carlo Likelihood Ratio Tests for Markov Switching Models.” \emph{Unpublished manuscript}.
"USRGDP"
# -----------------------------------------------------------------------------
#' @title US GNP data 1947Q2 - 2024Q2
#' 
#' @format 
#' \describe{This data is used in Rodriguez-Rondon & Dufour (2024). The series ranges from 
#' 1947Q2 to 2024Q2.
#'   \item{Date}{Vector of dates}
#'   \item{GNP}{US GNP series}
#'   \item{GNP_gr}{log difference of US GNP series}
#' }
#' @source \url{https://fred.stlouisfed.org/graph/?g=UKDQ}
#' @references U.S. Bureau of Economic Analysis, Gross National Product [GNP], retrieved from FRED, Federal Reserve Bank of St. Louis; \emph{\url{https://fred.stlouisfed.org/series/GNP}}, September, 2024.
#' @references Rodriguez-Rondon, Gabriel and Jean-Marie Dufour. 2024. “Monte Carlo Likelihood Ratio Tests for Markov Switching Models.” \emph{Unpublished manuscript}.
"USGNP"
# -----------------------------------------------------------------------------
#' @title Carrasco, Hu, & Ploberger 2010 GNP data
#' 
#' @format 
#' \describe{This data is the extension of the GNP series used in CHP (2014), Econometrica. This series ranges from 
#' 1951Q2 to 2010Q4. 
#'   \item{Date}{Vector of dates}
#'   \item{GNP}{US GNP series}
#'   \item{GNP_gr}{log difference of US GNP series}
#' }
#' @source \url{https://www.econometricsociety.org/content/supplement-optimal-test-markov-switching-parameters}
#' 
#' @references Carrasco, Marine, Liang Hu, and Werner Ploberger. 2014. “Optimal test for Markov switching parameters.” \emph{Econometrica} 82 (2): 765–784.
"chp10GNP"
# -----------------------------------------------------------------------------
#' @title Hamilton 1984 & Hansen 1992 GNP data
#' 
#' @format 
#' \describe{This data set is used in Hansen (1992) to test the US GNP model proposed by 
#' Hamilton (1989). This series ranges from 1951Q2 to 1984Q4. 
#'   \item{Date}{Vector of dates}
#'   \item{GNP}{US GNP series}
#'   \item{GNP_gr}{US GNP log difference}
#' }
#' @source \url{https://www.ssc.wisc.edu/~bhansen/progs/jae_92.html}
#' 
#' @references Hansen, Bruce E. 1992. “The likelihood ratio test under nonstandard conditions: testing the Markov switching model of GNP.” \emph{Journal of applied Econometrics} 7 (S1): S61–S82.
#' @references Hamilton, James D. 1989. “A new approach to the economic analysis of nonstationary time series and the business cycle.” \emph{Econometrica} 57 (2): 357–384.
"hamilton84GNP"
