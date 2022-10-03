#' @title Kasaharra 2018 US GDP data
#' 
#' @format 
#' \describe{ This data is used in Kasahara & Shimotsu (2018) working paper. The series ranges from 
#' 1960Q1 to 2014Q4.
#'   \item{DATE}{Vector of dates}
#'   \item{GNP_logdiff}{log difference of US GDP series}
#' }
#' @source \url{https://arxiv.org/abs/1801.06862}
#' @references Kasahara, H., & Shimotsu, K. (2018). Testing the number of regimes in Markov regime 
#' switching models. arXiv preprint arXiv:1801.06862.
"kasshi14GDP"
# -----------------------------------------------------------------------------
#' @title Carasco 2010 GNP
#' 
#' @format 
#' \describe{This dataset is the extension of the GNP series used in CHP (2014). Specifically, 
#' the series \emph{GNP_logdiff} from this dataset is used. This series ranges from 
#' 1951Q2 to 2010Q4. 
#'   \item{DATE}{Vector of dates}
#'   \item{GNP}{US GNP series}
#'   \item{GNP_logdiff}{log difference of US GNP series}
#' }
#' @source \url{https://www.econometricsociety.org/content/supplement-optimal-test-markov-switching-parameters}
#' 
#' @references Carrasco, Marine, Liang Hu, and Werner Ploberger. 2014. 
#' “Optimal test for Markov switch- ing parameters.” \emph{Econometrica} 
#' 
#' \bold{82 (2)}: 765–784.
#' 
"carhuplo10GNP"
# -----------------------------------------------------------------------------
#' @title Hamilton 1984 GNP
#' 
#' @format 
#' \describe{This data set is used in Hansen (1992) to test the US GNP model proposed by 
#' Hamilton (1989). Secifically, the series \emph{GNP_logdiff}. This series ranges from 
#' 1951Q2 to 1984Q4. 
#'   \item{DATE}{Vector of dates}
#'   \item{GNP}{US GNP seies}
#'   \item{GNP_logdiff}{US GNP log difference}
#' }
#' @source \url{https://www.ssc.wisc.edu/~bhansen/progs/jae_92.html}
#' 
#' @references Hansen, B. E. (1992). The likelihood ratio test under nonstandard 
#' conditions: testing the Markov switching model of GNP. Journal of applied Econometrics, 
#' 7(S1), S61-S82.
#' @references Hamilton, J. D. (1989). A new approach to the economic analysis of 
#' nonstationary time series and the business cycle. Econometrica: Journal of the 
#' Econometric Society, 357-384.
"hamilton84GNP"