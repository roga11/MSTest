#' -----------------------------------------------------------------------------
#' Rodriguez-Rondon & Dufour 2024 Real GDP data
rm(list=ls())
USRGDP <- read.table('inst/extdata/US_RGDP_19470401_20240401.csv', sep = ',', header = T)
colnames(USRGDP) <- c('Date','RGDP','RGDP_gr')
usethis::use_data(USRGDP, overwrite = TRUE)
#' -----------------------------------------------------------------------------
#' Rodriguez-Rondon & Dufour 2024 GNP data
rm(list=ls())
USGNP <- read.table('inst/extdata/US_GNP_19470401_20240401.csv', sep = ',', header = T)
colnames(USGNP) <- c('Date','GNP','GNP_gr')
usethis::use_data(USGNP, overwrite = TRUE)
#' -----------------------------------------------------------------------------
#' Carrasco, Hu and Ploberger 2010 GNP Data
rm(list=ls())
GNPC96 <- read.table('inst/extdata/GNPC96.txt', sep = ',', header = F)
date <- seq(as.Date("1951/1/1"), as.Date("2010/10/1"), "quarter")
GNP_change = 100*diff(log(GNPC96$V1))
chp10GNP <- data.frame(date[2:length(date)],GNP_change,GNPC96$V1[2:length(date)])
chp10GNP <- chp10GNP[,c(1,3,2)]
colnames(chp10GNP) <- c('DATE','GNP','GNP_gr')
usethis::use_data(chp10GNP, overwrite = TRUE)
#' -----------------------------------------------------------------------------
#' Hamilton 1984 GNP Data
rm(list=ls())
hamilton84GNP <- read.table('inst/extdata/GNP82.DAT', header = F)
date <- seq(as.Date("1951/1/1"), as.Date("1984/10/1"), "quarter")
GNP <- matrix(0, dim(hamilton84GNP)[1]*dim(hamilton84GNP)[2], 1 )
xk <- 1
for (xi in 1:dim(hamilton84GNP)[1] ){
  for (xj in 1:dim(hamilton84GNP)[2]){
    GNP[xk] <- hamilton84GNP[xi,xj]
    xk <- xk +1
  }
}
date <- date[2:length(date)]
Tsize <- length(GNP)
Y <- matrix(0, Tsize-1, 1)
Y <- 100*( log(GNP[2:Tsize]) -log(GNP[1:(Tsize-1)]) )
hamilton84GNP <- data.frame(date,Y,GNP[2:Tsize])
hamilton84GNP <- hamilton84GNP[,c(1,3,2)]
colnames(hamilton84GNP) <- c('DATE','GNP','GNP_gr')
usethis::use_data(hamilton84GNP, overwrite = TRUE)