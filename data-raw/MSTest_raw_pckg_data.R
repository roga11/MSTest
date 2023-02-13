#' -----------------------------------------------------------------------------
#' Rodriguez Rondon & Dufour 2022 GNP data
rm(list=ls())
USGNP <- read.table('inst/extdata/GNP22Q3.csv', sep = ',', header = T)
USGNP[2:nrow(USGNP),3] <- diff(log(USGNP[,2]))*100
USGNP <- USGNP[2:nrow(USGNP), c(1,3,2)]
colnames(USGNP) <- c('DATE','GNP_logdiff','GNP')
usethis::use_data(USGNP, overwrite = TRUE)
#' -----------------------------------------------------------------------------
#' Carrasco, Hu and Ploberger 2010 GNP Data
rm(list=ls())
GNPC96 <- read.table('inst/extdata/GNPC96.txt', sep = ',', header = F)
date <- seq(as.Date("1951/1/1"), as.Date("2010/10/1"), "quarter")
GNP_change = 100*diff(log(GNPC96$V1))
chp10GNP <- data.frame(date[2:length(date)],GNP_change,GNPC96$V1[2:length(date)])
colnames(chp10GNP) <- c('DATE','GNP_logdiff','GNP')
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
colnames(hamilton84GNP) <- c('DATE','GNP_logdiff','GNP')
usethis::use_data(hamilton84GNP, overwrite = TRUE)