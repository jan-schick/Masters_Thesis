#############################################
###                                       ###
###   Multivariate Time Series Analysis   ###
###        as a Tool for Studying         ###
###      Forest Business Development      ###
###          in Westphalia-Lippe          ###
###                                       ###
###          - Master's Thesis -          ###
###         Göttingen University          ###
###          Faculty of Forestry          ###
###                                       ###
###---------------------------------------###
###                                       ###
### Author:       Jan Schick              ###
### 1st examiner: Prof. Bernhard Möhring  ###
### 2nd examiner: Prof. Helmut Herwartz   ###
### Advisor:      Kai Husmann             ###
### Date:         20.03.2019              ###
###                                       ###
#############################################


###############################################################################
### White noise and AR(3) ###
#############################


### AR-process ----------------------------------------------------------------

y <- vector(length = 520, mode = "numeric")

set.seed(13)
for(i in seq(4, 520)){
  y[i] <- 0.9 * y[i-1] -0.5 * y[i-2] + 0.05 * y[i-3] + rnorm(1, 0, 1)
}
y <- as.ts(y[-1:-20])

plot(y)
# see fig. 1.1.1


### white noise ---------------------------------------------------------------

set.seed(123)
wn <- as.ts(rnorm(500, 0, 3))

plot(wn)
# see fig. 1.1.1


###############################################################################
### Stationarity ###
####################


### Non-Stationary data -------------------------------------------------------

x <- vector(length = 520, mode = "numeric")

set.seed(1337)
for(i in seq(2, 520)){
  x[i] <- 1.005 * x[i-1] + rnorm(1, 0, 3)
}
x <- as.ts(x[-1:-20])

plot(x)
# see fig. 1.3.1


### first order differences ---------------------------------------------------

x.diff <- diff(x)

plot(x.diff)
# see fig. 1.3.1


### testing for stationarity --------------------------------------------------

library(tseries)

adf.test(x)

adf.test(x.diff)

