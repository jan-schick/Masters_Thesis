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
### Author: Jan Schick                    ###
### 1st examiner: Prof. Bernhard Möhring  ###
### 2nd examiner: Prof. Helmut Herwartz   ###
### Advisor:      Kai Husmann             ###
### Date:         20.03.2019              ###
###                                       ###
#############################################


### Used packages -------------------------------------------------------------

library(tseries)
library(vars)
library(urca)
library(gtools)

### Data ----------------------------------------------------------------------

### Please set the working directory to the folder "R_Scripts/Data/" of the 
### data obtained from the CD / GitHub

setwd("")

###############################################################################
###                                                                         ###
###   Data obtained from table 61231-0002 of the GENESIS database           ###
###   https://www-genesis.destatis.de/genesis/online                        ###
###   Retrieved 19.02.2019 / 10:38:39                                       ###
###   (C)opyright Statistisches Bundesamt (Destatis), 2019                  ###
###                                                                         ###
###############################################################################

pi.data <- read.csv("pi_data.csv",
                    skip = 9,
                    row.names = 1,
                    fileEncoding = "utf8")
pi.data <- ts(pi.data,
              start = 1988, 
              frequency = 12)

pi.data.d <- read.csv("pi_data_d.csv",
                      skip = 9,
                      row.names = 1,
                      fileEncoding = "utf8")
pi.data.d <- ts(pi.data.d,
                start = 1988, 
                frequency = 12)

plot(pi.data)
# cp. fig. 4.11


### Storm dummies -------------------------------------------------------------

pi.storm <- read.csv("pi_storm.csv", 
                     row.names = 1, 
                     fileEncoding = "utf8")


###############################################################################
### Order of integration ###
############################


### Original data
apply(pi.data, 2, adf.test)
apply(pi.data, 2, kpss.test)
# cp. tab. 4.15

### Differences
apply(pi.data.d, 2, adf.test)
apply(pi.data.d, 2, kpss.test)
# cp. tab. 4.15


###############################################################################
### VAR analysis ###
####################


### Selecting lag order -------------------------------------------------------

VARselect(pi.data, 
          lag.max = 36, 
          season = 12,
          type = "const",
          exogen = pi.storm)
# cp. tab. 4.16


### Johansen test -------------------------------------------------------------

pi.vecm <- ca.jo(pi.data, 
                 K = 12, 
                 season = 12, 
                 ecdet = "const",
                 dumvar = pi.storm)
summary(pi.vecm)
# cp. tab. 4.17


### VAR estimation ------------------------------------------------------------

pi.var <- VAR(pi.data, 
           p = 12, 
           season = 12,
           type = "const", 
           exogen = pi.storm)
summary(pi.var)
# cp. R Output 4.2, tab. 4.18, fig. 4.12, fig. 4.13


### Residuals -----------------------------------------------------------------

plot(resid(pi.var)[, 1], type = "l")
plot(resid(pi.var)[, 2], type = "l")
plot(resid(pi.var)[, 3], type = "l")
plot(resid(pi.var)[, 4], type = "l")
# cp. fig. 4.14


### Correlation and lm of fitted and real -------------------------------------

for(i in c("Oak", "Beech", "Spruce", "Pine")){
  val <- c(round(cor(fitted(pi.var)[, i], pi.data[-1:-12, i]), 3), 
           round(lm(fitted(pi.var)[, i] ~ 
                      pi.data[-1:-12, i])$coefficients, 3))
  names(val) <- c("Correlation", "Intercept", "Slope")
  print(i)
  print(val)
}


### Display of fitted and real ------------------------------------------------

plot(pi.data[, 1])
lines(fitted(pi.var)[, 1] ~ names(fitted(pi.var)[, 1]), 
      col = "red")
plot(fitted(pi.var)[, 1] ~ pi.data[-1:-12, 1])
# cp. fig. 4.15, Oak

plot(pi.data[, 2])
lines(fitted(pi.var)[, 2] ~ names(fitted(pi.var)[, 2]), 
      col = "red")
plot(fitted(pi.var)[, 2] ~ pi.data[-1:-12, 2])
# cp. fig. 4.15, Beech

plot(pi.data[, 3])
lines(fitted(pi.var)[, 3] ~ names(fitted(pi.var)[, 3]), 
      col = "red")
plot(fitted(pi.var)[, 3] ~ pi.data[-1:-12, 3])
# cp. fig. 4.15, Spruce

plot(pi.data[, 4])
lines(fitted(pi.var)[, 4] ~ names(fitted(pi.var)[, 4]), 
      col = "red")
plot(fitted(pi.var)[, 4] ~ pi.data[-1:-12, 4])
# cp. fig. 4.15, Pine


### Checking model assumptions ------------------------------------------------

### Graphical
plot(serial.test(pi.var, 
                 lags.pt = 12))
# cp. fig. 4.14, fig. 4.16

### Autocorrelation
serial.test(pi.var,
            lags.pt = 12)
# cp. tab. 4.20

### Homoscedasticity
arch.test(pi.var,
          lags.single = 12,
          lags.multi = 12,
          multivariate.only = F)
# cp. tab. 4.20

### Normality
normality.test(pi.var,
               multivariate.only = F)
# cp. tab. 4.20


###############################################################################
### Impulse Responses ###
#########################


### Ordering variables --------------------------------------------------------

pi.var.sort <- VAR(pi.data[, c("Spruce", "Pine", "Beech", "Oak")], 
                   p = 12, 
                   season = 12,
                   type = "const", 
                   exogen = pi.storm)


### IRF estimation ------------------------------------------------------------

pi.irf <- irf(pi.var.sort, 
              n.ahead = 36, 
              ortho = T, 
              runs = 1000)
# This might take a while!

plot(pi.irf)
# cp. fig. 4.18


### Checking stability --------------------------------------------------------

pi.irf.stab <- irf(pi.var.sort, n.ahead = 600, boot = F)
plot(pi.irf.stab)


### Checking other orderings --------------------------------------------------

var.order <- permutations(4, 4, c("Oak", "Beech", "Spruce", "Pine"))
irf.list <- list()

for(i in seq(1, length(var.order[, 1]))){
  pi.var.sort.temp <- VAR(pi.data[, var.order[i, ]], 
                                  p = 12, 
                                  season = 12,
                                  type = "const", 
                                  exogen = pi.storm)
  x <- paste(t(var.order[i, ]), collapse = "-")
  irf.list[[x]] <- irf(pi.var.sort.temp,
                       n.ahead = 36,
                       boot = F,
                       ortho = T)
}

### E.g.:
plot(irf.list[[1]])
# cp. fig. 4.18


###############################################################################
### Forecast Error Variance Decomposition ###
#############################################


pi.fevd <- fevd(pi.var.sort, n.ahead = 36)

plot(pi.fevd)
# cp. fig. 4.18

