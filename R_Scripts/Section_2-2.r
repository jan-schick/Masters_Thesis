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


###############################################################################
### Generating and analyzing cointegrated data ###
##################################################


### Generating data -----------------------------------------------------------

set.seed(2048)
yg  <- cumsum(rnorm(100, 0, 3)) + 60 
ye1 <- 1.2 * yg + rnorm(100, 0, 10)
ye2 <- 0.9 * yg + rnorm(100, 0, 8)
y <- data.frame(ye1, ye2, yg)

plot(y[, 3], 
     type = "l",
     ylim = c(30, 130))
lines(y[, 1],
      col = "red")
lines(y[, 2],
      col = "blue")
# fig. 2.2.2


### VECM analysis -------------------------------------------------------------

library(urca)
vecm <- ca.jo(y,
              type = "eigen",
              ecdet = "none",
              K = 2,
              spec = "transitory")
summary(vecm)

(vecm.r2 <- cajorls(vecm, r = 2))
summary(vecm.r2$rlm)

beta <- vecm.r2$beta

### Cointegration relations ---------------------------------------------------

EC <- as.matrix(y) %*% beta
plot(EC[, 1], type = "l")
plot(EC[, 2], type = "l")
# fig. 2.2.3

R1t <- vecm@RK
EC.R1t <- R1t %*% beta
plot(EC.R1t[, 1], type = "l")
plot(EC.R1t[, 2], type = "l")
# fig. 2.2.3

