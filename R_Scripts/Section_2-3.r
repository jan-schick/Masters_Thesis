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
### Impulse Response Function ###
#################################


### Generating a VAR(2) -------------------------------------------------------

set.seed(5)
n <- 30
y <- matrix(vector(length = (n+20)*2, mode = "numeric"),
            ncol = 2)
A1 <- matrix(c(0.5, 0.2,
               0, 0), 
             nrow = 2, byrow = T)
A2 <- matrix(c(0, 0,
               0.95, 0), 
             nrow = 2, byrow = T)
for(i in 3:(n+20)){
  y[i, ] <- A1 %*% y[i-1, ] + A2 %*% y[i-2, ] + c(rnorm(1), rnorm(1, 0, 0.1))
}
y <- as.ts(y[-1:-20, ])


### VAR analysis --------------------------------------------------------------

library(vars)

var.1 <- VAR(y, p = 2, type = "none")
summary(var.1)


### Impulse response analysis  ------------------------------------------------

library(vars)
var.1.irf   <- irf(var.1, ortho = F, impulse = "Series.2")
var.1.irf.o <- irf(var.1, ortho = T, impulse = "Series.2")

var.2 <- VAR(y[, 2:1], type = "none", p = 2)
var.2.irf   <- irf(var.2, ortho = F, impulse = "Series.2")
var.2.irf.o <- irf(var.2, ortho = T, impulse = "Series.2")


plot(var.1.irf)
plot(var.1.irf.o)
plot(var.2.irf)
plot(var.2.irf.o)
# fig. 2.3.1

