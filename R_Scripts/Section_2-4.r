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
### Forecast Error Variance Decomposition ###
#############################################


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

var.1 <- VAR(y, 
             p = 2, 
             type = "none")
summary(var.1)


### Forecast Error Variance Decomposition  ------------------------------------

var.1.fevd <- fevd(var.1,
                   n.ahead = 10)

plot(var.1.fevd)
# fig. 2.4.1

