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
### Generating and analyzing VAR processes ###
##############################################


### Interdependent Time Series ------------------------------------------------

set.seed(13)
d.m <- matrix(rep(0, 240),
              ncol = 2, 
              dimnames = list(NULL, c("do", "dw")))
A1 <- matrix(c(0.9, 0,
               0.5, 0.8),
             nrow = 2, byrow = F)
for(i in seq(2, length(d.m[,1]))){
  d.m[i, ] <- d.m[i-1, ] %*% A1 + c(rnorm(1), rnorm(1, 0, 0.1))
}
d.m <- as.ts(d.m[-20:-1, ])

plot(d.m[, 2], col = "red")
lines(d.m[, 1])
# fig. 2.1.1

ar(d.m[,1])
ar(d.m[,2])

### Analyzing a VAR(1) --------------------------------------------------------

library(vars)

(var.model <- VAR(d.m, type = "none"))
Acoef(var.model)


###############################################################################
### Forecasting ###
###################


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


### Forecast ------------------------------------------------------------------

library(vars)
var.1 <- VAR(y, p = 2, type = "none")
summary(var.1)
var.1.pred <- predict(var.1, n.ahead = 5)

plot(y[, 1], 
     xlim = c(0, 35))
lines(y[, 2], 
      col = "red")
lines(var.1.pred$fcst$Series.1[, 1] ~ c(31:35), 
      lty = 2)
lines(var.1.pred$fcst$Series.2[, 1] ~ c(31:35), 
      lty = 2, 
      col = "red")
# fig. 2.1.2


###############################################################################
### Restricted VAR ###
######################


### Restricting VAR  ----------------------------------------------------------

Acoef(restrict(var.1, method = "ser", thresh = 2.0))
res.mat <- matrix(c(1, 0, 0, 0,
                    0, 0, 1, 0),
                  nrow = 2, byrow = T)
Acoef(restrict(var.1, method = "manual", resmat = res.mat))


###############################################################################
### Estimation within vars ###
##############################


### Estimation in VAR() -------------------------------------------------------

lm(y[3:30, 1] ~ y[2:29, 1] + y[2:29, 2] + 
     y[1:28, 1] + y[1:28, 2] - 1)$coefficients
var.1$varresult$Series.1$coefficients
