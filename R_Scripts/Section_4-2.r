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


### Data ----------------------------------------------------------------------

### Please set the working directory to the folder "R_Scripts/Data/" of the 
### data obtained from the CD / GitHub

setwd("")

###############################################################################
###                                                                         ###
###   Logging obtained from "Privatwald-Betriebsverlgeich Westfalen-        ###
###   Lippe" (BVGL)                                                         ###
###   Göttingen University, Department of Forest Economics                  ###
###   See Master's thesis for more details                                  ###
###                                                                         ###
###   -------------------------------------------------------------------   ###
###                                                                         ###
###   Wood price index obtained from table 61231-0003 of the GENESIS        ###
###   database                                                              ###
###   https://www-genesis.destatis.de/genesis/online                        ###
###   Retrieved 14.02.2019 / 16:27:54                                       ###
###   (C)opyright Statistisches Bundesamt (Destatis), 2019                  ###
###                                                                         ###
###   -------------------------------------------------------------------   ###
###                                                                         ###
###                     !!!IMPORTANT NOTES!!!                               ###
###   (i)  The data of the BVGL is highly sensitive. Therefore a            ###
###   randomized component was added to the data used here. For details,    ###
###   see the Master's thesis mentioned above.                              ###
###                                                                         ###
###############################################################################

warning("Data of the BVGL is distorted! See notes above!")

### Ring 1
pt.data.r1.d <- read.csv("pt_data_r1_d.csv",
                         fileEncoding = "utf8")
pt.data.r1.d <- ts(pt.data.r1.d,
                   start = 1971)

### Ring 2
pt.data.r2.d <- read.csv("pt_data_r2_d.csv",
                         fileEncoding = "utf8")
pt.data.r2.d <- ts(pt.data.r2.d,
                   start = 1971)

### Ring 3
pt.data.r3.d <- read.csv("pt_data_r3_d.csv",
                         fileEncoding = "utf8")
pt.data.r3.d <- ts(pt.data.r3.d,
                   start = 1971)


### Differenced wood price index and storm dummies ----------------------------

pt.exo.d <- read.csv("pt_exo_d.csv",
                     fileEncoding = "utf8",
                     skip = 9)
pt.exo.d <- ts(pt.exo.d,
               start = 1971)


###############################################################################
### Order of integration ###
############################


###-------------------------------------------------------------------------###
###   IMPORTANT: As mentioned, the BVGL data was changed due to its         ###
###   sensitivity. Thus all results will differ from the Master's thesis!   ###
###-------------------------------------------------------------------------###


### Ring 1
plot(pt.data.r1.d)
# cp. fig. 4.7, Ring 1

adf.test(pt.data.r1.d[, 1])
kpss.test(pt.data.r1.d[, 1])
adf.test(pt.data.r1.d[, 2])
kpss.test(pt.data.r1.d[, 2])
# cp. tab. 4.8, Ring 1

### Ring 2
plot(pt.data.r2.d)
# cp. fig. 4.7, Ring 2

adf.test(pt.data.r2.d[, 1])
kpss.test(pt.data.r2.d[, 1])
adf.test(pt.data.r2.d[, 2])
kpss.test(pt.data.r2.d[, 2])
# cp. tab. 4.8, Ring 2

### Ring 3
plot(pt.data.r3.d)
# cp. fig. 4.7, Ring 3

adf.test(pt.data.r3.d[, 1])
kpss.test(pt.data.r3.d[, 1])
adf.test(pt.data.r3.d[, 2])
kpss.test(pt.data.r3.d[, 2])
# cp. tab. 4.8, Ring 3


###############################################################################
### VAR analysis ###
####################


### Means of the data ---------------------------------------------------------

apply(pt.data.r1.d, 2, mean)
apply(pt.data.r2.d, 2, mean)
apply(pt.data.r3.d, 2, mean)
# cp. tab. 4.0


### Selecting lag order -------------------------------------------------------

### Ring 1
(pt.varselect.r1 <- VARselect(pt.data.r1.d, 
                              lag.max = 5, 
                              type = "none", 
                              exogen = pt.exo.d))

### Ring 2
(pt.varselect.r2 <- VARselect(pt.data.r2.d, 
                              lag.max = 5, 
                              type = "none", 
                              exogen = pt.exo.d))

### Ring 3
(pt.varselect.r3 <- VARselect(pt.data.r3.d, 
                              lag.max = 5, 
                              type = "none", 
                              exogen = pt.exo.d))


### VAR estimation ------------------------------------------------------------

### Ring 1
pt.var.r1 <- VAR(pt.data.r1.d, 
                 p = 3, 
                 type = "none", 
                 exogen = pt.exo.d)
summary(pt.var.r1)
# R Output 4.1

### Ring 2
pt.var.r2 <- VAR(pt.data.r2.d, 
                 p = 2, 
                 type = "none", 
                 exogen = pt.exo.d)
summary(pt.var.r2)
# R Output 8.2

### Ring 3
pt.var.r3 <- VAR(pt.data.r3.d, 
                 p = 3, 
                 type = "none", 
                 exogen = pt.exo.d)
summary(pt.var.r3)
# R Output 8.3


### Correlation of fitted and real --------------------------------------------
# cp. tab. 4.11

### Ring 1
# Deselecting storm years
sel <- seq(-1, (length(fitted(pt.var.r1)[, 1]) - length(pt.data.r1.d[, 1])))
storm.sel <- time(pt.data.r1.d)[sel] %in% c(1990:1992, 2007:2009)

# Beech
round(cor(fitted(pt.var.r1)[, 1][!storm.sel], 
          pt.data.r1.d[sel, 1][!storm.sel]), 3)
round(lm(fitted(pt.var.r1)[, 1][!storm.sel] ~ 
           pt.data.r1.d[sel, 1][!storm.sel])$coefficients, 3)

# Spruce
round(cor(fitted(pt.var.r1)[, 2][!storm.sel], 
          pt.data.r1.d[sel, 2][!storm.sel]), 3)
round(lm(fitted(pt.var.r1)[, 2][!storm.sel] ~ 
           pt.data.r1.d[sel, 2][!storm.sel])$coefficients, 3)

### Ring 2
# Deselecting storm years
sel <- seq(-1, (length(fitted(pt.var.r2)[, 1]) - length(pt.data.r2.d[, 1])))
storm.sel <- time(pt.data.r2.d)[sel] %in% c(1990:1992, 2007:2009)

# Beech
round(cor(fitted(pt.var.r2)[, 1][!storm.sel], 
          pt.data.r2.d[sel, 1][!storm.sel]), 3)
round(lm(fitted(pt.var.r2)[, 1][!storm.sel] ~ 
           pt.data.r2.d[sel, 1][!storm.sel])$coefficients, 3)

# Spruce
round(cor(fitted(pt.var.r2)[, 2][!storm.sel], 
          pt.data.r2.d[sel, 2][!storm.sel]), 3)
round(lm(fitted(pt.var.r2)[, 2][!storm.sel] ~ 
           pt.data.r2.d[sel, 2][!storm.sel])$coefficients, 3)

### Ring 3
# Deselecting storm years
sel <- seq(-1, (length(fitted(pt.var.r3)[, 1]) - length(pt.data.r3.d[, 1])))
storm.sel <- time(pt.data.r3.d)[sel] %in% c(1990:1992, 2007:2009)

# Beech
round(cor(fitted(pt.var.r3)[, 1][!storm.sel], 
          pt.data.r3.d[sel, 1][!storm.sel]), 3)
round(lm(fitted(pt.var.r3)[, 1][!storm.sel] ~ 
           pt.data.r3.d[sel, 1][!storm.sel])$coefficients, 3)

# Spruce
round(cor(fitted(pt.var.r3)[, 2][!storm.sel], 
          pt.data.r3.d[sel, 2][!storm.sel]), 3)
round(lm(fitted(pt.var.r3)[, 2][!storm.sel] ~ 
           pt.data.r3.d[sel, 2][!storm.sel])$coefficients, 3)


###############################################################################
### Checking model assumptions ###
##################################


### Graphical display ---------------------------------------------------------

plot(serial.test(pt.var.r1, type = "ES", lags.bg = 5))
# cp. fig. 4.9
plot(serial.test(pt.var.r2, type = "ES", lags.bg = 5))
# cp. fig. 8.4
plot(serial.test(pt.var.r3, type = "ES", lags.bg = 5))
# cp. fig. 8.5


### Autocorrelation -----------------------------------------------------------

serial.test(pt.var.r1, type = "PT.adjusted", lags.pt = 5)
serial.test(pt.var.r1, type = "ES", lags.bg = 5)
serial.test(pt.var.r2, type = "PT.adjusted", lags.pt = 5)
serial.test(pt.var.r2, type = "ES", lags.bg = 5)
serial.test(pt.var.r3, type = "PT.adjusted", lags.pt = 5)
serial.test(pt.var.r3, type = "ES", lags.bg = 5)
# cp. tab. 4.12


### Homoscedasticity ----------------------------------------------------------

arch.test(pt.var.r1, multivariate.only = F)
arch.test(pt.var.r2, multivariate.only = F)
arch.test(pt.var.r3, multivariate.only = F)
# cp. tab. 4.13

### Normality -----------------------------------------------------------------

### Ring 1
normality.test(pt.var.r1, multivariate.only = F)
# without outliers for the beech
jarque.bera.test(resid(pt.var.r1)[resid(pt.var.r1)[, 1] < 2, 1])

### Ring 2
normality.test(pt.var.r2, multivariate.only = F)
# without outliers for the beech
jarque.bera.test(resid(pt.var.r2)[resid(pt.var.r2)[, 1] < 2.5, 1])
# without outliers for the spruce
jarque.bera.test(resid(pt.var.r2)[resid(pt.var.r2)[, 2] > -5, 2])

### Ring 3
normality.test(pt.var.r3, multivariate.only = F)
# without outliers for the spruce
jarque.bera.test(resid(pt.var.r3)[resid(pt.var.r3)[, 2] < 5, 2])

# cp. tab. 4.14

