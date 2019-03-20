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
###   Revenue obtained from "Privatwald-Betriebsverlgeich Westfalen-        ###
###   Lippe" (BVGL)                                                         ###
###   Göttingen University, Department of Forest Economics                  ###
###   See Master's thesis for more details                                  ###
###                                                                         ###
###   -------------------------------------------------------------------   ###
###                                                                         ###
###   GDP was obtained from statista.com:                                   ###
###   https://de.statista.com/statistik/daten/studie/5023/umfrage/entwicklung-des-bruttoinlandsprodukts-von-nordrhein-westfalen-seit-1970/
###   Retrieved 24.10.2018 / 13:17:24                                       ###
###                                                                         ###
###   -------------------------------------------------------------------   ###
###                                                                         ###
###                     !!!IMPORTANT NOTES!!!                               ###
###   (i)  The data of the BVGL is highly sensitive. Therefore a            ###
###   randomized component was added to the data used here. For details,    ###
###   see the Master's thesis mentioned above.                              ###
###   (ii) As Statista GmbH does not allow transferring the data, it has    ###
###   to be obtained manually by the reader. When doing this, recall,       ###
###   that the GDP was transformed to billions for the analysis.            ###
###                                                                         ###
###############################################################################

gr.data <- read.csv("gr_data.csv",
                    fileEncoding = "utf8",
                    row.names = 1)

if(any(is.na(gr.data))){
  warning(paste0("GDP is missing! As Statista GmbH does not allow transferrin", 
                 "g the data, it has to be retrieved manually from\n\n",
                 "https://de.statista.com/statistik/daten/studie/5023/umfrage/entwicklung-des-bruttoinlandsprodukts-von-nordrhein-westfalen-seit-1970/\n\n",
                 "Please read notes above!"))
}

names(gr.data) <- c("revenue", "gdp.nrw")
gr.data <- ts(gr.data,
              start = min(as.numeric(rownames(gr.data))),
              end = max(as.numeric(rownames(gr.data))))

plot(gr.data[, 2],
     col = "red")
lines(gr.data[, 1])
# cp. fig. 4.1


### Storm dummies -------------------------------------------------------------

gr.storm <- read.csv("gr_storm.csv",
                     fileEncoding = "utf8",
                     row.names = 1)


###############################################################################
### Order of integration ###
############################


###-------------------------------------------------------------------------###
###   IMPORTANT: As mentioned, the BVGL data was changed due to its         ###
###   sensitivity. Thus some results will differ from the Master's thesis!  ###
###-------------------------------------------------------------------------###


### Graphical
plot(diff(gr.data[, "gdp.nrw"]))
plot(diff(gr.data[, "revenue"]))
# cp. fig. 4.1


### For the GDP ---------------------------------------------------------------

adf.test(gr.data[, "gdp.nrw"])
kpss.test(gr.data[, "gdp.nrw"])

### First difference
adf.test(diff(gr.data[, "gdp.nrw"]))
kpss.test(diff(gr.data[, "gdp.nrw"]))

kpss.test(gr.data[, "gdp.nrw"], 
          null = "Trend")


### For the revenue -----------------------------------------------------------

adf.test(gr.data[, "revenue"])
kpss.test(gr.data[, "revenue"])

### First differences
adf.test(diff(gr.data[, "revenue"]))
kpss.test(diff(gr.data[, "revenue"]))

kpss.test(gr.data[, "revenue"], 
          null = "Trend")


###############################################################################
### VECM analysis ###
#####################


### Lag order selection -------------------------------------------------------

VARselect(gr.data, 
          type = "both", 
          lag.max = 5, 
          exogen = gr.storm)


### Johansen test -------------------------------------------------------------

### Trace statistic
gr.vecm.trace <- ca.jo(gr.data, 
                       type = "trace",
                       ecdet = "none",
                       K = 2,
                       dumvar = gr.storm,
                       spec = "transitory")
summary(gr.vecm.trace)
# cp. tab. 4.3

### Maximal eigenvalue statistic
gr.vecm <- ca.jo(gr.data, 
                 type = "eigen",
                 ecdet = "none",
                 K = 2,
                 dumvar = gr.storm,
                 spec = "transitory")
summary(gr.vecm)
# cp. tab. 4.3


### RLS estimation ------------------------------------------------------------

(gr.vecm.r1 <- cajorls(gr.vecm, 
                       r = 1))
summary(gr.vecm.r1$rlm)


### Cointegration relations ---------------------------------------------------

gr.beta <- gr.vecm.r1$beta

### Multiplication with original data
gr.EC <- gr.data %*% gr.beta

plot(gr.EC, type = "l")
# fig. 4.2 a

adf.test(gr.EC)
kpss.test(gr.EC)

### Multiplication with R1t
gr.EC.R1t <- gr.vecm@RK %*% gr.beta

plot(gr.EC.R1t, type = "l")
# fig. 4.2 b

### Stationarity test
adf.test(gr.EC.R1t)
kpss.test(gr.EC.R1t)

### Stationarity test on reduced datasets
p.mat <- matrix(0, nrow = 5, 
                ncol = 2, 
                dimnames = list(c(1:5), c("ADF", "KPSS")))
for(gr.cut in 1:5){ 
  p.mat[gr.cut, 1] <- adf.test(gr.EC.R1t[c(seq(-1, -gr.cut), 
                                           seq(-length(gr.EC.R1t) - 1 + 
                                                 gr.cut, 
                                               -length(gr.EC.R1t)))])$p.value
  p.mat[gr.cut, 2] <- kpss.test(gr.EC.R1t[c(seq(-1, -gr.cut), 
                                            seq(-length(gr.EC.R1t) - 1 + 
                                                  gr.cut, 
                                                -length(gr.EC.R1t)))])$p.value
}
round(p.mat, 3)
# tab. 4.5


### Fitted values -------------------------------------------------------------

gr.vecm.coef <- gr.vecm.r1$rlm$coefficients[-1, ]
gr.alpha <- as.matrix(gr.vecm.r1$rlm$coefficients[1, ])

y.pred <- gr.vecm@ZK %*% (gr.alpha %*% t(gr.beta)) + 
            gr.vecm@Z1 %*% gr.vecm.coef
y.pred <- ts(y.pred,
             start = 2016 - length(y.pred[, 1]) + 1,
             end = 2016)


### Graphical display ---------------------------------------------------------

plot(diff(gr.data[, 1]))
lines(y.pred[, 1], 
      col = "red")
plot(as.numeric(y.pred[, 1]) ~ diff(gr.data[-1, 1]), 
     type = "p")
# cp. fig. 4.3, revenue

plot(diff(gr.data[, 2]))
lines(y.pred[, 2], 
      col = "red")
plot(as.numeric(y.pred[, 2]) ~ diff(gr.data[-1, 2]), 
     type = "p")
# cp. fig. 4.3, GDP


### Residuals -----------------------------------------------------------------

res <- gr.vecm@Z0 - y.pred
plot(res)
hist(res[, 1])
hist(res[, 2])
acf(res[, 1])
acf(res[, 2])
# cp. fig. 4.4

### Stationarity
adf.test(res[, 1])
kpss.test(res[, 1])

adf.test(res[, 2])
kpss.test(res[, 2])
# tab. 4.6

### Autocorrelation
Box.test(res[, 1], type = "Ljung-Box")

Box.test(res[, 2], type = "Ljung-Box")
# tab. 4.6

### Normality
jarque.bera.test(res[, 1])

jarque.bera.test(res[, 2])
# tab. 4.6


### Correlation of fitted and original (differences) --------------------------

### Deselecting storm years for correlation
sel <- !(time(y.pred) %in% c(1990:1992, 2007:2009))

### Revenue
cor(y.pred[sel, 1], diff(gr.data)[c(FALSE, sel), 1])
mean(diff(gr.data)[-1, 1])
mean(y.pred[, 1])

### GDP
cor(y.pred[sel, 2], diff(gr.data)[c(FALSE, sel), 2])
mean(diff(gr.data)[-1, 2])
mean(y.pred[, 2])


###############################################################################
### Conversion to VAR ###
#########################


(gr.vec2var <- vec2var(gr.vecm, r = 1))


### Correlation of fitted and original (level) --------------------------------

### Deselecting storm years for correlation
sel <- !(rownames(fitted(gr.vec2var)) %in% c(1990:1992, 2007:2009))

### Revenue
cor(fitted(gr.vec2var)[sel, "fit of revenue"], 
    gr.data[sel, "revenue"][seq(-1, -2)])
mean(gr.data[, "revenue"][seq(-1, -2)])
mean(fitted(gr.vec2var)[, "fit of revenue"])

### GDP
cor(fitted(gr.vec2var)[sel, "fit of gdp.nrw"], 
    gr.data[sel, "gdp.nrw"][seq(-1, -2)])
mean(gr.data[, "gdp.nrw"][seq(-1, -2)])
mean(fitted(gr.vec2var)[, "fit of gdp.nrw"])


### Graphical display ---------------------------------------------------------

plot(gr.data[, "revenue"][seq(-1, -2)], 
     type = "l")
lines(fitted(gr.vec2var)[, "fit of revenue"],
      col = "red")
plot(fitted(gr.vec2var)[, "fit of revenue"] ~ 
       gr.data[, "revenue"][seq(-1, -2)])
# cp. fig. 4.5, revenue

plot(gr.data[, "gdp.nrw"][seq(-1, -2)], 
     type = "l")
lines(fitted(gr.vec2var)[, "fit of gdp.nrw"],
      col = "red")
plot(fitted(gr.vec2var)[, "fit of gdp.nrw"] ~ 
       gr.data[, "gdp.nrw"][seq(-1, -2)])
# cp. fig. 4.5, GDP

