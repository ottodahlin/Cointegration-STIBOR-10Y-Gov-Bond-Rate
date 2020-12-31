
#############################################
# Cointegration - Otto Dahlin
#############################################

library(writexl)
library(ggplot2)
library(dplyr)
library(AER)
library(lmtest)
library(tseries)
library(urca)
library(dynlm)
library(sandwich)
library(readxl)
library(forecast)
library(xts)
library(vars)
library(zoo)
library(timeSeries)
library(quantmod)
library(mFilter)
library(seasonal)
library(lubridate)
library(CARS)
library(car)
library(fUnitRoots)


data <- read_excel("Stiborand10Yrates.xlsx")
data

Stibor1w<- data[,"STIBOR1W"]
Stibor1w

Gov10Y<- data[,"10Y"]
Gov10Y

# STIBOR 1W converting to time series
Stibor1w.ts <- ts(Stibor1w, frequency=12, start=c(1987,09), end =c(2020,07))
Stibor1w.ts
is.ts(Stibor1w.ts)
ts.plot(Stibor1w.ts, main="STIBOR 1W series", ylab="Percentage %")


# 10Y Gov Bond converting to time series
Gov10Y.ts <- ts(Gov10Y, frequency=12, start=c(1987,09), end =c(2020,07))
Gov10Y.ts
is.ts(Gov10Y.ts)
ts.plot(Gov10Y.ts, main="10Y Swedish Government Bond rate", ylab="Percentage %")


########################################################################
# Plot the time series.  Could these potentially be cointegrated
########################################################################

merged<- cbind(Gov10Y.ts, Stibor1w.ts)
plot(merged, main="Merged 10Y Gov Bond and STIBOR 1W rates")

# plotting both two series in one single figure:
ts.plot(Gov10Y.ts, Stibor1w.ts, main="Merged 10Y Gov Bond and STIBOR 1W rates",
        gpars=list(xlab="Time", col=c("black", "red")))
legend("topright", legend = c("10Y Gov Bond", "STIOBOR 1W"), 
       col = c("black", "red"), lty = 1, cex=0.8)

# Visually both series tend to decrease/fall with time.

# Checking for stationarity and Unit root tests:

#Starting with KPSS - Stationarity test - SERIES IN LEVELS

# KPSS-test Stationarity Test <Gov10Y.ts> in levels
ur.kpss(Gov10Y.ts, type = "tau")@teststat
ur.kpss(Gov10Y.ts, type ="tau")@cval
# 10 Y gov bond series is NOT stationary
# 0.66 > 0.146 at 5pct level - NOT Stationary.

# First Difference of Gov10Y.ts levels series
d.Gov10Y.ts <- diff(Gov10Y.ts)
ur.kpss(d.Gov10Y.ts, type = "tau")@teststat
ur.kpss(d.Gov10Y.ts, type ="tau")@cval
# 10Y Gov Bond rate is Stationary now in first diff transformation! 
# It is I(1) in other words since 0.02050925<0.146 at 5 pct level.

# KPSS-test Stationarity Test <Stibor1w.ts> in levels
ur.kpss(Stibor1w.ts, type = "tau")@teststat
ur.kpss(Stibor1w.ts, type ="tau")@cval
# STIBOR 1W series is NOT stationary either since 0.6163559>0.146 at 5 pct

# First Difference of Stibor1w.ts levels series
d.Stibor1w.ts<- diff(Stibor1w.ts)
ur.kpss(d.Stibor1w.ts, type = "tau")@teststat
ur.kpss(d.Stibor1w.ts, type ="tau")@cval
# STIBOR 1W rate is Stationary now in first diff transformation! 
# It is I(1) in other words as well since now 0.01688493<0.146  at 5 pct level.

# So, both series are I(1) from KPSS test.
# They are both seperately non-stationary once transformed.


# **** DF Test BELOW (using Urca-package)


# using "urca-package" because compared to other
# Unit root test/packages "urca" has automatic lag selection criterion.
# Therefore "urca" has certain benefits over other Unit root packages.


# Below follows Unit-Root Tests for the two series.
# In total four seperate unit root tests will be done:
# 2 Unit root tests on variables in levels.
# 2 unit root tests on variables in first differences.

# DF test in levels on Gov 10Y rate
Gov10y.unitroot.levels.none <-ur.df(Gov10Y.ts, type = "none", selectlags="AIC")
summary(Gov10y.unitroot.levels)

Gov10y.unitroot.levels.drift <-ur.df(Gov10Y.ts, type = "drift", selectlags="AIC")
summary(Gov10y.unitroot.levels.drift)


# DF test in first difference on Gov 10Y rate
Gov10y.unitroot.firstdiff <-ur.df(d.Gov10Y.ts, type = "drift", selectlags="AIC")
summary(Gov10y.unitroot.firstdiff)

####################   ###########################
# DF test in levels on STIBOR
STIBOR.unitroot.levels.none <-ur.df(Stibor1w.ts, type = "none", selectlags="AIC")
summary(STIBOR.unitroot.levels.none)

STIBOR.unitroot.levels.drift <-ur.df(Stibor1w.ts, type = "drift", selectlags="AIC")
summary(STIBOR.unitroot.levels.drift)


# DF test in first difference on Gov 10Y rate
STIBOR.unitroot.firstdiff.drift <-ur.df(d.Stibor1w.ts, type = "drift", selectlags="AIC")
summary(STIBOR.unitroot.firstdiff.drift)


############################################################
# Can test using no constant, constant (drift) and trend

###################################################
# CONSTANT = NO DRIFT, NO TREND
# DRIFT = With DRIFT, NO TREND
# TREND:with DRIFT and TREND
###################################################

# If our retrieved test-statistics is greater than (in absolute value) 
# to tabulated value, we then can reject the Null Hypothesis

# Hypothesis:
# H0: Unit root in time series/data series

# We will select 5% critical value going forward.


##############################################################
# UNIT ROOT TEST: 10Y Gov Bonb rate in Levels:
##############################################################
# Lags can be chosen by AIC or BIC if we desire.
Gov10y.unitroot.none <-ur.df(Gov10Y.ts, type = "none", selectlags="AIC")
summary(Gov10y.unitroot.none)
# Can reject H0 at 5% level with constant term

Gov10y.unitroot.drift<- ur.df(Gov10Y.ts, type = "drift", selectlags = "AIC")
summary(Gov10y.unitroot.drift)

# TEST STATISTICS:
#tau2: -1.3143
#phi1: 2.2707
# Cannnot reject H0 at 5% level.

Gov10y.unitroot.trend <- ur.df(Gov10Y.ts, type = "trend", selectlags = "AIC")
summary(Gov10y.unitroot.trend) 

# TEST STATISTICS:
# tau3: -3.2253
# Phi2: 4.4754
# Phi3: 5.2779 

# Cannot reject H0 at 5% level when in Levels.


##############################################################
# UNIT ROOT TEST: 10Y Gov Bonb rate in differenced data
##############################################################
# d.Gov10Y.ts = differenced data.

# Lags can be chosen by AIC or BIC if we desire.
diff.Gov10y.unitroot.none <-ur.df(d.Gov10Y.ts, type = "none", selectlags="AIC")
summary(diff.Gov10y.unitroot.none)
# Can reject H0 at 5% level with constant term. 
# Clearly the test-statistic is highly extreme and exceeds 
# 5% critical value of -1.95

d.Gov10y.unitroot.drift<- ur.df(d.Gov10Y.ts, type = "drift", selectlags = "AIC")
summary(d.Gov10y.unitroot.drift)

# TEST STATISTICS:
#tau2: -11.0845
#phi1: 61.4358
# H0 can be rejected at 5% level.

d.Gov10y.unitroot.trend <- ur.df(d.Gov10Y.ts, type = "trend", selectlags = "AIC")
summary(d.Gov10y.unitroot.trend) 

# TEST STATISTICS:
# tau3: -11.0728
# Phi2: 40.8779
# Phi3: 61.3138 

# Can reject H0 at 5% level!


##############################################################
# UNIT ROOT TEST: STIBOR 1W rate in levels data
##############################################################

# Lags can be chosen by AIC or BIC if we desire.
Stibor1w.unitroot.none <-ur.df(Stibor1w.ts, type = "none", selectlags="AIC")
summary(Stibor1w.unitroot.none)
# Can reject H0 at 5% level with constant term at 5% level.

Stibor1w.unitroot.drift<- ur.df(Stibor1w.ts, type = "drift", selectlags = "AIC")
summary(Stibor1w.unitroot.drift)

# TEST STATISTICS:
#tau2: -3.3414
#phi1: 5.6273
# Can reject at 5% level but not at 1% level.

Stibor1w.unitroot.trend <- ur.df(Stibor1w.ts, type = "trend", selectlags = "AIC")
summary(Stibor1w.unitroot.trend) 

# TEST STATISTICS:
# tau3: -6.7252
# Phi2: 15.1084
# Phi3: 22.614  

# Cant reject H0 at 5% level 


##############################################################
# UNIT ROOT TEST: STIBOR 1W rate in differenced data
##############################################################
# d.Stibor1w.ts  = differenced data

# Lags can be chosen by AIC or BIC if we desire.
d.Stibor1w.unitroot.none <-ur.df(d.Stibor1w.ts, type = "none", selectlags="AIC")
summary(Stibor1w.unitroot.none)
# Can reject H0 at 5% level with constant term at 5% level.

d.Stibor1w.unitroot.drift<- ur.df(d.Stibor1w.ts, type = "drift", selectlags = "AIC")
summary(d.Stibor1w.unitroot.drift)

# TEST STATISTICS:
#tau2: -21.0301
#phi1: 221.1322 
# Can reject at 5% level

d.Stibor1w.unitroot.trend <- ur.df(d.Stibor1w.ts, type = "trend", selectlags = "AIC")
summary(d.Stibor1w.unitroot.trend) 

# TEST STATISTICS:
# tau3: -21.003
# Phi2: 147.0427
# Phi3: 220.564   

# Cant reject H0 at 5% level 



####################################################

# UNIT ROOT TESTS using "adfTest" package

# adfTest is in fUnitRoots package.

# can test using no constant, constant or constant and trend

################################################################
# nc = no constant nor time trend
# c = constant but no time trend
# ct = constant and time trend.

# Default is "c" i.e. Constant but no time trend.
# Default lag is 1

# H0: existence of unit root
# IF p-value is less that a=0.05 we then reject H0 of unit root
################################################################


##############################################################
# UNIT ROOT TEST (adfTest): 10Y Gov Bonb rate in Levels:
#############################################################

adfTest(Gov10Y.ts, type = "nc", lags=12) # CAN Reject H0 
adfTest(Gov10Y.ts, type = "c", lags=12) # CANNOT reject H0
adfTest(Gov10Y.ts, type = "ct", lags=12) # CAN reject H0 / 
# (when lag=12 then CANNOT reject H0)

#(Tried increasing lag order to 12 due to monthly data, then "ct" CANNOT reject H0)
# Different lag orders may yield different results....
# have tried different lag lengths, received same results
# despite altered lag lengths.

# ***DIfferenced data of Gov 10y bond rate

adfTest(d.Gov10Y.ts, type = "nc") # CAN Reject H0 of unit root
adfTest(d.Gov10Y.ts, type = "c") # CAN reject H0 of unit root
adfTest(d.Gov10Y.ts, type = "ct") # CAN reject H0 of unit root

# Can reject in all cases when differenced data and despite components used.
# can reject in all cases despite lag length ofcourse.


###########################################################################

##############################################################
# UNIT ROOT TEST (adfTest): STIBOR1W rate in Levels:
#############################################################

adfTest(Stibor1w.ts, type = "nc", lags=4) # CANT Reject H0 
adfTest(Stibor1w.ts, type = "c", lags=4) # CANT reject H0
adfTest(Stibor1w.ts, type = "ct", lags=4) # CAN reject H0 
# Can reject in all above cases in levels STIBOR 1W data
# P-values are all less than a=0.05

# DIfferenced data of Stibor1w data series

adfTest(d.Stibor1w.ts, type = "nc", lags=12) # CAN Reject H0 of unit root
adfTest(d.Stibor1w.ts, type = "c") # CAN reject H0 of unit root
adfTest(d.Stibor1w.ts, type = "ct") # CAN reject H0 of unit root

# Can reject in all cases when differenced data and despite "types".
# Can reject in all above cases in levels STIBOR 1W differenced data
# P-values are all less than a=0.05

###########################################################################


#############################################################################
# Estimating the cointegrating rank of the system by a sequence of tests:
##############################################################################


# Cointegration 
# Need to bind the two series into a system.

binded.series.levels <- cbind(Stibor1w.ts, Gov10Y.ts)
binded.series.levels

# Lag selection criterion (just checking what lag length it suggests)

lagselect.const <- VARselect(binded.series.levels, lag.max=12, type="const")
lagselect.const$selection

lagselect.both <- VARselect(binded.series.levels, lag.max=12, type="both")
lagselect.both$selection

lagselect.trend <- VARselect(binded.series.levels, lag.max=12, type="trend")
lagselect.trend$selection

lagselect.none <- VARselect(binded.series.levels, lag.max=12, type="none")
lagselect.none$selection

# 4 most common otherwise BIC: 2

# Johansen Testing TRACE test:


trace.test.lag2 <- ca.jo(binded.series.levels, type ="trace", spec = "longrun", K = 2, ecdet="none")
# have tried both const and none... same results received.
summary(trace.test.lag2)

# when "r=0" then we see that 64,37 > 17.95. We therefore reject the
#Hypothesis (H0) of r=0
# implying that there should be AT LEAST 1 cointegrated relation between
# where as when increasing ranks to "r<=1" we CANNOT reject the H0 due to
#4.64 < 9.24 at 5pct level
# STIBOR1W and 10Y Gov bond rate when both data series are in levels.
# Same results received if choice of model formulation was set
#to "transitory" as in Hamilton. 

# Testing Maximum eigen value (benchmarking)

eigen.test.lag2 <- ca.jo(binded.series.levels, type ="eigen",spec = "longrun", 
                         K = 2, ecdet="none")
summary(eigen.test.lag2)

# same results received where hypothesis of r=0 rejected due 
#to (64.63 > 15.67)
# Proceeding with higher ranks when (r<=1), we cannot further 
#reject H0 due to (4.64<9.24) 
# implying that even the Eigen Value test suggest that there is 
#AT LEAST 1 cointegrated
# relationship between STIBOR 1W and 10Y GOV BOND rate when both held in levels.

# CONCLUSION: The cointegrated rank of the system is 1.


#############################################################################
#  Assuming the cointegrating rank is 1 and
# using the Johansen approach, test that (STIBOR1W - 10Y) is a cointegrating
# relation, i.e. the spread between the two rates is stationary
#############################################################################

# Spread:
spread <- (Stibor1w.ts-Gov10Y.ts)
is.ts(spread)
spread.ts <- ts(spread, frequency=12, start=c(1987,09), end =c(2020,07))
spread.ts
is.ts(spread.ts)
ts.plot(spread.ts, main="Spread: STIBOR1W minus 10Y Gov Bond rate", ylab="Percentage %")


trace.test.lag2.johansen <- ca.jo(binded.series.levels, type ="trace",
                                  spec = "longrun", K = 2, ecdet="none")
summary(trace.test.lag2.johansen)

# TESTING:

# Testing if spread is stationary.
B1 <- matrix(c(1, -1), nrow = 2)
Test <- blrtest(z = johansen, H=B1, r=1)

# RESULTS from LR test using urca-package:
summary(Test)

#############################################################################
#  Test the Null hypothesis using a simple unit root test
##############################################################################

# USING EG-Procedure. Tests performed on residuals, for stationarity

adfTest(spread)
s <- 1.000*ts_stibor -1.099679*ts_gvb
adf.test(s)

# EG-procedure: Using the "Urca" package.
# Both inserted series are non-stationary in levels.

reg_longrun = lm(Gov10Y.ts ~ Stibor1w.ts)
reg_longrun
# Extracting residuals
resid_reg_longrun = reg_longrun$residuals # storing residuals

# testing for stationarity on residuals below:
resid.test <- ur.df(resid_reg_longrun, type = "none", selectlags="AIC")
summary(resid.test)

resid.test <- ur.df(resid_reg_longrun, type = "drift", selectlags="AIC")
summary(resid.test)

# -6.9109 < -1.95, we can reject the null hypothesis of unit root.
# No unit root in residuals.

####################################################################
####################################################################
# END
####################################################################
####################################################################
