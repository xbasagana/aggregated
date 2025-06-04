############################################################################################
# Updated R CODE to fit time series models with temporally aggregated response data
#   and temporally disaggregated exposure data. 
# Method described in:
#   "Unbiased temperature-related mortality estimates using weekly and monthly health data:
#     a new method for environmental epidemiology and climate impact studies"
#   Basagaña X and Ballester J
#   Lancet Planetary Health 2024
#   https://www.thelancet.com/journals/lanplh/article/PIIS2542-5196(24)00212-2/fulltext
############################################################################################

# Load libraries

library(splines)
library(FluMoDL) # attrdl function

# Load functions
source("agr_functions_v36_cpp.r")


################################################
# -------------------------------------------- #
# 1. Example with a crossbasis for temperature #
# -------------------------------------------- #
################################################


# Reading the Input Data (time series from Paris)
data <- read.csv(file = "Paris20.txt", sep = ",")

# create Date variable
data$date <- as.Date(paste0(data$year, "-" , data$month, "-", data$day), "%Y-%m-%d")

########################
# Model specifications #
########################

# Exposure-Response Association
var_fun <- "ns"
var_prc <- c(0.10,0.75,0.90)
var_knots <- quantile(data$temp, var_prc, na.rm = TRUE)
boundary_knots <- range(data$temp, na.rm=T)

# Lag parameters
max_lag <- 21
lag_knots <- logknots( max_lag, 3 )

# Percentiles for the Predictions of the Cumulative Exposure-Response Association
pred_prc <- c(seq(0, 1, 0.1), 2:98, seq(99, 100, 0.1))/100
at_temp <- quantile(data$temp, pred_prc)

# CROSSBASIS MATRIX OF TEMPERATURE
cross_basis = crossbasis( data$temp,
                          lag = c(0, max_lag),
                          argvar = list( fun = var_fun, knots = var_knots,
                                         Boundary.knots = boundary_knots),
                          arglag = list( knots = lag_knots ) )

# Degrees of Freedom for Seasonality and Long-Term Trend
df_seas <- 8

# Model Formula (for daily data)
form_mod <- death ~ ns( date, df = round( df_seas * nrow(data) / 365.25 )) + cross_basis


#######################
#######################
# Daily response data #
#######################
#######################

# Fit glm model
model_daily <- glm( form_mod, data, family = "quasipoisson", na.action = "na.exclude" )

# Predict response variable based on the model
pred_daily <- predict(model_daily, type = "response")

# centering point. For now, at the median temperature = 17ºC
cen <- 17

# Predictions of crossbasis associations 
crosspred_daily <- crosspred( cross_basis, model_daily, cen = cen, at = at_temp, bylag = 1)

# Minimum mortality temperature (MMT) based on the model
mmt_daily <- crosspred_daily$predvar[which.min(crosspred_daily$allRRfit)[1]]

# Predictions of crossbasis associations, centering at the MMT 
crosspred_daily <- crosspred( cross_basis, model_daily, cen = mmt_daily, at = at_temp, bylag = 1)

# Reduced (cumulative) crossbasis predictions 
red_daily <- crossreduce(cross_basis, model_daily, cen = mmt_daily, at = at_temp)

# Attributable number

# Total

an_daily <- attrdl(x = data$temp, basis = cross_basis, cases = data$death, model = model_daily,
                         type = "an", dir = "forw", tot = TRUE, cen = mmt_daily, sim = FALSE)

an_daily_sim <- attrdl(x = data$temp, basis = cross_basis, cases = data$death, model = model_daily,
                   type = "an", dir = "forw", tot = TRUE, cen = mmt_daily, sim = TRUE, nsim = 1000)

# Cold

an_daily_cold <- attrdl(x = data$temp, basis = cross_basis, cases = data$death, model = model_daily,
                   type = "an", dir = "forw", tot = TRUE, cen = mmt_daily, sim = FALSE, range = c(min(data$temp, na.rm=T), mmt_daily))

an_daily_cold_sim <- attrdl(x = data$temp, basis = cross_basis, cases = data$death, model = model_daily,
                        type = "an", dir = "forw", tot = TRUE, cen = mmt_daily, sim = TRUE, nsim = 1000, range = c(min(data$temp, na.rm=T), mmt_daily))

# Heat

an_daily_heat <- attrdl(x = data$temp, basis = cross_basis, cases = data$death, model = model_daily,
                        type = "an", dir = "forw", tot = TRUE, cen = mmt_daily, sim = FALSE, range = c(mmt_daily, max(data$temp, na.rm=T)))

an_daily_heat_sim <- attrdl(x = data$temp, basis = cross_basis, cases = data$death, model = model_daily,
                        type = "an", dir = "forw", tot = TRUE, cen = mmt_daily, sim = TRUE, nsim = 1000, range = c(mmt_daily, max(data$temp, na.rm=T)))


########################
########################
# Weekly response data #
########################
########################

# Need to provide Y, a 2-column matrix with the period indicator and the counts
# Aggregate mortality counts by week - this is the outcome variable used in the analysis 
Y <- aggregate(data$death, by = list(week = data$week), FUN = sum)
names(Y)[2] <- "death"

# Variable that distributes the sum of cases per week equally across all 7 days
#   This variable is used when calculating attributable numbers
Ymean_week <- ave(data$death, week = data$week, FUN = mean)

# Fit model 
model_week <- fit_aggregate_cpp(Y = Y, X = data, CB = cross_basis, formula = form_mod, family = "quasipoisson",
                             tol = 1e-8, maxit = 500, seed = 1234, ntryini = 10, ntryfit = 0) 

# check convergence
model_week$converged

summary(model_week)

# Predict response variable based on the model
pred_week <- predict(model_week)

# centering point. For now, at the median temperature = 17ºC
cen <- 17

# Predictions of crossbasis associations
crosspred_week <- crosspred(cross_basis, model_week, cen = cen, at = at_temp, model.link = "log")

# calculate minimum mortality temperature (MMT) from model predictions 
mmt_week <- crosspred_week$predvar[which.min(crosspred_week$allRRfit)[1]]

# REDUCTION TO OVERALL CUMULATIVE CURVE (ALL LAGS), centering at MMT
red_week <- crossreduce(basis = cross_basis, model = model_week, cen = mmt_week, at = at_temp, model.link = "log")
# adding exponentiated coefficients (RR) to the object with an internal function
red_week <- add_exp(red_week)


# Attributable numbers

# Total

an_week <- attrdl.aggr(x = data$temp, basis = cross_basis, cases = Ymean_week,
                      model = model_week, type = "an", dir = "forw", tot = TRUE,
                      cen = mmt_week, sim = FALSE)

an_week_sim <- attrdl.aggr(x = data$temp, basis = cross_basis, cases = Ymean_week,
                       model = model_week, type = "an", dir = "forw", tot = TRUE,
                       cen = mmt_week, sim = TRUE, nsim = 1000)

# Cold

an_week_cold <- attrdl.aggr(x = data$temp, basis = cross_basis, cases = Ymean_week,
                            model = model_week, type = "an", dir = "forw", tot = TRUE,
                            cen = mmt_week, sim = FALSE, range = c(min(data$temp, na.rm = T), mmt_week))

an_week_cold_sim <- attrdl.aggr(x = data$temp, basis = cross_basis, cases = Ymean_week,
                            model = model_week, type = "an", dir = "forw", tot = TRUE,
                            cen = mmt_week, sim = TRUE, nsim = 1000, range = c(min(data$temp, na.rm = T), mmt_week))

# Heat

an_week_heat <- attrdl.aggr(x = data$temp, basis = cross_basis, cases = Ymean_week,
                            model = model_week, type = "an", dir = "forw", tot = TRUE,
                            cen = mmt_week, sim = FALSE, range = c(mmt_week, max(data$temp, na.rm = T)))

an_week_heat_sim <- attrdl.aggr(x = data$temp, basis = cross_basis, cases = Ymean_week,
                            model = model_week, type = "an", dir = "forw", tot = TRUE,
                            cen = mmt_week, sim = TRUE, nsim = 1000, range = c(mmt_week, max(data$temp, na.rm = T)))


####################################################
####################################################
# Monthly (irregular number of days) response data #
####################################################
####################################################

# Create indicator of year and month
data$ymonth <- data$year * 100 + data$month

# Aggregate data by year and month - this is the outcome variable used in the analysis
Y <- aggregate(data$death, by = list(ymonth = data$ymonth), FUN = sum)
names(Y)[2] <- "death"

# Variable that distributes the sum of cases per month equally across all days
#   This variable is used when calculating attributable numbers
Ymean_month <- ave(data$death, data$ymonth, FUN = mean)


# Fit model 
model_month <- fit_aggregate_cpp(Y = Y, X = data, CB = cross_basis, formula = form_mod, family = "quasipoisson",
                              tol = 1e-8, maxit = 5000, seed = 1234, ntryini = 10, ntryfit = 0) 
 

# check convergence
model_month$converged

summary(model_month)

# centering point. For now, at the median temperature = 17ºC
cen <- 17

# Predictions of crossbasis associations
crosspred_month <- crosspred(cross_basis, model_month, cen = cen, at = at_temp, model.link = "log")

# calculate minimum mortality temperature (MMT) from model predictions 
mmt_month <- crosspred_month$predvar[which.min(crosspred_month$allRRfit)[1]]

# REDUCTION TO OVERALL CUMULATIVE CURVE (ALL LAGS), centering at MMT
red_month <- crossreduce(basis = cross_basis, model = model_month, cen = mmt_month, at = at_temp, model.link = "log")
# adding exponentiated coefficients (RR) to the object with an internal function
red_month <- add_exp(red_month)


# Attributable number

# Total

an_month <- attrdl.aggr(x = data$temp, basis = cross_basis, cases = Ymean_month,
                       model = model_month, type = "an", dir = "forw", tot = TRUE, 
                       cen = mmt_month, sim = FALSE)

an_month_sim <- attrdl.aggr(x = data$temp, basis = cross_basis, cases = Ymean_month,
                        model = model_month, type = "an", dir = "forw", tot = TRUE, 
                        cen = mmt_month, sim = TRUE, nsim = 1000)

# Cold

an_month_cold <- attrdl.aggr(x = data$temp, basis = cross_basis, cases = Ymean_month,
                        model = model_month, type = "an", dir = "forw", tot = TRUE, 
                        cen = mmt_month, sim = FALSE, range = c(min(data$temp, na.rm=T), mmt_month))

an_month_cold_sim <- attrdl.aggr(x = data$temp, basis = cross_basis, cases = Ymean_month,
                             model = model_month, type = "an", dir = "forw", tot = TRUE, 
                             cen = mmt_month, sim = TRUE, nsim = 1000, range = c(min(data$temp, na.rm=T), mmt_month))

# Heat

an_month_heat <- attrdl.aggr(x = data$temp, basis = cross_basis, cases = Ymean_month,
                             model = model_month, type = "an", dir = "forw", tot = TRUE, 
                             cen = mmt_month, sim = FALSE, range = c(mmt_month, max(data$temp, na.rm = T)))

an_month_heat_sim <- attrdl.aggr(x = data$temp, basis = cross_basis, cases = Ymean_month,
                             model = model_month, type = "an", dir = "forw", tot = TRUE, 
                             cen = mmt_month, sim = TRUE, nsim = 1000, range = c(mmt_month, max(data$temp, na.rm = T)))


#################################################
#################################################
# Year, Month and day of the week response data #
#################################################
#################################################

# dow variable
data$dow <- as.POSIXlt(data$date)$wday

# Create indicator of year, month and dow
data$ymonthdow <- data$year * 1000 + data$month * 10 + data$dow

# Aggregate data by year and month - this is the outcome variable used in the analysis
Y <- aggregate(data$death, by = list(ymonthdow = data$ymonthdow), FUN = sum)
names(Y)[2] <- "death"

# Variable that distributes the sum of cases per month equally across all days
#   This variable is used when calculating attributable numbers
Ymean_dow <- ave(data$death, data$ymonthdow, FUN = mean)


# Fit model 
model_dow <- fit_aggregate_cpp(Y = Y, X = data, CB = cross_basis, formula = form_mod, family = "quasipoisson",
                                 tol = 1e-8, maxit = 5000, seed = 1234, ntryini = 10, ntryfit = 0) 


# check convergence
model_dow$converged

summary(model_dow)

# centering point. For now, at the median temperature = 17ºC
cen <- 17

# Predictions of crossbasis associations
crosspred_dow <- crosspred(cross_basis, model_dow, cen = cen, at = at_temp, model.link = "log")

# calculate minimum mortality temperature (MMT) from model predictions 
mmt_dow <- crosspred_dow$predvar[which.min(crosspred_dow$allRRfit)[1]]

# REDUCTION TO OVERALL CUMULATIVE CURVE (ALL LAGS), centering at MMT
red_dow <- crossreduce(basis = cross_basis, model = model_dow, cen = mmt_dow, at = at_temp, model.link = "log")
# adding exponentiated coefficients (RR) to the object with an internal function
red_dow <- add_exp(red_dow)


# Attributable number

# Total

an_dow <- attrdl.aggr(x = data$temp, basis = cross_basis, cases = Ymean_dow,
                        model = model_dow, type = "an", dir = "forw", tot = TRUE, 
                        cen = mmt_dow, sim = FALSE)

an_dow_sim <- attrdl.aggr(x = data$temp, basis = cross_basis, cases = Ymean_dow,
                      model = model_dow, type = "an", dir = "forw", tot = TRUE, 
                      cen = mmt_dow, sim = TRUE, nsim = 1000)

# Cold

an_dow_cold <- attrdl.aggr(x = data$temp, basis = cross_basis, cases = Ymean_dow,
                             model = model_dow, type = "an", dir = "forw", tot = TRUE, 
                             cen = mmt_dow, sim = FALSE, range = c(min(data$temp, na.rm=T), mmt_month))

an_dow_cold_sim <- attrdl.aggr(x = data$temp, basis = cross_basis, cases = Ymean_dow,
                           model = model_dow, type = "an", dir = "forw", tot = TRUE, 
                           cen = mmt_dow, sim = TRUE, nsim = 1000, range = c(min(data$temp, na.rm=T), mmt_month))

# Heat

an_dow_heat <- attrdl.aggr(x = data$temp, basis = cross_basis, cases = Ymean_dow,
                             model = model_dow, type = "an", dir = "forw", tot = TRUE, 
                             cen = mmt_dow, sim = FALSE, range = c(mmt_month, max(data$temp, na.rm = T)))

an_dow_heat_sim <- attrdl.aggr(x = data$temp, basis = cross_basis, cases = Ymean_dow,
                           model = model_dow, type = "an", dir = "forw", tot = TRUE, 
                           cen = mmt_dow, sim = TRUE, nsim = 1000, range = c(mmt_month, max(data$temp, na.rm = T)))

#########################
# Comparison of results #
#########################

# Plot cumulative exposure-response functions

plot(red_daily, xlab = "Temperature",ylab = "RR", ci = "n", lwd = 2)
lines(red_week, ci = "n", col = "red", lwd = 2)
lines(red_month, ci = "n", col = "blue", lwd = 2)
lines(red_dow, ci = "n", col="#44AA99", lwd = 2)
legend("top", c("Daily", "Weekly", "Monthly", "Dow"), lty = 1, lwd = 2, col = c("black", "red", "blue", "#44AA99"))

# Same plot with confidence intervals

plot(red_daily, xlab = "Temperature", ylab = "RR", ci = "area", lwd = 3, ylim = c(1,7),
     ci.arg = list(col = rgb(0, 0, 0, 0.2), border = FALSE, lwd = 1))
lines(red_week, col = "red", ci = "area", lwd = 3,
      ci.arg = list(col = rgb(1, 0, 0, 0.2), border = FALSE, lwd = 1))
lines(red_month, col = "blue", ci = "area", lwd = 3,
      ci.arg = list(col = rgb(0, 0, 1, 0.2), border = FALSE, lwd = 1))
lines(red_dow, col = "#44AA99", ci = "area", lwd = 3,
      ci.arg = list(col = rgb(0.27, 0.67, 0.6, 0.2), border = FALSE, lwd = 1))


# Compare ANs

cat(
  paste0("AN (D|D)  : ", round(an_daily), " (95% CI: ", round(quantile(an_daily_sim, 0.025)), ", ", round(quantile(an_daily_sim, 0.975)), ")\n",
         "AN (W|D)  : ", round(an_week), " (95% CI: ", round(quantile(an_week_sim, 0.025)), ", ", round(quantile(an_week_sim, 0.975)), ")\n",
         "AN (M|D)  : ", round(an_month), " (95% CI: ", round(quantile(an_month_sim, 0.025)), ", ", round(quantile(an_month_sim, 0.975)), ")\n",
         "AN (Dow|D): ", round(an_dow), " (95% CI: ", round(quantile(an_dow_sim, 0.025)), ", ", round(quantile(an_dow_sim, 0.975)), ")\n"
))

cat(
  paste0("AN (Cold, D|D)  : ", round(an_daily_cold), " (95% CI: ", round(quantile(an_daily_cold_sim, 0.025)), ", ", round(quantile(an_daily_cold_sim, 0.975)), ")\n",
         "AN (Cold, W|D)  : ", round(an_week_cold), " (95% CI: ", round(quantile(an_week_cold_sim, 0.025)), ", ", round(quantile(an_week_cold_sim, 0.975)), ")\n",
         "AN (Cold, M|D)  : ", round(an_month_cold), " (95% CI: ", round(quantile(an_month_cold_sim, 0.025)), ", ", round(quantile(an_month_cold_sim, 0.975)), ")\n",
         "AN (Cold, Dow|D): ", round(an_dow_cold), " (95% CI: ", round(quantile(an_dow_cold_sim, 0.025)), ", ", round(quantile(an_dow_cold_sim, 0.975)), ")\n"
  ))

cat(
  paste0("AN (Heat, D|D)  : ", round(an_daily_heat), " (95% CI: ", round(quantile(an_daily_heat_sim, 0.025)), ", ", round(quantile(an_daily_heat_sim, 0.975)), ")\n",
         "AN (Heat, W|D)  : ", round(an_week_heat), " (95% CI: ", round(quantile(an_week_heat_sim, 0.025)), ", ", round(quantile(an_week_heat_sim, 0.975)), ")\n",
         "AN (Heat, M|D)  : ", round(an_month_heat), " (95% CI: ", round(quantile(an_month_heat_sim, 0.025)), ", ", round(quantile(an_month_heat_sim, 0.975)), ")\n",
         "AN (Heat, Dow|D): ", round(an_dow_heat), " (95% CI: ", round(quantile(an_dow_heat_sim, 0.025)), ", ", round(quantile(an_dow_heat_sim, 0.975)), ")\n"
  ))




###################################################################
# --------------------------------------------------------------- #
# 2. Example with a single parameter of interest (air pollution)  #
# --------------------------------------------------------------- #
###################################################################

library(dlnm)  # to use the dataset chicagoNMMAPS

# function to create lags
lagpad <- function(x, k) {
  c(rep(NA, k), x)[1 : length(x)] 
}

chicagoNMMAPS$o3_iqr = chicagoNMMAPS$o3/IQR(chicagoNMMAPS$o3, na.rm=T)

# Adjust for temperature, separately for hot and cold temperatures, as in Stafoggia et al. EHP (2013)
# ´Create the necessary variables

chicagoNMMAPS$temp_lag1 = lagpad(chicagoNMMAPS$temp, 1)
chicagoNMMAPS$temp_lag2 = lagpad(chicagoNMMAPS$temp, 2)
chicagoNMMAPS$temp_lag3 = lagpad(chicagoNMMAPS$temp, 3)
chicagoNMMAPS$temp_lag4 = lagpad(chicagoNMMAPS$temp, 4)
chicagoNMMAPS$temp_lag5 = lagpad(chicagoNMMAPS$temp, 5)
chicagoNMMAPS$temp_lag6 = lagpad(chicagoNMMAPS$temp, 6)

chicagoNMMAPS$temp_lag01 = (chicagoNMMAPS$temp + chicagoNMMAPS$temp_lag1 )/2
chicagoNMMAPS$temp_lag16 = (chicagoNMMAPS$temp_lag1 + chicagoNMMAPS$temp_lag2 + chicagoNMMAPS$temp_lag3 + chicagoNMMAPS$temp_lag4 + chicagoNMMAPS$temp_lag5 + chicagoNMMAPS$temp_lag6)/6

chicagoNMMAPS$temp_lag01m = chicagoNMMAPS$temp_lag01
# replace values below the median of temp
chicagoNMMAPS$temp_lag01m[chicagoNMMAPS$temp_lag01 < quantile(chicagoNMMAPS$temp, .5, na.rm=T)] = quantile(chicagoNMMAPS$temp_lag01, .5, na.rm=T)

chicagoNMMAPS$temp_lag16m = chicagoNMMAPS$temp_lag16
# replace values above the median
chicagoNMMAPS$temp_lag16m[chicagoNMMAPS$temp_lag16 > quantile(chicagoNMMAPS$temp, .5, na.rm=T)] = quantile(chicagoNMMAPS$temp_lag16, .5, na.rm=T)


# Degrees of Freedom for Seasonality and Long-Term Trend
df_seas <- 8

# Model Formula (for daily data)

form_mod_o3 <- death ~ o3_iqr + dow + 
  ns(temp_lag01m, knots=c(quantile(chicagoNMMAPS$temp_lag01, c(.75,.9), na.rm=T))) +
  ns(temp_lag16m, knots=c(quantile(chicagoNMMAPS$temp_lag16, c(.25), na.rm=T))) +
  ns(time, df = round( df_seas * nrow(chicagoNMMAPS) / 365.25 ))


#######################
#######################
# Daily response data #
#######################
#######################

# Fit daily model

mod_o3_daily = glm(form_mod_o3, data = chicagoNMMAPS, family = "quasipoisson",
                     na.action = "na.exclude")


# estimated regression coefficient for "o3_iqr"
coef(mod_o3_daily)["o3_iqr"]

# standard error of "o3_iqr" coefficient
sqrt(vcov(mod_o3_daily)["o3_iqr", "o3_iqr"])


########################
########################
# Weekly response data #
########################
########################

# create week variable
#  - weeks start on Mondays

chicagoNMMAPS$week <- c( rep(1,4) , rep(2:((nrow(chicagoNMMAPS) - 4)/7 + 1), each = 7) )

# Need to provide Y, a 2-column matrix with the period indicator and the counts
# Aggregate mortality counts by week - this is the outcome variable used in the analysis 
Y <- aggregate(chicagoNMMAPS$death, by = list(week = chicagoNMMAPS$week), FUN = sum)
names(Y)[2] <- "death"

# Variable that distributes the sum of cases per week equally across all days of the week
#   This variable is used when calculating attributable numbers
Ymean_week <- ave(chicagoNMMAPS$death, week = chicagoNMMAPS$week, FUN = mean)

# Fit model 
mod_o3_week <- fit_aggregate_cpp(Y = Y, X = chicagoNMMAPS, name_exposure = "o3_iqr", 
                                   formula = form_mod_o3, family = "quasipoisson",
                                   tol = 1e-8, maxit = 500, seed = 1234, ntryini = 10, ntryfit = 0) 

mod_o3_week$converged

# estimated regression coefficient for "o3_iqr"
coef(mod_o3_week)["o3_iqr"]

mod_o3_week$error_code
# error_code indicates that variances are infinity. Check:
# standard error of "o3_iqr" coefficient
sqrt( vcov(mod_o3_week)["o3_iqr", "o3_iqr"] )

# we see the same with summary

summary(mod_o3_week)

# When error_code = -555, I can obtain the generalized covariance matrix, which provides appropriate values
sqrt( gvcov_modagr(mod_o3_week)["o3_iqr", "o3_iqr"] )

# we can also use summary_gcov to obtain the summary based on generalized covariance

summary_gcov(mod_o3_week)
# Warnings for 'NaNs produced' because the day of the week effects are not identifiable when using weekly outcome data


####################################################
####################################################
# Monthly (irregular number of days) response data #
####################################################
####################################################

# Create indicator of year and month
chicagoNMMAPS$ymonth <- chicagoNMMAPS$year*100 + chicagoNMMAPS$month

# Need to provide Y, a 2-column matrix with the period indicator and the counts
# Aggregate mortality counts by month - this is the outcome variable used in the analysis 
Y <- aggregate(chicagoNMMAPS$death, by = list(ymonth = chicagoNMMAPS$ymonth), FUN = sum)
names(Y)[2] <- "death"

# Variable that distributes the sum of cases per week equally across all days of the month
#   This variable is used when calculating attributable numbers
Ymean_month <- ave(chicagoNMMAPS$death, ymonth = chicagoNMMAPS$ymonth, FUN = mean)

# Fit model 
mod_o3_month <- fit_aggregate_cpp(Y = Y, X = chicagoNMMAPS, name_exposure = "o3_iqr", 
                                   formula = form_mod_o3, family = "quasipoisson",
                                   tol = 1e-8, maxit = 500, seed = 1234, ntryini = 10, ntryfit = 0) 

mod_o3_month$converged

# estimated regression coefficient for "o3_iqr"
coef(mod_o3_month)["o3_iqr"]

# standard error of "o3_iqr" coefficient
sqrt( vcov(mod_o3_month)["o3_iqr", "o3_iqr"] )


#################################################
#################################################
# Year, Month and day of the week response data #
#################################################
#################################################

# Create indicator of year and month
chicagoNMMAPS$ymonthdow <- chicagoNMMAPS$year*1000 + chicagoNMMAPS$month*10 + as.numeric(chicagoNMMAPS$dow)

# Need to provide Y, a 2-column matrix with the period indicator and the counts
# Aggregate mortality counts by day of the week, month and year - this is the outcome variable used in the analysis 
Y <- aggregate(chicagoNMMAPS$death, by = list(ymonthdow = chicagoNMMAPS$ymonthdow), FUN = sum)
names(Y)[2] <- "death"

# Variable that distributes the sum of cases equally across all days with same day of the week, month and year
#   This variable is used when calculating attributable numbers
Ymean_dow <- ave(chicagoNMMAPS$death, ymonthdow = chicagoNMMAPS$ymonthdow, FUN = mean)

# Fit model 
mod_o3_dow <- fit_aggregate_cpp(Y = Y, X = chicagoNMMAPS, name_exposure = "o3_iqr", 
                                  formula = form_mod_o3, family = "quasipoisson",
                                  tol = 1e-8, maxit = 500, seed = 1234, ntryini = 10, ntryfit = 0) 

mod_o3_dow$converged

# estimated regression coefficient for "o3_iqr"
coef(mod_o3_dow)["o3_iqr"]

# standard error of "o3_iqr" coefficient
sqrt( vcov(mod_o3_dow)["o3_iqr", "o3_iqr"] )



#########################
# Comparison of results #
#########################

coef(mod_o3_daily)["o3_iqr"]
coef(mod_o3_week)["o3_iqr"]
coef(mod_o3_month)["o3_iqr"]
coef(mod_o3_dow)["o3_iqr"]

sqrt( vcov(mod_o3_daily)["o3_iqr", "o3_iqr"] )
sqrt( gvcov_modagr(mod_o3_week)["o3_iqr", "o3_iqr"] )
sqrt( vcov(mod_o3_month)["o3_iqr", "o3_iqr"] )
sqrt( vcov(mod_o3_dow)["o3_iqr", "o3_iqr"] )

