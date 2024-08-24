################################################################################
# R CODE INCLUDED AS SUPPLEMENTARY ONLINE MATERIAL FOR:
#   "Unbiased temperature-related mortality estimates using weekly and monthly health data..."
#   Basagaña X and Ballester J
################################################################################

# Load functions
source("aggregated_functions.r")

# Reading the Input Data (time series from Paris)
data <- read.csv(file="Paris20.txt", sep=",")

# create Date variable
data$date <- as.Date(paste0(data$year,"-",data$month,"-",data$day), "%Y-%m-%d")

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
pred_prc <- c(seq(0,1,0.1), 2:98, seq(99,100,0.1))/100
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
pred_daily <- predict(model_daily, type="response")

# centering point. For now, at the median temperature = 17ºC
cen <- 17

# Predictions of crossbasis associations 
crosspred_daily <- crosspred( cross_basis, model_daily, cen=cen, at = at_temp, bylag = 1)

# Minimum mortality temperature (MMT) based on the model
mmt_daily <- crosspred_daily$predvar[which.min(crosspred_daily$allRRfit)[1]]

# Predictions of crossbasis associations, centering at the MMT 
crosspred_daily <- crosspred( cross_basis, model_daily, cen=mmt_daily, at = at_temp, bylag = 1)

# Reduced (cumulative) crossbasis predictions 
red_daily <- crossreduce(cross_basis, model_daily, cen=mmt_daily, at=at_temp)

# Attributable number

# Total

an_daily <- attrdl(x=data$temp, basis=cross_basis, cases=data$death, model=model_daily,
                         type="an", dir="forw", tot=TRUE, cen=mmt_daily, sim=FALSE)

# Cold

an_daily_cold <- attrdl(x=data$temp, basis=cross_basis, cases=data$death, model=model_daily,
                   type="an", dir="forw", tot=TRUE, cen=mmt_daily, sim=FALSE, range=c(min(data$temp, na.rm=T), mmt_daily))

# Heat

an_daily_heat <- attrdl(x=data$temp, basis=cross_basis, cases=data$death, model=model_daily,
                        type="an", dir="forw", tot=TRUE, cen=mmt_daily, sim=FALSE, range=c(mmt_daily, max(data$temp, na.rm=T)))


########################
########################
# Weekly response data #
########################
########################

# Need to provide Y, a 2-column matrix with the period indicator and the counts
# Aggregate mortality counts by week - this is the outcome variable used in the analysis 
Y <- aggregate(data$death, by=list(week=data$week), FUN=sum)
names(Y)[2] <- "death"

# Variable that distributes the sum of cases per week equally across all 7 days
#   This variable is used when calculating attributable numbers
Ymean_week <- ave(data$death, week=data$week, FUN=mean)

# Fit model 
model_week <- fit_aggregate(Y=Y, X=data, CB=cross_basis, formula=form_mod, family="quasipoisson",
                             tol = 1e-8, maxit = 500, start=NULL) 

# check convergence
model_week$converged

summary(model_week)

# Predict response variable based on the model
pred_week <- predict(model_week)

# centering point. For now, at the median temperature = 17ºC
cen <- 17

# Predictions of crossbasis associations
crosspred_week <- crosspred(cross_basis, model_week, cen=cen, at=at_temp, model.link = "log")

# calculate minimum mortality temperature (MMT) from model predictions 
mmt_week <- crosspred_week$predvar[which.min(crosspred_week$allRRfit)[1]]

# REDUCTION TO OVERALL CUMULATIVE CURVE (ALL LAGS), centering at MMT
red_week <- crossreduce(basis=cross_basis, model=model_week, cen=mmt_week, at=at_temp, model.link = "log")
# adding exponentiated coefficients (RR) to the object with an internal function
red_week <- add_exp(red_week)


# Attributable numbers

# Total
an_week <- attrdl.aggr(x=data$temp, basis=cross_basis, cases=Ymean_week,
                      model=model_week, type="an", dir="forw", tot=TRUE,
                      cen=mmt_week, sim=FALSE)

# Cold

an_week_cold <- attrdl.aggr(x=data$temp, basis=cross_basis, cases=Ymean_week,
                            model=model_week, type="an", dir="forw", tot=TRUE,
                            cen=mmt_week, sim=FALSE, range=c(min(data$temp, na.rm=T), mmt_week))

# Heat

an_week_heat <- attrdl.aggr(x=data$temp, basis=cross_basis, cases=Ymean_week,
                            model=model_week, type="an", dir="forw", tot=TRUE,
                            cen=mmt_week, sim=FALSE, range=c(mmt_week, max(data$temp, na.rm=T)))



####################################################
####################################################
# Monthly (irregular number of days) response data #
####################################################
####################################################

# Create indicator of year and month
data$ymonth <- data$year*100 + data$month

# Aggregate data by year and month - this is the outcome variable used in the analysis
Y <- aggregate(data$death, by=list(ymonth=data$ymonth),FUN=sum)
names(Y)[2] <- "death"

# Variable that distributes the sum of cases per month equally across all days
#   This variable is used when calculating attributable numbers
Ymean_month <- ave(data$death, data$ymonth, FUN=mean)


# Fit model 
model_month <- fit_aggregate(Y=Y, X=data, CB=cross_basis, formula=form_mod, family="quasipoisson",
                              tol = 1e-8, maxit = 5000, start=NULL) 

# check convergence
model_month$converged

summary(model_month)

# centering point. For now, at the median temperature = 17ºC
cen <- 17

# Predictions of crossbasis associations
crosspred_month <- crosspred(cross_basis, model_month, cen=cen, at=at_temp, model.link = "log")

# Calculate minimum mortality temperature (MMT) from model predictions 

# calculate minimum mortality temperature (MMT) from model predictions 
mmt_month <- crosspred_month$predvar[which.min(crosspred_month$allRRfit)[1]]

# REDUCTION TO OVERALL CUMULATIVE CURVE (ALL LAGS), centering at MMT
red_month <- crossreduce(basis=cross_basis, model=model_month, cen=mmt_month, at=at_temp, model.link = "log")
# adding exponentiated coefficients (RR) to the object with an internal function
red_month <- add_exp(red_month)


# Attributable number

# Total

an_month <- attrdl.aggr(x=data$temp, basis=cross_basis, cases=Ymean_month,
                       model=model_month, type="an", dir="forw", tot=TRUE, 
                       cen=mmt_month,sim=FALSE)

an_month_cold <- attrdl.aggr(x=data$temp, basis=cross_basis, cases=Ymean_month,
                        model=model_month, type="an", dir="forw", tot=TRUE, 
                        cen=mmt_month,sim=FALSE, range=c(min(data$temp, na.rm=T), mmt_month))

an_month_heat <- attrdl.aggr(x=data$temp, basis=cross_basis, cases=Ymean_month,
                             model=model_month, type="an", dir="forw", tot=TRUE, 
                             cen=mmt_month,sim=FALSE, range=c(mmt_month, max(data$temp, na.rm=T)))

##########################
# Comparision of results #
##########################

# Plot cumulative exposure-response functions

plot(red_daily, xlab="Temperature",ylab="RR", ci="n", lwd=2)
lines(red_week, ci="n", col="red", lwd=2)
lines(red_month, ci="n", col="blue", lwd=2)
legend("top",c("Daily", "Weekly", "Monthly"), lty=1, lwd=2, col=c("black", "red","blue"))


# Compare ANs

an_daily
an_week
an_month

an_daily_cold
an_week_cold
an_month_cold

an_daily_heat
an_week_heat
an_month_heat
