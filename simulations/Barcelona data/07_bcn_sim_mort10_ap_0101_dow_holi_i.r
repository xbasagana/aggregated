
source("agr_functions_v36_cpp.r")

library(splines)
library(Rfast)
library(MASS)


## Function to fit the air pollution models
## For memory issues, I call the function separately for each number of years

run_ap = function(nsim, i.years, dati, outc, rhs_simulate, form_simulate, rhs_fit, incl.dow, incl.holi, coef_true, disp_par, ind_expos) {  
  
  # Obtain matrix of predictors - simulate an arbitrary outcome variable
  dati[, outc] <- rpois(nrow(dati), 10)
  X <- model.matrix.lm(as.formula(form_simulate), data=dati, na.action = "na.pass")
  
  # Calculate linear predictor
  linpred <- X %*% coef_true
  
  # Formula that will be used to fit model using simulated data
  form_mod_fit <- paste0( "Y_sim", rhs_fit ) 
  
  # Initialiaze elements
  
  Ysim_mat <- matrix(NA, nrow=length(linpred), ncol=nsim)
  Ypred_daily <- matrix(NA, nrow=length(linpred), ncol=nsim)
  Ypred_week <- matrix(NA, nrow=length(linpred), ncol=nsim)
  Ypred_month <- matrix(NA, nrow=length(linpred), ncol=nsim)
  Ypred_dow <- matrix(NA, nrow=length(linpred), ncol=nsim)
  
  if (incl.dow & incl.holi) {
    coefs_daily <- coefs_week <- coefs_month <- coefs_dow <- matrix(NA, nrow=nsim, ncol=nrow(coef_true))
    ses_daily <- ses_week <- ses_month <- ses_dow <- matrix(NA, nrow=nsim, ncol=nrow(coef_true)) 
  } else if (!incl.dow & incl.holi) {
    coefs_daily <- coefs_week <- coefs_month <- coefs_dow <- matrix(NA, nrow=nsim, ncol=nrow(coef_true) - 6) # excluding dow
    ses_daily <- ses_week <- ses_month <- ses_dow <- matrix(NA, nrow=nsim, ncol=nrow(coef_true) - 6) # excluding dow
  } else if (incl.dow & !incl.holi) {
    coefs_daily <- coefs_week <- coefs_month <- coefs_dow <- matrix(NA, nrow=nsim, ncol=nrow(coef_true) - 1) # excluding holi
    ses_daily <- ses_week <- ses_month <- ses_dow <- matrix(NA, nrow=nsim, ncol=nrow(coef_true) - 1) # excluding holi
  } else if (!incl.dow & !incl.holi) {
    coefs_daily <- coefs_week <- coefs_month <- coefs_dow <- matrix(NA, nrow=nsim, ncol=nrow(coef_true) - 7) # excluding dow, holi
    ses_daily <- ses_week <- ses_month <- ses_dow <- matrix(NA, nrow=nsim, ncol=nrow(coef_true) - 7) # excluding dow, holi
  }

  disp_daily <- disp_week <- disp_month <- disp_dow <- vector()
  
  AF_daily <- AF_week <- AF_month <- AF_dow <- vector()
  lower_AF_daily <- lower_AF_week <- lower_AF_month <- lower_AF_dow <- vector()
  upper_AF_daily <- upper_AF_week <- upper_AF_month <- upper_AF_dow <- vector()
  
  converged_week <- converged_month <- converged_dow <- vector()
  
  # Loop over number of simulations
  for (i in 1:nsim) {
    if (i==1) cat("i= ")
    if ((i%%10)==0) cat(i,",\n") else cat(i,", ")
    
    # Simulate Y
    Y_sim <- rpois.od(n=nrow(dati), lambda=exp(linpred)/10, d=disp_par)
    Y_sim[is.nan(Y_sim)] <- NA
    dati$Y_sim <- Y_sim
    
    Ysim_mat[,i] <- Y_sim
    
    ################################################
    ### Fit daily model
    ################################################
    
    mod_daily <- glm(as.formula(form_mod_fit), data=dati, family = "quasipoisson", na.action = "na.exclude" )
    Ypred_daily[,i] <- predict(mod_daily, type="response")
    
    coefs_daily[i,] <- coef(mod_daily)
    ses_daily[i,] <- sqrt( diag(vcov(mod_daily)) )
    
    disp_daily[i] <- summary(mod_daily)$dispersion
    
    
    # Attributable number
    
    RRj = exp(coef(mod_daily)[ind_expos]*X[, ind_expos])
    
    AF_daily[i] <- sum( (Y_sim * (RRj-1)/RRj) / sum(Y_sim, na.rm=T), na.rm=T)
    coefsim = coef(mod_daily)[ind_expos] + rnorm(nrow(X), 0, sd=sqrt(vcov(mod_daily)[ind_expos, ind_expos]))
    RRjsim <- exp( t(t(matrix(rep( X[, ind_expos], nrow(X)), nrow = nrow(X), byrow=F )) * coefsim ))
    AFsim = colsums( Y_sim * ((RRjsim-1)/RRjsim)  / sum(Y_sim, na.rm=T) , na.rm=T )
    lower_AF_daily[i] = quantile(AFsim, 0.025)
    upper_AF_daily[i] = quantile(AFsim, 0.975)
    
    rm(mod_daily)
    
    
    ################################################
    ### Fit aggregated model (Y: weekly, X: daily)
    ################################################
    
    # Weekly counts
    
    Y <- aggregate(Y_sim, by=list(week=dati$week),FUN=sum)
    names(Y)[2] <- "Y_sim"
    
    avgcount <- ave(Y_sim, dati$week)
    
    # Fit model 
    
    try( mod_agr <- fit_aggregate_cpp(Y=Y, X=dati, name_exposure = "no2_lag01_20", formula=as.formula(form_mod_fit), family="quasipoisson",
                                      tol = 1e-8, maxit = 500, seed=NULL, ntryini = 10, ntryfit=0))
    if (class(mod_agr)[1] == "try-error" | (class(mod_agr)=="modagr" & !(mod_agr$error_code == 0 | mod_agr$error_code == -555))  )  {
      converged_week[i] <- FALSE
    } else {
      if (mod_agr$converged==TRUE) {
        converged_week[i] <- mod_agr$converged
        
        coefs_week[i,] <- coef(mod_agr)
        ses_week[i,] <- sqrt( diag(vcov(mod_agr)) )
        
        disp_week[i] <- mod_agr$disp
        
        Ypred_week[,i] <- predict(mod_agr)
        
        if (!(any(!(is.finite(vcov(mod_agr))) ) )) {
          
          # Attributable number
          
          RRj = exp(coef(mod_agr)[ind_expos]*X[, ind_expos])
          
          AF_week[i] <- sum( (avgcount * (RRj-1)/RRj) / sum(avgcount, na.rm=T), na.rm=T)
          coefsim = coef(mod_agr)[ind_expos] + rnorm(nrow(X), 0, sd=sqrt(vcov(mod_agr)[ind_expos, ind_expos]))
          RRjsim <- exp( t(t(matrix(rep( X[, ind_expos], nrow(X)), nrow = nrow(X), byrow=F )) * coefsim ))
          AFsim = colsums( Y_sim * ((RRjsim-1)/RRjsim)  / sum(Y_sim, na.rm=T) , na.rm=T )
          lower_AF_week[i] = quantile(AFsim, 0.025)
          upper_AF_week[i] = quantile(AFsim, 0.975)
          
        } else {
          converged_week[i] <- FALSE
        }
      } else if (mod_agr$converged == FALSE & mod_agr$error_code == -555) {
        converged_week[i] <- mod_agr$converged
        
        coefs_week[i,] <- coef(mod_agr)
        vcov1 = vcov(mod_agr)
        D2 = d2func_aggr_cpp(beta = coef(mod_agr), X = mod_agr$X, y = mod_agr$y, period = mod_agr$period)
        vcov2 = mod_agr$disp * ginv(-D2)
        ses_week[i,] <- sqrt( diag(vcov2) )
        
        disp_week[i] <- mod_agr$disp
        
        Ypred_week[,i] <- predict(mod_agr)
        
        if (!(any(!(is.finite(vcov2)) ) )) {
          
          # Attributable number
          
          RRj = exp(coef(mod_agr)[ind_expos]*X[, ind_expos])
          
          AF_week[i] <- sum( (avgcount * (RRj-1)/RRj) / sum(avgcount, na.rm=T), na.rm=T)
          coefsim = coef(mod_agr)[ind_expos] + rnorm(nrow(X), 0, sd=sqrt(vcov2[ind_expos, ind_expos]))
          RRjsim <- exp( t(t(matrix(rep( X[, ind_expos], nrow(X)), nrow = nrow(X), byrow=F )) * coefsim ))
          AFsim = colsums( Y_sim * ((RRjsim-1)/RRjsim)  / sum(Y_sim, na.rm=T) , na.rm=T )
          lower_AF_week[i] = quantile(AFsim, 0.025)
          upper_AF_week[i] = quantile(AFsim, 0.975)
          
        }
        
      } else {
        converged_week[i] <- FALSE
      }
      
    }
    rm(mod_agr)
    
    ################################################
    ### Fit aggregated model (Y: monthly, X: daily)
    ################################################
    
    # Monthly counts
    
    Y <- aggregate(Y_sim, by=list(ymonth=dati$ymonth), FUN=sum)
    names(Y)[2] <- "Y_sim"
    
    avgcount <- ave(Y_sim, dati$ymonth)
    
    
    # Fit model 
    
    try( mod_agr <- fit_aggregate_cpp(Y=Y, X=dati, name_exposure = "no2_lag01_20", formula=as.formula(form_mod_fit), family="quasipoisson",
                                      tol = 1e-8, maxit = 500, seed=NULL, ntryini = 10, ntryfit=0))
    if (class(mod_agr)[1] == "try-error" | (class(mod_agr)=="modagr" & !(mod_agr$error_code == 0 | mod_agr$error_code == -555))  )  {
        converged_month[i] <- FALSE
    } else {  
      if (mod_agr$converged==TRUE) {
        
        converged_month[i] <- mod_agr$converged
        
        coefs_month[i,] <- coef(mod_agr)
        ses_month[i,] <- sqrt( diag(vcov(mod_agr)) )
        
        disp_month[i] <- mod_agr$disp
        
        Ypred_month[,i] <- predict(mod_agr)
        
        if (!(any(!(is.finite(vcov(mod_agr))) ) )) {
          
          # Attributable number
          
          RRj = exp(coef(mod_agr)[ind_expos]*X[, ind_expos])
          AF_month[i] <- sum( (avgcount * (RRj-1)/RRj) / sum(avgcount, na.rm=T), na.rm=T)
          coefsim = coef(mod_agr)[ind_expos] + rnorm(nrow(X), 0, sd=sqrt(vcov(mod_agr)[ind_expos, ind_expos]))
          RRjsim <- exp( t(t(matrix(rep( X[, ind_expos], nrow(X)), nrow = nrow(X), byrow=F )) * coefsim ))
          AFsim = colsums( Y_sim * ((RRjsim-1)/RRjsim)  / sum(Y_sim, na.rm=T) , na.rm=T )
          lower_AF_month[i] = quantile(AFsim, 0.025)
          upper_AF_month[i] = quantile(AFsim, 0.975)
          
        } else {
          converged_month[i] <- FALSE
        }
      } else if (mod_agr$converged == FALSE & mod_agr$error_code == -555) {
        converged_month[i] <- mod_agr$converged
        
        coefs_month[i,] <- coef(mod_agr)
        vcov1 = vcov(mod_agr)
        D2 = d2func_aggr_cpp(beta = coef(mod_agr), X = mod_agr$X, y = mod_agr$y, period = mod_agr$period)
        vcov2 = mod_agr$disp * ginv(-D2)
        ses_month[i,] <- sqrt( diag(vcov2) )
        
        disp_month[i] <- mod_agr$disp
        
        Ypred_month[,i] <- predict(mod_agr)
        
        if (!(any(!(is.finite(vcov2)) ) )) {
          
          # Attributable number
          
          RRj = exp(coef(mod_agr)[ind_expos]*X[, ind_expos])
          
          AF_month[i] <- sum( (avgcount * (RRj-1)/RRj) / sum(avgcount, na.rm=T), na.rm=T)
          coefsim = coef(mod_agr)[ind_expos] + rnorm(nrow(X), 0, sd=sqrt(vcov2[ind_expos, ind_expos]))
          RRjsim <- exp( t(t(matrix(rep( X[, ind_expos], nrow(X)), nrow = nrow(X), byrow=F )) * coefsim ))
          AFsim = colsums( Y_sim * ((RRjsim-1)/RRjsim)  / sum(Y_sim, na.rm=T) , na.rm=T )
          lower_AF_month[i] = quantile(AFsim, 0.025)
          upper_AF_month[i] = quantile(AFsim, 0.975)
          
        }
        
      } else {
        converged_month[i] <- FALSE
      }
      
    }
    rm(mod_agr)
    
    ################################################
    ### Fit aggregated model (Y: dow, X: daily)
    ################################################
    
    # dow counts
    
    Y <- aggregate(Y_sim, by=list(ymonthdow=dati$ymonthdow), FUN=sum)
    names(Y)[2] <- "Y_sim"
    
    avgcount <- ave(Y_sim, dati$ymonthdow)
    
    # Fit model 
    
    try( mod_agr <- fit_aggregate_cpp(Y=Y, X=dati, name_exposure = "no2_lag01_20", formula=as.formula(form_mod_fit), family="quasipoisson",
                                      tol = 1e-8, maxit = 500, seed=NULL, ntryini = 10, ntryfit=0))
    
    if (class(mod_agr)[1] == "try-error" | (class(mod_agr)=="modagr" & mod_agr$error_code!=0))  {
      converged_dow[i.years,i] <- FALSE
    } else {
      if (mod_agr$converged==TRUE) {
        converged_dow[i] <- mod_agr$converged
        
        coefs_dow[i,] <- coef(mod_agr)
        ses_dow[i,] <- sqrt( diag(vcov(mod_agr)) )
        
        disp_dow[i] <- mod_agr$disp
        
        Ypred_dow[,i] <- predict(mod_agr)
        
        if (!(any(!(is.finite(vcov(mod_agr))) ) )) {
          
          # Attributable number
          
          RRj = exp(coef(mod_agr)[ind_expos]*X[, ind_expos])
          
          AF_dow[i] <- sum( (avgcount * (RRj-1)/RRj) / sum(avgcount, na.rm=T), na.rm=T)
          coefsim = coef(mod_agr)[ind_expos] + rnorm(nrow(X), 0, sd=sqrt(vcov(mod_agr)[ind_expos, ind_expos]))
          RRjsim <- exp( t(t(matrix(rep( X[, ind_expos], nrow(X)), nrow = nrow(X), byrow=F )) * coefsim ))
          AFsim = colsums( Y_sim * ((RRjsim-1)/RRjsim)  / sum(Y_sim, na.rm=T) , na.rm=T )
          lower_AF_dow[i] = quantile(AFsim, 0.025)
          upper_AF_dow[i] = quantile(AFsim, 0.975)
          
        } else {
          converged_dow[i] <- FALSE
        }
      } else {
        converged_dow[i] <- FALSE
      }  
    }
    rm(mod_agr)
    gc()
  }
  
  return(list(nsim, Ysim_mat, 
              converged_week, converged_month, converged_dow,
              coefs_daily, coefs_week, coefs_month, coefs_dow,
              ses_daily, ses_week, ses_month, ses_dow,       
              Ypred_daily, Ypred_week, Ypred_month, Ypred_dow,  
              AF_daily, AF_week, AF_month, AF_dow, 
              lower_AF_daily, lower_AF_week, lower_AF_month, lower_AF_dow, 
              upper_AF_daily, upper_AF_week, upper_AF_month, upper_AF_dow, 
              disp_daily, disp_week, disp_month, disp_dow))
  
}


# Function to generate random numbers following
# over-dispersed Poisson distribution
# based on simple cheat of using a standard negative binomial,
# but choosing the scale parameter to give the desired mean/variance
# ratio at the given value of the mean.
# Taken from: https://stat.ethz.ch/pipermail/r-help/2002-June/022425.html

rpois.od<-function (n, lambda,d=1) {
  if (d==1) {
    rpois(n, lambda)
  } else {
    suppressWarnings(rnbinom(n, size=(lambda/(d-1)), mu=lambda))
  }
}

##################
### Parameters
##################

outc = "mort"
anal = "mort_ap"
nsim <- 500

incl.dow = TRUE
incl.holi = TRUE

load(file = "bcn_ts_exposures.RData")

##################################################################
# Specify parameters of the model: outcome- & exposure-specific 
##################################################################

# Read parameters to simulate from

load(file=paste0("parameters_", anal, "_01_dow_holi.RData"))

# function to create lags
lagpad <- function(x, k) {
  c(rep(NA, k), x)[1 : length(x)] 
}

data$no2_lag1 = lagpad(data$no2interp, 1)
data$no2_lag2 = lagpad(data$no2interp, 2)
data$no2_lag3 = lagpad(data$no2interp, 3)
data$no2_lag4 = lagpad(data$no2interp, 4)
data$no2_lag5 = lagpad(data$no2interp, 5)
data$no2_lag6 = lagpad(data$no2interp, 6)

data$no2_lag01 = (data$no2interp + data$no2_lag1)/2

data$no2_lag06 = (data$no2interp + data$no2_lag1 + data$no2_lag2 + data$no2_lag3 +
                    data$no2_lag4 + data$no2_lag5 + data$no2_lag6)/7

# For 20 units (~IQR of no2_lag01)
data$no2_lag01_20 = data$no2_lag01/20
data$no2_lag06_20 = data$no2_lag06/20

data$tmean_lag1 = lagpad(data$tmean, 1)
data$tmean_lag2 = lagpad(data$tmean, 2)
data$tmean_lag3 = lagpad(data$tmean, 3)
data$tmean_lag4 = lagpad(data$tmean, 4)
data$tmean_lag5 = lagpad(data$tmean, 5)
data$tmean_lag6 = lagpad(data$tmean, 6)

data$tmean_lag01 = (data$tmean + data$tmean_lag1)/2
data$tmean_lag16 = (data$tmean_lag1 + data$tmean_lag2 + data$tmean_lag3 + data$tmean_lag4 + data$tmean_lag5 + data$tmean_lag6)/6


data$tmean_lag01m = data$tmean_lag01
# replace values below the median
data$tmean_lag01m[data$tmean_lag01 < quantile(data$tmean_lag01, .5, na.rm=T)] = quantile(data$tmean_lag01, .5, na.rm=T)

data$tmean_lag16m = data$tmean_lag16
# replace values above the median
data$tmean_lag16m[data$tmean_lag16 > quantile(data$tmean_lag16, .5, na.rm=T)] = quantile(data$tmean_lag16, .5, na.rm=T)



# Degrees of Freedom for Seasonality and Long-Term Trend
df_seas <- 8


# Model Formula to simulate the data
rhs_simulate = " ~ no2_lag01_20 + dow + holidays_nat_reg + ns(tmean_lag01m, knots=c(quantile(data$tmean_lag01, c(.75,.9), na.rm=T))) +
  ns(tmean_lag16m, knots=c(quantile(data$tmean_lag16, c(.25), na.rm=T))) + ns( Date, df = round( df_seas * nrow(dati) / 365.25 ))"

if (incl.dow == FALSE & incl.holi == FALSE) {
  rhs_fit = " ~ no2_lag01_20 + ns(tmean_lag01m, knots=c(quantile(data$tmean_lag01, c(.75,.9), na.rm=T))) +
  ns(tmean_lag16m, knots=c(quantile(data$tmean_lag16, c(.25), na.rm=T))) + ns( Date, df = round( df_seas * nrow(dati) / 365.25 ))"
  dow_holi_name = "nodow_noholi"
} else if (incl.dow == TRUE & incl.holi == FALSE) {
  rhs_fit = " ~ no2_lag01_20 + dow + ns(tmean_lag01m, knots=c(quantile(data$tmean_lag01, c(.75,.9), na.rm=T))) +
  ns(tmean_lag16m, knots=c(quantile(data$tmean_lag16, c(.25), na.rm=T))) + ns( Date, df = round( df_seas * nrow(dati) / 365.25 ))"
  dow_holi_name = "dow_noholi"
} else if (incl.dow == FALSE & incl.holi == TRUE) {
  rhs_fit = " ~ no2_lag01_20 + holidays_nat_reg + ns(tmean_lag01m, knots=c(quantile(data$tmean_lag01, c(.75,.9), na.rm=T))) +
  ns(tmean_lag16m, knots=c(quantile(data$tmean_lag16, c(.25), na.rm=T))) + ns( Date, df = round( df_seas * nrow(dati) / 365.25 ))"
  dow_holi_name = "nodow_holi"
} else if (incl.dow == TRUE & incl.holi == TRUE) {
  rhs_fit = " ~ no2_lag01_20 + dow +  holidays_nat_reg + ns(tmean_lag01m, knots=c(quantile(data$tmean_lag01, c(.75,.9), na.rm=T))) +
  ns(tmean_lag16m, knots=c(quantile(data$tmean_lag16, c(.25), na.rm=T))) + ns( Date, df = round( df_seas * nrow(dati) / 365.25 ))"
  dow_holi_name = "dow_holi"
}

form_simulate <- paste0( outc, rhs_simulate)

###################################################
## Datasets by year, and rest of the coefficients
###################################################

inddf <- data.frame(ind1 = (data$year>=2019), ind2 = (data$year>=2018),
                    ind3 = (data$year>=2017), ind4 = (data$year>=2016),
                    ind5 = (data$year>=2015), ind6 = (data$year>=2014),
                    ind7 = (data$year>=2013), ind8 = (data$year>=2012),
                    ind9 = (data$year>=2011), ind10 = (data$year>=2010),
                    ind11 = (data$year>=2009), ind12 = (data$year>=2008),
                    ind13 = (data$year>=2007), ind14 = (data$year>=2006),
                    ind15 = (data$year>=2005), ind16 = (data$year>=2004),
                    ind17 = (data$year>=2003), ind18 = (data$year>=2002),
                    ind19 = (data$year>=2001), ind20 = (data$year>=2000),
                    ind21 = (data$year>=1999), ind22 = (data$year>=1998),
                    ind23 = (data$year>=1997), ind24 = (data$year>=1996),
                    ind25 = (data$year>=1995), ind26 = (data$year>=1994),
                    ind27 = (data$year>=1993), ind28 = (data$year>=1992),
                    ind29 = (data$year>=1991), ind30 = (data$year>=1990),
                    ind31 = (data$year>=1989), ind32 = (data$year>=1988),
                    ind33 = (data$year>=1987), ind34 = (data$year>=1986),
                    ind35 = (data$year>=1985), ind36 = (data$year>=1984),
                    ind37 = (data$year>=1983), ind38 = (data$year>=1982),
                    ind39 = (data$year>=1981), ind40 = (data$year>=1980))

indyears <- c(1:10, seq(15, 25, 5))


######################
# Start simulation
######################


for (i.years in 1:length(indyears)) {
  cat(",\n", "Years: ", indyears[i.years], ",\n")

  # Select years
  dati <- data[inddf[,indyears[i.years]],]

  # Put together true coefficients 
  coef_true <- matrix(c(coef_int[i.years], coef_expos, coef_rest[[i.years]]), ncol=1)
  
  rlist = run_ap(nsim=nsim, i.years=i.years, dati=dati, outc=outc, rhs_simulate = rhs_simulate, form_simulate=form_simulate, 
                             rhs_fit = rhs_fit, incl.dow = incl.dow, incl.holi = incl.holi, coef_true=coef_true, disp_par=disp_par, ind_expos=2)
  
  save(rlist, file = paste0("results_sim_bcn10_v14_", anal, "_", dow_holi_name, "_0101_", indyears[i.years], ".RData"))
  rm(rlist, dati, coef_true)  
}

