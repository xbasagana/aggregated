
source("agr_functions_v36_cpp.r")

library(splines)
library(Rfast)
library(FluMoDL)

## Function to fit the temperature models
## For memory issues, I call the function separately for each number of years

run_temp = function(nsim, i.years, dati, crossbasis, outc, form_mod_i, rhs_i,
                    coef_true, disp_par,
                    pred_prc, at_temp, pos_p99, pos_p01) {  
  
  # Obtain matrix of predictors - simulate an arbitrary outcome variable
  dati[, outc] <- rpois(nrow(dati), 10)
  X <- model.matrix.lm(as.formula(form_mod_i), data=dati, na.action = "na.pass")
  
  # Calculate linear predictor
  linpred <- X %*% coef_true
  
  # Formula that will be used to fit model using simulated data
  form_mod_sim <- paste0( "Y_sim", rhs_i ) 
  
  # Initialiaze elements 
  
  Ysim_mat <- matrix(NA, nrow=length(linpred), ncol=nsim)
  Ypred_daily <- matrix(NA, nrow=length(linpred), ncol=nsim)
  Ypred_week <- matrix(NA, nrow=length(linpred), ncol=nsim)
  Ypred_month <- matrix(NA, nrow=length(linpred), ncol=nsim)
  Ypred_dow <- matrix(NA, nrow=length(linpred), ncol=nsim)
  
  coefs_daily <-  coefs_week <- coefs_month <- coefs_dow <- matrix(NA, nrow=nsim, ncol=nrow(coef_true))
  
  cumpred_daily <-  cumpred_week <- cumpred_month <- cumpred_dow <- matrix(NA, nrow=nsim, ncol=length(unique(at_temp)))
  cumpred_low_daily <-  cumpred_low_week <- cumpred_low_month <- cumpred_low_dow <- matrix(NA, nrow=nsim, ncol=length(unique(at_temp)))
  cumpred_high_daily <-  cumpred_high_week <- cumpred_high_month <- cumpred_high_dow <- matrix(NA, nrow=nsim, ncol=length(unique(at_temp)))
  
  cumpred_daily_fixedmmt <-  cumpred_week_fixedmmt <- cumpred_month_fixedmmt <- cumpred_dow_fixedmmt <- matrix(NA, nrow=nsim, ncol=length(unique(at_temp)))
  cumpred_low_daily_fixedmmt <-  cumpred_low_week_fixedmmt <- cumpred_low_month_fixedmmt <- cumpred_low_dow_fixedmmt <- matrix(NA, nrow=nsim, ncol=length(unique(at_temp)))
  cumpred_high_daily_fixedmmt <-  cumpred_high_week_fixedmmt <- cumpred_high_month_fixedmmt <- cumpred_high_dow_fixedmmt <- matrix(NA, nrow=nsim, ncol=length(unique(at_temp)))
  
  RRp99lag_daily_fixedmmt <- RRp99lag_week_fixedmmt <- RRp99lag_month_fixedmmt <- RRp99lag_dow_fixedmmt <- matrix(NA, nrow=nsim, ncol=max_lag+1)
  RRp99lag_low_daily_fixedmmt <- RRp99lag_low_week_fixedmmt <- RRp99lag_low_month_fixedmmt <- RRp99lag_low_dow_fixedmmt <- matrix(NA, nrow=nsim, ncol=max_lag+1)
  RRp99lag_high_daily_fixedmmt <- RRp99lag_high_week_fixedmmt <- RRp99lag_high_month_fixedmmt <- RRp99lag_high_dow_fixedmmt <- matrix(NA, nrow=nsim, ncol=max_lag+1)
  
  RRp01lag_daily_fixedmmt <- RRp01lag_week_fixedmmt <- RRp01lag_month_fixedmmt <- RRp01lag_dow_fixedmmt <- matrix(NA, nrow=nsim, ncol=max_lag+1)
  RRp01lag_low_daily_fixedmmt <- RRp01lag_low_week_fixedmmt <- RRp01lag_low_month_fixedmmt <- RRp01lag_low_dow_fixedmmt <- matrix(NA, nrow=nsim, ncol=max_lag+1)
  RRp01lag_high_daily_fixedmmt <- RRp01lag_high_week_fixedmmt <- RRp01lag_high_month_fixedmmt <- RRp01lag_high_dow_fixedmmt <- matrix(NA, nrow=nsim, ncol=max_lag+1)
  
  converged_week <- converged_month <- converged_dow <- vector()
  
  mmt_daily <- mmt_week <-  mmt_month <- mmt_dow <- vector()
  disp_daily <- disp_week <- disp_month <- disp_dow <- vector()
  
  # The same but centering at true MMT
  
  AN_daily_fixedmmt <- AN_week_fixedmmt <- AN_month_fixedmmt <- AN_dow_fixedmmt <- vector()
  AN_daily_cold_fixedmmt <- AN_week_cold_fixedmmt <- AN_month_cold_fixedmmt <- AN_dow_cold_fixedmmt <- vector()
  AN_daily_heat_fixedmmt <- AN_week_heat_fixedmmt <- AN_month_heat_fixedmmt <- AN_dow_heat_fixedmmt <- vector()
  
  lower_AN_daily_fixedmmt <- lower_AN_week_fixedmmt <- lower_AN_month_fixedmmt <- lower_AN_dow_fixedmmt <- vector()
  lower_AN_daily_cold_fixedmmt <- lower_AN_week_cold_fixedmmt <- lower_AN_month_cold_fixedmmt <- lower_AN_dow_cold_fixedmmt <- vector()
  lower_AN_daily_heat_fixedmmt <- lower_AN_week_heat_fixedmmt <- lower_AN_month_heat_fixedmmt <- lower_AN_dow_heat_fixedmmt <- vector()
  
  upper_AN_daily_fixedmmt <- upper_AN_week_fixedmmt <- upper_AN_month_fixedmmt <- upper_AN_dow_fixedmmt <- vector()
  upper_AN_daily_cold_fixedmmt <- upper_AN_week_cold_fixedmmt <- upper_AN_month_cold_fixedmmt <- upper_AN_dow_cold_fixedmmt <- vector()
  upper_AN_daily_heat_fixedmmt <- upper_AN_week_heat_fixedmmt <- upper_AN_month_heat_fixedmmt <- upper_AN_dow_heat_fixedmmt <- vector()
  
  
  # Loop over number of simulations
  for (i in 1:nsim) {
    if (i==1) cat("i= ")
    if ((i%%10)==0) cat(i,",\n") else cat(i,", ")
    
    # Simulate Y
    Y_sim <- rpois.od(n=nrow(dati), lambda=exp(linpred), d=disp_par)
    Y_sim[is.nan(Y_sim)] <- NA
    dati$Y_sim <- Y_sim
    
    Ysim_mat[,i] <- Y_sim
    
    ################################################
    ### Fit daily model
    ################################################
    
    mod_daily <- glm(form_mod_sim, data=dati, family = "quasipoisson", na.action = "na.exclude" )
    Ypred_daily[,i] <- predict(mod_daily, type="response")
    
    coefs_daily[i,] <- coef(mod_daily)
    disp_daily[i] <- summary(mod_daily)$dispersion
    
    crosspred_daily <- crosspred( cross_basis, mod_daily, cen = 20, at = at_temp, bylag = 1, model.link="log")
    
    iMMT <- 1 + which( diff( sign( diff( crosspred_daily$allRRfit ) ) ) == 2 ) 
    if( length( iMMT ) > 0 ) { 
      MMT <- crosspred_daily$predvar[iMMT[ which.min( crosspred_daily$allRRfit[ iMMT ] ) ]  ]
    } else {
      MMT <- crosspred_daily$predvar[ which( pred_prc == .1 ) - 1 + which.min( crosspred_daily$allRRfit[ which( pred_prc == 0.1 ) : which( pred_prc == 0.9 ) ] ) ]
    } 
    
    mmt_daily[i] <- MMT
    crosspred_daily <- crosspred( cross_basis, mod_daily, cen = MMT, at = at_temp, bylag = 1)
    cumpred_daily[i,] <- crosspred_daily$allRRfit
    cumpred_low_daily[i,] <- crosspred_daily$allRRlow
    cumpred_high_daily[i,] <- crosspred_daily$allRRhigh
    
    
    # Attributable number
    
    # Using fixed MMT (MMTtrue)
    
    crosspred <-  crosspred( cross_basis, coef=(coefs_daily[i,])[ind_expos],
                             vcov=(vcov(mod_daily))[ind_expos, ind_expos],
                             cen = MMTtrue, at = at_temp, bylag = 1, model.link="log")
    cumpred_daily_fixedmmt[i,] <- crosspred$allRRfit
    cumpred_low_daily_fixedmmt[i,] <- crosspred$allRRlow
    cumpred_high_daily_fixedmmt[i,] <- crosspred$allRRhigh
    
    RRp99lag_daily_fixedmmt[i,] <- crosspred$matRRfit[pos_p99,]
    RRp99lag_low_daily_fixedmmt[i,] <- crosspred$matRRlow[pos_p99,]
    RRp99lag_high_daily_fixedmmt[i,] <- crosspred$matRRhigh[pos_p99,]
    
    RRp01lag_daily_fixedmmt[i,] <- crosspred$matRRfit[pos_p01,]
    RRp01lag_low_daily_fixedmmt[i,] <- crosspred$matRRlow[pos_p01,]
    RRp01lag_high_daily_fixedmmt[i,] <- crosspred$matRRhigh[pos_p01,]
    
    AN_daily_fixedmmt[i] <- attrdl(x=dati$tmean, basis=cross_basis, cases=dati$Y_sim, 
                                           coef=(coefs_daily[i,])[ind_expos],
                                           vcov=(vcov(mod_daily))[ind_expos, ind_expos],
                                           type="an", dir="forw", tot=TRUE, cen=MMTtrue, sim=FALSE)
    ANsim <- attrdl(x=dati$tmean, basis=cross_basis, cases=dati$Y_sim, 
                    coef=(coefs_daily[i,])[ind_expos],
                    vcov=(vcov(mod_daily))[ind_expos, ind_expos],
                    type="an", dir="forw", tot=TRUE, cen=MMTtrue, sim=TRUE, nsim=1000)
    lower_AN_daily_fixedmmt[i] = quantile(ANsim, 0.025)
    upper_AN_daily_fixedmmt[i] = quantile(ANsim, 0.975)
    rm(ANsim)
    
    AN_daily_cold_fixedmmt[i] <- attrdl(x=dati$tmean, basis=cross_basis, cases=dati$Y_sim,
                                                coef=(coefs_daily[i,])[ind_expos],
                                                vcov=(vcov(mod_daily))[ind_expos, ind_expos],
                                                type="an", dir="forw", tot=TRUE, cen=MMTtrue, sim=FALSE, range=c(-100, MMTtrue))
    ANsim <- attrdl(x=dati$tmean, basis=cross_basis, cases=dati$Y_sim,
                    coef=(coefs_daily[i,])[ind_expos],
                    vcov=(vcov(mod_daily))[ind_expos, ind_expos],
                    type="an", dir="forw", tot=TRUE, cen=MMTtrue, sim=TRUE, nsim=1000, range=c(-100, MMTtrue))
    lower_AN_daily_cold_fixedmmt[i] = quantile(ANsim, 0.025)
    upper_AN_daily_cold_fixedmmt[i] = quantile(ANsim, 0.975)
    rm(ANsim)
    
    AN_daily_heat_fixedmmt[i] <- attrdl(x=dati$tmean, basis=cross_basis, cases=dati$Y_sim, 
                                                coef=(coefs_daily[i,])[ind_expos],
                                                vcov=(vcov(mod_daily))[ind_expos, ind_expos],
                                                type="an", dir="forw", tot=TRUE, cen=MMTtrue, sim=FALSE, range=c(MMTtrue, 100))
    ANsim <- attrdl(x=dati$tmean, basis=cross_basis, cases=dati$Y_sim, 
                    coef=(coefs_daily[i,])[ind_expos],
                    vcov=(vcov(mod_daily))[ind_expos, ind_expos],
                    type="an", dir="forw", tot=TRUE, cen=MMTtrue, sim=TRUE, nsim=1000, range=c(MMTtrue, 100))
    lower_AN_daily_heat_fixedmmt[i] = quantile(ANsim, 0.025)
    upper_AN_daily_heat_fixedmmt[i] = quantile(ANsim, 0.975)
    rm(ANsim)
    
    rm(mod_daily, crosspred)
    
    
    ################################################
    ### Fit aggregated model (Y: weekly, X: daily)
    ################################################
    
    # Weekly counts
    
    Y <- aggregate(Y_sim, by=list(week=dati$week),FUN=sum)
    names(Y)[2] <- "Y_sim"
    
    avgcount <- ave(Y_sim, dati$week)
    
    # Fit model 
    
    try( mod_agr <- fit_aggregate_cpp(Y=Y, X=dati, CB = cross_basis, formula=as.formula(form_mod_sim), family="quasipoisson",
                                      tol = 1e-8, maxit = 500, ntryini = 10, ntryfit = 0))
    
    if (class(mod_agr)[1] == "try-error" | (class(mod_agr)=="modagr" & mod_agr$error_code!=0))  {
      converged_week[i] <- FALSE
    } else {
      if (mod_agr$converged==TRUE) {
        converged_week[i] <- mod_agr$converged
        
        coefs_week[i,] <- coef(mod_agr)
        disp_week[i] <- mod_agr$disp
        
        Ypred_week[,i] <- predict(mod_agr)
        
        if (!(any(!(is.finite(vcov(mod_agr))) ) )) {
          
          crosspred_agr <- crosspred( cross_basis, mod_agr, cen = MMT, at = at_temp, bylag = 1, model.link = "log")
          
          iMMT <- 1 + which( diff( sign( diff( crosspred_agr$allRRfit ) ) ) == 2 ) 
          if( length( iMMT ) > 0 ) { 
            MMT_agr <- crosspred_agr$predvar[iMMT[ which.min( crosspred_agr$allRRfit[ iMMT ] ) ]  ]
          } else {
            MMT_agr <- crosspred_agr$predvar[ which( pred_prc == .1 ) - 1 + which.min( crosspred_agr$allRRfit[ which( pred_prc == 0.1 ) : which( pred_prc == 0.9 ) ] ) ]
          } 
          
          mmt_week[i] <- MMT_agr
          
          crosspred_agr <- crosspred( cross_basis, mod_agr, cen = MMT_agr, at = at_temp, bylag = 1, model.link = "log")
          
          cumpred_week[i,] <- crosspred_agr$allRRfit
          cumpred_low_week[i,] <- crosspred_agr$allRRlow
          cumpred_high_week[i,] <- crosspred_agr$allRRhigh
          
          # Attributable number
          
          # Using fixed MMT (MMTtrue)
          
          crosspred <-  crosspred( cross_basis, coef=(coefs_week[i,])[ind_expos],
                                   vcov=(vcov(mod_agr))[ind_expos, ind_expos],
                                   cen = MMTtrue, at = at_temp, bylag = 1, model.link="log")
          cumpred_week_fixedmmt[i,] <- crosspred$allRRfit
          cumpred_low_week_fixedmmt[i,] <- crosspred$allRRlow
          cumpred_high_week_fixedmmt[i,] <- crosspred$allRRhigh
          
          RRp99lag_week_fixedmmt[i,] <- crosspred$matRRfit[pos_p99,]
          RRp99lag_low_week_fixedmmt[i,] <- crosspred$matRRlow[pos_p99,]
          RRp99lag_high_week_fixedmmt[i,] <- crosspred$matRRhigh[pos_p99,]
          
          RRp01lag_week_fixedmmt[i,] <- crosspred$matRRfit[pos_p01,]
          RRp01lag_low_week_fixedmmt[i,] <- crosspred$matRRlow[pos_p01,]
          RRp01lag_high_week_fixedmmt[i,] <- crosspred$matRRhigh[pos_p01,]
          
          AN_week_fixedmmt[i] <- attrdl(x=dati$tmean, basis=cross_basis, cases=avgcount, 
                                                coef=(coefs_week[i,])[ind_expos],
                                                vcov=(vcov(mod_agr))[ind_expos, ind_expos],
                                                type="an", dir="forw", tot=TRUE, cen=MMTtrue, sim=FALSE)
          ANsim <- attrdl(x=dati$tmean, basis=cross_basis, cases=avgcount, 
                          coef=(coefs_week[i,])[ind_expos],
                          vcov=(vcov(mod_agr))[ind_expos, ind_expos],
                          type="an", dir="forw", tot=TRUE, cen=MMTtrue, sim=TRUE, nsim=1000)
          lower_AN_week_fixedmmt[i] = quantile(ANsim, 0.025)
          upper_AN_week_fixedmmt[i] = quantile(ANsim, 0.975)
          rm(ANsim)
          
          
          AN_week_cold_fixedmmt[i] <- attrdl(x=dati$tmean, basis=cross_basis, cases=avgcount,
                                                     coef=(coefs_week[i,])[ind_expos],
                                                     vcov=(vcov(mod_agr))[ind_expos, ind_expos],
                                                     type="an", dir="forw", tot=TRUE, cen=MMTtrue, sim=FALSE, range=c(-100, MMTtrue))
          ANsim <- attrdl(x=dati$tmean, basis=cross_basis, cases=avgcount,
                          coef=(coefs_week[i,])[ind_expos],
                          vcov=(vcov(mod_agr))[ind_expos, ind_expos],
                          type="an", dir="forw", tot=TRUE, cen=MMTtrue, sim=TRUE, nsim=1000, range=c(-100, MMTtrue))
          lower_AN_week_cold_fixedmmt[i] = quantile(ANsim, 0.025)
          upper_AN_week_cold_fixedmmt[i] = quantile(ANsim, 0.975)
          rm(ANsim)
          
          
          AN_week_heat_fixedmmt[i] <- attrdl(x=dati$tmean, basis=cross_basis, cases=avgcount, 
                                                     coef=(coefs_week[i,])[ind_expos],
                                                     vcov=(vcov(mod_agr))[ind_expos, ind_expos],
                                                     type="an", dir="forw", tot=TRUE, cen=MMTtrue, sim=FALSE, range=c(MMTtrue, 100))
          ANsim <- attrdl(x=dati$tmean, basis=cross_basis, cases=avgcount, 
                          coef=(coefs_week[i,])[ind_expos],
                          vcov=(vcov(mod_agr))[ind_expos, ind_expos],
                          type="an", dir="forw", tot=TRUE, cen=MMTtrue, sim=TRUE, nsim=1000, range=c(MMTtrue, 100))
          lower_AN_week_heat_fixedmmt[i] = quantile(ANsim, 0.025)
          upper_AN_week_heat_fixedmmt[i] = quantile(ANsim, 0.975)
          rm(ANsim)
          
          rm(crosspred_agr, crosspred)  
          
        } else {
          converged_week[i] <- FALSE
        }
      } else {
        converged_week[i] <- FALSE
      }
      
    }
    rm(mod_agr)
    gc()
    
    ################################################
    ### Fit aggregated model (Y: monthly, X: daily)
    ################################################
    
    # Monthly counts
    
    Y <- aggregate(Y_sim, by=list(ymonth=dati$ymonth), FUN=sum)
    names(Y)[2] <- "Y_sim"
    
    avgcount <- ave(Y_sim, dati$ymonth)
    
    
    # Fit model 
    
    try( mod_agr <- fit_aggregate_cpp(Y=Y, X=dati, CB = cross_basis, formula=as.formula(form_mod_sim), family="quasipoisson",
                                      tol = 1e-8, maxit = 500, ntryini = 10, ntryfit = 0))
    if (class(mod_agr)[1] == "try-error" | (class(mod_agr)=="modagr" & mod_agr$error_code!=0))  {
      converged_month[i] <- FALSE
    } else {  
      if (mod_agr$converged==TRUE) {
        
        converged_month[i] <- mod_agr$converged
        
        coefs_month[i,] <- coef(mod_agr)
        disp_month[i] <- mod_agr$disp
        
        Ypred_month[,i] <- predict(mod_agr)
        
        if (!(any(!(is.finite(vcov(mod_agr))) ) )) {
          
          crosspred_agr <- crosspred( cross_basis, mod_agr, cen = MMT, at = at_temp, bylag = 1, model.link = "log")
          
          iMMT <- 1 + which( diff( sign( diff( crosspred_agr$allRRfit ) ) ) == 2 ) 
          if( length( iMMT ) > 0 ) { 
            MMT_agr <- crosspred_agr$predvar[iMMT[ which.min( crosspred_agr$allRRfit[ iMMT ] ) ]  ]
          } else {
            MMT_agr <- crosspred_agr$predvar[ which( pred_prc == .1 ) - 1 + which.min( crosspred_agr$allRRfit[ which( pred_prc == 0.1 ) : which( pred_prc == 0.9 ) ] ) ]
          } 
          
          mmt_month[i] <- MMT_agr
          
          crosspred_agr <- crosspred( cross_basis, mod_agr, cen = MMT_agr, at = at_temp, bylag = 1, model.link = "log")
          
          cumpred_month[i,] <- crosspred_agr$allRRfit
          cumpred_low_month[i,] <- crosspred_agr$allRRlow
          cumpred_high_month[i,] <- crosspred_agr$allRRhigh
          
          
          # Attributable number
          
          # Using fixed MMT (MMTtrue)
          
          crosspred <-  crosspred( cross_basis, coef=(coefs_month[i,])[ind_expos],
                                   vcov=(vcov(mod_agr))[ind_expos, ind_expos],
                                   cen = MMTtrue, at = at_temp, bylag = 1, model.link="log")
          cumpred_month_fixedmmt[i,] <- crosspred$allRRfit
          cumpred_low_month_fixedmmt[i,] <- crosspred$allRRlow
          cumpred_high_month_fixedmmt[i,] <- crosspred$allRRhigh
          
          RRp99lag_month_fixedmmt[i,] <- crosspred$matRRfit[pos_p99,]
          RRp99lag_low_month_fixedmmt[i,] <- crosspred$matRRlow[pos_p99,]
          RRp99lag_high_month_fixedmmt[i,] <- crosspred$matRRhigh[pos_p99,]
          
          RRp01lag_month_fixedmmt[i,] <- crosspred$matRRfit[pos_p01,]
          RRp01lag_low_month_fixedmmt[i,] <- crosspred$matRRlow[pos_p01,]
          RRp01lag_high_month_fixedmmt[i,] <- crosspred$matRRhigh[pos_p01,]
          
          AN_month_fixedmmt[i] <- attrdl(x=dati$tmean, basis=cross_basis, cases=avgcount, 
                                                 coef=(coefs_month[i,])[ind_expos],
                                                 vcov=(vcov(mod_agr))[ind_expos, ind_expos],
                                                 type="an", dir="forw", tot=TRUE, cen=MMTtrue, sim=FALSE)
          ANsim <- attrdl(x=dati$tmean, basis=cross_basis, cases=avgcount, 
                          coef=(coefs_month[i,])[ind_expos],
                          vcov=(vcov(mod_agr))[ind_expos, ind_expos],
                          type="an", dir="forw", tot=TRUE, cen=MMTtrue, sim=TRUE, nsim=1000)
          lower_AN_month_fixedmmt[i] = quantile(ANsim, 0.025)
          upper_AN_month_fixedmmt[i] = quantile(ANsim, 0.975)
          rm(ANsim)
          
          
          AN_month_cold_fixedmmt[i] <- attrdl(x=dati$tmean, basis=cross_basis, cases=avgcount,
                                                      coef=(coefs_month[i,])[ind_expos],
                                                      vcov=(vcov(mod_agr))[ind_expos, ind_expos],
                                                      type="an", dir="forw", tot=TRUE, cen=MMTtrue, sim=FALSE, range=c(-100, MMTtrue))
          ANsim <- attrdl(x=dati$tmean, basis=cross_basis, cases=avgcount,
                          coef=(coefs_month[i,])[ind_expos],
                          vcov=(vcov(mod_agr))[ind_expos, ind_expos],
                          type="an", dir="forw", tot=TRUE, cen=MMTtrue, sim=TRUE, nsim=1000, range=c(-100, MMTtrue))
          lower_AN_month_cold_fixedmmt[i] = quantile(ANsim, 0.025)
          upper_AN_month_cold_fixedmmt[i] = quantile(ANsim, 0.975)
          rm(ANsim)
          
          
          AN_month_heat_fixedmmt[i] <- attrdl(x=dati$tmean, basis=cross_basis, cases=avgcount, 
                                                      coef=(coefs_month[i,])[ind_expos],
                                                      vcov=(vcov(mod_agr))[ind_expos, ind_expos],
                                                      type="an", dir="forw", tot=TRUE, cen=MMTtrue, sim=FALSE, range=c(MMTtrue, 100))
          ANsim <- attrdl(x=dati$tmean, basis=cross_basis, cases=avgcount, 
                          coef=(coefs_month[i,])[ind_expos],
                          vcov=(vcov(mod_agr))[ind_expos, ind_expos],
                          type="an", dir="forw", tot=TRUE, cen=MMTtrue, sim=TRUE, nsim=1000, range=c(MMTtrue, 100))
          lower_AN_month_heat_fixedmmt[i] = quantile(ANsim, 0.025)
          upper_AN_month_heat_fixedmmt[i] = quantile(ANsim, 0.975)
          rm(ANsim)
          
          rm(crosspred_agr, crosspred)  
        } else {
          converged_month[i] <- FALSE
        }
      } else {
        converged_month[i] <- FALSE
      }
      
    }
    rm(mod_agr)
    gc()
    
    ################################################
    ### Fit aggregated model (Y: dow, X: daily)
    ################################################
    
    # dow counts
    
    Y <- aggregate(Y_sim, by=list(ymonthdow=dati$ymonthdow), FUN=sum)
    names(Y)[2] <- "Y_sim"
    
    avgcount <- ave(Y_sim, dati$ymonthdow)
    
    # Fit model 
    
    try( mod_agr <- fit_aggregate_cpp(Y=Y, X=dati, CB = cross_basis, formula=as.formula(form_mod_sim), family="quasipoisson",
                                      tol = 1e-8, maxit = 500, ntryini = 10, ntryfit = 0))
    if (class(mod_agr)[1] == "try-error" | (class(mod_agr)=="modagr" & mod_agr$error_code!=0))  {
      converged_dow[i] <- FALSE
    } else {
      if (mod_agr$converged==TRUE) {
        converged_dow[i] <- mod_agr$converged
        
        coefs_dow[i,] <- coef(mod_agr)
        disp_dow[i] <- mod_agr$disp
        
        Ypred_dow[,i] <- predict(mod_agr)
        
        if (!(any(!(is.finite(vcov(mod_agr))) ) )) {
          
          crosspred_agr <- crosspred( cross_basis, mod_agr, cen = MMT, at = at_temp, bylag = 1, model.link = "log")
          
          iMMT <- 1 + which( diff( sign( diff( crosspred_agr$allRRfit ) ) ) == 2 ) 
          if( length( iMMT ) > 0 ) { 
            MMT_agr <- crosspred_agr$predvar[iMMT[ which.min( crosspred_agr$allRRfit[ iMMT ] ) ]  ]
          } else {
            MMT_agr <- crosspred_agr$predvar[ which( pred_prc == .1 ) - 1 + which.min( crosspred_agr$allRRfit[ which( pred_prc == 0.1 ) : which( pred_prc == 0.9 ) ] ) ]
          } 
          
          mmt_dow[i] <- MMT_agr
          
          crosspred_agr <- crosspred( cross_basis, mod_agr, cen = MMT_agr, at = at_temp, bylag = 1, model.link = "log")
          
          cumpred_dow[i,] <- crosspred_agr$allRRfit
          cumpred_low_dow[i,] <- crosspred_agr$allRRlow
          cumpred_high_dow[i,] <- crosspred_agr$allRRhigh
          
          
          # Attributable number
          
          # Using fixed MMT (MMTtrue)
          
          crosspred <-  crosspred( cross_basis, coef=(coefs_dow[i,])[ind_expos],
                                   vcov=(vcov(mod_agr))[ind_expos, ind_expos],
                                   cen = MMTtrue, at = at_temp, bylag = 1, model.link="log")
          cumpred_dow_fixedmmt[i,] <- crosspred$allRRfit
          cumpred_low_dow_fixedmmt[i,] <- crosspred$allRRlow
          cumpred_high_dow_fixedmmt[i,] <- crosspred$allRRhigh
          
          RRp99lag_dow_fixedmmt[i,] <- crosspred$matRRfit[pos_p99,]
          RRp99lag_low_dow_fixedmmt[i,] <- crosspred$matRRlow[pos_p99,]
          RRp99lag_high_dow_fixedmmt[i,] <- crosspred$matRRhigh[pos_p99,]
          
          RRp01lag_dow_fixedmmt[i,] <- crosspred$matRRfit[pos_p01,]
          RRp01lag_low_dow_fixedmmt[i,] <- crosspred$matRRlow[pos_p01,]
          RRp01lag_high_dow_fixedmmt[i,] <- crosspred$matRRhigh[pos_p01,]
          
          AN_dow_fixedmmt[i] <- attrdl(x=dati$tmean, basis=cross_basis, cases=avgcount, 
                                               coef=(coefs_dow[i,])[ind_expos],
                                               vcov=(vcov(mod_agr))[ind_expos, ind_expos],
                                               type="an", dir="forw", tot=TRUE, cen=MMTtrue, sim=FALSE)
          ANsim <- attrdl(x=dati$tmean, basis=cross_basis, cases=avgcount, 
                          coef=(coefs_dow[i,])[ind_expos],
                          vcov=(vcov(mod_agr))[ind_expos, ind_expos],
                          type="an", dir="forw", tot=TRUE, cen=MMTtrue, sim=TRUE, nsim=1000)
          
          lower_AN_dow_fixedmmt[i] = quantile(ANsim, 0.025)
          upper_AN_dow_fixedmmt[i] = quantile(ANsim, 0.975)
          
          rm(ANsim)
          
          
          AN_dow_cold_fixedmmt[i] <- attrdl(x=dati$tmean, basis=cross_basis, cases=avgcount,
                                                    coef=(coefs_dow[i,])[ind_expos],
                                                    vcov=(vcov(mod_agr))[ind_expos, ind_expos],
                                                    type="an", dir="forw", tot=TRUE, cen=MMTtrue, sim=FALSE, range=c(-100, MMTtrue))
          ANsim <- attrdl(x=dati$tmean, basis=cross_basis, cases=avgcount,
                          coef=(coefs_dow[i,])[ind_expos],
                          vcov=(vcov(mod_agr))[ind_expos, ind_expos],
                          type="an", dir="forw", tot=TRUE, cen=MMTtrue, sim=TRUE, nsim=1000, range=c(-100, MMTtrue))
          lower_AN_dow_cold_fixedmmt[i] = quantile(ANsim, 0.025)
          upper_AN_dow_cold_fixedmmt[i] = quantile(ANsim, 0.975)
          rm(ANsim)
          
          
          AN_dow_heat_fixedmmt[i] <- attrdl(x=dati$tmean, basis=cross_basis, cases=avgcount, 
                                                    coef=(coefs_dow[i,])[ind_expos],
                                                    vcov=(vcov(mod_agr))[ind_expos, ind_expos],
                                                    type="an", dir="forw", tot=TRUE, cen=MMTtrue, sim=FALSE, range=c(MMTtrue, 100))
          ANsim <- attrdl(x=dati$tmean, basis=cross_basis, cases=avgcount, 
                          coef=(coefs_dow[i,])[ind_expos],
                          vcov=(vcov(mod_agr))[ind_expos, ind_expos],
                          type="an", dir="forw", tot=TRUE, cen=MMTtrue, sim=TRUE, nsim=1000, range=c(MMTtrue, 100))
          lower_AN_dow_heat_fixedmmt[i] = quantile(ANsim, 0.025)
          upper_AN_dow_heat_fixedmmt[i] = quantile(ANsim, 0.975)
          rm(ANsim)
          
          rm(crosspred_agr, crosspred)
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
  
  
  return(list(pred_prc, at_temp, nsim, pos_p99, pos_p01, Ysim_mat, 
       converged_week, converged_month, converged_dow,
       coefs_daily, coefs_week, coefs_month, coefs_dow,
       disp_daily, disp_week, disp_month, disp_dow,
       Ypred_daily, Ypred_week, Ypred_month, Ypred_dow, 
       mmt_daily, mmt_week, mmt_month, mmt_dow, 
       cumpred_daily, cumpred_week, cumpred_month, cumpred_dow,
       cumpred_low_daily, cumpred_low_week, cumpred_low_month, cumpred_low_dow,
       cumpred_high_daily, cumpred_high_week, cumpred_high_month, cumpred_high_dow,
       cumpred_daily_fixedmmt, cumpred_week_fixedmmt, cumpred_month_fixedmmt, cumpred_dow_fixedmmt, 
       cumpred_low_daily_fixedmmt, cumpred_low_week_fixedmmt, cumpred_low_month_fixedmmt, cumpred_low_dow_fixedmmt,
       cumpred_high_daily_fixedmmt, cumpred_high_week_fixedmmt, cumpred_high_month_fixedmmt, cumpred_high_dow_fixedmmt,
       RRp99lag_daily_fixedmmt, RRp99lag_week_fixedmmt, RRp99lag_month_fixedmmt, RRp99lag_dow_fixedmmt,
       RRp99lag_low_daily_fixedmmt, RRp99lag_low_week_fixedmmt, RRp99lag_low_month_fixedmmt, RRp99lag_low_dow_fixedmmt,
       RRp99lag_high_daily_fixedmmt, RRp99lag_high_week_fixedmmt, RRp99lag_high_month_fixedmmt, RRp99lag_high_dow_fixedmmt, 
       RRp01lag_daily_fixedmmt, RRp01lag_week_fixedmmt, RRp01lag_month_fixedmmt, RRp01lag_dow_fixedmmt,
       RRp01lag_low_daily_fixedmmt, RRp01lag_low_week_fixedmmt, RRp01lag_low_month_fixedmmt, RRp01lag_low_dow_fixedmmt,
       RRp01lag_high_daily_fixedmmt, RRp01lag_high_week_fixedmmt, RRp01lag_high_month_fixedmmt, RRp01lag_high_dow_fixedmmt,
       AN_daily_fixedmmt, AN_week_fixedmmt, AN_month_fixedmmt, AN_dow_fixedmmt,
       AN_daily_cold_fixedmmt, AN_week_cold_fixedmmt, AN_month_cold_fixedmmt, AN_dow_cold_fixedmmt, 
       AN_daily_heat_fixedmmt, AN_week_heat_fixedmmt, AN_month_heat_fixedmmt, AN_dow_heat_fixedmmt,
       lower_AN_daily_fixedmmt, lower_AN_week_fixedmmt, lower_AN_month_fixedmmt, lower_AN_dow_fixedmmt,
       lower_AN_daily_cold_fixedmmt, lower_AN_week_cold_fixedmmt, lower_AN_month_cold_fixedmmt, lower_AN_dow_cold_fixedmmt, 
       lower_AN_daily_heat_fixedmmt, lower_AN_week_heat_fixedmmt, lower_AN_month_heat_fixedmmt, lower_AN_dow_heat_fixedmmt,
       upper_AN_daily_fixedmmt, upper_AN_week_fixedmmt, upper_AN_month_fixedmmt, upper_AN_dow_fixedmmt,
       upper_AN_daily_cold_fixedmmt, upper_AN_week_cold_fixedmmt, upper_AN_month_cold_fixedmmt, upper_AN_dow_cold_fixedmmt, 
       upper_AN_daily_heat_fixedmmt, upper_AN_week_heat_fixedmmt, upper_AN_month_heat_fixedmmt, upper_AN_dow_heat_fixedmmt))

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
anal = "mort_temp"
nsim <- 500

load(file = "bcn_ts_exposures.RData")


##################################################################
# Specify parameters of the model: outcome- & exposure-specific 
##################################################################

# Read parameters to simulate from

load(file=paste0("parameters_", anal, ".RData"))

mean_tmean <- mean(data$tmean, na.rm=T)
sd_tmean <- sd(data$tmean, na.rm=T)

# Maximum Lag
max_lag <- 21

lag_knots <- logknots( max_lag, 3 )

# Exposure-Response Association
var_fun = "ns";
var_prc = c(0.10,0.75,0.90);
var_knots = quantile(data$tmean, var_prc, na.rm = TRUE )
boundary_knots <- c(min(data$tmean, na.rm=T), max(data$tmean, na.rm=T)) # from full dataset (all years)

# Percentiles for the Predictions of the Cumulative Exposure-Response Association
pred_prc <- c(seq(0, 1, 0.1), 2:98, seq(99, 100, 0.1))/100
at_temp <- quantile(data$tmean, pred_prc, na.rm=T)

pos_p99 <- which(unique(at_temp) == quantile(data$tmean, prob=.99, na.rm=T))
pos_p01 <- which(unique(at_temp) == quantile(data$tmean, prob=.01, na.rm=T))

df_seas <- 8

# Model Formula (for daily data)

rhs = " ~ cross_basis + ns( Date, df = round( df_seas * nrow(data) / 365.25 ))"
rhs_i = " ~ cross_basis + ns( Date, df = round( df_seas * nrow(dati) / 365.25 ))"

form_mod <- paste0( outc, rhs)
form_mod_i <- paste0( outc, rhs_i)

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
                    ind25 = (data$year>=1995))
indyears <- c(1:10, seq(15, 25, 5))


######################
# Start simulation
######################

for (i.years in 1:length(indyears)) {
  cat(",\n", "Years: ", indyears[i.years], ",\n")
  
  # Select years
  dati <- data[inddf[,indyears[i.years]],]
  
  # standardize to global mean and sd
  dati$tmean <-  (dati$tmean - mean(dati$tmean, na.rm=T)) * (sd_tmean/sd(dati$tmean, na.rm=T)) + mean_tmean
  
  # Rebuild crossbasis for the current number of years
  cross_basis <- crossbasis( dati$tmean,
                             lag = c(0,max_lag),
                             argvar = list( fun = var_fun, knots = var_knots,
                                            Boundary.knots = boundary_knots),
                             arglag = list( knots = lag_knots ) )
  
  
  # Put together true coefficients 
  
  coef_true <- matrix(c(coef_int[i.years], coef_expos, coef_rest[[i.years]]), ncol=1)
  
  rlist = run_temp(nsim=nsim, i.years=i.years, dati=dati, crossbasis=crossbasis,
                   outc=outc, form_mod_i=form_mod_i, rhs_i=rhs_i,
                              coef_true=coef_true, disp_par=disp_par,
                              pred_prc=pred_prc, at_temp=at_temp, pos_p99=pos_p99, pos_p01=pos_p01)
    
  save(rlist, file=paste0("results_sim_bcn_v14_", anal, "_", indyears[i.years], ".RData"))
  rm(rlist, dati, coef_true)  
}

