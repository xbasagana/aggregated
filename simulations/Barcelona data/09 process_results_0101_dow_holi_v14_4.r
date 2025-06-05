# Using new results for temperature

library(FluMoDL)
library(Rfast)

#############################################################
# Functions to produce results: bias, rmse, coverage, power, 
#############################################################

# For temperature: 

results_temp = function(data, inddf, indyears, iyearmax, 
                        max_lag, var_fun, var_knots, boundary_knots,
                        lag_knots, coef_expos, vcov_expos, MMTtrue, at_temp, 
                        cumpred, cumpred_low, cumpred_high,
                        mmt,   
                        AF ,
                        AN, AN_heat, AN_cold, 
                        lower_AN, lower_AN_heat, lower_AN_cold,
                        upper_AN, upper_AN_heat, upper_AN_cold,
                        Ysim_mat) {
  
  # Function to calculate the results needed for the plots
  
  #################
  ## Properties ##
  #################
  
  # Calculate true Cumulative E-R curve for all years
  
  truth.daily <- list()
  for (i in 1:iyearmax) {
    
    #dati <- data[inddf[,i],]
    dati <- data[inddf[,indyears[i]],]
    
    # standardize to global mean and sd
    dati$tmean <-  (dati$tmean - mean(dati$tmean, na.rm=T)) * (sd_tmean/sd(dati$tmean, na.rm=T)) + mean_tmean
    
    cross_basis <- crossbasis( dati$tmean,
                               lag = c(0,max_lag),
                               argvar = list( fun = var_fun, knots = var_knots,
                                              Boundary.knots = boundary_knots),
                               arglag = list( knots = lag_knots ) )
    
    crosspred_i <- crosspred( cross_basis, coef=coef_expos,
                              vcov=vcov_expos,
                              cen = MMTtrue, at = at_temp, bylag = 1, model.link="log")
    
    truth.daily[[i]] =  matrix(rep(log(crosspred_i$allRRfit), nrow(cumpred[[i]])),byrow=T,
                               ncol=ncol(cumpred[[i]]))
  }  
  
  pct100 <- c(1, 11:109, 119)
  
  # objects with fixedmmt mean that we used trueMMT for centering (instead of the MMT of each analysis)
  
  bias_cum <- vector()
  bias_perc <- matrix(NA, nrow = length(pct100), ncol = iyearmax )
  coverage_cum <- vector()
  coverage_perc <- matrix(NA, nrow = length(pct100) - 1, ncol = iyearmax )
  rmse_cum <- vector()
  rmse_perc <- matrix(NA, nrow = length(pct100), ncol = iyearmax )
  power01 <- power99 <- vector()
  
  # CIs
  lower_bias_cum <- upper_bias_cum <- vector()
  lower_bias_perc <- upper_bias_perc <- matrix(NA, nrow = length(pct100), ncol = iyearmax )
  
  lower_rmse_cum <- upper_rmse_cum <- vector()
  lower_coverage_cum <- upper_coverage_cum <- vector()
  lower_coverage_perc <- upper_coverage_perc <- matrix(NA, nrow = length(pct100) - 1, ncol = iyearmax )
  
  for (i in 1:iyearmax) {
    bias <- rowMeans(  log(cumpred[[i]][,pct100]) - truth.daily[[i]][,pct100] , na.rm=T)
    bias_cum[i] <- mean(bias, na.rm=T)
    lower_bias_cum[i] <- mean(bias, na.rm=T) - 1.96* sd(bias, na.rm=T)/sqrt(nrow(cumpred[[1]]))
    upper_bias_cum[i] <- mean(bias, na.rm=T) + 1.96* sd(bias, na.rm=T)/sqrt(nrow(cumpred[[1]]))
    
    bias <- log(cumpred[[i]][,pct100]) - truth.daily[[i]][,pct100]
    bias_perc[, i] <- apply(bias, 2, mean, na.rm=TRUE)
    lower_bias_perc[, i] <- bias_perc[, i] - 1.96 * apply(bias, 2, sd, na.rm=TRUE)/sqrt(nrow(bias) )
    upper_bias_perc[, i] <- bias_perc[, i] + 1.96 * apply(bias, 2, sd, na.rm=TRUE)/sqrt(nrow(bias) )
    
    rmse <- sqrt( rowMeans( (log(cumpred[[i]][,pct100]) - truth.daily[[i]][,pct100])^2 , na.rm=T))
    rmse_cum[i] = mean(rmse, na.rm=T)
    lower_rmse_cum[i] = mean(rmse, na.rm=T) - 1.96* sd(rmse, na.rm=T)/sqrt(nsim)
    upper_rmse_cum[i] = mean(rmse, na.rm=T) + 1.96* sd(rmse, na.rm=T)/sqrt(nsim)
  
    rmse_perc[, i] <- sqrt( colMeans( (log(cumpred[[i]][,pct100]) - truth.daily[[i]][,pct100])^2 , na.rm=T) )

    # remove the reference value (RR=1 and no CI)
    # 1) position of reference value
    pos_ref <- apply(cumpred[[i]][,pct100], 1, function(x) which(x==1)[1])
    # 2) exclude reference value from results
    includes_true <- matrix(NA, nrow=nrow(cumpred[[i]]), ncol=length(pct100)-1)
    for (jj in 1:nrow(cumpred[[i]])) {
      includes_true[jj, ] = (log(cumpred_low[[i]][jj, pct100[-pos_ref[[jj]]]]) <  truth.daily[[i]][jj, pct100[-pos_ref[[jj]]]] )  & 
        (log(cumpred_high[[i]][jj ,pct100[-pos_ref[[jj]]]]) > truth.daily[[i]][jj, pct100[-pos_ref[[jj]]]] ) 
    }
    # Calculate coverage in each percentile
    cover <- colMeans(includes_true, na.rm=T)
    # Average coverage
    coverage_cum[i] = mean(cover, na.rm=T) 
    lower_coverage_cum[i] = mean(cover, na.rm=T) - 1.96 * sqrt(mean(cover, na.rm=T)*(1-mean(cover, na.rm=T))/nsim)
    upper_coverage_cum[i] = mean(cover, na.rm=T) + 1.96 * sqrt(mean(cover, na.rm=T)*(1-mean(cover, na.rm=T))/nsim)
    
    coverage_perc[, i] <- cover
    lower_coverage_perc[, i] <- cover - 1.96 * sqrt( cover * (1 - cover) / nsim ) 
    upper_coverage_perc[, i] <- cover + 1.96 * sqrt( cover * (1 - cover) / nsim ) 
    
    power01[i] <-  mean( cumpred_low[[i]][,pos_p01] > 1, na.rm=T )
    
    power99[i] <-  mean( cumpred_low[[i]][,pos_p99] > 1, na.rm=T)
  }
  
  
  # MMT
  #########
  
  bias_mmt <- vector()
  lower_bias_mmt <- vector()
  upper_bias_mmt <- vector()
  
  rmse_mmt <- vector()
  lower_rmse_mmt <- vector()
  upper_rmse_mmt <- vector()
  
  for (i in 1:iyearmax) {
    
    bias = mmt[[i]] - MMTtrue
    bias_mmt[i] = mean(bias)
    lower_bias_mmt[i] = mean(bias) - 1.96* sd(bias, na.rm=T)/sqrt(nsim)
    upper_bias_mmt[i] = mean(bias) + 1.96* sd(bias, na.rm=T)/sqrt(nsim)
    
    rmse_mmt[i] = sqrt(mean( (mmt[[i]] - MMTtrue)^2 , na.rm=T))
    
    # https://stats.stackexchange.com/questions/78079/confidence-interval-of-rmse 
    
    # I need to add coverage here
  }
  
  
  # AF
  #######
  
  if (AF==FALSE) {
    AF <- AN
    AF_cold <- AN_cold
    AF_heat <- AN_heat
    
    lower_AF <- lower_AN
    lower_AF_cold <- lower_AN_cold
    lower_AF_heat <- lower_AN_heat
    upper_AF <- upper_AN
    upper_AF_cold <- upper_AN_cold
    upper_AF_heat <- upper_AN_heat
    
    AN_true <- vector()
    AF_true <- vector()
    AN_true_cold <- vector()
    AF_true_cold <- vector()
    AN_true_heat <- vector()
    AF_true_heat <- vector()
    
    for (i in 1:iyearmax) {
      dati <- data[inddf[,indyears[i]],]
      
      # standardize to global mean and sd
      dati$tmean <-  (dati$tmean - mean(dati$tmean, na.rm=T) ) * (sd_tmean/sd(dati$tmean, na.rm=T)) + mean_tmean
      
      cross_basis <- crossbasis( dati$tmean,
                                 lag = c(0,max_lag),
                                 argvar = list( fun = var_fun, knots = var_knots,
                                                Boundary.knots = boundary_knots),
                                 arglag = list( knots = lag_knots ) )
      
      # As I change the model equation used to simulate the data by forcing the same coefficients of 
      #   the crossbasis in each period, the observed mortality counts no longer represent the
      #   reality I simulate from. To obtain the "true" daily number of cases, I take the average
      #   number over all simulations
      
      daily_cases <- rowmeans(Ysim_mat[[i]])
      
      AN_true[i] <- attrdl(x=dati$tmean, basis=cross_basis, cases=daily_cases, coef=coef_expos,
                           vcov=vcov_expos,  type="an", dir="forw", tot=TRUE, cen=MMTtrue, 
                           sim=FALSE)
      AF_true[i] <- 100*AN_true[i]/sum(daily_cases,na.rm=T)
      
      AN_true_cold[i] <- attrdl(x=dati$tmean, basis=cross_basis, cases=daily_cases, coef=coef_expos,
                                vcov=vcov_expos, type="an", dir="forw", tot=TRUE, cen=MMTtrue,
                                sim=FALSE, range=c(-100, MMTtrue))
      AF_true_cold[i] <- 100*AN_true_cold[i]/sum(daily_cases,na.rm=T)
      
      AN_true_heat[i] <- attrdl(x=dati$tmean, basis=cross_basis, cases=daily_cases, coef=coef_expos,
                                vcov=vcov_expos, type="an", dir="forw", tot=TRUE, cen=MMTtrue,
                                sim=FALSE, range=c(MMTtrue, 100))
      AF_true_heat[i] <- 100*AN_true_heat[i]/sum(daily_cases,na.rm=T)
      
      for (j in 1:nsim) {
        
        AF[[i]][j] <- 100*AN[[i]][j]/sum(Ysim_mat[[i]][,j],na.rm=T)
        AF_cold[[i]][j] <- 100*AN_cold[[i]][j]/sum(Ysim_mat[[i]][,j],na.rm=T)
        AF_heat[[i]][j] <- 100*AN_heat[[i]][j]/sum(Ysim_mat[[i]][,j],na.rm=T)
        
        lower_AF[[i]][j] <- 100*lower_AN[[i]][j]/sum(Ysim_mat[[i]][,j],na.rm=T)
        lower_AF_cold[[i]][j] <- 100*lower_AN_cold[[i]][j]/sum(Ysim_mat[[i]][,j],na.rm=T)
        lower_AF_heat[[i]][j] <- 100*lower_AN_heat[[i]][j]/sum(Ysim_mat[[i]][,j],na.rm=T)
        
        upper_AF[[i]][j] <- 100*upper_AN[[i]][j]/sum(Ysim_mat[[i]][,j],na.rm=T)
        upper_AF_cold[[i]][j] <- 100*upper_AN_cold[[i]][j]/sum(Ysim_mat[[i]][,j],na.rm=T)
        upper_AF_heat[[i]][j] <- 100*upper_AN_heat[[i]][j]/sum(Ysim_mat[[i]][,j],na.rm=T)
        
      }
      
    }
  }
  
  bias_AF <- vector()
  lower_bias_AF <- vector()
  upper_bias_AF <- vector()
  
  rmse_AF <- vector()
  
  bias_AF_cold <- vector()
  lower_bias_AF_cold <- vector()
  upper_bias_AF_cold <- vector()
  
  rmse_AF_cold <- vector()
  
  bias_AF_heat <- vector()
  lower_bias_AF_heat <- vector()
  upper_bias_AF_heat <- vector()
  
  rmse_AF_heat <- vector()
  
  
  cover_AF <- vector()
  lower_cover_AF <- vector()
  upper_cover_AF <- vector()
  
  cover_AF_cold <- vector()
  lower_cover_AF_cold <- vector()
  upper_cover_AF_cold <- vector()
  
  cover_AF_heat <- vector()
  lower_cover_AF_heat <- vector()
  upper_cover_AF_heat <- vector()
  
  for (i in 1:iyearmax) {
    
    bias = AF[[i]] - AF_true[i]
    bias_AF[i] = mean(bias, na.rm=T)
    lower_bias_AF[i] = mean(bias, na.rm=T) - 1.96* sd(bias, na.rm=T)/sqrt(nsim)
    upper_bias_AF[i] = mean(bias, na.rm=T) + 1.96* sd(bias, na.rm=T)/sqrt(nsim)
    
    rmse_AF[i] = sqrt(mean( (AF[[i]] - AF_true[i])^2 , na.rm=T))
    
    bias =  AF_cold[[i]] - AF_true_cold[i]
    bias_AF_cold[i] = mean( bias, na.rm=T)
    lower_bias_AF_cold[i] = mean( bias, na.rm=T) - 1.96* sd(bias, na.rm=T)/sqrt(nsim)
    upper_bias_AF_cold[i] = mean( bias, na.rm=T) + 1.96* sd(bias, na.rm=T)/sqrt(nsim)
    
    rmse_AF_cold[i] = sqrt(mean( (AF_cold[[i]] - AF_true_cold[i])^2 , na.rm=T))
    
    bias = AF_heat[[i]] - AF_true_heat[i]
    bias_AF_heat[i] = mean(bias, na.rm=T)
    lower_bias_AF_heat[i] = mean(bias, na.rm=T) - 1.96* sd(bias, na.rm=T)/sqrt(nsim)
    upper_bias_AF_heat[i] = mean(bias, na.rm=T) + 1.96* sd(bias, na.rm=T)/sqrt(nsim)
    
    
    rmse_AF_heat[i] = sqrt(mean( (AF_heat[[i]] - AF_true_heat[i])^2 , na.rm=T))
    
    cover <- ((lower_AF[[i]] <= AF_true[i]) & (upper_AF[[i]] >= AF_true[i]))
    cover_AF[i] = mean(cover, na.rm=T) 
    lower_cover_AF[i] = mean(cover, na.rm=T) - 1.96 * sqrt(mean(cover, na.rm=T)*(1-mean(cover, na.rm=T))/nsim)
    upper_cover_AF[i] = mean(cover, na.rm=T) + 1.96 * sqrt(mean(cover, na.rm=T)*(1-mean(cover, na.rm=T))/nsim)
    
    cover = ((lower_AF_cold[[i]] <= AF_true_cold[i]) & (upper_AF_cold[[i]] >= AF_true_cold[i]))
    cover_AF_cold[i] = mean(cover, na.rm=T) 
    lower_cover_AF_cold[i] = mean(cover, na.rm=T) - 1.96 * sqrt(mean(cover, na.rm=T)*(1-mean(cover, na.rm=T))/nsim)
    upper_cover_AF_cold[i] = mean(cover, na.rm=T) + 1.96 * sqrt(mean(cover, na.rm=T)*(1-mean(cover, na.rm=T))/nsim)
    
    cover = ((lower_AF_heat[[i]] <= AF_true_heat[i]) & (upper_AF_heat[[i]] >= AF_true_heat[i]))
    cover_AF_heat[i] = mean(cover, na.rm=T) 
    lower_cover_AF_heat[i] = mean(cover, na.rm=T) - 1.96 * sqrt(mean(cover, na.rm=T)*(1-mean(cover, na.rm=T))/nsim)
    upper_cover_AF_heat[i] = mean(cover, na.rm=T) + 1.96 * sqrt(mean(cover, na.rm=T)*(1-mean(cover, na.rm=T))/nsim)
    
    
  }
  
  return(list(bias_cum = bias_cum, coverage_cum = coverage_cum, rmse_cum = rmse_cum,
              lower_bias_cum = lower_bias_cum, upper_bias_cum = upper_bias_cum, 
              lower_rmse_cum = lower_rmse_cum, upper_rmse_cum = upper_rmse_cum,
              bias_perc = bias_perc, lower_bias_perc = lower_bias_perc, upper_bias_perc = upper_bias_perc,
              rmse_perc = rmse_perc,
              lower_coverage_cum = lower_coverage_cum, upper_coverage_cum = upper_coverage_cum,
              coverage_perc = coverage_perc, lower_coverage_perc = lower_coverage_perc, upper_coverage_perc = upper_coverage_perc,
              power99 = power99, power01 = power01,
              truth.daily = truth.daily,
              bias_mmt = bias_mmt, lower_bias_mmt = lower_bias_mmt, upper_bias_mmt = upper_bias_mmt,
              rmse_mmt = rmse_mmt, lower_rmse_mmt = lower_rmse_mmt, upper_rmse_mmt = upper_rmse_mmt,
              AF = AF, AF_cold = AF_cold, AF_heat = AF_heat,
              AN_true = AN_true, AF_true = AF_true,
              AN_true_cold = AN_true_cold, AF_true_cold = AF_true_cold,
              AN_true_heat = AN_true_heat, AF_true_heat = AF_true_heat,
              bias_AF = bias_AF, lower_bias_AF = lower_bias_AF, upper_bias_AF = upper_bias_AF,
              rmse_AF = rmse_AF, bias_AF_cold = bias_AF_cold, 
              lower_bias_AF_cold = lower_bias_AF_cold, upper_bias_AF_cold = upper_bias_AF_cold,
              rmse_AF_cold = rmse_AF_cold, 
              bias_AF_heat = bias_AF_heat, lower_bias_AF_heat = lower_bias_AF_heat,
              upper_bias_AF_heat = upper_bias_AF_heat,
              rmse_AF_heat = rmse_AF_heat,
              cover_AF = cover_AF, lower_cover_AF=lower_cover_AF, upper_cover_AF=upper_cover_AF,
              cover_AF_cold=cover_AF_cold, lower_cover_AF_cold=lower_cover_AF_cold, 
              upper_cover_AF_cold=upper_cover_AF_cold,
              cover_AF_heat=cover_AF_heat, lower_cover_AF_heat=lower_cover_AF_heat, 
              upper_cover_AF_heat=upper_cover_AF_heat))
}


# For 1 parameter (NO2)

results_single = function(data, inddf, indyears, iyearmax, name_expos, ind_expos, coef_expos, vcov_expos, at_no2, 
                          coefs, ses, AF, lower_AF, upper_AF, Ysim_mat) {
  # function to calculate the results needed for the plots for the case where the exposure is
  #   a single parameter
  
  # objects with fixedmmt mean that we used trueMMT for centering (instead of the MMT of each analysis)
  
  bias_cum <- vector()
  coverage_cum <- vector()
  rmse_cum <- vector()
  power <- vector()
  
  # CIs
  lower_bias_cum <- upper_bias_cum <- vector()
  lower_rmse_cum <- upper_rmse_cum <- vector()
  lower_coverage_cum <- upper_coverage_cum <- vector()
  
  for (i in 1:length(coefs)) {
    bias <- coefs[[i]][,ind_expos] - coef_expos
    bias_cum[i] <- mean(bias, na.rm=T)
    lower_bias_cum[i] <- mean(bias, na.rm=T) - 1.96* sd(bias, na.rm=T)/sqrt(nrow(coefs[[1]]))
    upper_bias_cum[i] <- mean(bias, na.rm=T) + 1.96* sd(bias, na.rm=T)/sqrt(nrow(coefs[[1]]))
    
    rmse <- sqrt( mean((coefs[[i]][,ind_expos] - coef_expos )^2, na.rm=T) )
    rmse_cum[i] = mean(rmse, na.rm=T)
    lower_rmse_cum[i] = mean(rmse, na.rm=T) - 1.96* sd(rmse, na.rm=T)/sqrt(nsim)
    upper_rmse_cum[i] = mean(rmse, na.rm=T) + 1.96* sd(rmse, na.rm=T)/sqrt(nsim)
    
    lower = vector()
    upper = vector()
    
    #if (class(ses[[1]])[1]=="list") {
    #  for (i.sim in 1:nsim) {
    #    lower[i.sim] = coefs[[i]][i.sim,ind_expos] - 1.96*sqrt( (vcovs[[i]][[i.sim]])[ind_expos,ind_expos])
    #    upper[i.sim] = coefs[[i]][i.sim,ind_expos] + 1.96*sqrt( (vcovs[[i]][[i.sim]])[ind_expos,ind_expos])
    #  }
    #} else {
      for (i.sim in 1:nsim) {
        lower[i.sim] = coefs[[i]][i.sim,ind_expos] - 1.96*ses[[i]][i.sim,ind_expos]
        upper[i.sim] = coefs[[i]][i.sim,ind_expos] + 1.96*ses[[i]][i.sim,ind_expos]    
      }
    #}
    
    cover <- (lower < coef_expos) & (upper > coef_expos)
    coverage_cum[i] = mean(cover, na.rm=T) 
    lower_coverage_cum[i] = mean(cover, na.rm=T) - 1.96 * sqrt(mean(cover, na.rm=T)*(1-mean(cover, na.rm=T))/nsim)
    upper_coverage_cum[i] = mean(cover, na.rm=T) + 1.96 * sqrt(mean(cover, na.rm=T)*(1-mean(cover, na.rm=T))/nsim)
    
    power[i] <- mean( lower > 0, na.rm=T)
  }
  
  AF_true = vector()
  for (i in 1:iyearmax) {
    dati <- data[inddf[,indyears[i]],]
    
    # As I change the model equation used to simulate the data by forcing the same coefficients of 
    #   the crossbasis in each period, the observed mortality counts no longer represent the
    #   reality I simulate from. To obtain the "true" daily number of cases, I take the average
    #   number over all simulations
    
    daily_cases <- rowmeans(Ysim_mat[[i]])
    
    RRj = exp(coef_expos*dati[, name_expos])
    AF_true[i] <- sum( (daily_cases * (RRj-1)/RRj) / sum(daily_cases, na.rm=T), na.rm=T)
    
  }  
  
  bias_AF <- vector()
  lower_bias_AF <- vector()
  upper_bias_AF <- vector()
  
  rmse_AF <- vector()
  
  cover_AF <- vector()
  lower_cover_AF <- vector()
  upper_cover_AF <- vector()
  
  for (i in 1:iyearmax) {
    
    bias = AF[[i]] - AF_true[i]
    bias_AF[i] = mean(bias, na.rm=T)
    lower_bias_AF[i] = mean(bias, na.rm=T) - 1.96* sd(bias, na.rm=T)/sqrt(nsim)
    upper_bias_AF[i] = mean(bias, na.rm=T) + 1.96* sd(bias, na.rm=T)/sqrt(nsim)
    
    rmse_AF[i] = sqrt(mean( (AF[[i]] - AF_true[i])^2 , na.rm=T))
    
    cover <- ((lower_AF[[i]] <= AF_true[i]) & (upper_AF[[i]] >= AF_true[i]))
    cover_AF[i] = mean(cover, na.rm=T) 
    lower_cover_AF[i] = mean(cover, na.rm=T) - 1.96 * sqrt(mean(cover, na.rm=T)*(1-mean(cover, na.rm=T))/nsim)
    upper_cover_AF[i] = mean(cover, na.rm=T) + 1.96 * sqrt(mean(cover, na.rm=T)*(1-mean(cover, na.rm=T))/nsim)
    
  }
  
  
  return(list(bias_cum = bias_cum, coverage_cum = coverage_cum, rmse_cum = rmse_cum,
              lower_bias_cum = lower_bias_cum, upper_bias_cum = upper_bias_cum, 
              lower_rmse_cum = lower_rmse_cum, upper_rmse_cum = upper_rmse_cum,
              lower_coverage_cum = lower_coverage_cum, upper_coverage_cum = upper_coverage_cum,
              bias_AF=bias_AF, lower_bias_AF=lower_bias_AF, upper_bias_AF=upper_bias_AF,
              rmse_AF=rmse_AF, cover_AF=cover_AF, lower_cover_AF=lower_cover_AF,
              upper_cover_AF=upper_cover_AF,
              power=power))
  
}



indyears <- c(1:10, seq(15, 25, 5))

########################################################
# Append individual results files into a single list
########################################################

## Temperature results
##########################

## Mortality 

anal = "mort_temp"

k = 1
#for (i in 1:length(indyears)) {
for (i in 1:13) {
  load(file=paste0("results_sim_bcn_v14_", anal, "_", indyears[i], ".RData"))
  if (i==1) {
    lmort_temp = rlist
    names(lmort_temp) = c("pred_prc", "at_temp", "nsim", "pos_p99", "pos_p01", 
                          "Ysim_mat", "converged_week", "converged_month", "converged_dow",
                          "coefs_daily", "coefs_week", "coefs_month", "coefs_dow",
                          "disp_daily", "disp_week", "disp_month", "disp_dow",
                          "Ypred_daily", "Ypred_week", "Ypred_month", "Ypred_dow", 
                          "mmt_daily", "mmt_week", "mmt_month", "mmt_dow",
                          "cumpred_daily", "cumpred_week", "cumpred_month", "cumpred_dow",
                          "cumpred_low_daily", "cumpred_low_week", "cumpred_low_month", "cumpred_low_dow",
                          "cumpred_high_daily", "cumpred_high_week", "cumpred_high_month", "cumpred_high_dow",
                          "cumpred_daily_fixedmmt", "cumpred_week_fixedmmt", "cumpred_month_fixedmmt", "cumpred_dow_fixedmmt", 
                          "cumpred_low_daily_fixedmmt", "cumpred_low_week_fixedmmt", "cumpred_low_month_fixedmmt", "cumpred_low_dow_fixedmmt",
                          "cumpred_high_daily_fixedmmt", "cumpred_high_week_fixedmmt", "cumpred_high_month_fixedmmt", "cumpred_high_dow_fixedmmt",
                          "RRp99lag_daily_fixedmmt", "RRp99lag_week_fixedmmt", "RRp99lag_month_fixedmmt", "RRp99lag_dow_fixedmmt",
                          "RRp99lag_low_daily_fixedmmt", "RRp99lag_low_week_fixedmmt", "RRp99lag_low_month_fixedmmt", "RRp99lag_low_dow_fixedmmt",
                          "RRp99lag_high_daily_fixedmmt", "RRp99lag_high_week_fixedmmt", "RRp99lag_high_month_fixedmmt", "RRp99lag_high_dow_fixedmmt", 
                          "RRp01lag_daily_fixedmmt", "RRp01lag_week_fixedmmt", "RRp01lag_month_fixedmmt", "RRp01lag_dow_fixedmmt",
                          "RRp01lag_low_daily_fixedmmt", "RRp01lag_low_week_fixedmmt", "RRp01lag_low_month_fixedmmt", "RRp01lag_low_dow_fixedmmt",
                          "RRp01lag_high_daily_fixedmmt", "RRp01lag_high_week_fixedmmt", "RRp01lag_high_month_fixedmmt", "RRp01lag_high_dow_fixedmmt",
                          "AN_daily_fixedmmt", "AN_week_fixedmmt", "AN_month_fixedmmt", "AN_dow_fixedmmt",
                          "AN_daily_cold_fixedmmt", "AN_week_cold_fixedmmt", "AN_month_cold_fixedmmt", "AN_dow_cold_fixedmmt", 
                          "AN_daily_heat_fixedmmt", "AN_week_heat_fixedmmt", "AN_month_heat_fixedmmt", "AN_dow_heat_fixedmmt",
                          "lower_AN_daily_fixedmmt", "lower_AN_week_fixedmmt", "lower_AN_month_fixedmmt", "lower_AN_dow_fixedmmt",
                          "lower_AN_daily_cold_fixedmmt", "lower_AN_week_cold_fixedmmt", "lower_AN_month_cold_fixedmmt", "lower_AN_dow_cold_fixedmmt", 
                          "lower_AN_daily_heat_fixedmmt", "lower_AN_week_heat_fixedmmt", "lower_AN_month_heat_fixedmmt", "lower_AN_dow_heat_fixedmmt",
                          "upper_AN_daily_fixedmmt", "upper_AN_week_fixedmmt", "upper_AN_month_fixedmmt", "upper_AN_dow_fixedmmt",
                          "upper_AN_daily_cold_fixedmmt", "upper_AN_week_cold_fixedmmt", "upper_AN_month_cold_fixedmmt", "upper_AN_dow_cold_fixedmmt",
                          "upper_AN_daily_heat_fixedmmt", "upper_AN_week_heat_fixedmmt", "upper_AN_month_heat_fixedmmt", "upper_AN_dow_heat_fixedmmt")
      
    for (j in 2:length(lmort_temp)) {
      lmort_temp[[j]] = list(lmort_temp[[j]])
    }
  } else {
    lmort_temp[[1]] = c(lmort_temp[[1]], rlist[[1]])
    for (j in 2:length(lmort_temp)) {
      lmort_temp[[j]][[k]] = rlist[[j]] 
    }
  }
  k = k + 1
}

rm(rlist)



## Hospitalizations 

anal = "hosp_temp"

k = 1
#for (i in 1:length(indyears)) {
for (i in 1:13) {
  load(file=paste0("results_sim_bcn_v14_", anal, "_", indyears[i], ".RData"))
  if (i==1) {
    lhosp_temp = rlist
    names(lhosp_temp) = c("pred_prc", "at_temp", "nsim", "pos_p99", "pos_p01", 
                          "Ysim_mat", "converged_week", "converged_month", "converged_dow",
                          "coefs_daily", "coefs_week", "coefs_month", "coefs_dow",
                          "disp_daily", "disp_week", "disp_month", "disp_dow",
                          "Ypred_daily", "Ypred_week", "Ypred_month", "Ypred_dow", 
                          "mmt_daily", "mmt_week", "mmt_month", "mmt_dow",
                          "cumpred_daily", "cumpred_week", "cumpred_month", "cumpred_dow",
                          "cumpred_low_daily", "cumpred_low_week", "cumpred_low_month", "cumpred_low_dow",
                          "cumpred_high_daily", "cumpred_high_week", "cumpred_high_month", "cumpred_high_dow",
                          "cumpred_daily_fixedmmt", "cumpred_week_fixedmmt", "cumpred_month_fixedmmt", "cumpred_dow_fixedmmt", 
                          "cumpred_low_daily_fixedmmt", "cumpred_low_week_fixedmmt", "cumpred_low_month_fixedmmt", "cumpred_low_dow_fixedmmt",
                          "cumpred_high_daily_fixedmmt", "cumpred_high_week_fixedmmt", "cumpred_high_month_fixedmmt", "cumpred_high_dow_fixedmmt",
                          "RRp99lag_daily_fixedmmt", "RRp99lag_week_fixedmmt", "RRp99lag_month_fixedmmt", "RRp99lag_dow_fixedmmt",
                          "RRp99lag_low_daily_fixedmmt", "RRp99lag_low_week_fixedmmt", "RRp99lag_low_month_fixedmmt", "RRp99lag_low_dow_fixedmmt",
                          "RRp99lag_high_daily_fixedmmt", "RRp99lag_high_week_fixedmmt", "RRp99lag_high_month_fixedmmt", "RRp99lag_high_dow_fixedmmt", 
                          "RRp01lag_daily_fixedmmt", "RRp01lag_week_fixedmmt", "RRp01lag_month_fixedmmt", "RRp01lag_dow_fixedmmt",
                          "RRp01lag_low_daily_fixedmmt", "RRp01lag_low_week_fixedmmt", "RRp01lag_low_month_fixedmmt", "RRp01lag_low_dow_fixedmmt",
                          "RRp01lag_high_daily_fixedmmt", "RRp01lag_high_week_fixedmmt", "RRp01lag_high_month_fixedmmt", "RRp01lag_high_dow_fixedmmt",
                          "AN_daily_fixedmmt", "AN_week_fixedmmt", "AN_month_fixedmmt", "AN_dow_fixedmmt",
                          "AN_daily_cold_fixedmmt", "AN_week_cold_fixedmmt", "AN_month_cold_fixedmmt", "AN_dow_cold_fixedmmt", 
                          "AN_daily_heat_fixedmmt", "AN_week_heat_fixedmmt", "AN_month_heat_fixedmmt", "AN_dow_heat_fixedmmt",
                          "lower_AN_daily_fixedmmt", "lower_AN_week_fixedmmt", "lower_AN_month_fixedmmt", "lower_AN_dow_fixedmmt",
                          "lower_AN_daily_cold_fixedmmt", "lower_AN_week_cold_fixedmmt", "lower_AN_month_cold_fixedmmt", "lower_AN_dow_cold_fixedmmt", 
                          "lower_AN_daily_heat_fixedmmt", "lower_AN_week_heat_fixedmmt", "lower_AN_month_heat_fixedmmt", "lower_AN_dow_heat_fixedmmt",
                          "upper_AN_daily_fixedmmt", "upper_AN_week_fixedmmt", "upper_AN_month_fixedmmt", "upper_AN_dow_fixedmmt",
                          "upper_AN_daily_cold_fixedmmt", "upper_AN_week_cold_fixedmmt", "upper_AN_month_cold_fixedmmt", "upper_AN_dow_cold_fixedmmt",
                          "upper_AN_daily_heat_fixedmmt", "upper_AN_week_heat_fixedmmt", "upper_AN_month_heat_fixedmmt", "upper_AN_dow_heat_fixedmmt")
    for (j in 2:length(lhosp_temp)) {
      lhosp_temp[[j]] = list(lhosp_temp[[j]])
    }
  } else {
    lhosp_temp[[1]] = c(lhosp_temp[[1]], rlist[[1]])
    for (j in 2:length(lhosp_temp)) {
      lhosp_temp[[j]][[k]] = rlist[[j]] 
    }
  }
  k = k + 1
}

rm(rlist)


## Air pollution results
##########################

## Mortality 

anal = "mort_ap"

k = 1
#for (i in 1:length(indyears)) {
for (i in 1:13) {
  load(file=paste0("results_sim_bcn_v14_", anal, "_dow_holi_0101_", indyears[i], ".RData"))
  if (i==1) {
    lmort_ap = rlist
    names(lmort_ap) = c("nsim", "Ysim_mat", 
                        "converged_week", "converged_month", "converged_dow",
                        "coefs_daily", "coefs_week", "coefs_month", "coefs_dow",
                        "ses_daily", "ses_week", "ses_month", "ses_dow",       
                        "Ypred_daily", "Ypred_week", "Ypred_month", "Ypred_dow",  
                        "AF_daily", "AF_week", "AF_month", "AF_dow", 
                        "lower_AF_daily", "lower_AF_week", "lower_AF_month", "lower_AF_dow", 
                        "upper_AF_daily", "upper_AF_week", "upper_AF_month", "upper_AF_dow", 
                        "disp_daily", "disp_week", "disp_month", "disp_dow")
    for (j in 2:length(lmort_ap)) {
      lmort_ap[[j]] = list(lmort_ap[[j]])
    }
  } else {
    lmort_ap[[1]] = c(lmort_ap[[1]], rlist[[1]])
    for (j in 2:length(lmort_ap)) {
      lmort_ap[[j]][[k]] = rlist[[j]] 
    }
  }
  k = k + 1
}

rm(rlist)



## Hospitalizations 

anal = "hosp_ap"

k = 1
#for (i in 1:length(indyears)) {
for (i in 1:13) {
  load(file=paste0("results_sim_bcn_v14_", anal, "_dow_holi_0101_1pct_", indyears[i], ".RData"))
  if (i==1) {
    lhosp_ap = rlist
    names(lhosp_ap) = c("nsim", "Ysim_mat", 
                        "converged_week", "converged_month", "converged_dow",
                        "coefs_daily", "coefs_week", "coefs_month", "coefs_dow",
                        "ses_daily", "ses_week", "ses_month", "ses_dow",       
                        "Ypred_daily", "Ypred_week", "Ypred_month", "Ypred_dow",  
                        "AF_daily", "AF_week", "AF_month", "AF_dow", 
                        "lower_AF_daily", "lower_AF_week", "lower_AF_month", "lower_AF_dow", 
                        "upper_AF_daily", "upper_AF_week", "upper_AF_month", "upper_AF_dow", 
                        "disp_daily", "disp_week", "disp_month", "disp_dow")
    for (j in 2:length(lhosp_ap)) {
      lhosp_ap[[j]] = list(lhosp_ap[[j]])
    }
  } else {
    lhosp_ap[[1]] = c(lhosp_ap[[1]], rlist[[1]])
    for (j in 2:length(lhosp_ap)) {
      lhosp_ap[[j]][[k]] = rlist[[j]] 
    }
  }
  k = k + 1
}

rm(rlist)


###############################################
# Load original data and model specifications
###############################################

load(file = "bcn_ts_exposures.RData")

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

at_no2 <- seq(0, max(data$no2_lag06, na.rm=T), by=0.1)

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



###########################################################
# Obtain results from data: bias, RMSE, coverage, power, 
###########################################################

nsim <- 500

results_daily = list()

### Mortality - temperature
###############################

iyearmax = 13

anal = "mort_temp"
load(file=paste0("parameters_", anal, ".RData"))


results_daily[[1]] <- results_temp(data=data, inddf=inddf, indyears=indyears, 
                                   iyearmax=iyearmax, max_lag=max_lag, var_fun=var_fun,
                                   var_knots=var_knots, boundary_knots=boundary_knots,
                                   lag_knots=lag_knots, coef_expos=coef_expos, 
                                   vcov_expos=vcov_expos, MMTtrue=MMTtrue, at_temp=at_temp,
                                   cumpred= lmort_temp$cumpred_daily_fixedmmt , 
                                   cumpred_low=lmort_temp$cumpred_low_daily_fixedmmt,
                                   cumpred_high=lmort_temp$cumpred_high_daily_fixedmmt,
                                   mmt=lmort_temp$mmt_daily,
                                   AF=FALSE, AN = lmort_temp$AN_daily_fixedmmt, 
                                   AN_heat=lmort_temp$AN_daily_heat_fixedmmt,
                                   AN_cold=lmort_temp$AN_daily_cold_fixedmmt,
                                   lower_AN=lmort_temp$lower_AN_daily_fixedmmt,
                                   lower_AN_heat=lmort_temp$lower_AN_daily_heat_fixedmmt,
                                   lower_AN_cold=lmort_temp$lower_AN_daily_cold_fixedmmt,
                                   upper_AN=lmort_temp$upper_AN_daily_fixedmmt,
                                   upper_AN_heat=lmort_temp$upper_AN_daily_heat_fixedmmt,
                                   upper_AN_cold=lmort_temp$upper_AN_daily_cold_fixedmmt,
                                   Ysim_mat=lmort_temp$Ysim_mat) 

results_week = list()
results_week[[1]] = results_temp(data=data, inddf=inddf, indyears=indyears, 
                                 iyearmax=iyearmax, max_lag=max_lag, var_fun=var_fun,
                                 var_knots=var_knots, boundary_knots=boundary_knots,
                                 lag_knots=lag_knots, coef_expos=coef_expos, 
                                 vcov_expos=vcov_expos, MMTtrue=MMTtrue, at_temp=at_temp,
                                 cumpred= lmort_temp$cumpred_week_fixedmmt , 
                                 cumpred_low=lmort_temp$cumpred_low_week_fixedmmt,
                                 cumpred_high=lmort_temp$cumpred_high_week_fixedmmt,
                                 mmt=lmort_temp$mmt_week,
                                 AF=FALSE, AN = lmort_temp$AN_week_fixedmmt, 
                                 AN_heat=lmort_temp$AN_week_heat_fixedmmt,
                                 AN_cold=lmort_temp$AN_week_cold_fixedmmt,
                                 lower_AN=lmort_temp$lower_AN_week_fixedmmt,
                                 lower_AN_heat=lmort_temp$lower_AN_week_heat_fixedmmt,
                                 lower_AN_cold=lmort_temp$lower_AN_week_cold_fixedmmt,
                                 upper_AN=lmort_temp$upper_AN_week_fixedmmt,
                                 upper_AN_heat=lmort_temp$upper_AN_week_heat_fixedmmt,
                                 upper_AN_cold=lmort_temp$upper_AN_week_cold_fixedmmt,
                                 Ysim_mat=lmort_temp$Ysim_mat)

results_month = list()
results_month[[1]] = results_temp(data=data, inddf=inddf, indyears=indyears, 
                                  iyearmax=iyearmax, max_lag=max_lag, var_fun=var_fun,
                                  var_knots=var_knots, boundary_knots=boundary_knots,
                                  lag_knots=lag_knots, coef_expos=coef_expos, 
                                  vcov_expos=vcov_expos, MMTtrue=MMTtrue, at_temp=at_temp,
                                  cumpred= lmort_temp$cumpred_month_fixedmmt , 
                                  cumpred_low=lmort_temp$cumpred_low_month_fixedmmt,
                                  cumpred_high=lmort_temp$cumpred_high_month_fixedmmt,
                                  mmt=lmort_temp$mmt_month,
                                  AF=FALSE, AN = lmort_temp$AN_month_fixedmmt, 
                                  AN_heat=lmort_temp$AN_month_heat_fixedmmt,
                                  AN_cold=lmort_temp$AN_month_cold_fixedmmt,
                                  lower_AN=lmort_temp$lower_AN_month_fixedmmt,
                                  lower_AN_heat=lmort_temp$lower_AN_month_heat_fixedmmt,
                                  lower_AN_cold=lmort_temp$lower_AN_month_cold_fixedmmt,
                                  upper_AN=lmort_temp$upper_AN_month_fixedmmt,
                                  upper_AN_heat=lmort_temp$upper_AN_month_heat_fixedmmt,
                                  upper_AN_cold=lmort_temp$upper_AN_month_cold_fixedmmt,
                                  Ysim_mat=lmort_temp$Ysim_mat)

results_dow = list()
results_dow[[1]] = results_temp(data=data, inddf=inddf, indyears=indyears, 
                                iyearmax=iyearmax, max_lag=max_lag, var_fun=var_fun,
                                var_knots=var_knots, boundary_knots=boundary_knots,
                                lag_knots=lag_knots, coef_expos=coef_expos, 
                                vcov_expos=vcov_expos, MMTtrue=MMTtrue, at_temp=at_temp,
                                cumpred= lmort_temp$cumpred_dow_fixedmmt , 
                                cumpred_low=lmort_temp$cumpred_low_dow_fixedmmt,
                                cumpred_high=lmort_temp$cumpred_high_dow_fixedmmt,
                                mmt=lmort_temp$mmt_dow,
                                AF=FALSE, AN = lmort_temp$AN_dow_fixedmmt, 
                                AN_heat=lmort_temp$AN_dow_heat_fixedmmt,
                                AN_cold=lmort_temp$AN_dow_cold_fixedmmt,
                                lower_AN=lmort_temp$lower_AN_dow_fixedmmt,
                                lower_AN_heat=lmort_temp$lower_AN_dow_heat_fixedmmt,
                                lower_AN_cold=lmort_temp$lower_AN_dow_cold_fixedmmt,
                                upper_AN=lmort_temp$upper_AN_dow_fixedmmt,
                                upper_AN_heat=lmort_temp$upper_AN_dow_heat_fixedmmt,
                                upper_AN_cold=lmort_temp$upper_AN_dow_cold_fixedmmt,
                                Ysim_mat=lmort_temp$Ysim_mat)


### Hospitalizations - temperature
###############################

anal = "hosp_temp"
load(file=paste0("parameters_", anal, ".RData"))

iyearmax = 13

results_daily[[2]] <- results_temp(data=data, inddf=inddf, indyears=indyears, 
                                   iyearmax=iyearmax, max_lag=max_lag, var_fun=var_fun,
                                   var_knots=var_knots, boundary_knots=boundary_knots,
                                   lag_knots=lag_knots, coef_expos=coef_expos, 
                                   vcov_expos=vcov_expos, MMTtrue=MMTtrue, at_temp=at_temp,
                                   cumpred= lhosp_temp$cumpred_daily_fixedmmt , 
                                   cumpred_low=lhosp_temp$cumpred_low_daily_fixedmmt,
                                   cumpred_high=lhosp_temp$cumpred_high_daily_fixedmmt,
                                   mmt=lhosp_temp$mmt_daily,
                                   AF=FALSE, AN = lhosp_temp$AN_daily_fixedmmt, 
                                   AN_heat=lhosp_temp$AN_daily_heat_fixedmmt,
                                   AN_cold=lhosp_temp$AN_daily_cold_fixedmmt,
                                   lower_AN=lhosp_temp$lower_AN_daily_fixedmmt,
                                   lower_AN_heat=lhosp_temp$lower_AN_daily_heat_fixedmmt,
                                   lower_AN_cold=lhosp_temp$lower_AN_daily_cold_fixedmmt,
                                   upper_AN=lhosp_temp$upper_AN_daily_fixedmmt,
                                   upper_AN_heat=lhosp_temp$upper_AN_daily_heat_fixedmmt,
                                   upper_AN_cold=lhosp_temp$upper_AN_daily_cold_fixedmmt,
                                   Ysim_mat=lhosp_temp$Ysim_mat) 

results_week[[2]] = results_temp(data=data, inddf=inddf, indyears=indyears, 
                                 iyearmax=iyearmax, max_lag=max_lag, var_fun=var_fun,
                                 var_knots=var_knots, boundary_knots=boundary_knots,
                                 lag_knots=lag_knots, coef_expos=coef_expos, 
                                 vcov_expos=vcov_expos, MMTtrue=MMTtrue, at_temp=at_temp,
                                 cumpred= lhosp_temp$cumpred_week_fixedmmt , 
                                 cumpred_low=lhosp_temp$cumpred_low_week_fixedmmt,
                                 cumpred_high=lhosp_temp$cumpred_high_week_fixedmmt,
                                 mmt=lhosp_temp$mmt_week,
                                 AF=FALSE, AN = lhosp_temp$AN_week_fixedmmt, 
                                 AN_heat=lhosp_temp$AN_week_heat_fixedmmt,
                                 AN_cold=lhosp_temp$AN_week_cold_fixedmmt,
                                 lower_AN=lhosp_temp$lower_AN_week_fixedmmt,
                                 lower_AN_heat=lhosp_temp$lower_AN_week_heat_fixedmmt,
                                 lower_AN_cold=lhosp_temp$lower_AN_week_cold_fixedmmt,
                                 upper_AN=lhosp_temp$upper_AN_week_fixedmmt,
                                 upper_AN_heat=lhosp_temp$upper_AN_week_heat_fixedmmt,
                                 upper_AN_cold=lhosp_temp$upper_AN_week_cold_fixedmmt,
                                 Ysim_mat=lhosp_temp$Ysim_mat)

results_month[[2]] = results_temp(data=data, inddf=inddf, indyears=indyears, 
                                  iyearmax=iyearmax, max_lag=max_lag, var_fun=var_fun,
                                  var_knots=var_knots, boundary_knots=boundary_knots,
                                  lag_knots=lag_knots, coef_expos=coef_expos, 
                                  vcov_expos=vcov_expos, MMTtrue=MMTtrue, at_temp=at_temp,
                                  cumpred= lhosp_temp$cumpred_month_fixedmmt , 
                                  cumpred_low=lhosp_temp$cumpred_low_month_fixedmmt,
                                  cumpred_high=lhosp_temp$cumpred_high_month_fixedmmt,
                                  mmt=lhosp_temp$mmt_month,
                                  AF=FALSE, AN = lhosp_temp$AN_month_fixedmmt, 
                                  AN_heat=lhosp_temp$AN_month_heat_fixedmmt,
                                  AN_cold=lhosp_temp$AN_month_cold_fixedmmt,
                                  lower_AN=lhosp_temp$lower_AN_month_fixedmmt,
                                  lower_AN_heat=lhosp_temp$lower_AN_month_heat_fixedmmt,
                                  lower_AN_cold=lhosp_temp$lower_AN_month_cold_fixedmmt,
                                  upper_AN=lhosp_temp$upper_AN_month_fixedmmt,
                                  upper_AN_heat=lhosp_temp$upper_AN_month_heat_fixedmmt,
                                  upper_AN_cold=lhosp_temp$upper_AN_month_cold_fixedmmt,
                                  Ysim_mat=lhosp_temp$Ysim_mat)

results_dow[[2]] = results_temp(data=data, inddf=inddf, indyears=indyears, 
                                iyearmax=iyearmax, max_lag=max_lag, var_fun=var_fun,
                                var_knots=var_knots, boundary_knots=boundary_knots,
                                lag_knots=lag_knots, coef_expos=coef_expos, 
                                vcov_expos=vcov_expos, MMTtrue=MMTtrue, at_temp=at_temp,
                                cumpred= lhosp_temp$cumpred_dow_fixedmmt , 
                                cumpred_low=lhosp_temp$cumpred_low_dow_fixedmmt,
                                cumpred_high=lhosp_temp$cumpred_high_dow_fixedmmt,
                                mmt=lhosp_temp$mmt_dow,
                                AF=FALSE, AN = lhosp_temp$AN_dow_fixedmmt, 
                                AN_heat=lhosp_temp$AN_dow_heat_fixedmmt,
                                AN_cold=lhosp_temp$AN_dow_cold_fixedmmt,
                                lower_AN=lhosp_temp$lower_AN_dow_fixedmmt,
                                lower_AN_heat=lhosp_temp$lower_AN_dow_heat_fixedmmt,
                                lower_AN_cold=lhosp_temp$lower_AN_dow_cold_fixedmmt,
                                upper_AN=lhosp_temp$upper_AN_dow_fixedmmt,
                                upper_AN_heat=lhosp_temp$upper_AN_dow_heat_fixedmmt,
                                upper_AN_cold=lhosp_temp$upper_AN_dow_cold_fixedmmt,
                                Ysim_mat=lhosp_temp$Ysim_mat)



### Mortality - air pollution
###############################

anal = "mort_ap"
load(file=paste0("parameters_", anal, "_01_dow_holi.RData"))

iyearmax = 13

results_daily[[3]] = results_single(data = data, inddf = inddf, indyears=indyears,
                                    iyearmax = iyearmax, name_expos = "no2_lag01_20",
                                    ind_expos = ind_expos,
                                    coef_expos = coef_expos, vcov_expos = vcov_expos, at_no2 = at_no2,
                                    coefs = lmort_ap$coefs_daily, ses = lmort_ap$ses_daily, 
                                    AF = lmort_ap$AF_daily, lower_AF = lmort_ap$lower_AF_daily, 
                                    upper_AF = lmort_ap$upper_AF_daily,
                                    Ysim_mat = lmort_ap$Ysim_mat)


results_week[[3]] = results_single(data = data, inddf = inddf, indyears=indyears,
                                   iyearmax = iyearmax, name_expos = "no2_lag01_20",
                                   ind_expos = ind_expos,
                                   coef_expos = coef_expos, vcov_expos = vcov_expos, at_no2 = at_no2,
                                   coefs = lmort_ap$coefs_week, ses = lmort_ap$ses_week, 
                                   AF = lmort_ap$AF_week, lower_AF = lmort_ap$lower_AF_week, 
                                   upper_AF = lmort_ap$upper_AF_week,
                                   Ysim_mat = lmort_ap$Ysim_mat)

results_month[[3]] = results_single(data = data, inddf = inddf, indyears=indyears,
                                    iyearmax = iyearmax, name_expos = "no2_lag01_20",
                                    ind_expos = ind_expos,
                                    coef_expos = coef_expos, vcov_expos = vcov_expos, at_no2 = at_no2,
                                    coefs = lmort_ap$coefs_month, ses = lmort_ap$ses_month, 
                                    AF = lmort_ap$AF_month, lower_AF = lmort_ap$lower_AF_month, 
                                    upper_AF = lmort_ap$upper_AF_month,
                                    Ysim_mat = lmort_ap$Ysim_mat)

results_dow[[3]] = results_single(data = data, inddf = inddf, indyears=indyears,
                                  iyearmax = iyearmax, name_expos = "no2_lag01_20",
                                  ind_expos = ind_expos,
                                  coef_expos = coef_expos, vcov_expos = vcov_expos, at_no2 = at_no2,
                                  coefs = lmort_ap$coefs_dow, ses = lmort_ap$ses_dow, 
                                  AF = lmort_ap$AF_dow, lower_AF = lmort_ap$lower_AF_dow, 
                                  upper_AF = lmort_ap$upper_AF_dow,
                                  Ysim_mat = lmort_ap$Ysim_mat)



### Hospitalizations - air pollution
######################################

anal = "hosp_ap"
load(file=paste0("parameters_", anal, "_01_dow_holi.RData"))

iyearmax = 13

results_daily[[4]] = results_single(data = data, inddf = inddf, indyears=indyears,
                                    iyearmax = iyearmax, name_expos = "no2_lag06_20",
                                    ind_expos = ind_expos,
                                    coef_expos = coef_expos, vcov_expos = vcov_expos, at_no2 = at_no2,
                                    coefs = lhosp_ap$coefs_daily, ses = lhosp_ap$ses_daily, 
                                    AF = lhosp_ap$AF_daily, lower_AF = lhosp_ap$lower_AF_daily, 
                                    upper_AF = lhosp_ap$upper_AF_daily,
                                    Ysim_mat = lhosp_ap$Ysim_mat)


results_week[[4]] = results_single(data = data, inddf = inddf, indyears=indyears,
                                   iyearmax = iyearmax, name_expos = "no2_lag06_20",
                                   ind_expos = ind_expos,
                                   coef_expos = coef_expos, vcov_expos = vcov_expos, at_no2 = at_no2,
                                   coefs = lhosp_ap$coefs_week, ses = lhosp_ap$ses_week, 
                                   AF = lhosp_ap$AF_week, lower_AF = lhosp_ap$lower_AF_week, 
                                   upper_AF = lhosp_ap$upper_AF_week,
                                   Ysim_mat = lhosp_ap$Ysim_mat)

results_month[[4]] = results_single(data = data, inddf = inddf, indyears=indyears,
                                    iyearmax = iyearmax, name_expos = "no2_lag06_20",
                                    ind_expos = ind_expos,
                                    coef_expos = coef_expos, vcov_expos = vcov_expos, at_no2 = at_no2,
                                    coefs = lhosp_ap$coefs_month, ses = lhosp_ap$ses_month, 
                                    AF = lhosp_ap$AF_month, lower_AF = lhosp_ap$lower_AF_month, 
                                    upper_AF = lhosp_ap$upper_AF_month,
                                    Ysim_mat = lhosp_ap$Ysim_mat)

results_dow[[4]] = results_single(data = data, inddf = inddf, indyears=indyears,
                                  iyearmax = iyearmax, name_expos = "no2_lag06_20",
                                  ind_expos = ind_expos,
                                  coef_expos = coef_expos, vcov_expos = vcov_expos, at_no2 = at_no2,
                                  coefs = lhosp_ap$coefs_dow, ses = lhosp_ap$ses_dow, 
                                  AF = lhosp_ap$AF_dow, lower_AF = lhosp_ap$lower_AF_dow, 
                                  upper_AF = lhosp_ap$upper_AF_dow,
                                  Ysim_mat = lhosp_ap$Ysim_mat)




##################################################
##################################################
##################################################
##################################################
# Case where the number of cases is divided by 10 #
##################################################
##################################################
##################################################
##################################################

indyears <- c(1:10, seq(15, 25, 5))

########################################################
# Append individual results files into a single list
########################################################

## Temperature results
##########################

## Mortality 

anal = "mort_temp"

k = 1
#for (i in 1:length(indyears)) {
for (i in 1:13) {
  load(file=paste0("results_sim_bcn10_v14_", anal, "_", indyears[i], ".RData"))
  if (i==1) {
    lmort_temp = rlist
    names(lmort_temp) = c("pred_prc", "at_temp", "nsim", "pos_p99", "pos_p01", 
                          "Ysim_mat", "converged_week", "converged_month", "converged_dow",
                          "coefs_daily", "coefs_week", "coefs_month", "coefs_dow",
                          "disp_daily", "disp_week", "disp_month", "disp_dow",
                          "Ypred_daily", "Ypred_week", "Ypred_month", "Ypred_dow", 
                          "mmt_daily", "mmt_week", "mmt_month", "mmt_dow",
                          "cumpred_daily", "cumpred_week", "cumpred_month", "cumpred_dow",
                          "cumpred_low_daily", "cumpred_low_week", "cumpred_low_month", "cumpred_low_dow",
                          "cumpred_high_daily", "cumpred_high_week", "cumpred_high_month", "cumpred_high_dow",
                          "cumpred_daily_fixedmmt", "cumpred_week_fixedmmt", "cumpred_month_fixedmmt", "cumpred_dow_fixedmmt", 
                          "cumpred_low_daily_fixedmmt", "cumpred_low_week_fixedmmt", "cumpred_low_month_fixedmmt", "cumpred_low_dow_fixedmmt",
                          "cumpred_high_daily_fixedmmt", "cumpred_high_week_fixedmmt", "cumpred_high_month_fixedmmt", "cumpred_high_dow_fixedmmt",
                          "RRp99lag_daily_fixedmmt", "RRp99lag_week_fixedmmt", "RRp99lag_month_fixedmmt", "RRp99lag_dow_fixedmmt",
                          "RRp99lag_low_daily_fixedmmt", "RRp99lag_low_week_fixedmmt", "RRp99lag_low_month_fixedmmt", "RRp99lag_low_dow_fixedmmt",
                          "RRp99lag_high_daily_fixedmmt", "RRp99lag_high_week_fixedmmt", "RRp99lag_high_month_fixedmmt", "RRp99lag_high_dow_fixedmmt", 
                          "RRp01lag_daily_fixedmmt", "RRp01lag_week_fixedmmt", "RRp01lag_month_fixedmmt", "RRp01lag_dow_fixedmmt",
                          "RRp01lag_low_daily_fixedmmt", "RRp01lag_low_week_fixedmmt", "RRp01lag_low_month_fixedmmt", "RRp01lag_low_dow_fixedmmt",
                          "RRp01lag_high_daily_fixedmmt", "RRp01lag_high_week_fixedmmt", "RRp01lag_high_month_fixedmmt", "RRp01lag_high_dow_fixedmmt",
                          "AN_daily_fixedmmt", "AN_week_fixedmmt", "AN_month_fixedmmt", "AN_dow_fixedmmt",
                          "AN_daily_cold_fixedmmt", "AN_week_cold_fixedmmt", "AN_month_cold_fixedmmt", "AN_dow_cold_fixedmmt", 
                          "AN_daily_heat_fixedmmt", "AN_week_heat_fixedmmt", "AN_month_heat_fixedmmt", "AN_dow_heat_fixedmmt",
                          "lower_AN_daily_fixedmmt", "lower_AN_week_fixedmmt", "lower_AN_month_fixedmmt", "lower_AN_dow_fixedmmt",
                          "lower_AN_daily_cold_fixedmmt", "lower_AN_week_cold_fixedmmt", "lower_AN_month_cold_fixedmmt", "lower_AN_dow_cold_fixedmmt", 
                          "lower_AN_daily_heat_fixedmmt", "lower_AN_week_heat_fixedmmt", "lower_AN_month_heat_fixedmmt", "lower_AN_dow_heat_fixedmmt",
                          "upper_AN_daily_fixedmmt", "upper_AN_week_fixedmmt", "upper_AN_month_fixedmmt", "upper_AN_dow_fixedmmt",
                          "upper_AN_daily_cold_fixedmmt", "upper_AN_week_cold_fixedmmt", "upper_AN_month_cold_fixedmmt", "upper_AN_dow_cold_fixedmmt",
                          "upper_AN_daily_heat_fixedmmt", "upper_AN_week_heat_fixedmmt", "upper_AN_month_heat_fixedmmt", "upper_AN_dow_heat_fixedmmt")
    
    for (j in 2:length(lmort_temp)) {
      lmort_temp[[j]] = list(lmort_temp[[j]])
    }
  } else {
    lmort_temp[[1]] = c(lmort_temp[[1]], rlist[[1]])
    for (j in 2:length(lmort_temp)) {
      lmort_temp[[j]][[k]] = rlist[[j]] 
    }
  }
  k = k + 1
}

rm(rlist)



## Hospitalizations 

anal = "hosp_temp"

k = 1
#for (i in 1:length(indyears)) {
for (i in 1:13) {
  load(file=paste0("results_sim_bcn10_v14_", anal, "_", indyears[i], ".RData"))
  if (i==1) {
    lhosp_temp = rlist
    names(lhosp_temp) = c("pred_prc", "at_temp", "nsim", "pos_p99", "pos_p01", 
                          "Ysim_mat", "converged_week", "converged_month", "converged_dow",
                          "coefs_daily", "coefs_week", "coefs_month", "coefs_dow",
                          "disp_daily", "disp_week", "disp_month", "disp_dow",
                          "Ypred_daily", "Ypred_week", "Ypred_month", "Ypred_dow", 
                          "mmt_daily", "mmt_week", "mmt_month", "mmt_dow",
                          "cumpred_daily", "cumpred_week", "cumpred_month", "cumpred_dow",
                          "cumpred_low_daily", "cumpred_low_week", "cumpred_low_month", "cumpred_low_dow",
                          "cumpred_high_daily", "cumpred_high_week", "cumpred_high_month", "cumpred_high_dow",
                          "cumpred_daily_fixedmmt", "cumpred_week_fixedmmt", "cumpred_month_fixedmmt", "cumpred_dow_fixedmmt", 
                          "cumpred_low_daily_fixedmmt", "cumpred_low_week_fixedmmt", "cumpred_low_month_fixedmmt", "cumpred_low_dow_fixedmmt",
                          "cumpred_high_daily_fixedmmt", "cumpred_high_week_fixedmmt", "cumpred_high_month_fixedmmt", "cumpred_high_dow_fixedmmt",
                          "RRp99lag_daily_fixedmmt", "RRp99lag_week_fixedmmt", "RRp99lag_month_fixedmmt", "RRp99lag_dow_fixedmmt",
                          "RRp99lag_low_daily_fixedmmt", "RRp99lag_low_week_fixedmmt", "RRp99lag_low_month_fixedmmt", "RRp99lag_low_dow_fixedmmt",
                          "RRp99lag_high_daily_fixedmmt", "RRp99lag_high_week_fixedmmt", "RRp99lag_high_month_fixedmmt", "RRp99lag_high_dow_fixedmmt", 
                          "RRp01lag_daily_fixedmmt", "RRp01lag_week_fixedmmt", "RRp01lag_month_fixedmmt", "RRp01lag_dow_fixedmmt",
                          "RRp01lag_low_daily_fixedmmt", "RRp01lag_low_week_fixedmmt", "RRp01lag_low_month_fixedmmt", "RRp01lag_low_dow_fixedmmt",
                          "RRp01lag_high_daily_fixedmmt", "RRp01lag_high_week_fixedmmt", "RRp01lag_high_month_fixedmmt", "RRp01lag_high_dow_fixedmmt",
                          "AN_daily_fixedmmt", "AN_week_fixedmmt", "AN_month_fixedmmt", "AN_dow_fixedmmt",
                          "AN_daily_cold_fixedmmt", "AN_week_cold_fixedmmt", "AN_month_cold_fixedmmt", "AN_dow_cold_fixedmmt", 
                          "AN_daily_heat_fixedmmt", "AN_week_heat_fixedmmt", "AN_month_heat_fixedmmt", "AN_dow_heat_fixedmmt",
                          "lower_AN_daily_fixedmmt", "lower_AN_week_fixedmmt", "lower_AN_month_fixedmmt", "lower_AN_dow_fixedmmt",
                          "lower_AN_daily_cold_fixedmmt", "lower_AN_week_cold_fixedmmt", "lower_AN_month_cold_fixedmmt", "lower_AN_dow_cold_fixedmmt", 
                          "lower_AN_daily_heat_fixedmmt", "lower_AN_week_heat_fixedmmt", "lower_AN_month_heat_fixedmmt", "lower_AN_dow_heat_fixedmmt",
                          "upper_AN_daily_fixedmmt", "upper_AN_week_fixedmmt", "upper_AN_month_fixedmmt", "upper_AN_dow_fixedmmt",
                          "upper_AN_daily_cold_fixedmmt", "upper_AN_week_cold_fixedmmt", "upper_AN_month_cold_fixedmmt", "upper_AN_dow_cold_fixedmmt",
                          "upper_AN_daily_heat_fixedmmt", "upper_AN_week_heat_fixedmmt", "upper_AN_month_heat_fixedmmt", "upper_AN_dow_heat_fixedmmt")
    for (j in 2:length(lhosp_temp)) {
      lhosp_temp[[j]] = list(lhosp_temp[[j]])
    }
  } else {
    lhosp_temp[[1]] = c(lhosp_temp[[1]], rlist[[1]])
    for (j in 2:length(lhosp_temp)) {
      lhosp_temp[[j]][[k]] = rlist[[j]] 
    }
  }
  k = k + 1
}

rm(rlist)


## Air pollution results
##########################

## Mortality 

anal = "mort_ap"

k = 1
for (i in 1:13) {
  load(file=paste0("results_sim_bcn10_v14_", anal, "_dow_holi_0101_", indyears[i], ".RData"))
  if (i==1) {
    lmort_ap = rlist
    names(lmort_ap) = c("nsim", "Ysim_mat", 
                        "converged_week", "converged_month", "converged_dow",
                        "coefs_daily", "coefs_week", "coefs_month", "coefs_dow",
                        "ses_daily", "ses_week", "ses_month", "ses_dow",       
                        "Ypred_daily", "Ypred_week", "Ypred_month", "Ypred_dow",  
                        "AF_daily", "AF_week", "AF_month", "AF_dow", 
                        "lower_AF_daily", "lower_AF_week", "lower_AF_month", "lower_AF_dow", 
                        "upper_AF_daily", "upper_AF_week", "upper_AF_month", "upper_AF_dow", 
                        "disp_daily", "disp_week", "disp_month", "disp_dow")
    for (j in 2:length(lmort_ap)) {
      lmort_ap[[j]] = list(lmort_ap[[j]])
    }
  } else {
    lmort_ap[[1]] = c(lmort_ap[[1]], rlist[[1]])
    for (j in 2:length(lmort_ap)) {
      lmort_ap[[j]][[k]] = rlist[[j]] 
    }
  }
  k = k + 1
}

rm(rlist)



## Hospitalizations 

anal = "hosp_ap"

k = 1
#for (i in 1:length(indyears)) {
for (i in 1:13) {
  load(file=paste0("results_sim_bcn10_v14_", anal, "_dow_holi_0101_1pct_", indyears[i], ".RData"))
  if (i==1) {
    lhosp_ap = rlist
    names(lhosp_ap) = c("nsim", "Ysim_mat", 
                        "converged_week", "converged_month", "converged_dow",
                        "coefs_daily", "coefs_week", "coefs_month", "coefs_dow",
                        "ses_daily", "ses_week", "ses_month", "ses_dow",       
                        "Ypred_daily", "Ypred_week", "Ypred_month", "Ypred_dow",  
                        "AF_daily", "AF_week", "AF_month", "AF_dow", 
                        "lower_AF_daily", "lower_AF_week", "lower_AF_month", "lower_AF_dow", 
                        "upper_AF_daily", "upper_AF_week", "upper_AF_month", "upper_AF_dow", 
                        "disp_daily", "disp_week", "disp_month", "disp_dow")
    for (j in 2:length(lhosp_ap)) {
      lhosp_ap[[j]] = list(lhosp_ap[[j]])
    }
  } else {
    lhosp_ap[[1]] = c(lhosp_ap[[1]], rlist[[1]])
    for (j in 2:length(lhosp_ap)) {
      lhosp_ap[[j]][[k]] = rlist[[j]] 
    }
  }
  k = k + 1
}

rm(rlist)


###############################################
# Load original data and model specifications
###############################################

load(file = "bcn_ts_exposures.RData")

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

at_no2 <- seq(0, max(data$no2_lag01, na.rm=T), by=0.1)

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

###########################################################
# Obtain results from data: bias, RMSE, coverage, power, 
###########################################################

nsim <- 500

results_daily10 = list()

### Mortality - temperature
###############################

iyearmax = 13

anal = "mort_temp"
load(file=paste0("parameters_", anal, ".RData"))


results_daily10[[1]] <- results_temp(data=data, inddf=inddf, indyears=indyears, 
                                   iyearmax=iyearmax, max_lag=max_lag, var_fun=var_fun,
                                   var_knots=var_knots, boundary_knots=boundary_knots,
                                   lag_knots=lag_knots, coef_expos=coef_expos, 
                                   vcov_expos=vcov_expos, MMTtrue=MMTtrue, at_temp=at_temp,
                                   cumpred= lmort_temp$cumpred_daily_fixedmmt , 
                                   cumpred_low=lmort_temp$cumpred_low_daily_fixedmmt,
                                   cumpred_high=lmort_temp$cumpred_high_daily_fixedmmt,
                                   mmt=lmort_temp$mmt_daily,
                                   AF=FALSE, AN = lmort_temp$AN_daily_fixedmmt, 
                                   AN_heat=lmort_temp$AN_daily_heat_fixedmmt,
                                   AN_cold=lmort_temp$AN_daily_cold_fixedmmt,
                                   lower_AN=lmort_temp$lower_AN_daily_fixedmmt,
                                   lower_AN_heat=lmort_temp$lower_AN_daily_heat_fixedmmt,
                                   lower_AN_cold=lmort_temp$lower_AN_daily_cold_fixedmmt,
                                   upper_AN=lmort_temp$upper_AN_daily_fixedmmt,
                                   upper_AN_heat=lmort_temp$upper_AN_daily_heat_fixedmmt,
                                   upper_AN_cold=lmort_temp$upper_AN_daily_cold_fixedmmt,
                                   Ysim_mat=lmort_temp$Ysim_mat) 

results_week10 = list()
results_week10[[1]] = results_temp(data=data, inddf=inddf, indyears=indyears, 
                                 iyearmax=iyearmax, max_lag=max_lag, var_fun=var_fun,
                                 var_knots=var_knots, boundary_knots=boundary_knots,
                                 lag_knots=lag_knots, coef_expos=coef_expos, 
                                 vcov_expos=vcov_expos, MMTtrue=MMTtrue, at_temp=at_temp,
                                 cumpred= lmort_temp$cumpred_week_fixedmmt , 
                                 cumpred_low=lmort_temp$cumpred_low_week_fixedmmt,
                                 cumpred_high=lmort_temp$cumpred_high_week_fixedmmt,
                                 mmt=lmort_temp$mmt_week,
                                 AF=FALSE, AN = lmort_temp$AN_week_fixedmmt, 
                                 AN_heat=lmort_temp$AN_week_heat_fixedmmt,
                                 AN_cold=lmort_temp$AN_week_cold_fixedmmt,
                                 lower_AN=lmort_temp$lower_AN_week_fixedmmt,
                                 lower_AN_heat=lmort_temp$lower_AN_week_heat_fixedmmt,
                                 lower_AN_cold=lmort_temp$lower_AN_week_cold_fixedmmt,
                                 upper_AN=lmort_temp$upper_AN_week_fixedmmt,
                                 upper_AN_heat=lmort_temp$upper_AN_week_heat_fixedmmt,
                                 upper_AN_cold=lmort_temp$upper_AN_week_cold_fixedmmt,
                                 Ysim_mat=lmort_temp$Ysim_mat)

results_month10 = list()
results_month10[[1]] = results_temp(data=data, inddf=inddf, indyears=indyears, 
                                  iyearmax=iyearmax, max_lag=max_lag, var_fun=var_fun,
                                  var_knots=var_knots, boundary_knots=boundary_knots,
                                  lag_knots=lag_knots, coef_expos=coef_expos, 
                                  vcov_expos=vcov_expos, MMTtrue=MMTtrue, at_temp=at_temp,
                                  cumpred= lmort_temp$cumpred_month_fixedmmt , 
                                  cumpred_low=lmort_temp$cumpred_low_month_fixedmmt,
                                  cumpred_high=lmort_temp$cumpred_high_month_fixedmmt,
                                  mmt=lmort_temp$mmt_month,
                                  AF=FALSE, AN = lmort_temp$AN_month_fixedmmt, 
                                  AN_heat=lmort_temp$AN_month_heat_fixedmmt,
                                  AN_cold=lmort_temp$AN_month_cold_fixedmmt,
                                  lower_AN=lmort_temp$lower_AN_month_fixedmmt,
                                  lower_AN_heat=lmort_temp$lower_AN_month_heat_fixedmmt,
                                  lower_AN_cold=lmort_temp$lower_AN_month_cold_fixedmmt,
                                  upper_AN=lmort_temp$upper_AN_month_fixedmmt,
                                  upper_AN_heat=lmort_temp$upper_AN_month_heat_fixedmmt,
                                  upper_AN_cold=lmort_temp$upper_AN_month_cold_fixedmmt,
                                  Ysim_mat=lmort_temp$Ysim_mat)

results_dow10 = list()
results_dow10[[1]] = results_temp(data=data, inddf=inddf, indyears=indyears, 
                                iyearmax=iyearmax, max_lag=max_lag, var_fun=var_fun,
                                var_knots=var_knots, boundary_knots=boundary_knots,
                                lag_knots=lag_knots, coef_expos=coef_expos, 
                                vcov_expos=vcov_expos, MMTtrue=MMTtrue, at_temp=at_temp,
                                cumpred= lmort_temp$cumpred_dow_fixedmmt , 
                                cumpred_low=lmort_temp$cumpred_low_dow_fixedmmt,
                                cumpred_high=lmort_temp$cumpred_high_dow_fixedmmt,
                                mmt=lmort_temp$mmt_dow,
                                AF=FALSE, AN = lmort_temp$AN_dow_fixedmmt, 
                                AN_heat=lmort_temp$AN_dow_heat_fixedmmt,
                                AN_cold=lmort_temp$AN_dow_cold_fixedmmt,
                                lower_AN=lmort_temp$lower_AN_dow_fixedmmt,
                                lower_AN_heat=lmort_temp$lower_AN_dow_heat_fixedmmt,
                                lower_AN_cold=lmort_temp$lower_AN_dow_cold_fixedmmt,
                                upper_AN=lmort_temp$upper_AN_dow_fixedmmt,
                                upper_AN_heat=lmort_temp$upper_AN_dow_heat_fixedmmt,
                                upper_AN_cold=lmort_temp$upper_AN_dow_cold_fixedmmt,
                                Ysim_mat=lmort_temp$Ysim_mat)


### Hospitalizations - temperature
###############################

anal = "hosp_temp"
load(file=paste0("parameters_", anal, ".RData"))

iyearmax = 13

results_daily10[[2]] <- results_temp(data=data, inddf=inddf, indyears=indyears, 
                                   iyearmax=iyearmax, max_lag=max_lag, var_fun=var_fun,
                                   var_knots=var_knots, boundary_knots=boundary_knots,
                                   lag_knots=lag_knots, coef_expos=coef_expos, 
                                   vcov_expos=vcov_expos, MMTtrue=MMTtrue, at_temp=at_temp,
                                   cumpred= lhosp_temp$cumpred_daily_fixedmmt , 
                                   cumpred_low=lhosp_temp$cumpred_low_daily_fixedmmt,
                                   cumpred_high=lhosp_temp$cumpred_high_daily_fixedmmt,
                                   mmt=lhosp_temp$mmt_daily,
                                   AF=FALSE, AN = lhosp_temp$AN_daily_fixedmmt, 
                                   AN_heat=lhosp_temp$AN_daily_heat_fixedmmt,
                                   AN_cold=lhosp_temp$AN_daily_cold_fixedmmt,
                                   lower_AN=lhosp_temp$lower_AN_daily_fixedmmt,
                                   lower_AN_heat=lhosp_temp$lower_AN_daily_heat_fixedmmt,
                                   lower_AN_cold=lhosp_temp$lower_AN_daily_cold_fixedmmt,
                                   upper_AN=lhosp_temp$upper_AN_daily_fixedmmt,
                                   upper_AN_heat=lhosp_temp$upper_AN_daily_heat_fixedmmt,
                                   upper_AN_cold=lhosp_temp$upper_AN_daily_cold_fixedmmt,
                                   Ysim_mat=lhosp_temp$Ysim_mat) 

results_week10[[2]] = results_temp(data=data, inddf=inddf, indyears=indyears, 
                                 iyearmax=iyearmax, max_lag=max_lag, var_fun=var_fun,
                                 var_knots=var_knots, boundary_knots=boundary_knots,
                                 lag_knots=lag_knots, coef_expos=coef_expos, 
                                 vcov_expos=vcov_expos, MMTtrue=MMTtrue, at_temp=at_temp,
                                 cumpred= lhosp_temp$cumpred_week_fixedmmt , 
                                 cumpred_low=lhosp_temp$cumpred_low_week_fixedmmt,
                                 cumpred_high=lhosp_temp$cumpred_high_week_fixedmmt,
                                 mmt=lhosp_temp$mmt_week,
                                 AF=FALSE, AN = lhosp_temp$AN_week_fixedmmt, 
                                 AN_heat=lhosp_temp$AN_week_heat_fixedmmt,
                                 AN_cold=lhosp_temp$AN_week_cold_fixedmmt,
                                 lower_AN=lhosp_temp$lower_AN_week_fixedmmt,
                                 lower_AN_heat=lhosp_temp$lower_AN_week_heat_fixedmmt,
                                 lower_AN_cold=lhosp_temp$lower_AN_week_cold_fixedmmt,
                                 upper_AN=lhosp_temp$upper_AN_week_fixedmmt,
                                 upper_AN_heat=lhosp_temp$upper_AN_week_heat_fixedmmt,
                                 upper_AN_cold=lhosp_temp$upper_AN_week_cold_fixedmmt,
                                 Ysim_mat=lhosp_temp$Ysim_mat)

results_month10[[2]] = results_temp(data=data, inddf=inddf, indyears=indyears, 
                                  iyearmax=iyearmax, max_lag=max_lag, var_fun=var_fun,
                                  var_knots=var_knots, boundary_knots=boundary_knots,
                                  lag_knots=lag_knots, coef_expos=coef_expos, 
                                  vcov_expos=vcov_expos, MMTtrue=MMTtrue, at_temp=at_temp,
                                  cumpred= lhosp_temp$cumpred_month_fixedmmt , 
                                  cumpred_low=lhosp_temp$cumpred_low_month_fixedmmt,
                                  cumpred_high=lhosp_temp$cumpred_high_month_fixedmmt,
                                  mmt=lhosp_temp$mmt_month,
                                  AF=FALSE, AN = lhosp_temp$AN_month_fixedmmt, 
                                  AN_heat=lhosp_temp$AN_month_heat_fixedmmt,
                                  AN_cold=lhosp_temp$AN_month_cold_fixedmmt,
                                  lower_AN=lhosp_temp$lower_AN_month_fixedmmt,
                                  lower_AN_heat=lhosp_temp$lower_AN_month_heat_fixedmmt,
                                  lower_AN_cold=lhosp_temp$lower_AN_month_cold_fixedmmt,
                                  upper_AN=lhosp_temp$upper_AN_month_fixedmmt,
                                  upper_AN_heat=lhosp_temp$upper_AN_month_heat_fixedmmt,
                                  upper_AN_cold=lhosp_temp$upper_AN_month_cold_fixedmmt,
                                  Ysim_mat=lhosp_temp$Ysim_mat)

results_dow10[[2]] = results_temp(data=data, inddf=inddf, indyears=indyears, 
                                iyearmax=iyearmax, max_lag=max_lag, var_fun=var_fun,
                                var_knots=var_knots, boundary_knots=boundary_knots,
                                lag_knots=lag_knots, coef_expos=coef_expos, 
                                vcov_expos=vcov_expos, MMTtrue=MMTtrue, at_temp=at_temp,
                                cumpred= lhosp_temp$cumpred_dow_fixedmmt , 
                                cumpred_low=lhosp_temp$cumpred_low_dow_fixedmmt,
                                cumpred_high=lhosp_temp$cumpred_high_dow_fixedmmt,
                                mmt=lhosp_temp$mmt_dow,
                                AF=FALSE, AN = lhosp_temp$AN_dow_fixedmmt, 
                                AN_heat=lhosp_temp$AN_dow_heat_fixedmmt,
                                AN_cold=lhosp_temp$AN_dow_cold_fixedmmt,
                                lower_AN=lhosp_temp$lower_AN_dow_fixedmmt,
                                lower_AN_heat=lhosp_temp$lower_AN_dow_heat_fixedmmt,
                                lower_AN_cold=lhosp_temp$lower_AN_dow_cold_fixedmmt,
                                upper_AN=lhosp_temp$upper_AN_dow_fixedmmt,
                                upper_AN_heat=lhosp_temp$upper_AN_dow_heat_fixedmmt,
                                upper_AN_cold=lhosp_temp$upper_AN_dow_cold_fixedmmt,
                                Ysim_mat=lhosp_temp$Ysim_mat)



### Mortality - air pollution
###############################

anal = "mort_ap"
load(file=paste0("parameters_", anal, "_01_dow_holi.RData"))

iyearmax = 13

results_daily10[[3]] = results_single(data = data, inddf = inddf, indyears=indyears,
                                      iyearmax = iyearmax, name_expos = "no2_lag01_20",
                                      ind_expos = ind_expos,
                                      coef_expos = coef_expos, vcov_expos = vcov_expos, at_no2 = at_no2,
                                      coefs = lmort_ap$coefs_daily, ses = lmort_ap$ses_daily, 
                                      AF = lmort_ap$AF_daily, lower_AF = lmort_ap$lower_AF_daily, 
                                      upper_AF = lmort_ap$upper_AF_daily,
                                      Ysim_mat = lmort_ap$Ysim_mat)


results_week10[[3]] = results_single(data = data, inddf = inddf, indyears=indyears,
                                     iyearmax = iyearmax, name_expos = "no2_lag01_20",
                                     ind_expos = ind_expos,
                                     coef_expos = coef_expos, vcov_expos = vcov_expos, at_no2 = at_no2,
                                     coefs = lmort_ap$coefs_week, ses = lmort_ap$ses_week, 
                                     AF = lmort_ap$AF_week, lower_AF = lmort_ap$lower_AF_week, 
                                     upper_AF = lmort_ap$upper_AF_week,
                                     Ysim_mat = lmort_ap$Ysim_mat)

results_month10[[3]] = results_single(data = data, inddf = inddf, indyears=indyears,
                                      iyearmax = iyearmax, name_expos = "no2_lag01_20",
                                      ind_expos = ind_expos,
                                      coef_expos = coef_expos, vcov_expos = vcov_expos, at_no2 = at_no2,
                                      coefs = lmort_ap$coefs_month, ses = lmort_ap$ses_month, 
                                      AF = lmort_ap$AF_month, lower_AF = lmort_ap$lower_AF_month, 
                                      upper_AF = lmort_ap$upper_AF_month,
                                      Ysim_mat = lmort_ap$Ysim_mat)

results_dow10[[3]] = results_single(data = data, inddf = inddf, indyears=indyears,
                                    iyearmax = iyearmax, name_expos = "no2_lag01_20",
                                    ind_expos = ind_expos,
                                    coef_expos = coef_expos, vcov_expos = vcov_expos, at_no2 = at_no2,
                                    coefs = lmort_ap$coefs_dow, ses = lmort_ap$ses_dow, 
                                    AF = lmort_ap$AF_dow, lower_AF = lmort_ap$lower_AF_dow, 
                                    upper_AF = lmort_ap$upper_AF_dow,
                                    Ysim_mat = lmort_ap$Ysim_mat)



### Hospitalizations - air pollution
######################################

anal = "hosp_ap"
load(file=paste0("parameters_", anal, "_01_dow_holi.RData"))

iyearmax = 13

results_daily10[[4]] = results_single(data = data, inddf = inddf, indyears=indyears,
                                      iyearmax = iyearmax, name_expos = "no2_lag06_20",
                                      ind_expos = ind_expos,
                                      coef_expos = coef_expos, vcov_expos = vcov_expos, at_no2 = at_no2,
                                      coefs = lhosp_ap$coefs_daily, ses = lhosp_ap$ses_daily, 
                                      AF = lhosp_ap$AF_daily, lower_AF = lhosp_ap$lower_AF_daily, 
                                      upper_AF = lhosp_ap$upper_AF_daily,
                                      Ysim_mat = lhosp_ap$Ysim_mat)


results_week10[[4]] = results_single(data = data, inddf = inddf, indyears=indyears,
                                     iyearmax = iyearmax, name_expos = "no2_lag06_20",
                                     ind_expos = ind_expos,
                                     coef_expos = coef_expos, vcov_expos = vcov_expos, at_no2 = at_no2,
                                     coefs = lhosp_ap$coefs_week, ses = lhosp_ap$ses_week, 
                                     AF = lhosp_ap$AF_week, lower_AF = lhosp_ap$lower_AF_week, 
                                     upper_AF = lhosp_ap$upper_AF_week,
                                     Ysim_mat = lhosp_ap$Ysim_mat)

results_month10[[4]] = results_single(data = data, inddf = inddf, indyears=indyears,
                                      iyearmax = iyearmax, name_expos = "no2_lag06_20",
                                      ind_expos = ind_expos,
                                      coef_expos = coef_expos, vcov_expos = vcov_expos, at_no2 = at_no2,
                                      coefs = lhosp_ap$coefs_month, ses = lhosp_ap$ses_month, 
                                      AF = lhosp_ap$AF_month, lower_AF = lhosp_ap$lower_AF_month, 
                                      upper_AF = lhosp_ap$upper_AF_month,
                                      Ysim_mat = lhosp_ap$Ysim_mat)

results_dow10[[4]] = results_single(data = data, inddf = inddf, indyears=indyears,
                                    iyearmax = iyearmax, name_expos = "no2_lag06_20",
                                    ind_expos = ind_expos,
                                    coef_expos = coef_expos, vcov_expos = vcov_expos, at_no2 = at_no2,
                                    coefs = lhosp_ap$coefs_dow, ses = lhosp_ap$ses_dow, 
                                    AF = lhosp_ap$AF_dow, lower_AF = lhosp_ap$lower_AF_dow, 
                                    upper_AF = lhosp_ap$upper_AF_dow,
                                    Ysim_mat = lhosp_ap$Ysim_mat)


