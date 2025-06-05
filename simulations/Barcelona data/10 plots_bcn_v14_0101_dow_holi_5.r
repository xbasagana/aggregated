# Using new temp results

#########
#########
# Plots #
#########
#########

##########################################################
# Figure bias-mse-coverage for exposure response function 
##########################################################

# Excluding nyears = 1
##########################

# Add a column with power
##########################

aa = c(results_daily[[1]]$lower_bias_cum[-1], results_daily[[1]]$upper_bias_cum[-1],
       results_week[[1]]$lower_bias_cum[-1], results_week[[1]]$upper_bias_cum[-1],
       results_month[[1]]$lower_bias_cum[-1], results_month[[1]]$upper_bias_cum[-1],
       results_dow[[1]]$lower_bias_cum[-1], results_dow[[1]]$upper_bias_cum[-1])
ylim_a = c( min(aa, na.rm=T), max(aa, na.rm=T))

aa = c(results_daily[[1]]$rmse_cum[-1], results_week[[1]]$rmse_cum[-1], 
       results_month[[1]]$rmse_cum[-1], results_dow[[1]]$rmse_cum[-1])
ylim_b = c( 0, max(aa, na.rm=T))

aa = c(results_daily[[1]]$lower_coverage_cum[-1], results_week[[1]]$lower_coverage_cum[-1], 
       results_month[[1]]$lower_coverage_cum[-1], results_dow[[1]]$lower_coverage_cum[-1])
ylim_c = c( min(aa, na.rm=T), 1)

aa = c(results_daily[[1]]$power01[-1], results_week[[1]]$power01[-1], 
       results_month[[1]]$power01[-1], results_dow[[1]]$power01[-1],
       results_daily[[1]]$power99[-1], results_week[[1]]$power99[-1], 
       results_month[[1]]$power99[-1], results_dow[[1]]$power99[-1])
ylim_d = c( min(aa, na.rm=T), 1)


aa = c(results_daily[[2]]$lower_bias_cum[-1], results_daily[[2]]$upper_bias_cum[-1],
       results_week[[2]]$lower_bias_cum[-1], results_week[[2]]$upper_bias_cum[-1],
       results_month[[2]]$lower_bias_cum[-1], results_month[[2]]$upper_bias_cum[-1],
       results_dow[[2]]$lower_bias_cum[-1], results_dow[[2]]$upper_bias_cum[-1])
ylim_e = c( min(aa, na.rm=T), max(aa, na.rm=T))

aa = c(results_daily[[2]]$rmse_cum[-1], results_week[[2]]$rmse_cum[-1],
       results_month[[2]]$rmse_cum[-1], results_dow[[2]]$rmse_cum[-1])
ylim_f = c( 0, max(aa, na.rm=T))

aa = c(results_daily[[2]]$lower_coverage_cum[-1], results_week[[2]]$lower_coverage_cum[-1], 
       results_month[[2]]$lower_coverage_cum[-1], results_dow[[2]]$lower_coverage_cum[-1])
ylim_g = c( min(aa, na.rm=T), 1)

aa = c(results_daily[[2]]$power01[-1], results_week[[2]]$power01[-1], 
       results_month[[2]]$power01[-1], results_dow[[2]]$power01[-1],
       results_daily[[2]]$power99[-1], results_week[[2]]$power99[-1], 
       results_month[[2]]$power99[-1], results_dow[[2]]$power99[-1])
ylim_h = c( min(aa, na.rm=T), 1)



aa = c(results_daily[[3]]$lower_bias_cum[-1], results_daily[[3]]$upper_bias_cum[-1],
       results_week[[3]]$lower_bias_cum[-1], results_week[[3]]$upper_bias_cum[-1],
       results_month[[3]]$lower_bias_cum[-1], results_month[[3]]$upper_bias_cum[-1],
       results_dow[[3]]$lower_bias_cum[-1], results_dow[[3]]$upper_bias_cum[-1])
ylim_i = c( min(aa, na.rm=T), max(aa, na.rm=T))

aa = c(results_daily[[3]]$rmse_cum[-1], results_week[[3]]$rmse_cum[-1], results_month[[3]]$rmse_cum[-1],
       results_dow[[3]]$rmse_cum[-1])
ylim_j = c( 0, max(aa, na.rm=T))

aa = c(results_daily[[3]]$lower_coverage_cum[-1], results_week[[3]]$lower_coverage_cum[-1], 
       results_month[[3]]$lower_coverage_cum[-1], results_dow[[3]]$lower_coverage_cum[-1])
ylim_k = c( min(aa, na.rm=T), 1)

aa = c(results_daily[[3]]$power[-1], results_week[[3]]$power[-1], 
       results_month[[3]]$power[-1], results_dow[[3]]$power[-1])
ylim_l = c( min(aa, na.rm=T), 1)


aa = c(results_daily[[4]]$lower_bias_cum[-1], results_daily[[4]]$upper_bias_cum[-1],
       results_week[[4]]$lower_bias_cum[-1], results_week[[4]]$upper_bias_cum[-1],
       results_month[[4]]$lower_bias_cum[-1], results_month[[4]]$upper_bias_cum[-1],
       results_dow[[4]]$lower_bias_cum[-1], results_dow[[4]]$upper_bias_cum[-1])
ylim_m = c( min(aa, na.rm=T), max(aa, na.rm=T))

aa = c(results_daily[[4]]$rmse_cum[-1], results_week[[4]]$rmse_cum[-1],
       results_month[[4]]$rmse_cum[-1], results_dow[[4]]$rmse_cum[-1])
ylim_n = c( 0, max(aa, na.rm=T))

aa = c(results_daily[[4]]$lower_coverage_cum[-1], results_week[[4]]$lower_coverage_cum[-1], 
       results_month[[4]]$lower_coverage_cum[-1], results_dow[[4]]$lower_coverage_cum[-1])
ylim_o = c( min(aa, na.rm=T), 1)

aa = c(results_daily[[4]]$power[-1], results_week[[4]]$power[-1], 
       results_month[[4]]$power[-1], results_dow[[4]]$power[-1])
ylim_p = c( min(aa, na.rm=T), 1)



ylim1 = matrix( c(ylim_a, ylim_e, ylim_i, ylim_m),
                byrow=T, ncol=2)

ylim2 = matrix( c(ylim_b, ylim_f, ylim_j, ylim_n),
                byrow=T, ncol=2)

ylim3 = matrix( c(ylim_c, ylim_g, ylim_k, ylim_o),
                byrow=T, ncol=2)

ylim4 = matrix( c(ylim_d, ylim_h, ylim_l, ylim_p),
                byrow=T, ncol=2)

lets = matrix(LETTERS[1:16], byrow=T, ncol=4)
at_lets = -6
cex_lets = 1.2

titles = c("Mortality | Temperature", "Hospit. | Temperature", 
           expression('Mortality | NO'[2]), expression('Hospit. | NO'[2]))

# imax = c(13-1, 13-1, 13-1, 13-1)

imax = c(13-1, 13-1, 13-1, 13-1)

indyears = c(2:10, 15, 20, 25)

pdf(file = "Figure bias-mse-coverage-power E-R no1year 0101_1pct dow_holi new temp.pdf",
    width = 8.2, height = 8.2 + 8.2/3)

par(mfrow=c(4,4))

for (i in 1:4) {
  
  # Bias
  
  plot(indyears[1:imax[i]], results_daily[[i]]$bias_cum[-1], type="l", ylab="Bias", 
       xlab="Number of years", lwd=5, cex.lab=1.2, ylim=ylim1[i,], xlim=c(1,25))
  polygon(c(indyears[1:imax[i]], rev(indyears[1:imax[i]])),
          c(results_daily[[i]]$lower_bias_cum[-1], 
            rev(results_daily[[i]]$upper_bias_cum[-1])),
          col=rgb(0, 0, 0,0.3), border=FALSE, lwd=1)
  lines(indyears[1:imax[i]], results_week[[i]]$bias_cum[-1], col="red", lwd=5)
  polygon(c(indyears[1:imax[i]], rev(indyears[1:imax[i]])),
          c(results_week[[i]]$lower_bias_cum[-1], 
            rev(results_week[[i]]$upper_bias_cum[-1])),
          col=rgb(1, 0, 0,0.3), border=FALSE, lwd=1)
  lines(indyears[1:imax[i]], results_dow[[i]]$bias_cum[-1], col="#44AA99", lwd=5)
  polygon(c(indyears[1:imax[i]], rev(indyears[1:imax[i]])),
          c(results_dow[[i]]$lower_bias_cum[-1], 
            rev(results_dow[[i]]$upper_bias_cum[-1])),
          col=rgb(0.27, 0.67, 0.6, 0.3), border=FALSE, lwd=1)
  
  lines(indyears[1:imax[i]], results_month[[i]]$bias_cum[-1], col="blue", lwd=5)
  polygon(c(indyears[1:imax[i]], rev(indyears[1:imax[i]])),
          c(results_month[[i]]$lower_bias_cum[-1], 
            rev(results_month[[i]]$upper_bias_cum[-1])),
          col=rgb(0, 0, 1,0.3), border=FALSE, lwd=1)
  
  abline(h=0,lty=1)
  mtext(titles[i], side=3, line=0.5, cex=.8, font=2)
  if (i==1) {
    legend("bottomright",c("D|D", "W|D","Dow|D", "M|D"), col=c("black","red", "#44AA99", "blue"), lty=1,
           lwd=5, cex=1.2, bty="n")
  }
  mtext(lets[i,1], side=3, at=at_lets, cex=cex_lets, line=1, font=2)
  
  
  # RMSE
  
  plot(indyears[1:imax[i]], results_daily[[i]]$rmse_cum[-1], type="l", ylab="RMSE", 
       xlab="Number of years",  lwd=5, cex.lab=1.2, ylim=ylim2[i,], xlim=c(1,25))
  lines(indyears[1:imax[i]], results_week[[i]]$rmse_cum[-1], col="red", lwd=5)
  lines(indyears[1:imax[i]], results_dow[[i]]$rmse_cum[-1], col="#44AA99", lwd=5)
  lines(indyears[1:imax[i]], results_month[[i]]$rmse_cum[-1], col="blue", lwd=5)
  
  abline(h=0,lty=1)
  #legend("topright",c("D|D", "W|D","Dow|D", "M|D"), col=c("black","red", "#44AA99", "blue"), lty=1,
  #       lwd=5, cex=1.2)
  
  mtext(titles[i], side=3, line=0.5, cex=.8, font=2)
  mtext(lets[i,2], side=3, at=at_lets, cex=cex_lets, line=1, font=2)
  
  
  # Coverage
  
  plot(indyears[1:imax[i]], results_daily[[i]]$coverage_cum[-1], type="l", ylab="Coverage", 
       xlab="Number of years", lwd=5, cex.lab=1.2, ylim=ylim3[i,], xlim=c(1,25))
  lines(indyears[1:imax[i]], results_week[[i]]$coverage_cum[-1], col="red", lwd=5)
  lines(indyears[1:imax[i]], results_dow[[i]]$coverage_cum[-1], col="#44AA99", lwd=5)
  polygon(c(indyears[1:imax[i]], rev(indyears[1:imax[i]])),
          c(results_week[[i]]$lower_coverage_cum[-1], 
            rev(results_week[[i]]$upper_coverage_cum[-1])),
          col=rgb(1, 0, 0,0.3), border=FALSE, lwd=1)
  polygon(c(indyears[1:imax[i]], rev(indyears[1:imax[i]])),
          c(results_dow[[i]]$lower_coverage_cum[-1], 
            rev(results_dow[[i]]$upper_coverage_cum[-1])),
          col=rgb(0.27, 0.67, 0.6,0.3), border=FALSE, lwd=1)
  lines(indyears[1:imax[i]], results_month[[i]]$coverage_cum[-1], col="blue", lwd=5)
  polygon(c(indyears[1:imax[i]], rev(indyears[1:imax[i]])),
          c(results_month[[i]]$lower_coverage_cum[-1], 
            rev(results_month[[i]]$upper_coverage_cum[-1])),
          col=rgb(0, 0, 1,0.3), border=FALSE, lwd=1)
  polygon(c(indyears[1:imax[i]], rev(indyears[1:imax[i]])),
          c(results_daily[[i]]$lower_coverage_cum[-1], 
            rev(results_daily[[i]]$upper_coverage_cum[-1])),
          col=rgb(0, 0, 0,0.3), border=FALSE, lwd=1)
  abline(h=0.95,lty=1)
  #legend("bottomright",c("D|D", "W|D","Dow|D", "M|D"), col=c("black","red", "#44AA99", "blue"), lty=1,
  #       lwd=5, cex=1.2)
  mtext(titles[i], side=3, line=0.5, cex=.8, font=2)
  mtext(lets[i,3], side=3, at=at_lets, cex=cex_lets, line=1, font=2)

  # Power
  
  if (i==1 | i==2) {
    plot(indyears[1:imax[i]], results_daily[[i]]$power99[-1], type="l", ylab="Power", 
       xlab="Number of years", lwd=5, cex.lab=1.2, ylim=ylim4[i,], xlim=c(1,25))
    lines(indyears[1:imax[i]], results_week[[i]]$power99[-1], col="red", lwd=5)
    lines(indyears[1:imax[i]], results_dow[[i]]$power99[-1], col="#44AA99", lwd=5)
    lines(indyears[1:imax[i]], results_month[[i]]$power99[-1], col="blue", lwd=5)
    
    lines(indyears[1:imax[i]], results_daily[[i]]$power01[-1], col="black", lwd=2)
    lines(indyears[1:imax[i]], results_week[[i]]$power01[-1], col="red", lwd=2)
    lines(indyears[1:imax[i]], results_dow[[i]]$power01[-1], col="#44AA99", lwd=2)
    lines(indyears[1:imax[i]], results_month[[i]]$power01[-1], col="blue", lwd=2)
    
    #abline(h=0.95,lty=1)
    mtext(titles[i], side=3, line=0.5, cex=.8, font=2)
    mtext(lets[i,4], side=3, at=at_lets, cex=cex_lets, line=1, font=2)
  } else if (i==3 | i==4) {
    plot(indyears[1:imax[i]], results_daily[[i]]$power[-1], type="l", ylab="Power", 
         xlab="Number of years", lwd=5, cex.lab=1.2, ylim=ylim4[i,], xlim=c(1,25))
    lines(indyears[1:imax[i]], results_week[[i]]$power[-1], col="red", lwd=5)
    lines(indyears[1:imax[i]], results_dow[[i]]$power[-1], col="#44AA99", lwd=5)
    lines(indyears[1:imax[i]], results_month[[i]]$power[-1], col="blue", lwd=5)
    #abline(h=0.95,lty=1)
    mtext(titles[i], side=3, line=0.5, cex=.8, font=2)
    mtext(lets[i,4], side=3, at=at_lets, cex=cex_lets, line=1, font=2)
    
  }
}

dev.off()


##########################################
##########################################
### Figure for every percentile
##########################################
##########################################

########
########
# Bias #
########
########

### Shown using true and estimated exposure response curves (and their average)
################################################################################

#################
# Using trueMMT #
#################

## Mortality - temperature
############################

ylim_1 = c( 
  min(c( log(min(unlist(lmort_temp$cumpred_daily_fixedmmt[10:13]))),
         log(min(unlist(lmort_temp$cumpred_week_fixedmmt[10:13] ))),
         log(min(unlist(lmort_temp$cumpred_month_fixedmmt[10:13]),na.rm=T)),
         log(min(unlist(lmort_temp$cumpred_dow_fixedmmt[10:13]),na.rm=T)))),
  max(c( log(max(unlist(lmort_temp$cumpred_daily_fixedmmt[10:13]))),
         log(max(unlist(lmort_temp$cumpred_week_fixedmmt[10:13] ))),
         log(max(unlist(lmort_temp$cumpred_month_fixedmmt[10:13]),na.rm=T)),
         log(max(unlist(lmort_temp$cumpred_dow_fixedmmt[10:13]),na.rm=T))))   )

lets = matrix(LETTERS[1:16], byrow=T, ncol=4)
at_lets = -10
cex_lets = 1.2

indices = 10:13
years = seq(10,25,5)


# Calculate true Cumulative E-R curve for all years

pdf(file = "Figure lines mort-temp_trueMMT.pdf",
    width = 8.2, height = 8.2 + 8.2/3)

par(mfrow = c(4,4))

for (i in 1:length(indices)) {
  
  plot(unique(at_temp), results_daily[[1]]$truth.daily[[indices[i]]][1,], col="coral", 
       lwd=6, type="l", xlab="Temperature (ºC)", ylab="log(RR)",ylim=ylim_1, cex.lab=1.2)
  for (j in 1:nrow(lmort_temp$cumpred_daily_fixedmmt[[i]])) {
    lines(unique(at_temp), log(lmort_temp$cumpred_daily_fixedmmt[[indices[i]]][j,]), col="gray80")
  }
  lines(unique(at_temp),  results_daily[[1]]$truth.daily[[indices[i]]][1,], col="coral", 
        lwd=6)
  lines(unique(at_temp), apply(log(lmort_temp$cumpred_daily_fixedmmt[[indices[i]]]),2,mean), 
        col="black", lwd=2)
  mtext(paste0("D|D - ", years[i], " years"), side=3, line=0.5, cex=.8, font=2)
  
  abline(v=unique(at_temp)[c(11,108)],lty=2,col="gray50")
  abline(h=0,lty=1)
  
  #legend("top", c("Truth", "D|D"), lty=1, col=c("coral","black"),
  #     lwd=c(5,2), bg="white", cex=1.2, bty = "n")
  
  mtext(lets[1,i], side=3, at=at_lets, cex=cex_lets, line=1, font=2)
  
}

for (i in 1:length(indices)) {
  
  plot(unique(at_temp), results_week[[1]]$truth.daily[[indices[i]]][1,], col="coral", 
       lwd=6, type="l", xlab="Temperature (ºC)", ylab="log(RR)",ylim=ylim_1, cex.lab=1.2)
  for (j in 1:nrow(lmort_temp$cumpred_week_fixedmmt[[i]])) {
    lines(unique(at_temp), log(lmort_temp$cumpred_week_fixedmmt[[indices[i]]][j,]), col="gray80")
  }
  lines(unique(at_temp),  results_week[[1]]$truth.daily[[indices[i]]][1,], col="coral", 
        lwd=6)
  lines(unique(at_temp), apply(log(lmort_temp$cumpred_week_fixedmmt[[indices[i]]]),2,mean), 
        col="black", lwd=2)
  mtext(paste0("W|D - ", years[i], " years"), side=3, line=0.5, cex=.8, font=2)
  
  abline(v=unique(at_temp)[c(11,108)],lty=2,col="gray50")
  abline(h=0,lty=1)
  
  #legend("top", c("Truth", "D|D"), lty=1, col=c("coral","black"),
  #     lwd=c(5,2), bg="white", cex=1.2, bty = "n")
  
  mtext(lets[2,i], side=3, at=at_lets, cex=cex_lets, line=1, font=2)
  
}


for (i in 1:length(indices)) {
  
  plot(unique(at_temp), results_dow[[1]]$truth.daily[[indices[i]]][1,], col="coral", 
       lwd=6, type="l", xlab="Temperature (ºC)", ylab="log(RR)",ylim=ylim_1, cex.lab=1.2)
  for (j in 1:nrow(lmort_temp$cumpred_dow_fixedmmt[[i]])) {
    lines(unique(at_temp), log(lmort_temp$cumpred_dow_fixedmmt[[indices[i]]][j,]), col="gray80")
  }
  lines(unique(at_temp),  results_dow[[1]]$truth.daily[[indices[i]]][1,], col="coral", 
        lwd=6)
  lines(unique(at_temp), apply(log(lmort_temp$cumpred_dow_fixedmmt[[indices[i]]]),2,mean), 
        col="black", lwd=2)
  mtext(paste0("Dow|D - ", years[i], " years"), side=3, line=0.5, cex=.8, font=2)
  
  abline(v=unique(at_temp)[c(11,108)],lty=2,col="gray50")
  abline(h=0,lty=1)
  
  #legend("top", c("Truth", "D|D"), lty=1, col=c("coral","black"),
  #     lwd=c(5,2), bg="white", cex=1.2, bty = "n")
  
  mtext(lets[3,i], side=3, at=at_lets, cex=cex_lets, line=1, font=2)
  
}

for (i in 1:length(indices)) {
  
  plot(unique(at_temp), results_month[[1]]$truth.daily[[indices[i]]][1,], col="coral", 
       lwd=6, type="l", xlab="Temperature (ºC)", ylab="log(RR)",ylim=ylim_1, cex.lab=1.2)
  for (j in 1:nrow(lmort_temp$cumpred_month_fixedmmt[[i]])) {
    lines(unique(at_temp), log(lmort_temp$cumpred_month_fixedmmt[[indices[i]]][j,]), col="gray80")
  }
  lines(unique(at_temp),  results_month[[1]]$truth.daily[[indices[i]]][1,], col="coral", 
        lwd=6)
  lines(unique(at_temp), apply(log(lmort_temp$cumpred_month_fixedmmt[[indices[i]]]),2,mean, na.rm=T), 
        col="black", lwd=2)
  mtext(paste0("M|D - ", years[i], " years"), side=3, line=0.5, cex=.8, font=2)
  
  abline(v=unique(at_temp)[c(11,108)],lty=2,col="gray50")
  abline(h=0,lty=1)
  
  #legend("top", c("Truth", "D|D"), lty=1, col=c("coral","black"),
  #     lwd=c(5,2), bg="white", cex=1.2, bty = "n")
  
  mtext(lets[4,i], side=3, at=at_lets, cex=cex_lets, line=1, font=2)
  
}


dev.off()




## Hosp - temperature
############################

ylim_1 = c( 
  min(c( log(min(unlist(lhosp_temp$cumpred_daily_fixedmmt[10:13]))),
         log(min(unlist(lhosp_temp$cumpred_week_fixedmmt[10:13] ))),
         log(min(unlist(lhosp_temp$cumpred_month_fixedmmt[10:13]),na.rm=T)),
         log(min(unlist(lhosp_temp$cumpred_dow_fixedmmt[10:13]),na.rm=T)))),
  max(c( log(max(unlist(lhosp_temp$cumpred_daily_fixedmmt[10:13]))),
         log(max(unlist(lhosp_temp$cumpred_week_fixedmmt[10:13] ))),
         log(max(unlist(lhosp_temp$cumpred_month_fixedmmt[10:13]),na.rm=T)),
         log(max(unlist(lhosp_temp$cumpred_dow_fixedmmt[10:13]),na.rm=T))))   )

lets = matrix(LETTERS[1:16], byrow=T, ncol=4)
at_lets = -10
cex_lets = 1.2

indices = 10:13
years = seq(10,25,5)


# Calculate true Cumulative E-R curve for all years

pdf(file = "Figure lines hosp-temp_trueMMT.pdf",
    width = 8.2, height = 8.2 + 8.2/3)

par(mfrow = c(4,4))

for (i in 1:length(indices)) {
  
  plot(unique(at_temp), results_daily[[2]]$truth.daily[[indices[i]]][1,], col="coral", 
       lwd=6, type="l", xlab="Temperature (ºC)", ylab="log(RR)",ylim=ylim_1, cex.lab=1.2)
  for (j in 1:nrow(lhosp_temp$cumpred_daily_fixedmmt[[i]])) {
    lines(unique(at_temp), log(lhosp_temp$cumpred_daily_fixedmmt[[indices[i]]][j,]), col="gray80")
  }
  lines(unique(at_temp),  results_daily[[2]]$truth.daily[[indices[i]]][1,], col="coral", 
        lwd=6)
  lines(unique(at_temp), apply(log(lhosp_temp$cumpred_daily_fixedmmt[[indices[i]]]),2,mean), 
        col="black", lwd=2)
  mtext(paste0("D|D - ", years[i], " years"), side=3, line=0.5, cex=.8, font=2)
  
  abline(v=unique(at_temp)[c(11,108)],lty=2,col="gray50")
  abline(h=0,lty=1)
  
  #legend("top", c("Truth", "D|D"), lty=1, col=c("coral","black"),
  #     lwd=c(5,2), bg="white", cex=1.2, bty = "n")
  
  mtext(lets[1,i], side=3, at=at_lets, cex=cex_lets, line=1, font=2)
  
}

for (i in 1:length(indices)) {
  
  plot(unique(at_temp), results_week[[2]]$truth.daily[[indices[i]]][1,], col="coral", 
       lwd=6, type="l", xlab="Temperature (ºC)", ylab="log(RR)",ylim=ylim_1, cex.lab=1.2)
  for (j in 1:nrow(lhosp_temp$cumpred_week_fixedmmt[[i]])) {
    lines(unique(at_temp), log(lhosp_temp$cumpred_week_fixedmmt[[indices[i]]][j,]), col="gray80")
  }
  lines(unique(at_temp),  results_week[[2]]$truth.daily[[indices[i]]][1,], col="coral", 
        lwd=6)
  lines(unique(at_temp), apply(log(lhosp_temp$cumpred_week_fixedmmt[[indices[i]]]),2,mean), 
        col="black", lwd=2)
  mtext(paste0("W|D - ", years[i], " years"), side=3, line=0.5, cex=.8, font=2)
  
  abline(v=unique(at_temp)[c(11,108)],lty=2,col="gray50")
  abline(h=0,lty=1)
  
  #legend("top", c("Truth", "D|D"), lty=1, col=c("coral","black"),
  #     lwd=c(5,2), bg="white", cex=1.2, bty = "n")
  
  mtext(lets[2,i], side=3, at=at_lets, cex=cex_lets, line=1, font=2)
  
}


for (i in 1:length(indices)) {
  
  plot(unique(at_temp), results_dow[[2]]$truth.daily[[indices[i]]][1,], col="coral", 
       lwd=6, type="l", xlab="Temperature (ºC)", ylab="log(RR)",ylim=ylim_1, cex.lab=1.2)
  for (j in 1:nrow(lhosp_temp$cumpred_dow_fixedmmt[[i]])) {
    lines(unique(at_temp), log(lhosp_temp$cumpred_dow_fixedmmt[[indices[i]]][j,]), col="gray80")
  }
  lines(unique(at_temp),  results_dow[[2]]$truth.daily[[indices[i]]][1,], col="coral", 
        lwd=6)
  lines(unique(at_temp), apply(log(lhosp_temp$cumpred_dow_fixedmmt[[indices[i]]]),2,mean), 
        col="black", lwd=2)
  mtext(paste0("Dow|D - ", years[i], " years"), side=3, line=0.5, cex=.8, font=2)
  
  abline(v=unique(at_temp)[c(11,108)],lty=2,col="gray50")
  abline(h=0,lty=1)
  
  #legend("top", c("Truth", "D|D"), lty=1, col=c("coral","black"),
  #     lwd=c(5,2), bg="white", cex=1.2, bty = "n")
  
  mtext(lets[3,i], side=3, at=at_lets, cex=cex_lets, line=1, font=2)
  
}

for (i in 1:length(indices)) {
  
  plot(unique(at_temp), results_month[[2]]$truth.daily[[indices[i]]][1,], col="coral", 
       lwd=6, type="l", xlab="Temperature (ºC)", ylab="log(RR)",ylim=ylim_1, cex.lab=1.2)
  for (j in 1:nrow(lhosp_temp$cumpred_month_fixedmmt[[i]])) {
    lines(unique(at_temp), log(lhosp_temp$cumpred_month_fixedmmt[[indices[i]]][j,]), col="gray80")
  }
  lines(unique(at_temp),  results_month[[2]]$truth.daily[[indices[i]]][1,], col="coral", 
        lwd=6)
  lines(unique(at_temp), apply(log(lhosp_temp$cumpred_month_fixedmmt[[indices[i]]]),2,mean, na.rm=T), 
        col="black", lwd=2)
  mtext(paste0("M|D - ", years[i], " years"), side=3, line=0.5, cex=.8, font=2)
  
  abline(v=unique(at_temp)[c(11,108)],lty=2,col="gray50")
  abline(h=0,lty=1)
  
  #legend("top", c("Truth", "D|D"), lty=1, col=c("coral","black"),
  #     lwd=c(5,2), bg="white", cex=1.2, bty = "n")
  
  mtext(lets[4,i], side=3, at=at_lets, cex=cex_lets, line=1, font=2)
  
}


dev.off()




########
########
# RMSE #
########
########

aa = c(results_month[[1]]$rmse_perc[,10],
       results_month[[1]]$rmse_perc[,11],
       results_month[[1]]$rmse_perc[,12],
       results_month[[1]]$rmse_perc[,13],
       results_month[[2]]$rmse_perc[,10],
       results_month[[2]]$rmse_perc[,11],
       results_month[[2]]$rmse_perc[,12],
       results_month[[2]]$rmse_perc[,13])

ylim1 = c( min(aa, na.rm=T), max(aa, na.rm=T))

lets = matrix(LETTERS[1:8], byrow=T, ncol=2)
at_lets = -10
cex_lets = 1.2

titles = c("Mort, 10 years", "Hosp, 10 years", 
           "Mort, 15 years", "Hosp, 15 years",
           "Mort, 20 years", "Hosp, 20 years", 
           "Mort, 25 years", "Hosp, 25 years")


pdf(file = "Figure RMSE E-R percentiles.pdf",
    width = 8.2*2/3, height = 8.2 + 8.2/3)


par(mfrow = c(4,2))

ylab1 = "RMSE"

plot(0:100, results_daily[[1]]$rmse_perc[,10], type="l", lwd=3, ylim = c(0, max(ylim1)), ylab = ylab1,
     xlab = "Percentile", main = titles[1])
lines(0:100, results_week[[1]]$rmse_perc[,10], col="red", lwd=3)
lines(0:100, results_month[[1]]$rmse_perc[,10], col="blue", lwd=3)
lines(0:100, results_dow[[1]]$rmse_perc[,10], col="#44AA99", lwd=3)
abline(h=0, lty = 2)
mtext(lets[1,1], side=3, at=at_lets, cex=cex_lets, line=1, font=2)
legend("top",c("D|D", "W|D","Dow|D", "M|D"), col=c("black","red", "#44AA99", "blue"), lty=1,
       lwd=3, cex=1.2, bty="n")

plot(0:100, results_daily[[2]]$rmse_perc[,10], type="l", lwd=3, ylim = c(0, max(ylim1)), ylab = ylab1,
     xlab = "Percentile", main = titles[2])
lines(0:100, results_week[[2]]$rmse_perc[,10], col="red", lwd=3)
lines(0:100, results_month[[2]]$rmse_perc[,10], col="blue", lwd=3)
lines(0:100, results_dow[[2]]$rmse_perc[,10], col="#44AA99", lwd=3)
abline(h=0, lty = 2)
mtext(lets[1,2], side=3, at=at_lets, cex=cex_lets, line=1, font=2)

plot(0:100, results_daily[[1]]$rmse_perc[,11], type="l", lwd=3, ylim = c(0, max(ylim1)), ylab = ylab1,
     xlab = "Percentile", main = titles[3])
lines(0:100, results_week[[1]]$rmse_perc[,11], col="red", lwd=3)
lines(0:100, results_month[[1]]$rmse_perc[,11], col="blue", lwd=3)
lines(0:100, results_dow[[1]]$rmse_perc[,11], col="#44AA99", lwd=3)
abline(h=0, lty = 2)
mtext(lets[2,1], side=3, at=at_lets, cex=cex_lets, line=1, font=2)

plot(0:100, results_daily[[2]]$rmse_perc[,11], type="l", lwd=3, ylim = c(0, max(ylim1)), ylab = ylab1,
     xlab = "Percentile", main = titles[4])
lines(0:100, results_week[[2]]$rmse_perc[,11], col="red", lwd=3)
lines(0:100, results_month[[2]]$rmse_perc[,11], col="blue", lwd=3)
lines(0:100, results_dow[[2]]$rmse_perc[,11], col="#44AA99", lwd=3)
abline(h=0, lty = 2)
mtext(lets[2,2], side=3, at=at_lets, cex=cex_lets, line=1, font=2)


plot(0:100, results_daily[[1]]$rmse_perc[,12], type="l", lwd=3, ylim = c(0, max(ylim1)), ylab = ylab1,
     xlab = "Percentile", main = titles[5])
lines(0:100, results_week[[1]]$rmse_perc[,12], col="red", lwd=3)
lines(0:100, results_month[[1]]$rmse_perc[,12], col="blue", lwd=3)
lines(0:100, results_dow[[1]]$rmse_perc[,12], col="#44AA99", lwd=3)
abline(h=0, lty = 2)
mtext(lets[3,1], side=3, at=at_lets, cex=cex_lets, line=1, font=2)

plot(0:100, results_daily[[2]]$rmse_perc[,12], type="l", lwd=3, ylim = c(0, max(ylim1)), ylab = ylab1,
     xlab = "Percentile", main = titles[6])
lines(0:100, results_week[[2]]$rmse_perc[,12], col="red", lwd=3)
lines(0:100, results_month[[2]]$rmse_perc[,12], col="blue", lwd=3)
lines(0:100, results_dow[[2]]$rmse_perc[,12], col="#44AA99", lwd=3)
abline(h=0, lty = 2)
mtext(lets[3,2], side=3, at=at_lets, cex=cex_lets, line=1, font=2)

plot(0:100, results_daily[[1]]$rmse_perc[,13], type="l", lwd=3, ylim = c(0, max(ylim1)), ylab = ylab1,
     xlab = "Percentile", main = titles[7])
lines(0:100, results_week[[1]]$rmse_perc[,13], col="red", lwd=3)
lines(0:100, results_month[[1]]$rmse_perc[,13], col="blue", lwd=3)
lines(0:100, results_dow[[1]]$rmse_perc[,13], col="#44AA99", lwd=3)
abline(h=0, lty = 2)
mtext(lets[4,1], side=3, at=at_lets, cex=cex_lets, line=1, font=2)

plot(0:100, results_daily[[2]]$rmse_perc[,13], type="l", lwd=3, ylim = c(0, max(ylim1)), ylab = ylab1,
     xlab = "Percentile", main = titles[8])
lines(0:100, results_week[[2]]$rmse_perc[,13], col="red", lwd=3)
lines(0:100, results_month[[2]]$rmse_perc[,13], col="blue", lwd=3)
lines(0:100, results_dow[[2]]$rmse_perc[,13], col="#44AA99", lwd=3)
abline(h=0, lty = 2)
mtext(lets[4,2], side=3, at=at_lets, cex=cex_lets, line=1, font=2)

dev.off()



############
############
# Coverage #
############
############

#######
# D|D
#######

anal = "mort_temp"
load(file=paste0("parameters_", anal, ".RData"))
excl = as.numeric(substr(names(which(at_temp == MMTtrue)), 1, 2))
percs_mort = c(0:(excl-1), (excl+1):100 )

anal = "hosp_temp"
load(file=paste0("parameters_", anal, ".RData"))
excl = as.numeric(substr(names(which(at_temp == MMTtrue)), 1, 2))
percs_hosp = c(0:(excl-1), (excl+1):100 )

aa = c(results_daily[[1]]$lower_coverage_perc[,10], results_daily[[1]]$upper_coverage_perc[,10],
       results_daily[[1]]$lower_coverage_perc[,11], results_daily[[1]]$upper_coverage_perc[,11],
       results_daily[[1]]$lower_coverage_perc[,12], results_daily[[1]]$upper_coverage_perc[,12],
       results_daily[[1]]$lower_coverage_perc[,13], results_daily[[1]]$upper_coverage_perc[,13],
       results_daily[[2]]$lower_coverage_perc[,10], results_daily[[2]]$upper_coverage_perc[,10],
       results_daily[[2]]$lower_coverage_perc[,11], results_daily[[2]]$upper_coverage_perc[,11],
       results_daily[[2]]$lower_coverage_perc[,12], results_daily[[2]]$upper_coverage_perc[,12],
       results_daily[[2]]$lower_coverage_perc[,13], results_daily[[2]]$upper_coverage_perc[,13])

ylim1 = c( min(aa, na.rm=T), max(aa, na.rm=T))
ylim1 = c( min(aa, na.rm=T), 1)

lets = matrix(LETTERS[1:8], byrow=T, ncol=2)
at_lets = -10
cex_lets = 1.2

titles = c("Mort, D|D, 10 years", "Hosp, D|D, 10 years", 
           "Mort, D|D, 15 years", "Hosp, D|D, 15 years",
           "Mort, D|D, 20 years", "Hosp, D|D, 20 years", 
           "Mort, D|D, 25 years", "Hosp, D|D, 25 years")

ylab1 = "Coverage (%)"

pdf(file = "Figure Coverage E-R percentiles DD.pdf",
    width = 8.2*2/3, height = 8.2 + 8.2/3)


par(mfrow = c(4,2))


plot(percs_mort, results_daily[[1]]$coverage_perc[,10], type="l", lwd=3, ylim = ylim1, ylab = ylab1,
     xlab = "Percentile", main = titles[1])
polygon(c(percs_mort, rev(percs_mort)), 
        c(results_daily[[1]]$lower_coverage_perc[,10], rev(results_daily[[1]]$upper_coverage_perc[,10])),
        col=rgb(0, 0, 0,0.1), border=FALSE)
abline(h=0.95, lty = 2)
mtext(lets[1,1], side=3, at=at_lets, cex=cex_lets, line=1, font=2)

plot(percs_hosp, results_daily[[2]]$coverage_perc[,10], type="l", lwd=3, ylim = ylim1, ylab = ylab1,
     xlab = "Percentile", main = titles[2])
polygon(c(percs_hosp, rev(percs_hosp)), 
        c(results_daily[[2]]$lower_coverage_perc[,10], rev(results_daily[[2]]$upper_coverage_perc[,10])),
        col=rgb(0, 0, 0,0.1), border=FALSE)
abline(h=0.95, lty = 2)
mtext(lets[1,2], side=3, at=at_lets, cex=cex_lets, line=1, font=2)

plot(percs_mort, results_daily[[1]]$coverage_perc[,11], type="l", lwd=3, ylim = ylim1, ylab = ylab1,
     xlab = "Percentile", main = titles[3])
polygon(c(percs_mort, rev(percs_mort)), 
        c(results_daily[[1]]$lower_coverage_perc[,11], rev(results_daily[[1]]$upper_coverage_perc[,11])),
        col=rgb(0, 0, 0,0.1), border=FALSE)
abline(h=0.95, lty = 2)
mtext(lets[2,1], side=3, at=at_lets, cex=cex_lets, line=1, font=2)

plot(percs_hosp, results_daily[[2]]$coverage_perc[,11], type="l", lwd=3, ylim = ylim1, ylab = ylab1,
     xlab = "Percentile", main = titles[4])
polygon(c(percs_hosp, rev(percs_hosp)), 
        c(results_daily[[2]]$lower_coverage_perc[,11], rev(results_daily[[2]]$upper_coverage_perc[,11])),
        col=rgb(0, 0, 0,0.1), border=FALSE)
abline(h=0.95, lty = 2)
mtext(lets[2,2], side=3, at=at_lets, cex=cex_lets, line=1, font=2)


plot(percs_mort, results_daily[[1]]$coverage_perc[,12], type="l", lwd=3, ylim = ylim1, ylab = ylab1,
     xlab = "Percentile", main = titles[5])
polygon(c(percs_mort, rev(percs_mort)), 
        c(results_daily[[1]]$lower_coverage_perc[,12], rev(results_daily[[1]]$upper_coverage_perc[,12])),
        col=rgb(0, 0, 0,0.1), border=FALSE)
abline(h=0.95, lty = 2)
mtext(lets[3,1], side=3, at=at_lets, cex=cex_lets, line=1, font=2)

plot(percs_hosp, results_daily[[2]]$coverage_perc[,12], type="l", lwd=3, ylim = ylim1, ylab = ylab1,
     xlab = "Percentile", main = titles[6])
polygon(c(percs_hosp, rev(percs_hosp)), 
        c(results_daily[[2]]$lower_coverage_perc[,12], rev(results_daily[[2]]$upper_coverage_perc[,12])),
        col=rgb(0, 0, 0,0.1), border=FALSE)
abline(h=0.95, lty = 2)
mtext(lets[3,2], side=3, at=at_lets, cex=cex_lets, line=1, font=2)

plot(percs_mort, results_daily[[1]]$coverage_perc[,13], type="l", lwd=3, ylim = ylim1, ylab = ylab1,
     xlab = "Percentile", main = titles[7])
polygon(c(percs_mort, rev(percs_mort)), 
        c(results_daily[[1]]$lower_coverage_perc[,13], rev(results_daily[[1]]$upper_coverage_perc[,13])),
        col=rgb(0, 0, 0,0.1), border=FALSE)
abline(h=0.95, lty = 2)
mtext(lets[4,1], side=3, at=at_lets, cex=cex_lets, line=1, font=2)

plot(percs_hosp, results_daily[[2]]$coverage_perc[,13], type="l", lwd=3, ylim = ylim1, ylab = ylab1,
     xlab = "Percentile", main = titles[8])
polygon(c(percs_hosp, rev(percs_hosp)), 
        c(results_daily[[2]]$lower_coverage_perc[,13], rev(results_daily[[2]]$upper_coverage_perc[,13])),
        col=rgb(0, 0, 0,0.1), border=FALSE)
abline(h=.95, lty = 2)
mtext(lets[4,2], side=3, at=at_lets, cex=cex_lets, line=1, font=2)

dev.off()



#######
# W|D
#######

aa = c(results_week[[1]]$lower_coverage_perc[,10], results_week[[1]]$upper_coverage_perc[,10],
       results_week[[1]]$lower_coverage_perc[,11], results_week[[1]]$upper_coverage_perc[,11],
       results_week[[1]]$lower_coverage_perc[,12], results_week[[1]]$upper_coverage_perc[,12],
       results_week[[1]]$lower_coverage_perc[,13], results_week[[1]]$upper_coverage_perc[,13],
       results_week[[2]]$lower_coverage_perc[,10], results_week[[2]]$upper_coverage_perc[,10],
       results_week[[2]]$lower_coverage_perc[,11], results_week[[2]]$upper_coverage_perc[,11],
       results_week[[2]]$lower_coverage_perc[,12], results_week[[2]]$upper_coverage_perc[,12],
       results_week[[2]]$lower_coverage_perc[,13], results_week[[2]]$upper_coverage_perc[,13])

ylim1 = c( min(aa, na.rm=T), max(aa, na.rm=T))
ylim1 = c( min(aa, na.rm=T), 1)

lets = matrix(LETTERS[1:8], byrow=T, ncol=2)
at_lets = -10
cex_lets = 1.2

titles = c("Mort, W|D, 10 years", "Hosp, W|D, 10 years", 
           "Mort, W|D, 15 years", "Hosp, W|D, 15 years",
           "Mort, W|D, 20 years", "Hosp, W|D, 20 years", 
           "Mort, W|D, 25 years", "Hosp, W|D, 25 years")

ylab1 = "Coverage (%)"

pdf(file = "Figure Coverage E-R percentiles WD.pdf",
    width = 8.2*2/3, height = 8.2 + 8.2/3)


par(mfrow = c(4,2))


plot(percs_mort, results_week[[1]]$coverage_perc[,10], type="l", lwd=3, ylim = ylim1, ylab = ylab1,
     xlab = "Percentile", main = titles[1])
polygon(c(percs_mort, rev(percs_mort)), 
        c(results_week[[1]]$lower_coverage_perc[,10], rev(results_week[[1]]$upper_coverage_perc[,10])),
        col=rgb(0, 0, 0,0.1), border=FALSE)
abline(h=0.95, lty = 2)
mtext(lets[1,1], side=3, at=at_lets, cex=cex_lets, line=1, font=2)

plot(percs_hosp, results_week[[2]]$coverage_perc[,10], type="l", lwd=3, ylim = ylim1, ylab = ylab1,
     xlab = "Percentile", main = titles[2])
polygon(c(percs_hosp, rev(percs_hosp)), 
        c(results_week[[2]]$lower_coverage_perc[,10], rev(results_week[[2]]$upper_coverage_perc[,10])),
        col=rgb(0, 0, 0,0.1), border=FALSE)
abline(h=0.95, lty = 2)
mtext(lets[1,2], side=3, at=at_lets, cex=cex_lets, line=1, font=2)

plot(percs_mort, results_week[[1]]$coverage_perc[,11], type="l", lwd=3, ylim = ylim1, ylab = ylab1,
     xlab = "Percentile", main = titles[3])
polygon(c(percs_mort, rev(percs_mort)), 
        c(results_week[[1]]$lower_coverage_perc[,11], rev(results_week[[1]]$upper_coverage_perc[,11])),
        col=rgb(0, 0, 0,0.1), border=FALSE)
abline(h=0.95, lty = 2)
mtext(lets[2,1], side=3, at=at_lets, cex=cex_lets, line=1, font=2)

plot(percs_hosp, results_week[[2]]$coverage_perc[,11], type="l", lwd=3, ylim = ylim1, ylab = ylab1,
     xlab = "Percentile", main = titles[4])
polygon(c(percs_hosp, rev(percs_hosp)), 
        c(results_week[[2]]$lower_coverage_perc[,11], rev(results_week[[2]]$upper_coverage_perc[,11])),
        col=rgb(0, 0, 0,0.1), border=FALSE)
abline(h=0.95, lty = 2)
mtext(lets[2,2], side=3, at=at_lets, cex=cex_lets, line=1, font=2)


plot(percs_mort, results_week[[1]]$coverage_perc[,12], type="l", lwd=3, ylim = ylim1, ylab = ylab1,
     xlab = "Percentile", main = titles[5])
polygon(c(percs_mort, rev(percs_mort)), 
        c(results_week[[1]]$lower_coverage_perc[,12], rev(results_week[[1]]$upper_coverage_perc[,12])),
        col=rgb(0, 0, 0,0.1), border=FALSE)
abline(h=0.95, lty = 2)
mtext(lets[3,1], side=3, at=at_lets, cex=cex_lets, line=1, font=2)

plot(percs_hosp, results_week[[2]]$coverage_perc[,12], type="l", lwd=3, ylim = ylim1, ylab = ylab1,
     xlab = "Percentile", main = titles[6])
polygon(c(percs_hosp, rev(percs_hosp)), 
        c(results_week[[2]]$lower_coverage_perc[,12], rev(results_week[[2]]$upper_coverage_perc[,12])),
        col=rgb(0, 0, 0,0.1), border=FALSE)
abline(h=0.95, lty = 2)
mtext(lets[3,2], side=3, at=at_lets, cex=cex_lets, line=1, font=2)

plot(percs_mort, results_week[[1]]$coverage_perc[,13], type="l", lwd=3, ylim = ylim1, ylab = ylab1,
     xlab = "Percentile", main = titles[7])
polygon(c(percs_mort, rev(percs_mort)), 
        c(results_week[[1]]$lower_coverage_perc[,13], rev(results_week[[1]]$upper_coverage_perc[,13])),
        col=rgb(0, 0, 0,0.1), border=FALSE)
abline(h=0.95, lty = 2)
mtext(lets[4,1], side=3, at=at_lets, cex=cex_lets, line=1, font=2)

plot(percs_hosp, results_week[[2]]$coverage_perc[,13], type="l", lwd=3, ylim = ylim1, ylab = ylab1,
     xlab = "Percentile", main = titles[8])
polygon(c(percs_hosp, rev(percs_hosp)), 
        c(results_week[[2]]$lower_coverage_perc[,13], rev(results_week[[2]]$upper_coverage_perc[,13])),
        col=rgb(0, 0, 0,0.1), border=FALSE)
abline(h=.95, lty = 2)
mtext(lets[4,2], side=3, at=at_lets, cex=cex_lets, line=1, font=2)

dev.off()


#######
# M|D
#######


aa = c(results_month[[1]]$lower_coverage_perc[,10], results_month[[1]]$upper_coverage_perc[,10],
       results_month[[1]]$lower_coverage_perc[,11], results_month[[1]]$upper_coverage_perc[,11],
       results_month[[1]]$lower_coverage_perc[,12], results_month[[1]]$upper_coverage_perc[,12],
       results_month[[1]]$lower_coverage_perc[,13], results_month[[1]]$upper_coverage_perc[,13],
       results_month[[2]]$lower_coverage_perc[,10], results_month[[2]]$upper_coverage_perc[,10],
       results_month[[2]]$lower_coverage_perc[,11], results_month[[2]]$upper_coverage_perc[,11],
       results_month[[2]]$lower_coverage_perc[,12], results_month[[2]]$upper_coverage_perc[,12],
       results_month[[2]]$lower_coverage_perc[,13], results_month[[2]]$upper_coverage_perc[,13])

ylim1 = c( min(aa, na.rm=T), max(aa, na.rm=T))
ylim1 = c( min(aa, na.rm=T), 1)

lets = matrix(LETTERS[1:8], byrow=T, ncol=2)
at_lets = -10
cex_lets = 1.2

titles = c("Mort, M|D, 10 years", "Hosp, M|D, 10 years", 
           "Mort, M|D, 15 years", "Hosp, M|D, 15 years",
           "Mort, M|D, 20 years", "Hosp, M|D, 20 years", 
           "Mort, M|D, 25 years", "Hosp, M|D, 25 years")

ylab1 = "Coverage (%)"

pdf(file = "Figure Coverage E-R percentiles MD.pdf",
    width = 8.2*2/3, height = 8.2 + 8.2/3)


par(mfrow = c(4,2))


plot(percs_mort, results_month[[1]]$coverage_perc[,10], type="l", lwd=3, ylim = ylim1, ylab = ylab1,
     xlab = "Percentile", main = titles[1])
polygon(c(percs_mort, rev(percs_mort)), 
        c(results_month[[1]]$lower_coverage_perc[,10], rev(results_month[[1]]$upper_coverage_perc[,10])),
        col=rgb(0, 0, 0,0.1), border=FALSE)
abline(h=0.95, lty = 2)
mtext(lets[1,1], side=3, at=at_lets, cex=cex_lets, line=1, font=2)

plot(percs_hosp, results_month[[2]]$coverage_perc[,10], type="l", lwd=3, ylim = ylim1, ylab = ylab1,
     xlab = "Percentile", main = titles[2])
polygon(c(percs_hosp, rev(percs_hosp)), 
        c(results_month[[2]]$lower_coverage_perc[,10], rev(results_month[[2]]$upper_coverage_perc[,10])),
        col=rgb(0, 0, 0,0.1), border=FALSE)
abline(h=0.95, lty = 2)
mtext(lets[1,2], side=3, at=at_lets, cex=cex_lets, line=1, font=2)

plot(percs_mort, results_month[[1]]$coverage_perc[,11], type="l", lwd=3, ylim = ylim1, ylab = ylab1,
     xlab = "Percentile", main = titles[3])
polygon(c(percs_mort, rev(percs_mort)), 
        c(results_month[[1]]$lower_coverage_perc[,11], rev(results_month[[1]]$upper_coverage_perc[,11])),
        col=rgb(0, 0, 0,0.1), border=FALSE)
abline(h=0.95, lty = 2)
mtext(lets[2,1], side=3, at=at_lets, cex=cex_lets, line=1, font=2)

plot(percs_hosp, results_month[[2]]$coverage_perc[,11], type="l", lwd=3, ylim = ylim1, ylab = ylab1,
     xlab = "Percentile", main = titles[4])
polygon(c(percs_hosp, rev(percs_hosp)), 
        c(results_month[[2]]$lower_coverage_perc[,11], rev(results_month[[2]]$upper_coverage_perc[,11])),
        col=rgb(0, 0, 0,0.1), border=FALSE)
abline(h=0.95, lty = 2)
mtext(lets[2,2], side=3, at=at_lets, cex=cex_lets, line=1, font=2)


plot(percs_mort, results_month[[1]]$coverage_perc[,12], type="l", lwd=3, ylim = ylim1, ylab = ylab1,
     xlab = "Percentile", main = titles[5])
polygon(c(percs_mort, rev(percs_mort)), 
        c(results_month[[1]]$lower_coverage_perc[,12], rev(results_month[[1]]$upper_coverage_perc[,12])),
        col=rgb(0, 0, 0,0.1), border=FALSE)
abline(h=0.95, lty = 2)
mtext(lets[3,1], side=3, at=at_lets, cex=cex_lets, line=1, font=2)

plot(percs_hosp, results_month[[2]]$coverage_perc[,12], type="l", lwd=3, ylim = ylim1, ylab = ylab1,
     xlab = "Percentile", main = titles[6])
polygon(c(percs_hosp, rev(percs_hosp)), 
        c(results_month[[2]]$lower_coverage_perc[,12], rev(results_month[[2]]$upper_coverage_perc[,12])),
        col=rgb(0, 0, 0,0.1), border=FALSE)
abline(h=0.95, lty = 2)
mtext(lets[3,2], side=3, at=at_lets, cex=cex_lets, line=1, font=2)

plot(percs_mort, results_month[[1]]$coverage_perc[,13], type="l", lwd=3, ylim = ylim1, ylab = ylab1,
     xlab = "Percentile", main = titles[7])
polygon(c(percs_mort, rev(percs_mort)), 
        c(results_month[[1]]$lower_coverage_perc[,13], rev(results_month[[1]]$upper_coverage_perc[,13])),
        col=rgb(0, 0, 0,0.1), border=FALSE)
abline(h=0.95, lty = 2)
mtext(lets[4,1], side=3, at=at_lets, cex=cex_lets, line=1, font=2)

plot(percs_hosp, results_month[[2]]$coverage_perc[,13], type="l", lwd=3, ylim = ylim1, ylab = ylab1,
     xlab = "Percentile", main = titles[8])
polygon(c(percs_hosp, rev(percs_hosp)), 
        c(results_month[[2]]$lower_coverage_perc[,13], rev(results_month[[2]]$upper_coverage_perc[,13])),
        col=rgb(0, 0, 0,0.1), border=FALSE)
abline(h=.95, lty = 2)
mtext(lets[4,2], side=3, at=at_lets, cex=cex_lets, line=1, font=2)

dev.off()

#######
# Dow|D
#######

aa = c(results_dow[[1]]$lower_coverage_perc[,10], results_dow[[1]]$upper_coverage_perc[,10],
       results_dow[[1]]$lower_coverage_perc[,11], results_dow[[1]]$upper_coverage_perc[,11],
       results_dow[[1]]$lower_coverage_perc[,12], results_dow[[1]]$upper_coverage_perc[,12],
       results_dow[[1]]$lower_coverage_perc[,13], results_dow[[1]]$upper_coverage_perc[,13],
       results_dow[[2]]$lower_coverage_perc[,10], results_dow[[2]]$upper_coverage_perc[,10],
       results_dow[[2]]$lower_coverage_perc[,11], results_dow[[2]]$upper_coverage_perc[,11],
       results_dow[[2]]$lower_coverage_perc[,12], results_dow[[2]]$upper_coverage_perc[,12],
       results_dow[[2]]$lower_coverage_perc[,13], results_dow[[2]]$upper_coverage_perc[,13])

ylim1 = c( min(aa, na.rm=T), max(aa, na.rm=T))
ylim1 = c( min(aa, na.rm=T), 1)

lets = matrix(LETTERS[1:8], byrow=T, ncol=2)
at_lets = -10
cex_lets = 1.2

titles = c("Mort, Dow|D, 10 years", "Hosp, Dow|D, 10 years", 
           "Mort, Dow|D, 15 years", "Hosp, Dow|D, 15 years",
           "Mort, Dow|D, 20 years", "Hosp, Dow|D, 20 years", 
           "Mort, Dow|D, 25 years", "Hosp, Dow|D, 25 years")

ylab1 = "Coverage (%)"

pdf(file = "Figure Coverage E-R percentiles DowD.pdf",
    width = 8.2*2/3, height = 8.2 + 8.2/3)


par(mfrow = c(4,2))


plot(percs_mort, results_dow[[1]]$coverage_perc[,10], type="l", lwd=3, ylim = ylim1, ylab = ylab1,
     xlab = "Percentile", main = titles[1])
polygon(c(percs_mort, rev(percs_mort)), 
        c(results_dow[[1]]$lower_coverage_perc[,10], rev(results_dow[[1]]$upper_coverage_perc[,10])),
        col=rgb(0, 0, 0,0.1), border=FALSE)
abline(h=0.95, lty = 2)
mtext(lets[1,1], side=3, at=at_lets, cex=cex_lets, line=1, font=2)

plot(percs_hosp, results_dow[[2]]$coverage_perc[,10], type="l", lwd=3, ylim = ylim1, ylab = ylab1,
     xlab = "Percentile", main = titles[2])
polygon(c(percs_hosp, rev(percs_hosp)), 
        c(results_dow[[2]]$lower_coverage_perc[,10], rev(results_dow[[2]]$upper_coverage_perc[,10])),
        col=rgb(0, 0, 0,0.1), border=FALSE)
abline(h=0.95, lty = 2)
mtext(lets[1,2], side=3, at=at_lets, cex=cex_lets, line=1, font=2)

plot(percs_mort, results_dow[[1]]$coverage_perc[,11], type="l", lwd=3, ylim = ylim1, ylab = ylab1,
     xlab = "Percentile", main = titles[3])
polygon(c(percs_mort, rev(percs_mort)), 
        c(results_dow[[1]]$lower_coverage_perc[,11], rev(results_dow[[1]]$upper_coverage_perc[,11])),
        col=rgb(0, 0, 0,0.1), border=FALSE)
abline(h=0.95, lty = 2)
mtext(lets[2,1], side=3, at=at_lets, cex=cex_lets, line=1, font=2)

plot(percs_hosp, results_dow[[2]]$coverage_perc[,11], type="l", lwd=3, ylim = ylim1, ylab = ylab1,
     xlab = "Percentile", main = titles[4])
polygon(c(percs_hosp, rev(percs_hosp)), 
        c(results_dow[[2]]$lower_coverage_perc[,11], rev(results_dow[[2]]$upper_coverage_perc[,11])),
        col=rgb(0, 0, 0,0.1), border=FALSE)
abline(h=0.95, lty = 2)
mtext(lets[2,2], side=3, at=at_lets, cex=cex_lets, line=1, font=2)


plot(percs_mort, results_dow[[1]]$coverage_perc[,12], type="l", lwd=3, ylim = ylim1, ylab = ylab1,
     xlab = "Percentile", main = titles[5])
polygon(c(percs_mort, rev(percs_mort)), 
        c(results_dow[[1]]$lower_coverage_perc[,12], rev(results_dow[[1]]$upper_coverage_perc[,12])),
        col=rgb(0, 0, 0,0.1), border=FALSE)
abline(h=0.95, lty = 2)
mtext(lets[3,1], side=3, at=at_lets, cex=cex_lets, line=1, font=2)

plot(percs_hosp, results_dow[[2]]$coverage_perc[,12], type="l", lwd=3, ylim = ylim1, ylab = ylab1,
     xlab = "Percentile", main = titles[6])
polygon(c(percs_hosp, rev(percs_hosp)), 
        c(results_dow[[2]]$lower_coverage_perc[,12], rev(results_dow[[2]]$upper_coverage_perc[,12])),
        col=rgb(0, 0, 0,0.1), border=FALSE)
abline(h=0.95, lty = 2)
mtext(lets[3,2], side=3, at=at_lets, cex=cex_lets, line=1, font=2)

plot(percs_mort, results_dow[[1]]$coverage_perc[,13], type="l", lwd=3, ylim = ylim1, ylab = ylab1,
     xlab = "Percentile", main = titles[7])
polygon(c(percs_mort, rev(percs_mort)), 
        c(results_dow[[1]]$lower_coverage_perc[,13], rev(results_dow[[1]]$upper_coverage_perc[,13])),
        col=rgb(0, 0, 0,0.1), border=FALSE)
abline(h=0.95, lty = 2)
mtext(lets[4,1], side=3, at=at_lets, cex=cex_lets, line=1, font=2)

plot(percs_hosp, results_dow[[2]]$coverage_perc[,13], type="l", lwd=3, ylim = ylim1, ylab = ylab1,
     xlab = "Percentile", main = titles[8])
polygon(c(percs_hosp, rev(percs_hosp)), 
        c(results_dow[[2]]$lower_coverage_perc[,13], rev(results_dow[[2]]$upper_coverage_perc[,13])),
        col=rgb(0, 0, 0,0.1), border=FALSE)
abline(h=.95, lty = 2)
mtext(lets[4,2], side=3, at=at_lets, cex=cex_lets, line=1, font=2)

dev.off()





################################################################################
# Figure bias-mse-coverage for exposure response function - Number of cases/10 #
################################################################################

# Excluding nyears = 1
##########################

# Add a column with power
##########################

aa = c(results_daily10[[1]]$lower_bias_cum[-1], results_daily10[[1]]$upper_bias_cum[-1],
       results_week10[[1]]$lower_bias_cum[-1], results_week10[[1]]$upper_bias_cum[-1],
       results_month10[[1]]$lower_bias_cum[-1], results_month10[[1]]$upper_bias_cum[-1],
       results_dow10[[1]]$lower_bias_cum[-1], results_dow10[[1]]$upper_bias_cum[-1])
ylim_a = c( min(aa, na.rm=T), max(aa, na.rm=T))

aa = c(results_daily10[[1]]$rmse_cum[-1], results_week10[[1]]$rmse_cum[-1], 
       results_month10[[1]]$rmse_cum[-1], results_dow10[[1]]$rmse_cum[-1])
ylim_b = c( 0, max(aa, na.rm=T))

aa = c(results_daily10[[1]]$lower_coverage_cum[-1], results_week10[[1]]$lower_coverage_cum[-1], 
       results_month10[[1]]$lower_coverage_cum[-1], results_dow10[[1]]$lower_coverage_cum[-1])
ylim_c = c( min(aa, na.rm=T), 1)

aa = c(results_daily10[[1]]$power01[-1], results_week10[[1]]$power01[-1], 
       results_month10[[1]]$power01[-1], results_dow10[[1]]$power01[-1],
       results_daily10[[1]]$power99[-1], results_week10[[1]]$power99[-1], 
       results_month10[[1]]$power99[-1], results_dow10[[1]]$power99[-1])
ylim_d = c( min(aa, na.rm=T), 1)


aa = c(results_daily10[[2]]$lower_bias_cum[-1], results_daily10[[2]]$upper_bias_cum[-1],
       results_week10[[2]]$lower_bias_cum[-1], results_week10[[2]]$upper_bias_cum[-1],
       results_month10[[2]]$lower_bias_cum[-1], results_month10[[2]]$upper_bias_cum[-1],
       results_dow10[[2]]$lower_bias_cum[-1], results_dow10[[2]]$upper_bias_cum[-1])
ylim_e = c( min(aa, na.rm=T), max(aa, na.rm=T))

aa = c(results_daily10[[2]]$rmse_cum[-1], results_week10[[2]]$rmse_cum[-1],
       results_month10[[2]]$rmse_cum[-1], results_dow10[[2]]$rmse_cum[-1])
ylim_f = c( 0, max(aa, na.rm=T))

aa = c(results_daily10[[2]]$lower_coverage_cum[-1], results_week10[[2]]$lower_coverage_cum[-1], 
       results_month10[[2]]$lower_coverage_cum[-1], results_dow10[[2]]$lower_coverage_cum[-1])
ylim_g = c( min(aa, na.rm=T), 1)

aa = c(results_daily10[[2]]$power01[-1], results_week10[[2]]$power01[-1], 
       results_month10[[2]]$power01[-1], results_dow10[[2]]$power01[-1],
       results_daily10[[2]]$power99[-1], results_week10[[2]]$power99[-1], 
       results_month10[[2]]$power99[-1], results_dow10[[2]]$power99[-1])
ylim_h = c( min(aa, na.rm=T), 1)



aa = c(results_daily10[[3]]$lower_bias_cum[-1], results_daily10[[3]]$upper_bias_cum[-1],
       results_week10[[3]]$lower_bias_cum[-1], results_week10[[3]]$upper_bias_cum[-1],
       results_month10[[3]]$lower_bias_cum[-1], results_month10[[3]]$upper_bias_cum[-1],
       results_dow10[[3]]$lower_bias_cum[-1], results_dow10[[3]]$upper_bias_cum[-1])
ylim_i = c( min(aa, na.rm=T), max(aa, na.rm=T))

aa = c(results_daily10[[3]]$rmse_cum[-1], results_week10[[3]]$rmse_cum[-1], results_month10[[3]]$rmse_cum[-1],
       results_dow10[[3]]$rmse_cum[-1])
ylim_j = c( 0, max(aa, na.rm=T))

aa = c(results_daily10[[3]]$lower_coverage_cum[-1], results_week10[[3]]$lower_coverage_cum[-1], 
       results_month10[[3]]$lower_coverage_cum[-1], results_dow10[[3]]$lower_coverage_cum[-1])
ylim_k = c( min(aa, na.rm=T), 1)

aa = c(results_daily10[[3]]$power[-1], results_week10[[3]]$power[-1], 
       results_month10[[3]]$power[-1], results_dow10[[3]]$power[-1])
ylim_l = c( min(aa, na.rm=T), 1)


aa = c(results_daily10[[4]]$lower_bias_cum[-1], results_daily10[[4]]$upper_bias_cum[-1],
       results_week10[[4]]$lower_bias_cum[-1], results_week10[[4]]$upper_bias_cum[-1],
       results_month10[[4]]$lower_bias_cum[-1], results_month10[[4]]$upper_bias_cum[-1],
       results_dow10[[4]]$lower_bias_cum[-1], results_dow10[[4]]$upper_bias_cum[-1])
ylim_m = c( min(aa, na.rm=T), max(aa, na.rm=T))

aa = c(results_daily10[[4]]$rmse_cum[-1], results_week10[[4]]$rmse_cum[-1],
       results_month10[[4]]$rmse_cum[-1], results_dow10[[4]]$rmse_cum[-1])
ylim_n = c( 0, max(aa, na.rm=T))

aa = c(results_daily10[[4]]$lower_coverage_cum[-1], results_week10[[4]]$lower_coverage_cum[-1], 
       results_month10[[4]]$lower_coverage_cum[-1], results_dow10[[4]]$lower_coverage_cum[-1])
ylim_o = c( min(aa, na.rm=T), 1)

aa = c(results_daily10[[4]]$power[-1], results_week10[[4]]$power[-1], 
       results_month10[[4]]$power[-1], results_dow10[[4]]$power[-1])
ylim_p = c( min(aa, na.rm=T), 1)



ylim1 = matrix( c(ylim_a, ylim_e, ylim_i, ylim_m),
                byrow=T, ncol=2)

ylim2 = matrix( c(ylim_b, ylim_f, ylim_j, ylim_n),
                byrow=T, ncol=2)

ylim3 = matrix( c(ylim_c, ylim_g, ylim_k, ylim_o),
                byrow=T, ncol=2)

ylim4 = matrix( c(ylim_d, ylim_h, ylim_l, ylim_p),
                byrow=T, ncol=2)

lets = matrix(LETTERS[1:16], byrow=T, ncol=4)
at_lets = -6
cex_lets = 1.2

titles = c("Mortality | Temperature", "Hospit. | Temperature", 
           expression('Mortality | NO'[2]), expression('Hospit. | NO'[2]))

imax = c(13-1, 13-1, 13-1, 13-1)


indyears = c(2:10, 15, 20, 25)

pdf(file = "Figure bias-mse-coverage-power E-R no1year 0101_1pct new temp Ndiv10.pdf",
    width = 8.2, height = 8.2 + 8.2/3)

par(mfrow=c(4,4))

for (i in 1:4) {
  
  # Bias
  
  plot(indyears[1:imax[i]], results_daily10[[i]]$bias_cum[-1], type="l", ylab="Bias", 
       xlab="Number of years", lwd=5, cex.lab=1.2, ylim=ylim1[i,], xlim=c(1,25))
  polygon(c(indyears[1:imax[i]], rev(indyears[1:imax[i]])),
          c(results_daily10[[i]]$lower_bias_cum[-1], 
            rev(results_daily10[[i]]$upper_bias_cum[-1])),
          col=rgb(0, 0, 0,0.3), border=FALSE, lwd=1)
  lines(indyears[1:imax[i]], results_week10[[i]]$bias_cum[-1], col="red", lwd=5)
  polygon(c(indyears[1:imax[i]], rev(indyears[1:imax[i]])),
          c(results_week10[[i]]$lower_bias_cum[-1], 
            rev(results_week10[[i]]$upper_bias_cum[-1])),
          col=rgb(1, 0, 0,0.3), border=FALSE, lwd=1)
  lines(indyears[1:imax[i]], results_dow10[[i]]$bias_cum[-1], col="#44AA99", lwd=5)
  polygon(c(indyears[1:imax[i]], rev(indyears[1:imax[i]])),
          c(results_dow10[[i]]$lower_bias_cum[-1], 
            rev(results_dow10[[i]]$upper_bias_cum[-1])),
          col=rgb(0.27, 0.67, 0.6, 0.3), border=FALSE, lwd=1)
  
  lines(indyears[1:imax[i]], results_month10[[i]]$bias_cum[-1], col="blue", lwd=5)
  polygon(c(indyears[1:imax[i]], rev(indyears[1:imax[i]])),
          c(results_month10[[i]]$lower_bias_cum[-1], 
            rev(results_month10[[i]]$upper_bias_cum[-1])),
          col=rgb(0, 0, 1,0.3), border=FALSE, lwd=1)
  
  abline(h=0,lty=1)
  mtext(titles[i], side=3, line=0.5, cex=.8, font=2)
  if (i==1) {
    legend("bottomright",c("D|D", "W|D","Dow|D", "M|D"), col=c("black","red", "#44AA99", "blue"), lty=1,
           lwd=5, cex=1.2, bty="n")
  }
  mtext(lets[i,1], side=3, at=at_lets, cex=cex_lets, line=1, font=2)
  
  
  # RMSE
  
  plot(indyears[1:imax[i]], results_daily10[[i]]$rmse_cum[-1], type="l", ylab="RMSE", 
       xlab="Number of years",  lwd=5, cex.lab=1.2, ylim=ylim2[i,], xlim=c(1,25))
  lines(indyears[1:imax[i]], results_week10[[i]]$rmse_cum[-1], col="red", lwd=5)
  lines(indyears[1:imax[i]], results_dow10[[i]]$rmse_cum[-1], col="#44AA99", lwd=5)
  lines(indyears[1:imax[i]], results_month10[[i]]$rmse_cum[-1], col="blue", lwd=5)
  
  abline(h=0,lty=1)
  #legend("topright",c("D|D", "W|D","Dow|D", "M|D"), col=c("black","red", "#44AA99", "blue"), lty=1,
  #       lwd=5, cex=1.2)
  
  mtext(titles[i], side=3, line=0.5, cex=.8, font=2)
  mtext(lets[i,2], side=3, at=at_lets, cex=cex_lets, line=1, font=2)
  
  
  # Coverage
  
  plot(indyears[1:imax[i]], results_daily10[[i]]$coverage_cum[-1], type="l", ylab="Coverage", 
       xlab="Number of years", lwd=5, cex.lab=1.2, ylim=ylim3[i,], xlim=c(1,25))
  lines(indyears[1:imax[i]], results_week10[[i]]$coverage_cum[-1], col="red", lwd=5)
  lines(indyears[1:imax[i]], results_dow10[[i]]$coverage_cum[-1], col="#44AA99", lwd=5)
  polygon(c(indyears[1:imax[i]], rev(indyears[1:imax[i]])),
          c(results_week10[[i]]$lower_coverage_cum[-1], 
            rev(results_week10[[i]]$upper_coverage_cum[-1])),
          col=rgb(1, 0, 0,0.3), border=FALSE, lwd=1)
  polygon(c(indyears[1:imax[i]], rev(indyears[1:imax[i]])),
          c(results_dow10[[i]]$lower_coverage_cum[-1], 
            rev(results_dow10[[i]]$upper_coverage_cum[-1])),
          col=rgb(0.27, 0.67, 0.6,0.3), border=FALSE, lwd=1)
  lines(indyears[1:imax[i]], results_month10[[i]]$coverage_cum[-1], col="blue", lwd=5)
  polygon(c(indyears[1:imax[i]], rev(indyears[1:imax[i]])),
          c(results_month10[[i]]$lower_coverage_cum[-1], 
            rev(results_month10[[i]]$upper_coverage_cum[-1])),
          col=rgb(0, 0, 1,0.3), border=FALSE, lwd=1)
  polygon(c(indyears[1:imax[i]], rev(indyears[1:imax[i]])),
          c(results_daily10[[i]]$lower_coverage_cum[-1], 
            rev(results_daily10[[i]]$upper_coverage_cum[-1])),
          col=rgb(0, 0, 0,0.3), border=FALSE, lwd=1)
  abline(h=0.95,lty=1)
  #legend("bottomright",c("D|D", "W|D","Dow|D", "M|D"), col=c("black","red", "#44AA99", "blue"), lty=1,
  #       lwd=5, cex=1.2)
  mtext(titles[i], side=3, line=0.5, cex=.8, font=2)
  mtext(lets[i,3], side=3, at=at_lets, cex=cex_lets, line=1, font=2)
  
  # Power
  
  if (i==1 | i==2) {
    plot(indyears[1:imax[i]], results_daily10[[i]]$power99[-1], type="l", ylab="Power", 
         xlab="Number of years", lwd=5, cex.lab=1.2, ylim=ylim4[i,], xlim=c(1,25))
    lines(indyears[1:imax[i]], results_week10[[i]]$power99[-1], col="red", lwd=5)
    lines(indyears[1:imax[i]], results_dow10[[i]]$power99[-1], col="#44AA99", lwd=5)
    lines(indyears[1:imax[i]], results_month10[[i]]$power99[-1], col="blue", lwd=5)
    
    lines(indyears[1:imax[i]], results_daily10[[i]]$power01[-1], col="black", lwd=2)
    lines(indyears[1:imax[i]], results_week10[[i]]$power01[-1], col="red", lwd=2)
    lines(indyears[1:imax[i]], results_dow10[[i]]$power01[-1], col="#44AA99", lwd=2)
    lines(indyears[1:imax[i]], results_month10[[i]]$power01[-1], col="blue", lwd=2)
    
    #abline(h=0.95,lty=1)
    mtext(titles[i], side=3, line=0.5, cex=.8, font=2)
    mtext(lets[i,4], side=3, at=at_lets, cex=cex_lets, line=1, font=2)
  } else if (i==3 | i==4) {
    plot(indyears[1:imax[i]], results_daily10[[i]]$power[-1], type="l", ylab="Power", 
         xlab="Number of years", lwd=5, cex.lab=1.2, ylim=ylim4[i,], xlim=c(1,25))
    lines(indyears[1:imax[i]], results_week10[[i]]$power[-1], col="red", lwd=5)
    lines(indyears[1:imax[i]], results_dow10[[i]]$power[-1], col="#44AA99", lwd=5)
    lines(indyears[1:imax[i]], results_month10[[i]]$power[-1], col="blue", lwd=5)
    #abline(h=0.95,lty=1)
    mtext(titles[i], side=3, line=0.5, cex=.8, font=2)
    mtext(lets[i,4], side=3, at=at_lets, cex=cex_lets, line=1, font=2)
    
  }
}

dev.off()


