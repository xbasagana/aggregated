# Written by Xavier Basagaña
# Version 0.27
#  2024/10/10

# Declaring Libraries
suppressMessages( library(dlnm) );
suppressMessages( library(splines) );
suppressMessages( library(tsModel))
suppressMessages( library(FluMoDL))  # attrdl
suppressMessages( library(data.table))
suppressMessages( library(Rfast))
suppressMessages( library(maxLik))
suppressMessages(library(matrixcalc))

# Function returning the log likelihood

agr_logLik <- function(beta, X, y, period) {
  # Auxiliary stuff
  sum_ex <- group(exp(X %*% beta), period, method="sum")
  
  sum(y * log(sum_ex) - sum_ex)   
}

# Function returning the score function 

d1score_ind <- function(beta, X, y, period) {
  
  ex <- as.vector(exp(X %*% beta))
  # Auxiliary stuff
  sum_ex <- group(exp(X %*% beta), period, method="sum")
  
  sum_xex <- rowsum(X*ex, period)
  
  # Individual score vectors
  d1i <- sum_xex*(y/sum_ex) - sum_xex
  
  return(d1i)
}

# Function returning the gradient

d1func_aggr <- function(beta, X, y, period) {
  
  # Auxiliary stuff
  ex <- as.vector(exp(X %*% beta))
  
  sum_ex <- group(exp(X %*% beta), period, method="sum")
  sum_xex <- rowsum(X*ex, period)
  
  # Gradient
  d1 <- crossprod(y, sum_xex/sum_ex) - colsums(sum_xex)
  
  return(d1)
}

# Function returning the information function

d2func_aggr <- function(beta, X, y, period) {
  periods <- unique(period)
  ex <- as.vector(exp(X %*% beta))
  sum_ex <- group(exp(X %*% beta), period, method="sum")
  sum_xex <- rowsum(X*ex, period)
  hessi <- matrix(0,nrow=ncol(X),ncol=ncol(X))
  
  for (i in 1:length(periods)) {
    xi <- X[period==periods[i],]
    if (is.vector(xi)) xi <- matrix(xi,nrow=1)
    xxex <- array(NA,dim=c(ncol(X), ncol(X), nrow(xi) ))
    for (j in 1:nrow(xi)) {
      xxex[,,j] <- tcrossprod(xi[j,],xi[j,]) * ex[period==periods[i]][j]
    }
    sum_xxexi <- rowSums(xxex,dims=2)
    if ((i%% 1000)==0) {
      rm(xxex)
      gc()
    }
    hessi <- hessi + (sum_xxexi * sum_ex[i] - tcrossprod(sum_xex[i,],sum_xex[i,])) * (y[i]/(sum_ex[i]^2)) - sum_xxexi
  }
  return(hessi)
}

return_null_model <- function(error_code, iter, logLik) {
  res <- list(coef = NULL, vcov = NULL, mle = NULL, converged = FALSE,
              X = NULL, period = NULL, disp = NULL, iter = iter,
              logLik = logLik,
              AIC = NULL,
              family = list(family=family, link="log", linkfun=function (mu) log(mu)),
              pred_day_NA = NULL, error_code = error_code)
  class(res) <- "modagr"
  return(res)
}

# Function to fit the model, using aggregated mortality and daily temperature

fit_aggregate <- function(Y, X, CB = NULL, name_exposure = NULL, 
                          formula, family = "quasipoisson", tol = 1e-8, maxit = 500,
                          start=NULL, seed = NULL, ntry=3) {
  # 'Y' is a data.frame with two variables: the period indicator and the aggregated 
  #     counts for that period (in this order)
  # 'X' is a data.frame with the exposures and covariates at daily level, and the period indicator
  # The period variable should have the same name in both Y and X 
  # 'formula' is the formula to be used if daily data were available. The name of the outcome
  #     variable in the formula should be the same than the second column in Y
  # 'CB' is a crossbasis object to be used as predictor in the models. If 'name_exposure'
  #      is provided, then 'CB' should be NULL
  # 'name_exposure' is the name of the main exposure variable, included in 'X'. If
  #     'CB' is provided, then 'name_exposure' should be NULL
  
  # Return error_code:
  #   -999: number of parameters > number of observations
  #   -888: family is not correct
  #   -777: covariance is positive semidefinite - more than one solution
  #   -666: covariance is not positive semidefinite 
  #   -555: covariance has infinite values
  
  
  error_code <- 0
  if (family!="quasipoisson" & family!="poisson") {
    cat("Family can only be 'poisson' or 'quasipoisson' \n")
    error_code <- -888
  }
  
  if (error_code == 0) {  
    name_period <- names(Y)[1]
    name_response <- names(Y)[2]
    tt <- table(X[,name_period])
    nperiods <- data.frame(x1 = as.numeric(names(tt)), ndays = as.numeric(tt))
    names(nperiods)[1] = name_period
    
    if (!is.null(seed)) set.seed(seed)
    
    # Distribute the counts randomly across the days, then fit a daily model to obtain initial values
    Yaux <- list()
    for (i in 1:nrow(Y)) {
      if (is.na(Y[i,2])) {
        Yaux[[i]] <- rep(NA, nperiods$ndays[i])
      } else {
        Yaux[[i]] <- rmultinom(1, Y[i,2], prob = rep(1/nperiods$ndays[i], nperiods$ndays[i]))
      }
    }
    dat <- data.frame(yaux = unlist(Yaux), X)
    names(dat)[which(names(dat)=="yaux")] <- name_response
    
    model_ini <- glm( formula, data=dat, family = family, na.action = "na.exclude" )
    
    # Extract model matrix (matrix of predictors) from daily model 
    X <- model.matrix(model_ini)
    # The same but including missing values
    X_NA <- model.matrix.lm(formula, data=dat, na.action = "na.pass")
    
    if (ncol(X) > nrow(Y)) {
      cat("Error: Number of parameters is greater than number of data points \n")
      error_code <- -999
    }
    
    if (error_code == 0) {
      
      # Check we have complete periods after excluding missing values
      
      dat2 <- dat[, c(name_period, name_response)]
      aa <- merge(X, dat2, by = 'row.names')
      aa <- aa[order(as.numeric(aa$Row.names)),]
      tt2 <- table(aa[name_period])
      nperiods2 <- data.frame(x1 = as.numeric(names(tt2)), ndays = as.numeric(tt2))
      names(nperiods2)[1] <- name_period
      bb <- merge(nperiods, nperiods2, by = name_period)
      
      if (min(bb$ndays.x == bb$ndays.y) == 0) {
        cat("After removing missing values, some periods are not complete \n")
        to_exclude <- bb[which(bb$ndays.x != bb$ndays.y), 1]
        cat(paste0("Removing periods: ", to_exclude,"\n"))
        # set response variable to missing for incomplete periods
        dat[,name_response][(dat[,name_period] %in% to_exclude)] <- NA
        # set exposure variable to missing for incomplete periods
        if (!is.null(CB)) {
          CB[(dat[,name_period] %in% to_exclude),1] <- NA
        } else if (!is.null(name_exposure)) {
          dat[(dat[,name_period] %in% to_exclude) , name_exposure] <- NA
        }
        # refit initial model with dataset with incomplete periods as NA
        model_ini <- glm( formula, data=dat, family = family, na.action = "na.exclude" )
        X <- model.matrix(model_ini)
        dat2 <- dat[, c(name_period, name_response)]
        aa <- merge(X, dat2, by = 'row.names')
        aa <- aa[order(as.numeric(aa$Row.names)), ]
        y <- Y[!(Y[,1] %in%  to_exclude), 2]
      }
      
      # Keep outcome data only for the periods without missing data
      
      y <- Y[Y[,1] %in% as.numeric(names(table(aa[,name_period]))),2]
      
      if (ncol(X) >= length(y)) {
        cat("Error: Number of parameters is equal or greater than number of data points (after excluding missing values) \n")
        error_code <- -999
        res <- return_null_model(error_code, iter=0, logLik=NULL)
      } else {
        # initial values
        if (is.null(start)) beta <- coef(model_ini) else beta <- start
        
        mleBHHH <- maxLik(logLik = agr_logLik, grad = d1score_ind, 
                          start = c(beta), method = "BHHH", 
                          control = list(tol=tol, iterlim = maxit),
                          X = X, y = y, period = aa[, name_period])
        
        if (is.infinite(max(vcov(mleBHHH)))) {
          error_code <- -555
        }
        
        if ((mleBHHH$code %in% c(1,2,8)) & error_code==0 ) { # converged
          
          beta <- coef(mleBHHH)
          # refine with Newton-Raphson
          mle <- maxLik(logLik = agr_logLik, grad = d1func_aggr,
                        hess = d2func_aggr,
                        start = c(beta), 
                        control=list(tol=tol,iterlim=maxit),
                        X=X, y=y, period=aa[,name_period])
          
          beta <- coef(mle)
          V <- vcov(mle)
          
          if (is.infinite(max(V))) {
            error_code <- -555
          }
          
          if (error_code == 0) {
            # Try to find convergence if V is not positive definite. Try try_iter times
            try_iter <- 1
            
            while ( error_code == 0 & (try_iter<=ntry))  {
              # refit (different random seed can fix stopping at a saddle point)
              
              if (is.null(start)) {
                
                # Distribute the counts randomly across the days, then fit a daily model to obtain initial values
                Yaux <- list()
                for (i in 1:nrow(Y)) {
                  if (is.na(Y[i,2])) {
                    Yaux[[i]] <- rep(NA, nperiods$ndays[i])
                  } else {
                    Yaux[[i]] <- rmultinom(1, Y[i,2], prob = rep(1/nperiods$ndays[i], nperiods$ndays[i]))
                  }
                }
                dat$yaux <- unlist(Yaux)
                names(dat)[which(names(dat)=="yaux")] <- name_response
                
                model_ini <- glm( formula, data=dat, family = family, na.action = "na.exclude" )
                beta <- coef(model_ini)            
                
                mle <- maxLik(logLik = agr_logLik, grad = d1func_aggr,
                              hess = d2func_aggr,
                              start = c(beta), 
                              control=list(tol=tol,iterlim=maxit),
                              X=X, y=y, period=aa[,name_period])
                
                beta <- coef(mle)
                V <- vcov(mle)

                if (is.infinite(max(V))) {
                  error_code <- -555
                } else {
                  if (!is.positive.definite(V)) {
                    if (is.positive.semi.definite(V)) (error_code <- -777) else (error_code <- -666)
                  }
                }

              }
              try_iter <- try_iter + 1
              
            }
            
            if (is.infinite(max(V))) {
              error_code <- -555
            }
            
            if (error_code==0) {
              # daily prediction
              pred_day <- exp(X %*% beta)
              
              # prediction for period
              pred_period <- aggregate(pred_day, by=list(period = aa[, name_period]),FUN=sum)$V1
              pearson_resid <- (y - pred_period)/sqrt(pred_period)
              
              pred_day_NA <- exp(X_NA %*% beta)
              
              if (family=="quasipoisson") {
                # Dispersion parameter
                dp <- sum(pearson_resid^2)/(length(y) - length(beta))
              } else {
                dp <- 1
              }
              
              # Covariance with dispersion
              Vdisp <- V*dp
              
              if (!is.positive.definite(Vdisp)) {
                if (is.positive.semi.definite(Vdisp)) (error_code <- -777) else (error_code <- -666)
              }
              
              converged <- ((mle$code %in% c(1,2,8)) & is.positive.semi.definite(Vdisp) & error_code==0)
              
              res <- list(coef=beta, vcov=Vdisp, mle=mle, converged=converged,
                          X=X, period=unlist(aa[name_period]), disp=dp, iter = mleBHHH$iterations + mle$iterations,
                          logLik = maxValue(mle), 
                          AIC = -2*maxValue(mle)/dp + 2*(length(beta)+1),
                          family = list(family=family,link="log",linkfun=function (mu) log(mu)),
                          pred_day_NA = pred_day_NA, error_code = error_code)
              class(res) <- "modagr"
            } else {
              res <- return_null_model(error_code, iter=mleBHHH$iterations, logLik = maxValue(mleBHHH))
            }
          } else {
            res <- return_null_model(error_code, iter=mleBHHH$iterations, logLik = maxValue(mleBHHH))
          }
        } else {
          res <- return_null_model(error_code, iter=mleBHHH$iterations, logLik = maxValue(mleBHHH))
        }
        
      }
      
    } else {
      res <- return_null_model(error_code, iter=0, logLik=NULL)
    }
  } else {
    res <- return_null_model(error_code, iter=0, logLik=NULL)
  }
  
  return(res)
}

coef.modagr <- function(x,...) {
  x$coef
}

vcov.modagr <- function(x,...) {
  x$vcov
}

model.matrix.modagr <- function(x,...) {
  x$X
}

predict.modagr <- function(mod,type="response") {
  if (type == "response") mod$pred_day_NA else if (type == "link") log(mod$pred_day_NA)
}

summary.modagr <- function(mod) {
  summary(mod$mle)
}

# Function to add RRfit, RRlow and RRhigh to a crossreduce object

add_exp <- function(cr) {
  cr$RRfit <- exp(cr$fit)
  z <- qnorm(1 - (1 - cr$ci.level)/2)
  cr$RRlow <- exp(cr$fit - z * cr$se) 
  cr$RRhigh <- exp(cr$fit + z * cr$se) 
  cr$model.link = "log"
  return(cr)
}

attrdl.aggr <- function (x, basis, cases, model = NULL, coef = NULL, vcov = NULL, 
                         type = "af", dir = "back", tot = TRUE, cen, range = NULL, 
                         sim = FALSE, nsim = 5000, sub = 1:length(cases)) 
{
  # Code taken from attrdl in library FluMoDL, modified to accomodate 
  #    objects of class modagr
  
  .getcoef <- getFromNamespace("getcoef", "dlnm")
  .getvcov <- getFromNamespace("getvcov", "dlnm")
  .getlink <- getFromNamespace("getlink", "dlnm")
  .seqlag <- getFromNamespace("seqlag", "dlnm")
  .mkXpred <- getFromNamespace("mkXpred", "dlnm")
  if (packageVersion("dlnm") < "2.2.0") 
    stop("update dlnm package to version >= 2.2.0")
  name <- deparse(substitute(basis))
  type <- match.arg(type, c("an", "af"))
  dir <- match.arg(dir, c("back", "forw"))
  if (missing(cen) && is.null(cen <- attr(basis, "argvar")$cen)) 
    stop("'cen' must be provided")
  if (!is.numeric(cen) && length(cen) > 1L) 
    stop("'cen' must be a numeric scalar")
  attributes(basis)$argvar$cen <- NULL
  if (!is.null(range)) 
    x[x < range[1] | x > range[2]] <- cen
  lag <- attr(basis, "lag")
  if (NCOL(x) == 1L) {
    at <- if (dir == "back") 
      tsModel::Lag(x, seq(lag[1], lag[2]))
    else matrix(rep(x, diff(lag) + 1), length(x))
  }
  else {
    if (dir == "forw") 
      stop("'x' must be a vector when dir='forw'")
    if (ncol(at <- x) != diff(lag) + 1) 
      stop("dimension of 'x' not compatible with 'basis'")
  }
  if (NROW(cases) != NROW(at)) 
    stop("'x' and 'cases' not consistent")
  if (NCOL(cases) > 1L) {
    if (dir == "back") 
      stop("'cases' must be a vector if dir='back'")
    if (ncol(cases) != diff(lag) + 1) 
      stop("dimension of 'cases' not compatible")
    den <- sum(rowMeans(cases, na.rm = TRUE), na.rm = TRUE)
    cases <- rowMeans(cases)
  }
  else {
    den <- sum(cases[sub], na.rm = TRUE)
    if (dir == "forw") 
      cases <- rowMeans(as.matrix(tsModel::Lag(cases, -seq(lag[1], 
                                                           lag[2]))))
  }
  if (!is.null(model)) {
    cond <- paste0(name, "[[:print:]]*v[0-9]{1,2}\\.l[0-9]{1,2}")
    if (ncol(basis) == 1L) 
      cond <- name
    model.class <- class(model)
    coef <- .getcoef(model, model.class)
    ind <- grep(cond, names(coef))
    coef <- coef[ind]
    vcov <- .getvcov(model, model.class)[ind, ind, drop = FALSE]
    #model.link <- .getlink(model, model.class)
    #if (model.link != "log") 
    #  stop("'model' must have a log link function")
    model.link = "log"
  }
  typebasis <- ifelse(length(coef) != ncol(basis), "one", 
                      "cb")
  predvar <- if (typebasis == "one") 
    x
  else seq(NROW(at))
  predlag <- if (typebasis == "one") 
    0
  else .seqlag(lag)
  if (typebasis == "cb") {
    Xpred <- .mkXpred(typebasis, basis, at, predvar, predlag, 
                      cen)
    Xpredall <- 0
    for (i in seq(length(predlag))) {
      ind <- seq(length(predvar)) + length(predvar) * (i - 
                                                         1)
      Xpredall <- Xpredall + Xpred[ind, , drop = FALSE]
    }
  }
  else {
    basis <- do.call(onebasis, c(list(x = x), attr(basis, 
                                                   "argvar")))
    Xpredall <- .mkXpred(typebasis, basis, x, predvar, predlag, 
                         cen)
  }
  if (length(coef) != ncol(Xpredall)) 
    stop("arguments 'basis' do not match 'model' or 'coef'-'vcov'")
  if (any(dim(vcov) != c(length(coef), length(coef)))) 
    stop("arguments 'coef' and 'vcov' do no match")
  if (typebasis == "one" && dir == "back") 
    stop("only dir='forw' allowed for reduced estimates")
  af <- 1 - exp(-drop(as.matrix(Xpredall %*% coef)))
  an <- (af * cases)
  if (tot) {
    isna <- is.na(an[sub])
    af <- sum(an[sub][!isna])/sum(cases[sub][!isna])
    an <- af * sum(cases[sub], na.rm = TRUE)
  }
  else {
    af <- af[sub]
    an <- an[sub]
  }
  if (sim) {
    k <- length(coef)
    eigen <- eigen(vcov)
    X <- matrix(rnorm(length(coef) * nsim), nsim)
    coefsim <- coef + eigen$vectors %*% diag(sqrt(eigen$values), 
                                             k) %*% t(X)
    if (tot) {
      afsim <- apply(coefsim, 2, function(coefi) {
        ani <- ((1 - exp(-drop(Xpredall %*% coefi))) * 
                  cases)
        sum(ani[!is.na(ani)])/sum(cases[!is.na(ani)])
      })
      ansim <- afsim * den
    }
    else {
      afsim <- apply(coefsim, 2, function(coefi) {
        ani <- ((1 - exp(-drop(Xpredall %*% coefi))) * 
                  cases)[sub]
        sum(ani[!is.na(ani)])/sum(cases[sub][!is.na(ani)])
      })
      ansim <- afsim * sum(cases[sub], na.rm = TRUE)
    }
  }
  res <- if (sim) {
    if (type == "an") 
      ansim
    else afsim
  }
  else {
    if (type == "an") 
      an
    else af
  }
  return(res)
}
