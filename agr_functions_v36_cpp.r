# Changes introduced in v36: 
#      - added line in 211. ifelse in get_ini for seed did not work
#      - fixing set.seed in get_ini()
#      - when V has Infinities, it returns the estimated coefficients
#      - fit_aggregate function now also returns  'y'
#      - fit_aggregate function now also returns gvcov, the generalized covariance matrix (non-NULL only when error_code == -555)
#      - gvcov_modagr function created to extract the the generalized covariance matrix
#      - summary.modagr modified, now it provides covariance matrix that integrates the dispersion parameter
#      - added summary_gcov() functions, which returns the summary results based on the generalized covariance matrix

# Declaring Libraries
suppressMessages( library(dlnm) )
suppressMessages( library(maxLik) ) # maxLik
suppressMessages(library(Rcpp))
suppressMessages(library(MASS)) # ginv

# C++ Functions
Rcpp::sourceCpp('agr_logLik_cpp.cpp')  # Loglikelihood Function
Rcpp::sourceCpp('d1score_ind_cpp.cpp') # Score Function for BHHHH
Rcpp::sourceCpp('d1func_aggr_cpp.cpp') # Gradient Function for NR
Rcpp::sourceCpp('d2func_aggr_cpp.cpp') # Hessian Function for NR

# Function that returns an empty Output Object
return_null_model <- function( error_code, iter, logLik ){
  res <- list( coef = NULL,
               vcov = NULL,
               mle = NULL,
               converged = FALSE,
               X = NULL,
               period = NULL,
               disp = NULL,
               iter = iter,
               logLik = logLik,
               QAIC = NULL,
               family = list( family = family, link = "log", linkfun = function(mu) log(mu) ),
               pred_disag_NA = NULL,
               error_code = error_code,
               gvcov = NULL)
  class(res) <- "modagr"
  return(res)
}



# Fitting a Disaggregated|Disaggregated Model to Obtain Initial Parameter Values After Randomly Distributing Counts Across Disaggregated Time Steps 
get_ini = function( Y, X, obs_in_period, name_response, formula, family, seed ){
  
  # Randomly Distributing Counts Across Disaggregated Time Steps
  Yaux <- list()
  for ( i in 1:nrow(Y) ) {
    if ( is.na( Y[i,2] ) ) { 
      Yaux[[i]] <- rep( NA, obs_in_period$nobs[i] ) 
    } else { 
      if ( !is.null(seed) ) set.seed(seed + i)
      Yaux[[i]] <- rmultinom( 1, Y[i,2], prob = rep( 1 / obs_in_period$nobs[i], obs_in_period$nobs[i] ) ) }
  }
  disag_model_data <- data.frame( yaux = unlist(Yaux), X )
  names(disag_model_data)[ which( names(disag_model_data) == "yaux" ) ] <- name_response
  
  # Fitting a Daily|Daily Model to Obtain Initial Parameter Values
  model_ini <- try( glm( formula, data = disag_model_data, family = family, na.action = "na.exclude" ), silent = FALSE )
  if ( "try-error" %in% class(model_ini) ) {
    error_code <- -444
    model_ini <- NULL
    model_ini$converged <- FALSE
  } else {
    error_code <- 0
  }
  
  return( list( error_code = error_code, model_ini = model_ini, disag_model_data = disag_model_data ) )
  
}



# Function to Fit a Temporally Aggregated Model
#   'Y' is a data.frame with two variables in order: (i) period indicator and (ii) period aggregated counts
#   'X' is a data.frame at disaggregated resolution with (i) period indicator and (ii) exposures/covariates 
#       The period indicator should have the same name in 'X' and 'Y'
#   'CB' is a crossbasis object to be used as predictor in the models. If 'name_exposure' is provided, then 'CB' should be NULL
#   'name_exposure' is the name of the main exposure variable, included in 'X'. If 'CB' is provided, then 'name_exposure' should be NULL
#   'formula' is the formula to be used if disaggregated data were available. The name of the outcome variable in the formula should be the same than the second column in 'Y' 
fit_aggregate_cpp <- function( Y,
                               X,
                               CB = NULL,
                               name_exposure = NULL,
                               formula,
                               family = "quasipoisson",
                               tol = 1e-8,
                               maxit = 5000,
                               seed = NULL,
                               ntryini = 10,
                               ntryfit = 0 ) {
  
  # Code Nomenclature for 'error_code'
  #  -9999: Invalid input parameters
  #   -999: Number of Parameters >= Number of Observations
  #   -777: Covariance Matrix Is Positive Semidefinite But Not Positive Definite, i.e. More than One Solution 
  #   -666: Covariance Matrix Is Not Positive Semidefinite 
  #   -555: Covariance Matrix Has Infinite Values
  #   -444: Non-Convergent Disaggregated|Disaggregated Model
  #   -333: Error Found in BHHH or NR Fitting
  #      0: Convergent Temporally Aggregated Model

  error_code <- 0
  
  # Checking that Y Is a Data Frame with 2 Columns: Temporal Aggregation, Temporally Aggregated Response
  if ( !is.data.frame(Y) | dim(Y)[2] != 2 ) {
    print( "fit_aggregate_cpp :: Invalid Y" )
    error_code <- -9999
  } 
  
  # Checking that X Is a Data Frame 
  if ( !is.data.frame(X) ) {
    print( "fit_aggregate_cpp :: X is not a data.frame" )
    error_code <- -9999
  }
  
  # Checking that CB Is a crossbasis 
  if (!is.null(CB) &  max(class(CB) == "crossbasis") != 1) {
    print( "fit_aggregate_cpp :: Invalid CB" )
    error_code <- -9999
  }
  
  # Checking that name_exposure is character
  if (!is.null(name_exposure) & (length(name_exposure) != 1 | !is.character(name_exposure) )) {
    print( "fit_aggregate_cpp :: Invalid name_exposure" )
    error_code <- -9999 
  }

  # Checking that one of name_exposure and CB are provided
  if (is.null(name_exposure) & is.null(CB)) {
    print( "fit_aggregate_cpp :: One of 'CB' or 'name_exposure' should be provided" )
    error_code <- -9999 
  }
  
  if (!is.null(name_exposure) & !is.null(CB)) {
    print( "fit_aggregate_cpp :: Only one of 'CB' or 'name_exposure' should be provided" )
    error_code <- -9999 
  }
  
  # Checking that formula is formula or character
  if (class(formula) != "formula" & class(formula) != "character") {
    print( "fit_aggregate_cpp :: Invalid formula" )
    error_code <- -9999 
  }
  
  if ( family != "quasipoisson" & family != "poisson" ) {
    cat( "Family can only be 'poisson' or 'quasipoisson'\n" )
    error_code <- -9999
  }
  
  # Checking tol
  if ( length(tol) != 1 | !is.numeric(tol) | tol <= 0 ) {
    print( "fit_aggregate_cpp :: Invalid tol" )
    error_code <- -9999 
  }
  
  # Checking maxit 
  if ( length(maxit) != 1 | !is.numeric(maxit) | maxit <= 0 ) {
    print( "fit_aggregate_cpp :: Invalid maxit" )
    error_code <- -9999 
  }
  
  # Checking seed
  if ( !is.null(seed) & (length(seed) != 1 | !is.numeric(seed) ) ) {
    print( "fit_aggregate_cpp :: Invalid seed" )
    error_code <- -9999 
  }
  
  # Checking ntryini
  if ( length(ntryini) != 1 | !is.numeric(ntryini) | ntryini < 0 ) {
    print( "fit_aggregate_cpp :: Invalid ntryini" )
    error_code <- -9999 
  }
  
  # Checking ntryfit
  if ( length(ntryfit) != 1 | !is.numeric(ntryfit) | ntryfit < 0 ) {
    print( "fit_aggregate_cpp :: Invalid ntryfit" )
    error_code <- -9999 
  }

  
  if ( error_code == 0 ) {
    
    name_period <- names(Y)[1] # Period Indicator Name
    name_response <- names(Y)[2] # Outcome Variable Name
    
    # If Period Indicators Are Dates, Transform Them into Numbers
    if ( class(Y[, name_period]) == "Date" ) {
      Y[, name_period] <- as.numeric(Y[, name_period])
      X[, name_period] <- as.numeric(X[, name_period])
    }
    
    # Number of Disaggregated Time Steps of Each Temporal Aggregation 
    tt <- table( X[, name_period] )
    obs_in_period_original <- data.frame( x1 = as.numeric(names(tt)), nobs = as.numeric(tt) )
    names(obs_in_period_original)[1] <- name_period
    rm(tt)
    
    # Fitting a Disaggregated|Disaggregated Model by Randomly Distributing Aggregated Counts across Disaggregated Time Steps
    ini <- get_ini( Y = Y, X = X, obs_in_period = obs_in_period_original, name_response = name_response, formula = formula, family = family, seed = seed )
    error_code <- ini$error_code
    model_ini <- ini$model_ini
    disag_model_data <- ini$disag_model_data
    
    # If the Disaggregated|Disaggregated Model Did Not Converge, Fit it Again until Convergence up to 'ntryini' Times
    try_iter <- 1
    while ( !model_ini$converged & try_iter <= ntryini ) {
      ini <- get_ini( Y = Y, X = X, obs_in_period = obs_in_period_original, name_response = name_response, formula = formula, family = family, seed = ifelse( is.null(seed), NULL, seed + max(ntryini,ntryfit) * (try_iter-1) + 1 ) )
      if ( ini$error_code != -444 ) {
        error_code <- ini$error_code
        model_ini <- ini$model_ini
        disag_model_data <- ini$disag_model_data
      }
      try_iter <- try_iter + 1
    }
    rm(try_iter)
    
  }
  
  # Extracting 'y', 'X' and 'period', Excluding Missing Values
  if ( error_code == 0 ) {
    
    # Extract Model Matrix (i.e. Predictor Matrix) from Disaggregated|Disaggregated Model (Excluding Missing Values) 
    X <- model.matrix(model_ini)
    
    # Extract Model Matrix (i.e. Predictor Matrix) from Disaggregated|Disaggregated Model (Not Excluding Missing Values) 
    X_NA <- model.matrix.lm( formula, data = disag_model_data, na.action = "na.pass" )
    
    # Outcome Dataframe (After Excluding Missing Values)
    disag_model_data_aux <- disag_model_data[ , c(name_period, name_response) ]
    
    # If Period Indicators Are Dates, Transform Them into Numbers
    if ( class( disag_model_data_aux[, name_period] ) == "Date" ) disag_model_data_aux[, name_period] = as.numeric(disag_model_data_aux[, name_period])
    
    # Outcome (disaggregated) and Covariate Dataframe - (After Excluding Missing Values) 
    disag_model_matrix <- merge( X, disag_model_data_aux, by = 'row.names' )
    disag_model_matrix <- disag_model_matrix[ order( as.numeric( disag_model_matrix$Row.names ) ), ]
    
    # Number of Disaggregated Time Steps of Each Temporal Aggregation (After Excluding Missing Values) 
    tt2 <- table( disag_model_matrix[name_period] )
    obs_in_period_noNA <- data.frame( x1 = as.numeric(names(tt2)), nobs = as.numeric(tt2) )
    names(obs_in_period_noNA)[1] <- name_period
    rm(tt2)
    
    # Checking that All Aggregation Periods Are Complete (After Excluding Missing Values)
    obs_in_period_merge <- merge( obs_in_period_original, obs_in_period_noNA, by = name_period )
    if ( any( obs_in_period_merge$nobs.x != obs_in_period_merge$nobs.y ) ) {

      # Excluding Incomplete Aggregation Periods (After Excluding Missing Values)
      to_exclude <- obs_in_period_merge[ which( obs_in_period_merge$nobs.x != obs_in_period_merge$nobs.y ), 1 ]
      disag_model_data[, name_response][ ( disag_model_data[, name_period] %in% to_exclude ) ] <- NA
      if     ( !is.null(CB)            ) {  CB[ ( disag_model_data[, name_period] %in% to_exclude ), 1 ] <- NA 
      } else if ( !is.null(name_exposure) ) { disag_model_data[ ( disag_model_data[, name_period] %in% to_exclude ), name_exposure ] <- NA }
      
      # Re-Fitting the Disaggregated|Disaggregated Model by Randomly Distributing Aggregated Counts across Disaggregated Time Steps 
      model_ini <- try( glm( formula, data = disag_model_data, family = family, na.action = "na.exclude" ), silent = FALSE )

      if ( "try-error" %in% class(model_ini) ) {
        error_code <- -444
        model_ini <- NULL
      } else {
        X <- model.matrix(model_ini)
        disag_model_data_aux <- disag_model_data[ , c(name_period, name_response) ]
        disag_model_matrix <- merge( X, disag_model_data_aux, by = 'row.names' )
        disag_model_matrix <- disag_model_matrix[ order( as.numeric( disag_model_matrix$Row.names ) ), ]
      }
    }
  }
  rm(disag_model_data_aux)
  
  # Outcome Data for the Periods without Missing Values
  if ( error_code == 0 ) {
    response_variable <- Y[ Y[,1] %in% as.numeric( names( table( disag_model_matrix[,name_period] ) ) ), name_response ]
    if ( ncol(X) >= length(response_variable) ) error_code <- -999
  }
  
  # Fitting the Aggregation Model, Trying 'ntryfit' Times Until Convergence
  if ( error_code == 0 ) {
    
    converged <- FALSE
    try_iter <- 0
    while ( !converged & try_iter <= ntryfit ) {
      
      if ( try_iter == 0 ) {
        if ( model_ini$converged ) {
          beta <- coef(model_ini)
        } else {
          beta <- rep( 0, length( coef(model_ini) ) )
          names(beta) <- names( coef(model_ini) )
          beta[ names(beta) == "(Intercept)" ] = log( mean( disag_model_data[,name_response] ) )
        }
      } else {
        if (!is.null(seed)) seed = seed - max(ntryini,ntryfit) * (try_iter-1) + 1
        ini <- get_ini( Y = Y, X = disag_model_data, obs_in_period = obs_in_period_original, name_response = name_response, formula = formula, family = family, seed = seed )
        if ( ini$model_ini$converged ) {
          beta <- coef(ini$model_ini)
        } else {
          if ( !is.null(seed) ) set.seed( seed - max(ntryini,ntryfit) * (try_iter-1) + 1 )
          beta <- rnorm( length( coef(model_ini) ), 0, 2^try_iter / 10 )
          names(beta) <- names( coef(model_ini) )
          beta[ names(beta) == "(Intercept)" ] <- log( mean( disag_model_data[, name_response] ) )
        }
      }
      
      mleBHHH <- try( maxLik::maxLik( logLik = agr_logLik_cpp,
                                      grad = d1score_ind_cpp,
                                      start = beta,
                                      method = "BHHH",
                                      control = list( tol = tol, iterlim = maxit ),
                                      X = X,
                                      y = response_variable,
                                      period = disag_model_matrix[, name_period] ),
                      silent = FALSE )
      if ( "try-error" %in% class(mleBHHH) ) {
        V <- NA
      } else {
        beta <- coef(mleBHHH)
        mle <- try( maxLik::maxLik( logLik = agr_logLik_cpp,
                                    grad = d1func_aggr_cpp,
                                    hess = d2func_aggr_cpp,
                                    start = beta,
                                    method = "NR",
                                    control = list( tol = tol, iterlim = maxit ),
                                    X = X,
                                    y = response_variable,
                                    period = disag_model_matrix[,name_period] ),
                    silent = FALSE )
        if ( "try-error" %in% class(mle) ) {
          V <- NA
        } else {
          beta <- coef(mle)
          V <- vcov(mle)
          if ( all( is.finite(V) ) & ( mle$code %in% c(1,2,8) ) ) {
            if ( all( eigen(V)$values > 0 ) ) converged = TRUE
          } 
        }
      }
      
      try_iter <- try_iter + 1
    }
    rm(converged, try_iter)
  }
  
  # Checking If 'V' Is Positive Definite
  if ( error_code == 0 ) {
    if ( any( is.na(V) ) ) {
      error_code <- -333
    } else if ( any( is.infinite(V) ) ) {
      error_code <- -555
    } else {
      if ( !all( eigen(V)$values > 0 ) ) {
        if ( all( eigen(V)$values >= 0 ) ) {
          error_code <- -777
        } else {
          error_code <- -666
        }
      }
    }
  }
  
  # Returning Values     
  if ( error_code != 0 & error_code != -555 ) {
    res <- return_null_model( error_code, iter = 0, logLik = NULL )
  } else {
    
    # Disaggregated Prediction (Excluding Or Not Missing Values)
    pred_disag <- exp( X %*% beta )
    pred_disag_NA <- exp( X_NA %*% beta )
    
    # Period Prediction
    pred_period <- aggregate( pred_disag, by = list( period = disag_model_matrix[, name_period] ), FUN = sum )$V1
    
    # Pearson Residuals
    pearson_resid <- ( response_variable - pred_period ) / sqrt(pred_period)
    
    # Dispersion Parameter
    if ( family == "quasipoisson" ) { dp <- sum( pearson_resid^2 ) / ( length(response_variable) - length(beta) ) }
    else if ( family == "poisson" ) { dp <- 1 }
    else                            { error_code <- -888 }
    
    # Covariance Matrix with Dispersion
    Vdisp <- V * dp
    
    if (error_code == 0) {
      # Did the Aggregation Model Converge?
      converged <- ( mle$code %in% c(1,2,8) ) & all( eigen(Vdisp)$values > 0 )
      gVdisp <- NULL
    } else if (error_code == -555) { 
      # In this specific case I allow labeling as converged, since a generalized covariance matrix
      #     can be calculated and provides results that are valid
      converged <- ( mle$code %in% c(1,2,8) )  
      
      # if converged, calculate generalized covariance matrix
      if (converged == TRUE) {
        D2 <- d2func_aggr_cpp(beta = beta, X = X, y = response_variable, period = unlist( disag_model_matrix[name_period] ))
        gVdisp <- dp * ginv(-D2)
      } else {
        gVdisp <- NULL
      }
    }
    
    res <- list( coef = beta,
                 vcov = Vdisp,
                 mle = mle,
                 converged = converged,
                 y = response_variable, 
                 X = X,
                 period = unlist( disag_model_matrix[name_period] ),
                 disp = dp,
                 iter = mleBHHH$iterations + mle$iterations,
                 iter_bhhh = mleBHHH$iterations,
                 iter_nr = mle$iterations,
                 logLik = maxLik::maxValue(mle),
                 QAIC = -2 * maxLik::maxValue(mle) + 2 * length(beta) * dp,
                 family = list( family = family, link = "log", linkfun = function(mu) log(mu) ),
                 pred_disag_NA = pred_disag_NA,
                 error_code = error_code, 
                 gvcov = gVdisp)
      class(res) <- "modagr"

  }
  
  return(res)
  
}

# Extraction of Temporally Aggregated Model Coefficients, Covariance Matrix and Generalized Covariance Matrix
coef.modagr <- function( x, ... ){ x$coef }
vcov.modagr <- function( x, ... ){ x$vcov }

gvcov_modagr <- function( x, ...) {
  cnames <- colnames(x$vcov)
  res <- x$gvcov
  colnames(res) <- cnames
  rownames(res) <- cnames
  res
}

# Extraction of model matrix for class modagr
model.matrix.modagr <- function(x,...) { x$X }

# Predict function for class modagr 
predict.modagr <- function(mod,type="response") {
  if (type == "response") mod$pred_disag_NA else if (type == "link") log(mod$pred_disag_NA)
}

# Summary function for class modagr

summary.modagr <- function(object,... ) {
  ## Modification of code for summary.maxLik to integrate the dispersion parameter
  
  ## object      object of class "modagr"
  ## 
  ## RESULTS:
  ## list of class "summary.maxLik" with following components:
  ## maximum    : function value at optimum
  ## estimate   : estimated parameter values at optimum
  ## gradient   :           gradient at optimum
  ## code       : code of convergence
  ## message    : message, description of the code
  ## iterations : number of iterations
  ## type       : type of optimisation
  ##
  if(!inherits(object, "modagr"))
    stop("'summary.modagr' called on a non-'modagr' object")

  result <- object$mle$maxim
  nParam <- length(coef.modagr(object))
  activePar <- activePar( object$mle )
  if((object$mle$code < 100) & !is.null(coef.modagr(object))) {
    # in case of infinity at initial values, the coefs are not provided
    t <- coef( object ) / sqrt(diag(vcov(object)))
    p <- 2*pnorm( -abs( t))
    t[!activePar(object)] <- NA
    p[!activePar(object)] <- NA
    results <- cbind("Estimate" = coef( object ),
                     "Std. error" = sqrt(diag(vcov(object))),
                     "t value" = t, "Pr(> t)" = p )
  } else {
    results <- NULL
  }
  summary <- list(maximType=object$mle$type,
                  iterations=object$iter,
                  returnCode=object$mle$code,
                  returnMessage=object$mle$message,
                  loglik=object$mle$maximum,
                  estimate=results,
                  fixed=!activePar,
                  NActivePar=sum(activePar),
                  constraints=object$mle$constraints)
  class(summary) <- "summary.modagr"
  summary
}

print.summary.modagr <- function( x,
                                  digits = max( 3L, getOption("digits") - 3L ), ... ) {
  
  cat("--------------------------------------------\n")
  cat("Maximum Likelihood estimation\n")
  cat(maximType(x), ", ", nIter(x), " iterations\n", sep="")
  cat("Return code ", returnCode(x), ": ", returnMessage(x), "\n", sep="")
  if(!is.null(x$estimate)) {
    cat("Log-Likelihood:", x$loglik, "\n")
    cat(x$NActivePar, " free parameters\n")
    cat("Estimates:\n")
    printCoefmat( x$estimate, digits = digits )
  }
  if(!is.null(x$constraints)) {
    cat("\nWarning: constrained likelihood estimation.",
        "Inference is probably wrong\n")
    cat("Constrained optimization based on", x$constraints$type,
        "\n")
    if(!is.null(x$constraints$code))
      cat("Return code:", x$constraints$code, "\n")
    # note: this is missing for 'constrOptim'
    if(!is.null(x$constraints$message))
      cat(x$constraints$message, "\n")
    # note: this is missing for 'constrOptim'
    cat(x$constraints$outer.iterations,
        " outer iterations, barrier value",
        x$constraints$barrier.value, "\n")
  }
  cat("--------------------------------------------\n")
}


# Summary that uses generalized covariance instead of vcov

summary_gcov <- function(object,...) {

  ## Modification of code for summary.maxLik to integrate the dispersion parameter
    
  ## object      object of class "modagr"
  ## 
  ## RESULTS:
  ## list of class "summary.maxLik" with following components:
  ## maximum    : function value at optimum
  ## estimate   : estimated parameter values at optimum
  ## gradient   :           gradient at optimum
  ## code       : code of convergence
  ## message    : message, description of the code
  ## iterations : number of iterations
  ## type       : type of optimisation
  ##
  if(!inherits(object, "modagr"))
    stop("'summary.modagr' called on a non-'modagr' object")
    
  result <- object$mle$maxim
  nParam <- length(coef.modagr(object))
  activePar <- activePar( object$mle )
  if((object$mle$code < 100) & !is.null(coef.modagr(object))) {
    # in case of infinity at initial values, the coefs are not provided
    t <- coef( object ) / sqrt(diag( object$gvcov))
    p <- 2*pnorm( -abs( t))
    t[!activePar(object)] <- NA
    p[!activePar(object)] <- NA
    results <- cbind("Estimate" = coef( object ),
                     "Std. error" = sqrt(diag( object$gvcov)),
                     "t value" = t, "Pr(> t)" = p )
  } else {
    results <- NULL
  }
  summary <- list(maximType=object$mle$type,
                  iterations=object$iter,
                  returnCode=object$mle$code,
                  returnMessage=object$mle$message,
                  loglik=object$mle$maximum,
                  estimate=results,
                  fixed=!activePar,
                  NActivePar=sum(activePar),
                  constraints=object$mle$constraints)
  class(summary) <- "summary.modagr"
  summary
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

# Function to calculate attributable fractions and numbers using object of class modagr

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


