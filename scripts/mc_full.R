########################################################################
# MONTE CARLO EXPERIMENTS
# FULL MULTIVARIATE SPATIO-TEMPORAL COUNT & BINARY MODELS
########################################################################

library(dimstar)
library(foreach)
library(doMC)
registerDoMC(10)

#################################################
# FUNCTION DEFS
#################################################

###################
# Data creation function
###################

sample_data <- function(config) {
  
  dp <- config$dep_params
  set.seed(config$seed)
  beta <- c(config$beta0, config$beta1)
  data <- simulate_data(N = config$N, TT = config$TT, G = config$G, count = config$count, 
                        rho_vec = dp$rho_vec, lambda_vec = dp$lambda_vec, gamma_vec = dp$gamma_vec,
                        beta.ls = rep(list(beta), config$G), sigma2_vec = rep(1, config$G), 
                        X_params = c(0, 1))
  return(data)
}

###################
# Data estimation functions
###################

fit_dimstar <- function(data) {
  
  model <- DISTAR$new(X.ls = data$X.ls, y.ls = data$y.ls, W_t = data$W_t, 
                      N = data$N, G = data$G, TT = data$TT,
                      count = data$count, 
                      spatial_dep = TRUE, 
                      temporal_dep = TRUE, 
                      outcome_dep = TRUE)
  
  timing <- system.time({
    res <- try({
      set.seed(0)
      model$train(maxiter = 25, M = 50, abs_tol = 1e-3, burnin = 0, thinning = 1)
      model$compute_vcov(M = 100, thinning = 1)
    }, silent = TRUE)
  })
  
  if(inherits(res, "try-error")) {
    
    # 'Gracefully' handle crash
    results <- as.list(rep(NA, length(model$pack_theta())))
    timing <- rep(NA, 3)
    
  } else {
    
    estim <- model$coeftest(print_out = F)
    theta_est <- estim$coef
    se <- estim$se
    
    theta_lo <- theta_est - qnorm(0.95)*se
    theta_hi <- theta_est +  qnorm(0.95)*se
    
    results <- as.list(c(theta_est, theta_lo, theta_hi))
  }
  
  nms <- names(theta_est)
  nms <- c(paste0(nms, '_est'), paste0(nms, '_lo'), paste0(nms, '_hi'))
  names(results) <- nms
  
  results$elapsed <- timing[3]
  
  return(results)
}


#################################################
# I. EVALUATE BIAS / RMSE FOR COUNT DATA
#################################################

###################
# Set up configurations
###################

## Constants
M <- 1

## Configs
config.df <- expand.grid(N = 64,
                         TT = 2,
                         G = 2,
                         dep_params = list(list(rho_vec = c(0,0), gamma_vec = c(0,0), lambda_vec = 0)),
                         beta0 = c(2),
                         beta1 = c(1),
                         count = TRUE,
                         seed = 1:M,
                         model = c("dimstar"),
                         stringsAsFactors = FALSE)
config.ls <- split(config.df, seq(nrow(config.df)))