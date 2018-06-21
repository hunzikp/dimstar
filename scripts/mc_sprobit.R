########################################################################
# MONTE CARLO EXPERIMENTS
# CROSS-SETIONAL UNIVERIATE SPATIAL PROBIT
########################################################################

library(dimstar)
library(spatialprobit)
library(McSpatial)
library(foreach)
library(doMC)
registerDoMC(6)

###################
# Set up configurations
###################

## Constants
M <- 30
BAYES_MAX <- 256
GMM_MAX <- 256
config.df <- expand.grid(N = c(64, 256),
                         rho = c(0, 0.5),
                         beta0 = c(-0.5),
                         beta1 = c(1),
                         seed = 1:M,
                         model = c("dimstar", "bayes", "gmm"),
                         stringsAsFactors = FALSE)
config.ls <- split(config.df, seq(nrow(config.df)))


###################
# Data creation function
###################

sample_data <- function(config) {
  
  set.seed(config$seed)
  beta <- c(config$beta0, config$beta1)
  data <- simulate_data(N = config$N, TT = 1, G = 1, count = FALSE, 
                        rho_vec = config$rho, lambda_vec = 0, gamma_vec = 0,
                        beta.ls = list(beta), sigma2_vec = 1)
  return(data)
}

###################
# Data estimation functions
###################

fit_dimstar <- function(data) {
  
  model <- DISTAR$new(X.ls = data$X.ls, y.ls = data$y.ls, W_t = data$W_t, 
                      N = data$N, G = 1, TT = 1,
                      count = FALSE, 
                      spatial_dep = TRUE, 
                      temporal_dep = FALSE, 
                      outcome_dep = FALSE)

  timing <- system.time({
    res <- try({
      set.seed(0)
      model$train(maxiter = 50, M = 50, abs_tol = 1e-5, burnin = 0, thinning = 1)
      model$compute_vcov(M = 100, thinning = 1)
    }, silent = TRUE)
  })
  
  if(inherits(res, "try-error")) {
    
    # 'Gracefully' handle crash
    results <- as.list(rep(NA, 9))
    timing <- rep(NA, 3)
    
  } else {
    
    estim <- model$coeftest(print_out = F)
    theta_est <- estim$coef
    se <- estim$se
    
    theta_lo <- theta_est - qnorm(0.95)*se
    theta_hi <- theta_est +  qnorm(0.95)*se
    
    results <- as.list(c(theta_est, theta_lo, theta_hi))
  }
  
  nms <- c('beta0', 'beta1', 'rho')
  nms <- c(paste0(nms, '_est'), paste0(nms, '_lo'), paste0(nms, '_hi'))
  names(results) <- nms
  
  results$elapsed <- timing[3]
  
  return(results)
}

fit_bayes <- function(data) {
  
  X <- data$X.ls[[1]]
  y <- data$y.ls[[1]]
  W <- data$W_t
  N <- data$N
  
  if (N > BAYES_MAX) {
    
    # Terminate gracefully
    results <- as.list(rep(NA, 9))
    timing <- rep(NA, 3)
    
  } else {
    
    df <- as.data.frame(X[,-1,drop=FALSE])
    df$y <- y
    timing <- system.time({
      g.fit <- sarprobit(formula = y ~ ., W = W, data = df)
    })
    
    theta_est <- g.fit$coefficients
    
    se <- c(sqrt(diag(var(g.fit$bdraw))), sd(g.fit$pdraw))
    theta_lo <- theta_est - qnorm(0.95)*se
    theta_hi <- theta_est +  qnorm(0.95)*se
    
    results <- as.list(c(theta_est, theta_lo, theta_hi))
  }
  
  nms <- c('beta0', 'beta1', 'rho')
  nms <- c(paste0(nms, '_est'), paste0(nms, '_lo'), paste0(nms, '_hi'))
  names(results) <- nms
  
  results$elapsed <- timing[3]
  
  return(results)
}

fit_gmm <- function(data) {
  
  X <- data$X.ls[[1]]
  y <- data$y.ls[[1]]
  W <- data$W_t
  N <- data$N
  
  if (N > GMM_MAX) {
    
    # Terminate
    results <- as.list(rep(NA, 9))
    timing <- rep(NA, 3)
    
  } else {
    
    df <- as.data.frame(X[,-1,drop=FALSE])
    df$y <- y
    timing <- system.time({
      res <- try({
        g.fit <- gmmprobit(form = y ~ ., shpfile = NULL, wmat = W, data = df, silent = TRUE,
                           cvcrit = sqrt(.Machine$double.eps))  ## level playing field...
      }, silent = TRUE)
    })
    
    if(inherits(res, "try-error")) {
      # 'Gracefully' handle crash
      
      results <- as.list(rep(NA, 9))
      timing <- rep(NA, 3)
      
    } else {
      
      theta_est <- g.fit$coef
      
      theta_lo <- theta_est - qnorm(0.95)*g.fit$se
      theta_hi <- theta_est +  qnorm(0.95)*g.fit$se
      
      results <- as.list(c(theta_est, theta_lo, theta_hi))
    }
  }
  
  nms <- c('beta0', 'beta1', 'rho')
  nms <- c(paste0(nms, '_est'), paste0(nms, '_lo'), paste0(nms, '_hi'))
  names(results) <- nms
  
  results$elapsed <- timing[3]
  
  return(results)
}


###################
# Main MC Loop
###################

results.ls <- foreach(config = iter(config.ls)) %dopar% {
  
  ## Generate data
  data <- sample_data(config)
  
  ## Fit approrpiate model
  if (config$model == 'dimstar') {
    results <- fit_dimstar(data)
  } else if (config$model == 'bayes') {
    results <- fit_bayes(data)
  } else if (config$model == 'gmm') {
    results <- fit_gmm(data)
  }
  
  ## Return result
  results
}

###################
# Concatenate results
###################

results.df <- do.call('rbind', lapply(results.ls, function(x) as.data.frame(x)))
results.df <- cbind(config.df, results.df)


###################
# Plot bias / RMSE
###################

library(ggplot2)

## Rho bias
rho_bias.df <- results.df
rho_bias.df$rb <- rho_bias.df$rho_est - rho_bias.df$rho

rho_bias.df$rho_label <- paste0('rho = ', rho_bias.df$rho)
rho_bias.df$N_label <- paste0('N = ', rho_bias.df$N)

p <- ggplot(rho_bias.df, aes(x=model, y=rb)) +  geom_violin()
p <- p + geom_hline(yintercept = 0) + stat_summary(fun.y=mean, geom="point", shape=23, size=2)
p <- p + facet_grid(rho_label ~ N_label) + xlab("") + ylab("rho bias")
p
