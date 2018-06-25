########################################################################
# MONTE CARLO EXPERIMENTS
# CROSS-SETIONAL UNIVERIATE SPATIAL PROBIT
########################################################################

library(dimstar)
library(spatialprobit)
library(McSpatial)
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
  
  set.seed(config$seed)
  beta <- c(config$beta0, config$beta1)
  data <- simulate_data(N = config$N, TT = 1, G = 1, count = FALSE, 
                        rho_vec = config$rho, lambda_vec = 0, gamma_vec = 0,
                        beta.ls = list(beta), sigma2_vec = 1, 
                        X_params = c(0, 1))
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
      model$train(maxiter = 75, M = 50, abs_tol = 1e-3, burnin = 0, thinning = 1)
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
        g.fit <- gmmprobit(form = y ~ ., shpfile = NULL, wmat = as.matrix(W), data = df, silent = TRUE,
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

fit_em <- function(data) {
  
  X <- data$X.ls[[1]]
  y <- data$y.ls[[1]]
  W <- data$W_t
  N <- data$N
  
  if (N > EM_MAX) {
    
    # Terminate
    results <- as.list(rep(NA, 9))
    timing <- rep(NA, 3)
    
  } else {
    
    df <- as.data.frame(X[,-1,drop=FALSE])
    df$y <- y
    timing <- system.time({
      res <- try({
        g.fit <- spprobitml(form = y ~ ., shpfile = NULL, wmat = as.matrix(W), data = df, stdprobit = FALSE)
      }, silent = TRUE)
    })
    
    if(inherits(res, "try-error")) {
      # 'Gracefully' handle crash
      
      results <- as.list(rep(NA, 9))
      timing <- rep(NA, 3)
      
    } else {
      
      theta_est <- g.fit$coef
      
      se <- sqrt(diag(g.fit$vmat2))
      theta_lo <- theta_est - qnorm(0.95)*se
      theta_hi <- theta_est +  qnorm(0.95)*se
      
      results <- as.list(c(theta_est, theta_lo, theta_hi))
    }
  }
  
  nms <- c('beta0', 'beta1', 'rho')
  nms <- c(paste0(nms, '_est'), paste0(nms, '_lo'), paste0(nms, '_hi'))
  names(results) <- nms
  
  results$elapsed <- timing[3]
  
  return(results)
}

#################################################
# I. EVALUATE BIAS / RMSE
#################################################

###################
# Set up configurations
###################

## Constants
M <- 30
BAYES_MAX <- 1024
GMM_MAX <- 1024
EM_MAX <- 1024

## Configs
config.df <- expand.grid(N = c(64, 1024),
                         rho = c(0, 0.5, 0.8),
                         beta0 = c(0),
                         beta1 = c(2),
                         seed = 1:M,
                         model = c("dimstar", "bayes", "gmm"),
                         stringsAsFactors = FALSE)
config.ls <- split(config.df, seq(nrow(config.df)))

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
  } else if (config$model == 'em') {
    results <- fit_em(data)
  }
  
  ## Return result
  results
}

###################
# Concatenate & save results
###################

results.df <- do.call('rbind', lapply(results.ls, function(x) as.data.frame(x)))
results.df <- cbind(config.df, results.df)
saveRDS(results.df, file = 'results/mc_sprobit_rmse.rds')

# results.df <- readRDS('results/mc_sprobit_rmse.rds')


###################
# Plot bias / RMSE
###################

library(ggplot2)
library(data.table)
library(gridExtra)

## Calculate Rho bias
rho_bias.df <- results.df[results.df$model != 'em',]
rho_bias.df$model[rho_bias.df$model == 'dimstar'] <- 'mcem'
rho_bias.df$value <- rho_bias.df$rho_est - rho_bias.df$rho

rho_bias.df$rho_label <- paste0('rho = ', rho_bias.df$rho)
rho_bias.df$N_label <- paste0('N = ', rho_bias.df$N)
rho_bias.df$N_label = factor(rho_bias.df$N_label, levels=unique(rho_bias.df$N_label)[order(unique(rho_bias.df$N_label), decreasing = T)])

## Prepare aggregates
agg.dt <- data.table(rho_bias.df)
rmse <- function(x) {sqrt(mean(x^2, na.rm=TRUE))}
agg.df <- data.frame(agg.dt[, j=list(bias = mean(value, na.rm = TRUE), rm = rmse(value)),by = list(rho_label, N_label, model)])

## Determine y limits
ylim <- c(min(rho_bias.df$value, na.rm = T) - 0.3, max(rho_bias.df$value, na.rm = T) + 0.3)

p <- ggplot(rho_bias.df, aes(x=model, y=value)) +  geom_violin()
p <- p + geom_hline(yintercept = 0) + stat_summary(fun.y=mean, geom="point", shape=23, size=2)
p <- p + facet_grid(rho_label ~ N_label) + xlab("") + ylab("rho bias")
p <- p + ylim(ylim[1], ylim[2])
p <- p + geom_text(data=agg.df, aes(x = model, y = ylim[1] + 0.1, 
                                    label=sprintf("%0.2f", round(bias, digits = 3))), size=4)
p1 <- p + geom_text(data=agg.df, aes(x = model, y = ylim[2] - 0.1, 
                                    label=sprintf("%0.2f", round(rm, digits = 3))), size=4)

## Calculate Beta bias
beta_bias.df <- results.df[results.df$model != 'em',]
beta_bias.df$model[beta_bias.df$model == 'dimstar'] <- 'mcem'
beta_bias.df$value <- beta_bias.df$beta1_est - beta_bias.df$beta1

beta_bias.df$rho_label <- paste0('rho = ', beta_bias.df$rho)
beta_bias.df$N_label <- paste0('N = ', beta_bias.df$N)
beta_bias.df$N_label = factor(beta_bias.df$N_label, levels=unique(beta_bias.df$N_label)[order(unique(beta_bias.df$N_label), decreasing = T)])

## Prepare aggregates
agg.dt <- data.table(beta_bias.df)
rmse <- function(x) {sqrt(mean(x^2, na.rm=TRUE))}
agg2.df <- data.frame(agg.dt[, j=list(bias = mean(value, na.rm = TRUE), rm = rmse(value)),by = list(rho_label, N_label, model)])

## Determine y limits
ylim2 <- c(min(beta_bias.df$value, na.rm = T) - 0.3, max(beta_bias.df$value, na.rm = T) + 0.3)

p <- ggplot(beta_bias.df, aes(x=model, y=value)) +  geom_violin()
p <- p + geom_hline(yintercept = 0) + stat_summary(fun.y=mean, geom="point", shape=23, size=2)
p <- p + facet_grid(rho_label ~ N_label) + xlab("") + ylab("beta bias")
p <- p + ylim(ylim2[1], ylim2[2])
p <- p + geom_text(data=agg2.df, aes(x = model, y = ylim2[1] + 0.1, 
                                    label=sprintf("%0.2f", round(bias, digits = 3))), size=4)
p2 <- p + geom_text(data=agg2.df, aes(x = model, y = ylim2[2] - 0.1, 
                                     label=sprintf("%0.2f", round(rm, digits = 3))), size=4)

p <- arrangeGrob(p1, p2, nrow=1, ncol=2) 
ggsave(filename = "plots/mc_sprobit_bias.pdf", plot = p)



#################################################
# II. EVALUATE RUNTIME / COMPLEXITY
#################################################

###################
# Set up configurations
###################

## Constants
M <- 2
BAYES_MAX <- 4096
GMM_MAX <- 1600
EM_MAX <- 1600

## Configs
config.df <- expand.grid(N = seq(8, 64, by = 8)^2,
                         rho = c(0.5),
                         beta0 = c(0),
                         beta1 = c(2),
                         seed = 1:M,
                         model = c("dimstar", "bayes", "gmm"),
                         stringsAsFactors = FALSE)
config.ls <- split(config.df, seq(nrow(config.df)))

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
  } else if (config$model == 'em') {
    results <- fit_em(data)
  }
  
  ## Return result
  results
}

###################
# Concatenate & save results
###################

results.df <- do.call('rbind', lapply(results.ls, function(x) as.data.frame(x)))
results.df <- cbind(config.df, results.df)
saveRDS(results.df, file = 'results/mc_sprobit_complexity.rds')

# results.df <- readRDS('results/mc_sprobit_complexity.rds')

###################
# Plot Scalability
###################

scal.dt <- data.table(results.df[results.df$model != 'em',])
scal.dt$model[scal.dt$model == 'dimstar'] <- 'mcem'
scal.dt <- scal.dt[, j=list(elapsed = mean(elapsed, na.rm = TRUE)), by = list(N, model)]

p <- ggplot(data=scal.dt, aes(x=N, y=elapsed, group=model)) +
  geom_line(aes(linetype=model)) +
  geom_point(aes(shape=model)) + xlab("N") + ylab("Estimation Time") +
  theme(legend.position="bottom")
ggsave(filename = "plots/mc_sprobit_complexity.pdf", plot = p, scale = 0.5)



