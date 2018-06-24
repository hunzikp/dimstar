########################################################################
# MONTE CARLO EXPERIMENTS
# FULL MULTIVARIATE SPATIO-TEMPORAL COUNT & BINARY MODELS
########################################################################

library(dimstar)
library(foreach)
library(doMC)
registerDoMC(10)
mcoptions <- list(preschedule = FALSE)


#################################################
# FUNCTION DEFS
#################################################

###################
# Data creation function
###################

sample_data <- function(config) {
  
  dp <- config$dep_params[[1]]
  dd <- config$data_dim[[1]]
  set.seed(config$seed)
  beta <- c(config$beta0, config$beta1)
  data <- simulate_data(N = dd$N, TT = dd$TT, G = dd$G, count = config$count, 
                        rho_vec = dp$rho_vec, lambda_vec = dp$lambda_vec, gamma_vec = dp$gamma_vec,
                        beta.ls = rep(list(beta), dd$G), sigma2_vec = rep(1, dd$G), 
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
      model$train(maxiter = 50, M = 50, abs_tol = 1e-4, burnin = 0, thinning = 1)
      model$compute_vcov(M = 300, thinning = 1)
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
M <- 50

## Configs
config.df <- expand.grid(data_dim = list(list(N = 36, TT = 10, G = 2), list(N = 256, TT = 10, G = 2)),
                         dep_params = list(list(rho_vec = c(0.25,0.25), gamma_vec = c(0.25,0.25), lambda_vec = 0.25)),
                         beta0 = c(2),
                         beta1 = c(1),
                         count = TRUE,
                         seed = 1:M,
                         model = c("dimstar"),
                         stringsAsFactors = FALSE)
config.ls <- split(config.df, seq(nrow(config.df)))


###################
# Main MC Loop
###################

results.ls <- foreach(config = iter(config.ls), .options.multicore = mcoptions) %dopar% {
  
  ## Generate data
  data <- sample_data(config)
  
  ## Fit approrpiate model
  if (config$model == 'dimstar') {
    results <- fit_dimstar(data)
  }
  
  ## Return result
  results
}

###################
# Concatenate & save results
###################

results.df <- do.call('rbind', lapply(results.ls, function(x) as.data.frame(x)))
results.df <- cbind(config.df, results.df)
saveRDS(results.df, file = 'results/mc_full_rmse_count.rds')


###################
# Plot bias / RMSE
###################

library(ggplot2)
library(data.table)

## Calculate bias
bias.df <- results.df
bias.df$dim <- unlist(lapply(bias.df$data_dim, function(x) paste(paste0(c("N = ", "T = ", "G = "), x), collapse = ", ")))
est_param_names <- names(bias.df)[grep("est", names(bias.df))]
bias.df <- bias.df[,c("dim", est_param_names)]
bias.df$beta0_1 <- bias.df$beta_1_0_est - 2
bias.df$beta0_2 <- bias.df$beta_2_0_est - 2
bias.df$beta1_1 <- bias.df$beta_1_1_est - 1
bias.df$beta1_2 <- bias.df$beta_2_1_est - 1
bias.df$lambda <- bias.df$lambda1_est - 0.25
bias.df$rho_1 <- bias.df$rho_1_est - 0.25
bias.df$rho_2 <- bias.df$rho_2_est - 0.25
bias.df$gamma_1 <- bias.df$gamma_1_est - 0.25
bias.df$gamma_2 <- bias.df$gamma_2_est - 0.25
bias.df <- bias.df[,c("dim", "beta0_1", "beta0_2", "beta1_1", "beta1_2", 
                      "lambda", "rho_1", "rho_2", "gamma_1", "gamma_2")]

## Wide to long
library(reshape2)
bias.df <- melt(bias.df, id.vars = 'dim')

## Prepare faceting
bias.df$dd <- 1
bias.df$dim_f = factor(bias.df$dim, levels=unique(bias.df$dim)[order(unique(bias.df$dim), decreasing = T)])

## Prepare aggregates
agg.dt <- data.table(bias.df)
rmse <- function(x) {sqrt(mean(x^2, na.rm=TRUE))}
agg.df <- data.frame(agg.dt[, j=list(bias = mean(value, na.rm = TRUE), rm = rmse(value)),by = list(dim_f, variable, dd)])

## Determine y limits
ylim <- c(min(bias.df$value) - 0.2, max(bias.df$value) + 0.2)

## Plot errors, bottom row - bias, top row - rmse
p <- ggplot(bias.df, aes(x=dd, y=value)) +  geom_violin()
p <- p + geom_hline(yintercept = 0) + stat_summary(fun.y=mean, geom="point", shape=23, size=2)
p <- p + facet_grid(dim_f ~ variable) + xlab("") + ylab("Error") 
p <- p + theme(axis.title.x=element_blank(),
           axis.text.x=element_blank(),
           axis.ticks.x=element_blank())
p <- p + ylim(ylim[1], ylim[2])
p <- p + geom_text(data=agg.df, aes(x = dd, y = ylim[1] + 0.1, 
                                    label=sprintf("%0.2f", round(bias, digits = 3))), size=4)
p <- p + geom_text(data=agg.df, aes(x = dd, y = ylim[2] - 0.1, 
                                    label=sprintf("%0.2f", round(rm, digits = 3))), size=4)
ggsave(filename = "plots/mc_full_bias_count.pdf", plot = p)


###################
# Coverage prob
###################
library(data.table)

## Coverage data frame
coverage.df <- results.df
coverage.df$dim <- unlist(lapply(coverage.df$data_dim, function(x) paste(paste0(c("N = ", "T = ", "G = "), x), collapse = ", ")))
coverage.df$beta0_1_hit <- 2 < coverage.df$beta_1_0_hi & 2 > coverage.df$beta_1_0_lo
coverage.df$beta0_2_hit <- 2 < coverage.df$beta_2_0_hi & 2 > coverage.df$beta_2_0_lo
coverage.df$beta1_1_hit <- 1 < coverage.df$beta_1_1_hi & 1 > coverage.df$beta_1_1_lo
coverage.df$beta1_2_hit <- 1 < coverage.df$beta_2_1_hi & 1 > coverage.df$beta_2_1_lo
coverage.df$lambda_hit <- 0.25 < coverage.df$lambda1_hi & 0.25 > coverage.df$lambda1_lo
coverage.df$rho_1_hit <- 0.25 < coverage.df$rho_1_hi & 0.25 > coverage.df$rho_1_lo
coverage.df$rho_2_hit <- 0.25 < coverage.df$rho_2_hi & 0.25 > coverage.df$rho_2_lo
coverage.df$gamma_1_hit <- 0.25 < coverage.df$gamma_1_hi & 0.25 > coverage.df$gamma_1_lo
coverage.df$gamma_2_hit <- 0.25 < coverage.df$gamma_2_hi & 0.25 > coverage.df$gamma_2_lo
coverage_var_names <- names(coverage.df)[grep("hit", names(coverage.df))]
coverage.df <- coverage.df[,c("dim", coverage_var_names)]

## Wide to long
coverage.df <- melt(coverage.df, id.vars = 'dim')
cov.dt <- data.table(coverage.df)

## Renampe
cov.dt$variable <- gsub(pattern = "_hit", replacement = "", cov.dt$variable)

## Aggregate
pt_lo <- function(x) {prop.test(x = sum(x), n = length(x))$conf.int[1]}
pt_hi <- function(x) {prop.test(x = sum(x), n = length(x))$conf.int[2]}
cov.dt <- cov.dt[, .(hit_mean = mean(value), hit_hi = pt_hi(value), hit_lo = pt_lo(value)), by = .(dim, variable)]

## Prepare faceting
cov.dt$dd <- 1
cov.dt$dim_f = factor(cov.dt$dim, levels=unique(cov.dt$dim)[order(unique(cov.dt$dim), decreasing = T)])

## Bar plot
p <- ggplot(data=cov.dt, aes(x=dd))
p <- p + geom_bar(aes(weight = hit_mean))
p <- p + geom_errorbar(aes(ymin=hit_lo, ymax=hit_hi), width = 0.2)
p <- p + facet_grid(dim_f ~ variable) + geom_hline(yintercept = 0.9, linetype = 2, size = 0.25)
p <- p + theme(axis.title.x=element_blank(),
              axis.text.x=element_blank(),
              axis.ticks.x=element_blank())
p <- p + ylab("Coverage")
ggsave(filename = "plots/mc_full_coverage_count.pdf", plot = p)



#################################################
# II. EVALUATE BIAS / RMSE FOR BINARY DATA
# Note: Needs M = 300 for SE estimation
#################################################

###################
# Set up configurations
###################

## Constants
M <- 30

## Configs
config.df <- expand.grid(data_dim = list(list(N = 36, TT = 10, G = 2), list(N = 256, TT = 10, G = 2)),
                         dep_params = list(list(rho_vec = c(0.25,0.25), gamma_vec = c(0.25,0.25), lambda_vec = 0.25)),
                         beta0 = c(-0.5),
                         beta1 = c(1),
                         count = FALSE,
                         seed = 1:M,
                         model = c("dimstar"),
                         stringsAsFactors = FALSE)
config.ls <- split(config.df, seq(nrow(config.df)))


###################
# Main MC Loop
###################

results.ls <- foreach(config = iter(config.ls), .options.multicore = mcoptions) %dopar% {
  
  ## Generate data
  data <- sample_data(config)
  
  ## Fit approrpiate model
  if (config$model == 'dimstar') {
    results <- fit_dimstar(data)
  }
  
  ## Return result
  results
}

###################
# Concatenate & save results
###################

results.df <- do.call('rbind', lapply(results.ls, function(x) as.data.frame(x)))
results.df <- cbind(config.df, results.df)
saveRDS(results.df, file = 'results/mc_full_rmse_binary.rds')


###################
# Plot bias / RMSE
###################

library(ggplot2)
library(data.table)

## Calculate bias
bias.df <- results.df
bias.df$dim <- unlist(lapply(bias.df$data_dim, function(x) paste(paste0(c("N = ", "T = ", "G = "), x), collapse = ", ")))
est_param_names <- names(bias.df)[grep("est", names(bias.df))]
bias.df <- bias.df[,c("dim", est_param_names)]
bias.df$beta0_1 <- bias.df$beta_1_0_est - -0.5
bias.df$beta0_2 <- bias.df$beta_2_0_est - -0.5
bias.df$beta1_1 <- bias.df$beta_1_1_est - 1
bias.df$beta1_2 <- bias.df$beta_2_1_est - 1
bias.df$lambda <- bias.df$lambda1_est - 0.25
bias.df$rho_1 <- bias.df$rho_1_est - 0.25
bias.df$rho_2 <- bias.df$rho_2_est - 0.25
bias.df$gamma_1 <- bias.df$gamma_1_est - 0.25
bias.df$gamma_2 <- bias.df$gamma_2_est - 0.25
bias.df <- bias.df[,c("dim", "beta0_1", "beta0_2", "beta1_1", "beta1_2", 
                      "lambda", "rho_1", "rho_2", "gamma_1", "gamma_2")]

## Wide to long
library(reshape2)
bias.df <- melt(bias.df, id.vars = 'dim')

## Prepare faceting
bias.df$dd <- 1
bias.df$dim_f = factor(bias.df$dim, levels=unique(bias.df$dim)[order(unique(bias.df$dim), decreasing = T)])

## Prepare aggregates
agg.dt <- data.table(bias.df)
rmse <- function(x) {sqrt(mean(x^2, na.rm=TRUE))}
agg.df <- data.frame(agg.dt[, j=list(bias = mean(value, na.rm = TRUE), rm = rmse(value)),by = list(dim_f, variable, dd)])

## Determine y limits
ylim <- c(min(bias.df$value) - 0.2, max(bias.df$value) + 0.2)

## Plot errors, bottom row - bias, top row - rmse
p <- ggplot(bias.df, aes(x=dd, y=value)) +  geom_violin()
p <- p + geom_hline(yintercept = 0) + stat_summary(fun.y=mean, geom="point", shape=23, size=2)
p <- p + facet_grid(dim_f ~ variable) + xlab("") + ylab("Error") 
p <- p + theme(axis.title.x=element_blank(),
               axis.text.x=element_blank(),
               axis.ticks.x=element_blank())
p <- p + ylim(ylim[1], ylim[2])
p <- p + geom_text(data=agg.df, aes(x = dd, y = ylim[1] + 0.1, 
                                    label=sprintf("%0.2f", round(bias, digits = 3))), size=4)
p <- p + geom_text(data=agg.df, aes(x = dd, y = ylim[2] - 0.1, 
                                    label=sprintf("%0.2f", round(rm, digits = 3))), size=4)
ggsave(filename = "plots/mc_full_bias_binary.pdf", plot = p)


###################
# Coverage prob
###################
library(data.table)

## Coverage data frame
coverage.df <- results.df
coverage.df$dim <- unlist(lapply(coverage.df$data_dim, function(x) paste(paste0(c("N = ", "T = ", "G = "), x), collapse = ", ")))
coverage.df$beta0_1_hit <- -0.5 < coverage.df$beta_1_0_hi & -0.5 > coverage.df$beta_1_0_lo
coverage.df$beta0_2_hit <- -0.5 < coverage.df$beta_2_0_hi & -0.5 > coverage.df$beta_2_0_lo
coverage.df$beta1_1_hit <- 1 < coverage.df$beta_1_1_hi & 1 > coverage.df$beta_1_1_lo
coverage.df$beta1_2_hit <- 1 < coverage.df$beta_2_1_hi & 1 > coverage.df$beta_2_1_lo
coverage.df$lambda_hit <- 0.25 < coverage.df$lambda1_hi & 0.25 > coverage.df$lambda1_lo
coverage.df$rho_1_hit <- 0.25 < coverage.df$rho_1_hi & 0.25 > coverage.df$rho_1_lo
coverage.df$rho_2_hit <- 0.25 < coverage.df$rho_2_hi & 0.25 > coverage.df$rho_2_lo
coverage.df$gamma_1_hit <- 0.25 < coverage.df$gamma_1_hi & 0.25 > coverage.df$gamma_1_lo
coverage.df$gamma_2_hit <- 0.25 < coverage.df$gamma_2_hi & 0.25 > coverage.df$gamma_2_lo
coverage_var_names <- names(coverage.df)[grep("hit", names(coverage.df))]
coverage.df <- coverage.df[,c("dim", coverage_var_names)]

## Wide to long
coverage.df <- melt(coverage.df, id.vars = 'dim')
cov.dt <- data.table(coverage.df)

## Renampe
cov.dt$variable <- gsub(pattern = "_hit", replacement = "", cov.dt$variable)

## Aggregate
pt_lo <- function(x) {prop.test(x = sum(x), n = length(x))$conf.int[1]}
pt_hi <- function(x) {prop.test(x = sum(x), n = length(x))$conf.int[2]}
cov.dt <- cov.dt[, .(hit_mean = mean(value), hit_hi = pt_hi(value), hit_lo = pt_lo(value)), by = .(dim, variable)]

## Prepare faceting
cov.dt$dd <- 1
cov.dt$dim_f = factor(cov.dt$dim, levels=unique(cov.dt$dim)[order(unique(cov.dt$dim), decreasing = T)])

## Bar plot
p <- ggplot(data=cov.dt, aes(x=dd))
p <- p + geom_bar(aes(weight = hit_mean))
p <- p + geom_errorbar(aes(ymin=hit_lo, ymax=hit_hi), width = 0.2)
p <- p + facet_grid(dim_f ~ variable) + geom_hline(yintercept = 0.9, linetype = 2, size = 0.25)
p <- p + theme(axis.title.x=element_blank(),
               axis.text.x=element_blank(),
               axis.ticks.x=element_blank())
p <- p + ylab("Coverage")
ggsave(filename = "plots/mc_full_coverage_binary.pdf", plot = p)



#################################################
# III. EVALUATE COMPLEXITY (FOR COUNT DATA)
#################################################


###################
# Data creation function
###################

sample_data <- function(config) {

  dd <- config$data_dim[[1]]
  
  if (dd$G == 1) {
    lambda_vec <- 0
  } else {
    n_pairs <- choose(dd$G, 2)
    lambda_vec <- rep(0.25, n_pairs)
  }
  
  if (dd$TT == 1) {
    gamma_vec <- rep(0, dd$G)
  } else {
    gamma_vec <- rep(0.25, dd$G)
  }
  
  rho_vec <- rep(0.25, dd$G)
  
  set.seed(config$seed)
  beta <- c(config$beta0, config$beta1)
  data <- simulate_data(N = dd$N, TT = dd$TT, G = dd$G, count = config$count, 
                        rho_vec = rho_vec, lambda_vec = lambda_vec, gamma_vec = gamma_vec,
                        beta.ls = rep(list(beta), dd$G), sigma2_vec = rep(1, dd$G), 
                        X_params = c(0, 1))
  return(data)
}

###################
# Data estimation function
###################

fit_dimstar <- function(data) {
  
  if (data$G == 1) {
    outcome_dep <- FALSE
  } else {
    outcome_dep <- TRUE
  }
  if (data$TT == 1) {
    temporal_dep <- FALSE
  } else {
    temporal_dep <- TRUE
  }
  
  model <- DISTAR$new(X.ls = data$X.ls, y.ls = data$y.ls, W_t = data$W_t, 
                      N = data$N, G = data$G, TT = data$TT,
                      count = data$count, 
                      spatial_dep = TRUE, 
                      temporal_dep = temporal_dep, 
                      outcome_dep = outcome_dep)
  model$sample_ystar(M = 50, ystar_init = colMeans(model$ystar_sample))
  
  timing <- system.time({
    res <- try({
      set.seed(0)
      # model$train(maxiter = 1, M = 50, abs_tol = 1e-4, burnin = 0, thinning = 1)
      ll <- model$E_llik(theta = model$pack_theta())
    }, silent = TRUE)
  })
  
  results <- list()
  results$elapsed <- timing[3]
  
  return(results)
}

###################
# Set up configurations
###################

## Constants
M <- 5

## Configs
config.df <- expand.grid(data_dim = list(list(N = 64, TT = 1, G = 1),
                                         list(N = 144, TT = 1, G = 1),
                                         list(N = 256, TT = 1, G = 1),
                                         list(N = 400, TT = 1, G = 1), 
                                         list(N = 576, TT = 1, G = 1), 
                                         list(N = 784, TT = 1, G = 1), 
                                         list(N = 1024, TT = 1, G = 1), 
                                         
                                         list(N = 64, TT = 1, G = 4),
                                         list(N = 144, TT = 1, G = 4),
                                         list(N = 256, TT = 1, G = 4),
                                         list(N = 400, TT = 1, G = 4), 
                                         list(N = 576, TT = 1, G = 4), 
                                         list(N = 784, TT = 1, G = 4), 
                                         list(N = 1024, TT = 1, G = 4), 
                                         
                                         list(N = 64, TT = 4, G = 1), 
                                         list(N = 64, TT = 8, G = 1),
                                         list(N = 64, TT = 12, G = 1),
                                         list(N = 64, TT = 16, G = 1),
                                         list(N = 64, TT = 20, G = 1),
                                         list(N = 64, TT = 24, G = 1),
                                         list(N = 64, TT = 48, G = 1),
                                         
                                         list(N = 64, TT = 4, G = 4), 
                                         list(N = 64, TT = 8, G = 4),
                                         list(N = 64, TT = 12, G = 4),
                                         list(N = 64, TT = 16, G = 4),
                                         list(N = 64, TT = 20, G = 4),
                                         list(N = 64, TT = 24, G = 4),
                                         list(N = 64, TT = 48, G = 4),
                                         
                                         list(N = 64, TT = 1, G = 4),
                                         list(N = 64, TT = 1, G = 8),
                                         list(N = 64, TT = 1, G = 12),
                                         list(N = 64, TT = 1, G = 16),
                                         list(N = 64, TT = 1, G = 20)),
                         beta0 = c(2),
                         beta1 = c(1),
                         count = TRUE,
                         seed = 1:M,
                         model = c("dimstar"),
                         stringsAsFactors = FALSE)
config.df$iter <- 1:nrow(config.df)
config.ls <- split(config.df, seq(nrow(config.df)))

###################
# Main MC Loop
###################

results.ls <- foreach(config = iter(config.ls), .options.multicore = mcoptions) %dopar% {
  
  # print(paste(config$iter, '/', nrow(config.df)))
  
  ## Generate data
  data <- sample_data(config)
  
  ## Fit approrpiate model
  if (config$model == 'dimstar') {
    results <- fit_dimstar(data)
  }
  
  ## Return result
  results
}

###################
# Concatenate & save results
###################

results.df <- do.call('rbind', lapply(results.ls, function(x) as.data.frame(x[length(x)])))
results.df <- cbind(config.df, results.df)
saveRDS(results.df, file = 'results/mc_full_complexity_count.rds')

###################
# Plot Scalability in N
###################
library(data.table)
library(ggplot2)

scal.df <- results.df
scal.df$N <- unlist(lapply(scal.df$data_dim, function(x) x$N))
scal.df$TT <- unlist(lapply(scal.df$data_dim, function(x) x$TT))
scal.df$G <- unlist(lapply(scal.df$data_dim, function(x) x$G))
scal.dt <- data.table(scal.df[scal.df$TT == 1 & scal.df$G <= 3,])
scal.dt <- scal.dt[, j=list(elapsed = mean(elapsed, na.rm = TRUE)), by = list(N, G)]

scal.dt$G <- as.factor(scal.dt$G)

p <- ggplot(data=scal.dt, aes(x=N, y=elapsed, group = G)) +
  geom_line(aes(linetype=G)) +
  geom_point(aes(shape=G)) + xlab("N") +
  theme(legend.position="bottom") + 
  xlab("N") + ylab("Evaluation Time (seconds)") +
  theme(legend.position="bottom")
p
ggsave(filename = "plots/mc_full_time_N.pdf", plot = p, scale = 0.5)


###################
# Plot Scalability in TT
###################

scal.df <- results.df
scal.df$N <- unlist(lapply(scal.df$data_dim, function(x) x$N))
scal.df$TT <- unlist(lapply(scal.df$data_dim, function(x) x$TT))
scal.df$G <- unlist(lapply(scal.df$data_dim, function(x) x$G))
scal.dt <- data.table(scal.df[scal.df$N == 64 & scal.df$G <= 2,])
scal.dt <- scal.dt[, j=list(elapsed = mean(elapsed, na.rm = TRUE)), by = list(TT, G)]

scal.dt$G <- as.factor(scal.dt$G)

p <- ggplot(data=scal.dt, aes(x=TT, y=elapsed, group = G)) +
  geom_line(aes(linetype=G)) +
  geom_point(aes(shape=G)) + xlab("N") +
  theme(legend.position="bottom") + 
  xlab("T") + ylab("Evaluation Time (seconds)") +
  theme(legend.position="bottom")
p
ggsave(filename = "plots/mc_full_time_T.pdf", plot = p, scale = 0.5)

###################
# Plot Scalability in G
###################

scal.df <- results.df
scal.df$N <- unlist(lapply(scal.df$data_dim, function(x) x$N))
scal.df$TT <- unlist(lapply(scal.df$data_dim, function(x) x$TT))
scal.df$G <- unlist(lapply(scal.df$data_dim, function(x) x$G))
scal.dt <- data.table(scal.df[scal.df$N == 64 & scal.df$TT == 1,])
scal.dt <- scal.dt[, j=list(elapsed = mean(elapsed, na.rm = TRUE)), by = list(G)]

p <- ggplot(data=scal.dt, aes(x=G, y=elapsed)) +
  geom_line() + 
  geom_point() + 
  xlab("G") + ylab("Evaluation Time (seconds)")
p
ggsave(filename = "plots/mc_full_time_G.pdf", plot = p, scale = 0.5)

