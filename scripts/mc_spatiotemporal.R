########################################################################
# MONTE CARLO EXPERIMENTS
# FULL UNIVARIATE SPATIO-TEMPORAL COUNT & BINARY MODELS
########################################################################

library(dimstar)
library(reticulate)
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
  data <- simulate_data(N = dd$N, TT = dd$TT, G = 1L, count = config$count, 
                        rho_vec = dp$rho_vec, lambda_vec = 0, gamma_vec = dp$gamma_vec,
                        beta.ls = list(beta), sigma2_vec = 1, 
                        X_params = c(0, 1))
  return(data)
}

###################
# Data estimation functions
###################

fit_dimstar <- function(data) {
  
  model <- DISTAR$new(X.ls = data$X.ls, y.ls = data$y.ls, W_t = data$W_t, 
                      N = data$N, G = 1, TT = data$TT,
                      count = data$count, 
                      spatial_dep = TRUE, 
                      temporal_dep = TRUE, 
                      outcome_dep = FALSE)
  
  timing <- system.time({
    res <- try({
      set.seed(0)
      model$train(maxiter = 200, M = 50, abs_tol = 1e-5, burnin = 0, thinning = 1, verbose = TRUE, soft_init = FALSE)
      model$compute_vcov(M = 500, thinning = 5)
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
    
    nms <- names(theta_est)
    nms <- c(paste0(nms, '_est'), paste0(nms, '_lo'), paste0(nms, '_hi'))
    names(results) <- nms
  }
  
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
config.df <- expand.grid(data_dim = list(list(N = 64, TT = 10, G = 1), list(N = 256, TT = 10, G = 1)),
                         dep_params = list(list(rho_vec = 0, gamma_vec = 0),
                                           list(rho_vec = 0.25, gamma_vec = 0.25),
                                           list(rho_vec = 0.45, gamma_vec = 0.45)),
                         beta0 = c(1),
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
saveRDS(results.df, file = 'results/mc_spatiotemporal_rmse_count.rds')


###################
# Print average run-times
###################

runtime.df <- results.df
runtime.df$dim <- unlist(lapply(runtime.df$data_dim, function(x) paste(paste0(c("N = ", "T = ", "G = "), x), collapse = ", ")))
runtime.dt <- data.table(runtime.df)
runtime.dt <- runtime.dt[,list(runtime = mean(elapsed)), by = list(dim)]
print(runtime.dt$runtime)
print(runtime.dt$runtime/60)

###################
# Plot bias / RMSE
###################

library(ggplot2)
library(data.table)

## Calculate bias
dimvec <- unlist(lapply(results.df$data_dim, function(x) paste(paste0(c("N = ", "T = "), x[-3]), collapse = ", ")))
depvec <- unlist(lapply(results.df$dep_params, function(x) paste(paste0(c("rho = ", "gamma = "), x), collapse = ", ")))
depnum <- unlist(lapply(results.df$dep_params, function(x) x[1]))
bias.df <- data.frame(dim = dimvec, dep = depvec, depnum = depnum)
bias.df$beta0 <- bias.df$beta1 <- bias.df$rho <- bias.df$gamma <- NA
for (i in 1:nrow(results.df)) {
  bias.df$beta0[i] <- results.df$beta_1_0_est[i] - results.df$beta0[i]
  bias.df$beta1[i] <- results.df$beta_1_1_est[i] - results.df$beta1[i]
  bias.df$gamma[i] <- results.df$gamma_1_est[i] - results.df$dep_params[[i]]$gamma_vec
  bias.df$rho[i] <- results.df$rho_1_est[i] - results.df$dep_params[[i]]$rho_vec
}

## Wide to long
library(reshape2)
bias.df <- melt(bias.df, id.vars = list('dim', "dep", "depnum"))

## Prepare faceting
bias.df$dd <- 1
bias.df$dim_f = factor(bias.df$dim, levels=unique(bias.df$dim)[order(unique(bias.df$dim), decreasing = T)])
bias.df$dep_f = factor(bias.df$dep, levels=unique(bias.df$dep)[order(unique(bias.df$depnum), decreasing = F)])
bias.df$variable = factor(bias.df$variable, levels=unique(bias.df$variable)[order(unique(as.character(bias.df$variable)), decreasing = F)])

## Prepare aggregates
agg.dt <- data.table(bias.df)
rmse <- function(x) {sqrt(mean(x^2, na.rm=TRUE))}
agg.df <- data.frame(agg.dt[, j=list(bias = mean(value, na.rm = TRUE), rm = rmse(value)),by = list(dim_f, dep_f, variable, dd)])

### Create separate plots for separate sample sizes
counter <- 1
for (dim in unique(bias.df$dim_f)) {
  
  plot.df <- bias.df[bias.df$dim_f == dim,]
  
  ## Determine y limits
  ylim <- c(min(plot.df$value) - 0.2, max(plot.df$value) + 0.2)
  
  ## Plot errors, bottom row - bias, top row - rmse
  p <- ggplot(plot.df, aes(x=dd, y=value)) +  geom_violin()
  p <- p + geom_hline(yintercept = 0) + stat_summary(fun.y=mean, geom="point", shape=23, size=2)
  p <- p + facet_grid(dep_f ~ variable) + xlab("") + ylab("Error")
  p <- p + theme(axis.title.x=element_blank(),
                 axis.text.x=element_blank(),
                 axis.ticks.x=element_blank(),
                 plot.title = element_text(face = "bold"))
  p <- p + ylim(ylim[1], ylim[2])
  p <- p + ggtitle(as.character(dim))
  p <- p + geom_text(data=agg.df[agg.df$dim_f == dim,], aes(x = dd, y = ylim[1] + 0.1, 
                                      label=sprintf("%0.2f", round(bias, digits = 3))), size=4)
  p <- p + geom_text(data=agg.df[agg.df$dim_f == dim,], aes(x = dd, y = ylim[2] - 0.1, 
                                      label=sprintf("%0.2f", round(rm, digits = 3))), size=4)
  
  filename <- paste0("plots/mc_bias_count_", counter, ".pdf")
  ggsave(filename = filename, plot = p)
  
  counter <- counter + 1
}


###################
# Coverage prob
###################
library(data.table)

## Coverage data frame
coverage.df <- results.df
coverage.df$dim <- unlist(lapply(coverage.df$data_dim, function(x) paste(paste0(c("N = ", "T = "), x[-3]), collapse = ", ")))
coverage.df$dep <- unlist(lapply(results.df$dep_params, function(x) paste(paste0(c("rho = ", "gamma = "), x), collapse = ", ")))
coverage.df$depnum <- unlist(lapply(results.df$dep_params, function(x) x[1]))
coverage.df$rho <- unlist(lapply(coverage.df$dep_params, function(x) x[1]))
coverage.df$gamma <- unlist(lapply(coverage.df$dep_params, function(x) x[2]))

coverage.df$beta0_hit <- coverage.df$beta0 < coverage.df$beta_1_0_hi & coverage.df$beta0 > coverage.df$beta_1_0_lo
coverage.df$beta1_hit <- coverage.df$beta1 < coverage.df$beta_1_1_hi & coverage.df$beta1 > coverage.df$beta_1_1_lo
coverage.df$gamma_hit <- coverage.df$gamma < coverage.df$gamma_1_hi & coverage.df$gamma > coverage.df$gamma_1_lo
coverage.df$rho_hit <- coverage.df$rho < coverage.df$rho_1_hi & coverage.df$rho > coverage.df$rho_1_lo
coverage_var_names <- names(coverage.df)[grep("hit", names(coverage.df))]
coverage.df <- coverage.df[,c("dim", "dep", "depnum", coverage_var_names)]

## Wide to long
coverage.df <- melt(coverage.df, id.vars = c('dim', "dep", "depnum"))
cov.dt <- data.table(coverage.df)

## Renampe
cov.dt$variable <- gsub(pattern = "_hit", replacement = "", cov.dt$variable)

## Aggregate
pt_lo <- function(x) {prop.test(x = sum(x), n = length(x))$conf.int[1]}
pt_hi <- function(x) {prop.test(x = sum(x), n = length(x))$conf.int[2]}
cov.dt <- cov.dt[, .(hit_mean = mean(value), hit_hi = pt_hi(value), hit_lo = pt_lo(value)), by = .(dim, dep, depnum, variable)]

## Prepare faceting
cov.dt$dd <- 1
cov.dt$dim_f = factor(cov.dt$dim, levels=unique(cov.dt$dim)[order(unique(cov.dt$dim), decreasing = T)])
cov.dt$dep_f = factor(cov.dt$dep, levels=unique(cov.dt$dep)[order(unique(cov.dt$depnum), decreasing = F)])

## Bar plot for each sample size
cov.df <- as.data.frame(cov.dt)
counter <- 1
for (dim in unique(cov.df$dim_f)) {
  
  plot.df <- cov.df[cov.df$dim_f == dim,]
  
  p <- ggplot(data=plot.df, aes(x=dd))
  p <- p + geom_bar(aes(weight = hit_mean))
  p <- p + geom_errorbar(aes(ymin=hit_lo, ymax=hit_hi), width = 0.2)
  p <- p + facet_grid(dep_f ~ variable) + geom_hline(yintercept = 0.9, linetype = 2, size = 0.25)
  p <- p + theme(axis.title.x=element_blank(),
                 axis.text.x=element_blank(),
                 axis.ticks.x=element_blank(),
                 plot.title = element_text(face = "bold"))
  p <- p + ggtitle(as.character(dim))
  p <- p + ylab("Coverage")
  
  filename <- paste0("plots/mc_coverage_count_", counter, ".pdf")
  ggsave(filename = filename, plot = p)
  
  counter <- counter + 1
}


#################################################
# II. EVALUATE BIAS / RMSE FOR BINARY DATA
# Note: Needs M = 300 for SE estimation
#################################################

# BAUSTELLE
# PROBLEM: STANDARD ERRORS TOO NARROW IF RHO/GAMMA = 0.45
# POSSIBLE EXPLANATION: maxiter NOT HIGH ENOUG DURING ESTIMATION

###################
# Set up configurations
###################

## Constants
M <- 50

## Configs
config.df <- expand.grid(data_dim = list(list(N = 36, TT = 10, G = 1)),
                         dep_params = list(list(rho_vec = 0.45, gamma_vec = 0.45)),
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
saveRDS(results.df, file = 'results/mc_spatiotemporal_rmse_binary.rds')


###################
# Print average run-times
###################

runtime.df <- results.df
runtime.df$dim <- unlist(lapply(runtime.df$data_dim, function(x) paste(paste0(c("N = ", "T = ", "G = "), x), collapse = ", ")))
runtime.dt <- data.table(runtime.df)
runtime.dt <- runtime.dt[,list(runtime = mean(elapsed)), by = list(dim)]
print(runtime.dt$runtime)
print(runtime.dt$runtime/60)

###################
# Plot bias / RMSE
###################

library(ggplot2)
library(data.table)

## Calculate bias
dimvec <- unlist(lapply(results.df$data_dim, function(x) paste(paste0(c("N = ", "T = "), x[-3]), collapse = ", ")))
depvec <- unlist(lapply(results.df$dep_params, function(x) paste(paste0(c("rho = ", "gamma = "), x), collapse = ", ")))
depnum <- unlist(lapply(results.df$dep_params, function(x) x[1]))
bias.df <- data.frame(dim = dimvec, dep = depvec, depnum = depnum)
bias.df$beta0 <- bias.df$beta1 <- bias.df$rho <- bias.df$gamma <- NA
for (i in 1:nrow(results.df)) {
  bias.df$beta0[i] <- results.df$beta_1_0_est[i] - results.df$beta0[i]
  bias.df$beta1[i] <- results.df$beta_1_1_est[i] - results.df$beta1[i]
  bias.df$gamma[i] <- results.df$gamma_1_est[i] - results.df$dep_params[[i]]$gamma_vec
  bias.df$rho[i] <- results.df$rho_1_est[i] - results.df$dep_params[[i]]$rho_vec
}

## Wide to long
library(reshape2)
bias.df <- melt(bias.df, id.vars = list('dim', "dep", "depnum"))

## Prepare faceting
bias.df$dd <- 1
bias.df$dim_f = factor(bias.df$dim, levels=unique(bias.df$dim)[order(unique(bias.df$dim), decreasing = T)])
bias.df$dep_f = factor(bias.df$dep, levels=unique(bias.df$dep)[order(unique(bias.df$depnum), decreasing = F)])
bias.df$variable = factor(bias.df$variable, levels=unique(bias.df$variable)[order(unique(as.character(bias.df$variable)), decreasing = F)])

## Prepare aggregates
agg.dt <- data.table(bias.df)
rmse <- function(x) {sqrt(mean(x^2, na.rm=TRUE))}
agg.df <- data.frame(agg.dt[, j=list(bias = mean(value, na.rm = TRUE), rm = rmse(value)),by = list(dim_f, dep_f, variable, dd)])

### Create separate plots for separate sample sizes
counter <- 1
for (dim in unique(bias.df$dim_f)) {
  
  plot.df <- bias.df[bias.df$dim_f == dim,]
  
  ## Determine y limits
  ylim <- c(min(plot.df$value) - 0.2, max(plot.df$value) + 0.2)
  
  ## Plot errors, bottom row - bias, top row - rmse
  p <- ggplot(plot.df, aes(x=dd, y=value)) +  geom_violin()
  p <- p + geom_hline(yintercept = 0) + stat_summary(fun.y=mean, geom="point", shape=23, size=2)
  p <- p + facet_grid(dep_f ~ variable) + xlab("") + ylab("Error")
  p <- p + theme(axis.title.x=element_blank(),
                 axis.text.x=element_blank(),
                 axis.ticks.x=element_blank(),
                 plot.title = element_text(face = "bold"))
  p <- p + ylim(ylim[1], ylim[2])
  p <- p + ggtitle(as.character(dim))
  p <- p + geom_text(data=agg.df[agg.df$dim_f == dim,], aes(x = dd, y = ylim[1] + 0.1, 
                                                            label=sprintf("%0.2f", round(bias, digits = 3))), size=4)
  p <- p + geom_text(data=agg.df[agg.df$dim_f == dim,], aes(x = dd, y = ylim[2] - 0.1, 
                                                            label=sprintf("%0.2f", round(rm, digits = 3))), size=4)
  
  filename <- paste0("plots/mc_bias_binary_", counter, ".pdf")
  ggsave(filename = filename, plot = p)
  
  counter <- counter + 1
}


###################
# Coverage prob
###################
library(data.table)

## Coverage data frame
coverage.df <- results.df
coverage.df$dim <- unlist(lapply(coverage.df$data_dim, function(x) paste(paste0(c("N = ", "T = "), x[-3]), collapse = ", ")))
coverage.df$dep <- unlist(lapply(results.df$dep_params, function(x) paste(paste0(c("rho = ", "gamma = "), x), collapse = ", ")))
coverage.df$depnum <- unlist(lapply(results.df$dep_params, function(x) x[1]))
coverage.df$rho <- unlist(lapply(coverage.df$dep_params, function(x) x[1]))
coverage.df$gamma <- unlist(lapply(coverage.df$dep_params, function(x) x[2]))

coverage.df$beta0_hit <- coverage.df$beta0 < coverage.df$beta_1_0_hi & coverage.df$beta0 > coverage.df$beta_1_0_lo
coverage.df$beta1_hit <- coverage.df$beta1 < coverage.df$beta_1_1_hi & coverage.df$beta1 > coverage.df$beta_1_1_lo
coverage.df$gamma_hit <- coverage.df$gamma < coverage.df$gamma_1_hi & coverage.df$gamma > coverage.df$gamma_1_lo
coverage.df$rho_hit <- coverage.df$rho < coverage.df$rho_1_hi & coverage.df$rho > coverage.df$rho_1_lo
coverage_var_names <- names(coverage.df)[grep("hit", names(coverage.df))]
coverage.df <- coverage.df[,c("dim", "dep", "depnum", coverage_var_names)]

## Wide to long
coverage.df <- melt(coverage.df, id.vars = c('dim', "dep", "depnum"))
cov.dt <- data.table(coverage.df)

## Renampe
cov.dt$variable <- gsub(pattern = "_hit", replacement = "", cov.dt$variable)

## Aggregate
pt_lo <- function(x) {prop.test(x = sum(x), n = length(x))$conf.int[1]}
pt_hi <- function(x) {prop.test(x = sum(x), n = length(x))$conf.int[2]}
cov.dt <- cov.dt[, .(hit_mean = mean(value), hit_hi = pt_hi(value), hit_lo = pt_lo(value)), by = .(dim, dep, depnum, variable)]

## Prepare faceting
cov.dt$dd <- 1
cov.dt$dim_f = factor(cov.dt$dim, levels=unique(cov.dt$dim)[order(unique(cov.dt$dim), decreasing = T)])
cov.dt$dep_f = factor(cov.dt$dep, levels=unique(cov.dt$dep)[order(unique(cov.dt$depnum), decreasing = F)])

## Bar plot for each sample size
cov.df <- as.data.frame(cov.dt)
counter <- 1
for (dim in unique(cov.df$dim_f)) {
  
  plot.df <- cov.df[cov.df$dim_f == dim,]
  
  p <- ggplot(data=plot.df, aes(x=dd))
  p <- p + geom_bar(aes(weight = hit_mean))
  p <- p + geom_errorbar(aes(ymin=hit_lo, ymax=hit_hi), width = 0.2)
  p <- p + facet_grid(dep_f ~ variable) + geom_hline(yintercept = 0.9, linetype = 2, size = 0.25)
  p <- p + theme(axis.title.x=element_blank(),
                 axis.text.x=element_blank(),
                 axis.ticks.x=element_blank(),
                 plot.title = element_text(face = "bold"))
  p <- p + ggtitle(as.character(dim))
  p <- p + ylab("Coverage")
  
  filename <- paste0("plots/mc_coverage_binary_", counter, ".pdf")
  ggsave(filename = filename, plot = p)
  
  counter <- counter + 1
}




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
  
  model$rho_vec[] <- 0.25
  if (temporal_dep) {
    model$gamma_vec[] <- 0.25
  }
  if (outcome_dep) {
    model$lambda_vec[] <- 0.25
  }
  model$beta.ls <- rep(list(c(2 ,1)), data$G)
  
  timing <- system.time({
    res <- try({
      set.seed(0)
      for (i in 1:5) {
        ll <- model$E_llik(theta = model$pack_theta())
      }
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
M <- 1

## Configs
config.df <- expand.grid(data_dim = list(list(N = 64, TT = 1, G = 1),
                                         list(N = 256, TT = 1, G = 1),
                                         list(N = 400, TT = 1, G = 1), 
                                         list(N = 576, TT = 1, G = 1), 
                                         list(N = 784, TT = 1, G = 1), 
                                         list(N = 1024, TT = 1, G = 1), 
                                         
                                         list(N = 64, TT = 1, G = 4),
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

results.ls <- foreach(config = iter(config.ls), .options.multicore = mcoptions) %do% {
  
  print(paste(config$iter, '/', nrow(config.df)))
  
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
library(gridExtra)

scal.df <- results.df
scal.df$N <- unlist(lapply(scal.df$data_dim, function(x) x$N))
scal.df$TT <- unlist(lapply(scal.df$data_dim, function(x) x$TT))
scal.df$G <- unlist(lapply(scal.df$data_dim, function(x) x$G))
scal.dt <- data.table(scal.df[scal.df$TT == 1 & scal.df$G <= 4,])
scal.dt <- scal.dt[, j=list(elapsed = mean(elapsed, na.rm = TRUE)), by = list(N, G)]

scal.dt$Glab <- paste('G =', scal.dt$G)

## Plot for G = 1 and G > 1
p1 <- ggplot(data=scal.dt[scal.dt$G == 1,], aes(x=N, y=elapsed)) +
  geom_line() +
  geom_point() + xlab("N") +
  facet_grid( ~ Glab) +
  xlab("N") + ylab("Evaluation Time (seconds)")
p2 <- ggplot(data=scal.dt[scal.dt$G == 4,], aes(x=N, y=elapsed)) +
  geom_line() +
  geom_point() + xlab("N") +
  facet_grid( ~ Glab) +
  xlab("N") + ylab("Evaluation Time (seconds)")
p <- arrangeGrob(p1, p2, nrow=1, ncol=2) 
ggsave(filename = "plots/mc_full_time_N.pdf", plot = p, scale = 0.75)


###################
# Plot Scalability in TT
###################

scal.df <- results.df
scal.df$N <- unlist(lapply(scal.df$data_dim, function(x) x$N))
scal.df$TT <- unlist(lapply(scal.df$data_dim, function(x) x$TT))
scal.df$G <- unlist(lapply(scal.df$data_dim, function(x) x$G))
scal.dt <- data.table(scal.df[scal.df$N == 64 & scal.df$G <= 4,])
scal.dt <- scal.dt[, j=list(elapsed = mean(elapsed, na.rm = TRUE)), by = list(TT, G)]

scal.dt$Glab <- paste('G =', scal.dt$G)

## Plot for G = 1 and G > 1
p1 <- ggplot(data=scal.dt[scal.dt$G == 1,], aes(x=TT, y=elapsed)) +
  geom_line() +
  geom_point() + 
  facet_grid( ~ Glab) +
  xlab("T") + ylab("Evaluation Time (seconds)")
p2 <- ggplot(data=scal.dt[scal.dt$G == 4,], aes(x=TT, y=elapsed)) +
  geom_line() +
  geom_point() + 
  facet_grid( ~ Glab) +
  xlab("T") + ylab("Evaluation Time (seconds)")
p <- arrangeGrob(p1, p2, nrow=1, ncol=2) 
ggsave(filename = "plots/mc_full_time_T.pdf", plot = p, scale = 0.75)

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

