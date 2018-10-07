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
      model$train(maxiter = 100, M = 50, abs_tol = 1e-5, burnin = 0, thinning = 1, verbose = TRUE, soft_init = FALSE)
      model$compute_vcov(M = 300, thinning = 3)
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
# Note: Needs M = 300, thinning = 3 for SE estimation
#       Needs maxiter = 100 
#################################################

# BAUSTELLE
# PROBLEM: STANDARD ERRORS TOO NARROW IF RHO/GAMMA = 0.45
# POSSIBLE EXPLANATION: maxiter NOT HIGH ENOUGH DURING ESTIMATION
# ACTUAL EXPLANATION: RMSE TOO LARGE (AND SE TOO LOW) IF OUTCOME IS HIGHLY UNBALANCED BETWEEN 0/1

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
                         beta0 = c(0),
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

  if (config$TT == 1) {
    gamma_vec <- 0
  } else {
    gamma_vec <- 0.25
  }
  rho_vec <- 0.25
  
  set.seed(1)
  beta <- c(1, 1)
  data <- simulate_data(N = config$N, TT = config$TT, G = 1L, count = TRUE, 
                        rho_vec = rho_vec, lambda_vec = 0, gamma_vec = gamma_vec,
                        beta.ls = list(beta), sigma2_vec = 1, 
                        X_params = c(0, 1))
  data$beta.ls <- list(beta)
  return(data)
}

###################
# Data estimation function
###################

fit_dimstar <- function(data) {
  
  outcome_dep <- FALSE
  if (data$TT == 1) {
    temporal_dep <- FALSE
  } else {
    temporal_dep <- TRUE
  }
  
  model <- DISTAR$new(X.ls = data$X.ls, y.ls = data$y.ls, W_t = data$W_t, 
                      N = data$N, G = 1, TT = data$TT,
                      count = data$count, 
                      spatial_dep = TRUE, 
                      temporal_dep = temporal_dep, 
                      outcome_dep = outcome_dep)
  model$sample_ystar(M = 50, ystar_init = colMeans(model$ystar_sample))
  
  model$rho_vec[] <- 0.25
  if (temporal_dep) {
    model$gamma_vec[] <- 0.25
  }

  model$beta.ls <- data$beta.ls
  
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
Nconfig.df <- expand.grid(N = seq(16, 64, 4)^2, TT = c(1, 4), experiment = "N")
Tconfig.df <- expand.grid(N = c(256, 1024), TT = seq(5, 50, 5), experiment = "T")
config.df <- rbind(Nconfig.df, Tconfig.df)
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
  results <- fit_dimstar(data)

  ## Return result
  results
}

###################
# Concatenate & save results
###################

results.df <- do.call('rbind', lapply(results.ls, function(x) as.data.frame(x[length(x)])))
results.df <- cbind(config.df, results.df)
saveRDS(results.df, file = 'results/mc_full_count.rds')


###################
# Plot Scalability in N
###################
library(data.table)
library(ggplot2)
library(gridExtra)

scal.df <- results.df
scal.dt <- data.table(scal.df[scal.df$experiment == "N",])
scal.dt <- scal.dt[, j=list(elapsed = mean(elapsed, na.rm = TRUE)), by = list(N, TT, experiment)]
scal.dt$TTlab <- paste('T =', scal.dt$TT)

## Plot for TT = 1 and TT > 1
p <- ggplot(data=scal.dt, aes(x=N, y=elapsed)) +
  geom_line() +
  geom_point() + xlab("N") +
  facet_grid( ~ TTlab) +
  xlab("N") + ylab("Evaluation Time (seconds)")
ggsave(filename = "plots/mc_time_N.pdf", plot = p, scale = 0.5)


###################
# Plot Scalability in TT
###################

scal.df <- results.df
scal.dt <- data.table(scal.df[scal.df$experiment == "T",])
scal.dt <- scal.dt[, j=list(elapsed = mean(elapsed, na.rm = TRUE)), by = list(N, TT, experiment)]
scal.dt$Nlab <- paste('N =', scal.dt$N)
scal.dt$Nlab <- factor(x = scal.dt$Nlab, levels = unique(scal.dt$Nlab)[order(unique(scal.dt$N))])

## Plot
p <- ggplot(data=scal.dt, aes(x=TT, y=elapsed)) +
  geom_line() +
  geom_point() + xlab("T") +
  facet_grid( ~ Nlab) +
  xlab("T") + ylab("Evaluation Time (seconds)")
ggsave(filename = "plots/mc_time_T.pdf", plot = p, scale = 0.5)

