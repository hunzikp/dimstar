##############################################################
# Tests with ADMIN UNIT / ACLED data
##############################################################

library(pgrid)
library(cshapes)
library(raster)
library(velox)
library(dimstar)
library(rgeos)
library(Matrix)
library(data.table)
library(plyr)

###########################
# Prepare the full data set
###########################

## Load the data
prio.ls <- readRDS('data/admin_acled.rds')

## Get data
panel.df <- prio.ls$panel
panel.df$gid <- panel.df$ADM1_CODE
gid.spdf <- prio.ls$spdf

## Plot the polygons
plot(gid.spdf)

## Get panel data, number of units, time-periods
N <- length(unique(panel.df$gid))
TT <- nrow(panel.df)/N

## Make lagged covariates
panel.df <- panel.df[order(panel.df$ADM1_CODE, panel.df$year),]
gid.spdf <- gid.spdf[order(gid.spdf$ADM1_CODE),]
vars <- c("battle_event", "civilian_event", "protest_event")
for (v in vars) {
  units <- unique(panel.df$gid)
  panel.df$newvar <- NA
  for (u in units) {
    target_var <- panel.df[panel.df$gid == u,v]
    panel.df$newvar[panel.df$gid == u] <- c(NA, target_var[1:(TT-1)])
  }
  names(panel.df)[ncol(panel.df)] <- paste0(v, "_lag")
}
panel.df <- panel.df[order(panel.df$year, panel.df$gid),]
panel.df <- na.omit(panel.df)
TT <- nrow(panel.df)/N

## Make the spatial weights matrix
B <- gIntersects(gid.spdf, byid = TRUE)
diag(B) <- FALSE
B <- B*1
rs <- rowSums(B)
rs <- ifelse(rs == 0, 1, rs)
W_t <- B / rs
W_t <- Matrix(W_t, sparse = TRUE)

## Make the panel model matrix, normalize variables
formula <- battle_event ~ log(Shape_Area+1) + log(pop+1) + log(ttime) + mountains # + log(battle_event_lag + 1)
X <- model.matrix(formula, data = panel.df)
for (k in 2:ncol(X)) {
  X[,k] <- (X[,k] - mean(X[,k]))/sd(X[,k])
}

## Summarize data
print(N)
print(TT)


###########################
# Fit DIMSTAR to battle events, no dependence
###########################

X.ls <- list(X)
y.ls <- list(panel.df$battle_event)

model1 <- DISTAR$new(X.ls = X.ls, y.ls = y.ls, W_t = W_t, 
                    N = N, G = 1, TT = TT,
                    count = TRUE, 
                    spatial_dep = FALSE, 
                    temporal_dep = FALSE, 
                    outcome_dep = FALSE)

ch_vec <- model1$train(maxiter = 15, M = 50, abs_tol = 1e-4, soft_init = FALSE)
plot(ch_vec)
model1$compute_vcov(M = 100)
model1$coeftest()


###########################
# Fit DIMSTAR to battle events, spatio-temporal dependence
###########################

X.ls <- list(X)
y.ls <- list(panel.df$battle_event)

model2 <- DISTAR$new(X.ls = X.ls, y.ls = y.ls, W_t = W_t, 
                     N = N, G = 1, TT = TT,
                     count = TRUE, 
                     spatial_dep = TRUE, 
                     temporal_dep = TRUE, 
                     outcome_dep = FALSE)

ch_vec <- model2$train(maxiter = 15, M = 50, abs_tol = 1e-4, soft_init = FALSE)
plot(ch_vec)
model2$compute_vcov(M = 100)
model2$coeftest()


###########################
# Predict battle events during last 5 years
###########################

years <- unique(panel.df$year)
predict_years <- tail(years, 5)

predictions_ls <- vector('list', length(predict_years))
for (i in 1:length(predictions_ls)) {
  
  predict_year <- predict_years[i]

  ## Forecast with spatio-temporal dimstar
  # Subset data
  panel_idx <- panel.df$year <= predict_year
  X_this <- X[panel_idx,]
  y_this <- panel.df$battle_event[panel_idx]
  TT_this <- length(y_this)/N
  years_this <- panel.df$year[panel_idx]
  
  # Prepare prediction
  predict_idx <- years_this == predict_year
  ytrue <- y_this[predict_idx]
  y_this[predict_idx] <- NA
  
  # Fit
  this_model <- DISTAR$new(X.ls = list(X_this), y.ls = list(y_this), W_t = W_t, 
                       N = N, G = 1, TT = TT_this,
                       count = TRUE, 
                       spatial_dep = TRUE, 
                       temporal_dep = TRUE, 
                       outcome_dep = FALSE)
  ch_vec <- this_model$train(maxiter = 15, M = 50, abs_tol = 1e-4, soft_init = TRUE)
  
  # Predict
  this_model$sample_ystar(M = 500, ystar_init = colMeans(this_model$ystar_sample))
  ypred <- colMeans(exp(this_model$ystar_sample))[predict_idx]
  
  
  ## Forecast with lagged-y regression
  train.df <- panel.df[panel.df$year < predict_year,]
  glm.fit <- glm.nb(formula = battle_event ~ log(Shape_Area+1) + log(pop+1) + log(ttime) + mountains + log(battle_event_lag+1), 
                    data = train.df)
  predict.df <- panel.df[panel.df$year == predict_year,]
  ypred_glm <- predict(glm.fit, newdata = predict.df, type = 'response')
  
  ## Store results
  res <- list(ytrue = ytrue, ypred = ypred, ypred_glm = ypred_glm)
  predictions_ls[[i]] <- res
  
  print(predict_year)
  flush.console()
}

## Calc stats
rmse <- function(x, y) {sqrt(mean((x-y)^2))}
mae <- function(x, y) {mean(abs(x-y))}
stats.ls <- vector('list', length(predictions_ls))
for (i in 1:length(predictions_ls)) {
  df = data.frame(year = predict_years[i])
  res <- predictions_ls[[i]]
  df$rmse <- rmse(res$ytrue, res$ypred)
  df$rmse_glm <- rmse(res$ytrue, res$ypred_glm)
  df$mae <- mae(res$ytrue, res$ypred)
  df$mae_glm <- mae(res$ytrue, res$ypred_glm)
  stats.ls[[i]] <- df
}
stats.df <- do.call('rbind', stats.ls)









###########################
# Fit DIMSTAR to full dataset, no dependence!
###########################

X.ls <- list(X, X)
y.ls <- list(panel.df$battle_event, panel.df$civilian_event)

model <- DISTAR$new(X.ls = X.ls, y.ls = y.ls, W_t = W_t, 
                    N = N, G = 2, TT = TT,
                    count = TRUE, 
                    spatial_dep = FALSE, 
                    temporal_dep = FALSE, 
                    outcome_dep = FALSE)

ch_vec <- model$train(maxiter = 15, M = 50, abs_tol = 1e-4, soft_init = FALSE)
plot(ch_vec)
model$compute_vcov(M = 100)
model$coeftest()


###########################
# Fit DIMSTAR to full dataset, with dep
###########################

X.ls <- list(X, X)
y.ls <- list(panel.df$battle_event, panel.df$civilian_event)

model1 <- DISTAR$new(X.ls = X.ls, y.ls = y.ls, W_t = W_t, 
                     N = N, G = 2, TT = TT,
                     count = TRUE, 
                     spatial_dep = TRUE, 
                     temporal_dep = TRUE, 
                     outcome_dep = TRUE)

ch_vec <- model1$train(maxiter = 15, M = 50, abs_tol = 1e-4)
model1$compute_vcov(M = 100)
model1$coeftest()


###########################
# Forecast last period, only battle
###########################

X.ls <- list(X)
y.ls <- list(panel.df$battle_event)
last_idx <- panel.df$year == max(panel.df$year)
y.ls[[1]][last_idx] <- NA

model <- DISTAR$new(X.ls = X.ls, y.ls = y.ls, W_t = W_t, 
                    N = N, G = 1, TT = TT,
                    count = TRUE, 
                    spatial_dep = TRUE, 
                    temporal_dep = TRUE, 
                    outcome_dep = FALSE)

ch_vec <- model$train(maxiter = 15, M = 50, abs_tol = 1e-4, soft_init = TRUE)
plot(ch_vec)
model$compute_vcov(M = 100)
model$coeftest()

model$sample_ystar(M = 1000, ystar_init = colMeans(model$ystar_sample))
ypred <- colMeans(exp(model$ystar_sample))[last_idx]
ytrue <- panel.df$battle_event[last_idx]
plot(ypred, ytrue)
mean(abs(ypred - ytrue))
# [1] 5.28039
# [1] 5.175352



###########################
# Forecast last period
###########################

X.ls <- list(X, X)
y.ls <- list(panel.df$battle_event, panel.df$civilian_event)
last_idx <- panel.df$year == max(panel.df$year)
y.ls[[1]][last_idx] <- NA
y.ls[[2]][last_idx] <- NA

model3 <- DISTAR$new(X.ls = X.ls, y.ls = y.ls, W_t = W_t, 
                     N = N, G = 2, TT = TT,
                     count = TRUE, 
                     spatial_dep = TRUE, 
                     temporal_dep = TRUE, 
                     outcome_dep = FALSE)

ch_vec <- model3$train(maxiter = 15, M = 50, abs_tol = 1e-4)
model3$compute_vcov(M = 100)
model3$coeftest()

model3$sample_ystar(M = 100, ystar_init = colMeans(model3$ystar_sample))
ystar_sample <- model3$ystar_sample
ystar_predict <- tail(colMeans(exp(ystar_sample)), N*2)

battle_true <- panel.df$battle_event[last_idx]
civilian_true <- panel.df$civilian_event[last_idx]
battle_predict <- ystar_predict[1:N]
civilian_predict <- ystar_predict[(N+1):(2*N)]

plot(battle_true, battle_predict)
plot(civilian_true, civilian_predict)

mean(abs(battle_predict - battle_true))
# 6.374384
mean(abs(civilian_true - civilian_predict))
# 5.995947


