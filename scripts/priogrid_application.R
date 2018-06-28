##############################################################
# Tests with PRIO-GRID / ACLED data
##############################################################

library(pgrid)
library(cshapes)
library(raster)
library(velox)
library(dimstar)
library(rgeos)
library(Matrix)
library(data.table)

###########################
# Prepare the full data set
###########################

## Load the data
prio.ls <- readRDS('data/prio_acled.rds')

## Get data
panel.df <- prio.ls$panel

## Determine target
panel.dt <- data.table(panel.df[panel.df$year == max(panel.df$year),])
agg.dt <- panel.dt[, j=list(cells = length(gid), bev = sum(battle_event), cev = sum(civilian_event)), by = list(gwno)]
agg.dt <- agg.dt[order(agg.dt$bev),]

## Set target country
TARGET <- 482

## Subset to target country

panel.df <- panel.df[panel.df$gwno == TARGET,]
gid.spdf <- prio.ls$spdf
gid.spdf <- gid.spdf[gid.spdf$gid %in% panel.df$gid,]

## Plot the prio-grid cells
plot(gid.spdf)

## Make the spatial weights matrix
B <- gIntersects(gid.spdf, byid = TRUE)
diag(B) <- FALSE
B <- B*1
rs <- rowSums(B)
rs <- ifelse(rs == 0, 1, rs)
W_t <- B / rs
W_t <- Matrix(W_t, sparse = TRUE)

## Get panel data, number of units, time-periods
N <- length(unique(panel.df$gid))
TT <- nrow(panel.df)/N

## Make the panel model matrix, normalize variables
formula <- battle_event ~ log(pop_gpw_sum+1) + mountains_mean + log(ttime_mean+1) + log(nlights_mean+1) + excluded
X <- model.matrix(formula, data = panel.df)
for (k in 2:ncol(X)) {
  X[,k] <- (X[,k] - mean(X[,k]))/sd(X[,k])
}

## Summarize data
print(N)
print(TT)

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

ypred <- colMeans(exp(model$ystar_sample))[last_idx]
ytrue <- panel.df$battle_event[last_idx]
plot(ypred, ytrue)
mean(abs(ypred - ytrue))
# [1] 1.40929


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
                    outcome_dep = TRUE)

ch_vec <- model3$train(maxiter = 15, M = 50, abs_tol = 1e-4)
model3$compute_vcov(M = 100)
model3$coeftest()

ypred <- colMeans(model3$ystar_sample)
y <- model3$y

oos_pred <- ypred[is.na(y)]
oos_y <- c(panel.df$battle_event[last_idx], panel.df$civilian_event[last_idx])

ystar_sample <- model3$ystar_sample
ystar_predict <- tail(colMeans(ystar_sample), N*2)

battle_true <- panel.df$battle_event[last_idx]
civilian_true <- panel.df$civilian_event[last_idx]
battle_predict <- exp(ystar_predict[1:N])
civilian_predict <- exp(ystar_predict[(N+1):(2*N)])

plot(battle_true, battle_predict)
plot(civilian_true, civilian_predict)



