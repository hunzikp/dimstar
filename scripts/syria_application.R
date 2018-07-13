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
library(blockmatrix)
library(data.table)

###########################
# Prepare the full data set
###########################

## Load the data
prio.ls <- readRDS('extdata/syria_acled.rds')

## Get data
panel.df <- prio.ls$panel
gid.spdf <- prio.ls$spdf

## Plot the polygons
plot(gid.spdf)

## Get panel data, number of units, time-periods
N <- length(unique(panel.df$gid))
TT <- nrow(panel.df)/N

## Make lagged covariates
panel.df <- panel.df[order(panel.df$gid, panel.df$startdate),]
gid.spdf <- gid.spdf[order(gid.spdf$gid),]
vars <- c("battle_event", "civilian_event", "protest_event", "isis_strike", "strike_event")
for (v in vars) {
  units <- unique(panel.df$gid)
  panel.df$newvar <- NA
  for (u in units) {
    target_var <- panel.df[panel.df$gid == u,v]
    panel.df$newvar[panel.df$gid == u] <- c(NA, target_var[1:(TT-1)])
  }
  names(panel.df)[ncol(panel.df)] <- paste0(v, "_lag")
}
panel.df <- panel.df[order(panel.df$startdate, panel.df$gid),]
gid.spdf <- gid.spdf[order(gid.spdf$gid),]
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
formula <- battle_event ~ log(Shape_Area+1) + log(pop+1) + log(ttime) + mountains + log(isis_strike_lag + 1)
X <- model.matrix(formula, data = panel.df)
# for (k in 2:ncol(X)) {
#   X[,k] <- (X[,k] - mean(X[,k]))/sd(X[,k])
# }

## Summarize data
print(N)
print(TT)


###########################
# Plot counts
###########################
library(broom)
library(ggplot2)
library(dtplyr)
library(dplyr)
library(magrittr)
library(viridis)

## Prepare counts
panel.dt <- data.table(panel.df)
agg.dt <- panel.dt[, j=list(battles = sum(isis_battle), civ = sum(isis_civilian), strike = sum(isis_strike)), by = list(gid)]
agg.df <- as.data.frame(agg.dt)

## Prepare polygons
spdf_fortified <- tidy(gid.spdf, region = "gid")
ggplot() +
  geom_polygon(data = spdf_fortified, aes( x = long, y = lat, group = group)) +
  theme_void() +
  coord_map()

## Merge
spdf_fortified =  merge(spdf_fortified , agg.df, by.x = "id", by.y = "gid", sort = FALSE) 
spdf_fortified$strike1 <- spdf_fortified$strike + 1
spdf_fortified$battles1 <- spdf_fortified$battles + 1
spdf_fortified$civ1 <- spdf_fortified$civ + 1

## Strike plot
p_strike <- ggplot() +
  geom_polygon(data = spdf_fortified, aes(fill = strike1, x = long, y = lat, group = group) , size=0, alpha=0.9) +
  theme_void() +
  scale_fill_viridis(trans = "log", breaks=c(5,10,20,50,100), 
                     name="Airstrikes against IS", 
                     guide = guide_legend( keyheight = unit(3, units = "mm"), 
                                           keywidth=unit(12, units = "mm"), 
                                           label.position = "bottom", 
                                           title.position = 'top', nrow=1), 
                     option = "inferno") +
  theme(
    text = element_text(color = "#22211d"), 
    plot.background = element_rect(fill = "#f5f5f2", color = NA), 
    panel.background = element_rect(fill = "#f5f5f2", color = NA), 
    legend.background = element_rect(fill = "#f5f5f2", color = NA),
    legend.position = c(0.7, 0.09)
  ) +
  coord_map()


## Battle plot
p_battle <- ggplot() +
  geom_polygon(data = spdf_fortified, aes(fill = battles1, x = long, y = lat, group = group) , size=0, alpha=0.9) +
  theme_void() +
  scale_fill_viridis(trans = "log", breaks=c(5,10,50,100,500), 
                     name="Battles involving IS", 
                     guide = guide_legend( keyheight = unit(3, units = "mm"), 
                                           keywidth=unit(12, units = "mm"), 
                                           label.position = "bottom", 
                                           title.position = 'top', nrow=1), 
                     option = "inferno") +
  theme(
    text = element_text(color = "#22211d"), 
    plot.background = element_rect(fill = "#f5f5f2", color = NA), 
    panel.background = element_rect(fill = "#f5f5f2", color = NA), 
    legend.background = element_rect(fill = "#f5f5f2", color = NA),
    legend.position = c(0.7, 0.09)
  ) +
  coord_map()


## Civ plot
p_civ <- ggplot() +
  geom_polygon(data = spdf_fortified, aes(fill = civ1, x = long, y = lat, group = group) , size=0, alpha=0.9) +
  theme_void() +
  scale_fill_viridis(trans = "log", breaks=c(5,10,20,50,100), 
                     name="IS violence against civilians", 
                     guide = guide_legend( keyheight = unit(3, units = "mm"), 
                                           keywidth=unit(12, units = "mm"), 
                                           label.position = "bottom", 
                                           title.position = 'top', nrow=1), 
                     option = "inferno") +
  theme(
    text = element_text(color = "#22211d"), 
    plot.background = element_rect(fill = "#f5f5f2", color = NA), 
    panel.background = element_rect(fill = "#f5f5f2", color = NA), 
    legend.background = element_rect(fill = "#f5f5f2", color = NA),
    legend.position = c(0.7, 0.09)
  ) +
  coord_map()

## Prepare time series
ts.dt <- panel.dt[, j=list(battles = sum(isis_battle), civ = sum(isis_civilian), strike = sum(isis_strike)), by = list(startdate)]
ts.df <- melt(as.data.frame(ts.dt), id="startdate")
levels(ts.df$variable)[levels(ts.df$variable)=="battles"] <- "Battles"
levels(ts.df$variable)[levels(ts.df$variable)=="strike"] <- "Airstrikes"
levels(ts.df$variable)[levels(ts.df$variable)=="civ"] <- "Violence agst. Civilians"
p_ts <- ggplot(ts.df, aes(x=startdate, y=value, group=variable, linetype=variable)) + geom_line() + ylab("Event Count") + xlab("") +
  theme(
    legend.position="bottom",
    text = element_text(color = "#22211d"), 
    plot.background = element_rect(fill = "#f5f5f2", color = NA), 
    panel.background = element_rect(fill = "#f5f5f2", color = NA), 
    legend.background = element_rect(fill = "#f5f5f2", color = NA),
    legend.text = element_text(size=15)
    ) + labs(linetype='') 


## Write out plot
require(gridExtra)
g <- arrangeGrob(p_strike, p_battle, p_civ, p_ts, layout_matrix = matrix(c(1,2,3,4,4,4), 2, 3, byrow = T))
ggsave(file="plots/application.pdf", g, height = 8.5*0.9, width = 14*0.9)




###########################
# Function to calculate elasticities
###########################

## Function for univariate outcome models
get_univ_elast <- function(model) {
  
  N <- model$N
  TT <- model$TT
  beta_est <- model$beta.ls[[1]]
  K <- length(beta_est)
  beta_names <- model$get_theta_names()[1:K]
  invA <- solve(model$get_A())
  
  # Direct elasticity; x_it on y_it
  multiplier <- diag(invA)
  elasticity_est <- rep(NA, K-1)
  for (i in 2:K) {
    elasticity_est[i-1] <- mean(multiplier*beta_est[i])
  }
  names(elasticity_est) <- beta_names[-1]
  
  # Indirect elasticity x_it on sum_j y_jt
  inv_blocks <- as.blockmatrix(as.matrix(invA), nrowe = N, ncole = N)
  indirect_elasticity_est <- rep(NA, K-1)
  for (k in 2:K) {
    beta_k <- beta_est[k]
    elast_k <- rep(NA, N)
    for (i in 1:N) {
      elast_k_i <- rep(NA, TT)
      for (t in 1:TT) {
        invA_block <- inv_blocks[t,t]
        multiplier <- invA_block[-i,i]
        elast_k_i[t] <- sum(multiplier)*beta_k
      }
      elast_k[i] <- mean(elast_k_i)
    }
    indirect_elasticity_est[k-1] <- mean(elast_k)
  }
  names(indirect_elasticity_est) <- beta_names[-1]
  
  out <- list(direct_elast = elasticity_est, indirect_elast = indirect_elasticity_est)
  return(out)
}


## Function for multivariate outcome models with 2 outcomes
get_multiv_elast <- function(model, outcome) {
    
  beta_est1 <- model$beta.ls[[1]]
  beta_est2 <- model$beta.ls[[2]]
  K <- length(beta_est1)
  G <- model$G
  N <- model$N
  TT <- model$TT
  invA <- solve(model$get_A())
  inv_blocks <- as.blockmatrix(as.matrix(invA), nrowe = N*G, ncole = N*G)
  
  ## Direct elasticities
  elasticity_est <- rep(NA, K-1)
  for (k in 2:K) {
    beta_k1 <- beta_est1[k]
    beta_k2 <- beta_est2[k]
    elast_k <- rep(NA, N)
    for (i in 1:N) {
      elast_k_i <- rep(NA, TT)
      for (t in 1:TT) {
        invA_block <- inv_blocks[t,t]
        if (outcome==1) {
          i_idx <- i
          multiplier_i <- invA_block[i_idx,i_idx]
          alt_idx <- i + N
          multiplier_alt <- invA_block[i_idx,alt_idx]
          elast_k_i[t] <- multiplier_i*beta_k1 + multiplier_alt*beta_k2
        } else {
          i_idx <- i + N
          multiplier_i <- invA_block[i_idx,i_idx]
          alt_idx <- i
          multiplier_alt <- invA_block[i_idx,alt_idx]
          elast_k_i[t] <- multiplier_i*beta_k2 + multiplier_alt*beta_k1
        }
      }
      elast_k[i] <- mean(elast_k_i)
    }
    elasticity_est[k-1] <- mean(elast_k)
  }
  
  ## Indirect elasticities
  indirect_elasticity_est <- rep(NA, K-1)
  for (k in 2:K) {
    beta_k1 <- beta_est1[k]
    beta_k2 <- beta_est2[k]
    elast_k <- rep(NA, N)
    for (i in 1:N) {
      elast_k_i <- rep(NA, TT)
      for (t in 1:TT) {
        invA_block <- inv_blocks[t,t]
        if (outcome==1) {
          i_idx <- i
          noti_idx <- (1:N)[-i]
          multiplier_i <- invA_block[noti_idx,i_idx]
          alt_idx <- i + N
          multiplier_alt <- invA_block[noti_idx,alt_idx]
          elast_k_i[t] <- sum(multiplier_i*beta_k1) + sum(multiplier_alt*beta_k2)
        } else {
          i_idx <- i + N
          noti_idx <- ((N+1):(2*N))[-i]
          multiplier_i <- invA_block[noti_idx,i_idx]
          alt_idx <- i
          multiplier_alt <- invA_block[noti_idx,alt_idx]
          elast_k_i[t] <- sum(multiplier_i)*beta_k2 + sum(multiplier_alt)*beta_k1
        }
      }
      elast_k[i] <- mean(elast_k_i)
    }
    indirect_elasticity_est[k-1] <- mean(elast_k)
  }
  out <- list(direct_elast = elasticity_est, indirect_elast = indirect_elasticity_est)
  return(out)
}




###########################
# Fit DIMSTAR to battle events and civilian events, no dependence
###########################

## Battles
X.ls <- list(X)
y.ls <- list(panel.df$isis_battle)

model1_bat <- DISTAR$new(X.ls = X.ls, y.ls = y.ls, W_t = W_t, 
                    N = N, G = 1, TT = TT,
                    count = TRUE, 
                    spatial_dep = FALSE, 
                    temporal_dep = FALSE, 
                    outcome_dep = FALSE)

# Train
set.seed(0)
ch_vec <- model1_bat$train(maxiter = 15, M = 50, abs_tol = 1e-4, soft_init = FALSE)
plot(ch_vec)

# Compute VCOV, print table
model1_bat$compute_vcov(M = 100)
model1_bat$coeftest()
coef1_bat <- model1_bat$coeftest(print_out = FALSE)

# Get elasticities
elast1_bat.ls <- get_univ_elast(model1_bat)


## Civilians
X.ls <- list(X)
y.ls <- list(panel.df$isis_civilian)

model1_civ <- DISTAR$new(X.ls = X.ls, y.ls = y.ls, W_t = W_t, 
                         N = N, G = 1, TT = TT,
                         count = TRUE, 
                         spatial_dep = FALSE, 
                         temporal_dep = FALSE, 
                         outcome_dep = FALSE)

# Train
set.seed(0)
ch_vec <- model1_civ$train(maxiter = 15, M = 50, abs_tol = 1e-4, soft_init = FALSE)
plot(ch_vec)

# Compute VCOV, print table
model1_civ$compute_vcov(M = 100)
model1_civ$coeftest()
coef1_civ <- model1_civ$coeftest(print_out = FALSE)

# Get elasticities
elast1_civ.ls <- get_univ_elast(model1_civ)


###########################
# Fit DIMSTAR to battle events and civilian events, with dependence
###########################

## Battles
X.ls <- list(X)
y.ls <- list(panel.df$isis_battle)

model2_bat <- DISTAR$new(X.ls = X.ls, y.ls = y.ls, W_t = W_t, 
                     N = N, G = 1, TT = TT,
                     count = TRUE, 
                     spatial_dep = TRUE, 
                     temporal_dep = TRUE, 
                     outcome_dep = FALSE)

# Train 
set.seed(0)
ch_vec <- model2_bat$train(maxiter = 15, M = 50, abs_tol = 1e-4, soft_init = FALSE)
plot(ch_vec)

# Compute VCOV, print table
model2_bat$compute_vcov(M = 100)
model2_bat$coeftest()
coef2_bat <- model2_bat$coeftest(print_out = FALSE)

# Get elasticities
elast2_bat.ls <- get_univ_elast(model2_bat)



## Civilians
X.ls <- list(X)
y.ls <- list(panel.df$isis_civilian)

model2_civ <- DISTAR$new(X.ls = X.ls, y.ls = y.ls, W_t = W_t, 
                         N = N, G = 1, TT = TT,
                         count = TRUE, 
                         spatial_dep = TRUE, 
                         temporal_dep = TRUE, 
                         outcome_dep = FALSE)

# Train
set.seed(0)
ch_vec <- model2_civ$train(maxiter = 15, M = 50, abs_tol = 1e-4, soft_init = FALSE)
plot(ch_vec)

# Compute VCOV, print table
model2_civ$compute_vcov(M = 100)
model2_civ$coeftest()
coef2_civ <- model2_civ$coeftest(print_out = FALSE)

# Get elasticities
elast2_civ.ls <- get_univ_elast(model2_civ)



###########################
# Fit DIMSTAR to battle and civilian events, spatio-temporal dependence
###########################

X.ls <- list(X, X)
y.ls <- list(panel.df$isis_battle, panel.df$isis_civilian)

model3 <- DISTAR$new(X.ls = X.ls, y.ls = y.ls, W_t = W_t, 
                     N = N, G = 2, TT = TT,
                     count = TRUE, 
                     spatial_dep = TRUE, 
                     temporal_dep = TRUE, 
                     outcome_dep = TRUE)

# Train
ch_vec <- model3$train(maxiter = 15, M = 50, abs_tol = 1e-4, soft_init = TRUE)
plot(ch_vec)

# Compute VCOV, print table
model3$compute_vcov(M = 100)
model3$coeftest()
coef3 <- model3$coeftest(print_out = FALSE)

# Get elasticities
elast3_out1.ls <- get_multiv_elast(model3, 1)
elast3_out2.ls <- get_multiv_elast(model3, 2)


###########################
# Set up out-of-saple prediction
###########################

set.seed(0)
N_missing <- round(0.33*nrow(panel.df))
missing_bat <- sample(nrow(panel.df), size = N_missing, replace = FALSE)
missing_civ <- sample(nrow(panel.df), size = N_missing, replace = FALSE)

rmse <- function(x, y) {sqrt(mean((x - y)^2))}
mae <- function(x, y) {mean(abs(x-y))}

###########################
# No dependence prediction
###########################

## Battles
X.ls <- list(X)
y.ls <- list(panel.df$isis_battle)
y.ls[[1]][missing_bat] <- NA

pred1_bat <- DISTAR$new(X.ls = X.ls, y.ls = y.ls, W_t = W_t, 
                         N = N, G = 1, TT = TT,
                         count = TRUE, 
                         spatial_dep = FALSE, 
                         temporal_dep = FALSE, 
                         outcome_dep = FALSE)

# Train
set.seed(0)
ch_vec <- pred1_bat$train(maxiter = 15, M = 50, abs_tol = 1e-4, soft_init = FALSE)
plot(ch_vec)

# Predict out of sample
pred1_bat$sample_ystar(M = 100, ystar_init = colMeans(pred1_bat$ystar_sample))
pred1_bat_pred <- colMeans(exp(pred1_bat$ystar_sample))
pred1_bat_pred <- pred1_bat_pred[missing_bat]

# Evaluate predictions
pred1_bat_true <- panel.df$isis_battle[missing_bat]
plot(pred1_bat_pred, pred1_bat_true)
rmse(pred1_bat_pred, pred1_bat_true)
mae(pred1_bat_pred, pred1_bat_true)

## Civ
X.ls <- list(X)
y.ls <- list(panel.df$isis_civilian)
y.ls[[1]][missing_civ] <- NA

pred1_civ <- DISTAR$new(X.ls = X.ls, y.ls = y.ls, W_t = W_t, 
                        N = N, G = 1, TT = TT,
                        count = TRUE, 
                        spatial_dep = FALSE, 
                        temporal_dep = FALSE, 
                        outcome_dep = FALSE)

# Train
set.seed(0)
ch_vec <- pred1_civ$train(maxiter = 15, M = 50, abs_tol = 1e-4, soft_init = FALSE)
plot(ch_vec)

# Predict out of sample
pred1_civ$sample_ystar(M = 100, ystar_init = colMeans(pred1_civ$ystar_sample))
pred1_civ_pred <- colMeans(exp(pred1_civ$ystar_sample))
pred1_civ_pred <- pred1_civ_pred[missing_civ]

# Evaluate predictions
pred1_civ_true <- panel.df$isis_civilian[missing_civ]
plot(pred1_civ_pred, pred1_civ_true)
rmse(pred1_civ_pred, pred1_civ_true)
mae(pred1_civ_pred, pred1_civ_true)

###########################
# Spatio-temporal dependence prediction
###########################

## Battles
X.ls <- list(X)
y.ls <- list(panel.df$isis_battle)
y.ls[[1]][missing_bat] <- NA

pred2_bat <- DISTAR$new(X.ls = X.ls, y.ls = y.ls, W_t = W_t, 
                        N = N, G = 1, TT = TT,
                        count = TRUE, 
                        spatial_dep = TRUE, 
                        temporal_dep = TRUE, 
                        outcome_dep = FALSE)

# Train
set.seed(0)
ch_vec <- pred2_bat$train(maxiter = 15, M = 50, abs_tol = 1e-4, soft_init = FALSE)
plot(ch_vec)

# Predict out of sample
pred2_bat$sample_ystar(M = 100, ystar_init = colMeans(pred2_bat$ystar_sample))
pred2_bat_pred <- colMeans(exp(pred2_bat$ystar_sample))
pred2_bat_pred <- pred2_bat_pred[missing_bat]

# Evaluate predictions
pred2_bat_true <- panel.df$isis_battle[missing_bat]
plot(pred2_bat_pred, pred2_bat_true)
rmse(pred2_bat_pred, pred2_bat_true)
mae(pred2_bat_pred, pred2_bat_true)

## Civ
X.ls <- list(X)
y.ls <- list(panel.df$isis_civilian)
y.ls[[1]][missing_civ] <- NA

pred2_civ <- DISTAR$new(X.ls = X.ls, y.ls = y.ls, W_t = W_t, 
                        N = N, G = 1, TT = TT,
                        count = TRUE, 
                        spatial_dep = TRUE, 
                        temporal_dep = TRUE, 
                        outcome_dep = FALSE)

# Train
set.seed(0)
ch_vec <- pred2_civ$train(maxiter = 15, M = 50, abs_tol = 1e-4, soft_init = FALSE)
plot(ch_vec)

# Predict out of sample
pred2_civ$sample_ystar(M = 100, ystar_init = colMeans(pred1_civ$ystar_sample))
pred2_civ_pred <- colMeans(exp(pred2_civ$ystar_sample))
pred2_civ_pred <- pred2_civ_pred[missing_civ]

# Evaluate predictions
pred2_civ_true <- panel.df$isis_civilian[missing_civ]
plot(pred2_civ_pred, pred2_civ_true)
rmse(pred2_civ_pred, pred2_civ_true)
mae(pred2_civ_pred, pred2_civ_true)

###########################
# Multivariate dependence prediction
###########################

## Battles
X.ls <- list(X, X)
y.ls <- list(panel.df$isis_battle, panel.df$isis_civilian)
y.ls[[1]][missing_bat] <- NA
y.ls[[2]][missing_civ] <- NA

pred3 <- DISTAR$new(X.ls = X.ls, y.ls = y.ls, W_t = W_t, 
                        N = N, G = 2, TT = TT,
                        count = TRUE, 
                        spatial_dep = TRUE, 
                        temporal_dep = TRUE, 
                        outcome_dep = TRUE)

# Train
set.seed(0)
ch_vec <- pred3$train(maxiter = 15, M = 50, abs_tol = 1e-4, soft_init = TRUE)
plot(ch_vec)

# Predict out of sample
pred3$sample_ystar(M = 100, ystar_init = colMeans(pred3$ystar_sample))
pred3_pred <- colMeans(exp(pred3$ystar_sample))

bat_idx <- rep(c(rep(TRUE, N), rep(FALSE, N)), TT)
pred3_pred_bat <- pred3_pred[bat_idx]
pred3_pred_bat <- pred3_pred_bat[missing_bat]

civ_idx <- rep(c(rep(F, N), rep(T, N)), TT)
pred3_pred_civ <- pred3_pred[civ_idx]
pred3_pred_civ <- pred3_pred_civ[missing_civ]

# Evaluate predictions
pred3_bat_true <- panel.df$isis_battle[missing_bat]
plot(pred3_pred_bat, pred3_bat_true)
rmse(pred3_pred_bat, pred3_bat_true)
mae(pred3_pred_bat, pred3_bat_true)

pred3_civ_true <- panel.df$isis_civilian[missing_civ]
plot(pred3_pred_civ, pred3_civ_true)
rmse(pred3_pred_civ, pred3_civ_true)
mae(pred3_pred_civ, pred3_civ_true)






