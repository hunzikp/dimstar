##################################################
# Test
##################################################

library(dimstar)
reticulate::use_condaenv("py36")

#########################
# Set parameters

N <- 10*10
TT <- 1
G <- 1
n_pairs <- ifelse(G > 1, ncol(combn(G, 2)), 1)

temporal_dep <- FALSE
spatial_dep <- TRUE
outcome_dep <- FALSE
count <- FALSE

rho <- ifelse(spatial_dep, 0.75, 0)
gamma <- ifelse(temporal_dep, 0.2, 0)
lambda <- ifelse(outcome_dep, 0.5, 0)

if (count) {
  beta <- c(1,1)
} else {
  beta <- c(-0.5,1)
}
sigma2 <- 1


#########################
# Simulate data

set.seed(0)
out <- simulate_data(N, TT, G, count, rho_vec = rep(rho, G), lambda_vec = rep(lambda, n_pairs), gamma_vec = rep(gamma, G),
                     beta.ls = rep(list(beta), G), sigma2_vec = rep(sigma2, G))
X.ls <- out$X.ls
y.ls <- out$y.ls
W_t <- out$W_t
ystar_full <- out$ystar


#########################
# Make model

model <- DISTAR$new(X.ls = X.ls, y.ls = y.ls, W_t = W_t, N = N, G = G, TT = TT,
                    count = count, spatial_dep = spatial_dep, temporal_dep = temporal_dep, outcome_dep = outcome_dep)


#########################
# Train model

chng_vec <- model$train(maxiter = 20, M = 50, abs_tol = 1e-5, burnin = 0, thinning = 1, verbose = TRUE)
plot(chng_vec)
model$beta.ls
model$sigma2_vec
model$lambda_vec
model$gamma_vec
model$rho_vec

#########################
# Compute SEs

model$compute_vcov(M = 100, thinning = 4)
model$coeftest()


#########################
# Interpretation

interObj <- DSI$new(model)

## Get expected values for all obs
ev <- rep(NA, model$NGT)
iter <- 1
for (j in 1:G) {
  for (tt in 1:TT) {
    for (i in 1:N) {
      ev[iter] <- interObj$expected_value(i = i, t = tt, j = j)
      iter <- iter + 1
    }
  }
}
plot(ev, unlist(y.ls))
sqrt(mean(((ev - unlist(y.ls))^2)))
abline(a = 0, b = 1)


#########################
# Get marginal effects

interObj$marginal_effect(iy = 1, ty = 1, jy = 1, 
                         ix = 1, tx = 1, jx = 2, k = 2)

#########################
# Simulate conditional expectation

interObj$expected_value(i = 1, t = 1, j = 1)

Ey_1 <- interObj$conditional_expectation(i = 1, j = 1, t = 1, 
                                         icond = 2, jcond = 1, tcond = 1, 
                                         ycond = 0, M = 5000)
Ey_1

Ey_2 <- interObj$conditional_expectation(i = 1, j = 1, t = 1, 
                                         icond = 2, jcond = 1, tcond = 1, 
                                         ycond = 1, M = 5000)
Ey_2

