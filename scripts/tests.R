##################################################
# Test
##################################################

library(dimstar)
reticulate::use_condaenv("py36")

#########################
# Set parameters

N <- 256
TT <- 20
G <- 1
n_pairs <- ifelse(G > 1, ncol(combn(G, 2)), 1)

temporal_dep <- TRUE
spatial_dep <- TRUE
outcome_dep <- FALSE
count <- FALSE

rho <- ifelse(spatial_dep, 0.45, 0)
gamma <- ifelse(temporal_dep, 0.45, 0)
lambda <- ifelse(outcome_dep, 0.25, 0)

if (count) {
  beta <- c(3,1)
} else {
  beta <- c(-0.5, 1)
}
sigma2 <- 1


#########################
# Simulate data

set.seed(6)
out <- simulate_data(N, TT, G, count, rho_vec = rep(rho, G), lambda_vec = rep(lambda, n_pairs), gamma_vec = rep(gamma, G),
                     beta.ls = rep(list(beta), G), sigma2_vec = rep(sigma2, G))
X.ls <- out$X.ls
y.ls <- out$y.ls
W_t <- out$W_t
ystar_full <- out$ystar

# Randomly censor outcome
# y.ls[[1]][sample(length(y.ls[[1]]), size = 10, replace = FALSE)] <- NA

#########################
# Make model

model <- DISTAR$new(X.ls = X.ls, y.ls = y.ls, W_t = W_t, N = N, G = G, TT = TT,
                    count = count, spatial_dep = spatial_dep, temporal_dep = temporal_dep, outcome_dep = outcome_dep)

#########################
# Benchmarking
# library(rbenchmark)
# 
# set.seed(0)
# model$sample_ystar(M = 50, ystar_init = colMeans(model$ystar_sample))
# plot(ystar_full, colMeans(model$ystar_sample))
# 
# theta <- model$pack_theta()
# rbenchmark::benchmark(test = {model$E_llik(theta)})
# test replications elapsed relative user.self sys.self user.child sys.child
# 1 test          100   0.102        1       0.1        0          0         0


#########################
# Train test, iterative

model$sample_ystar(M = 50, ystar_init = colMeans(model$ystar_sample))
system.time(model$update_theta())
theta <- model$pack_theta()
system.time(model$E_llik(theta))
# model$beta.ls
# model$sigma2_vec
# model$lambda_vec
# model$gamma_vec
# model$rho_vec

#########################
# Train test

chng_vec <- model$train(maxiter = 50, M = 50, abs_tol = 1e-5, burnin = 0, thinning = 1, soft_init = FALSE, verbose = TRUE)
plot(chng_vec)
model$beta.ls
model$sigma2_vec
model$lambda_vec
model$gamma_vec
model$rho_vec

#########################
# VCOV test
model$compute_vcov(M = 1000, thinning = 10, verbose = TRUE)
model$coeftest()

estim <- model$coeftest(print_out = F)
theta_est <- estim$coef
se <- estim$se
se
# 0.05516480 0.03882166 0.02401263 0.02210312

theta_lo <- theta_est - qnorm(0.95)*se
theta_hi <- theta_est +  qnorm(0.95)*se


### TODO
# - [DONE] Move code to own project/repo; package; write basic tests.
# - [DONE] Check out whether we can use constrained optim on lambda/gamma/rho instead of tanh/atanh transform
# - [NOT POSSIBLE ATM] Speed up likelihood eval, esp for VC estimation
# - [DONE] Ensure that the model works with G = T = 1 
# - [NOT EFFICIENT] Ensure that we use sparse cholesky for det eval when G > 1
# - [DONE] Ensure that we use efficient log-det eval if G = 1