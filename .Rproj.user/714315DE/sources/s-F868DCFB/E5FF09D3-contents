##################################################
# Monte Carlo simulation of multivariate discrete spatio-temporally autoregressive data
##################################################

## Make spatial weights matrix from lattice
# N must have natural root
make_W_t <- function(N) {
  
  ras <- raster(matrix(1, sqrt(N), sqrt(N)))
  spdf <- rasterToPolygons(ras)
  B <- gIntersects(spdf, byid = TRUE)
  diag(B) <- FALSE
  B <- B*1
  W_t <- B / rowSums(B)
  W_t <- Matrix(W_t, sparse = TRUE)
  
  return(W_t)
}

## Monte Carlo simulation of DIM-STAR data
# N must have natural root
simulate_data <- function(N, TT, G, count,
                          rho_vec, lambda_vec, gamma_vec,
                          beta.ls, sigma2_vec) {
  
  if (!count) {
    sigma2_vec <- rep(1, G)
  }
  
  ## Make NxN spatial weights matrix
  W_t <- make_W_t(N)
  
  ## Make NGT x NGT spatial weights matrix list (one per outcome)
  W.ls <- make_spatial_weights(W_t, G, TT)
  
  ## Make NGT x NGT temporal weights matrix list (one per outcome)
  TM.ls <- make_temporal_weights(N, G, TT)
  
  ## Make NGT x NGT outcome weights matrix list (one per outcome pair)
  O.ls <- make_outcome_weights(N, G, TT)
  
  ## Compute A
  A <- .sparseDiagonal(N*TT*G)
  for (k in 1:G) {
    A <- A - rho_vec[k]*W.ls[[k]] - gamma_vec[k]*TM.ls[[k]]
  }
  for (k in 1:length(O.ls)) {
    A <- A - O.ls[[k]]*lambda_vec[k]
  }
  
  ## Make predictors
  X.ls <- vector('list', G)
  for (k in 1:G) {
    K <- length(beta.ls[[k]])
    X <- cbind(1, matrix(rnorm((K-1)*N*TT), N*TT, K-1))
    X.ls[[k]] <- X
  }
  X <- make_X(X_list = X.ls, N = N, TT = TT)
  Xbeta <- X%*%unlist(beta.ls)
  
  ## Make sigma vector
  outcome_index_vec <- rep(rep(1:G, each = N), TT)
  sigma2_long <- sigma2_vec[outcome_index_vec]
  Sigma <- .sparseDiagonal(N*G*TT)*sigma2_long
  sd <- sqrt(sigma2_long)
  
  ## Sample ystar
  eps <- rnorm(N*TT*G, 0, sd)
  ystar <- as.vector(scipy$sparse$linalg$bicgstab(r_to_py(A), Xbeta + eps)[[1]])
  
  ## Sample y
  if (count) {
    y <- rpois(N*TT*G, exp(ystar))
  } else {
    y <- (ystar>0)*1
  }
  
  ## Pack y into list
  y.ls <- vector('list', G)
  for (k in 1:G) {
    inner <- rep(FALSE, N*G)
    st <- ((k-1)*N + 1)
    en <- st + N - 1
    inner[st:en] <- TRUE
    idx <- (1:(N*G*TT))[inner]
    y.ls[[k]] <- y[idx]
  }
  
  ## Prepare return list
  out <- list(y.ls = y.ls, X.ls = X.ls, W_t = W_t, ystar = ystar, A = A, N = N, G = G, TT = TT)
  
  return(out)
}
