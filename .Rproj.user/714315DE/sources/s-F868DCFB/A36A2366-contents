##################################################
# Utility Functions
##################################################

## Softplus and inverse softplus
softplus <- function(x) {
  log(exp(x) + 1)
}

asoftplus <- function(x) {
  log(exp(x)-1)
}

## Make full NGTxNGT panel spatial weights matrices; one for each outcome
make_spatial_weights <- function(W_t, G, TT) {
  
  matrix_list <- vector('list', G)
  N <- nrow(W_t)
  sparseZeroMatrix <- sparseMatrix(i = {}, j = {}, dims = c(N, N))
  
  for (k in 1:G) {
    bmat_list <- rep(list(sparseZeroMatrix), G)
    bmat_list[[k]] <- W_t
    W_period <- bdiag(bmat_list)
    matrix_list[[k]] <- bdiag(rep(list(W_period), TT))
    matrix_list[[k]] <- as(matrix_list[[k]], 'dgCMatrix')
  }
  
  return(matrix_list)
}

## Make full NGTxNGT panel temporal weights matrices; one for each outcome
make_temporal_weights <- function(N, G, TT) {
  
  matrix_list <- vector('list', G)
  
  NG <- N*G
  for (k in 1:G) {
    if (TT > 1) {
      rows <- c()
      cols <- c()
      for (tt in 2:TT) {
        start_row <- (tt-1)*NG + (k-1)*N + 1
        end_row <- start_row +  N - 1
        start_col <- start_row - NG
        end_col <- end_row - NG
        rows <- c(rows, start_row:end_row)
        cols <- c(cols, start_col:end_col)
      }
      TM <- sparseMatrix(i = rows, j = cols, x = TRUE, dims = c(NG*TT, NG*TT))
    } else {
      TM <- sparseMatrix(i = {}, j = {}, dims = c(NG*TT, NG*TT))
    }
    matrix_list[[k]] <- TM
    matrix_list[[k]] <- as(matrix_list[[k]], 'dgCMatrix')
  }
  
  return(matrix_list)
}

## Make full NGTxNGT panel outcome weights matrices; one for each outcome-pair
make_outcome_weights <- function(N, G, TT) {
  
  if (G == 1) {
    matrix_list <- list(sparseMatrix(i = {}, j = {}, dims = c(N*TT, N*TT)))
    matrix_list[[1]] <- as(matrix_list[[1]], 'dgCMatrix')
    return(matrix_list)
  }
  
  pair_mat <- combn(G, 2)
  n_pairs <- ncol(pair_mat)
  matrix_list <- vector('list', n_pairs)
  
  NG <- N*G
  for (k in 1:n_pairs) {
    
    ## Make the inner matrix
    outcome_1 <- pair_mat[1,k]
    outcome_2 <- pair_mat[2,k]
    
    row_start <- (outcome_1-1)*N + 1
    row_end <- row_start + N - 1
    col_start <- (outcome_2-1)*N + 1
    col_end <- col_start + N - 1
    idx1 <- c(row_start:row_end, col_start:col_end)
    idx2 <- c(col_start:col_end, row_start:row_end)
    mat <- sparseMatrix(i = idx1, j = idx2, x = TRUE, dims = c(NG, NG))
    
    matrix_list[[k]] <- bdiag(rep(list(mat), TT))
    matrix_list[[k]] <- as(matrix_list[[k]], 'dgCMatrix')
  }
  
  return(matrix_list)
}

## Make covariate matrix
make_X <- function(X_list, N, TT) {
  ## Creates a covariate matrix for multiple outcomes
  # X_list: list containing one covariate matrix per outcome
  # N: Number of units
  # TT: Number of periods
  
  G <- length(X_list)  # Number of outcomes
  K_vec <- unlist(lapply(X_list, ncol))  # Number of regressors per outcome
  K <- sum(K_vec)  # Total number of regressors
  X <- matrix(0, N*TT*G, K)
  for (t in 1:TT) {
    # Determine which rows of X_list to pick
    t_start <- (t-1)*N + 1
    t_end <- t_start + N - 1
    for (k in 1:G) {
      # Determine which columns of X to fill in
      if (k == 1) {
        col_start <- 1
        col_end <- K_vec[1]
      } else {
        col_start <- sum(K_vec[1:(k-1)]) + 1
        col_end <- col_start + K_vec[k] - 1
      }
      # Determine which rows of X to fill in
      row_start <- (t-1)*N*G + (k-1)*N + 1
      row_end <- row_start + N - 1
      # Fill in
      X[row_start:row_end, col_start:col_end] <- X_list[[k]][t_start:t_end,]
    }
  }
  
  return(X)
}

## Dummy class so we can use lmtest::coeftest
make_dummyclass <- function(x) {
  x <- list(coef = x)
  class(x) <- 'dummyclass'
  return(x)
}
coef.dummyclass <- function(x) {
  x$coef
}