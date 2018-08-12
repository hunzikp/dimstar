

##################################################
# Estimation Class
##################################################

DISTAR <- R6Class("DISTAR",
                  public = list(

                    count = NULL,
                    spatial_dep = NULL,
                    temporal_dep = NULL,
                    outcome_dep = NULL,

                    N = NULL,
                    TT = NULL,
                    G = NULL,
                    NGT = NULL,

                    X = NULL,
                    X.ls = NULL,
                    y = NULL,
                    W_t = NULL,

                    W.ls = NULL,
                    TM.ls = NULL,
                    O.ls = NULL,
                    n_pairs = NULL,
                    lambda_map.ls = NULL,

                    rho_vec = NULL,
                    gamma_vec = NULL,
                    lambda_vec = NULL,

                    sigma2_vec = NULL,
                    beta.ls = NULL,

                    ystar_sample = NULL,
                    kern_obj = NULL,
                    ldA_fun = NULL,

                    outcome_index_vec = NULL,

                    VC_mat = NULL,

                    initialize = function(X.ls, y.ls, W_t, N, G, TT, count,
                                          spatial_dep = FALSE, temporal_dep = FALSE, outcome_dep = FALSE) {

                      ## Check inputs
                      if (TT == 1 & temporal_dep) {
                        warning("You requested temporal dependence estimation, but the data only cover a single period. 
                                Temporal dependence will be set to FALSE.")
                        temporal_dep <- FALSE
                      }
                      if (G == 1 & outcome_dep) {
                        warning("You requested outcome dependence estimation, but the data only cover a single outcome. 
                                Outcome dependence will be set to FALSE.")
                        outcome_dep <- FALSE
                      }
                      
                      ## Assign data
                      self$X <- make_X(X_list = X.ls, N = N, TT = TT)
                      self$X.ls <- X.ls
                      self$W_t <- W_t

                      ## Make y vector
                      y <- c()
                      for (t in 1:TT) {
                        start <- (t-1)*N + 1
                        end <- start + N - 1
                        for (k in 1:G) {
                          y <- c(y, y.ls[[k]][start:end])
                        }
                      }
                      self$y <- y

                      ## Assign outcome switch
                      self$count <- count

                      ## Assign dependence switches
                      self$spatial_dep <- spatial_dep
                      self$temporal_dep <- temporal_dep
                      self$outcome_dep <- outcome_dep

                      ## Assign data dimensions
                      self$N <- N
                      self$TT <- TT
                      self$G <- G
                      self$NGT <- N*G*TT

                      ## Make W.ls: One NGTxNGT matrix per outcome
                      self$W.ls <- make_spatial_weights(W_t, G, TT)

                      ## Make TM.ls: One NGTxNGT matrix per outcome
                      # NOTE: If TT = 1, returns G NGxNG matrices of zeros
                      self$TM.ls <- make_temporal_weights(N, G, TT)

                      ## Make O.ls: One NGTxNGT matrix per outcome-pair
                      # NOTE: If G = 1, returns a single NTxNT matrix of zeros
                      self$O.ls <- make_outcome_weights(N, G, TT)
                      self$n_pairs <- length(self$O.ls) # set to 1 if G = 1

                      ## Map lambdas/Os onto outcomes
                      if (self$G > 1) {
                        pair_mat <- combn(self$G, 2)
                        self$lambda_map.ls <- vector('list', self$G)
                        for (i in 1:self$G) {
                          self$lambda_map.ls[[i]] <- c(which(pair_mat[1,]==i), which(pair_mat[2,]==i))
                        }
                      } else {
                        self$lambda_map.ls <- rep(list(1), self$G)
                      }

                      ## Initialize interdependence parameters
                      self$rho_vec <- rep(0, self$G)
                      self$gamma_vec <- rep(0, self$G)
                      self$lambda_vec <- rep(0, self$n_pairs)

                      ## Initialize variance
                      self$sigma2_vec <- rep(1, self$G)

                      ## Initialize beta
                      beta.ls <- vector('list', self$G)
                      for (k in 1:self$G) {
                        beta.ls[[k]] <- rep(0, ncol(X.ls[[k]]))
                      }
                      self$beta.ls <- beta.ls

                      ## Initialize ystar
                      self$ystar_sample <- matrix(0, 1, N*TT*G)

                      ## Make the kernel compute object
                      self$kern_obj <- Kernel$new(self$W.ls, self$TM.ls, self$O.ls, self$X, self$N*self$G)

                      ## Make the outcome index vector
                      # Indicating outcome per y position
                      self$outcome_index_vec <- rep(rep(1:self$G, each = self$N), self$TT)
                      
                      ## Setup of log|A| function
                      # If G = 1, we can use lookup to eval log|A|
                      # Here we prepare the lookup function
                      if (G == 1) {
                        print(paste('Setting up log-determinants...'))
                        flush.console()
                        rho_seq <- seq(-0.99, 0.99, length = 20)
                        ldA_seq <- rep(0, length(rho_seq))
                        for (j in 1:length(rho_seq)) {
                          ldA_seq[j] <- self$TT * Matrix::determinant(.sparseDiagonal(self$N) - rho_seq[j]*self$W_t)$modulus
                        }
                        self$ldA_fun <- splinefun(x = rho_seq, y = ldA_seq)
                        print(paste('Done.'))
                        flush.console()
                      }

                    },

                    get_A = function() {
                      # Returns the NGTxNGT inverse dependence multiplier with class params

                      A <- .sparseDiagonal(self$NGT)
                      for (k in 1:self$G) {
                        A <- A - self$rho_vec[k]*self$W.ls[[k]] - self$gamma_vec[k]*self$TM.ls[[k]]
                      }
                      for (k in 1:self$n_pairs) {
                        A <- A - self$O.ls[[k]]*self$lambda_vec[k]
                      }

                      return(A)
                    },

                    sample_ystar = function(M, ystar_init) {
                      if (self$count) {
                        self$sample_ystar_count(M, ystar_init)
                      } else {
                        self$sample_ystar_binary(M, ystar_init)
                      }

                      self$kern_obj$update_ystar_sample(self$ystar_sample)
                    },

                    sample_ystar_count = function(M, ystar_init) {

                      # Get parameters
                      A <- self$get_A()
                      Xbeta <- self$X%*%unlist(self$beta.ls)
                      mu <- scipy$sparse$linalg$bicgstab(r_to_py(A), Xbeta)[[1]]
                      inv_sigma2_long <- (1/self$sigma2_vec)[self$outcome_index_vec]
                      invSigma <- .sparseDiagonal(self$NGT)*inv_sigma2_long
                      H <- t(A)%*%invSigma%*%A

                      # Sample
                      sampler <- arscpp::mplnposterior(self$y, as.vector(mu), H)
                      smp <- sampler$sample(M, ystar_init)

                      # Assign
                      self$ystar_sample <- smp
                    },

                    sample_ystar_binary = function(M, ystar_init) {

                      # Get parameters
                      A <- self$get_A()
                      Xbeta <- self$X%*%unlist(self$beta.ls)
                      mu <- scipy$sparse$linalg$bicgstab(r_to_py(A), Xbeta)[[1]]
                      inv_sigma2_long <- (1/self$sigma2_vec)[self$outcome_index_vec]
                      invSigma <- .sparseDiagonal(self$NGT)*inv_sigma2_long
                      H <- t(A)%*%invSigma%*%A

                      # Set limits
                      lower <- rep(0, length(self$y))
                      lower[is.na(self$y) | self$y == 0] <- -Inf
                      upper <- rep(0, length(self$y))
                      upper[is.na(self$y) | self$y == 1] <- Inf

                      # Sample
                      smp <- rtmvnorm.sparseMatrix(M, mean = as.vector(mu), H = H,
                                                   lower = lower, upper = upper,
                                                   start.value = ystar_init)

                      # Assign
                      self$ystar_sample <- smp
                    },

                    pack_theta = function() {
                      ##  Get vector containing all parameters to be estimated
                      
                      theta <- c(unlist(self$beta.ls))
                      if (self$outcome_dep) {
                        theta <- c(theta, atanh(self$lambda_vec))
                      }
                      if (self$count) {
                        theta <- c(theta, asoftplus(self$sigma2_vec))
                      }
                      if (self$spatial_dep) {
                        theta <- c(theta, atanh(self$rho_vec))
                      }
                      if (self$temporal_dep) {
                        theta <- c(theta, atanh(self$gamma_vec))
                      }

                      return(theta)
                    },

                    unpack_theta = function(theta) {
                      ## Takes the packed theta vector and returns a list with all parameters as attribues

                      K <- length(unlist(self$beta.ls))
                      beta <- theta[1:K]
                      beta.ls <- vector('list', self$G)
                      idx <- 1
                      for (k in 1:self$G) {
                        K.k <- length(self$beta.ls[[k]])
                        beta.ls[[k]] <- beta[idx:(idx+K.k-1)]
                        idx <- idx + K.k
                      }

                      idx <- K+1
                      if (self$outcome_dep) {
                        idx_end <- idx + self$n_pairs - 1
                        lambda_vec <- tanh(theta[idx:idx_end])
                        idx <- idx_end + 1
                      } else {
                        lambda_vec <- rep(0, self$n_pairs)
                      }
                      if (self$count) {
                        idx_end <- idx + self$G - 1
                        sigma2_vec <- softplus(theta[idx:idx_end])
                        idx <- idx_end + 1
                      } else {
                        sigma2_vec <- rep(1, self$G)
                      }
                      if (self$spatial_dep) {
                        idx_end <- idx + self$G - 1
                        rho_vec <- tanh(theta[idx:idx_end])
                        idx <- idx_end + 1
                      } else {
                        rho_vec <- rep(0, self$G)
                      }
                      if (self$temporal_dep) {
                        idx_end <- idx + self$G - 1
                        gamma_vec <- tanh(theta[idx:idx_end])
                        idx <- idx_end + 1
                      } else {
                        gamma_vec <- rep(0, self$G)
                      }

                      out <- list(beta.ls = beta.ls, sigma2_vec = sigma2_vec,
                                  rho_vec = rho_vec, gamma_vec = gamma_vec,
                                  lambda_vec = lambda_vec)
                      return(out)
                    },
                    
                    transform_theta = function(theta) {
                      ## Takes the packed theta vector and returns a theta vector where all 
                      ## transformed params are returned to their original support
                      # i.e. from atanh(rho) to rho
                      
                      transformed_theta <- rep(NA, length(theta))
                      
                      K <- length(unlist(self$beta.ls))
                      transformed_theta[1:K] <- theta[1:K]
                      
                      idx <- K+1
                      if (self$outcome_dep) {
                        idx_end <- idx + self$n_pairs - 1
                        transformed_theta[idx:idx_end] <- tanh(theta[idx:idx_end])
                        idx <- idx_end + 1
                      }
                      if (self$count) {
                        idx_end <- idx + self$G - 1
                        transformed_theta[idx:idx_end] <- softplus(theta[idx:idx_end])
                        idx <- idx_end + 1
                      }
                      if (self$spatial_dep) {
                        idx_end <- idx + self$G - 1
                        transformed_theta[idx:idx_end] <- tanh(theta[idx:idx_end])
                        idx <- idx_end + 1
                      }
                      if (self$temporal_dep) {
                        idx_end <- idx + self$G - 1
                        transformed_theta[idx:idx_end] <- tanh(theta[idx:idx_end])
                      }

                      return(transformed_theta)
                    },
                    
                    get_ldA = function(rho_vec, lambda_vec) {
                      ## Get the log determinant of A
                      
                      if (self$G == 1) {
                        ldA <- self$ldA_fun(rho_vec[1])
                      } else {
                        A_t <- self$kern_obj$make_A_t(rho_vec, lambda_vec)
                        # lud <- lu(A_t)
                        # ldA <- self$TT*(sum(log(diag(lud@L))) + sum(log(diag(lud@U))))
                        ldA <- self$TT*determinant(A_t)$modulus
                      }
                      
                      return(ldA)
                    },
                    
                    E_llik = function(theta) {

                      par <- self$unpack_theta(theta)
                      beta.ls <- par$beta.ls
                      sigma2_vec <- par$sigma2_vec
                      rho_vec <- par$rho_vec
                      gamma_vec <- par$gamma_vec
                      lambda_vec <- par$lambda_vec

                      ## Safety check
                      if (any(!is.finite(1/sigma2_vec))) {
                        return(NaN)
                      }

                      ## Compute expected kernel
                      inv_sigma2_long <- (1/sigma2_vec)[self$outcome_index_vec]
                      kern <- self$kern_obj$get_E_kernel(rho_vec, gamma_vec, lambda_vec,
                                                         inv_sigma2_long, unlist(beta.ls))
                      
                      ## Compute log determinant
                      ldA <- self$get_ldA(rho_vec, lambda_vec)

                      ## Expected likelihood
                      sigma2_long <- (sigma2_vec)[self$outcome_index_vec]
                      ll <- ldA - sum(log(sqrt(sigma2_long))) - 0.5*kern

                      return(-1*ll)
                    },

                    check_constraint = function(theta) {
                      par <- self$unpack_theta(theta)

                      constraint_vec <- rep(NA, self$G)
                      for (k in 1:self$G) {
                        lambda_k <- sum(par$lambda_vec[self$lambda_map.ls[[k]]])
                        constraint_vec[k] <- par$rho_vec[k] + lambda_k # + par$gamma_vec[k] 
                      }

                      return(constraint_vec)
                    },

                    update_theta = function() {

                      theta_init <- self$pack_theta()
                      
                      ineqLB <- rep(-1, self$G)
                      ineqUB <- rep(1, self$G)

                      fit <- solnp(pars = theta_init, fun = self$E_llik,
                                   ineqfun = self$check_constraint,
                                   ineqLB = ineqLB, ineqUB = ineqUB,
                                   control = list(trace = 0))
                      par <- self$unpack_theta(fit$pars)

                      self$beta.ls <- par$beta.ls
                      self$sigma2_vec <- par$sigma2_vec
                      self$rho_vec <- par$rho_vec
                      self$gamma_vec <- par$gamma_vec
                      self$lambda_vec <- par$lambda_vec
                    },

                    train = function(maxiter = 20, M = 50, abs_tol = 1e-3, soft_init = TRUE, burnin = 50, thinning = 1, verbose = FALSE) {
                      # Learn the parameters

                      iter <- 0
                      change_vec <- c()
                      max_change <- Inf
                      while (iter < maxiter & max_change > abs_tol) {
                        
                        ## Record params from last iter
                        theta_old <- self$pack_theta()
                        
                        ## In first itertion: initialize using only spatial dep
                        # Otherwise coordinate ascent tends to get stuck in local minimum
                        if (iter == 0 & soft_init) {
                          
                          # Sample
                          Ms <- burnin + M*thinning
                          self$sample_ystar(M = Ms, ystar_init = colMeans(self$ystar_sample))
                          Mseq <- seq(burnin+1, Ms, by = thinning)
                          self$ystar_sample <- self$ystar_sample[Mseq,]
                          
                          # Shut off temporal and outcome dep
                          od <- self$outcome_dep
                          self$outcome_dep <- FALSE
                          td <- self$temporal_dep
                          self$temporal_dep <- FALSE
                          
                          # Update params
                          self$update_theta()
                          
                          # Set temporal and outcome flags back to original state
                          self$outcome_dep <- od
                          self$temporal_dep <- td
                        }

                        ## Train regularly and record parameter changes
                        Ms <- burnin + M*thinning
                        self$sample_ystar(M = Ms, ystar_init = colMeans(self$ystar_sample))
                        Mseq <- seq(burnin+1, Ms, by = thinning)
                        self$ystar_sample <- self$ystar_sample[Mseq,]
                        self$update_theta()

                        max_change <- max(abs(self$pack_theta() - theta_old))
                        change_vec <- c(change_vec, max_change)
                        iter <- iter + 1
                        
                        if (verbose) {
                          print(paste0("Iteration ", iter, "; Max. parameter change = ", round(max_change, 3)))
                          flush.console()
                        }
                      }
                      return(change_vec)
                    },

                    llik_i = function(theta, i) {

                      par <- self$unpack_theta(theta)
                      beta.ls <- par$beta.ls
                      sigma2_vec <- par$sigma2_vec
                      rho_vec <- par$rho_vec
                      gamma_vec <- par$gamma_vec
                      lambda_vec <- par$lambda_vec

                      ## Safety check
                      if (any(!is.finite(1/sigma2_vec))) {
                        return(NaN)
                      }
                      
                      ## Compute expected kernel
                      inv_sigma2_long <- (1/sigma2_vec)[self$outcome_index_vec]
                      kern <- self$kern_obj$get_kernel_i(rho_vec, gamma_vec, lambda_vec,
                                                         inv_sigma2_long, unlist(beta.ls), i-1)

                      ## Compute log determinant
                      ldA <- self$get_ldA(rho_vec, lambda_vec)
                      
                      ## Expected likelihood
                      sigma2_long <- (sigma2_vec)[self$outcome_index_vec]
                      ll <- ldA - sum(log(sqrt(sigma2_long))) - 0.5*kern

                      return(ll)
                    },

                    get_gradient = function(theta, i) {
                      # Get gradient at a single sample

                      gr <- grad(func = self$llik_i, x = theta, method = 'simple',
                                 i = i)
                      return(gr)
                    },

                    compute_vcov = function(M, thinning = 2, verbose = FALSE) {
                      ## See Louis (1982): "Finding the observed information matrix when using the EM algorithm"

                      self$sample_ystar(M*thinning, ystar_init = colMeans(self$ystar_sample))
                      Mseq <- seq(1, M*thinning, by = thinning)
                      self$ystar_sample <- self$ystar_sample[Mseq,]
                      
                      theta <- self$pack_theta()
                      
                      grad.mat <- matrix(NA, M, length(theta))
                      hess.ls <- vector('list', M)
                      for (i in 1:M) {

                        gr <- self$get_gradient(theta, i)
                        hs <- nl.jacobian(theta, self$get_gradient, i = i)

                        grad.mat[i,] <- gr
                        hess.ls[[i]] <- hs
                        
                        if (verbose) {
                          print(i)
                          flush.console()
                        }
                      }

                      # Get observed information matrix
                      IM <- -1*Reduce('+', hess.ls)/M - var(grad.mat)

                      # Get variance estimate for parameters
                      VC <- solve(Matrix(IM))  # VC is inverse of observed Fisher Information
                      VC_pd <- as.matrix(nearPD(VC)$mat) # Force positive definite

                      self$VC_mat <- VC_pd
                    },

                    get_theta_names = function() {

                      theta_names <- c()

                      for (k in 1:self$G) {
                        if (is.null(colnames(self$X.ls[[k]]))) {
                          beta_names <- paste0(paste0('beta_', k, '_', 0:(ncol(self$X.ls[[k]])-1)))
                        } else {
                          beta_names <- colnames(self$X.ls[[k]])
                        }
                        theta_names <- c(theta_names, beta_names)
                      }
                      if (self$outcome_dep) {
                        theta_names <- c(theta_names, paste0("lambda", 1:self$n_pairs))
                      }
                      if (self$count) {
                        theta_names <- c(theta_names, paste0("sigma2_", 1:self$G))
                      }
                      if (self$spatial_dep) {
                        theta_names <- c(theta_names, paste0("rho_", 1:self$G))
                      }
                      if (self$temporal_dep) {
                        theta_names <- c(theta_names, paste0("gamma_", 1:self$G))
                      }

                      return(theta_names)
                    },

                    coeftest = function(S = 1000, print_out = TRUE) {
                      ## Does two things:
                      # 1. Uses MC simulation to get VC matrix for transformed params
                      # 2. Prints pretty regression table (or returns list with coefs & SEs)
                      
                      ## MC sim
                      theta <- self$pack_theta()
                      theta_trans <- self$transform_theta(theta)
                      theta_sim <- mvtnorm::rmvnorm(n = S, mean = theta, sigma = self$VC_mat)
                      theta_trans_sim <- t(apply(theta_sim, MARGIN = 1, FUN = self$transform_theta))
                      theta_trans_vc <- var(theta_trans_sim)

                      ## Pack in dummy class and print coef test
                      names(theta_trans) <- self$get_theta_names()
                      dc <- make_dummyclass(theta_trans)
                      
                      if (print_out) {
                        lmtest::coeftest(x = dc, vcov. = theta_trans_vc)
                      } else {
                        out <- list(coef = theta_trans, se = sqrt(diag(theta_trans_vc)))
                        return(out)
                      }
                    }
                  ))
