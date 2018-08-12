
##################################################
# Interpretation Class
##################################################

DSI <- R6Class("DSI",
                  public = list(
                    
                    count = NULL,
                    
                    N = NULL,
                    TT = NULL,
                    G = NULL,
                    NGT = NULL,
                    
                    X = NULL,
                    
                    sigma2_vec = NULL,
                    beta.ls = NULL,
                    
                    A = NULL,
                    Ainv = NULL,
                    V_z = NULL,
                    mu_z = NULL,
                    
                    obs_index_df = NULL,
                    
                    initialize = function(obj) {
                      
                      ## Assign attributes from dimstar object
                      self$count <- obj$count
                      self$N <- obj$N
                      self$TT <- obj$TT
                      self$G <- obj$G
                      self$NGT <- obj$NGT
                      self$X <- obj$X
                      self$sigma2_vec <- obj$sigma2_vec
                      self$beta.ls <- obj$beta.ls
                      self$A <- obj$get_A()
                      
                      ## Create a matrix mapping model indices onto i,t,j indices (i.e. i=unit, t=period, j=outcome)
                      ivec <- rep(rep(1:self$N, self$G), self$TT)
                      tvec <- rep(1:self$TT, each = self$N * self$G)
                      jvec <- rep(rep(1:self$G, each = self$N), self$TT)
                      self$obs_index_df <- data.frame(i = ivec, t = tvec, j = jvec)
                      
                      ## Calculate mean (mu_z) and variance (V_z) of z (= ystar)
                      Xbeta <- self$X%*%unlist(self$beta.ls)
                      self$mu_z <- scipy$sparse$linalg$bicgstab(r_to_py(self$A), Xbeta)[[1]]
                      
                      sigma2_long <- (self$sigma2_vec)[jvec]
                      Sigma <- .sparseDiagonal(self$NGT)*sigma2_long
                      self$Ainv <- solve(self$A)
                      self$V_z <- self$Ainv%*%Sigma%*%t(self$Ainv)
                    }, 
                    
                    expected_value = function(i, t, j) {
                      
                      ## Get model index of requested observation
                      l <- which(self$obs_index_df$i == i & self$obs_index_df$j == j & self$obs_index_df$t == t)
                      
                      ## Get predictor probability or expected value
                      V_z_ll <- self$V_z[l,l]
                      mu_z_l <- self$mu_z[l]
                      if (self$count) {
                        ev <- exp(mu_z_l + V_z_ll/2)
                      } else {
                        ev <- pnorm(mu_z_l/V_z_ll)
                      }
                      
                      return(ev)
                    }, 
                    
                    marginal_effect = function(iy, ty, jy, 
                                               ix, tx, jx, 
                                               k, semi_elast = FALSE) {
                      
                      ## Get model index of target outcome & regressor
                      l <- which(self$obs_index_df$i == iy & self$obs_index_df$j == jy & self$obs_index_df$t == ty)
                      m <- which(self$obs_index_df$i == ix & self$obs_index_df$j == jx & self$obs_index_df$t == tx)
                      
                      ## Marginal effect
                      V_z_ll <- self$V_z[l,l]
                      mu_z_l <- self$mu_z[l]
                      beta_k <- self$beta.ls[[jx]][k]
                      abeta_k <- self$Ainv[l,m]*beta_k
                      if (self$count) {
                        if (semi_elast) {
                          me <- abeta_k
                        } else {
                          me <- exp(mu_z_l + V_z_ll/2)*abeta_k 
                        }
                      } else {
                        if (semi_elast) {
                          warning("Can't compute semi-elasticities for binary outcome models. Returning marginal effects instead.")
                        }
                        me <- dnorm(mu_z_l/V_z_ll)*(abeta_k/V_z_ll)
                      }
                      
                      return(me)
                    }, 
                    
                    conditional_expectation = function(i, j, t, 
                                                       icond, jcond, tcond, ycond, 
                                                       M = 100) {
                      
                      ## Get model index of target outcome & conditioning outcome
                      l <- which(self$obs_index_df$i == i & self$obs_index_df$j == j & self$obs_index_df$t == t)
                      m <- which(self$obs_index_df$i == icond & self$obs_index_df$j == jcond & self$obs_index_df$t == tcond)
                      y_m = ycond
                      
                      ## Get params of joint z distribution
                      mu <- self$mu_z[c(l,m)]
                      Sigma <- self$V_z[c(l,m), c(l,m)]
                      H <- solve(Sigma)
                      
                      ## Gibbs sampling
                      z_l <- 0
                      z_m <- 0
                      z_l_vec <- rep(NA, M)
                      for (iter in 1:M) {
                        
                        # Sample z_m | z_l, y_m
                        # propto p(z_m | z_l) * p(y_m | z_m)
                        mu_m <- mu[2] - (1/H[2,2])*H[2,1]*(z_l - mu[1])
                        Sigma_m <- 1/H[2,2]
                        if (self$count) {
                          sampler_m <- arscpp::plnposterior(y = y_m, mu = mu_m, sd = sqrt(Sigma_m))
                          z_m <- as.numeric(sampler_m$sample(1))
                        } else {
                          upper <- ifelse(y_m == 0, 0, Inf)
                          lower <- ifelse(y_m == 0, -Inf, 0)
                          z_m <- tmvtnorm::rtmvnorm(n = 1, mean = mu_m, sigma = Sigma_m, 
                                                    lower = lower, upper = upper, algorithm = "gibbs")
                        }
                        
                        # Sample z_l | z_m
                        mu_l <- mu[1] - (1/H[1,1])*H[1,2]*(z_m - mu[2])
                        Sigma_l <- 1/H[1,1]
                        z_l <- rnorm(1, mu_l, sqrt(Sigma_l))

                        # Store
                        z_l_vec[iter] <- z_l
                      }
                      
                      ## Sample y_l_vec
                      if (self$count) {
                        y_l_vec <- rpois(M, exp(z_l_vec))
                      } else {
                        y_l_vec <- ifelse(z_l_vec > 0, 1, 0)
                      }

                      ## Return expected value
                      return(mean(y_l_vec))
                    }
))



