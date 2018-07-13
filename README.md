dimstar: Discrete Multivariate Spatio-temporal Autoregressive Models
====================================================================

`dimstar` is an R package implementing a novel MCEM algorithm for
estimating latent-Guassian spatio-temporal autoregressive models for
potentially multivariate binary and count outcomes.

Installation
------------

Install dimstar straight from github:

    require(devtools)
    devtools::install_github('hunzikp/dimstar')

NOTE: This is a development version; documentation is pending. The
package assumes that Python (2.7 / 3.6) and
[scipy](https://www.scipy.org/) is installed on your machine.

Usage
-----

In the following, we replicate the county-level Spatial Probit analysis
of the 1996 presidential elections reported in Section 5.1 of
[Liesenfeld et al.
(2016)](https://ideas.repec.org/h/eme/aecozz/s0731-905320160000037009.html).
The outcome is a dummy indicating whether a county was won by Democratic
candidate Bill Clinton.

    ## Dependencies
    library(spdep)

    ## Loading required package: sp

    ## Loading required package: Matrix

    library(Matrix)
    library(dimstar)

    ## Prepare W
    coords <- as.matrix(elections.df[,c("lon" ,"lat")])
    el_knn <- knearneigh(coords, k = 6, longlat = TRUE)
    el_nb <- knn2nb(knn = el_knn, sym = TRUE)
    W <- nb2mat(el_nb, style = "W")
    W <- Matrix(W, sparse = TRUE)

    ## Prepare X
    X <- model.matrix(clinton ~ log_urban + prop_smcollege + prop_associate + prop_college + prop_gradprof, data = elections.df)

    ## Instantiate model object 
    X.ls <- list(X)
    y.ls <- list(elections.df$clinton)

    model <- DISTAR$new(X.ls = X.ls, y.ls = y.ls, W_t = W, 
                             N = nrow(X), G = 1, TT = 1,
                             count = FALSE, 
                             spatial_dep = TRUE, 
                             temporal_dep = FALSE, 
                             outcome_dep = FALSE)

    ## [1] "Setting up log-determinants..."
    ## [1] "Done."

    ## Train
    set.seed(0)
    ch_vec <- model$train(maxiter = 15, M = 50, abs_tol = 1e-4, verbose = TRUE)

    ## [1] "Iteration 1; Max. parameter change = 6.172"
    ## [1] "Iteration 2; Max. parameter change = 0.773"
    ## [1] "Iteration 3; Max. parameter change = 0.738"
    ## [1] "Iteration 4; Max. parameter change = 0.233"
    ## [1] "Iteration 5; Max. parameter change = 0.147"
    ## [1] "Iteration 6; Max. parameter change = 0.611"
    ## [1] "Iteration 7; Max. parameter change = 0.628"
    ## [1] "Iteration 8; Max. parameter change = 0.199"
    ## [1] "Iteration 9; Max. parameter change = 0.966"
    ## [1] "Iteration 10; Max. parameter change = 0.289"
    ## [1] "Iteration 11; Max. parameter change = 0.083"
    ## [1] "Iteration 12; Max. parameter change = 0.339"
    ## [1] "Iteration 13; Max. parameter change = 0.137"
    ## [1] "Iteration 14; Max. parameter change = 0.641"
    ## [1] "Iteration 15; Max. parameter change = 0.045"

    ## Standard errors
    model$compute_vcov(M = 100)

    ## Report
    model$coeftest()

    ## 
    ## z test of coefficients:
    ## 
    ##                 Estimate Std. Error z value  Pr(>|z|)    
    ## (Intercept)     0.605057   0.120679  5.0138 5.338e-07 ***
    ## log_urban       4.275470   5.346903  0.7996  0.423933    
    ## prop_smcollege -2.838460   0.442548 -6.4139 1.418e-10 ***
    ## prop_associate  0.795818   0.845979  0.9407  0.346855    
    ## prop_college   -1.721855   0.808532 -2.1296  0.033204 *  
    ## prop_gradprof   4.474180   1.427493  3.1343  0.001723 ** 
    ## rho_1           0.643665   0.023906 26.9249 < 2.2e-16 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

Note that the estimates are near-identical to the ones reported by
Liesenfeld et al. (2016).
