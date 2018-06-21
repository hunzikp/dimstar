Rcpp::loadModule("kernel", TRUE)

scipy <- NULL
.onLoad <- function(libname, pkgname) {
  scipy <<- reticulate::import("scipy", delay_load = TRUE)
}