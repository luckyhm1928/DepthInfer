#' @keywords internal
#' 
#' @title Depth-based inference for functional parameters
#' 
#' @author Hyemin Yeon \email{hyeon1@@kent.edu} 
#' 
#' @description
#' The package provides tools to test for functional (infinite-dimensional) parameters
#' focusing mainly on two-sample functional mean tests and mean respones inference for function-on-function regression
#' and using depth statistics proposed by Yeon (2026+).
#' 
#' \code{depthInfer2test} conducts two-sample functional mean tests using depth statistics.
#' 
#' \code{depthInferFoFR} conducts mean response inference for function-on-function regression.
#' 
#' \code{RegDepth} computes the regularized halfspace and projection depth values of functional data.
#' 
#' \code{KD} computes the kernel depth values of functional data.
#' 
#' 
#' 
#' @docType package
#' @name DepthInfer
#' @useDynLib DepthInfer, .registration = TRUE
#' @import Rcpp
#' 
#' @seealso
#' \code{\link{depthInfer2test}}
#' \code{\link{depthInferFoFR}}
#' \code{\link{RegDepth}}
#' \code{\link{KD}}
#' 
#' @references 
#' Yeon, H. (2026+) Effective and flexible depth-based inference for functional parameters. In preparation
#' 
"_PACKAGE"

## usethis namespace: start
## usethis namespace: end
NULL
