#' @title Kernel (h-mode) depth for functional data
#' 
#' @description Compute the kernel depth values of in-sample and out-of-sample functional data over different quantile levels
#' 
#' @param X An n by p matrix of in-sample functional data, where each row represent one observed curve. The depth is constructed by this dataset.
#' @param X0 An n0 by p matrix of in-sample functional data, where each row represent one observed curve. The depth is evaluated at these data points. If X0=NULL, only the in-sample depth values are provided.
#' @param u_vec A numerical vector of quantile levels to determine the bandwidths for kernel depth.
#' @param multi A logical vector. If TRUE, it is applied for functional data. Otherwise, it is applied for multivariate data.
#' 
#' @return A list of kernel depth values and selected bandwidths.
#' \item{KD}{The kernel depth values of both in-sample and out-of-sample functional data over the quantile levels in u_vec}
#' \item{h_vec}{The selected bandwidths given the quantile levels in u_vec}
#' 
#' @seealso 
#' \code{\link{RegDepth}}
#' \code{\link{depthInfer2test}}
#' \code{\link{depthInferFoFR}}
#' 
#' @references
#' Cuevas, Febrero, and Fraiman (2006) On the use of the bootstrap for estimating functions with functional data. Computational Statistics & Data Analysis, 51(2):1063--1074 
#' 
#' Wynne and Nagy (2025) Statistical depth meets machine learning: Kernel mean embeddings and depth in functional data analysis. International Statistical Review, 93(2):317--348
#' 
#' @examples
#' # example =========
#' 
#' # simulate the regressor curves from Karhunen-Loeve expansion and responses
#' J=15
#' 
#' # eigenvalues
#' c=2; a=3.5; dt = c*(1:J)^(-a)
#' g1 = c*VGAM::zeta(a)
#' ga = c(g1, sapply(1:(J-1), function(j){g1 - sum(dt[1:j])}))
#' 
#' # eigenfunctions
#' tt = 70; tGrid = 1:tt/tt; J=15
#' phi= t(fda::fourier(tGrid, J))
#' 
#' n=101
#' xi = matrix(rnorm(n*J), n, J)
#' X = xi %*% (phi * sqrt(ga))
#' n0=50
#' xi0 = matrix(rnorm(n0*J), n0, J)
#' X0 = xi0 %*% (phi * sqrt(ga))
#' 
#' u_vec = c(0.5, 0.3, 0.05, 0.01,0.001)
#' 
#' resKD = KD(X, X0, u_vec)
#' @export
KD = function(X, X0, u_vec, multi=FALSE){
  # Only for functions observed at dense equi-spaced time grid points on [0,1].
  
  if(is.vector(X)){X=matrix(X,nrow=1)}
  
  n = nrow(X); tt = ncol(X); 
  if(multi){scal=1}else{scal = 1/tt}
  
  rownames(X) = paste0("i",1:n)
  colnames(X) = paste0("t",1:tt)
  # u_vec must be within [0,1]
  
  Xdist = sqrt(Rfast::Dist(X,sq=T)*scal)
  Xdist_vec = Xdist[upper.tri(Xdist)]
  h_vec = quantile(Xdist_vec, u_vec)
  
  KDin = vapply(h_vec, function(h){
    colMeans(dnorm(Xdist/h))/h
  }, numeric(length = ncol(Xdist)))
  if(is.vector(KDin)){KDin=matrix(KDin,nrow=1)}
  rownames(KDin) = rownames(X)
  
  if(!is.null(X0)){
    if(is.vector(X0)){X0=matrix(X0,nrow=1)}
    n0 = nrow(X0)
    XdistX0 = sqrt(proxy::dist(X, X0)^2*scal)
    KDout = vapply(h_vec, function(h){
      colMeans(dnorm(XdistX0/h))/h
    }, numeric(length = ncol(XdistX0)))
    if(is.vector(KDout)){KDout=matrix(KDout,nrow=1)}
    rownames(KDout) = paste0("i0_",1:n0)
    KD = rbind(KDin, KDout)
  }else{
    KD = KD
  }
  return(list(
    KD = KD,
    h_vec = h_vec
  ))
}





