#' @title Depth-based inference for mean response in function-on-function regression models
#' 
#' @description
#' Conduct mean response inference for function-on-function regression (Yeon, 2026+)
#' with kernel (h-mode) depth, regularized halfspace depth, and regularized projection depth.
#' 
#' @param X An n by pX matrix of functional regressors. Each row represent one observed curve.
#' @param Y An n by pY matrix of functional responses. Each row represent one observed curve.
#' @param h_max An integer. The largest number of truncation levels that are used for test statistics based on functional principal components.
#' @param rho_vec A numeric vector of threshold values used for fraction of variance explained. These must be within zero and one.
#' @param u_vec A numerical vector of quantile levels to determine the regularizations for regularized halfspace and projection depths and the bandwidths for kernel depth.
#' @param Mproj An integer. The number of random projections used to approximate the regularized halfspace and projection depths.
#' @param Mbts An integer. The number of bootstrap resamples used for depth-based inference.
#' 
#' @return A list containing the p-values of the following statistics. All p-values are given across new functional predictors X0 and threshold levels in rho_vec.
#' \item{RD}{The depth p-values from the regularized halfspace depth (HD), the regularized projection depth (PD), and the tie-broken regularized halfspace depth (HDtb) over quantile levels in u_vec. The regularization was done by both standard deviation and median absoluted deviation of the projecitons.}
#' \item{KD}{The depth p-value from the kernel depth over quantile levels in u_vec.}
#' \item{OneD}{The p-values from the classical summary test statistics: the squared L2 norm (L2) and supremum norm (sup)}
#' 
#' @seealso 
#' \code{\link{KD}}
#' \code{\link{RegDepth}}
#' \code{\link{depthInfer2test}}
#' 
#' @references
#' Yeon, H. (2026) Inference for function-on-function regression: central limit theorem and residual bootstrap. Statistics and Its Interface. To appear
#' 
#' Yeon, H. (2026+) Effective and flexible depth-based inference for functional parameters. In preparation
#' 
#' @examples
#' 
#' ev_poly = function(a, c, J){
#'   dt = c*(1:J)^(-a)
#'   ld1 = c*VGAM::zeta(a)
#'   ld = c(ld1, sapply(1:(J-1), function(j){ld1 - sum(dt[1:j])}))
#'   return(ld)
#' }
#' 
#' #'set.seed(20260212)
#' tt = 50; tGrid_ttt = seq(0, 1, len=tt+1);
#' tGrid = tGrid_ttt[1:tt]
#' diffrange = diff(tGrid)[1] + diff(range(tGrid))       # when using left points
#' scal = diffrange / tt   # scaling factor for integration
#' 
#' Jtrue=20 # rows are functions in phi
#' 
#' # orthonormal systems
#' 
#' # trigonometric functions
#' phi_base = t(fda::fourier(tGrid_ttt, 2*Jtrue))[,1:tt]
#' phi = phi_base[(1:Jtrue),]
#' # orthonormalized monomials
#' mono = t(sapply(1:Jtrue, function(j)tGrid^j))
#' phi_mono = DepthInfer:::orthoL2Equidense(mono, tGrid)
#' # orthonormalized Chebyshev functions
#' cheb = t(sapply(1:Jtrue, function(j)pracma::chebPoly(j, 2*tGrid-1)))
#' phi_cheb = DepthInfer:::orthoL2Equidense(cheb, tGrid)
#' 
#' # eigenvalues
#' a1 = 2.5; a2 = 3.5
#' ga1  = ev_poly(a1, 2, Jtrue) # a=2.5
#' ga2  = ev_poly(a2, 2, Jtrue) # a=3.5
#' 
#' # slope operators
#' 
#' J0 = 5; Bscal = 2;
#' Wbetaj = rbinom(2*J0, 1, 0.5)*2-1
#' 
#' Bnull1 = Reduce("+",lapply(J0+1:J0, function(j){
#'   Bscal*j^(-(1.5)) *Wbetaj[j]* outer(phi_mono[j,], phi_mono[j,])
#' }))
#' Balter1 = Reduce("+",lapply(1:(2*J0), function(j){
#'   Bscal*j^(-(1.5)) *Wbetaj[j]* outer(phi[j,], phi[j,])
#' }))
#' Bnull2 = Reduce("+",lapply(J0+1:J0, function(j){
#'   Bscal*j^(-(2.5)) *Wbetaj[j]* outer(phi_mono[j,], phi_mono[j,])
#' }))
#' Balter2 = Reduce("+",lapply(1:(2*J0), function(j){
#'   Bscal*j^(-(2.5)) *Wbetaj[j]* outer(phi[j,], phi[j,])
#' }))
#' c_vec = seq(0,1,by=0.2)   # degrees of alternative
#' 
#' rexp_cent1 = function(n){rexp(n)-1} # centered exponential distribution
#' simX = function(n, Jtrue, rdist_Wj, rdist_xi, ga, phi){
#'   WjX = matrix(rdist_Wj(n*Jtrue), n, Jtrue)
#'     # xiX = rexp(n) - 1
#'   xiX = rdist_xi(n)
#'   Xfpc = t(sqrt(ga) * t(xiX*WjX))
#'   X = Xfpc %*% phi
#'   colnames(X) = paste0("t", 1:ncol(phi))
#'   rownames(X) = paste0("i", 1:n)
#'   return(X)
#' }
#' 
#' # ==== set scenarios ====
#' n_vec = c(50, 200, 1000)
#' n = n_vec[1]; iiB=1; ii_c=6
#' if(iiB==1){ Bnull=Bnull1;Balter=Balter1
#' }else if(iiB==2){ Bnull=Bnull2;Balter=Balter2
#' }
#' c = c_vec[ii_c]
#' B = (1-c)*Bnull + c*Balter
#' colnames(B) = paste0("t",1:ncol(B))
#' rownames(B) = paste0("t",1:nrow(B))
#' 
#' evX=xiX=evEr=xiEr=1
#' if(evX==1){ gaX=ga1 # a=2.5
#' }else{ gaX=ga2 # a=3.5
#' }
#' if(xiX==1){ rX=rnorm # \xi \sim \nd(0,1)
#' }else{ rX=rexp_cent1 # \xi \sim \Exp(1)-1
#' }
#' if(evEr==1){ gaEr=ga1 # a=2.5
#' }else{ gaEr=ga2 # a=3.5
#' }
#' if(xiEr==1){ rEr=rnorm # \xi \sim \nd(0,1)
#' }else{ rEr=rexp_cent1 # \xi \sim \Exp(1)-1
#' }
#' 
#' 
#' # ========= data generation ========
#' X = simX(n, Jtrue, rnorm, rX, gaX, phi_mono)
#' er = simX(n, Jtrue, rnorm, rEr, gaEr, phi_cheb)
#' BX = t(B %*% t(X) * scal)
#' Y = BX + er
#' 
#' # new predictors
#' n0=1 # It can be more than 1
#' X0 = simX(n0, J0, rnorm, rX, gaX[1:J0], phi_mono[1:J0,])
#' 
#' res = depthInferFoFR(X, Y, X0, Mproj=15)
#' 
#' @export
depthInferFoFR = function(
    X, Y, X0, 
    h_max = 20, rho_vec = c(seq(0.85, 0.95, by=0.05), 0.99),
    u_vec = c(0.5, 0.4, 0.3, 0.2, 0.1, 0.05, 0.01, 0.001), 
    Mproj = 2e4, Mbts=1000
    ){
  rho_names = paste0("rho=", rho_vec)
  names(rho_vec) = rho_names
  Hnames = rho_names
  if(is.vector(X0)){X0 = matrix(X0, nrow=1)}
  
  Xbar = colMeans(X); Xcent=t(t(X)-Xbar)
  Ybar = colMeans(Y); Ycent=t(t(Y)-Ybar)
  X0cent=t(t(X0)-Xbar)
  n=nrow(X); pX = ncol(X); pY=ncol(Y)
  
  GaHat = cov(X)*(n-1)/n
  
  resFPCA = DepthInfer:::FPCAcov(GaHat, rho_vec, h_max)
  gaHat = resFPCA$ev; phiHat = resFPCA$ef
  hFVE = c(resFPCA$hFVE)
  
  # (leave-one-out) cross validation ====
  resCV = DepthInfer:::kFoldCVcpp(X,Y,h_max,5)
  gHat = which.min(resCV)
  tune = c(gHat, hFVE)
  names(tune) = c("CV", rho_names)
  hFVE = pmax(hFVE, gHat)
  
  # Inference ====
  XprojCent = Xcent %*% phiHat * scal
  X0projCent = X0cent %*% phiHat * scal
  
  resDataVer = DepthInfer:::fofrMRtest(
    Y, X, Ybar, XprojCent, X0projCent, gaHat, phiHat, h_max
  )
  
  TstatArr=resDataVer$TstatArr
  YhatArr=resDataVer$YhatArr
  erHatArr=resDataVer$erHatArr
  tHatMat = resDataVer$tHatMat
  BhatX0cent = resDataVer$BhatX0cent
  BhatX0cent_cent = array(BhatX0cent[,,gHat], dim=dim(BhatX0cent))
  
  TstatArr = TstatArr[,,hFVE,drop=FALSE]
  dimnames(TstatArr) = list(
    dimnames(X0)[[1]], 
    dimnames(B)[[1]],
    Hnames
  )
  
  # ==========resample=========
  
  TstatStarListOrigin = DepthInfer:::fofrMRtestBTS(
    erHatArr[,,gHat], YhatArr[,,gHat],
    X, X0projCent, gaHat, phiHat, tHatMat, BhatX0cent_cent,
    Mbts
  )
  TstatStarList = DepthInfer:::AtoB(TstatStarListOrigin)
  
  for(i0 in 1:n0){
    TstatStarArr = TstatStarList[[i0]]
    TstatStarList[[i0]] = TstatStarArr[,hFVE,,drop=F]
  }
  
  H = dim(TstatArr)[[3]]
  u_len = length(u_vec)
  pvalRegDepth = array(
    0, dim=c(n0, H, u_len, 3, 2),
    dimnames = list(
      rownames(X0),Hnames,
      paste0("u=",u_vec),
      c("HD", "PD", "HDtb"),
      c("SD", "MAD")
    )
  )
  pvalKD = array(
    0, dim=c(n0, H, u_len),
    dimnames = list(
      rownames(X0),Hnames,
      paste0("u=",u_vec)
    )
  )
  pvalBTS1d = array(
    0, dim=c(n0, H, 2),
    dimnames = list(
      rownames(X0),Hnames,
      c("L2","sup")
    )
  )
  
  for(i0 in 1:n0){
    for(h in 1:H){
      # i0=1
      # h=1
      TnInner = matrix(TstatArr[i0, , h],nrow=1)
      TnStarInner = t(TstatStarList[[i0]][,h,])
      
      L2Star = rowSums(TnStarInner^2)*scal
      supStar = matrixStats::rowMaxs(TnStarInner)
      pvalBTS1d[i0,h,1] = mean(L2Star>sum(TnInner^2)*scal)
      pvalBTS1d[i0,h,2] = mean(supStar>max(TnInner))
      
      resKD = KD(TnStarInner,TnInner,u_vec)
      resRegDepth = RegDepth(
        TnStarInner,
        TnInner,
        u_vec,
        Mproj,
        verbose = TRUE
      )
      
      # tie breaks in RHD
      PO_mat = resRegDepth$PO
      HP_mat = resRegDepth$HP
      
      rkHDnew = rkHDold = resRegDepth$depth[,,,"HD"]
      RHD_BTS_pval = matrix(0,length(u_vec),2)
      for(ii_reg in 1:2){
        # ii_reg=1
        dRank = matrixStats::colRanks(
          resRegDepth$depth[,,ii_reg,"HD"],
          ties.method="min"
        )
        
        for(ii_u in 1:length(u_vec)){
          
          
          HPminIndices = DepthInfer:::rowMins(
            HP_mat[,as.logical(resRegDepth$whichReg[,ii_u,ii_reg]),drop=FALSE]
          )$indices
          
          # break ties in RHD rankings
          rk_new = rk = dRank[ii_u,]
          is_tie_vec = duplicated(rk)|duplicated(rk, fromLast = T)
          if(any(is_tie_vec)){
            
            rk_tied = unique(rk[is_tie_vec])
            for(ii_rk_tied in 1:length(rk_tied)){
              
              # ii_rk_tied = 2
              which_rk_tied = which(rk==rk_tied[ii_rk_tied])
              
              # Find outlyingness for projections of tied elements
              out_proj_vec = rep(0,length(which_rk_tied))
              for(ii_which_rk_tied in 1:length(which_rk_tied)){
                # ii_which_rk_tied=3
                idx_inner = which_rk_tied[ii_which_rk_tied]
                
                # HP_mat[idx_inner,c(HPminIndices[[idx_inner]])]
                out_proj_vec[ii_which_rk_tied] = 
                  mean(PO_mat[idx_inner,c(HPminIndices[[idx_inner]])])
                
              }
              # ranking tied elements based on 1/(1+outlyingness)
              detied_rk = rank(out_proj_vec, ties.method = "min")
              rk_new[which_rk_tied] = detied_rk+rk_tied[ii_rk_tied]-1
            }
          }
          
          rkHDnew[,ii_u,ii_reg]= rk_new
          rkHDold[,ii_u,ii_reg]= rk
          
          RHD_BTS_pval[ii_u,ii_reg] = mean(rk_new[Mbts+1]>rk_new[1:Mbts])
        }
      }
      
      tempp = rep(
        resRegDepth$depth[Mbts+1,,,,drop=F], each = Mbts
      )
      dim(tempp) = c(Mbts,dim(resRegDepth$depth)[-1])
      RegDepthBTSpval = apply(
        resRegDepth$depth[-(Mbts+1),,,]<tempp, 2:4, mean
      )
      
      pvalRegDepth[i0,h,,,] = abind::abind(
        RegDepthBTSpval,
        HDtb = RHD_BTS_pval
      )
      
      pvalKD[i0,h,] = rowMeans(t(resKD$KD[-(Mbts+1),])<resKD$KD[Mbts+1,])
    }
  }
  
  
  
  # results
  
  list(
    RD = pvalRegDepth,
    KD = pvalKD,
    OneD = pvalBTS1d,
    tune=tune
  )
  
}










