#' @title Depth-based inference: two-sample functional mean tests
#' 
#' @description
#' Conduct two-sample functional mean tests using depth statistics (Yeon, 2026+)
#' with kernel (h-mode) depth, regularized halfspace depth, regularized projection depth,
#' integrated depth, and inifimal depth
#' 
#' @param X1 An n1 by p matrix of functional data of the first group. Each row represent one observed curve.
#' @param X2 An n2 by p matrix of functional data of the second group. Each row represent one observed curve.
#' @param h_max An integer. The largest number of truncation levels that are used for test statistics based on functional principal components.
#' @param rho_vec A numeric vector of threshold values used for fraction of variance explained. These must be within zero and one.
#' @param u_vec A numerical vector of quantile levels to determine the regularizations for regularized halfspace and projection depths and the bandwidths for kernel depth.
#' @param Mproj An integer. The number of random projections used to approximate the regularized halfspace and projection depths.
#' @param B An integer. The number of bootstrap resamples used for depth-based inference.
#' 
#' @return A list containing the p-values of the following statistics.
#' \item{ITD}{The depth p-value from the integrated depth}
#' \item{IFD}{The depth p-value from the infimal depth}
#' \item{RD}{The depth p-values from the regularized halfspace depth (HD), the regularized projection depth (PD), and the tie-broken regularized halfspace depth (HDtb) over quantile levels in u_vec. The regularization was done by both standard deviation and median absoluted deviation of the projecitons.}
#' \item{KD}{The depth p-value from the kernel depth over quantile levels in u_vec.}
#' \item{BTS1d}{The p-values from the classical summary test statistics: the squared L2 norm (L2normSq), supremum norm (SupNorm), integrated pointwise F statistic (Fint), and maximal pointwise F statistic (Fmax).}
#' \item{BTSfpc}{The p-values from the test statistics based on functional principal components analysis: QuadInv (Horvath et al., 2013) and QuadMul (Sharghi Ghale-Joogh and Hosseini-Nasab, 2018).}
#' 
#' @seealso 
#' \code{\link{KD}}
#' \code{\link{RegDepth}}
#' \code{\link{depthInferFoFR}}
#' 
#' @references
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
#' set.seed(20260212)
#' tt = 50; tGrid_ttt = seq(0, 1, len=tt+1); 
#' tGrid = tGrid_ttt[1:tt]
#' diffrange = diff(tGrid)[1] + diff(range(tGrid))       # when using left points
#' scal = diffrange / tt   # scaling factor for integration
#' 
#' Jtrue=20 # rows are functions in phi
#' 
#' # eigenvalues
#' ga_rough  = ev_poly(2.5, 2, Jtrue) # a=2.5
#' ga_smooth = ev_poly(5, 2, Jtrue)   # a=5
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
#' # orthonormalized spline functions
#' rangeval = c(0, 1); nbasis = Jtrue
#' bs_basis <- fda::create.bspline.basis(
#'   rangeval = rangeval, nbasis = nbasis, norder = 4
#' )  # norder = degree + 1
#' phi_splines = DepthInfer:::orthoL2Equidense(
#'   t(fda::eval.basis(tGrid, bs_basis)), tGrid
#' )
#' 
#' # alternative mean functions
#' mu1Alter_mat = 1*rbind(
#'   mag = rep(1, tt), # mag
#'   jump = sapply(tGrid, function(t){ifelse(t<=0.2, -1, 1)}), # jump
#'   peak = sapply(tGrid, function(t){ifelse(t>0.2&t<=0.4, 1, -1)}), # peak
#'   lin = 2*tGrid-1, # lin
#'   quad = 8*(tGrid-1/2)^2-1, # quad
#'   cubic = 12*sqrt(3)*tGrid*(tGrid-1/2)*(tGrid-1), # cubic
#'   wiggle = sin(10*pi*(tGrid-0.05)) # wiggle
#' )
#' 
#' c_vec = seq(0,1,by=0.2)
#' mu1 = rep(0,tt); c_alter=c_vec[6]; Htype=4 # (linear alternative)
#' mu2 = (1-c_alter)*mu1 + c_alter*mu1Alter_mat[Htype,]
#' 
#' eveq = TRUE # TRUE: eigenvalues are equal
#' efeq = TRUE # TRUE: eigenfunctions are equal
#' if(eveq){ ga1 = ga2 = ga_rough # equal eigenvalues 
#' }else{ ga1 = ga_smooth; ga2 = ga_rough # unequal eigenvalues
#' }
#' if(efeq){ phi1 = phi; phi2 = phi # equal eigenfunctions
#' }else{ phi1 = phi_mono; phi2 = phi_cheb # unequal eigenfunctions
#' }
#' 
#' n=50; n1 = n2 = n/2 
#' idx1 = as.logical(c(rep(1,n1), rep(0,n2))); idx2 = !idx1
#' 
#' 
#' # # # # # # # # # # # 
#' # generate data
#' 
#' Wj = matrix(rnorm(n*Jtrue), n, Jtrue); xi_type=1
#' if(xi_type==2){ xi = rexp(n) - 1 # normal*exp FPC scores
#' }else if(xi_type==1){ xi = rnorm(n) # normal*normal FPC scores
#' }else{ xi = rep(1,n) # independent FPC scores
#' }
#' 
#' Xfpc = matrix(0, n, Jtrue)
#' Xfpc[idx1,] = t(sqrt(ga1) * t(xi[idx1]*Wj[idx1,]))
#' Xfpc[idx2,] = t(sqrt(ga2) * t(xi[idx2]*Wj[idx2,]))
#' 
#' X = matrix(0, n, tt)
#' X[idx1,] = t(mu1 + t(Xfpc[idx1,] %*% phi1))
#' X[idx2,] = t(mu2 + t(Xfpc[idx2,] %*% phi2))
#' 
#' X1 = X[idx1,]; X2 = X[idx2,]
#' 
#' res = depthInfer2test(X1, X2, Mproj=15)
#' 
#' @export
depthInfer2test = function(
    X1, X2,
    h_max = 20, rho_vec = c(seq(0.85, 0.95, by=0.05), 0.99),
    u_vec = c(0.5, 0.4, 0.3, 0.2, 0.1, 0.05, 0.01, 0.001), 
    Mproj = 2e4, B=1000
    ){
  X = rbind(X1, X2)
  n1 = nrow(X1); n2 = nrow(X2)
  gp_vec = c(rep(1, n1), rep(2,n2))
  tt = ncol(X)
  idx1 = as.logical(c(rep(1,n1), rep(0,n2))); idx2 = !idx1
  
  # statistics
  
  DtHat = colMeans(X1) - colMeans(X2)
  
  ftn2stat = function(X, gp_vec){
    
    idx1 = ifelse(gp_vec==unique(gp_vec)[1], TRUE, FALSE)
    idx2 = !idx1; n=nrow(X)
    X1 = X[idx1,]; X2 = X[idx2,]
    X1bar = colMeans(X1); X2bar = colMeans(X2)
    
    Xcent=X
    Xcent[idx1,]=t(t(X1)-X1bar)
    Xcent[idx2,]=t(t(X2)-X2bar)
    DtHat = X1bar-X2bar
    
    # pooled covariance estimation
    GaHat_pool = (cov(X1)*(n1-1)+cov(X2)*(n2-1)) / n
    eig_pool = eigen(GaHat_pool)
    gaHat_pool = eig_pool$values[1:h_max] * scal
    phiHat_pool = eig_pool$vectors[,1:h_max] / sqrt(scal)
    
    DtHat_j_vec = c(DtHat %*% phiHat_pool * scal)
    
    gaHatPDsum = cumsum(gaHat_pool)
    FVE_PD = gaHatPDsum/gaHatPDsum[h_max]
    # XprojCent_pool = Xcent %*% phiHat_pool * scal
    
    
    # ============stat=======
    
    # bootstrap tests (6 stats)
    
    L2normSq = mean(DtHat^2)                   # L2 norm
    SupNorm = max(abs(DtHat))                  # sup norm
    QuadInv = cumsum(DtHat_j_vec^2/gaHat_pool) # T2
    QuadMul = cumsum(DtHat_j_vec^2*gaHat_pool) # \|\ga^{1/2}v\|^2
    
    Xbar = colMeans(X)
    SSEvec = colSums(t(t(X1) - X1bar)^2)+colSums(t(t(X2) - X2bar)^2)
    SSRvec = n1*(X1bar-Xbar)^2 + n2*(X2bar-Xbar)^2
    Fvec = (SSRvec / 1) / (SSEvec / (n-2))
    Fint = sum(Fvec)*scal # Fint, Zhang: fdANOVA::fanova.tests(t(X), gp_vec, test="GPF")
    Fmax = max(Fvec)      # Fmax, Zhang: fdANOVA::fanova.tests(t(X), gp_vec, test="Fmaxb")
    
    res = list(
      statFPC = rbind(FVE_PD = FVE_PD, QuadInv=QuadInv, QuadMul=QuadMul),
      stat1d = c(L2normSq=L2normSq, SupNorm=SupNorm, Fint=Fint, Fmax = Fmax)
    )
    return(res)
  }
  res = ftn2stat(X, gp_vec)
  
  hSelFVE = sapply(rho_vec, function(rho){
    min(which(res$statFPC["FVE_PD",]>rho))
  })
  resStatFPC = cbind(res$statFPC[2:3,],res$statFPC[2:3,hSelFVE])
  colnames(resStatFPC) = c(
    paste0("h=",1:h_max),
    paste0("rho=", rho_vec)
  )
  
  # bootstrap 
  Xbar_pool = colMeans(X); XStar = X
  
  erHat1=t(t(X[idx1,])-colMeans(X[idx1,]))
  erHat2=t(t(X[idx2,])-colMeans(X[idx2,]))
  DtHatStar_mat = matrix(0, B, tt)
  
  resFPC = array(
    0,c(2, h_max+length(rho_vec), B),
    dimnames=list(
      c("QuadInv", "QuadMul"),
      c(paste0('h=',1:h_max), paste0("rho=", rho_vec)),
      NULL
    ))
  res1d = matrix(0,4, B)
  rownames(res1d) = c("L2normSq", "SupNorm", "Fint", "Fmax")
  for(b in 1:B){
    idxBTS1 = ceiling(n1*runif(n1))
    idxBTS2 = ceiling(n2*runif(n2))
    
    Xcent1Star = erHat1[idxBTS1,]
    Xcent2Star = erHat2[idxBTS2,]
    XStar[idx1,] = t(Xbar_pool+t(Xcent1Star))
    XStar[idx2,] = t(Xbar_pool+t(Xcent2Star))
    X1Star = XStar[idx1,]; X2Star = XStar[idx2,] 
    
    DtHatStar_mat[b,] = colMeans(X1Star) - colMeans(X2Star)
    
    # 
    resStar = ftn2stat(XStar, gp_vec)
    res1d[,b] = resStar$stat1d
    resFPC[,,b] = cbind(
      resStar$statFPC[2:3,],
      resStar$statFPC[2:3,hSelFVE]
    )
  }
  
  resBTS1d = rowMeans(res1d>res$stat1d)
  resBTSfpc = t(sapply(rownames(resFPC), function(aa){rowMeans(resFPC[aa,,]>resStatFPC[aa,])}))
  
  # depth-based inference
  
  # 1) ITD/IFD
  DtHatStar_list = apply(
    rbind(DtHatStar_mat, DtHat), 
    1, function(x){list(args = tGrid, vals = x)}
  )
  depthStar = ddalpha::depthf.fd1(
    DtHatStar_list, 
    DtHatStar_list[1:B]
  )
  
  depth_pval_ITD = mean(
    depthStar$Half_FD[B+1]>depthStar$Half_FD[1:B]
  )
  depthStar_IFD_rank = ddalpha::infimalRank(
    depthStar$Half_ID, # tie breaking for IFD
    depthStar$Half_IA
  )
  depth_pval_IFD = mean(
    depthStar_IFD_rank[B+1]<depthStar_IFD_rank[1:B]
    )
  
  # Regularized halfspace/projection depths
  
  resRegDepth = RegDepth(
    DtHatStar_mat,
    DtHat,
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
      
      RHD_BTS_pval[ii_u,ii_reg] = mean(rk_new[B+1]>rk_new[1:B])
    }
  }
  
  
  tempp = rep(
    resRegDepth$depth[B+1,,,,drop=F], each = B
  )
  dim(tempp) = c(B,dim(resRegDepth$depth)[-1])
  RegDepthBTSpval = apply(
    resRegDepth$depth[-(B+1),,,]<tempp, 2:4, mean
  )
  RDpval = abind::abind(
    RegDepthBTSpval,
    HDtb = RHD_BTS_pval,
    along=3
  )
  
  # kernel depth
  
  resKD = KD(DtHatStar_mat,DtHat,u_vec)
  KDpval = rowMeans(t(resKD$KD[-(B+1),])<resKD$KD[B+1,])
  
  # results
  list(
    ITD = depth_pval_ITD,
    IFD = depth_pval_IFD,
    RD = RDpval,
    KD = KDpval,
    BTS1d = resBTS1d,
    BTSfpc = resBTSfpc
  )
}





