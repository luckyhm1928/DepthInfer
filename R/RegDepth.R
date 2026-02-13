#' @title Regarulized halfspace and projection depths for functional data
#' 
#' @description Compute the regularized halfspace and projection depth values of in-sample and out-of-sample functional data over different quantile levels
#' 
#' @param X An n by p matrix of in-sample functional data, where each row represent one observed curve. The depth is constructed by this dataset.
#' @param X0 An n0 by p matrix of in-sample functional data, where each row represent one observed curve. The depth is evaluated at these data points. If X0=NULL, only the in-sample depth values are provided.
#' @param u_vec A numerical vector of quantile levels to determine the regularizations.
#' @param Mproj An integer. The total number of random projections used to approximate the regularized halfspace and projection depths. For each quantile level u in u_vec, (1-u)*Mproj directions will be used for constructing the depth functions.
#' @param verbose A logical vector. If FALSE, only depth and RegPara appear in the result. Otherwise, it also shows shows regularized indices, projections of the functional data, halfspace probabilties, projection outlyingness, generated projection directions, standard deviations and median absolute deviations of the projected data.
#' 
#' @return 
#' \item{depth}{The regularized halfspace and projection depth values of both in-sample and out-of-sample functional data over the quantile levesl in u_vec by regularizing either standard deviations or median absolute deviations of the projected data}
#' \item{RegPara}{The selected regularizations by the quantile levels in u_vec.}
#' \item{whichReg}{The indices of the regularized projections.}
#' \item{Xproj}{The projected functional data onto the generated random directions.}
#' \item{HP}{The halfspace probabilities.}
#' \item{PO}{The projection outlyingness.}
#' \item{ProjDir}{The randomly generated directions.}
#' \item{ProjDispersion}{The standard deviations and median absolute deviations of the projected data.}
#' 
#' @seealso 
#' \code{\link{KD}}
#' \code{\link{depthInfer2test}}
#' \code{\link{depthInferFoFR}}
#' 
#' @references
#' Yeon, H., Dai, X., and Lopez-Pintado, S. (2025) Regularized halfspace depth for functional data. Journal of the Royal Statistical Society Series B: Statistical Methodology, 87(5):1553--1575
#' 
#' Bočinec, F., Nagy, S., and Yeon, H. (2026+) Projection detph for functional data: Theoretical properties. Under major revision
#' 
#' Yeon, H. (2026+) Effective and flexible depth-based inference for functional parameters. In preparation
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
#' resRD = RegDepth(X, X0, u_vec)
#' 
#' @export
RegDepth = function(X, X0, u_vec, Mproj=1e4, verbose=FALSE){
  # Only for functions observed at dense equi-spaced time grid points on [0,1].
  
  if(is.vector(X)){X=matrix(X,nrow=1)}
  if(is.vector(X0)){X0=matrix(X0,nrow=1)}
  
  n = nrow(X); tt = ncol(X); scal = 1/tt
  rownames(X) = paste0("i",1:n)
  colnames(X) = paste0("t",1:tt)
  
  if(!is.null(X0)){ # Outliers ====
    n0 = nrow(X0)
    rownames(X0) = paste0("i0_",1:n0)
    colnames(X0) = paste0("t",1:tt)
  }
  
  # Generate projections =====
  # generate from non-isotropic Gaussian distribution
  z_mat = mvtnorm::rmvnorm(Mproj, sigma=cov(X)) 
  
  v_mat = z_mat / sqrt(rowSums(z_mat^2)*scal) # normalize
  rownames(v_mat) = paste0("v", 1:Mproj)
  
  # Compute halfspace probabilities and projection outlyingness =====
  Xproj_mat = X %*% t(v_mat) * scal
  if(!is.null(X0)){ # Outliers ====
    X0proj_mat = X0 %*% t(v_mat) * scal
  }
  SDproj_vec = matrixStats::colSds(Xproj_mat)
  MADproj_vec = matrixStats::colMads(Xproj_mat)/1.4826
  MedProj_vec = matrixStats::colMedians(Xproj_mat)
  POin_mat = 1/(1+abs(t((t(Xproj_mat) - MedProj_vec) / MADproj_vec)))
  HPin_mat = (n-t(matrixStats::colRanks(Xproj_mat, ties.method = "max"))+1)/n
  if(!is.null(X0)){ # Outliers ====
    POout_mat = 1/(1+abs(t((t(X0proj_mat) - MedProj_vec) / MADproj_vec)))
    HPout_mat = X0proj_mat
    for(ii_v in 1:Mproj){
      HPout_mat[, ii_v] = colMeans(outer(
        Xproj_mat[, ii_v], X0proj_mat[, ii_v], FUN = ">="
      ))
    }
  }
  
  u_vec = sort(u_vec, dec=T)
  u_len = length(u_vec)
  ldSDvec = quantile(SDproj_vec, u_vec)
  ldMADvec = quantile(MADproj_vec, u_vec)
  
  # compute depth values ========
  resInArr = array(
    0, dim=c(n,u_len, 2, 2),
    dimnames=list(
      paste0("i",1:n),
      paste0("u=",u_vec),
      c("SD", "MAD"),
      c("HD", "PD")
    )
  )
  resOutArr = array(
    0, dim=c(n0,u_len, 2, 2),
    dimnames=list(
      paste0("i0_",1:n0),
      paste0("u=",u_vec),
      c("SD", "MAD"),
      c("HD", "PD")
    )
  )
  whichSDreg_mat = whichMADreg_mat = matrix(
    0,Mproj,u_len,
    dimnames = list(
      paste0("v",1:Mproj),
      paste0("u=", u_vec)
    )
    )
  for(ii_u in 1:u_len){
    if(ii_u==1){
      # Inliers =====
      
      # SDreg
      whichSDreg_mat[,ii_u] = whichSDreg_past = SDproj_vec>ldSDvec[ii_u]
      HDinSDreg_vec_past = matrixStats::rowMins(
        HPin_mat[,whichSDreg_past,drop=FALSE]
        )
      PDinSDreg_vec_past = matrixStats::rowMins(
        POin_mat[,whichSDreg_past,drop=FALSE]
        )
      
      resInArr[,ii_u,"SD","HD"] = HDinSDreg_vec_past
      resInArr[,ii_u,"SD","PD"] = PDinSDreg_vec_past
      
      # MADreg
      whichMADreg_mat[,ii_u] = whichMADreg_past = MADproj_vec>ldMADvec[ii_u]
      HDinMADreg_vec_past = matrixStats::rowMins(
        HPin_mat[,whichMADreg_past,drop=FALSE]
        )
      PDinMADreg_vec_past = matrixStats::rowMins(
        POin_mat[,whichMADreg_past,drop=FALSE]
        )
      
      resInArr[,ii_u,"MAD","HD"] = HDinMADreg_vec_past
      resInArr[,ii_u,"MAD","PD"] = PDinMADreg_vec_past
      
      if(!is.null(X0)){ # Outliers ====
        # SDreg
        HDoutSDreg_vec_past = matrixStats::rowMins(
          HPout_mat[,whichSDreg_past,drop=FALSE]
          )
        PDoutSDreg_vec_past = matrixStats::rowMins(
          POout_mat[,whichSDreg_past,drop=FALSE]
          )
        
        resOutArr[,ii_u,"SD","HD"] = HDoutSDreg_vec_past
        resOutArr[,ii_u,"SD","PD"] = PDoutSDreg_vec_past
        # MADreg
        HDoutMADreg_vec_past = matrixStats::rowMins(
          HPout_mat[,whichMADreg_past,drop=FALSE]
          )
        PDoutMADreg_vec_past = matrixStats::rowMins(
          POout_mat[,whichMADreg_past,drop=FALSE]
          )
        
        resOutArr[,ii_u,"MAD","HD"] = HDoutMADreg_vec_past
        resOutArr[,ii_u,"MAD","PD"] = PDoutMADreg_vec_past
      }
      
    }else{
      # SDreg
      whichSDreg_mat[,ii_u] = whichSDreg_current = SDproj_vec>ldSDvec[ii_u]
      whichSDreg_subtraction = whichSDreg_current & !whichSDreg_past
      if(all(!whichSDreg_subtraction)){ # reg sets are equal 
        # Inliners
        HDinSDreg_vec_current = HDinSDreg_vec_past
        PDinSDreg_vec_current = PDinSDreg_vec_past
        if(!is.null(X0)){ # Outliers ====
          HDoutSDreg_vec_current = HDoutSDreg_vec_past
          PDoutSDreg_vec_current = PDoutSDreg_vec_past
        }
      }else{
        # Inliers
        HDinSDreg_vec_current = pmin(
          matrixStats::rowMins(
            HPin_mat[,whichSDreg_subtraction,drop=F]
            ),
          HDinSDreg_vec_past
        )
        PDinSDreg_vec_current = pmin(
          matrixStats::rowMins(
            POin_mat[,whichSDreg_subtraction,drop=F]
            ),
          PDinSDreg_vec_past
        )
        if(!is.null(X0)){ # Outliers ====
          HDoutSDreg_vec_current = pmin(
            matrixStats::rowMins(
              HPout_mat[,whichSDreg_subtraction,drop=F]
              ),
            HDoutSDreg_vec_past
          )
          PDoutSDreg_vec_current = pmin(
            matrixStats::rowMins(
              POout_mat[,whichSDreg_subtraction,drop=F]
              ),
            PDoutSDreg_vec_past
          )
        }
      }
      
      # MADreg
      whichMADreg_mat[,ii_u] = whichMADreg_current = MADproj_vec>ldMADvec[ii_u]
      whichMADreg_subtraction = whichMADreg_current & !whichMADreg_past
      if(all(!whichMADreg_subtraction)){# reg sets are equal 
        # Inliers
        HDinMADreg_vec_current = HDinMADreg_vec_past
        PDinMADreg_vec_current = PDinMADreg_vec_past
        if(!is.null(X0)){ # Outliers ====
          HDoutMADreg_vec_current = HDoutMADreg_vec_past
          PDoutMADreg_vec_current = PDoutMADreg_vec_past
        }
      }else{
        # Inliers
        HDinMADreg_vec_current = pmin(
          matrixStats::rowMins(
            HPin_mat[,whichMADreg_subtraction,drop=F]
            ),
          HDinMADreg_vec_past
        )
        PDinMADreg_vec_current = pmin(
          matrixStats::rowMins(
            POin_mat[,whichMADreg_subtraction,drop=F]
            ),
          PDinMADreg_vec_past
        )
        if(!is.null(X0)){ # Outliers ====
          HDoutMADreg_vec_current = pmin(
            matrixStats::rowMins(
              HPout_mat[,whichMADreg_subtraction,drop=F]
              ),
            HDoutMADreg_vec_past
          )
          PDoutMADreg_vec_current = pmin(
            matrixStats::rowMins(
              POout_mat[,whichMADreg_subtraction,drop=F]
              ),
            PDoutMADreg_vec_past
          )
        }
      }
      
      # store
      # Inliers
      resInArr[,ii_u,"SD","HD"] = HDinSDreg_vec_current
      resInArr[,ii_u,"SD","PD"] = PDinSDreg_vec_current
      resInArr[,ii_u,"MAD","HD"] = HDinMADreg_vec_current
      resInArr[,ii_u,"MAD","PD"] = PDinMADreg_vec_current
      if(!is.null(X0)){ # Outliers ====
        resOutArr[,ii_u,"SD","HD"] = HDoutSDreg_vec_current
        resOutArr[,ii_u,"SD","PD"] = PDoutSDreg_vec_current
        resOutArr[,ii_u,"MAD","HD"] = HDoutMADreg_vec_current
        resOutArr[,ii_u,"MAD","PD"] = PDoutMADreg_vec_current
      }
      
      # update
      # Inliers
      HDinSDreg_vec_past = HDinSDreg_vec_current
      PDinSDreg_vec_past = PDinSDreg_vec_current
      HDinMADreg_vec_past = HDinMADreg_vec_current
      PDinMADreg_vec_past = PDinMADreg_vec_current
      if(!is.null(X0)){ # Outliers ====
        HDoutSDreg_vec_past = HDoutSDreg_vec_current
        PDoutSDreg_vec_past = PDoutSDreg_vec_current
        HDoutMADreg_vec_past = HDoutMADreg_vec_current
        PDoutMADreg_vec_past = PDoutMADreg_vec_current
      }
    }
  }
  
  # results ====
  if(!is.null(X0)){
    res = abind::abind(
      resInArr,
      resOutArr,
      along=1
    )
    Xproj = rbind(Xproj_mat,X0proj_mat)
    PO_mat = rbind(POin_mat, POout_mat)
    HP_mat = rbind(HPin_mat, HPout_mat)
  }else{
    res = resInArr
    Xproj = Xproj_mat
    PO_mat = POin_mat
    HP_mat = HPin_mat
  }
  if(verbose){
    return(
      list(
        depth=res, # depth values
        RegPara = rbind( # regularization parameters
          SD = ldSDvec,
          MAD = ldMADvec
        ),
        whichReg = abind::abind( # indices of projections that are regularized
          SD = whichSDreg_mat,
          MAD = whichMADreg_mat,
          along=3
        ),
        Xproj = Xproj, # projections
        HP = HP_mat, # halfspace probabilities
        PO = PO_mat, # projection outlyingness
        ProjDir = v_mat, # All projection directions
        ProjDispersion = rbind( # SD/MAD of projections
          SD = SDproj_vec,
          MAD = MADproj_vec
        )
      )
    )
  }else{
    return(
      list(
        depth=res, # depth values
        RegPara = rbind( # regularization parameters
          SD = ldSDvec,
          MAD = ldMADvec
        )
      )
    )
  }
}






