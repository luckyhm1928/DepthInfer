#include <RcppArmadillo.h>
#include <map>
#include <vector>
#include <algorithm>

using namespace Rcpp;
using namespace arma;

// [[Rcpp::depends(RcppArmadillo)]]


// [[Rcpp::export]]
arma::mat FourierBasis(const arma::vec& pts, int K) {
  int nGrid = pts.n_elem;
  arma::mat res(K, nGrid);
  
  for (int k = 1; k <= K; ++k) {
    if (k == 1) {
      res.row(k - 1).fill(1.0);  // f_1(x) = 1
    } else if (k % 2 == 0) {
      res.row(k - 1) = sqrt(2.0) * arma::trans(arma::sin(k * M_PI * pts));
    } else {
      res.row(k - 1) = sqrt(2.0) * arma::trans(arma::cos((k - 1) * M_PI * pts));
    }
  }
  
  return res;
}

//[[Rcpp::export]]
arma::mat orthoL2Equidense(arma::mat& X, arma::vec& tGrid){
  int n = X.n_rows;
  int tt = X.n_cols;
  double scal = (range(tGrid) + tGrid(1) - tGrid(0)) / tt;  // left 
  
  arma::mat Xon(n, tt, fill::zeros);
  
  arma::vec x1 = trans(X.row(0));
  Xon.row(0) = trans(x1/sqrt(sum(square(x1))*scal));
  
  for(int i=1; i<n; ++i){
    arma::vec x_i = trans(X.row(i));
    arma::mat X_i1 = Xon.rows(0, i-1);
    arma::vec v_i = x_i - trans(X_i1) * (X_i1 * x_i * scal);
    Xon.row(i) = trans(v_i/sqrt(sum(square(v_i))*scal));
    
  }
  return(Xon);
  
}

// [[Rcpp::export]]
List rowMins(const arma::mat& x) {
  int n_rows = x.n_rows;
  NumericVector min_values(n_rows);
  List min_indices(n_rows);
  
  for (int i = 0; i < n_rows; ++i) {
    arma::rowvec row = x.row(i);
    double min_val = row.min();
    arma::uvec indices = arma::find(row == min_val);
    
    min_values[i] = min_val;
    min_indices[i] = indices + 1;  // Convert to 1-based indexing
  }
  
  return List::create(
    _["min_values"] = min_values,
    _["indices"] = min_indices
  );
}



// functions for FoFR


// [[Rcpp::export]]
List FPCAcov(const arma::mat& GaHat, const arma::vec& rho_vec, int h_max) {
  int tt = GaHat.n_cols;
  double scal = 1.0 / tt;
  
  // Eigen decomposition
  arma::vec gaHat;
  arma::mat phiHat;
  eig_sym(gaHat, phiHat, GaHat);
  
  // Flip and scale
  gaHat = reverse(gaHat) * scal;
  phiHat = fliplr(phiHat) / std::sqrt(scal);
  
  // FVE and hFVE
  arma::vec FVE = cumsum(gaHat) / sum(gaHat);
  arma::mat FVE_mat = repmat(FVE, 1, rho_vec.n_elem);
  arma::mat rho_mat = repmat(rho_vec.t(), tt, 1);
  arma::umat mask = FVE_mat < rho_mat;
  arma::rowvec hFVE = conv_to<rowvec>::from(sum(mask, 0)) + 1;
  hFVE = clamp(hFVE, 1, h_max);
  
  return List::create(
    Named("ev") = gaHat,
    Named("ef") = phiHat,
    Named("FVE") = FVE,
    Named("hFVE") = hFVE
  );
}

// [[Rcpp::export]]
arma::vec kFoldCVcpp(const arma::mat& X, const arma::mat& Y, int h_max, int K) {
  int n = X.n_rows;
  int tt = X.n_cols;
  double scal = 1.0 / tt;
  arma::vec resCVkfold(h_max, fill::zeros);
  
  // Randomly shuffle indices
  arma::uvec indices = arma::randperm(n);
  
  // Determine fold sizes
  int fold_size = n / K;
  int remainder = n % K;
  
  for (int k = 0; k < K; ++k) {
    // Determine fold boundaries
    int start_idx = k * fold_size + std::min(k, remainder);
    int end_idx = start_idx + fold_size - 1;
    if (k < remainder) end_idx += 1;
    
    arma::uvec test_idx = indices.subvec(start_idx, end_idx);
    arma::uvec train_idx = arma::regspace<arma::uvec>(0, n - 1);
    train_idx.shed_rows(test_idx);
    
    arma::mat X_train = X.rows(train_idx);
    arma::mat Y_train = Y.rows(train_idx);
    arma::mat X_test = X.rows(test_idx);
    arma::mat Y_test = Y.rows(test_idx);
    
    // Centering
    arma::rowvec X_train_bar = mean(X_train, 0);
    arma::rowvec Y_train_bar = mean(Y_train, 0);
    arma::mat X_test_cent = X_test.each_row() - X_train_bar;
    
    int n_train = X_train.n_rows;
    
    // Covariance and eigen decomposition
    arma::mat GaHat_train = cov(X_train) * (n_train - 1.0) / n_train;
    arma::vec eigval;
    arma::mat eigvec;
    eig_sym(eigval, eigvec, GaHat_train);
    
    eigval = reverse(eigval);
    eigvec = fliplr(eigvec);
    
    arma::vec gaHat_train = eigval * scal;
    arma::mat phiHat_train = eigvec / std::sqrt(scal);
    
    arma::mat DtHat_train = cov(Y_train, X_train) * (n_train - 1.0) / n_train;
    arma::mat DtHat_phiHat_train = DtHat_train * phiHat_train * scal;
    
    arma::mat X_test_proj = X_test_cent * phiHat_train * scal;
    
    for (int i = 0; i < X_test.n_rows; ++i) {
      arma::rowvec BHatXicent(tt, fill::zeros);
      arma::vec PEiSq_vec(h_max, fill::zeros);
      
      for (int j = 0; j < h_max; ++j) {
        arma::rowvec tempBHatXicent = (X_test_proj(i, j) / gaHat_train[j]) * DtHat_phiHat_train.col(j).t();
        BHatXicent += tempBHatXicent;
        arma::rowvec muHatXi = Y_train_bar + BHatXicent;
        PEiSq_vec[j] = accu(square(Y_test.row(i) - muHatXi)) * scal;
      }
      
      resCVkfold += PEiSq_vec;
    }
  }
  
  return resCVkfold / n;
}

// [[Rcpp::export]]
List fofrMRtest(
    const arma::mat& Y, const arma::mat& X, const arma::rowvec& Ybar,
    const arma::mat& XprojCent, const arma::mat& X0projCent,
    const arma::vec& gaHat, const arma::mat& phiHat,
    int h_max
    ) {
  int n = Y.n_rows;
  int n0 = X0projCent.n_rows;
  int tt = Y.n_cols;
  double scal = 1.0 / Y.n_cols;
  
  // DtHat and DtHat_phiHat
  arma::mat DtHat = cov(Y, X) * (n - 1.0) / n;
  arma::mat DtHat_phiHat = DtHat * phiHat * scal;
  
  // Initialize containers
  arma::cube YhatArr(n, tt, h_max, fill::zeros);
  arma::cube erHatArr(n, tt, h_max, fill::zeros);
  arma::cube TstatArr(n0, tt, h_max, fill::zeros);
  
  arma::mat BhatXcent(n, tt, fill::zeros);
  arma::mat BhatX0cent(n0, tt, fill::zeros);
  arma::cube BhatX0centCube(n0, tt, h_max, fill::zeros);
  arma::vec tHat(n0, fill::zeros);
  arma::mat tHatMat(n0, h_max, fill::zeros);
  
  for (int h = 0; h < h_max; ++h) {
    int j = h;
    
    // BhatXcent
    arma::mat tempBhatXcent = (XprojCent.col(j) / gaHat[j]) * DtHat_phiHat.col(j).t();
    BhatXcent += tempBhatXcent;
    
    // Yhat and error
    arma::mat Yhat = repmat(Ybar, n, 1) + BhatXcent;
    YhatArr.slice(h) = Yhat;
    erHatArr.slice(h) = Y - Yhat;
    
    // BhatX0cent
    arma::mat tempBhatX0cent = (X0projCent.col(j) / gaHat[j]) * DtHat_phiHat.col(j).t();
    BhatX0cent += tempBhatX0cent;
    BhatX0centCube.slice(h) = BhatX0cent;
    
    // tHat and Tstat
    arma::vec temp_tHat = square(X0projCent.col(j)) / gaHat[j];
    tHat += temp_tHat;
    tHatMat.col(h) = tHat;
    
    arma::mat Tstat = BhatX0cent;
    Tstat.each_col() /= sqrt(tHat);
    Tstat *= std::sqrt(n);
    TstatArr.slice(h) = Tstat;
  }
  
  return List::create(
    Named("YhatArr") = YhatArr,
    Named("erHatArr") = erHatArr,
    Named("tHatMat") = tHatMat,
    Named("TstatArr") = TstatArr,
    Named("BhatX0cent") = BhatX0centCube
  );
}

arma::cube fofrMRtest_short(
    const arma::mat& Y,
    const arma::mat& X,
    const arma::mat& X0projCent,
    const arma::vec& gaHat,
    const arma::mat& phiHat,
    const arma::mat& tHatMat,
    const arma::cube& BhatX0cent_cent
) {
  int n = X.n_rows;
  int tt = X.n_cols;
  int n0 = X0projCent.n_rows;
  int h_max = tHatMat.n_cols;
  double scal = 1.0 / tt;
  
  // Compute projection matrix
  arma::mat DtHat = arma::cov(Y, X) * (n - 1.0) / n;
  arma::mat DtHat_phiHat = DtHat * phiHat * scal;
  
  // Initialize output cube
  arma::cube TstatCube(n0, tt, h_max, arma::fill::zeros);
  
  for (int i0 = 0; i0 < n0; ++i0) {
    arma::rowvec BhatX0cent(tt, arma::fill::zeros);
    
    for (int h = 0; h < h_max; ++h) {
      arma::rowvec tempBhat = (X0projCent(i0, h) / gaHat[h]) * DtHat_phiHat.col(h).t();
      BhatX0cent += tempBhat;
      
      // arma::rowvec TstatRow = std::sqrt(n / tHatMat(i0, h)) * BhatX0cent;
      for (int t = 0; t < tt; ++t) {
        // TstatCube(i0, t, h) = TstatRow[t];
        TstatCube(i0, t, h) = std::sqrt(n / tHatMat(i0, h)) *(
          BhatX0cent[t] - BhatX0cent_cent(i0, t, h)
        );
      }
    }
  }
  
  return TstatCube;
}

// [[Rcpp::export]]
Rcpp::List fofrMRtestBTS(
    const arma::mat& erHat,              // n × tt residual matrix
    const arma::mat& Yhat,               // n × tt fitted values
    const arma::mat& X,                  // n × tt predictors
    const arma::mat& phiHat,             // tt × h_max
    const arma::vec& gaHat,              // h_max × 1
    const arma::mat& X0projCent,         // n0 × h_max
    const arma::mat& tHatMat,            // n0 × h_max
    const arma::cube& BhatX0cent_cent,   // n0 x tt x h_max (for bootstrap centering)
    int Mbts
) {
  int n = X.n_rows;
  int tt = X.n_cols;
  int n0 = X0projCent.n_rows;
  int h_max = tHatMat.n_cols;
  
  Rcpp::List TstatCubeList(Mbts);
  
  for (int iBTS = 0; iBTS < Mbts; ++iBTS) {
    // Bootstrap indices
    arma::uvec idx = arma::conv_to<arma::uvec>::from(n * arma::randu<arma::vec>(n));
    arma::mat erStar = erHat.rows(idx);
    arma::mat YStar = Yhat + erStar;
    //  n0 × tt × h_max
    arma::cube TstatCube = fofrMRtest_short(
      YStar, X, phiHat, gaHat, X0projCent, tHatMat, BhatX0cent_cent
    );
    TstatCubeList[iBTS] = TstatCube;
  }
  
  return TstatCubeList;
}


// [[Rcpp::export]]
Rcpp::List AtoB(Rcpp::List A) {
  int Mbts = A.size();
  if (Mbts == 0) {
    Rcpp::stop("Input list A is empty.");
  }
  
  arma::cube firstCube = Rcpp::as<arma::cube>(A[0]);
  int n0 = firstCube.n_rows;
  int tt = firstCube.n_cols;
  int h_max = firstCube.n_slices;
  
  std::vector<arma::cube> B_cubes(n0);  // store actual cubes
  
  for (int i0 = 0; i0 < n0; ++i0) {
    B_cubes[i0] = arma::cube(tt, h_max, Mbts, arma::fill::none);
  }
  
  for (int iBTS = 0; iBTS < Mbts; ++iBTS) {
    arma::cube cubeA = Rcpp::as<arma::cube>(A[iBTS]);
    for (int i0 = 0; i0 < n0; ++i0) {
      B_cubes[i0].slice(iBTS) = cubeA.row(i0);
    }
  }
  
  Rcpp::List B(n0);
  for (int i0 = 0; i0 < n0; ++i0) {
    B[i0] = B_cubes[i0];  // assign once
  }
  
  return B;
}




