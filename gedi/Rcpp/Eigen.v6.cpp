// [[Rcpp::depends(RcppEigen)]]

#include <Rcpp.h>
#include <RcppEigen.h>

// Returns A %*% B in the form of a matrix
// [[Rcpp::export]]
SEXP eigenMatProduct(Eigen::MatrixXd A, Eigen::MatrixXd B){
  return Rcpp::wrap(A * B);
}

// Returns A %*% t(B) in the form of a matrix
// [[Rcpp::export]]
SEXP eigenMatTcrossprod(Eigen::MatrixXd A, Eigen::MatrixXd B){
  return Rcpp::wrap(A * B.transpose());
}

// Returns t(A) %*% B in the form of a matrix
// [[Rcpp::export]]
SEXP eigenMatCrossprod(Eigen::MatrixXd A, Eigen::MatrixXd B){
  return Rcpp::wrap(A.transpose() * B);
}

// Returns A %*% b in the form of a column vector
// [[Rcpp::export]]
SEXP eigenMatVecProduct(Eigen::MatrixXd A, Eigen::VectorXd b){
  return Rcpp::wrap(A * b);
}

// Returns the transpose of t(a) %*% B in the form of a column vector
// [[Rcpp::export]]
SEXP eigenVecMatProduct(Eigen::VectorXd a, Eigen::MatrixXd B){
  return Rcpp::wrap(B.transpose() * a);
}

// Returns the transpose of a %*% t(b) in the form of a matrix
// [[Rcpp::export]]
SEXP eigenVecVecProduct(Eigen::VectorXd a, Eigen::VectorXd b){
  return Rcpp::wrap(a * b.transpose());
}

// Returns the row-wise L2 norm of a matrix
// [[Rcpp::export]]
SEXP eigenRowL2(Eigen::MatrixXd A){
  return Rcpp::wrap(A.rowwise().squaredNorm());
}

// Returns the column-wise L2 norm of a matrix
// [[Rcpp::export]]
SEXP eigenColL2(Eigen::MatrixXd A){
  return Rcpp::wrap(A.colwise().squaredNorm());
}

// Returns the column-wise L1 norm of a matrix
// [[Rcpp::export]]
SEXP eigenColL1(Eigen::MatrixXd A){
  return Rcpp::wrap(A.colwise().lpNorm<1>());
}

// Multiplies each row of a matrix by a vector (coefficient-wise)
// [[Rcpp::export]]
SEXP eigenRowMult( Eigen::MatrixXd A, Eigen::ArrayXd b ){
  return Rcpp::wrap( (A.array().rowwise()*b.transpose()).matrix() );
}

/*******************************************************/
// Returns the transpose of a %*% t(b) in the form of a matrix
// [[Rcpp::export]]
SEXP solveBi(
    Eigen::MatrixXd Yi,
    Eigen::VectorXd D,
    Eigen::MatrixXd Z, Eigen::MatrixXd Qi,
    Eigen::MatrixXd diagK,
    Eigen::VectorXd si, Eigen::VectorXd o, Eigen::VectorXd oi ) {
  
  Eigen::MatrixXd Wi = (Z+Qi)*D.asDiagonal();
  return Rcpp::wrap(
    (Wi.transpose()*Wi).ldlt().solve(diagK) *
      ( Wi.transpose() * ((Yi.rowwise()-si.transpose()).colwise()-(o+oi)) )
    );
  
  /* The following are other alternatives that could be used, although they are slightly slower*/
  //return Rcpp::wrap( WiWi.colPivHouseholderQr().solve(diagK) * (Wi.transpose()*Yi_res) );
  //return Rcpp::wrap( (Wi.transpose()*Wi).ldlt().solve(Wi.transpose()*Yi_res) );
}



// Returns the residual of Yi after taking into account everything other than ZDBi
// [[Rcpp::export]]
SEXP Yi_resZ(
    Eigen::MatrixXd Yi,
    Eigen::MatrixXd QiDBi,
    Eigen::VectorXd si, Eigen::VectorXd o, Eigen::VectorXd oi ) {
  
  return Rcpp::wrap( ((Yi-QiDBi).rowwise()-si.transpose()).colwise()-(o+oi) );
}


// Solve Z without C
// [[Rcpp::export]]
SEXP solveZ_noC(
    Eigen::MatrixXd Y_res,
    Eigen::VectorXd D,
    Eigen::MatrixXd B,
    Eigen::MatrixXd diagK,
    double lambda ) {
  
  Eigen::MatrixXd DB = D.asDiagonal()*B;
  return Rcpp::wrap( (Y_res*DB.transpose()) * (DB*DB.transpose()+lambda*diagK).ldlt().solve(diagK) );
}


// Solve Z without C
// [[Rcpp::export]]
SEXP solveZ_wC(
    Eigen::MatrixXd Y_res,
    Eigen::VectorXd D,
    Eigen::MatrixXd B,
    Eigen::MatrixXd C, Eigen::MatrixXd A,
    Eigen::MatrixXd diagK,
    double lambda ) {
  
  Eigen::MatrixXd DB = D.asDiagonal()*B;
  return Rcpp::wrap( (Y_res*DB.transpose()+lambda*C*A) * (DB*DB.transpose()+lambda*diagK).ldlt().solve(diagK) );
}


// Solve U so as to minimize ||Y-USB||, with you being orthonormal
// [[Rcpp::export]]
SEXP solveU(
    Eigen::MatrixXd Y_res,
    Eigen::VectorXd S,
    Eigen::MatrixXd B ) {
  
  Eigen::BDCSVD<Eigen::MatrixXd> svd(
      Y_res * (S.asDiagonal()*B).transpose(),
      Eigen::ComputeThinU | Eigen::ComputeThinV );
  return Rcpp::wrap( svd.matrixU()*svd.matrixV().transpose() );
}


// Solve Qi without H
// [[Rcpp::export]]
SEXP solveQi_noH(
    Eigen::MatrixXd Yi,
    Eigen::MatrixXd ZDBi,
    Eigen::VectorXd si, Eigen::VectorXd o, Eigen::VectorXd oi,
    Eigen::VectorXd D,
    Eigen::MatrixXd Bi,
    Eigen::MatrixXd diagK,
    double lambda ) {

  Eigen::MatrixXd DBi = D.asDiagonal()*Bi;
  return Rcpp::wrap(
    ( (((Yi-ZDBi).rowwise()-si.transpose()).colwise()-(o+oi)) * DBi.transpose() ) *
      (DBi*DBi.transpose()+lambda*diagK).ldlt().solve(diagK) );
}


// Solve Qi with H
// [[Rcpp::export]]
SEXP solveQi_wH(
    Eigen::MatrixXd Yi,
    Eigen::MatrixXd ZDBi,
    Eigen::VectorXd si, Eigen::VectorXd o, Eigen::VectorXd oi,
    Eigen::VectorXd D,
    Eigen::MatrixXd Bi,
    Eigen::MatrixXd H, Eigen::MatrixXd Ri,
    Eigen::MatrixXd diagK,
    double lambda ) {
  
  Eigen::MatrixXd DBi = D.asDiagonal()*Bi;
  return Rcpp::wrap(
    ( (((Yi-ZDBi).rowwise()-si.transpose()).colwise()-(o+oi)) * DBi.transpose() +
      lambda*H*Ri ) *
      (DBi*DBi.transpose()+lambda*diagK).ldlt().solve(diagK) );
}


// Solve oi
// [[Rcpp::export]]
SEXP solveOi(
    Eigen::MatrixXd Yi,
    Eigen::MatrixXd ZDBi, Eigen::MatrixXd QiDBi,
    Eigen::VectorXd ni,
    Eigen::VectorXd si, Eigen::VectorXd o,
    double Ni, double lambda ) {
  
  return Rcpp::wrap(
    (((Yi-ZDBi-QiDBi).rowwise()-si.transpose()).colwise()-o) * ni / (Ni+lambda) );
}


// Solve si
// [[Rcpp::export]]
SEXP solveSi(
    Eigen::MatrixXd Yi,
    Eigen::MatrixXd ZDBi, Eigen::MatrixXd QiDBi,
    Eigen::VectorXd j,
    Eigen::VectorXd o, Eigen::VectorXd oi,
    Eigen::VectorXd si_0,
    double J, double lambda ) {
  
  return Rcpp::wrap(
    ( ((Yi-ZDBi-QiDBi).colwise()-(o+oi)).transpose()*j + lambda*si_0 ) / (J+lambda) );
}


// Returns the row sum of the residual of Yi after taking into account everything other than o
// [[Rcpp::export]]
SEXP Yi_resO_rowSum(
    Eigen::MatrixXd Yi,
    Eigen::MatrixXd ZDBi, Eigen::MatrixXd QiDBi,
    Eigen::VectorXd ni,
    Eigen::VectorXd si, Eigen::VectorXd oi ) {
  
  return Rcpp::wrap( (((Yi-ZDBi-QiDBi).rowwise()-si.transpose()).colwise()-oi)*ni );
}



// Solves Yi as log of unobserved expression
// [[Rcpp::export]]
SEXP solveYi(
    Eigen::SparseMatrix<double> Mi,
    Eigen::ArrayXXd Yi,
    Eigen::ArrayXXd ZDBi, Eigen::ArrayXXd QiDBi,
    Eigen::ArrayXd si, Eigen::ArrayXd o, Eigen::ArrayXd oi,
    double sigma2 ) {
  
  // Calculate the estimate of Yi based on model parameters
  Eigen::ArrayXXd Yi_hat = ((ZDBi+QiDBi).rowwise()+si.transpose()).colwise()+(o+oi);
  // Use Halley's method to obtain an approximation of Yi.
  //   Note that with a reasonable starting point, Halley's method
  //   leads to reasonable estimates in one iteration, which should improve
  //   in subsequent iterations
  Eigen::ArrayXXd fpp = sigma2 * Yi.exp();
  Eigen::ArrayXXd fp = fpp + 1;
  Eigen::ArrayXXd f = fpp + Yi - Yi_hat - Eigen::ArrayXXd(sigma2*Mi);
  
  return Rcpp::wrap( Yi - 2*f*fp/(2*fp.square()-f*fpp) );

  // TODO: This function can become numerically unstable if the current
  // estimate of Yi is too large for some entries, and therefore exp(Yi)
  // becomes unfeasible to compute. This should not happen in principle, since
  // Yi is meant to stay close to values that lead to Mi~exp(Yi). But in case
  // this happens, we can use an alternative solution for large Yi values, in
  // which f, fp, and fpp are all multiplied by exp(-Yi), similar to what is used in the
  // solveYi_paired function below.
  
}


// Solves Yi as log-ratio of the unobserved paired abundances
// [[Rcpp::export]]
SEXP solveYi_paired(
    Eigen::SparseMatrix<double> M1i, Eigen::SparseMatrix<double> M2i,
    Eigen::ArrayXXd Yi,
    Eigen::ArrayXXd ZDBi, Eigen::ArrayXXd QiDBi,
    Eigen::ArrayXd si, Eigen::ArrayXd o, Eigen::ArrayXd oi,
    double sigma2 ) {
  
  // Calculate the estimate of Yi based on model parameters
  Eigen::ArrayXXd Yi_hat = ((ZDBi+QiDBi).rowwise()+si.transpose()).colwise()+(o+oi);
  
  // Use Halley's method to obtain an approximation of Yi.
  //   Note that with a reasonable starting point, Halley's method
  //   leads to reasonable estimates in one iteration, which should improve
  //   in subsequent iterations
  Eigen::ArrayXXd solution = Yi; // start with the current Yi, and obtain a better estimate via Halley's method
  Eigen::ArrayXXd alpha = Yi_hat - Eigen::ArrayXXd(sigma2*M2i);
  Eigen::ArrayXXd beta = Yi_hat + Eigen::ArrayXXd(sigma2*M1i);
  
  double exp_Yi, f, fp, fpp;
  for( int i = 0; i < solution.size(); i ++ ) {
    if( Yi(i) > 0 ) {
      exp_Yi = exp(-Yi(i));
      f = Yi(i)-alpha(i)+(Yi(i)-beta(i))*exp_Yi;
      fp = 1+exp_Yi*(beta(i)-Yi(i)+1);
      fpp = exp_Yi*(Yi(i)-beta(i)-2);
    } else {
      exp_Yi = exp(Yi(i));
      f = exp_Yi*(Yi(i)-alpha(i))+Yi(i)-beta(i);
      fp = exp_Yi+beta(i)-Yi(i)+1;
      fpp = Yi(i)-beta(i)-2;
    }
    solution(i) -= 2*f*fp/(2*fp*fp-f*fpp);
    // check for potential overshoot; if so, set the solution to the respective limit
    if( solution(i) < alpha(i) ) {
      solution(i) = alpha(i);
    } else if( solution(i) > beta(i) ) {
      solution(i) = beta(i);
    }
  }
  
  return Rcpp::wrap( solution );
}

// Below is an older version which could have resulted in numerical instability for
// large negative values of the current estimate of Yi.

// // Solves Yi as log-ratio of the unobserved paired abundances
// // [[Rcpp::export]]
// SEXP solveYi_paired(
//     Eigen::SparseMatrix<double> M1i, Eigen::SparseMatrix<double> M2i,
//     Eigen::ArrayXXd Yi,
//     Eigen::ArrayXXd ZDBi, Eigen::ArrayXXd QiDBi,
//     Eigen::ArrayXd si, Eigen::ArrayXd o, Eigen::ArrayXd oi,
//     double sigma2 ) {
//   
//   // Calculate the estimate of Yi based on model parameters
//   Eigen::ArrayXXd Yi_hat = ((ZDBi+QiDBi).rowwise()+si.transpose()).colwise()+(o+oi);
//   // Use Halley's method to obtain an approximation of Yi.
//   //   Note that with a reasonable starting point, Halley's method
//   //   leads to reasonable estimates in one iteration, which should improve
//   //   in subsequent iterations
//   Eigen::ArrayXXd inverse_exp_Yi = (-Yi).exp();
//   Eigen::ArrayXXd alpha = Yi_hat - Eigen::ArrayXXd(sigma2*M2i);
//   Eigen::ArrayXXd beta = Yi_hat + Eigen::ArrayXXd(sigma2*M1i);
// 
//   Eigen::ArrayXXd f = -alpha + Yi + (Yi-beta) * inverse_exp_Yi;
//   Eigen::ArrayXXd fp = 1 + inverse_exp_Yi * (beta-Yi+1);
//   Eigen::ArrayXXd fpp = inverse_exp_Yi * (Yi-beta-2);
//   
//   Eigen::ArrayXXd solution = Yi - 2*f*fp/(2*fp.square()-f*fpp);
//   
//   // check for potential overshoot; if so, set the solution to the respective limit
//   for( int i = 0; i < solution.size(); i ++ ) {
//     if( solution(i) < alpha(i) ) {
//       solution(i) = alpha(i);
//     } else if( solution(i) > beta(i) ) {
//       solution(i) = beta(i);
//     }
//   }
//   
//   return Rcpp::wrap( solution );
// }


// Returns the L2 norm of a vector
// [[Rcpp::export]]
double vecL2_noPrior( Eigen::VectorXd x ) {
  return x.squaredNorm();
}
// Returns the L2 norm of a vector after subtracting a prior (squared Euclidean distance)
// [[Rcpp::export]]
double vecL2_wPrior( Eigen::VectorXd x, Eigen::VectorXd prior ) {
  return (x-prior).squaredNorm();
}
// Returns the L2 norm of a matrix
// [[Rcpp::export]]
double matL2_noPrior( Eigen::MatrixXd X ) {
  return X.squaredNorm();
}

// Returns the L2 norm of a matrix after subtracting a prior 
// [[Rcpp::export]]
double matL2_wPrior( Eigen::MatrixXd X, Eigen::MatrixXd A, Eigen::MatrixXd B ) {
  return (X-A*B).squaredNorm();
}

// // Returns the SSE (sum-squared-error) of Yi
// // [[Rcpp::export]]
// double Yi_SSE(
//     Eigen::MatrixXd Yi,
//     Eigen::MatrixXd ZDBi, Eigen::MatrixXd QiDBi,
//     Eigen::VectorXd si, Eigen::VectorXd o, Eigen::VectorXd oi ) {
//   
//   return ( ((Yi-ZDBi-QiDBi).rowwise()-si.transpose()).colwise()-(o+oi) ).squaredNorm();
// }

// Returns the variance of Yi-Yhat
// [[Rcpp::export]]
double Yi_SSE(
    Eigen::MatrixXd Yi,
    Eigen::MatrixXd ZDBi, Eigen::MatrixXd QiDBi,
    Eigen::VectorXd si, Eigen::VectorXd o, Eigen::VectorXd oi,
    double sigma2 ) {
  
  return ( ((Yi-ZDBi-QiDBi).rowwise()-si.transpose()).colwise()-(o+oi) ).squaredNorm() +
    (1/(Eigen::ArrayXXd(Yi).exp()+1/sigma2)).sum();
}


// Returns the variance of Yi-Yhat
// [[Rcpp::export]]
double Yi_SSE_paired(
    Eigen::MatrixXd Yi,
    Eigen::SparseMatrix<double> M1i, Eigen::SparseMatrix<double> M2i,
    Eigen::MatrixXd ZDBi, Eigen::MatrixXd QiDBi,
    Eigen::VectorXd si, Eigen::VectorXd o, Eigen::VectorXd oi,
    double sigma2 ) {

  Eigen::ArrayXXd expYi = (-Eigen::ArrayXXd(Yi).abs()).exp();
  Eigen::ArrayXXd M = Eigen::ArrayXXd(M1i+M2i);
  

  return ( ((Yi-ZDBi-QiDBi).rowwise()-si.transpose()).colwise()-(o+oi) ).squaredNorm() +
    (1/(M*expYi/(1+expYi).square()+1/sigma2)).sum();
}


// Returns the root mean squared difference of two matrices
// [[Rcpp::export]]
double matRMSD( Eigen::MatrixXd A, Eigen::MatrixXd B ) {
  return sqrt( (A-B).squaredNorm()/A.size() );
}
// Returns the root mean squared difference of two vectors
// [[Rcpp::export]]
double vecRMSD( Eigen::VectorXd A, Eigen::VectorXd B ) {
  return sqrt( (A-B).squaredNorm()/A.size() );
}


// Performs cbind of two matrices
// Benchmarking shows that R cbind is twice as fast. So, let's continue using that.
// SEXP eigenCBind( Eigen::MatrixXd A, Eigen::MatrixXd B ) {
//   
//   Eigen::MatrixXd AB(A.rows(),A.cols()+B.cols()); AB << A,B;
//   return Rcpp::wrap( AB );
// }

// Performs cbind of two matrices
// Benchmarking shows that R cbind is still slightly faster. So, let's continue using that.
// Rcpp::NumericMatrix RcppCBind(Rcpp::NumericMatrix a, Rcpp::NumericMatrix b) {
//   int acoln = a.ncol();
//   int bcoln = b.ncol();
//   Rcpp::NumericMatrix out = Rcpp::no_init_matrix(a.nrow(), acoln + bcoln);
//   for (int j = 0; j < acoln + bcoln; j++) {
//     if (j < acoln) {
//       out(Rcpp::_, j) = a(Rcpp::_, j);
//     } else {
//       out(Rcpp::_, j) = b(Rcpp::_, j - acoln);
//     }
//   }
//   return out;
// }
