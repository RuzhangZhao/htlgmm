#include <Rcpp.h>
#include <RcppEigen.h>
using namespace Rcpp;
using namespace RcppEigen;
// [[Rcpp::depends(RcppEigen)]]

// [[Rcpp::export]]
Eigen::MatrixXd crossprod_rcpp(
    Eigen::MatrixXd A, Eigen::MatrixXd B) {
  // matrix-matrix multiplication
  Eigen::MatrixXd AB = A.transpose() * B;
  // return
  return(AB);
}

// [[Rcpp::export]]
Eigen::MatrixXd prod_rcpp(
    Eigen::MatrixXd A, Eigen::MatrixXd B) {
  // matrix-matrix multiplication
  Eigen::MatrixXd AB = A * B;
  // return
  return(AB);
}

// [[Rcpp::export]]
Eigen::VectorXd crossprodv_rcpp(
    Eigen::MatrixXd A, Eigen::VectorXd v) {
  // matrix-vector multiplication
  Eigen::VectorXd Av = A.transpose() * v;
  // return
  return(Av);
}

// [[Rcpp::export]]
Eigen::VectorXd prodv_rcpp(
    Eigen::MatrixXd A, Eigen::VectorXd v) {
  // matrix-vector multiplication
  Eigen::VectorXd Av = A * v;
  // return
  return(Av);
}

// [[Rcpp::export]]
Eigen::MatrixXd self_crossprod_rcpp(
    Eigen::MatrixXd A) {
  // matrix-matrix multiplication
  Eigen::MatrixXd AA = crossprod_rcpp(A,A);
  // return
  return(AA);
}

// [[Rcpp::export]]
Eigen::VectorXd expit_rcpp(
    Eigen::VectorXd v) {
  // expit(s) = exp(s)/(1+exp(s))
  Eigen::VectorXd negv = -v;
  Eigen::VectorXd expitv = 1.0/(negv.array().exp()+1);
  // return
  return(expitv);
}

// [[Rcpp::export]]
Eigen::VectorXd dexpit_rcpp(
    Eigen::VectorXd v) {
  // dexpit(s) = expit(s)*(1-expit(s))
  Eigen::VectorXd expitv = expit_rcpp(v);
  Eigen::VectorXd dexpitv = expitv.array()*(1.0-expitv.array());
  // return
  return(dexpitv);
}


// [[Rcpp::export]]
Eigen::VectorXd timesv_rcpp(Eigen::VectorXd u,
                         Eigen::VectorXd v) {
  // matrix-matrix multiplication
  Eigen::VectorXd w = u.array()*v.array();
  // return
  return(w);
}

// [[Rcpp::export]]
Eigen::MatrixXd timesm_rcpp(Eigen::MatrixXd A,
                           Eigen::MatrixXd B) {
  // matrix-matrix multiplication
  Eigen::MatrixXd AB = A.array()*B.array();
  // return
  return(AB);
}


// [[Rcpp::export]]
Eigen::MatrixXd square_rcpp(const Eigen::MatrixXd mat) {
  // Calculate the square of each element
  Eigen::MatrixXd mat2 = mat.array().square();;
  return mat2;
}


// [[Rcpp::export]]
Eigen::MatrixXd expitm_rcpp(
    Eigen::MatrixXd A) {
  // expit(s) = exp(s)/(1+exp(s))
  Eigen::MatrixXd negA = -A;
  Eigen::MatrixXd expitA = 1.0/(negA.array().exp()+1);
  // return
  return(expitA);
}

// [[Rcpp::export]]
Eigen::MatrixXd sqrtchoinv_rcpp(const Eigen::MatrixXd& matrix) {
  // Check if the matrix is square
  if (matrix.rows() != matrix.cols()) {
    Rcpp::stop("Matrix must be square.");
  }

  // Compute the Cholesky decomposition
  Eigen::LLT<Eigen::MatrixXd> llt(matrix);
  if (llt.info() != Eigen::Success) {
    Rcpp::stop("Cholesky decomposition failed.");
  }
  // Compute the inverse of L
  Eigen::MatrixXd L_inv = llt.matrixL().solve(Eigen::MatrixXd::Identity(matrix.rows(), matrix.cols()));

  // The square root of the inverse of A is L^{-T}
  return L_inv;
}

// [[Rcpp::export]]
Eigen::MatrixXd choinv_rcpp(const Eigen::MatrixXd& matrix) {
  // Ensure the matrix is square
  if (matrix.rows() != matrix.cols()) {
    Rcpp::stop("Matrix must be square.");
  }
  // Compute the Cholesky decomposition (LLT)
  Eigen::LLT<Eigen::MatrixXd> llt(matrix);

  // Check if the matrix is positive definite
  if(llt.info() != Eigen::Success) {
    Rcpp::stop("Matrix is not positive definite.");
  }

  // Compute the inverse based on the Cholesky decomposition
  return llt.solve(Eigen::MatrixXd::Identity(matrix.rows(), matrix.cols()));
}


// [[Rcpp::export]]
Eigen::MatrixXd sqrtcho_rcpp(const Eigen::MatrixXd& matrix) {
  // Check if the matrix is square
  if (matrix.rows() != matrix.cols()) {
    Rcpp::stop("Matrix must be square.");
  }

  // Compute the Cholesky decomposition
  Eigen::LLT<Eigen::MatrixXd> llt(matrix);
  if (llt.info() != Eigen::Success) {
    Rcpp::stop("Cholesky decomposition failed.");
  }
  // Compute the inverse of L
  Eigen::MatrixXd L_inv = llt.matrixL().transpose();

  // The square root of A is L^{T}
  return L_inv;
}

