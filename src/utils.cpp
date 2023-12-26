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
  Eigen::MatrixXd Av = A.transpose() * v;
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
Eigen::MatrixXd addm_rcpp(Eigen::MatrixXd A,
    Eigen::MatrixXd B) {
  // vector-vector addition
  Eigen::MatrixXd C = A+B;
  // return
  return(C);
}

// [[Rcpp::export]]
Eigen::VectorXd addv_rcpp(Eigen::VectorXd u,
                         Eigen::VectorXd v) {
  // vector-vector addition
  Eigen::VectorXd w = u+v;
  // return
  return(w);
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


