# include <RcppArmadillo.h>

// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;

// This is a simple example of exporting a C++ function to R. You can
// source this function into an R session using the Rcpp::sourceCpp 
// function (or via the Source button on the editor toolbar). Learn
// more about Rcpp at:
//
//   http://www.rcpp.org/
//   http://adv-r.had.co.nz/Rcpp.html
//   http://gallery.rcpp.org/
//

// [[Rcpp::export]]
NumericVector timesTwo(NumericVector x) {
  return x * 2;
}


// You can include R code blocks in C++ files processed with sourceCpp
// (useful for testing and development). The R code will be automatically 
// run after the compilation.
//
// 
// /*** R
// timesTwo(42)
// */


// [[Rcpp::export]]
arma::mat cpp_solve(arma::mat x) {
    arma::mat y = inv(x);
    return(y);
}


// [[Rcpp::export]]
double cpp_log_det(arma::mat x) {
    double y = arma::det(x);
    double z = std::abs(y);
    return(std::log(z));
}


// [[Rcpp::export]]
arma::mat cpp_kron(arma::mat x, arma::mat y) {
    arma::mat z = arma::kron(x, y);
    return(z);
}