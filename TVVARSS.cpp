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

// [[Rcpp::export]]
arma::mat cpp_mmult(arma::mat x, arma::mat y) {
    arma::mat z = x * y;
    return(z);
}



arma::mat cpp_TVVARSS_ml(arma::vec par, arma::mat X, arma::mat U, arma::vec par_fixed) {
    // Note: X and U are transposed, so time runs through columns
    int n = X.n_rows;
    int Tmax = X.n_cols;
    
    arma::vec par_full = par_fixed;
    arma::uvec par_nas = arma::find_nonfinite(par_fixed);
    par_full.elem(par_nas) = par;
    
    // set up coefficient matrices
    
//         B0 <- matrix(par.full[1:n], nrow=n, ncol=1)
//             B <- matrix(par.full[(n+1):(n+n^2)], nrow=n, ncol=n, byrow = TRUE)
//             Se <- diag(par.full[(n+n^2+1):(n+n^2+n)]^2)
//             Su <- diag(par.full[(n+n^2+n+1):(n+n^2+n+n)]^2)	
//             Sb <- diag(par.full[(n+n^2+n+n+1):(n+n^2+n+n+n*(n+1))]^2)
//             
// # set up independent variable
//             if(!is.null(U)){
//                 nu <- dim(U)[1]
//                 C <- matrix(par.full[(n+n^2+n+n+n*(n+1)+1):(n+n^2+n+n+n*(n+1)+nu*n)], 
//                             nrow=n, ncol=nu, byrow = TRUE)
//             }		
//             
//             S <- as.matrix(rbind(cbind(Se, matrix(0, n, n*(n+1))), 
//                                  cbind(matrix(0, n*(n+1), n), Sb)))		
//                 Z <- as.matrix(cbind(diag(n), matrix(0, n, n*(n+1))))
//                 
// # Initial unconditional values
//                 x <- X[,1]
//                 
// # If the initial parameter values imply a stationary distribution, then the initial
// # covariance matrix of the process error, PP, is computed from the coefficients. If 
// # not, the initial value of PP given by the covariance matrix of the process error
// # variation (effectively assuming that the dominant eigenvalue of the system is zero).
//                 
//                 if (max(abs(eigen(B)$values)) < 1) {
//                     PP <- solve(diag(n*n)-kronecker(B,B)) %*% matrix(Se, nrow = n*n)		
//                     PP <- matrix(PP, n, n)
//                 } else {
//                     PP <- Se
//                 }
//                 PP <- as.matrix(rbind(cbind(PP, matrix(0, n, n*(n+1))), 
//                                       cbind(matrix(0, n*(n+1), n), Sb)))
//                     
//                     logFt <- 0
//                 vFv <- 0
//                 
//                 for(t in 2:Tmax) {
//                     
// # PREDICTION EQUATIONS
//                     
//                     B12 <- diag(n) - B
//                         B13 <- kronecker(t(x-B0),diag(n))
//                         
//                         BB <- as.matrix(rbind(cbind(B, B12, B13), 
//                                               cbind(matrix(0, n, n), diag(n), matrix(0, n, n^2)), 
//                                               cbind(matrix(0, n^2, 2*n), diag(n^2))))
//                         PP <- BB %*% PP %*% t(BB) + S
//                         
//                         if (is.null(U)){
//                             x <- B0 + B %*% (x-B0)
//                         } else {				
//                             if (nu == 1) {
//                                 x <- B0 + B %*% (x-B0) + C * U[t]
//                             } else {
//                                 x <- B0 + B %*% (x-B0) + C %*% U[,t]
//                             }
//                         }
//                         
// # UPDATING EQUATIONS
//                         if(!any(is.na(X[,t]))){
//                             
//                             FF <- Z %*% PP %*% t(Z) + Su
//                             invF <- solve(FF)
//                             
//                             y <- matrix(c(x, B0, t(B)), ncol=1)		
//                             v <- X[,t] - Z %*% y
//                             
//                             y <- y + PP %*% t(Z) %*% invF %*% v
//                             PP <- PP - PP %*% t(Z) %*% invF %*% Z %*% PP
//                             
//                             x <- y[1:n]
//                             B0 <- y[(n+1):(2*n)]
//                             B <- matrix(y[(2*n+1):length(y)], nrow=n, ncol=n, byrow = TRUE)
//                             
// # TERMS OF LIKELIHOOD FUNCTION
//                             
// # "determinant(FF)$modulus[1]" gives the log of the determinant by default
//                             logdetFF <- determinant(FF)$modulus[1]
//                             logFt <- logFt + logdetFF
//                                 
//                                 vFv <- vFv + t(v) %*% invF %*% v
//                         }
//                 }
//                 
//                 LL <- logFt + vFv
//                     if (is.complex(LL)) {
//                         LL <- 10^10
//                     }
// # show(c(LL,par.full[1:6]))
// # browser()
//                     return(LL)
}

