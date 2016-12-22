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


// [[Rcpp::export]]
arma::mat cpp_test(arma::vec par, arma::mat X, arma::mat U, arma::vec par_fixed) {
    
    arma::uword n = X.n_rows;
    arma::uword Tmax = X.n_cols;
    
    arma::vec par_full = par_fixed;
    arma::uvec par_nas = arma::find_nonfinite(par_fixed);
    par_full.elem(par_nas) = par;
    
    // =============
    // set up coefficient matrices
    // =============
    
    arma::mat B0;
    B0.insert_cols(0, par_full.subvec(0,n-1));
    B0.reshape(n, 1);
    
    arma::mat B = arma::zeros(n, n);
    arma::uword j = n;
    arma::mat row;
    for(arma::uword i=0; i<n; ++i) {
        row = par_full.subvec(j, (j + n - 1));
        row.reshape(1,n);
        B.row(i) = row;
        j += n;
    }
    
    arma::uword N;
    arma::sword n2 = std::pow(n,2);
    
    arma::vec Se_vec = par_full.subvec((n+n2),(n+n2+n-1));
    N = Se_vec.n_elem;
    for(arma::uword i=0; i<N; ++i) {
        Se_vec(i) = std::pow(Se_vec(i), 2);
    }
    arma::mat Se = arma::diagmat(Se_vec);
    
    arma::vec Su_vec = par_full.subvec((n+n2+n),(n+n2+n+n-1));
    N = Su_vec.n_elem;
    for(arma::uword i=0; i<N; ++i) {
        Su_vec(i) = std::pow(Su_vec(i), 2);
    }
    arma::mat Su = arma::diagmat(Su_vec);
    
    arma::vec Sb_vec = par_full.subvec((n+n2+n+n), (n+n2+n+n+n*(n+1))-1);
    N = Sb_vec.n_elem;
    for(arma::uword i=0; i<N; ++i) {
        Sb_vec(i) = std::pow(Sb_vec(i), 2);
    }
    arma::mat Sb = arma::diagmat(Sb_vec);
    
    
    // =============
    // set up independent variable
    // =============
    arma::uword nu;
    arma::mat C;
    
    if(U.n_cols > 1) {
        nu = U.n_rows;
        C = arma::zeros(n, nu);
        j = (n+n2+n+n+n*(n+1));
        for(arma::uword i=0; i<n; ++i) {
            row = par_full.subvec(j, (j + nu - 1));
            row.reshape(1,nu);
            C.row(i) = row;
            j += nu;
        }
    }
    
    arma::mat S;
    arma::mat S_temp = arma::zeros(n, n*(n+1));
    S = arma::join_rows(Se, S_temp);
    S_temp = arma::zeros(n*(n+1), n);
    S_temp = arma::join_rows(S_temp, Sb);
    S = arma::join_cols(S, S_temp);
    
    arma::mat Z;
    Z = arma::join_rows(arma::eye<arma::mat>(n,n), arma::zeros(n, n*(n+1)));
    
    // =============
    // Initial unconditional values
    // =============
    arma::vec x = X.col(0);
    
    
    
    // If the initial parameter values imply a stationary distribution, then the initial
    // covariance matrix of the process error, PP, is computed from the coefficients. If 
    // not, the initial value of PP given by the covariance matrix of the process error
    // variation (effectively assuming that the dominant eigenvalue of the system is zero).
    
    // double std::abs(y)

    arma::cx_vec eigval = arma::eig_gen(B);
    
    arma::vec eig_reals = arma::zeros(eigval.n_elem);
    double eig_real_i;
    for(arma::uword i=0; i<(eigval.n_elem); ++i){
        eig_real_i = eigval(i).real();
        eig_reals(i) = std::abs(eig_real_i);
    }
    
    arma::mat PP;
    arma::mat PP_k;
    arma::mat PP_Se = Se;
    arma::mat PP_dm;
    if(max(eig_reals) < 1){
        PP_k = cpp_kron(B,B);
        PP_Se.reshape(n*n,1);
        PP_dm = arma::eye<arma::mat>(n2,n2);
        PP = cpp_solve(PP_dm - PP_k);
        PP = PP * PP_Se;
        PP.reshape(n, n);
    } else {
        PP = Se;
    }
    
    PP = join_rows(PP, arma::zeros(n, n*(n+1)));
    arma::mat PP_Sb = join_rows(arma::zeros(n*(n+1), n), Sb);
    PP = join_cols(PP, PP_Sb);
    
    double logFt = 0;
    double vFv = 0;
    
    
    
//     for(t in 2:Tmax) {
//         
// # PREDICTION EQUATIONS
//         
//         B12 <- diag(n) - B
//             B13 <- kronecker(t(x-B0),diag(n))
//             
//             BB <- as.matrix(rbind(cbind(B, B12, B13), 
//                                   cbind(matrix(0, n, n), diag(n), matrix(0, n, n^2)), 
//                                   cbind(matrix(0, n^2, 2*n), diag(n^2))))
//             PP <- BB %*% PP %*% t(BB) + S
//             
//             if (is.null(U)){
//                 x <- B0 + B %*% (x-B0)
//             } else {				
//                 if (nu == 1) {
//                     x <- B0 + B %*% (x-B0) + C * U[t]
//                 } else {
//                     x <- B0 + B %*% (x-B0) + C %*% U[,t]
//                 }
//             }
//             
// # UPDATING EQUATIONS
//             if(!any(is.na(X[,t]))){
//                 
//                 FF <- Z %*% PP %*% t(Z) + Su
//                 invF <- solve(FF)
//                 
//                 y <- matrix(c(x, B0, t(B)), ncol=1)		
//                 v <- X[,t] - Z %*% y
//                 
//                 y <- y + PP %*% t(Z) %*% invF %*% v
//                 PP <- PP - PP %*% t(Z) %*% invF %*% Z %*% PP
//                 
//                 x <- y[1:n]
//                 B0 <- y[(n+1):(2*n)]
//                 B <- matrix(y[(2*n+1):length(y)], nrow=n, ncol=n, byrow = TRUE)
//                 
// # TERMS OF LIKELIHOOD FUNCTION
//                 
// # "determinant(FF)$modulus[1]" gives the log of the determinant by default
//                 logdetFF <- determinant(FF)$modulus[1]
//                 logFt <- logFt + logdetFF
//                     
//                     vFv <- vFv + t(v) %*% invF %*% v
//             }
//     }
    
    return(PP);
}

