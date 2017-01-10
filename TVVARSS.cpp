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


// Same as `determinant(FF)$modulus[1]` in R
// [[Rcpp::export]]
double cpp_log_det(arma::mat x) {
    double y = arma::det(x);
    double z = std::abs(y);
    return(std::log(z));
}


// [[Rcpp::export]]
arma::mat cpp_TVVARSS_ml(arma::vec par, arma::mat X, arma::mat U, arma::vec par_fixed) {
    
    // Note: X and U are transposed, so time runs through columns
    arma::uword n = X.n_rows;
    arma::uword Tmax = X.n_cols;
    
    // Identity matrix of size n by n
    arma::mat n_ident = arma::eye<arma::mat>(n,n);
    // And for n^2 by n^2
    arma::sword n2 = std::pow(n,2);
    arma::mat n2_ident = arma::eye<arma::mat>(n2,n2);
    
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
    arma::mat row_iter;
    for(arma::uword i=0; i<n; ++i) {
        row_iter = par_full.subvec(j, (j + n - 1));
        row_iter.reshape(1,n);
        B.row(i) = row_iter;
        j += n;
    }
    
    arma::uword N;
    
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
    arma::uword nu = U.n_rows;
    arma::mat C;
    
    if(U.n_cols > 1) {
        C = arma::zeros(n, nu);
        j = (n+n2+n+n+n*(n+1));
        for(arma::uword i=0; i<n; ++i) {
            row_iter = par_full.subvec(j, (j + nu - 1));
            row_iter.reshape(1,nu);
            C.row(i) = row_iter;
            j += nu;
        }
    }
    
    arma::mat S;
    arma::mat S_temp = arma::zeros(n, n*(n+1));
    S = arma::join_rows(Se, S_temp);
    S_temp = arma::zeros(n*(n+1), n);
    S_temp = arma::join_rows(S_temp, Sb);
    S = arma::join_cols(S, S_temp);
    
    arma::mat Z = arma::join_rows(n_ident, arma::zeros(n, n*(n+1)));
    
    // =============
    // Initial unconditional values
    // =============
    arma::vec x = X.col(0);
    
    
    
    // If the initial parameter values imply a stationary distribution, then the initial
    // covariance matrix of the process error, PP, is computed from the coefficients. If 
    // not, the initial value of PP given by the covariance matrix of the process error
    // variation (effectively assuming that the dominant eigenvalue of the system is 
    // zero).
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
        PP_k = arma::kron(B,B);
        PP_Se.reshape(n*n,1);
        PP_dm = n2_ident;
        PP = arma::inv(PP_dm - PP_k);
        PP = PP * PP_Se;
        PP.reshape(n, n);
    } else {
        PP = Se;
    }
    
    PP = arma::join_rows(PP, arma::zeros(n, n*(n+1)));
    arma::mat PP_Sb = arma::join_rows(arma::zeros(n*(n+1), n), Sb);
    PP = arma::join_cols(PP, PP_Sb);
    
    double logFt = 0;
    arma::mat vFv = arma::zeros(1,1);
    
    arma::mat BB;
    arma::mat B12;
    arma::mat B13;
    arma::mat BB_temp;
    arma::mat FF;
    arma::mat invF;
    arma::mat y;
    arma::mat v;
    for(arma::uword t=1; t<Tmax; ++t){
        // PREDICTION EQUATIONS
        B12 = n_ident - B;
        B13 = x - B0;
        B13 = arma::kron(B13.t(), n_ident);
        BB = arma::join_rows(B, B12);
        BB = arma::join_rows(BB, B13);
        BB_temp = arma::join_rows(arma::zeros(n, n), n_ident);
        BB_temp = arma::join_rows(BB_temp, arma::zeros(n, n2));
        BB = arma::join_cols(BB, BB_temp);
        BB = arma::join_cols(BB, arma::join_rows(arma::zeros(n2, 2*n), n2_ident));
        PP = BB * PP * BB.t() + S;
        if(U.n_cols < 2) {
            x = B0 + B * (x - B0);
        } else {
            if (nu == 1) {
                x = B0 + B * (x - B0) + C * U(0,t);
            } else {
                x = B0 + B * (x - B0) + C * U.col(t);
            }
        }
        
        // UPDATING EQUATIONS
        
        if(!X.col(t).has_nan()){
            FF = Z * PP * Z.t() + Su;
            invF = arma::inv(FF);
            y = arma::join_cols(arma::vectorise(x), arma::vectorise(B0));
            y = arma::join_cols(y, arma::vectorise(B.t()));
            v = X.col(t) - Z * y;
            
            y = y + PP * Z.t() * invF * v;
            PP = PP - PP * Z.t() * invF * Z * PP;
            
            x = y.rows(0,(n-1));
            B0 = y.rows(n, (2*n-1));
            B = arma::zeros(n,n);
            j = 2*n;
            for(arma::uword i=0; i<n; ++i) {
                row_iter = y.rows(j,(j+n-1));
                row_iter.reshape(1,n);
                B.row(i) = row_iter;
                j += n;
            }

            // TERMS OF LIKELIHOOD FUNCTION

            double logdetFF = cpp_log_det(FF);
            logFt += logdetFF;
            vFv = vFv + v.t() * invF * v;
        }
    }
    
    arma::mat LL = logFt + vFv;
    
    // My attempt at including check for imaginary numbers
    // arma::cx_mat LL_i = logFt + vFv;
    // arma::mat LL = arma::zeros(1,1);
    // if(LL_i(0,0).imag() != 0){
    //     LL(0,0) = std::pow(10,10);
    // } else {
    //     LL = arma::real(LL_i);
    // }
    
    return(LL);
}

