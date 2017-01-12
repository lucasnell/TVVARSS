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
double cpp_log_det(arma::mat x) {
    double y = arma::det(x);
    double z = std::abs(y);
    double log_z = std::log(z);
    return(log_z);
}

// Same as `determinant(FF)$modulus[1]` in R
arma::cx_double cx_cpp_log_det(arma::cx_mat x) {
    arma::cx_double y = arma::det(x);
    arma::cx_double z = std::abs(y);
    arma::cx_double log_z = std::log(z);
    return(log_z);
}

// Same as `matrix(x, byrow = TRUE)` in R
arma::mat mat_byrow(arma::vec V, arma::uword nrow, arma::uword ncol) {
    if(V.n_elem != (nrow * ncol)){
        Rcpp::stop("Length of V != nrow * ncol");
    }
    arma::mat M = arma::zeros(nrow, ncol);
    arma::mat row_iter;
    arma::uword j = 0;
    for(arma::uword i=0; i<nrow; ++i) {
        row_iter = V.subvec(j, (j + ncol - 1));
        row_iter.reshape(1,ncol);
        M.row(i) = row_iter;
        j += ncol;
    }
    return(M);
}

// Same as `matrix(x, byrow = TRUE)` in R
arma::cx_mat cx_mat_byrow(arma::cx_vec V, arma::uword nrow, arma::uword ncol) {
    if(V.n_elem != (nrow * ncol)){
        Rcpp::stop("Length of V != nrow * ncol");
    }
    arma::cx_mat M = arma::zeros<arma::cx_mat>(nrow, ncol);
    arma::cx_mat row_iter;
    arma::uword j = 0;
    for(arma::uword i=0; i<nrow; ++i) {
        row_iter = V.subvec(j, (j + ncol - 1));
        row_iter.reshape(1,ncol);
        M.row(i) = row_iter;
        j += ncol;
    }
    return(M);
}


// Checks if there's a complex number in a matrix
bool cx_present(arma::cx_mat M){
    const arma::cx_mat& cM = M;
    arma::uword nrow = M.n_rows;
    arma::uword ncol = M.n_cols;
    double im;
    bool out = false;
    for (arma::uword i=0; i<nrow; ++i){
        for (arma::uword j=0; j<ncol; ++j){
            im = cM(i,j).imag();
            if(im != 0){
                out = true;
                return(out);
            }
        }
    }
    return(out);
}





// [[Rcpp::export]]
arma::mat cpp_TVVARSS_ml(arma::vec par, arma::mat X, arma::mat U, arma::vec par_fixed) {
    
    // Defining output type now, in case of complex numbers
    arma::mat LL = arma::zeros(1,1);
    LL(0,0) = std::pow(10,10);
    // This may be used more than once, but definitely will be used at the end
    bool is_cx;
    
    // Note: X and U are transposed, so time runs through columns
    arma::uword n = X.n_rows;
    arma::uword Tmax = X.n_cols;
    
    // Identity matrix of size n by n
    arma::cx_mat n_ident = arma::eye<arma::cx_mat>(n,n);
    // And for n^2 by n^2
    arma::sword n2 = std::pow(n,2);
    arma::cx_mat n2_ident = arma::eye<arma::cx_mat>(n2,n2);
    
    arma::cx_vec par_full = arma::conv_to<arma::cx_vec>::from(par_fixed);
    arma::uvec par_nas = arma::find_nonfinite(par_fixed);
    par_full.elem(par_nas) = arma::conv_to<arma::cx_vec>::from(par);
    
    // =============
    // set up coefficient matrices
    // =============
    
    arma::cx_mat B0;
    B0.insert_cols(0, par_full.subvec(0,n-1));
    B0.reshape(n, 1);
    
    arma::cx_mat B = cx_mat_byrow(par_full.subvec(n, (n + n2 - 1)), n, n);
    
    arma::cx_vec Se_vec = par_full.subvec((n+n2),(n+n2+n-1));
    Se_vec = arma::pow(Se_vec, 2);
    arma::cx_mat Se = arma::diagmat(Se_vec);
    
    arma::cx_vec Su_vec = par_full.subvec((n+n2+n),(n+n2+n+n-1));
    Su_vec = arma::pow(Su_vec, 2);
    arma::cx_mat Su = arma::diagmat(Su_vec);
    
    arma::cx_vec Sb_vec = par_full.subvec((n+n2+n+n), (n+n2+n+n+n*(n+1))-1);
    Sb_vec = arma::pow(Sb_vec, 2);
    arma::cx_mat Sb = arma::diagmat(Sb_vec);
    
    
    // =============
    // set up independent variable
    // =============
    arma::uword nu = U.n_rows;
    arma::cx_mat C;
    
    if(U.n_cols > 1) {
        C = cx_mat_byrow(par_full.subvec((n+n2+n+n+n*(n+1)), (n+n2+n+n+n*(n+1)+nu*n-1)), 
                      n, nu);
    }
    
    arma::cx_mat S;
    arma::cx_mat S_temp = arma::zeros<arma::cx_mat>(n, n*(n+1));
    S = arma::join_rows(Se, S_temp);
    S_temp = arma::zeros<arma::cx_mat>(n*(n+1), n);
    S_temp = arma::join_rows(S_temp, Sb);
    S = arma::join_cols(S, S_temp);
    
    arma::cx_mat Z = arma::join_rows(n_ident, arma::zeros<arma::cx_mat>(n, n*(n+1)));
    
    // =============
    // Initial unconditional values
    // =============
    arma::cx_vec x = arma::conv_to<arma::cx_vec>::from(X.col(0));
    
    
    // If the initial parameter values imply a stationary distribution, then the initial
    // covariance matrix of the process error, PP, is computed from the coefficients. If 
    // not, the initial value of PP given by the covariance matrix of the process error
    // variation (effectively assuming that the dominant eigenvalue of the system is 
    // zero).
    arma::cx_vec eigvals = arma::eig_gen(B);
    arma::vec eig_abs = abs(eigvals);
    // is_cx = cx_present(arma::conv_to<arma::cx_mat>::from(eigvals));
    // if(is_cx){
    //     return(LL);
    // }

    arma::cx_mat PP;
    arma::cx_mat PP_cx;
    // arma::cx_mat B_cx = arma::conv_to<arma::cx_mat>::from(B);
    arma::cx_mat B_cx = B;
    arma::cx_mat PP_k;
    arma::cx_mat PP_Se = Se;
    arma::cx_mat PP_dm;
    if(max(eig_abs) < 1){
        PP_k = arma::kron(B_cx, B_cx);
        PP_Se.reshape(n*n,1);
        PP_dm = n2_ident;
        PP_cx = arma::inv(PP_dm - PP_k);
        PP_cx = PP_cx * PP_Se;
        PP_cx.reshape(n, n);
        // is_cx = cx_present(PP_cx);
        // if (is_cx){
        //     return(LL);
        // }
        // // PP = arma::conv_to<arma::cx_mat>::from(PP_cx);
        PP = PP_cx;
    } else {
        PP = Se;
    }
    
    PP = arma::join_rows(PP, arma::zeros<arma::cx_mat>(n, n*(n+1)));
    arma::cx_mat PP_Sb = arma::join_rows(arma::zeros<arma::cx_mat>(n*(n+1), n), Sb);
    PP = arma::join_cols(PP, PP_Sb);
    
    arma::cx_double logFt = 0;
    arma::cx_mat vFv = arma::zeros<arma::cx_mat>(1,1);
    
    arma::cx_mat BB;
    arma::cx_mat B12;
    arma::cx_mat B13;
    arma::cx_mat B13_cx;
    arma::cx_mat BB_temp;
    arma::cx_mat FF;
    arma::cx_mat invF;
    arma::cx_mat y;
    arma::cx_mat v;
    arma::cx_double logdetFF;
    for(arma::uword t=1; t<Tmax; ++t){
        // PREDICTION EQUATIONS
        B12 = n_ident - B;
        B13 = x - B0;
        B13_cx = arma::conv_to<arma::cx_mat>::from(B13);
        B13_cx = arma::kron(B13_cx.t(), n_ident);
        // is_cx = cx_present(B13_cx);
        // if (is_cx){
        //     return(LL);
        // }
        B13 = arma::conv_to<arma::cx_mat>::from(B13_cx);
        BB = arma::join_rows(B, B12);
        BB = arma::join_rows(BB, arma::conv_to<arma::cx_mat>::from(B13));
        BB_temp = arma::join_rows(arma::zeros<arma::cx_mat>(n, n), n_ident);
        BB_temp = arma::join_rows(BB_temp, arma::zeros<arma::cx_mat>(n, n2));
        BB = arma::join_cols(BB, BB_temp);
        BB = arma::join_cols(BB, arma::join_rows(arma::zeros<arma::cx_mat>(n2, 2*n), 
                                                 n2_ident));
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

            B = cx_mat_byrow(y(span(2*n, (y.n_rows - 1)), 0), n, n);
            
            // TERMS OF LIKELIHOOD FUNCTION
            logdetFF = cx_cpp_log_det(FF);
            is_cx = logdetFF.imag() != 0;
            if (is_cx){
                return(LL);
            }
            logFt += logdetFF;
            vFv = vFv + v.t() * invF * v;
        }
    }
    
    // LL = logFt + vFv;
    arma::cx_mat LL_cx = logFt + vFv;
    is_cx = cx_present(LL_cx);
    if (is_cx){
        return(LL);
    }
    LL = arma::conv_to<arma::mat>::from(LL_cx);
    
    return(LL);
}

