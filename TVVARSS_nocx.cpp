# include <RcppArmadillo.h>

// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;


// Same as `determinant(FF)$modulus[1]` in R
double cpp_log_det(mat x) {
    double y = det(x);
    double z = std::abs(y);
    double log_z = std::log(z);
    return(log_z);
}

// Same as `matrix(x, byrow = TRUE)` in R
mat mat_byrow(vec V, uword nrow, uword ncol) {
    if(V.n_elem != (nrow * ncol)){
        Rcpp::stop("Length of V != nrow * ncol");
    }
    mat M = zeros(nrow, ncol);
    mat row_iter;
    uword j = 0;
    for(uword i=0; i<nrow; ++i) {
        row_iter = V.subvec(j, (j + ncol - 1));
        row_iter.reshape(1,ncol);
        M.row(i) = row_iter;
        j += ncol;
    }
    return(M);
}



// Checks if there's a complex number in a matrix
short int cx_present(cx_mat M){
    const cx_mat& cM = M;
    uword nrow = M.n_rows;
    uword ncol = M.n_cols;
    double im;
    short int out = 0;
    for (uword i=0; i<nrow; ++i){
        for (uword j=0; j<ncol; ++j){
            im = cM(i,j).imag();
            if (im != 0.0){
                out = 1;
                break;
            }
        }
    }
    return out;
}





// [[Rcpp::export]]
double cpp_TVVARSS_ml(vec par, mat X, mat U, vec par_fixed) {
    
    // Defining output type now, in case of complex numbers
    double LL = std::pow(10,10);
    // This may be used more than once, but definitely will be used at the end
    short int is_cx;
    
    // Note: X and U are transposed, so time runs through columns
    uword n = X.n_rows;
    uword Tmax = X.n_cols;
    
    // Identity matrix of size n by n
    mat n_ident = eye<mat>(n,n);
    // And for n^2 by n^2
    sword n2 = n * n;
    mat n2_ident = eye<mat>(n2,n2);
    
    vec par_full = conv_to<vec>::from(par_fixed);
    uvec par_nas = find_nonfinite(par_fixed);
    par_full.elem(par_nas) = conv_to<vec>::from(par);
    
    // =============
    // set up coefficient matrices
    // =============
    
    mat B0;
    B0.insert_cols(0, par_full.subvec(0,n-1));
    B0.reshape(n, 1);
    
    mat B = mat_byrow(par_full.subvec(n, (n + n2 - 1)), n, n);
    
    vec Se_vec = par_full.subvec((n+n2),(n+n2+n-1));
    Se_vec = pow(Se_vec, 2);
    mat Se = diagmat(Se_vec);
    
    vec Su_vec = par_full.subvec((n+n2+n),(n+n2+n+n-1));
    Su_vec = pow(Su_vec, 2);
    mat Su = diagmat(Su_vec);
    
    vec Sb_vec = par_full.subvec((n+n2+n+n), (n+n2+n+n+n*(n+1))-1);
    Sb_vec = pow(Sb_vec, 2);
    mat Sb = diagmat(Sb_vec);
    
    
    // =============
    // set up independent variable
    // =============
    uword nu = U.n_rows;
    mat C;
    
    if(U.n_cols > 1) {
        C = mat_byrow(par_full.subvec((n+n2+n+n+n*(n+1)), (n+n2+n+n+n*(n+1)+nu*n-1)), 
                      n, nu);
    }
    
    mat S;
    mat S_temp = zeros<mat>(n, n*(n+1));
    S = join_rows(Se, S_temp);
    S_temp = zeros<mat>(n*(n+1), n);
    S_temp = join_rows(S_temp, Sb);
    S = join_cols(S, S_temp);
    
    mat Z = join_rows(n_ident, zeros<mat>(n, n*(n+1)));
    
    // =============
    // Initial unconditional values
    // =============
    vec x = conv_to<vec>::from(X.col(0));
    
    
    // If the initial parameter values imply a stationary distribution, then the initial
    // covariance matrix of the process error, PP, is computed from the coefficients. If 
    // not, the initial value of PP given by the covariance matrix of the process error
    // variation (effectively assuming that the dominant eigenvalue of the system is 
    // zero).
    vec eigvals = eig_gen(B);
    vec eig_abs = abs(eigvals);
    // is_cx = present(conv_to<mat>::from(eigvals));
    // if(is_cx){
    //     return(LL);
    // }

    mat PP;
    mat PP_cx;
    mat B_cx = B;
    mat PP_k;
    mat PP_Se = Se;
    mat PP_dm;
    if(max(eig_abs) < 1){
        PP_k = kron(B_cx, B_cx);
        PP_Se.reshape(n*n,1);
        PP_dm = n2_ident;
        PP_cx = inv(PP_dm - PP_k);
        PP_cx = PP_cx * PP_Se;
        PP_cx.reshape(n, n);
        // is_cx = present(PP_cx);
        // if (is_cx){
        //     return(LL);
        // }
        PP = PP_cx;
    } else {
        PP = Se;
    }
    
    PP = join_rows(PP, zeros<mat>(n, n*(n+1)));
    mat PP_Sb = join_rows(zeros<mat>(n*(n+1), n), Sb);
    PP = join_cols(PP, PP_Sb);
    
    double logFt = 0;
    mat vFv = zeros<mat>(1,1);
    
    mat BB;
    mat B12;
    mat B13;
    mat B13_cx;
    mat BB_temp;
    mat FF;
    mat invF;
    mat y;
    mat v;
    double logdetFF;
    for(uword t=1; t<Tmax; ++t){
        // PREDICTION EQUATIONS
        B12 = n_ident - B;
        B13 = x - B0;
        B13_cx = B13;
        B13_cx = kron(B13_cx.t(), n_ident);
        // is_cx = present(B13_cx);
        // if (is_cx){
        //     return(LL);
        // }
        B13 = B13_cx;
        BB = join_rows(B, B12);
        BB = join_rows(BB, B13);
        BB_temp = join_rows(zeros<mat>(n, n), n_ident);
        BB_temp = join_rows(BB_temp, zeros<mat>(n, n2));
        BB = join_cols(BB, BB_temp);
        BB = join_cols(BB, join_rows(zeros<mat>(n2, 2*n), 
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
            
            logdetFF = cpp_log_det(FF);
            // is_cx = logdetFF.imag() != 0;
            
            invF = inv(FF);
            y = join_cols(vectorise(x), vectorise(B0));
            y = join_cols(y, vectorise(B.t()));
            v = X.col(t) - Z * y;
            
            y = y + PP * Z.t() * invF * v;
            PP = PP - PP * Z.t() * invF * Z * PP;
            
            x = y.rows(0,(n-1));
            B0 = y.rows(n, (2*n-1));

            B = mat_byrow(y(span(2*n, (y.n_rows - 1)), 0), n, n);
            
            // TERMS OF LIKELIHOOD FUNCTION
            
            // if (is_cx){
            //     return(LL);
            // }
            logFt += logdetFF;
            vFv = vFv + v.t() * invF * v;
        }
    }
    
    // LL = logFt + vFv;
    mat LL_cx = logFt + vFv;
    is_cx = cx_present(LL_cx);
    // std::cout << is_cx << endl;
    
    if (is_cx == 1){
        std::cout << "here 1" << endl;
        LL = pow(10,10);
    } else {
        double LL_i = LL_cx(0,0);
        LL = real(LL_i);
        // std::cout << LL << endl;
    }
    // double LL_i = LL_cx(0,0);
    // LL = real(LL_i);
    return LL;
}

