
# R version of TVVARSS.ml

TVVARSS.ml <- function(par, X, U, par.fixed) {
    
    # Note: X and U are transposed, so time runs through columns
    n <- dim(X)[1]
    Tmax <- dim(X)[2]
    
    largenumber <-10.0^10
    
    par.full <- par.fixed
    par.full[is.na(par.fixed)] <- par
    
    # set up coefficient matrices
    B0 <- matrix(par.full[1:n], nrow=n, ncol=1)
    B <- matrix(par.full[(n+1):(n+n^2)], nrow=n, ncol=n, byrow = TRUE)
    Se <- diag(par.full[(n+n^2+1):(n+n^2+n)]^2)
    Su <- diag(par.full[(n+n^2+n+1):(n+n^2+n+n)]^2)
    Sb <- diag(par.full[(n+n^2+n+n+1):(n+n^2+n+n+n*(n+1))]^2)
    
    # set up independent variable
    if(!is.null(U)){
        nu <- dim(U)[1]
        C <- matrix(par.full[(n+n^2+n+n+n*(n+1)+1):(n+n^2+n+n+n*(n+1)+nu*n)], 
                    nrow=n, ncol=nu, byrow = TRUE)
    }		
    
    S <- as.matrix(rbind(cbind(Se, matrix(0, n, n*(n+1))), 
                         cbind(matrix(0, n*(n+1), n), Sb)))		
    Z <- as.matrix(cbind(diag(n), matrix(0, n, n*(n+1))))
    
    # Initial unconditional values
    x <- X[,1]
    
    # If the initial parameter values imply a stationary distribution, then the initial
    # covariance matrix of the process error, PP, is computed from the coefficients. If 
    # not, the initial value of PP given by the covariance matrix of the process error
    # variation (effectively assuming that the dominant eigenvalue of the system is zero).
    B_eigen <- eigen(B)$values
    # if (is.complex(B_eigen)) {
    #     return(matrix(10^10))
    # }
    if (max(abs(B_eigen)) < 1) {
        PP <- solve(diag(n*n)-kronecker(B,B)) %*% matrix(Se, nrow = n*n)		
        PP <- matrix(PP, n, n)
        # if (is.complex(PP)) {
        #     return(matrix(10^10))
        # }
    } else {
        PP <- Se
    }
    PP <- as.matrix(rbind(cbind(PP, matrix(0, n, n*(n+1))), 
                          cbind(matrix(0, n*(n+1), n), Sb)))
    
    logFt <- 0
    vFv <- 0
    
    for(t in 2:Tmax) {
        
        # PREDICTION EQUATIONS
        
        B12 <- diag(n) - B
        B13 <- kronecker(t(x-B0),diag(n))
        # if (is.complex(B13)) {
        #     return(matrix(10^10))
        # }
        BB <- as.matrix(rbind(cbind(B, B12, B13), 
                              cbind(matrix(0, n, n), diag(n), matrix(0, n, n^2)), 
                              cbind(matrix(0, n^2, 2*n), diag(n^2))))
        PP <- BB %*% PP %*% t(BB) + S
        
        if (is.null(U)){
            x <- B0 + B %*% (x-B0)
        } else {				
            if (nu == 1) {
                x <- B0 + B %*% (x-B0) + C * U[t]
            } else {
                x <- B0 + B %*% (x-B0) + C %*% U[,t]
            }
        }
        
        # UPDATING EQUATIONS
        if(!any(is.na(X[,t]))){
            
            FF <- Z %*% PP %*% t(Z) + Su
            # "determinant(FF)$modulus[1]" gives the log of the determinant by default
            logdetFF <- determinant(FF)$modulus[1]
            if (logdetFF == -Inf) return(largenumber)
            invF <- solve(FF)
            
            y <- matrix(c(x, B0, t(B)), ncol=1)		
            v <- X[,t] - Z %*% y
            
            y <- y + PP %*% t(Z) %*% invF %*% v
            PP <- PP - PP %*% t(Z) %*% invF %*% Z %*% PP
            
            x <- y[1:n]
            B0 <- y[(n+1):(2*n)]
            B <- matrix(y[(2*n+1):length(y)], nrow=n, ncol=n, byrow = TRUE)
            
            # TERMS OF LIKELIHOOD FUNCTION
            

            # if (is.complex(logdetFF)){
            #     return(matrix(10^10))
            # }
            logFt <- logFt + logdetFF
            
            vFv <- vFv + t(v) %*% invF %*% v
        }
    }
    
    LL <- as.numeric(logFt + vFv)
    if (is.complex(LL)) {
        cat(LL, ' large \n')
        # system(paste('echo', LL, '>> LLs.txt'))
        LL <- largenumber
    } else{
        # cat(LL, '\n')
    }
    # show(c(LL,par.full[1:6]))
    # browser()
    return(LL)
}