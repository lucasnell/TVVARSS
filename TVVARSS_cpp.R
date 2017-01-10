library(Rcpp)
library(RcppArmadillo)
library(GenSA)

sourceCpp('TVVARSS.cpp')




cpp_TVVARSS <- function(X, U = NULL, B0.start = matrix(NA, nrow=n, ncol=1), B.start = matrix(NA, nrow=n, ncol=n), se.start = matrix(NA, nrow=n, ncol=1), su.start = matrix(0.001, nrow=n, ncol=1), Sb0.start = matrix(0.1, nrow=n, ncol=1), Sb.start = matrix(0.1, nrow=n, ncol=n), C.start = if(is.null(U)) NULL else matrix(0.1, nrow=n, ncol=ncol(U)), B0.fixed = matrix(NA, nrow=n, ncol=1), B.fixed = matrix(NA, nrow=n, ncol=n), se.fixed = matrix(NA, nrow=n, ncol=1), su.fixed = matrix(NA, nrow=n, ncol=1), Sb0.fixed = matrix(NA, nrow=n, ncol=1), Sb.fixed = matrix(NA, nrow=n, ncol=n), C.fixed = if(is.null(U)) NULL else matrix(NA, nrow=n, ncol=ncol(U)), Tsamplefract = .5, annealing = T, method = 'Nelder-Mead', show.fig = T, maxit.optim = 10^5, maxit.SANN = 10^2, optim.control = NULL, optim.bobyqa.control = NULL, optim.SANN.control = NULL) {
    
    
    # TVVARSS.ml
    
    ####################################################
    # The function below mirrors the one above but saves output and makes plots. It should stay in R.
    ####################################################
    ####################################################
    # Begin TVVARSS.fit
    TVVARSS.fit <- function(par, X, U, par.fixed, show.fig) {
        
        if (all(is.na(U))){
            U <- NULL
        }
        
        # Note: X and U are transposed, so time runs through columns
        n <- dim(X)[1]
        Tmax <- dim(X)[2]
        
        par.full <- par.fixed
        par.full[is.na(par.fixed)] <- par
        
        B0 <- matrix(par.full[1:n], nrow=n, ncol=1)
        B <- matrix(par.full[(n+1):(n+n^2)], nrow=n, ncol=n, byrow = TRUE)
        Se <- diag(par.full[(n+n^2+1):(n+n^2+n)]^2)
        Su <- diag(par.full[(n+n^2+n+1):(n+n^2+n+n)]^2)	
        Sb <- diag(par.full[(n+n^2+n+n+1):(n+n^2+n+n+n*(n+1))]^2)
        
        if(!is.null(U)){
            nu <- dim(U)[1]
            C <- matrix(par.full[(n+n^2+n+n+n*(n+1)+1):(n+n^2+n+n+n*(n+1)+nu*n)], nrow=n, ncol=nu, byrow = TRUE)
        }		
        
        S <- as.matrix(rbind(cbind(Se, matrix(0, n, n*(n+1))), cbind(matrix(0, n*(n+1), n), Sb)))		
        Z <- as.matrix(cbind(diag(n),matrix(0, n, n*(n+1))))
        
        # Initial unconditional values
        x <- X[,1]
        
        # If the initial parameter values imply a stationary distribution, then the initial covariance matrix of the process error, PP, is computed from the coefficients. If not, the initial value of PP given by the covariance matrix of the process error variation (effectively assuming that the dominant eigenvalue of the system is zero).
        if(max(abs(eigen(B)$values)) < 1) {
            PP <- solve(diag(n*n)-kronecker(B,B)) %*% matrix(Se, nrow = n*n)		
            PP <- matrix(PP, n, n)
        }else{
            PP <- Se
        }
        PP <- as.matrix(rbind(cbind(PP, matrix(0, n, n*(n+1))), matrix(0, n*(n+1), n*(2+n))))
        
        logFt <- 0
        vFv <- 0
        
        # Information used for plotting
        eiglist <- matrix(NA, nrow=n, ncol=Tmax)
        eiglist[,1] <- eigen(B)$values
        
        xlist <- matrix(NA, nrow=n, ncol=Tmax)
        xlist[,1] <- x
        
        b0list <- matrix(NA, nrow=n, ncol=Tmax)
        b0list[,1] <- B0
        
        blist <- matrix(NA, nrow=n^2, ncol=Tmax)
        blist[,1] <- matrix(t(B), nrow=n^2, ncol=1)
        
        PPlist <- matrix(NA, nrow=(2*n+n^2)^2, ncol=Tmax)
        PPlist[,1] <- matrix(t(PP), nrow=(2*n+n^2)^2, ncol=1)
        
        for(t in 2:Tmax) {
            
            # PREDICTION EQUATIONS
            
            B12 <- diag(n) - B
            B13 <- kronecker(t(x-B0),diag(n))
            
            BB <- as.matrix(rbind(cbind(B, B12, B13), cbind(matrix(0, n, n), diag(n), matrix(0, n, n^2)), cbind(matrix(0, n^2, 2*n), diag(n^2))))
            PP <- BB %*% PP %*% t(BB) + S
            
            if(is.null(U)){
                x <- B0 + B %*% (x-B0)
            }else{				
                if(nu == 1) {
                    x <- B0 + B %*% (x-B0) + C * U[t]
                }else{
                    x <- B0 + B %*% (x-B0) + C %*% U[,t]
                }
            }		
            
            # UPDATING EQUATIONS
            if(!any(is.na(X[,t]))){
                FF <- Z %*% PP %*% t(Z) + Su
                invF <- solve(FF)
                
                y <- matrix(c(x, B0, t(B)), ncol=1)		
                v <- X[,t] - Z %*% y
                
                y <- y + PP %*% t(Z) %*% invF %*% v
                PP <- PP - PP %*% t(Z) %*% invF %*% Z %*% PP
                
                x <- y[1:n]
                B0 <- y[(n+1):(2*n)]
                B <- matrix(y[(2*n+1):length(y)], nrow=n, ncol=n, byrow = TRUE)
                
                # TERMS OF LIKELIHOOD FUNCTION
                logdetFF <- determinant(FF)$modulus[1]
                logFt <- logFt + logdetFF
                
                vFv <- vFv + t(v) %*% invF %*% v
                
                eiglist[,t] <- eigen(B)$values
                xlist[,t] <- x
                b0list[,t] <- B0
                blist[,t] <- matrix(t(B), nrow=n^2, ncol=1)
                PPlist[,t] <- matrix(t(PP), nrow=(2*n+n^2)^2, ncol=1)
            }
        }
        
        LL <- logFt + vFv
        
        if(is.complex(LL)) LL <- 10^10
        
        if(show.fig == T){
            
            npar <- length(par)
            logLik <- -(n * (Tmax - 1)/2) * log(2*pi) - LL/2
            
            par(mfrow=c(2,1), mar=c(4, 4, 2, .5), cex.axis=1, cex.lab=1, cex.main=1)
            matplot(1:Tmax, t(X), typ = "p", pch = 1, col = "black", xlab="", ylab="X (o), fitted.X (b), and fitted.mean (g)", main=paste('logLik=', .001*round(1000*logLik)))
            matlines(1:Tmax, t(xlist), typ = "l", lty = 1, col = "blue")
            if(any(is.na(X))) matlines(1:Tmax, t(xlist), typ = "p", pch = 20, cex = 0.5, col = "blue")
            
            if(is.null(U)){
                matlines(1:Tmax, t(b0list), col = "green", lty = 1)
                if(any(is.na(X))) matlines(1:Tmax, t(b0list), typ = "p", col = "green", pch = 20, cex = 0.5)
            }else{				
                matlines(1:Tmax, t(b0list + C %*% U), col = "green", lty = 1)
                if(any(is.na(X))) matlines(1:Tmax, t(b0list + C %*% U), typ = "p", col = "green", pch = 20, cex = 0.5)
            }		
            
            maxEigs <- apply(abs(eiglist), 2, pmax)[1,]
            plot(1:Tmax, maxEigs, typ="l", col = "red", xlab = "Time", ylab = expression(lambda), ylim = c(0,1.3), main=paste('max.Eig=', .001*round(1000*max(maxEigs))))
            if(any(is.na(X))) points(1:Tmax, maxEigs, col = "red", pch = 20, cex = 0.5)
            lines(c(0, Tmax), c(1, 1), col="black")
        }		
        
        X.fitted <<- xlist
        B0.fitted <<- b0list
        B.fitted <<- blist
        PP.fitted <<- PPlist
        eigen.fitted <<- eiglist		
        
        return(LL)
    }			
    # End TVVARSS.fit
    ####################################################
    
    if(!is.matrix(X)){
        stop("X should be a matrix.")
    }
    if(!is.null(U)) if(!is.matrix(U)){
        stop("U should be a matrix.")
    }
    if (any(var(X, na.rm = TRUE) == 0)) {
        stop("The response (dependent variable) has no variation.")
    }
    
    X <- as.matrix(X)
    if (dim(X)[2] > dim(X)[1]) {
        stop("Time should run in the vertical direction.")
    }
    n <- dim(X)[2]
    
    if(!is.null(U)){
        if(nrow(U) != nrow(X)) stop("Variables X and U have different lengths.")
        U <- as.matrix(U)
        nu <- dim(U)[2]
    }else{
        nu <- NULL
    }
    
    Tmax <- dim(X)[1]	
    Tsample <- floor(Tsamplefract * Tmax)
    
    # set up initial values of parameters
    B0.init <- matrix(NA, nrow = n, ncol = 1)
    B.init <- matrix(NA, nrow = n, ncol = n)
    if(!is.null(U)) {
        C.init <- matrix(NA, nrow = n, ncol = nu)
    }else{
        C.init <- NULL
    }
    se.init <- matrix(NA, nrow = n, ncol = 1)
    for(i in 1:n){
        y <- X[2:Tsample, i]
        x <- X[1:(Tsample-1),]
        if(!is.null(U)) {
            if(nu == 1){
                u <- U[1:(Tsample-1)]
                df <- as.data.frame(cbind(y, x, u))
                z.lm <- lm(y ~ ., df, na.action = na.omit)
                B0.init[i] <- z.lm$coef[1]
                B.init[i,] <- z.lm$coef[2:(n+1)]
                C.init[i] <- z.lm$coef[(n+2):(n+1+nu)]
                se.init[i] <- var(z.lm$resid)^.5
            }else{
                u <- U[1:(Tsample-1),]
                df <- as.data.frame(cbind(y, x, u))
                z.lm <- lm(y ~ ., df, na.action = na.omit)
                B0.init[i] <- z.lm$coef[1]
                B.init[i,] <- z.lm$coef[2:(n+1)]
                C.init[i,] <- z.lm$coef[(n+2):(n+1+nu)]
                se.init[i] <- var(z.lm$resid)^.5			}
        }else{
            df <- as.data.frame(cbind(y, x))
            z.lm <- lm(y ~ ., df, na.action = na.omit)
            B0.init[i] <- z.lm$coef[1]
            B.init[i,] <- z.lm$coef[2:(n+1)]
            se.init[i] <- var(z.lm$resid)^.5
        }
    }
    
    # variables.init
    b0.init <- array(B0.init, dim=n)
    b.init <- array(t(B.init), dim=n^2)
    se.init <- array(se.init, dim=n)
    su.init <- array(su.start, dim=n)
    sb0.init <- array(t(Sb0.start), dim=n)
    sb.init <- array(t(Sb.start), dim=n*n)
    if(!is.null(U)) {
        c.init <- array(t(C.init), dim=n*nu)
    }else{
        c.init <- NULL
    }
    par.init <- c(b0.init, b.init, se.init, su.init, sb0.init, sb.init, c.init) 
    
    # variables.start	
    b0.start <- array(B0.start, dim=n)
    b.start <- array(t(B.start), dim=n^2)
    se.start <- array(se.start, dim=n)
    su.start <- array(su.start, dim=n)
    sb0.start <- array(t(Sb0.start), dim=n)
    sb.start <- array(t(Sb.start), dim=n*n)
    if(!is.null(U)) {
        c.start <- array(t(C.start), dim=n*nu)
    }else{
        c.start <- NULL
    }
    par.start <- c(b0.start, b.start, se.start, su.start, sb0.start, sb.start, c.start)
    
    # variables.fixed	
    b0.fixed <- array(B0.fixed, dim=n)
    b.fixed <- array(t(B.fixed), dim=n^2)
    se.fixed <- array(se.fixed, dim=n)
    su.fixed <- array(su.fixed, dim=n)
    sb0.fixed <- array(t(Sb0.fixed), dim=n)
    sb.fixed <- array(t(Sb.fixed), dim=n*n)
    if(!is.null(U)) {
        c.fixed <- array(t(C.fixed), dim=n*nu)
    }else{
        c.fixed <- NULL
    }
    par.fixed <- c(b0.fixed, b.fixed, se.fixed, su.fixed, sb0.fixed, sb.fixed, c.fixed)	
    
    # set up variables for fitting
    par.full <- par.init
    par.full[!is.na(par.start)] <- par.start[!is.na(par.start)]
    
    par.full[!is.na(par.fixed)] <- par.fixed[!is.na(par.fixed)]	
    par <- par.full[is.na(par.fixed)]
    
    counter <- 0
    tX <- t(X)
    if(!is.null(U)) {
        tU <- as.matrix(t(U))
    }else{
        tU <- matrix()
    }
    
    if(!is.element(method, c("Nelder-Mead", "bobyqa", "BFGS"))) stop("Acceptable methods are Nelder-Mead {optim},  BFGS {optim}, and bobyqa {minqa}.")
    
    if(annealing == T) require(GenSA)
    if(method == "bobyqa") require(minqa)
    
    if(is.null(optim.control)) optim.control = list(maxit = maxit.optim)
    if(is.null(optim.bobyqa.control)) optim.bobyqa.control = list(maxfun = maxit.optim)
    if(is.null(optim.SANN.control)) optim.SANN.control = list(temp = 0.01, nb.stop.improvement = 40, maxit = maxit.SANN)	
    
    if(annealing == T){
        
        if(!is.null(U)) {
            par.upper <- c(array(1, dim=n), array(2, dim=n^2), array(10, dim=n), array(10, dim=n), array(10, dim=n*(n+1)), array(100, dim=n*nu))
            par.lower <- c(array(-1, dim=n), array(-2, dim=n^2), array(0, dim=n), array(0, dim=n), array(0, dim=n*(n+1)), array(-100, dim=n*nu))
        }else{
            par.upper <- c(array(1, dim=n), array(2, dim=n^2), array(10, dim=n), array(10, dim=n), array(10, dim=n*(n+1)))
            par.lower <- c(array(-1, dim=n), array(-2, dim=n^2), array(0, dim=n), array(0, dim=n), array(0, dim=n*(n+1)))
        }
        
        par.upper <- par.upper[is.na(par.fixed)]
        par.lower <- par.lower[is.na(par.fixed)]
        
        optSANN <- GenSA(fn = cpp_TVVARSS_ml, par = par, lower = par.lower, upper = par.upper, X = tX, U = tU, par_fixed = par.fixed, control=list(smooth = F, maxit = maxit.SANN))	
        par <- optSANN$par
    }
    
    if(method == 'Nelder-Mead') {
        opt <- optim(fn = cpp_TVVARSS_ml, par = par, X = tX, U = tU, par_fixed = par.fixed, method = "Nelder-Mead", control = optim.control)
        if(opt$convergence != 0) cat("/nNelder-Mead optimization failed to converge")
        opt.convergence <- opt$convergence
    }
    
    if(method == 'BFGS') {
        opt <- optim(fn = cpp_TVVARSS_ml, par = par, X = tX, U = tU, par_fixed = par.fixed, method = "BFGS", control = optim.control)
        if(opt$convergence != 0) cat("/nBFGS optimization failed to converge")
        opt.convergence <- opt$convergence
    }
    
    if(method == 'bobyqa') {
        opt <- bobyqa(fn = cpp_TVVARSS_ml, par = par, X = tX, U = tU, par_fixed = par.fixed, control = optim.bobyqa.control)
        if(opt$ierr != 0) cat("/nbobyqa optimization failed to converge")
        opt.convergence <- opt$ierr
        opt$value <- opt$fval
    }
    
    # retrieve final fitted values
    par <- opt$par
    TVVARSS.fit(par = par, X = tX, U = tU, par.fixed = par.fixed, show.fig = show.fig)
    
    par.full <- par.fixed
    par.full[is.na(par.fixed)] <- par
    
    B0 <- matrix(par.full[1:n], nrow=n, ncol=1)
    B <- matrix(par.full[(n+1):(n+n^2)], nrow=n, ncol=n, byrow = TRUE)
    se <- abs(matrix(par.full[(n+n^2+1):(n+n^2+n)], nrow=n, ncol=1))
    su <- abs(matrix(par.full[(n+n^2+n+1):(n+n^2+n+n)], nrow=n, ncol=1))
    Sb0 <- abs(matrix(par.full[(n+n^2+n+n+1):(n+n^2+n+n+n)], nrow=n, ncol=1))
    Sb <- abs(matrix(par.full[(n+n^2+n+n+n+1):(n+n^2+n+n+n+n*n)], nrow=n, ncol=n, byrow = TRUE))
    
    if(!is.null(U)){
        C <- matrix(par.full[(n+n^2+n+n+n*(n+1)+1):(n+n^2+n+n+n*(n+1)+nu*n)], nrow=n, ncol=nu, byrow = TRUE)
    }else{
        C <- NULL
    }
    
    LL <- opt$value
    npar <- length(par)
    logLik <- -(n * (Tmax - 1)/2) * log(2*pi) - LL/2
    AIC <- -2*logLik + 2*npar;
    
    results <- list(X = X, U = U, n = n, se = se, su = su, Sb0 = Sb0, Sb = Sb, C = C, B0 = B0, B = B, logLik = logLik, AIC = AIC, npar = npar, B0.fitted = B0.fitted, B.fitted = B.fitted, X.fitted = X.fitted, PP.fitted = PP.fitted, eigen.fitted = eigen.fitted, B0.init = B0.init, B.init = B.init, se.init = se.init, C.init = C.init, B0.start = B0.start, B.start = B.start, se.start = se.start, su.start = su.start, Sb0.start = Sb0.start, Sb.start = Sb.start, C.start = C.start, B0.fixed = B0.fixed, B.fixed = B.fixed, se.fixed = se.fixed, su.fixed = su.fixed, Sb0.fixed = Sb0.fixed, Sb.fixed = Sb.fixed, C.fixed = C.fixed, par.fixed = par.fixed, opt.par = par, Tsamplefract = Tsamplefract, annealing = annealing, optim.control = optim.control, optim.SANN.control = optim.SANN.control, opt.converge = opt.convergence)
    
    class(results) <- "TVVARSS"
    
    rm(list=c("X.fitted", "B0.fitted", "B.fitted", "PP.fitted", "eigen.fitted"), envir = .GlobalEnv)
    
    return(results)
}





