## TVARSS.R

##   Time-Varying Autoregressive State Space model

## Copyright 2016 Anthony R. Ives


library(Rcpp)
library(RcppArmadillo)
library(GenSA)


set.seed(9)
Tmax <- 100
n <- 2
B0 <- matrix(c(1,1), nrow = n, ncol = 1)
B <- matrix(c(.7,-.25,-.1,.6), nrow = n, ncol = n, byrow = TRUE)
SE <- matrix(c(1,1), nrow = n, ncol = 1)
C <- matrix(c(.0,-.0), nrow = n, ncol=1)

X <- matrix(0, nrow=n, ncol=Tmax)

U <- (1:Tmax) + rnorm(n=Tmax)

for(t in 2:Tmax){
    X[,t] <- B0 + B %*% X[,t-1] + C %*% U[t-1] + SE * rnorm(n)
}
X <- t(X)
X[,1] <- (X[,1] - mean(X[,1]))/sd(X[,1])
X[,2] <- (X[,2] - mean(X[,2]))/sd(X[,2])


set.seed(8); U <- cbind(U, U * rnorm(100))
nu <- if (is.null(ncol(U))) 1 else ncol(U)
# plot(X ~ U)




# =================
# Defined as arguments
# =================
B0.start = matrix(NA, nrow=n, ncol=1)
B.start = matrix(NA, nrow=n, ncol=n)
se.start = matrix(NA, nrow=n, ncol=1)
su.start = matrix(0.001, nrow=n, ncol=1)
Sb0.start = matrix(0.1, nrow=n, ncol=1)
Sb.start = matrix(0.1, nrow=n, ncol=n)
C.start = if(is.null(U)) NULL else matrix(0.1, nrow=n, ncol=ncol(U))

B0.fixed = matrix(NA, nrow=n, ncol=1)
B.fixed = matrix(NA, nrow=n, ncol=n)
se.fixed = matrix(NA, nrow=n, ncol=1)
su.fixed = matrix(NA, nrow=n, ncol=1)
Sb0.fixed = matrix(NA, nrow=n, ncol=1)
Sb.fixed = matrix(NA, nrow=n, ncol=n)
C.fixed = if(is.null(U)) NULL else matrix(NA, nrow=n, ncol=ncol(U))


# =================
# variables.init
# =================
# set up initial values of parameters
Tmax <- nrow(X)
Tsamplefract = 0.5
Tsample <- floor(Tsamplefract * Tmax)
B0.init <- matrix(NA, nrow = n, ncol = 1)
B.init <- matrix(NA, nrow = n, ncol = n)
if(!is.null(U)) {
    C.init <- matrix(NA, nrow = n, ncol = nu)
}else{
    C.init <- NULL
}
se.init <- matrix(NA, nrow = n, ncol = 1)
for(i in 1:dim(X)[2]){
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

# =================
# variables.start	
# =================

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



# =================
# variables.fixed
# =================

b0.fixed <- array(B0.fixed, dim = n)
b.fixed <- array(t(B.fixed), dim=n^2)
se.fixed <- array(se.fixed, dim=n)
su.fixed <- array(su.fixed, dim=n)
sb0.fixed <- array(t(Sb0.fixed), dim=n)
sb.fixed <- array(t(Sb.fixed), dim=n*n)
c.fixed <- array(t(C.fixed), dim=n*nu)

par.fixed <- c(b0.fixed, b.fixed, se.fixed, su.fixed, sb0.fixed, sb.fixed, c.fixed)



# =================
# set up variables for fitting
# =================
par.full <- par.init
par.full[!is.na(par.start)] <- par.start[!is.na(par.start)]

par.full[!is.na(par.fixed)] <- par.fixed[!is.na(par.fixed)]	
par <- par.full[is.na(par.fixed)]


# library(Rcpp)
# library(RcppArmadillo)
# sourceCpp('TVVARSS.cpp')
# optim(fn = cpp_TVVARSS_ml, par = par, X = t(X), U = t(U), par_fixed = par.fixed, method = "Nelder-Mead", control = list(maxit = 10e5))






TVVARSS.ml <- function(par, X, U, par.fixed) {
    
    # Note: X and U are transposed, so time runs through columns
    n <- dim(X)[1]
    Tmax <- dim(X)[2]
    
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
    
    if (max(abs(eigen(B)$values)) < 1) {
        PP <- solve(diag(n*n)-kronecker(B,B)) %*% matrix(Se, nrow = n*n)		
        PP <- matrix(PP, n, n)
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
            invF <- solve(FF)
            
            y <- matrix(c(x, B0, t(B)), ncol=1)		
            v <- X[,t] - Z %*% y
            
            y <- y + PP %*% t(Z) %*% invF %*% v
            PP <- PP - PP %*% t(Z) %*% invF %*% Z %*% PP
            
            x <- y[1:n]
            B0 <- y[(n+1):(2*n)]
            B <- matrix(y[(2*n+1):length(y)], nrow=n, ncol=n, byrow = TRUE)
            
            # TERMS OF LIKELIHOOD FUNCTION
            
            # "determinant(FF)$modulus[1]" gives the log of the determinant by default
            logdetFF <- determinant(FF)$modulus[1]
            logFt <- logFt + logdetFF
            
            vFv <- vFv + t(v) %*% invF %*% v
        }
    }
    
    LL <- logFt + vFv
    if (is.complex(LL)) {
        LL <- 10^10
    }
    # show(c(LL,par.full[1:6]))
    # browser()
    return(LL)
}



# system.time(replicate(100, TVVARSS.ml(par = par, X = t(X), U = t(U), 
#                                       par.fixed = par.fixed)))
# #  user  system elapsed 
# # 1.860   0.015   1.880
# TVVARSS.ml(par = par, X = t(X), U = t(U), par.fixed = par.fixed)
# #          [,1]
# # [1,] 10960.63
# TVVARSS.ml(par = par, X = t(X), U = NULL, par.fixed = par.fixed)
# #          [,1]
# # [1,] 1618.454



# 
# set.seed(0)
# x <- matrix(rnorm(25), nrow = 5)
# y <- matrix(rnorm(5), nrow = 5)
# system.time(replicate(1e3, solve(x)))
# system.time(replicate(1e3, cpp_solve(x)))
# 
# system.time(replicate(10000, determinant(x)$modulus[1]))
# system.time(replicate(10000, cpp_log_det(x)))
# 
# system.time(replicate(1000, kronecker(x, y)))
# system.time(replicate(1000, cpp_kron(x, y)))
# 
# cpp_mmult(x,y)
# x %*% y

# ~~~~~~~~~~~~~~~~
# Things changed:
# ~~~~~~~~~~~~~~~~
# Input U as matrix() instead of NULL


# sourceCpp('TVVARSS.cpp')
# cpp_test(par, t(X), t(U), par.fixed)

par.full <- par.fixed
par.full[is.na(par.fixed)] <- par


# Se
diag(par.full[(n+n^2+1):(n+n^2+n)]^2)
# Su
diag(par.full[(n+n^2+n+1):(n+n^2+n+n)]^2)
# Sb
diag(par.full[(n+n^2+n+n+1):(n+n^2+n+n+n*(n+1))]^2)

# C
C = matrix(par.full[(n+n^2+n+n+n*(n+1)+1):(n+n^2+n+n+n*(n+1)+nu*n)], nrow=n, ncol=nu, byrow = TRUE); C

# S
# as.matrix(rbind(cbind(Se, matrix(0, n, n*(n+1))), 
#                 cbind(matrix(0, n*(n+1), n), Sb)))
S = as.matrix(rbind(cbind(diag(par.full[(n+n^2+1):(n+n^2+n)]^2), 
                      matrix(0, n, n*(n+1))),
                cbind(matrix(0, n*(n+1), n), 
                      diag(par.full[(n+n^2+n+n+1):(n+n^2+n+n+n*(n+1))]^2))))

# Z
Z = as.matrix(cbind(diag(n), matrix(0, n, n*(n+1)))); Z

# x
t(X)[,1]

# Check "If the initial parameter values imply a stationary distribution..."
max(abs(eigen(matrix(par.full[(n+1):(n+n^2)], nrow=n, ncol=n, byrow = TRUE))$values))

# PP
matrix(solve(diag(n*n)-kronecker(matrix(par.full[(n+1):(n+n^2)], nrow=n, ncol=n, byrow = TRUE),matrix(par.full[(n+1):(n+n^2)], nrow=n, ncol=n, byrow = TRUE))) %*% matrix(diag(par.full[(n+n^2+1):(n+n^2+n)]^2), nrow = n*n), n, n)

# PP after manipulation
PP0 = as.matrix(rbind(cbind(matrix(solve(diag(n*n)-kronecker(matrix(par.full[(n+1):(n+n^2)], nrow=n, ncol=n, byrow = TRUE),matrix(par.full[(n+1):(n+n^2)], nrow=n, ncol=n, byrow = TRUE))) %*% matrix(diag(par.full[(n+n^2+1):(n+n^2+n)]^2), nrow = n*n), n, n), matrix(0, n, n*(n+1))), 
                cbind(matrix(0, n*(n+1), n), diag(par.full[(n+n^2+n+n+1):(n+n^2+n+n+n*(n+1))]^2))))


# BB
B12 <- diag(n) - matrix(par.full[(n+1):(n+n^2)], nrow=n, ncol=n, byrow = TRUE)
B13 <- kronecker(t(t(X)[,1]-matrix(par.full[1:n], nrow=n, ncol=1)),diag(n))

BB = as.matrix(rbind(cbind(matrix(par.full[(n+1):(n+n^2)], nrow=n, ncol=n, byrow = TRUE), B12, B13), 
                      cbind(matrix(0, n, n), diag(n), matrix(0, n, n^2)), 
                      cbind(matrix(0, n^2, 2*n), diag(n^2))))


# PP after BB creation
PP = BB %*% PP0 %*% t(BB) + S

# x after first iteration of for loop
# B0 + B %*% (x-B0) + C * U[t]
# # B0 + B %*% (x-B0) + C %*% U[,t]
x0 = matrix(par.full[1:n], nrow=n, ncol=1) + matrix(par.full[(n+1):(n+n^2)], nrow=n, ncol=n, byrow = TRUE) %*% (as.matrix(t(X)[,1]) - matrix(par.full[1:n], nrow=n, ncol=1)) + C %*% as.matrix(t(U)[,2])

B0_0 <- matrix(par.full[1:n], nrow=n, ncol=1)
B_0 <- matrix(par.full[(n+1):(n+n^2)], nrow=n, ncol=n, byrow = TRUE)

FF <- Z %*% PP %*% t(Z) + diag(par.full[(n+n^2+n+1):(n+n^2+n+n)]^2)
invF <- solve(FF)

y0 <- matrix(c(x0, B0_0, t(B_0)), ncol=1)
v <- as.matrix(t(X)[,2]) - Z %*% y0
y <- y0 + PP %*% t(Z) %*% invF %*% v



sourceCpp('TVVARSS.cpp')
# cpp_test(par, t(X), t(U), par.fixed)
# system.time(replicate(100, cpp_TVVARSS_ml(par, t(X), t(U), par.fixed)))
# system.time(replicate(100, TVVARSS.ml(par = par, X = tX, U = tU, par.fixed = par.fixed)))

summary_output <- function(out1, out2, t_vec, name_vec = c('R', 'Rcpp')){
    cat(paste(
        sprintf('%s version took %s minutes',
                name_vec[1], 
                round(as.numeric(t_vec[2] - t_vec[1], units = 'mins'), 2)),
        sprintf('%s version took %s minutes',
                name_vec[2],
                round(as.numeric(t_vec[3] - t_vec[2], units = 'mins'), 2)),
        '\nAre they equal?',
        paste(all.equal(out1, out2), collapse = '\n'),
        sep = '\n'))
}

# cpp_TVVARSS_ml(par = par, X = t(X), U = t(U), par_fixed = par.fixed)
# #          [,1]
# # [1,] 10960.63
# TVVARSS.ml(par = par, X = t(X), U = t(U), par.fixed = par.fixed)
# #          [,1]
# # [1,] 10960.63
# 
# 
# cpp_TVVARSS_ml(par = par, X = t(X), U = matrix(), par_fixed = par.fixed)
# #          [,1]
# # [1,] 1618.454
# TVVARSS.ml(par = par, X = t(X), U = NULL, par.fixed = par.fixed)
# #          [,1]
# # [1,] 1618.454
# 
# 
# t1 <- Sys.time()
# R_optim <- optim(fn = TVVARSS.ml, par = par, X = t(X), U = t(U), par.fixed = par.fixed,
#                  method = "Nelder-Mead", control = list(maxit = 10^5))
# t2 <- Sys.time()
# Rcpp_optim <- optim(fn = cpp_TVVARSS_ml, par = par, X = t(X), U = t(U),
#                     par_fixed = par.fixed, method = "Nelder-Mead",
#                     control = list(maxit = 10^5))
# t3 <- Sys.time()
# summary_output(R_optim, Rcpp_optim, c(t1, t2, t3))
# system2('say', 'R script finished')
# # R version took 2.36 minutes
# # Rcpp version took 0.07 minutes
# # Are they equal? --> TRUE


n <- ncol(X)
maxit.SANN <- 10^2

par.upper <- c(array(1, dim=n), array(2, dim=n^2), array(10, dim=n), array(10, dim=n),
               array(10, dim=n*(n+1)))
par.lower <- c(array(-1, dim=n), array(-2, dim=n^2), array(0, dim=n), array(0, dim=n),
               array(0, dim=n*(n+1)))
par.upper <- par.upper[is.na(par.fixed)]
par.lower <- par.lower[is.na(par.fixed)]

# Added by me to account for using U = NULL
# Same reason I did head(par.fixed, -4) and head(par, -4) below
par.upper <- head(par.upper, -4)
par.lower <- head(par.lower, -4)


# t1 <- Sys.time()
# R_optSANN <- GenSA(fn = TVVARSS.ml, par = head(par, -4), lower = par.lower,
#                    upper = par.upper,
#                    X = t(X), U = NULL, par.fixed = head(par.fixed, -4),
#                    control=list(smooth = F, maxit = maxit.SANN))
# # par <- optSANN$par
# t2 <- Sys.time()
# Rcpp_optSANN <- GenSA(fn = cpp_TVVARSS_ml, par = head(par, -4), lower = par.lower,
#                       upper = par.upper,
#                       X = t(X), U = matrix(), par_fixed = head(par.fixed, -4),
#                       control=list(smooth = F, maxit = maxit.SANN))
# t3 <- Sys.time()
# summary_output(R_optSANN, Rcpp_optSANN, c(t1, t2, t3))
# system2('say', 'R script finished')
# system2('terminal-notifier', "-message Done -title Script")
# save(Rcpp_optSANN, R_optSANN, file = 'GenSA.RData')

load('GenSA.RData')





# t1 <- Sys.time()
# Rcpp_optSANN2 <- GenSA(fn = cpp_TVVARSS_ml, par = head(par, -4), lower = par.lower,
#                        upper = par.upper,
#                        X = t(X), U = matrix(), par_fixed = head(par.fixed, -4),
#                        control=list(smooth = F, maxit = maxit.SANN))
# # par <- optSANN$par
# t2 <- Sys.time()
# R_optSANN2 <- GenSA(fn = TVVARSS.ml, par = head(par, -4), lower = par.lower,
#                     upper = par.upper,
#                     X = t(X), U = NULL, par.fixed = head(par.fixed, -4),
#                     control=list(smooth = F, maxit = maxit.SANN))
# t3 <- Sys.time()
# summary_output(Rcpp_optSANN2, R_optSANN2, c(t1, t2, t3))
# system2('say', 'R script finished')
# system2('terminal-notifier', "-message Done -title Script")
# save(Rcpp_optSANN2, R_optSANN2, file = 'GenSA2.RData')


# load('GenSA2.RData')

identical(Rcpp_optSANN, R_optSANN)
# identical(Rcpp_optSANN2, R_optSANN2)
# identical(R_optSANN, R_optSANN2)
# identical(Rcpp_optSANN, Rcpp_optSANN2)



# t1 <- Sys.time()
# set.seed(9)
# R_optSANN_ss <- GenSA(fn = TVVARSS.ml, par = head(par, -4), lower = par.lower,
#                       upper = par.upper,
#                       X = t(X), U = NULL, par.fixed = head(par.fixed, -4),
#                       control=list(smooth = F, maxit = 10))
# t2 <- Sys.time()
# set.seed(9)
# Rcpp_optSANN_ss <- GenSA(fn = cpp_TVVARSS_ml, par = head(par, -4), lower = par.lower,
#                          upper = par.upper,
#                          X = t(X), U = matrix(), par_fixed = head(par.fixed, -4),
#                          control=list(smooth = F, maxit = 10))
# t3 <- Sys.time()
# summary_output(R_optSANN_ss, Rcpp_optSANN_ss, c(t1, t2, t3))
# system2('say', 'R script finished')
# system2('terminal-notifier',"-message Done -title Script")
# save(Rcpp_optSANN_ss, R_optSANN_ss, file = 'GenSA_ss.RData')

length(R_optSANN_ss$par)
length(Rcpp_optSANN_ss$par)


set.seed(9)
gensa1 <- GenSA(fn = cpp_TVVARSS_ml, par = head(par, -4), lower = par.lower,
                upper = par.upper, X = t(X), U = matrix(), 
                par_fixed = head(par.fixed, -4),
                control=list(smooth = F, maxit = 10, max.call=10, verbose=T))
set.seed(9)
# gensa2 <- GenSA(fn = cpp_TVVARSS_ml, par = head(par, -4), lower = par.lower,
#                 upper = par.upper, X = t(X), U = matrix(), 
#                 par_fixed = head(par.fixed, -4),
#                 control=list(smooth = F, maxit = 10, max.call=1, verbose=T))
gensa2 <- GenSA(fn = TVVARSS.ml, par = head(par, -4), lower = par.lower,
                upper = par.upper, X = t(X), U = NULL, 
                par.fixed = head(par.fixed, -4),
                control=list(smooth = F, maxit = 10, max.call=10, verbose=T))

identical(gensa1$par, gensa2$par)
identical(gensa1$par, head(par, -4))
gensa1$par; gensa2$par


