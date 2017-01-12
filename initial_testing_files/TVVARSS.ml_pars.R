# Setting up parameters necessary to test TVVARSS.ml (both R and C++ versions)

set.seed(79)
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


set.seed(8)
U <- cbind(U, U * rnorm(100))
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












# --------------------------
# =================

# For GenSA testing...

# =================
# --------------------------

maxit.SANN <- 10 # 10^2

par.upper <- c(array(1, dim=n), array(2, dim=n^2), array(10, dim=n), array(10, dim=n),
               array(10, dim=n*(n+1)), array(100, dim=n*nu))
par.lower <- c(array(-1, dim=n), array(-2, dim=n^2), array(0, dim=n), array(0, dim=n),
               array(0, dim=n*(n+1)), array(-100, dim=n*nu))









# =====================================================================================
# =====================================================================================
# =====================================================================================
# =====================================================================================

#       Intermediate objects within TVVARSS.ml function
#       (not necessary for initial setup, but left here in case they're useful for 
#        for troubleshooting)

# =====================================================================================
# =====================================================================================
# =====================================================================================
# =====================================================================================


# par.full <- par.fixed
# par.full[is.na(par.fixed)] <- par
# 
# # Se
# diag(par.full[(n+n^2+1):(n+n^2+n)]^2)
# # Su
# diag(par.full[(n+n^2+n+1):(n+n^2+n+n)]^2)
# # Sb
# diag(par.full[(n+n^2+n+n+1):(n+n^2+n+n+n*(n+1))]^2)
# 
# # C
# C = matrix(par.full[(n+n^2+n+n+n*(n+1)+1):(n+n^2+n+n+n*(n+1)+nu*n)], nrow=n, ncol=nu, byrow = TRUE); C
# 
# # S
# # as.matrix(rbind(cbind(Se, matrix(0, n, n*(n+1))), 
# #                 cbind(matrix(0, n*(n+1), n), Sb)))
# S = as.matrix(rbind(cbind(diag(par.full[(n+n^2+1):(n+n^2+n)]^2), 
#                           matrix(0, n, n*(n+1))),
#                     cbind(matrix(0, n*(n+1), n), 
#                           diag(par.full[(n+n^2+n+n+1):(n+n^2+n+n+n*(n+1))]^2))))
# 
# # Z
# Z = as.matrix(cbind(diag(n), matrix(0, n, n*(n+1)))); Z
# 
# # x
# t(X)[,1]
# 
# # Check "If the initial parameter values imply a stationary distribution..."
# max(abs(eigen(matrix(par.full[(n+1):(n+n^2)], nrow=n, ncol=n, byrow = TRUE))$values))
# 
# # PP
# matrix(solve(diag(n*n)-kronecker(matrix(par.full[(n+1):(n+n^2)], nrow=n, ncol=n, byrow = TRUE),matrix(par.full[(n+1):(n+n^2)], nrow=n, ncol=n, byrow = TRUE))) %*% matrix(diag(par.full[(n+n^2+1):(n+n^2+n)]^2), nrow = n*n), n, n)
# 
# # PP after manipulation
# PP0 = as.matrix(rbind(cbind(matrix(solve(diag(n*n)-kronecker(matrix(par.full[(n+1):(n+n^2)], nrow=n, ncol=n, byrow = TRUE),matrix(par.full[(n+1):(n+n^2)], nrow=n, ncol=n, byrow = TRUE))) %*% matrix(diag(par.full[(n+n^2+1):(n+n^2+n)]^2), nrow = n*n), n, n), matrix(0, n, n*(n+1))), 
#                       cbind(matrix(0, n*(n+1), n), diag(par.full[(n+n^2+n+n+1):(n+n^2+n+n+n*(n+1))]^2))))
# 
# 
# # BB
# B12 <- diag(n) - matrix(par.full[(n+1):(n+n^2)], nrow=n, ncol=n, byrow = TRUE)
# B13 <- kronecker(t(t(X)[,1]-matrix(par.full[1:n], nrow=n, ncol=1)),diag(n))
# 
# BB = as.matrix(rbind(cbind(matrix(par.full[(n+1):(n+n^2)], nrow=n, ncol=n, byrow = TRUE), B12, B13), 
#                      cbind(matrix(0, n, n), diag(n), matrix(0, n, n^2)), 
#                      cbind(matrix(0, n^2, 2*n), diag(n^2))))
# 
# 
# # PP after BB creation
# PP = BB %*% PP0 %*% t(BB) + S
# 
# # x after first iteration of for loop
# # B0 + B %*% (x-B0) + C * U[t]
# # # B0 + B %*% (x-B0) + C %*% U[,t]
# x0 = matrix(par.full[1:n], nrow=n, ncol=1) + matrix(par.full[(n+1):(n+n^2)], nrow=n, ncol=n, byrow = TRUE) %*% (as.matrix(t(X)[,1]) - matrix(par.full[1:n], nrow=n, ncol=1)) + C %*% as.matrix(t(U)[,2])
# 
# B0_0 <- matrix(par.full[1:n], nrow=n, ncol=1)
# B_0 <- matrix(par.full[(n+1):(n+n^2)], nrow=n, ncol=n, byrow = TRUE)
# 
# FF <- Z %*% PP %*% t(Z) + diag(par.full[(n+n^2+n+1):(n+n^2+n+n)]^2)
# invF <- solve(FF)
# 
# y0 <- matrix(c(x0, B0_0, t(B_0)), ncol=1)
# v <- as.matrix(t(X)[,2]) - Z %*% y0
# y <- y0 + PP %*% t(Z) %*% invF %*% v

