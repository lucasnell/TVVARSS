## TVARSS.R

##   Time-Varying Autoregressive State Space model

## Copyright 2016 Anthony R. Ives


library(Rcpp)
library(RcppArmadillo)
library(GenSA)
library(minqa)

# File creating parameters necessary to test TVVARSS.ml and cpp_TVVARSS_ml
source('./initial_testing_files/TVVARSS.ml_pars.R')

source('TVVARSS.ml.R')




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



source('TVVARSS.ml.R')
sourceCpp('TVVARSS.cpp')
# system.time(replicate(100, cpp_TVVARSS_ml(par, t(X), t(U), par.fixed)))
# system.time(replicate(100, TVVARSS.ml(par = par, t(X), t(U), par.fixed)))

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

source('TVVARSS.ml.R')
sourceCpp('TVVARSS.cpp')

t1 <- Sys.time()
R_optim <- optim(fn = TVVARSS.ml, par = par, X = t(X), U = t(U), par.fixed = par.fixed,
                 method = "Nelder-Mead", control = list(maxit = 10^5))
# R_optim <- optim(fn = TVVARSS.ml, par = par, X = t(X), U = t(U), par.fixed = par.fixed,
#                  method = "BFGS", control = list(maxit = 10^5))
# R_optim <- bobyqa(fn = TVVARSS.ml, par = par, X = t(X), U = t(U), par.fixed = par.fixed,
#                   control = list(maxfun = 10^5))
t2 <- Sys.time()
Rcpp_optim <- optim(fn = cpp_TVVARSS_ml, par = par, X = t(X), U = t(U),
                    par_fixed = par.fixed, method = "Nelder-Mead",
                    control = list(maxit = 10^5))
# Rcpp_optim <- optim(fn = cpp_TVVARSS_ml, par = par, X = t(X), U = t(U),
#                     par_fixed = par.fixed, method = "BFGS",
#                     control = list(maxit = 10^5))
# Rcpp_optim <- bobyqa(fn = cpp_TVVARSS_ml, par = par, X = t(X), U = t(U),
#                      par_fixed = par.fixed, control = list(maxfun = 10^5))
t3 <- Sys.time()
summary_output(R_optim, Rcpp_optim, c(t1, t2, t3))
system2('say', 'R script finished')
# R version took 2.14 minutes
# Rcpp version took 0.08 minutes
# 
# Are they equal?
# Component “par”: Mean relative difference: 0.4850069
# Component “value”: Mean relative difference: 0.03822511
# Component “counts”: Mean relative difference: 0.2229251


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


# source('TVVARSS.ml.R')
# sourceCpp('TVVARSS.cpp')
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
# # save(Rcpp_optSANN_ss, R_optSANN_ss, file = 'GenSA_ss.RData')

length(R_optSANN_ss$par)
length(Rcpp_optSANN_ss$par)

t1 <- Sys.time()
set.seed(9)
gensa1 <- GenSA(fn = cpp_TVVARSS_ml, par = head(par, -4), lower = par.lower,
                upper = par.upper, X = t(X), U = matrix(), 
                par_fixed = head(par.fixed, -4),
                control=list(smooth = F, maxit = 10, verbose=T))
set.seed(9)
# gensa2 <- GenSA(fn = cpp_TVVARSS_ml, par = head(par, -4), lower = par.lower,
#                 upper = par.upper, X = t(X), U = matrix(), 
#                 par_fixed = head(par.fixed, -4),
#                 control=list(smooth = F, maxit = 10, max.call=1, verbose=T))
gensa2 <- GenSA(fn = TVVARSS.ml, par = head(par, -4), lower = par.lower,
                upper = par.upper, X = t(X), U = NULL, 
                par.fixed = head(par.fixed, -4),
                control=list(smooth = F, maxit = 10, verbose=T))
t2 <- Sys.time()
t2 - t1
system2("terminal-notifier", 
        c("-message", "Done", "-title", "Script", #"-sender", "org.rstudio.RStudio", 
          "-activate", "org.rstudio.RStudio"),
        stdout = NULL, stderr = NULL, wait = FALSE)
all.equal(gensa1$par, gensa2$par)
gensa1$par; gensa2$par






