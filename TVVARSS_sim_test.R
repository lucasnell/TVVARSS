source('./initial_testing_files/TVVARSS_forC_17Dec16.R')
source('TVVARSS_cpp.R')

# simulation test of TVVARSS

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

par(mfrow=c(1,1))
matplot(X, typ="l")

# quartz()
# Takes ~30 min with annealing, ~1.5 min without
# set.seed(1)
# zNM <- TVVARSS(X, Tsamplefract = 0.9, method="Nelder-Mead", show.fig = FALSE,
#                #annealing = FALSE)
#                annealing = TRUE)
# 
# # Takes ~1.5 min with annealing, ~2 min without
# set.seed(1)
cpp_zNM <- cpp_TVVARSS(X, Tsamplefract = 0.9, method="Nelder-Mead", show.fig = FALSE,
                       #annealing = FALSE)
                       annealing = TRUE)
# system2('say', 'R script finished')
# system2('terminal-notifier',"-message Done -title Script")
# 
# 
# 
# save(zNM, cpp_zNM, file = 'TVVARSS.RData')
load('TVVARSS.RData')

# save(zNM, cpp_zNM, file = 'TVVARSS_no_annealing.RData')
# load('TVVARSS_no_annealing.RData')


all.equal(zNM, cpp_zNM)








cppFunction(depends = 'RcppArmadillo',
            code = 
'arma::mat mat_byrow(arma::vec x, arma::uword nrow, arma::uword ncol) {
    if(x.n_elem != (nrow * ncol)){
        arma::mat M = arma::zeros(0);
        return(M);
    }
    arma::mat M = arma::zeros(nrow, ncol);
    arma::mat row_iter;
    arma::uword j = 0;
    for(arma::uword i=0; i<nrow; ++i) {
        row_iter = x.subvec(j, (j + ncol - 1));
        row_iter.reshape(1,ncol);
        M.row(i) = row_iter;
        j += ncol;
    }
    return(M);
    }')




# Looking at how C++ manages complex numbers
# library(Rcpp)
# library(RcppArmadillo)
# cppFunction(depends = 'RcppArmadillo', 
#             code = 
# "arma::cx_mat imag_test(arma::mat X){
#     arma::cx_mat M = arma::randu<arma::cx_mat>(1,2);
#     M(0,1) = X(0,0);
#     arma::mat I = arma::imag(M);
#     // arma::mat C = M;
#     return(M);
# }")
# 
# set.seed(9); imag_test(matrix(4+1i^9))




# quartz()
# zBFGS <- TVVARSS(X, Tsamplefract = .9, method="BFGS", show.fig = T, annealing = T)
# summary(zBFGS)
# 
# quartz()
# zbobyqa <- TVVARSS(X, Tsamplefract = .9, method="bobyqa", show.fig = T, annealing = T)
# summary(zbobyqa)
# 
# quartz()
# zNM <- TVVARSS(X, Sb0.fixed = matrix(.01, nrow = n, ncol = 1), Tsamplefract = .9, method="Nelder-Mead", show.fig = T, annealing = T)
# summary(zNM)
# 
# quartz()
# zBFGS <- TVVARSS(X, Sb0.fixed = matrix(.01, nrow = n, ncol = 1), Tsamplefract = .9, method="BFGS", show.fig = T, annealing = T)
# summary(zBFGS)
# 
# quartz()
# zbobyqa <- TVVARSS(X, Sb0.fixed = matrix(.01, nrow = n, ncol = 1), Tsamplefract = .9, method="bobyqa", show.fig = T, annealing = T)
# summary(zbobyqa)
