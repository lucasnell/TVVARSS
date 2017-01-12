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
set.seed(1)
zNM <- TVVARSS(X, Tsamplefract = 0.9, method="Nelder-Mead", show.fig = FALSE,
               annealing = FALSE)
               # annealing = TRUE)

# Takes ~1.5 min with annealing, ~2 min without
set.seed(1)
cpp_zNM <- cpp_TVVARSS(X, Tsamplefract = 0.9, method="Nelder-Mead", show.fig = FALSE,
                       annealing = FALSE)
                       # annealing = TRUE)
system2('say', 'R script finished')
system2('terminal-notifier',"-message Done -title Script")
all.equal(zNM, cpp_zNM)

# With annealing...
# > all.equal(zNM, cpp_zNM)
# [1] "Component “se”: Mean relative difference: 7.369233e-06"             
# [2] "Component “su”: Mean relative difference: 3.267949e-05"             
# [3] "Component “Sb0”: Mean relative difference: 3.174057e-06"            
# [4] "Component “Sb”: Mean relative difference: 7.188581e-06"             
# [5] "Component “B0”: Mean relative difference: 1.60482e-05"              
# [6] "Component “B”: Mean relative difference: 7.586266e-06"              
# [7] "Component “logLik”: Mean relative difference: 2.343661e-06"         
# [8] "Component “AIC”: Mean relative difference: 2.228445e-06"            
# [9] "Component “B0.fitted”: Mean relative difference: 0.002025549"       
# [10] "Component “B.fitted”: Mean relative difference: 0.002272341"        
# [11] "Component “X.fitted”: Mean relative difference: 0.000153545"        
# [12] "Component “PP.fitted”: Mean relative difference: 0.001931707"       
# [13] "Component “eigen.fitted”: Mean relative Mod difference: 0.002331226"
# [14] "Component “opt.par”: Mean relative difference: 9.404639e-06"

# Without annealing...
# > all.equal(zNM, cpp_zNM)
# [1] "Component “se”: Mean relative difference: 0.08068179"            
# [2] "Component “su”: Mean relative difference: 1.238354"              
# [3] "Component “Sb0”: Mean relative difference: 1.615259"             
# [4] "Component “Sb”: Mean relative difference: 0.5464784"             
# [5] "Component “B0”: Mean relative difference: 0.186959"              
# [6] "Component “B”: Mean relative difference: 0.1338447"              
# [7] "Component “logLik”: Mean relative difference: 0.005695631"       
# [8] "Component “AIC”: Mean relative difference: 0.005253712"          
# [9] "Component “B0.fitted”: Mean relative difference: 1.199959"       
# [10] "Component “B.fitted”: Mean relative difference: 0.2313852"       
# [11] "Component “X.fitted”: Mean relative difference: 0.05634668"      
# [12] "Component “PP.fitted”: Mean relative difference: 1.284283"       
# [13] "Component “eigen.fitted”: Mean relative Mod difference: 0.151117"
# [14] "Component “opt.par”: Mean relative difference: 0.2581309"
# 
# 
# save(zNM, cpp_zNM, file = 'TVVARSS.RData')
# load('TVVARSS.RData')

save(zNM, cpp_zNM, file = 'TVVARSS_no_annealing.RData')
# load('TVVARSS_no_annealing.RData')


all.equal(cpp_zNM, cpp_zNM2)



summary(cpp_zNM2)
summary(zNM)

all.equal(zNM, cpp_zNM)








cppFunction(depends = 'RcppArmadillo',
            code = 
'bool cx_present(arma::cx_mat M){
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
}')



cx_present(matrix(c(1:9,1i), 2))
cx_present(matrix(1:10, 2))





# library(Rcpp)
# library(RcppArmadillo)
# library(GenSA)

sourceCpp('TVVARSS.cpp')




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
