# Testing R and C++ versions of TVVARSS.ml 

library(Rcpp)
library(RcppArmadillo)
library(GenSA)
library(minqa)

# File creating parameters necessary to test TVVARSS.ml and cpp_TVVARSS_ml
source('./initial_testing_files/TVVARSS.ml_pars.R')

# R and C++ versions of TVVARSS.ml (functions `TVVARSS.ml` and `cpp_TVVARSS_ml`, resp.)
source('./initial_testing_files/TVVARSS.ml.R')
sourceCpp('TVVARSS.cpp')




# ------------------------
# Testing initial output
# ------------------------

cpp_TVVARSS_ml(par = par, X = t(X), U = t(U), par_fixed = par.fixed)
TVVARSS.ml(par = par, X = t(X), U = t(U), par.fixed = par.fixed)

cpp_TVVARSS_ml(par = par, X = t(X), U = matrix(), par_fixed = par.fixed)
TVVARSS.ml(par = par, X = t(X), U = NULL, par.fixed = par.fixed)


# ------------------------
# Testing optim results
# ------------------------

# List of outputs
R_optim <- list()
Rcpp_optim <- list()
# Vector of times
tv.ml_t <- c(Sys.time())

# -------
# Using optim(..., method = 'Nelder-Mead')
# -------
tv.ml_t[1] <- Sys.time()
set.seed(1)
R_optim[['NM']] <- optim(fn = TVVARSS.ml, par = par, X = t(X), U = t(U), 
                         par.fixed = par.fixed, method = "Nelder-Mead", 
                         control = list(maxit = 10^5))
tv.ml_t[2] <- Sys.time()
set.seed(1)
Rcpp_optim[['NM']] <- optim(fn = cpp_TVVARSS_ml, par = par, X = t(X), U = t(U),
                       par_fixed = par.fixed, method = "Nelder-Mead",
                       control = list(maxit = 10^5))
tv.ml_t[3] <- Sys.time()


# -------
# Using optim(..., method = 'BFGS')
# -------
tv.ml_t[4] <- Sys.time()
set.seed(1)
R_optim[['BF']] <- optim(fn = TVVARSS.ml, par = par, X = t(X), U = t(U), 
                    par.fixed = par.fixed, method = "BFGS", 
                    control = list(maxit = 10^5))
tv.ml_t[5] <- Sys.time()
set.seed(1)
Rcpp_optim[['BF']] <- optim(fn = cpp_TVVARSS_ml, par = par, X = t(X), U = t(U),
                       par_fixed = par.fixed, method = "BFGS",
                       control = list(maxit = 10^5))
tv.ml_t[6] <- Sys.time()

# -------
# Using bobyqa
# -------
tv.ml_t[7] <- Sys.time()
set.seed(1)
R_optim[['bq']] <- bobyqa(fn = TVVARSS.ml, par = par, X = t(X), U = t(U), 
                     par.fixed = par.fixed, control = list(maxfun = 10^5))
tv.ml_t[8] <- Sys.time()
set.seed(1)
Rcpp_optim[['bq']] <- bobyqa(fn = cpp_TVVARSS_ml, par = par, X = t(X), U = t(U),
                        par_fixed = par.fixed, control = list(maxfun = 10^5))
tv.ml_t[9] <- Sys.time()







# ------------------------
# Testing GenSA
# ------------------------


# `maxit.SANN` reduced from 100 to 10 to reduce run time

# When U = NULL (U = matrix() for C++ version), you must do head(<par*>, -4) , 
# where <par*> is par.lower, par.upper, par, and par.fixed

tv.ml_t[10] <- Sys.time()
set.seed(1)
R_GenSA <- GenSA(fn = TVVARSS.ml, par = par, lower = par.lower, upper = par.upper,
                   X = t(X), U = t(U), par.fixed = par.fixed,
                   control=list(smooth = F, maxit = maxit.SANN))
tv.ml_t[11] <- Sys.time()
set.seed(1)
Rcpp_GenSA <- GenSA(fn = cpp_TVVARSS_ml, par = par, lower = par.lower, upper = par.upper,
                      X = t(X), U = t(U), par_fixed = par.fixed,
                      control=list(smooth = F, maxit = maxit.SANN))
tv.ml_t[12] <- Sys.time()


# Total time taken
message(sprintf('Total time taken: %s minutes',
                round(as.numeric(tail(tv.ml_t, 1) - tv.ml_t[1], units = 'mins'), 2)))

save(R_optim, Rcpp_optim, R_GenSA, Rcpp_GenSA, tv.ml_t, summary_output,
     file = './initial_testing_files/TVVARSS.ml_testing.RData')






# ------------------------
# Summarizing output
# ------------------------

# Function to return summary output for 2 objects and 3 times
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


summary_output(R_optim[['NM']], Rcpp_optim[['NM']], tv.ml_t[1:3])
summary_output(R_optim[['BF']], Rcpp_optim[['BF']], tv.ml_t[4:6])
summary_output(R_optim[['bq']], Rcpp_optim[['bq']], tv.ml_t[7:9])

summary_output(R_GenSA, Rcpp_GenSA, tv.ml_t[10:12])








