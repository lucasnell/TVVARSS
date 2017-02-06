# Simulation test of TVVARSS final output
# Output used in ../TVVARSS_testing.Rmd


source('./initial_testing_files/TVVARSS_forC_17Dec16.R')
source('TVVARSS_cpp.R')


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


run_method <- function(method_type, language, seed = 1){
    # Standardizing language
    lang <- c('R', 'cpp')[grepl(language, c('R', 'cpp'), ignore.case = TRUE)]
    if (length(lang) == 0){
        stop(paste('Function takes only "R" or "cpp" for the language argument.',
                   '(Capitalization does not matter.)'))
    }
    # Error function
    e_fun <- function(e, ta) {
        mess <- paste('Error on', lang, 'version of TVVARSS with method =', 
                      method_type, 'and annealing =', ta, '\n')
        warning(mess)
        return(NULL)
    }
    if (lang == 'R'){
        z <- lapply(
            c(TRUE, FALSE), 
            function(to_anneal) {
                set.seed(seed)
                tryCatch(TVVARSS(X, Tsamplefract = 0.9, method = method_type, 
                                 show.fig = FALSE, annealing = to_anneal), 
                         error = function(e) e_fun(e, to_anneal))})
    } else if (lang == 'cpp'){
        z <- lapply(
            c(TRUE, FALSE), 
            function(to_anneal) {
                set.seed(seed)
                tryCatch(cpp_TVVARSS(X, Tsamplefract = 0.9, method = method_type, 
                                     show.fig = FALSE, annealing = to_anneal), 
                         error = function(e) e_fun(e, to_anneal))})
    }
    names(z) <- c('anneal', 'no_anneal')
    return(z)
}


tv_t <- c(Sys.time())
cpp_out <- lapply(c('Nelder-Mead', 'BFGS', 'bobyqa'), 
                  function(x){run_method(x, 'cpp')})
tv_t[2] <- Sys.time()
R_out <- lapply(c('Nelder-Mead', 'BFGS', 'bobyqa'), 
                function(x){run_method(x, 'R')})
tv_t[3] <- Sys.time()


names(cpp_out) <- c('Nelder-Mead', 'BFGS', 'bobyqa')
names(R_out) <- c('Nelder-Mead', 'BFGS', 'bobyqa')

# X is saved so I can show output from failed run of TVVARSS
save(cpp_out, R_out, tv_t, X, file = './initial_testing_files/TVVARSS_testing.RData')

