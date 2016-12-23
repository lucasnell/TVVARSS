source('TVVARSS_forC_17Dec16.R')
source('TVVARSS_cpp.R')

# simulation test of TVVARSS

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
# zNM <- TVVARSS(X, Tsamplefract = .9, method="Nelder-Mead", show.fig = T, annealing = T)
# 12.5 minutes


t1 <- Sys.time()
cpp_zNM <- cpp_TVVARSS(X, Tsamplefract = .9, method="Nelder-Mead", show.fig = F, annealing = T)
t2 <- Sys.time()
t2 - t1
# summary(zNM)
# summary(cpp_zNM)












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
