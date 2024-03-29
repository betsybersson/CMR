source("/Users/betsybersson/Documents/betsy_R_functions.R")
source("functions.R")
library(parallel)

# sample size
n = 20
# feature dimension
p = 8
# latent dimension
k = 4

# X = matrix(rnorm(n*p),ncol = p)

# true cov
cov.true = cor.mat(p,.2)

# simulate fake data
set.seed(12)
Y = rmatnorm(matrix(0,ncol = p, nrow = n),V=cov.true)
set.seed(Sys.time())

# run GS
out = CMR_GS(Y,
             k = 4,
             S = 6100,
             burnin = 100)

# plot trace plots of cov
plot.trace.ind = sample(1:(p^2),5)
matplot(out$cov[,plot.trace.ind],type="l")

cov.hat.stein = qr.solve(matrix(colMeans(out$cov.inv),ncol = p)) ## stein estimator
cov.hat.frob = matrix(colMeans(out$cov),ncol = p) ## frobenius estimator

loss_stein(cov.hat.stein,cov.true)
