##########################################################
## Matrix Normal Density
# X: nxp matrix you want to evaluate density at
# M: nxp mean of normal density
# U: nxn row covariance matrix
# V: pxp column covariance matrix
##########################################################
dmatnorm = function(X,
                    M = matrix(0,ncol=ncol(X),nrow=nrow(X)),
                    U = eye(nrow(X)),
                    V = eye(ncol(X)),
                    U.inv = qr.solve(U),
                    V.inv = qr.solve(V),
                    if_log = FALSE){
  # function for density of X_{ n\times p} \sim N(M, V kron U)
  N = nrow(X)
  P = ncol(X)
  out = (2*pi)^(-N*P/2) * det(U.inv)^(P/2) * det(V.inv)^(N/2)*
    exp(-mat_trace( V.inv %*% t(X-M) %*% U.inv %*% (X-M) )/2)
  if (if_log == TRUE){
    out = -N*P*log(2*pi)/2 - P*determinant(U)$mod[1]/2 - N*determinant(V)$mod[1]/2 - 
      mat_trace( V.inv %*% t(X-M) %*% U.inv %*% (X-M) )/2
  }
  return(out)
}

##########################################################
## Random Multivariate or Matrix-variate Normal Sample
# M: nxp mean of normal density
# U: nxn row covariance matrix
# V: pxp column covariance matrix
##########################################################
rmatnorm = function(M,U = eye(nrow(M)),V = eye(ncol(M))) {
  # function for density of Y_{ n\times p} \sim N(M, V kron U)
  N = nrow(M)
  P = ncol(M)
  
  E = matrix(rnorm(N*P),ncol=P)
  
  out = M + t(chol(U)) %*% E %*% chol(V)
  
  return(out)
}
rmvnorm = function(M,V) {
  # function for density of Y_{px1} \sim N_p(M, V) (ind rows)
  P = length(M)
  
  E = rnorm(P)
  
  out = c(M) + t(chol(V)) %*% E  
  
  return(c(out))
}
##########################################################
## Wishart Density
# X: pxp matrix you want to evaluate density at
# M: pxp scale matrix
# nu: scaler degrees of freedom
##########################################################
dwish = function(X,M,nu, if_log = FALSE) {
  # function for density of X \sim W(M, nu)
  P = ncol(X)
  gamma.p = pi^(P*(P-1)/4) * prod(gamma( (nu+1-c(1:P))/2 ))
  out = 1/(2^(nu*P/2) * gamma.p * det(M)^(nu/2) ) * det(X)^((nu-P-1)/2) * 
    exp(mat_trace(-qr.solve(M)%*%X/2))
  if ( if_log == TRUE ){
    out = (-1) * ((nu*P/2)*log(2) + (P*(P-1)/4)*log(pi) + sum(lgamma( (nu+1-c(1:P))/2 )) +
                    (nu/2)*determinant(M)$mod) + ((nu-P-1)/2)*determinant(X)$mod +
      mat_trace(-qr.solve(M)%*%X/2)
    out = out[1]
  }
  return(out)
}
##########################################################
## Multivariate Normal Density
# X: nx1 matrix you want to evaluate density at
# M: nx1 mean of normal density
# U: nxn row covariance matrix
# dmatnorm could work, but this may be faster because one less det() call
##########################################################
dmvnorm = function(X,M,U, if_log = FALSE) {
  # function for density of X_{ n\times p} \sim N(0, V kron U)
  N = nrow(X)
  P = 1
  if (is.null(N)){N = 1}
  out = (2*pi)^(-N*P/2) * det(U)^(-P/2) *
    exp(-( t(X-M) %*% qr.solve(U) %*% (X-M) )/2)
  if (if_log == TRUE){
    out = -N*P*log(2*pi)/2 - P*log(det(U))/2  - 
      ( t(X-M) %*% qr.solve(U) %*% (X-M) )/2
  }
  return(c(out))
}
##########################################################
## Inverse Wishart Density
# x: pxp matrix you want to evaluate density at
# nu: degrees of freedom
# s: scale matrix
##########################################################
dinvwish = function(x,nu,s,log.out=F){
  ## do not invert s
  
  p = dim(s)[1]
  gprod = prod(sapply(1:p,function(kk)gamma((nu+1-kk)/2)))
  
  etr.term = exp(-sum(diag(qr.solve(s)%*%qr.solve(x)))/2)
  out = 1/(2^(nu*p/2)*pi^(choose(p,2)/2)*gprod*det(s)^(nu/2))*
    det(x)^(-(nu+p+1)/2)*etr.term

  if(log.out==T){
    gprod = sum(sapply(1:p,function(kk)lgamma((nu+1-kk)/2)))
    etr.term = -sum(diag(qr.solve(s)%*%qr.solve(x)))/2
    out = -(nu*p/2) * log(2) -
      (choose(p,2)/2) * log(pi) -
      gprod -
      (nu/2) * determinant(s)$mod -
      (nu+p+1)/2 * determinant(x)$mod +
      etr.term
    out = out[1]
  }
  
  return(out)
  
}
##########################################################
## Inverse Wishart Density- prop to nu
# x: pxp matrix you want to evaluate density at
# nu: degrees of freedom
# s: scale matrix
##########################################################
dinvwish.proptonu = function(x,nu,s,log.out=F){
  ## do not invert s
  
  p = dim(s)[1]

  ### log
  if ( log.out == T){
    gprod = sum(log(sapply(1:p,function(kk)gamma((nu+1-kk)/2))))
    
    etr.term = (-sum(diag(qr.solve(s)%*%solve(x)))/2)
    out = (-nu*p/2)*log(2) - gprod - (nu/2)*determinant(s)$mod +
      (-nu/2)*determinant(x)$mod + etr.term
    # out = exp(out)
    # if(out==0){
    #   out = 1e-15
    # }
  } else{
    gprod = prod(sapply(1:p,function(kk)gamma((nu+1-kk)/2)))
    
    etr.term = exp(-sum(diag(qr.solve(s)%*%solve(x)))/2)
    out = 1/(2^(nu*p/2)*gprod*det(s)^(nu/2))*
      det(x)^(-nu/2)*etr.term
  }
  
  
  return(out)
  
}
##########################################################
## Box's M multivariate test for homoskedasticity
# Y: nxp response variable
# group: vector of length n, factor, group membership of each entry of Y
# code based off of Friendly's heplot boxM code
##########################################################
boxM_robust = function(Y,group){
  
  p <- ncol(Y)
  nlev <- nlevels(group)
  lev <- levels(group)
  dfs <- tapply(group, group, length) - 1
  if (any(dfs < p)) 
    warning("there are one or more levels with less observations than variables!")
  mats <- aux <- list()
  for(i in 1:nlev) {
    mats[[i]] <- cov(Y[group == lev[i], ])
    aux[[i]] <- mats[[i]] * dfs[i]
  }
  names(mats) <- lev
  pooled <- Reduce("+", aux) / sum(dfs)
  logdet <- log(unlist(lapply(mats, det)))
  minus2logM <- sum(dfs) * log(det(pooled)) - sum(logdet * dfs)
  sum1 <- sum(1 / dfs) 
  sum2 <- sum((1 / dfs)^2)
  A1 <- (((2 * p^2) + (3 * p) - 1) / (6 * (p + 1) *
                                        (nlev - 1))) * (sum1 - (1 / sum(dfs)))
  dfchi <- (choose(p, 2) + p) * (nlev - 1)
  
  if ( any(dfs<19) ){
    A2 <- (p-1) * (p+2) / (6 * (nlev - 1) ) * (sum2 - (1 / sum(dfs)^2))
    if (A2-A1^2 > 0){
      dfF = (dfchi + 2)/(A2-A1^2)
      b = dfchi/(1-A1-dfchi/dfF)
      X2 <- minus2logM / b
      pval <- pf(X2, dfchi, dfF, lower.tail=FALSE)
    } else if (A2-A1^2<0){
      dfF = (dfchi + 2)/(A1^2-A2)
      b = dfF / (1-A1+2/dfF)
      X2 <- dfF * minus2logM / (dfchi * (b - minus2logM))
      pval <- pf(X2, dfchi, dfF, lower.tail=FALSE)
    } else{
      warning("error")
    }
  } else{
    X2 <- minus2logM * (1 - A1)
    pval <- pchisq(X2, dfchi, lower.tail = FALSE)
  }
  
  return(pval)
  
}
##########################################################
## Identity
# N: length of diagonal
##########################################################
eye = function(N){
  diag(rep(1,N))
}
##########################################################
## Trace of matrix
# X: square matrix
##########################################################
mat_trace = function(X){
  sum(diag(X))
}
##########################################################
## Covariance Stein Loss
# A: Estimator
# B: Truth
##########################################################
loss_stein = function(A,B){
  helper = A %*% qr.solve(B)
  ld_help = determinant(helper,logarithm = TRUE)$modulus[1]
  mat_trace(helper) - ld_help - nrow(A)
}
##########################################################
## Covariance Quadratic Loss
# A: Estimator
# B: Truth
##########################################################
loss_quad = function(A,B){
  helper = A %*% qr.solve(B) - eye(nrow(A))
  mat_trace(helper %*% helper) 
}
##########################################################
## Covariance Modified Quadratic Loss
# A: Estimator
# B: Truth
##########################################################
loss_mod_quad = function(A,B){
  helper = A - B
  mat_trace(helper %*% helper) /nrow(A)
}
##########################################################
## Covariance Modified Squared Error Loss
# A: Estimator
# B: Truth
##########################################################
loss_mod_sq = function(A,B){
  sum((A[upper.tri(A)] - B[upper.tri(B)])^2) + 
    sum((diag(A) - diag(B))^2)
}
##########################################################
## Squared Error Loss
# A: Estimator
# B: Truth
##########################################################
loss_sq = function(A,B){
  mat_trace(crossprod(A-B))
}
##########################################################
## Correlation Loss
# A: Estimator
# B: Truth
##########################################################
loss_corr = function(A,B){
  1 - mat_trace(A%*%B)/norm(A,type="F")/norm(B,type="F")
}
##########################################################
## PH z scores transformation- from copula.pdf notes for STA 832
# y: vec to have zcores found of
##########################################################
zscores <- function(y,ties.method="average"){
  u<-rank(y,na.last="keep",ties.method=ties.method)/(sum(!is.na(y))+1)
  z<-qnorm(u)
  names(z)<-names(y)
  m<-dim(y)
  if(length(m)==2){
    z<-matrix(z,nrow=m[1],ncol=m[2])
    dimnames(z)<-dimnames(y)
  }
  if(length(m)>=3){
    z<-array(z,dim=m)
    dimnames(z)<-dimnames(y)
  }
  z
}
##########################################################
## K Nearest Neighbors- obtained from a distance matrix
# D: distance metric
##########################################################
k_nn <- function(D,k){
  out = apply(D,1,function(j)order(j)[2:(k+1)]) %>% t()
  if (length(rownames(D))>0){
    rownames(out) = rownames(D)
  }
  return(out)
}
##########################################################
## Make a matrix symmetric: copy the upper triangle to the lower triangle
# m: matrix with desired upper triangle
##########################################################
makeSymm <- function(m) {
  m[lower.tri(m)] <- t(m)[lower.tri(m)]
  return(m)
}
##########################################################
## get mode of vector
# v: vector
##########################################################
getmode <- function(v) {
  uniqv <- unique(v)
  uniqv[which.max(tabulate(match(v, uniqv)))]
}
##########################################################
## noncentral t density evaluated at a point X
# X: point to evaluate density at
# nu: number of degrees of freedom
# mu: mean of t distribution
# sig2: variance of t distribution
## density is written down on STA 721 midterms in my classes
##########################################################
noncent_dt = function(X,nu,mu,sig2,propto = TRUE){
  # density at X for t distribution with 
  # nu degrees of freedom
  # mean mu
  # variance sig2
  # if propto = true, then return density up to normalizing constant
  
  out = gamma((nu+1)/2) * (nu * sig2 * pi)^(-1/2) / gamma(nu/2) *
    (1 + (X - mu)^2/(sig2 * nu))^(-(nu+1)/2) 
  
  if (propto == TRUE){
    out = out/(gamma((nu+1)/2)/ gamma(nu/2) /sqrt(nu*pi))
  }
  
  return(out)
}
##########################################################
## shifted poisson density
# x: point to evaluate density at
# lambda: poisson parameter
# C: shift parameter
##########################################################
dpois_shift = function(x,lambda,C,log.ind=FALSE){
  
  out = exp(-lambda)*lambda^(x-C)/factorial(x-C)
  
  if ( log.ind == TRUE) {
    out = (-lambda) + (x-C) * log(lambda) - log(factorial(x-C))
  }
  
  return(out)
  
}
##########################################################
## sample from truncated normal
# n: number of samples
# mu: mean
# sig2: variance
# a: lower bound
# b: upper bound
##########################################################
rtnorm = function(n,a=-Inf,b=Inf,mu=0,sig2=1){
  
  sig = sqrt(sig2)
  
  a.star = pnorm(a,mu,sig)
  b.star = pnorm(b,mu,sig)
  
  X = runif(n,a.star,b.star)
  
  out = qnorm(X,mu,sig)
  
  return(out)
  
}
##########################################################
## col and row medians of a matrix
# X: matrix
##########################################################
colMedians = function(X){
  
  apply(X,2,median)
  
}
rowMedians = function(X){
  
  apply(X,1,median)
  
}
##########################################################
## Matrix T Density
# X: matrix to evaluate density at (nxp)
# nu: degrees of freedom
# M: location (nxp matrix)
# Omega: pxp column covariance matrix
# Sigma: nxn row covariance matrix
##########################################################
dmatT = function(X,nu,M,Omega,Sigma,
                 Omega.inv = qr.solve(Omega),
                 Sigma.inv = qr.solve(Sigma),
                 log=FALSE){
  
  p = ncol(X)
  n = nrow(X)
  
  # help with huge numbers
  gamma.div = gamma.mv(p,(nu+n+p-1)/2,T) - gamma.mv(p,(nu+p-1)/2,T)

  if (log == TRUE){
    main.det = determinant(eye(n) + Sigma.inv %*% (X-M) %*% Omega.inv %*% t(X-M))
    
    out = gamma.div - (n*p/2) * pi + 
      (n/2) * determinant(Omega.inv) + 
      (p/2) * determinant(Sigma.inv) + 
      (-(nu+n+p-1)/2) * main.det
  } else {
    main.det = det(eye(n) + Sigma.inv %*% (X-M) %*% Omega.inv %*% t(X-M))
    out = exp(gamma.div) / pi^(n*p/2) * 
      det(Omega.inv)^(n/2) * det(Sigma.inv)^(p/2) * 
      main.det^(-(nu+n+p-1)/2)
  }
  
  
  return(out)
  
}
##########################################################
## Multivariate T Density
# X: matrix to evaluate density at (nxp)
# nu: degrees of freedom
# M: location (nxp matrix)
# Omega: pxp column covariance matrix
# Sigma: eye(n) row covariance matrix
##########################################################
dmvT = function(X,nu,M,Omega,
                 Omega.inv = qr.solve(Omega),
                 log=FALSE){
  
  p = ncol(X)
  n = nrow(X)
  
  # help with huge numbers
  gamma.div = gamma.mv(p,(nu+n+p-1)/2,T) - gamma.mv(p,(nu+p-1)/2,T)
  
  if (log == TRUE){
    main.det = determinant(eye(n) + (X-M) %*% Omega.inv %*% t(X-M))$mod
    
    out = gamma.div - (n*p/2) * pi + 
      (n/2) * determinant(Omega.inv)$mod + 
      (-(nu+n+p-1)/2) * main.det
  } else {
    main.det = det(eye(n) + (X-M) %*% Omega.inv %*% t(X-M))
    out = exp(gamma.div) / pi^(n*p/2) * 
      det(Omega.inv)^(n/2) * 
      main.det^(-(nu+n+p-1)/2)
  }
  
  
  return(out)
  
}
##########################################################
## Multivariate Gamma Function
# p: dimension
# a: object
# output: Gamma_p(a)
##########################################################
gamma.mv = function(p,a,log.out = FALSE){
  out = pi^(p*(p-1)/4) * prod(gamma(a+(1-c(1:p))/2))
  if (log.out == TRUE){
    out = (p*(p-1)/4) * log(pi) + sum(lgamma(a+(1-c(1:p))/2))
  }
  return(out)
}
##########################################################
## Random covariance matrix
# p: dimension
# sd: sd of the underlying normal matrix
##########################################################
rcovmat = function(p,sd = 1){
  
  out.temp = matrix(rnorm(p*p,sd = sd),ncol = p)
  out = out.temp %*% t(out.temp)
  
  return(out)
}
##########################################################
## Zero matrix
# p1:  number of rows 
# p2: numbr of columns (optional, return matrix if empty)
##########################################################
zeros = function(p1,p2 = NA){
  
  if (is.na(p2)){
    p2 = p1
  }
  
  out = matrix(0,nrow = p1, ncol = p2)
  
  return(out)
}
##########################################################
## Morans I metric
# XX: vector of response variable
# WW: row standardized adjacency matrix
# output: morans I value: near 0 means no/little correlation
##########################################################
moransI = function(XX,WW){

  NN = length(XX)
  
  mean.dif = (XX - mean(XX))
  
  num = c(t(mean.dif) %*% WW %*% (mean.dif))
  
  NN/sum(WW) * num/sum((mean.dif)^2)
  
}
##########################################################
## Geary's C metric
# XX: vector of response variable
# WW: row standardized adjacency matrix
# output: geary's C value: near 1 means no/little autocorrelation
# difference from morans I: more sensitive to local autocorrelation
##########################################################
gearysC = function(XX,WW){
  ## near 1 means no/little correlation
  
  NN = length(XX)
  
  num = 0
  for (ii in 1:NN){
    for (jj in 1:NN){
      num = num + WW[ii,jj] * (XX[ii]-XX[jj])^2
    }
  }
  
  (NN-1)/(2*sum(WW)) * num/sum((XX - mean(XX))^2)
  
}
##########################################################
## Symmetric proposal on domain (0,1) for Metropolis algorithm
# given in Hoff page 190
# center: the previous value, where the proposal uniform distribution is centered at
# delta: neighborhood around "center" to draw uniformly from
# out: draw from symmetric proposal on domain (0,1)
##########################################################
MH_sym_proposal_01 = function(center, delta){
  out = runif(1,center-delta,center+delta)
  if (out<0){
    out = abs(out)
  } else if (out>1){
    out = 2-out
  }
  return(out)
}
##########################################################
## Box M Test for Homogeneity of  (pxp) Covariances
# updated from install_github("friendly/heplots"). Forked and edited on my github
# Y: matrix of all data (N x p)
# group: vector (N) of group membership for each row in Y 
# f.indicator: if TRUE, use F test over chisq test. Set TRUE if sample size is small.
# out: object of output
# --------- 
### recall: if p<epsilon, reject null that cov mats are homogeneous
###  See Rencher (1998) for info on F test
##########################################################
boxM = function(Y,group,f.indicator = FALSE){
  ## See Rencher (1998) for info on F test
  
  dname <- deparse(substitute(Y))
  if (!inherits(Y, c("data.frame", "matrix")))
    stop(paste(dname, "must be a numeric data.frame or matrix!"))
  if (length(group) != nrow(Y))
    stop("incompatible dimensions in Y and group!")
  
  Y <- as.matrix(Y)
  gname <- deparse(substitute(group))
  if (!is.factor(group)) group <- as.factor(as.character(group))
  
  valid <- complete.cases(Y, group)
  if (nrow(Y) > sum(valid)) 
    warning(paste(nrow(Y) - sum(valid)), " cases with missing data have been removed.")
  Y <- Y[valid,]
  group <- group[valid]
  
  p <- ncol(Y)
  nlev <- nlevels(group)
  lev <- levels(group)
  dfs <- tapply(group, group, length) - 1
  if (any(dfs < p)) 
    warning("there are one or more levels with less observations than variables!")
  mats <- aux <- list()
  for(i in 1:nlev) {
    mats[[i]] <- cov(Y[group == lev[i], ])
    aux[[i]] <- mats[[i]] * dfs[i]
  }
  names(mats) <- lev
  pooled <- Reduce("+", aux) / sum(dfs)
  logdet <- log(unlist(lapply(mats, det)))
  minus2logM <- sum(dfs) * log(det(pooled)) - sum(logdet * dfs)
  sum1 <- sum(1 / dfs) 
  sum2 <- sum((1 / dfs)^2)
  A1 <- (((2 * p^2) + (3 * p) - 1) / (6 * (p + 1) *
                                        (nlev - 1))) * (sum1 - (1 / sum(dfs)))
  dfchi <- (choose(p, 2) + p) * (nlev - 1)
  
  ## if group sample sizes are small, use F test
  if (any(dfs<15) & f.indicator == F)
    warning("due to small sample size, we recommend the use of an F test!")
  if (!any(dfs<15) & f.indicator == T)
    warning("due to large sample sizes, we do not recommend the use of an F test!")
  
  if ( f.indicator == TRUE ){
    A2 <- (p-1) * (p+2) / (6 * (nlev - 1) ) * (sum2 - (1 / sum(dfs)^2))
    if (A2-A1^2 > 0){
      dfF = (dfchi + 2)/(A2-A1^2)
      b = dfchi/(1-A1-dfchi/dfF)
      X2 <- minus2logM / b
      pval <- pf(X2, dfchi, dfF, lower.tail=FALSE)
    } else if (A2-A1^2<0){
      dfF = (dfchi + 2)/(A1^2-A2)
      b = dfF / (1-A1+2/dfF)
      X2 <- dfF * minus2logM / (dfchi * (b - minus2logM))
      pval <- pf(X2, dfchi, dfF, lower.tail=FALSE)
    } else {
      warning("error")
    }
  } else {
    X2 <- minus2logM * (1 - A1)
    pval <- pchisq(X2, dfchi, lower.tail = FALSE)
  }
  
  means <- aggregate(Y, list(group), mean)
  rn <- as.character(means[,1])
  means <- means[,-1]
  means <- rbind( means, colMeans(Y) )
  rownames(means) <- c(rn, "pooled")
  
  logdet <- c(logdet, pooled=log(det(pooled)))
  df <- c(dfs, pooled=sum(dfs))
  out <- structure(
    list(statistic = c("Chi-Sq (approx.)" = X2),
         parameter = c(df = dfchi),
         p.value = pval,
         cov = mats, pooled = pooled, logDet = logdet, means = means, df=df,
         data.name = dname, group = gname,
         method = "Box's M-test for Homogeneity of Covariance Matrices"
    ),
    class = c("htest", "boxM")
  )
  return(out)
}
##########################################################
## tic/ toc
##########################################################
tic = function(){
  Sys.time()
}
toc = function(tic){
  Sys.time()-tic
}
##########################################################
## sq_exp_distance_matrix
# vector of lat locs 
# vector of lon locs
# spat.names optional unique names corresponding to vectors of lat/lon
#
# output is a squared exponential distance matrix
##########################################################
sq_exp_distance_matrix = function(lat,lon,spat.names = 1:length(lat)){
  J = length(lat)
  W = matrix(data=0,nrow=J,ncol=J)
  for ( j in 1:J ){
    
    lon.j = lon[j]
    lat.j = lat[j]
    
    dist.from.j = sapply(1:J,function(k)sqrt((lon.j-lon[k])^2+
                                               (lat.j-lat[k])^2) )
    stand.dist = exp(-dist.from.j^2)
    
    W[j,] = stand.dist
    
  }
  diag(W) = 0
  colnames(W) = rownames(W) = spat.names
  
  return(W)
}
##########################################################
## row_standardize
# matrix W
#
# output is a row standardized matrix W
##########################################################
row_standardize = function(W){
  diag(1/rowSums(W)) %*% W
}
##########################################################
## lbind
# input L: list of length g of objects of same dimension
#
# output row bind the entries in L
##########################################################
lbind = function(L){
  # stack (rbind) matrices of dim nxp in list L of length g
  
  if(is.matrix(L[[1]])){ ## if matrix- put in vec and return group indicating membership of rows
    p = ncol(L[[1]])
    ns = unlist(lapply(L,nrow)); N = sum(ns)
    g = length(ns)
    mat.out = matrix(NA,ncol = p,nrow = N)
    group = rep(1:g,times =  ns)
    
    for ( j in 1:g){
      mat.out[group == j,] = as.matrix(L[[j]])
    }
    out = list("mat" = mat.out,"group" = group)
  } else if (is.vector(L[[1]])){
    p = length(L[[1]])

    mat.out = matrix(unlist(L),ncol = p, byrow = T)
    out = list("mat" = mat.out)
  }
  

  return(out)
}
draw_line = function(X,Y,color){
  if ( is.matrix(Y)){
    Y.plot = rbind(Y,NA)
  } else {
    Y.plot = c(Y,NA)
  }
  lines(rep(X,each = 3),Y.plot,col = color)
  # points(rep(X,each = 2),rbind(Y),col = color,pch="-")
}
##########################################################
## rwish
## copied directly from hoff amen
# input S0: centrality parameter
# input nu: degrees of freedom parameter
#
# output one random sample from wishart with mean S0*nu
##########################################################
rwish <-function(S0,nu=dim(S0)[1]+2){
  # sample from a Wishart distribution 
  # with expected value nu*S0 
  sS0<-(chol(S0))
  Z<-matrix(rnorm(nu*dim(S0)[1]),nu,dim(S0)[1])%*% sS0
  t(Z)%*%Z
}
##########################################################
## correlation matrix
#
# p : matrix dimension
# rho: off-diag value in (-1,1)
#
# create a p-dimensional correlation matrix with correlation rho
##########################################################
cor.mat <-function(p,RHO){
  eye(p)*(1-RHO) + RHO
}
##########################################################
## Auto-regressive AR(1) correlation matrix
#
# p : matrix dimension
# rho: off-diag value in (-1,1); first off-diagonal is rho^1
#
# create a p-dimensional correlation matrix with AR(1) structure
##########################################################
ar1.cor.mat <- function(p,RHO) {
  exponent <- abs(matrix(1:p - 1, nrow = p, ncol = p, byrow = TRUE) - 
                    (1:p - 1))
  RHO^exponent
}
cov.mle = function(X){
  ## MLE covariance matrix
  ## X: n x p matrix that we want a pxp column sample covariance matrix from
  
  X = as.matrix(X)
  N = nrow(X)
  
  X.temp = t(apply(X,1,function(KK) KK - colMeans(X)))
  
  t(X.temp) %*% X.temp / (N)
}
##########################################################
## ones
#
# p: length of vector
#
# create a length p vector filled with 1s
##########################################################
ones = function(p){
  matrix(rep(1,p),ncol=1)
}
##########################################################
## Reflect value to positive side of A
##########################################################
reflect = function(X,A){
  if (X<A) {
    out = A + (A-X)
  } else {
    out = X 
  }
  return(out)
}
##########################################################
## Reflecting random walk: sample a value from a delta interval around center and reflect to the right side of L
##########################################################
RRW = function(center,delta,L){
  reflect(sample((center - delta):(center + delta),1),L)
}
##########################################################
## Extract matrix dimension indices from indieces of the vectorized matrix
## where the vectorized matrix stacks the columns of the matrix 
## on top of each other into a vector
##########################################################
vec_ind_to_mat_ind <- function(ind, p) {
  # ind : vectorized indices 
  # p : dimension of the p x p matrix
  
  row = ((ind - 1) %% p) + 1
  col = ((ind - 1) %/% p) + 1
  
  data.frame(row = row, col = col)
}
##########################################################
## Extract vectorized dimension indices from indices of a matrix
## where the vectorized version stacks the columns of the matrix 
## on top of each other into a vector
##########################################################
mat_ind_to_vec_ind <- function(row, col, p) {
  # row, col: indices of the matrix 
  # p : dimension of the p x p matrix
  
  ind = row + (col - 1) * p
  return(ind)
}
