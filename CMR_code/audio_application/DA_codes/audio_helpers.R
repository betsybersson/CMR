#############################################
## general
#############################################

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

dmatT.faster = function(eig.helper,
                        nu,nu0,
                        p,n){
  # output is log output
  # helper = (X-M) %*% Omega.inv %*% t(X-M) * (nu0-p-1)

  # if X is a vector not a matrix
  if (is.null(p)){
    p = length(X)
    n=1
  }
  
  # help with huge numbers
  main.det = determinant(eye(n) + eig.helper/(nu0-p-1))$modulus[1]
  # equal to main.det, probably faster but introduces complex numbers in case of numerical error: 
  # main.det = sum(log(eig.helper/(nu0-p-1) + 1))
  gamma.div = gamma.mv(p,(nu+n+p-1)/2,T) - gamma.mv(p,(nu+p-1)/2,T)
  
  out = gamma.div  + (-p*n/2) * log(nu0-p-1) + 
    (-(nu+n+p-1)/2) * (main.det)
  
  return(out)
  
}
vec.3d = function(XX){
  ## input: XX of dim (n,p1,p2)
  ## output: matrix of dim (n,p1*p2); 
  ## (columns p2 stacked on top of each other)
  t(apply(XX,1,c))
}
rwish.wrapper <-function(S0,nu=dim(S0)[1]+2){
  # sample from a Wishart distribution 
  # with expected value nu*S0 
  Z = matrix(rnorm(nu*dim(S0)[1]),nu,dim(S0)[1])
  crwish(S0,Z)
}
rmvnorm.wrapper = function(M,V = eye(ncol(M))) {
  # function for density of Y_{ n\times p} \sim N(M, V kron U) where U = eye(nrow(M))
  N = nrow(M)
  P = ncol(M)
  
  E = matrix(rnorm(N*P),ncol=P)
  
  out = crmvnorm(M,V,E)
  
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
dmvnorm.faster = function(S,N,V,V.inv = csolve(V)){
  P = ncol(S)

  out = - N*determinant(V)$mod[1]/2 - mat_trace( multAB(V.inv,S) )/2
  
  return(out)
}
dwish.const.prop = function(nu,P) {

  out = (-1) * ((nu*P/2)*log(2) + sum(lgamma( (nu+1-c(1:P))/2 )) +
                  (nu/2)*(-P)*log(nu))
  
  return(out)
}

dwish_log = function(X,M,nu){
  ## log = true
  p = ncol(M)
  (-1) * ((nu*p/2)*log(2) + 
            (p*(p-1)/4)*log(pi) +
            sum(lgamma( (nu+1-c(1:p))/2 )) +
            (nu/2)*determinant(M)$mod[1]) +
    ((nu-p-1)/2) * determinant(X)$mod[1] -
    mat_trace(-csolve(M)*X/2)
    
}

sampleNu0 = function(j){
  # helper = eigen(multABAprime(U[[j]],  V0.inv),
  #                only.values = T,symmetric=T)$val
  helper = multABAprime(U[[j]],  V0.inv)
  
  dmatT.faster(helper,nu0.star-p+1,nu0.star,p,ns[j]) - 
    dmatT.faster(helper,nu0-p+1,nu0,p,ns[j])
}
samplePsi = function(j){
  M = csolve(V0*(nu0 - p - 1) + multAprimeA(U[[j]]) )
  temp.Psi.inv = rwish.wrapper(M,nu0 + ns[j] - 1)
  temp.Psi = csolve(temp.Psi.inv)
  return(list("Psi"=temp.Psi,"Psi.inv"=temp.Psi.inv))
}
sampleNu = function(j){
  Vj.inv = kronecker(V2.inv[[j]],V1.inv[[j]])
  helper = multABAprime(Y.list[[j]] - pis^(1/2)*U[[j]],  Vj.inv) / (1-pis)
  
  dmatT.faster(helper,nu.star-p+1,nu.star,p,ns[j]) - 
    dmatT.faster(helper,nu-p+1,nu,p,ns[j])
}
sampleSig = function(j){
  V = kronecker(V2[[j]],V1[[j]])
  Y.tilde = Y.list[[j]] - pis^(1/2)*U[[j]]
  M = csolve(V * (nu - p - 1) + multAprimeA(Y.tilde) / (1-pis))
  temp.Sig.inv = rwish.wrapper(M,nu+ns[j]-1)
  temp.Sig = csolve(temp.Sig.inv)
  return(list("Sig"=temp.Sig,"Sig.inv"=temp.Sig.inv))
}
wDensity = function(j){
  dmvnorm.faster(S.list[[j]],ns[j],pis.star * Psi[[j]] + (1-pis.star) * Sig[[j]]) - 
    dmvnorm.faster(S.list[[j]],ns[j],pi.hat * Psi[[j]] + (1-pi.hat) * Sig[[j]])
}
sampleU = function(j){
  W = csolve(pis/(1-pis)*Sig.inv[[j]] + Psi.inv[[j]])
  M = pis^(1/2)/(1-pis) * multABC(Y.list[[j]], Sig.inv[[j]], W)
  rmvnorm.wrapper(M,W)
}
sampleV1 = function(j){
  L = array(Sig.inv.chol[[j]],dim=c(p1,p2,p)) # checked- should be correct
  helper = lapply(1:p,function(k) multABAprime(L[,,k], (V2[[j]])) )
  M = csolve(Reduce('+',helper)*(nu-p-1) + S1.inv * eta1)
  temp.V1 = rwish.wrapper(M,eta1 + nu * p2)
  temp.V1.inv = csolve(temp.V1)
  return(list("V1"=temp.V1,"V1.inv"=temp.V1.inv))
}
sampleV2 = function(j){
  L = array(Sig.inv.chol[[j]],dim=c(p1,p2,p))
  helper = lapply(1:p,function(k) multAprimeBA(L[,,k], (V1[[j]])) )
  M = csolve(Reduce('+',helper)*(nu-p-1) + S2.inv * eta2)
  temp.V2 = rwish.wrapper(M,eta2 + nu * p1)
  temp.V2.inv = csolve(temp.V2)
  return(list("V2"=temp.V2,"V2.inv"=temp.V2.inv))
}
SigInvChol = function(j){
  ccholL(Sig.inv[[j]])
}
lbind = function(L){
  # stack (rbind) matrices of dim njxp in list L of length g
  
  p = ncol(L[[1]])
  ns = unlist(lapply(L,nrow)); N = sum(ns)
  g = length(ns)
  
  mat.out = matrix(NA,ncol = p,nrow = N)
  group = rep(1:g,times =  ns)
  
  for ( j in 1:g){
    mat.out[group == j,] = as.matrix(L[[j]])
  }
  
  return(list("mat" = mat.out,"group" = group))
}

vec.inv.array = function(x,p1,p2){
  ## input: x: n x p1p2 matrix
  ## output: X: p1 x p2 x n array 
  
  N = nrow(x)
  XX = array(NA,dim=c(p1,p2,N))
  
  for ( nn in 1:N ){
    XX[,,nn] = matrix(x[nn,],ncol=p2,nrow=p1)
  }
  
  return(XX)
  
}
list.to.3d.array = function(L){
  # input: list of length J containing matrices of dimension p1 x p2
  # output: array of dimension p1 x p2 x J 
  
  J = length(L)
  p1 = nrow(L[[1]])
  p2 = ncol(L[[1]])
  
  array.out = array(NA,dim=c(p1,p2,J))
  for ( jj in 1:J ){
    array.out[,,jj] = L[[jj]]
  }
  
  return(array.out)
  
}
rep.array = function(A,p){
  ## repeat matrix A p times
  ## output: X : array of dimension (dim(A),p)
  p1 = nrow(A)
  p2 = ncol(A)
  out.array = array(NA,dim=c(p1,p2,p))
  for(j in seq_len(p)){
    out.array[,,j] = A
  }
  return(out.array)
}