da.function = function(x,mu,cov.hat){
  ## x is test data- vector of length p
  ## mu is list of length g of vectors of length p- mean estimate
  ## cov is list of length g of matrices of dim n_jxp- cov estimate
  
  g = length(mu)
  # he lper = lapply(mu,function(k)x-k)
  d.hat = sapply(1:g,function(j)
    multAprimeBA((x-mu[[j]]),csolve(cov.hat[[j]])) + determinant(cov.hat[[j]])$mod[1]
    )
  
  return(which(d.hat==min(d.hat)))

}
da.function.1word.manymu = function(x,mu,cov.hat){
  ## x is test data- vector of length p
  ## mu is list mean estimate
  ## cov is cov estimate
  
  g = length(mu)
  d.hat = sapply(1:g,function(j)
    t(x-mu[[j]]) %*% csolve(cov.hat) %*% (x-mu[[j]]) + determinant(cov.hat)$mod[1]
    #multAprimeBA((x-mu[[j]]),csolve(cov.hat)) + determinant(cov.hat)$mod[1]
  )
  
  return(which(d.hat==min(d.hat)))
  
}
da.mean.est = function(X){
  ## X is list of length g of n_jxp matrices
  
  MU = lapply(X,colMeans)
  
  return(MU)
}
pooled.cov = function(X){
  ## X is list of length g of n_jxp matrices
  
  E = lapply(X,function(j)t(apply(j,1,function(k)k-colMeans(j))))

  Si = lapply(E,function(j)t(j)%*%(j))
  
  ns = unlist(lapply(X, nrow))
  
  Sp = Reduce('+',Si)/(sum(ns)-length(X))
  
  return(Sp)
}
rep.list = function(L,p){
  out = list()
  for(k in 1:p){
    out[[k]] = L
  }
  return(out)
}
