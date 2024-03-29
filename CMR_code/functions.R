library(LaplacesDemon) # for dmvt (multivariate Student's t) for cusp

##########################################################
### Covariance Meta Regression Model Gibbs Sampler- CMR_GS
## data
# Y: n x p response matrix
# X: p x q meta covariate matrix (optional)
# k: dimension of latent factor
## classic GS parameters
# S: number of iterations for GS
# burnin: number of iterations to burnin
# thin: thinning mechanism
# save.all: logical to determine what parameters to return. If 0, save subset 
#       of params post thin/burnin. If 1, save and return all iterations 
#       of all parameters.
## optional parameters
# 
##########################################################
CMR_GS = function(Y,X = NA,
                  k = round(ncol(Y)/2),
                  S = 5100,
                  burnin = 100,
                  thin = 1,
                  save.all = 0,
                  my.seed = Sys.time()){
  
  set.seed(my.seed)
  
  ###########################
  ## set hyper params
  p = ncol(Y)
  n = nrow(Y)

  ## handle some things for meta covariates
  ## incl some warnings to make sure dimensions are right
  if (is.vector(X)){
    X = matrix(X,ncol=1)
  }
   
  if (all(is.na(X))){
    # if X is NA, make it the intercept
    X = matrix(1,ncol=1,nrow = p)
  }
  
  if (nrow(X)!=p){
    stop("Y,X dimension mismatch!")
  }
  
  q = ncol(X)
  
  ## IG shape/rate parameters
  a.d = b.d = 1
  a.xi = b.xi = 1
  a.tau = b.tau = 1
  ###########################
  
  ###########################
  ## initialize values
  # sampling model
  Lambda = matrix(0,ncol = k, nrow = p)
  eta = matrix(0,ncol = k, nrow = n)
  ds = rep(1,p); D = D.inv = diag(ds)
  # factor model
  Gamma = matrix(0,nrow = q, ncol = k)
  tau2 = 1; Tau = Tau.inv = eye(k)*tau2
  # variable selection scale
  xi2 = 1
  ###########################
  
  ###########################
  ## create storage 
  if (save.all == 1){
    n.save.out = S
  } else if (save.all == 0) {
    len.out = length(seq(from = burnin+1, to = S, by = thin))
    n.save.out = len.out
    index = 1  
  }
  cov.out = cov.inv.out = array(NA,dim = c(n.save.out,p*p))
  Lambda.out = array(NA,dim = c(n.save.out,p*k))
  eta.out = array(NA,dim=c(n.save.out, k*n))
  ds.out = array(NA,dim=c(n.save.out,p))
  Gamma.out = array(NA,dim=c(n.save.out,(q)*k))
  vars.out = array(NA,dim=c(n.save.out,2)) #for tau2 and xi2
  ###########################
  
  ###########################
  ## global helpers
  
  ###########################
  
  ###########################
  ## GS
  tic = Sys.time()
  for ( s in 1:S ){
    
    ## sample etais
    R = t(Lambda) %*% D.inv
    S.eta.inv = qr.solve(R %*% Lambda + eye(k))
    M.eta.help = S.eta.inv %*% R
    
    eta.list = mclapply(1:n,function(i) rmvnorm(M.eta.help %*% Y[i,],S.eta.inv))
    eta = rbind.list(eta.list)
    
    ## sample djs
    for ( j in 1:p ){
      lambda.j = Lambda[j,]
      helper = (lambda.j - t(Gamma) %*% X[j,])
      S.dj = sum((Y[,j] -  eta %*% lambda.j)^2) + t(helper) %*% helper / tau2 + b.d
      ds[j] = 1/rgamma(1,(n+k+a.d)/2,c(S.dj)/2)
    }
    # update diagonal matrices
    D = diag(ds); D.inv = diag(1/ds)
    
    ## sample Lambda
    S.gamma.inv = qr.solve(Tau.inv + t(eta) %*% eta)
    M.gamma = X %*% Gamma / tau2 + t(Y) %*% eta
    Lambda = rmatnorm(M.gamma %*% S.gamma.inv,
                      V = S.gamma.inv, U = D)
    
    # ## sample gamma0
    # s.gamma.inv = tau2/(sum(1/ds) + tau2)
    # S.gamma.inv = eye(k) * s.gamma.inv
    # M.gamma.helper = sapply(1:p,function(j)(Lambda[j,] - t(Gamma) %*% X[j,])/ds[j]) %>% 
    #   rowSums()
    # gamma0 = rmvnorm(M.gamma.helper/tau2*s.gamma.inv,S.gamma.inv)
    # Gamma0 = matrix(rep(gamma0,p), nrow = p, byrow = T)

    ## sample Gamma
    R = t(X) %*% D.inv
    S.gamma.inv = qr.solve(R %*% X + eye(q)/xi2)
    M.gamma = R %*% Lambda
    Gamma = rmatnorm(S.gamma.inv %*% M.gamma, V = Tau, U = S.gamma.inv)
  
    
    ## sample xi2
    GtG_trace = mat_trace(Gamma %*% t(Gamma))
    xi2 = 1/rgamma(1, (q*k + a.xi)/2, 
                 (GtG_trace/tau2 + b.xi)/2)
 
    ## sample tau2
    helper = Lambda - X %*% Gamma
    S.tau = mat_trace(t(helper) %*% D.inv %*% helper) + GtG_trace/xi2 + b.tau
    tau2 = 1/rgamma(1, (p*k + q*k + a.tau)/2, S.tau/2)
    Tau = eye(k)*tau2; Tau.inv = eye(k)/tau2
    
    ## save output
    if (save.all == 0){
      if ((s>burnin)&((s %% thin)==0)){
        
        cov.temp = Lambda %*% t(Lambda) + D
        cov.inv.temp = qr.solve(cov.temp)
        
        cov.out[index,] = c(cov.temp)
        cov.inv.out[index,] = c(cov.inv.temp)
        Lambda.out[index,] = c(Lambda)
        eta.out[index,] = c(eta)
        ds.out[index,] = c(ds)
        Gamma.out[index,] = c(Gamma)
        vars.out[index,] = c(tau2,xi2) #for tau2 and xi2
        
        index = index + 1
      }
    } else {
      cov.temp = Lambda %*% t(Lambda) + D
      cov.inv.temp = qr.solve(cov.temp)
      
      cov.out[s,] = c(cov.temp)
      cov.inv.out[s,] = c(cov.inv.temp)
      Lambda.out[s,] = c(Lambda)
      eta.out[s,] = c(eta)
      ds.out[s,] = c(ds)
      Gamma.out[s,] = c(Gamma)
      vars.out[s,] = c(tau2,xi2) #for tau2 and xi2
    }
    
  } # end GS iteration
  ###########################
  runtime = difftime(Sys.time(),tic,units = "secs")
  
  ###########################
  ## put output in list and save to function
  colnames(vars.out) = c("tau2","xi2")
  
  func.out = list("cov" = cov.out,"cov.inv" = cov.inv.out,
                  "Lambda" = Lambda.out,
                  "eta" = eta.out,
                  "D" = ds.out,
                  "Gamma" = Gamma.out,
                  "vars" = vars.out,
                  "runtime" = runtime)
  if (save.all == 1){
    func.out = func.out #c(func.out,list("out2" = out2))
  }
  
  set.seed(Sys.time())
  
  return(func.out)
  ###########################
  
} # end function

###########################
## rbind.list - rbind a list of vectors into a single matrix 
###########################
rbind.list = function(L){
  n.col = length(L[[1]])
  n.row = length(L)
  matrix(unlist(L),ncol = n.col, nrow = n.row,
         byrow = T)
}

###########################
## Latent factor model (basic)
# latent factor drawn from mean 0 normal distribution
###########################
LFM_GS = function(Y,
                  k = round(ncol(Y)/2),
                  S = 5100,
                  burnin = 100,
                  thin = 1,
                  save.all = 0,
                  my.seed = Sys.time()){
  
  set.seed(my.seed)
  
  ###########################
  ## set hyper params
  p = ncol(Y)
  n = nrow(Y)
  
  ## IG shape/rate parameters
  a.d = b.d = 1
  a.tau = b.tau = 1
  ###########################
  
  ###########################
  ## initialize values
  # sampling model
  Lambda = matrix(0,ncol = k, nrow = p)
  eta = matrix(0,ncol = k, nrow = n)
  ds = rep(1,p); D = D.inv = diag(ds)
  # factor model
  tau2 = 1; Tau = Tau.inv = eye(k)*tau2
  ###########################
  
  ###########################
  ## create storage 
  if (save.all == 1){
    n.save.out = S
  } else if (save.all == 0) {
    len.out = length(seq(from = burnin+1, to = S, by = thin))
    n.save.out = len.out
    index = 1  
  }
  cov.out = cov.inv.out = array(NA,dim = c(n.save.out,p*p))
  Lambda.out = array(NA,dim = c(n.save.out,p*k))
  eta.out = array(NA,dim=c(n.save.out, k*n))
  ds.out = array(NA,dim=c(n.save.out,p))
  tau2.out = array(NA,dim=c(n.save.out,1)) 
  ###########################
  
  ###########################
  ## global helpers
  
  ###########################
  
  ###########################
  ## GS
  tic = Sys.time()
  for ( s in 1:S ){
    
    ## sample etais
    R = t(Lambda) %*% D.inv
    S.eta.inv = qr.solve(R %*% Lambda + eye(k))
    M.eta.help = S.eta.inv %*% R
    
    eta.list = mclapply(1:n,function(i) rmvnorm(M.eta.help %*% Y[i,],S.eta.inv))
    eta = rbind.list(eta.list)
    
    ## sample djs
    for ( j in 1:p ){
      lambda.j = Lambda[j,]
      S.dj = sum((Y[,j] -  eta %*% lambda.j)^2) + t(lambda.j) %*% lambda.j / tau2 + b.d
      ds[j] = 1/rgamma(1,(n+k+a.d)/2,c(S.dj)/2)
    }
    # update diagonal matrices
    D = diag(ds); D.inv = diag(1/ds)
    
    ## sample Lambda
    S.gamma.inv = qr.solve(Tau.inv + t(eta) %*% eta)
    M.gamma = t(Y) %*% eta
    Lambda = rmatnorm(M.gamma %*% S.gamma.inv,
                      V = S.gamma.inv, U = D)
    
    ## sample tau2
    S.tau = mat_trace(t(Lambda) %*% D.inv %*% Lambda) + b.tau
    tau2 = 1/rgamma(1, (p*k + a.tau)/2, S.tau/2)
    Tau = eye(k)*tau2; Tau.inv = eye(k)/tau2
    
    ## save output
    if (save.all == 0){
      if ((s>burnin)&((s %% thin)==0)){
        
        cov.temp = Lambda %*% t(Lambda) + D
        cov.inv.temp = qr.solve(cov.temp)
        
        cov.out[index,] = c(cov.temp)
        cov.inv.out[index,] = c(cov.inv.temp)
        Lambda.out[index,] = c(Lambda)
        eta.out[index,] = c(eta)
        ds.out[index,] = c(ds)
        tau2.out[index,] = c(tau2) 
        
        index = index + 1
      }
    } else {
      cov.temp = Lambda %*% t(Lambda) + D
      cov.inv.temp = qr.solve(cov.temp)
      
      cov.out[s,] = c(cov.temp)
      cov.inv.out[s,] = c(cov.inv.temp)
      Lambda.out[s,] = c(Lambda)
      eta.out[s,] = c(eta)
      ds.out[s,] = c(ds)
      tau2.out[s,] = c(tau2) 
    }
    
  } # end GS iteration
  ###########################
  runtime = difftime(Sys.time(),tic,units = "secs")
  
  ###########################
  ## put output in list and save to function

  func.out = list("cov" = cov.out,"cov.inv" = cov.inv.out,
                  "Lambda" = Lambda.out,
                  "eta" = eta.out,
                  "D" = ds.out,
                  "tau2" = tau2.out,
                  "runtime" = runtime)

  
  set.seed(Sys.time())
  
  return(func.out)
  ###########################
  
} # end function

##########################################################################
# Adaptive Gibbs sampler for Gaussian factor models under the CUSP prior
# source:https://github.com/siriolegramanti/CUSP/blob/master/cusp.R
##########################################################################
cusp_factor_adapt <- function(y,my_seed,N_sampl,alpha,a_sig,b_sig,a_theta,b_theta,
                              theta_inf,start_adapt,Hmax,alpha0,alpha1,
                              save.all = 0,
                              burnin = round(N_sampl/10),
                              thin = 1){
  set.seed(my_seed)
  # statistical units
  n<-dim(y)[1]
  # number of variables
  p<-dim(y)[2]
  # uniform rvs for adaptation
  u<-runif(N_sampl)
  ########################################################################
  # VARIABLES TO UPDATE
  ########################################################################
  # number of factors
  H<-rep(NA,N_sampl)
  # number of active factors
  Hstar<-rep(NA,N_sampl)
  # factor loadings matrix
  Lambda_post<-vector("list",N_sampl)
  # factor scores matrix
  eta_post<-vector("list",N_sampl)
  # column-specific loadings precisions
  theta_inv_post<-vector("list",N_sampl)
  # stick-breaking weights
  w_post<-vector("list",N_sampl)
  # augmented data
  z_post<-vector("list",N_sampl)
  # error variances
  inv_sigma_sq_post<-matrix(NA,p,N_sampl)
  ########################################################################
  #INITIALIZATION
  ########################################################################
  H[1]<-p+1
  Hstar[1]<-p
  Lambda_post[[1]]<-matrix(rnorm(p*H[1]),p,H[1])
  eta_post[[1]]<-matrix(rnorm(n*H[1]),n,H[1])
  theta_inv_post[[1]]<-rep(1,H[1])
  w_post[[1]]<-rep(1/H[1],H[1])
  z_post[[1]]<-rep(H[1],H[1])
  inv_sigma_sq_post[,1]<-1
  ########################################################################
  # GIBBS SAMPLING #######################################################
  ########################################################################
  t0 <- Sys.time()
  for (t in 2:N_sampl){
    Lambda_post[[t]]<-matrix(NA,p,H[t-1])
    eta_post[[t]]<-matrix(NA,n,H[t-1])
    theta_inv_post[[t]]<-rep(NA,H[t-1])
    w_post[[t]]<-rep(NA,H[t-1])
    z_post[[t]]<-rep(NA,H[t-1])
    # 1) sample the j-th row of Lambda
    for (j in 1:p){
      V_j<-chol2inv(chol(diag(theta_inv_post[[t-1]],nrow=H[t-1])+inv_sigma_sq_post[j,t-1]*(t(eta_post[[t-1]])%*%eta_post[[t-1]])))
      mu_j<-V_j%*%t(eta_post[[t-1]])%*%(y[,j])*inv_sigma_sq_post[j,t-1]
      Lambda_post[[t]][j,]<-rmvnorm(mu_j,V_j) #mvrnorm(1,mu_j,V_j)
    }
    # 2) sample sigma^{-2}
    diff_y<-apply((y-(eta_post[[t-1]]%*%t(Lambda_post[[t]])))^2,2,sum)
    for (j in 1:p){
      inv_sigma_sq_post[j,t]<-rgamma(n=1,shape=a_sig+0.5*n,rate=b_sig+0.5*diff_y[j])
    }		
    # 3) sample eta
    V_eta<-chol2inv(chol(diag(H[t-1])+t(Lambda_post[[t]])%*%diag(inv_sigma_sq_post[,t])%*%Lambda_post[[t]]))
    for (i in 1:n){
      mu_eta_i<-V_eta%*%t(Lambda_post[[t]])%*%diag(inv_sigma_sq_post[,t])%*%y[i,]
      eta_post[[t]][i,]<-rmvnorm(mu_eta_i,V_eta) #mvrnorm(1,mu_eta_i,V_eta)
    }	
    # 4) sample z
    lhd_spike<-rep(0,H[t-1])
    lhd_slab<-rep(0,H[t-1])
    for (h in 1:H[t-1]){
      lhd_spike[h]<-exp(sum(log(dnorm(Lambda_post[[t]][,h], mean = 0, sd = theta_inf^(1/2), log = FALSE))))
      lhd_slab[h]<-dmvt(x=Lambda_post[[t]][,h], mu=rep(0,p), S=(b_theta/a_theta)*diag(p), df=2*a_theta)
      prob_h<-w_post[[t-1]]*c(rep(lhd_spike[h],h),rep(lhd_slab[h],H[t-1]-h))
      if (sum(prob_h)==0){
        prob_h<-c(rep(0,H[t-1]-1),1)
      }
      else{
        prob_h<-prob_h/sum(prob_h)
      }
      z_post[[t]][h]<-c(1:H[t-1])%*%rmultinom(n=1, size=1, prob=prob_h)
    }
    # 5) sample v and update w
    v<-rep(NA,H[t-1])
    for (h in 1:(H[t-1]-1)){
      v[h]<-rbeta(1, shape1 = 1+sum(z_post[[t]]==h), shape2 = alpha+sum(z_post[[t]]>h))
    }
    v[H[t-1]]<-1
    w_post[[t]][1]<-v[1]
    for (h in 2:H[t-1]){
      w_post[[t]][h]<-v[h]*prod(1-v[1:(h-1)])  
    }
    # 6) sample theta^{-1}
    for (h in 1:H[t-1]){
      if (z_post[[t]][h]<=h){
        theta_inv_post[[t]][h]<-theta_inf^(-1)
      }
      else{
        theta_inv_post[[t]][h]<-rgamma(n=1,shape=a_theta+0.5*p,rate=b_theta+0.5*t(Lambda_post[[t]][,h])%*%Lambda_post[[t]][,h])
      }
    }
    # 7) update H[t]
    active <- which(z_post[[t]]>c(1:H[t-1]))
    Hstar[t] <- length(active)
    H[t]<-H[t-1]
    # by default, keep the same truncation, unless...
    if (t>=start_adapt & u[t]<=exp(alpha0+alpha1*t)){
      if (Hstar[t]<H[t-1]-1){
        # set truncation to Hstar[t] and subset all variables, keeping only active columns
        H[t]<-Hstar[t]+1
        eta_post[[t]]<-cbind(eta_post[[t]][,active],rnorm(n))
        theta_inv_post[[t]]<-c(theta_inv_post[[t]][active],theta_inf^(-1))
        w_post[[t]]<-c(w_post[[t]][active],1-sum(w_post[[t]][active]))
        Lambda_post[[t]]<-cbind(Lambda_post[[t]][,active],rnorm(p,mean=0,sd=sqrt(theta_inf)))
      } else if (H[t-1]<Hmax) {
        # increase truncation by 1 and extend all variables, sampling from the prior/model
        H[t]<-H[t-1]+1
        eta_post[[t]]<-cbind(eta_post[[t]],rnorm(n))
        v[H[t-1]]<-rbeta(1,shape1=1,shape2=alpha)
        v<-c(v,1)
        w_post[[t]]<-rep(NA,H[t])
        w_post[[t]][1]<-v[1]
        for (h in 2:H[t]){
          w_post[[t]][h]<-v[h]*prod(1-v[1:(h-1)])  
        }
        theta_inv_post[[t]]<-c(theta_inv_post[[t]],theta_inf^(-1))
        Lambda_post[[t]]<-cbind(Lambda_post[[t]],rnorm(p,mean=0,sd=sqrt(theta_inf)))
      }
    }
    ##print(t)
  }
  runtime <- difftime(Sys.time(),t0,units = "secs")
  
  ## get thinned cov.inv - Betsy added
  thinning <- seq(burnin+1,N_sampl,by=thin)
  
  cov.inv <- array(NA,c(length(thinning),p*p))
  
  for (i in 1:length(thinning)){
    cov.inv[i,] <- c(qr.solve(Lambda_post[[thinning[i]]]%*%t(Lambda_post[[thinning[i]]]) +
                                diag(1/inv_sigma_sq_post[,thinning[i]]))
    )
  }
  
  ## save output
  output = list("runtime"=runtime,"cov.inv" = cov.inv)
  if (save.all == 1){
    output = c(output,c("y"=y,"my_seed"=my_seed,"N_sampl"=N_sampl,
                        "alpha"=alpha,"a_sig"=a_sig,"b_sig"=b_sig,"a_theta"=a_theta,"b_theta"=b_theta,
                        "theta_inf"=theta_inf,"start_adapt"=start_adapt,"Hmax"=Hmax,
                        "alpha0"=alpha0,"alpha1"=alpha1,
                        "H"=H,"Hstar"=Hstar,"Lambda_post"=Lambda_post,"eta_post"=eta_post,
                        "theta_inv_post"=theta_inv_post,"w_post"=w_post,"z_post"=z_post,
                        "inv_sigma_sq_post"=inv_sigma_sq_post)
    )
  }
  
  return(output)
}