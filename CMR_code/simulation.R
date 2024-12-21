source("functions.R")
source("Schiavon_SISGaussian.R")
library(parallel)
library(doParallel)
####################################

####################################
## helpers
on.server = TRUE
cov.method = "kron" ## options: eye, cor9, comSym3groups, kron
identifier = "p50"
####################################

####################################
## problem dimension parameters

# number of variables
p = 50
# sample sizes to loop through
Ns =  c(p+1,round(p*1.5),round(p*3))#,p*3 
Ns.names = c("1","1.5","3") #,"3"
# low dimension
Ks = "" # c(1,3,5)
####################################

print(paste0("Running the following scenario: ",
             "cov: ", cov.method,
             "; p: ", p,
             "; Ns: (", paste0(Ns, collapse = ","),
             "); Ks: (", paste0(Ks, collapse = ","),
             ") !!!!!!!!!!!!"))

X = NA

####################################
## gibbs sampler variables
S.fancy = 20000
burnin.fancy = 10000

S.simple = 11000
burnin.simple = 1000
# number of simulation replicates
sim = 25
####################################

print(paste0("Using the following GS values:",
		    " S: ",S.fancy,
		    "; burnin: ", burnin.fancy,
		    "; sims: ",sim,
		    "!!!!!!!!!!!!!"))

####################################
## have correct path to read in R functions
if (on.server == T){
  source("betsy_R_functions.R")
} else {
  source("/Users/betsybersson/Documents/Betsy_R_functions.R")
}
####################################


####################################
## true covariance matrix
if (cov.method == "eye"){
  true.cov = eye(p)
} else if (cov.method == "cor9"){
  true.cov = cor.mat(p,0.9)
} else if (cov.method == "comSym3groups"){
  
  p_l = round(c(4,6,8)/sum(4,6,8)*p)
  L = length(p_l)
  p = sum(p_l)
  which_group = rep(1:L,times = p_l)
  
  lam = c(1,0,1,1)
  Lam = matrix(rep(lam,p_l[1]),nrow=p_l[1],byrow=T)
  lam2 = c(0,0,0,1)
  Lam2 = matrix(rep(lam2,p_l[2]),nrow=p_l[2],byrow=T)
  lam3 = c(1,4,0,0)
  Lam3 = matrix(rep(lam3,p_l[3]),nrow=p_l[3],byrow=T)
  Lam = rbind(Lam,Lam2,Lam3)
  cov = Lam %*% t(Lam)
  cov = cov + eye(p)
  true.cov = cov2cor(cov)
  eigen(true.cov)$val ## check invertible

  # get design matrix
  X = matrix(0, nrow = p, ncol = length(p_l))
  for ( pl.ind in 1:length(p_l)){
    X[which_group==pl.ind,pl.ind] = 1
  }
    
} else if (cov.method == "kron"){
  p1 = round(p/3)
  p2 = round(p/p1)
  p = p1*p2
  
  R = cor.mat(p1,.7)
  C = cor.mat(p2,.2)
  true.cov = kronecker(C,R)
  eigen(true.cov)$val ## check invertible
  
  # get design matrix
  X = getMatDesignMat(p1,p2)
  
}
# propagate filename suffix
suffix = cov.method
if (identifier != ""){
  suffix = paste0(suffix,"_",identifier)
}
####################################

####################################
## run simulation
cores = detectCores()
loss.avg = toc.avg = c() #matrix(NA,ncol = length(Ks)*2 + 2, nrow = length(Ns))
for ( n.ind in 1:length(Ns) ){
    
  n = Ns[n.ind]

  if (on.server == T){
    cl <- makeCluster(cores[1]-5)  
  } else {
    cl <- makeCluster(cores[1] - 1)  # dont overload your computer
  }
  registerDoParallel(cl)
  
  ###
  parallel.out <- foreach(sim.ind=1:sim, .combine=cbind, .packages = c("parallel","LaplacesDemon","pgdraw","calculus","matrixStats")) %dopar% {
    
    ## output storage
    output = toc = list()

    ####################################
    ## sample dataset
    set.seed(sim.ind)
    Y = rmatnorm(matrix(0,ncol = p, nrow = n),V=true.cov)
    set.seed(Sys.time())
    ####################################
    
    # ####################################   
    # ## be stupid and use pca to determin num of factors
    # baing_cr <- baingcriterion(tibble(Y), rmax = p)$IC[1:2]
    # baing_cr
    # ####################################
    
    # ####################################
    # ## run CMR GS
    # for ( k.ind in 1:length(Ks) ){
    #   out.cmr  = CMR_GS(Y,X,
    #                     k = Ks[k.ind],
    #                     S = S.simple,
    #                     burnin = burnin.simple,
    #                     my.seed = sim.ind + 100 + Ks[k.ind])
    #   output[[k.ind]] = qr.solve(matrix(colMeans(out.cmr$cov.inv),ncol = p)) ## stein estimator
    #   toc[[k.ind]]  = out.cmr$runtime
    # }
    # # name based on value of K
    # names(output) = names(toc) = paste0("cmr.K",Ks)
    # ####################################
    
    ####################################
    ## run CMR GS with CUSP
    DUNSON_ALPHA = floor(p/3) ## prior expected number of active factors- used for all cusp and sis
    
    out.cmr.cusp  = CMR_cusp_GS(Y,X,
                                k = floor((p-1)/2),
                                S = S.fancy,
                                burnin = burnin.fancy,
                                my.seed = sim.ind + 200,
                                alpha = DUNSON_ALPHA,
                                a.theta = 1/2, b.theta = 1/2)
    output$cmr.cusp = qr.solve(matrix(colMeans(out.cmr.cusp$cov.inv),ncol = p)) ## stein estimator
    toc$cmr.cusp  = out.cmr.cusp$runtime
    ####################################
    
    ####################################
    ## run intercept CMR GS if you haven't
    if (!all(is.na(X))){
      # for ( k.ind in 1:length(Ks) ){
      #   out.cmr  = CMR_GS(Y,X = NA,
      #                     k = Ks[k.ind],
      #                     S = S.simple,
      #                     burnin = burnin.simple,
      #                     my.seed = sim.ind + 300 + k.ind)
      #   output[[length(output)+1]] = qr.solve(matrix(colMeans(out.cmr$cov.inv),ncol = p)) ## stein estimator
      #   toc[[length(toc)+1]]  = out.cmr$runtime
      #   
      #   # name based on value of K
      #   names(output)[length(output)+1] = names(toc)[length(output)+1] = paste0("cmr.K",Ks[k.ind],".intercept")
      # }
      
      ## run CMR GS with CUSP
      out.cmr.cusp  = CMR_cusp_GS(Y,X = NA,
                                  k = ceiling((p-1)/2),
                                  S = S.fancy,
                                  burnin = burnin.fancy,
                                  my.seed = sim.ind + 400,
                                  alpha = DUNSON_ALPHA,
                                  a.theta = 1/2, b.theta = 1/2)
      
      output$cmr.cusp.intercept = qr.solve(matrix(colMeans(out.cmr.cusp$cov.inv),ncol = p)) ## stein estimator
      toc$cmr.cusp.intercept  = out.cmr.cusp$runtime
    }
    ####################################
    
    ####################################
    ## run intercept kron stuff if matrix-variate data
    if (cov.method == "kron"){
      # get kron MLE from data
      cov.kron.temp = cov.kron.mle(vec.inv.array(Y,p1,p2),
                                   de.meaned = T,my.seed = sim.ind + 500)
      output$kron.mle = cov.kron.temp$Cov
      toc$kron.mle = cov.kron.temp$runtime
      
      # # shrink to kron
      # kron.shrink.temp = ShrinkSep_GS(Y,p1,p2,
      #                                 S = S.simple,
      #                                 burnin = burnin.simple,
      #                                 my.seed = sim.ind + 500)
      # output$kron = qr.solve(matrix(colMeans(kron.shrink.temp$cov.inv),ncol = p)) ## stein estimator
      # toc$kron = kron.shrink.temp$runtime
    }
    ####################################
    
    ####################################
    ## sample covariance
    if ( p <= n){
      output$mle = cov.mle(Y)
      toc$mle = 0
    }
    ####################################


    ####################################
    ## run competitor GS- CUSP
    out.cusp = cusp_factor_adapt(Y,
                                 my_seed =  sim.ind + 600,
                                 N_sampl = S.fancy,
                                 burnin = burnin.fancy,
                                 alpha = DUNSON_ALPHA,
                                 a_sig = 1, b_sig = 1,
                                 a_theta = 2, b_theta = 2, theta_inf = 0.05,
                                 start_adapt = S.fancy, Hmax = p + 1, # don't adapt
                                 alpha0 = -1, alpha1 = -5*10^(-4))
    output$cusp = qr.solve(matrix(colMeans(out.cusp$cov.inv),ncol = p)) ## stein estimator
    toc$cusp = out.cusp$runtime
    ####################################

    ####################################
    ## run competitor GS- SIS
    out.sis = Mcmc_SIS(Y,X,
                       as = 1, bs = 1, alpha = DUNSON_ALPHA,
                       a_theta = 2, b_theta = 2,
                       nrun = S.fancy, burn = burnin.fancy,
                       thin = 1,
                       start_adapt = S.fancy, kmax = p + 1, # don't adapt
                       my_seed =  sim.ind + 700,
                       output = c("covSamples", "covSamplesInv","time"))
    output$sis = qr.solve(matrix(colMeans(out.sis$covSamplesInv),ncol = p)) ## stein estimator #
    # output$sis =  out.sis$covMean
    toc$sis = out.sis$time
    ####################################
    
    ####################################
    ## get distance from each output and the truth
    loss = unlist(lapply(output,function(k)
      loss_stein(k,true.cov)))
    ####################################
    
    ####################################
    ## save output
    output = list("loss" = loss,
                  "toc" = unlist(toc) #as.numeric(toc.swag,units="secs")
                  )
    ####################################
    
  }
  
  #stop cluster
  stopCluster(cl)
  
  temp.loss = Reduce("+",parallel.out[1,])/sim
  loss.avg = rbind(loss.avg,temp.loss)
  toc.avg = rbind(toc.avg,Reduce("+",parallel.out[2,])/sim)

  print(paste0("Finished running the scenario for N = ", n,"!!!!!!!!!"))

}

colnames(loss.avg) = colnames(toc.avg) = names(parallel.out[1,][[1]])
rownames(loss.avg) = rownames(toc.avg) = paste0("N",Ns.names)

####################################


print("Saving output now !!!")


####################################
## save output
output.filename = paste0("./output/SIM_",suffix,".Rdata")
save(loss.avg,toc.avg,file = output.filename)
####################################

print(paste0("Saved output to: ", output.filename))

apply(loss.avg,1,function(j)j/min(j)) |> t()
