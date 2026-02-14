source("functions.R")
source("Schiavon_SISGaussian.R")
library(parallel)
library(doParallel)
####################################

####################################
## helpers
on.server = F
cov.method = "updatedGroup" ## options: continuous, eye, cor9, comSym3groups, kron, continuous
identifier = "feb"
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
S.fancy = 20#20000
burnin.fancy = 10#10000

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
} else if (cov.method == "comSym3groups" || cov.method == "continuous"){
  
  # ### original approach
  # p_l = round(c(4,6,8)/sum(4,6,8)*p)
  # L = length(p_l)
  # p = sum(p_l)
  # which_group = rep(1:L,times = p_l)
  # 
  # lam = c(1,0,1,1)
  # Lam = matrix(rep(lam,p_l[1]),nrow=p_l[1],byrow=T)
  # lam2 = c(0,0,0,1)
  # Lam2 = matrix(rep(lam2,p_l[2]),nrow=p_l[2],byrow=T)
  # lam3 = c(1,4,0,0)
  # Lam3 = matrix(rep(lam3,p_l[3]),nrow=p_l[3],byrow=T)
  # Lam = rbind(Lam,Lam2,Lam3)
  # cov = Lam %*% t(Lam)
  # cov = cov + eye(p)
  # true.cov = cov2cor(cov)
  # eigen(true.cov)$val ## check invertible
  
  # get block correlation value map
  val_map = matrix(c(
    0.75, 0.3535534, 0.1178511,  # Group 1 correlations (Within G1, G1-G2, G1-G3)
    0.3535534, 0.5, 0,  # Group 2 correlations (G2-G1, Within G2, G2-G3)
    0.1178511, 0, 0.9444444   # Group 3 correlations (G3-G1, G3-G2, Within G3)
  ), nrow = 3, byrow = TRUE)
  
  # get number of elements in group
  p_l = round(c(4,6,8)/sum(4,6,8)*p)
  L = length(p_l)
  p = sum(p_l)
  which_group = rep(1:L,times = p_l); group_sizes = table(which_group)
  true.cov = build_blocked_cov(group_sizes, val_map)
  eigen(true.cov)$val ## check invertible
  
  # get design matrix
  X = matrix(0, nrow = p, ncol = length(p_l))
  for ( pl.ind in 1:length(p_l)){
    X[which_group==pl.ind,pl.ind] = 1
  }
} else if (cov.method == "kron"){
  p1 = round(p/5)
  p2 = round(p/p1)
  p = p1*p2
  
  R = cor.mat(p1,.9)
  C = cor.mat(p2,.6)
  true.cov = kronecker(C,R)
  eigen(true.cov)$val ## check invertible
  
  # get design matrix
  X = getMatDesignMat(p1,p2)
  
} else if (cov.method == "updatedGroup"){
  val_map = matrix(c(
    0.75, -0.20,  0.14,  # Group 1 (Modified off-diagonals for safety)
    -0.20,  0.50,  0.60,  # Group 2 (Lowered G2-G3 from 0.7 to 0.4)
    0.14,  0.60,  0.94
  ), nrow = 3, byrow = TRUE)
  # get number of elements in group
  p_l = round(c(4,6,8)/sum(4,6,8)*p)
  L = length(p_l)
  p = sum(p_l)
  which_group = rep(1:L,times = p_l); group_sizes = table(which_group)
  true.cov = build_blocked_cov(group_sizes, val_map)
  eigen(true.cov)$val ## check invertible
  
  # get design matrix
  X = matrix(0, nrow = p, ncol = length(p_l))
  for ( pl.ind in 1:length(p_l)){
    X[which_group==pl.ind,pl.ind] = 1
  }
}
### if in continuous case, use true cov same as with three groups, but change the meta covariate used
if (cov.method == "continuous"){ 
  set.seed(123)
  X = matrix(rnorm(p,which_group,sd = 1/4),ncol=1)
  set.seed(Sys.time())
  X = cbind(1,X) # add an intercept
}
# propagate filename suffix
suffix = paste0(cov.method,"_p",p)
if (identifier != ""){
  suffix = paste0(suffix,"_",identifier)
}
####################################

####################################
## run simulation
if (on.server == T){
  cores = detectCores()
} else {
  cores = detectCores() - 1
}
loss.avg = toc.avg = c() #matrix(NA,ncol = length(Ks)*2 + 2, nrow = length(Ns))
loss.mat = list()
for ( n.ind in 1:length(Ns)){  #1:length(Ns)
  
  n = Ns[n.ind]
  
  if (on.server == T){
    cl <- makeCluster(cores[1])  
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
    output$MR.O = qr.solve(matrix(colMeans(out.cmr.cusp$cov.inv),ncol = p)) ## stein estimator
    toc$MR.O  = out.cmr.cusp$runtime
    ####################################
    
    ####################################
    ## run intercept CMR GS if you haven't
    if (!all(is.na(X))){

      ## run CMR GS with CUSP
      out.cmr.cusp  = CMR_cusp_GS(Y,X = NA,
                                  k = ceiling((p-1)/2),
                                  S = S.fancy,
                                  burnin = burnin.fancy,
                                  my.seed = sim.ind + 400,
                                  alpha = DUNSON_ALPHA,
                                  a.theta = 1/2, b.theta = 1/2)
      
      output$MR.I = qr.solve(matrix(colMeans(out.cmr.cusp$cov.inv),ncol = p)) ## stein estimator
      toc$MR.I  = out.cmr.cusp$runtime
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
    
    }
    ####################################
    
    ####################################
    ## sample covariance
    if ( p <= n){
      output$MLE = cov.mle(Y)
      toc$MLE = 0
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
    output$CUSP = qr.solve(matrix(colMeans(out.cusp$cov.inv),ncol = p)) ## stein estimator
    toc$CUSP = out.cusp$runtime
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
    output$SIS = qr.solve(matrix(colMeans(out.sis$covSamplesInv),ncol = p)) ## stein estimator #
    # output$sis =  out.sis$covMean
    toc$SIS = out.sis$time
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
    
    
    ####################################
    ## save output
    output.filename = paste0("./output_temp/SIM_",suffix,"_nInd",n,"_simInd",sim.ind,".Rdata")
    save(loss.avg,toc.avg,loss.mat,file = output.filename)
    ####################################
    
    
    cat(print(paste0("Finished running and saving:",
                 " Sim ind : ",sim.ind,
                 "; n ind: ", n,
                 "!!!!!!!!!!!!!")), file = "progress_log.txt", append = TRUE)
    
    
  }
  
  #stop cluster
  stopCluster(cl)
  
  print(paste0("Finished running and saving all of:",
               "n ind: ", n,
               "!!!!!!!!!!!!!"))
  
  
}



print(paste0("Finished everything !!!!!!!!!!!!!"))
