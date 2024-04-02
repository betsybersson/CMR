source("functions.R")
library(parallel)
library(doParallel)
####################################

####################################
## helpers
on.server = TRUE
cov.method = "eye" ## options: eye, cor9, comSym3groups
identifier = "p8"
####################################

####################################
## problem dimension parameters

# number of variables
p = 8
# sample sizes to loop through
Ns = p+2 # c(round(p/2),p+2,round(p*1.5),p*3)
Ns.names = "1" #c("0.5","1","1.5","3")
# low dimension
Ks = c(1,3,5,7)
####################################

print(paste0("Running the following scenario: ",
             "cov: ", cov.method,
             "; p: ", p,
             "; Ns: (", paste0(Ns, collapse = ","),
             "); Ks: (", paste0(Ks, collapse = ","),
             ") !!!!!!!!!!!!"))

####################################
## gibbs sampler variables
S = 600# 10000
burnin = 500
# number of simulation replicates
sim = 2 #25
####################################

####################################
## have correct path to read in R functions
if (on.server == T){
  source("betsy_R_functions.R")
} else {
  source("/Users/betsybersson/Documents/betsy_R_functions.R")
}
####################################


####################################
## true covariance matrix
if (cov.method == "eye"){
  true.cov = eye(p)
} else if (cov.method == "cor9"){
  true.cov = cor.mat(0.9)
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
loss.avg = toc.avg = matrix(NA,ncol = length(Ks)*2 + 2, nrow = length(Ns))
for ( n.ind in 1:length(Ns) ){
    
  n = Ns[n.ind]

  if (on.server == T){
    cl <- makeCluster(cores[1])  
  } else {
    cl <- makeCluster(cores[1] - 1)  # dont overload your computer
  }
  registerDoParallel(cl)
  
  ###
  parallel.out <- foreach(sim.ind=1:sim, .combine=cbind, .packages = c("parallel","LaplacesDemon")) %dopar% {
    
    ## output storage
    output = toc = list()

    ####################################
    ## sample dataset
    set.seed(sim.ind)
    Y = rmatnorm(matrix(0,ncol = p, nrow = n),V=true.cov)
    set.seed(Sys.time())
    ####################################
    
    
    ####################################
    ## run CMR GS
    for ( k.ind in 1:length(Ks) ){
      out.cmr  = CMR_GS(Y,
                        k = Ks[k.ind],
                        S = S,
                        burnin = burnin,
                        my.seed = sim.ind + 100)
      output[[k.ind]] = qr.solve(matrix(colMeans(out.cmr$cov.inv),ncol = p)) ## stein estimator
      toc[[k.ind]]  = out.cmr$runtime
    }
    # name based on value of K
    names(output) = names(toc) = paste0("cmr.K",Ks)
    ####################################
    
    ####################################
    ## run LFM GS
    for ( k.ind in 1:length(Ks) ){
      out.lfm  = LFM_GS(Y,
                        k = Ks[k.ind],
                        S = S,
                        burnin = burnin,
                        my.seed = sim.ind + 200)
      output[[k.ind + length(Ks)]] = qr.solve(matrix(colMeans(out.lfm$cov.inv),ncol = p)) ## stein estimator
      toc[[k.ind + length(Ks)]]  = out.lfm$runtime
    }
    # name based on value of K
    names(output)[-c(1:length(Ks))] = names(toc)[-c(1:length(Ks))] = paste0("lfm.K",Ks)
    ####################################
    
    ####################################
    ## sample covariance
    output$mle = cov(Y)
    toc$mle = 0
    ####################################
    
    ####################################
    ## run competitor GS- CUSP
    out.cusp = cusp_factor_adapt(Y,
                                 my_seed =  sim.ind + 300,
                                 N_sampl = S,
                                 alpha = 5, a_sig = 1, b_sig = 0.3,
                                 a_theta = 2, b_theta = 2, theta_inf = 0.05,
                                 start_adapt = S, Hmax = p + 1, # don't adapt
                                 alpha0 = -1, alpha1 = -5*10^(-4))
    output$cusp = qr.solve(matrix(colMeans(out.cusp$cov.inv),ncol = p)) ## stein estimator
    toc$cusp = out.cusp$runtime
    ####################################
    
    ###########################
    ## get distance from each output and the truth
    loss = unlist(lapply(output,function(k)
      loss_stein(k,true.cov)))
    ###########################
    
    output = list("loss" = loss,
                  "toc" = unlist(toc) #as.numeric(toc.swag,units="secs")
                  )
  }
  
  #stop cluster
  stopCluster(cl)
  
  temp.loss = Reduce("+",parallel.out[1,])/S
  loss.avg[n.ind,] = temp.loss/min(temp.loss)
  toc.avg[n.ind,] = Reduce("+",parallel.out[2,])/S

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