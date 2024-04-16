source("./DA_codes/cov_functions.R")
source("./DA_codes/audio_helpers.R")
source("./DA_codes/discrim_anal_functions.R")
library(Rcpp)
library(parallel)
sourceCpp("./DA_codes/fast_matrix_ops.cpp")
### get data and params
## read data
ddir<-"./word_matrices/" 
df = readRDS(paste0(ddir,"speech-MFCCs.rds"))
word.names = names(df)

## dimensions
ns = unlist(lapply(df,function(j)dim(j)[1]))
p1 = dim(df[[1]])[2]
p2 = dim(df[[1]])[3]

p = p1*p2
J = length(df)

n.test = 100
set.seed(123)
tosave = sapply(ns,function(j)sample(j,n.test))
set.seed(Sys.time())

Y.list = Y.test = list()
means = scale.attributes = list()

for( j in 1:J){
  
  df.temp = vec.3d(df[[j]][-tosave[,j],,])
  temp = scale(df.temp); scale.attributes[[j]] = attributes(temp)
  Y.list[[j]] = as.matrix(temp)
  means[[j]] = scale.attributes[[j]]$`scaled:center`
  
  Y.test[[j]] = vec.3d(df[[j]][tosave[,j],,])
}

# # save y.list (training data) for use in julia
# for ( j in 1:J){
#   write.table(Y.list[[j]],
#               file = paste0(HOMEDIR,"SWAG/SWAG-Bryant/word_matrices/",word.names[j],".txt"),
#               row.names = F,col.names = F)
# }

### get test data
ns = unlist(lapply(Y.list,nrow))
group = rep(1:J,times=ns)
group.compare = rep(1:J,each = n.test)

Y.da = c()
Sig.true= list()
for( j in 1:J){
  da.temp = Y.test[[j]]
  Y.da = rbind(Y.da,da.temp)
  Sig.true[[j]] = cov(scale(da.temp))
}


# update
g = J


###########-------------------------------
output = list()

## params
###### edit this part if running for more than 1 word
words = "no"
k = 5

## read in data
word = words[1]
out.filename = paste0("./output/output_CMR_",word,"_k",k)
out = readRDS(paste0(out.filename,"_COVHAT.RDS")) #".RDS"
output$cov.hat.stein = list() #qr.solve(matrix(out$cov.inv.mean,ncol = p)) ## stein estimator
output$cov.hat.stein[[1]] = out

###########------------------------------
# rescale mles
for ( j.ind in 1:length(words)){ ## loop through words 1:g
  
  j = which(word.names == words[j.ind])
  
  R = scale.attributes[[j]]$`scaled:scale`
  R = diag(R)
  
  output$cov.hat.stein[[j.ind]] = t(R) %*% output[[1]][[j.ind]] %*% R #multAprimeBA(R,output[[l]][[j]])
     
}
###########-------------------------------

print("getting covs done")

###########################
## DA Comparison
###########################

###########################
## reformat each output from 3d array to list
cov.hat.stein.list = list()
for(ZZ in 1:length(words)){
  cov.hat.stein.list[[ZZ]] = output$cov.hat.stein[[ZZ]]
}

###########################
## run DA with these estimates
DA = list()

DA$CMR = unlist(mclapply(1:nrow(Y.da),function(K)
  #da.function.1word(Y.da[K,],means[[which(word.names==word)]],cov.hat.stein.list[[1]])))
  da.function.1word.manymu(Y.da[K,],means,cov.hat.stein.list[[1]])))

print("DA done")


###########################

# saveRDS(DA,file="./output/audio_classrate_CMR_no_K5.RDS")


###########################
## get misclass rates
misclass = lapply(DA,function(ZZ)
  mean((ZZ-group.compare)==0))

missclass = unlist(misclass)
# temp[order(temp)]
round(missclass,2)

###########-------------------------------
# group specific misclass rates
group.misclass = lapply(DA,function(ZZ)
  tapply((ZZ-group.compare)==0,factor(group.compare),mean))

group.misclass = matrix(unlist(group.misclass),ncol = length(group.misclass),byrow = F)
colnames(group.misclass) = names(DA)
rownames(group.misclass) = paste0("group",1:g)
round(group.misclass,2)
###########-------------------------------



###########-------------------------------
# contrast tables

CA = array(NA, dim = c(J,J,1))
dimnames(CA)[[1]] = word.names ## target word
dimnames(CA)[[2]] = word.names ## predicted word
dimnames(CA)[[3]] = c("CMR")


for ( est in c(1:1) ){
  for ( target in 1:J ){
    for ( pred in 1:J ){
      CA[target,pred,est] = sum(DA[[est]][group.compare == target]==pred)
    }
  }
}

# plot
# pdf("./plots/audio_CA_99.pdf", family="Times",height=6)
par(mfrow=c(1,1),mar=c(3,1,1.35,1),mgp=c(1.75,.75,0),oma=c(.1,2,.1,.1) )
for ( j in 1:1){
  image(t(CA[J:1,,j]),xaxt="n",yaxt="n")
  axis(side=1,at=seq(0,1,length=10),
       labels=dimnames(CA)[[1]],tick=FALSE,las=3 ,cex.axis=1.25)  
  axis(side=2,at=seq(1,0,length=10),
       labels=dimnames(CA)[[1]],tick=FALSE,las=1 ,cex.axis=1.25)   
  mtext(side=3,dimnames(CA)[[3]][j],line=.2 )
}
# dev.off()









