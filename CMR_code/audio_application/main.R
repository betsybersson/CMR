source("../betsy_R_functions.R")
source("../functions.R")
library(parallel)

## for one word at a time
word = "down"
Y = read.table(paste0("./word_matrices/",word,".txt"),
               quote="\"", comment.char="")
Y = as.matrix(Y)

# sample size
n = nrow(Y)
# feature dimension
p = ncol(Y)

# use separable meta covariates
p1 = 99
p2 = 13
X = getMatDesignMat(p1,p2)

# latent dimension
n.element.uppertri = function(p){
  p*(p-1)/2
}
n.elem.kron = n.element.uppertri(p1) + n.element.uppertri(p2) ## num open params for kron
n.elem.kron / (p1+p2) ## 

k = 20

print(paste0("Running the following scenario: ",
             "word: ", word,
             "; Ks: (", paste0(k, collapse = ","),
             ") !!!!!!!!!!!!"))

# run GS
out = CMR_GS(Y,X,
             k = k,
             S = 5250,
             burnin = 250,
             my.seed = 1,
             save.all = "ICO")

print("Saving output now !!!!!!!!!!!!")

out.filename = paste0("./output/output_CMR_",word,"_k",k)
saveRDS(out,file = paste0(out.filename,".RDS"))

print("All done !!!!!!!!!!!!")


cov.hat.stein = qr.solve(matrix(out$cov.inv.mean,ncol = p)) ## stein estimator

saveRDS(cov.hat.stein,file = paste0(out.filename,"_COVHAT.RDS"))


# # plot trace plots of cov
# set.seed(1)
# plot.trace.ind = sample(1:(p^2),5)
# set.seed(Sys.time())
# pdf(paste0(out.filename,".pdf"),
#     family="Times",height = 7,width = 6.6)
# matplot(out$cov[,plot.trace.ind],type="l")
# dev.off()


