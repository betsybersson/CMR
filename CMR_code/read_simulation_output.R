library(xtable)
cov.method = "kron"  #"comSym3groups"
ps = c("9","16") #,"64"


loss.out = c()
for ( p.ind in 1:length(ps)){
  file.name = paste0("output/SIM_",cov.method,"_p",ps[p.ind],".Rdata")
  load(file.name)
  loss.temp = t(round(apply(loss.avg,1,function(j)j/min(j)),2))
  rownames(loss.temp) = paste0("p = ",ps[p.ind],
                               ", N = ",
                               substring(rownames(loss.temp),2), "p")
  
  loss.out = rbind(loss.out, 
                   loss.temp)
}

if (cov.method == "comSym3groups" | cov.method == "kron"){
  xtable(loss.out[,1:which(colnames(loss.out)=="lfm.K7")])
  xtable(loss.out[,-c(1:which(colnames(loss.out)=="lfm.K7"))])
} else {
  xtable(loss.out[,1:5])
  xtable(loss.out[,-c(1:5)])
}

