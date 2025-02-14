library(xtable)  
library(dplyr)
library(ggplot2)

cov.method = "cor9"  # "kron"  # "cor9"# "comSym3groups" 
plottitle = expression("(a) " ~ Sigma[cor]) # expression(Sigma[kron]) # expression("(b) " ~ Sigma[block]) #  
legendposit =  ""  # "bottom" #
plotidentifier = ""

ps = c("9","16","50") 
if (cov.method == "kron"){
  ps = c("8","15","50")
}

grand.df = data.frame()
for ( p.ind in 1:length(ps)){
  
  file.name = paste0("output/SIM_",cov.method,"_p",ps[p.ind],"_saveall.Rdata")
  load(file.name)
  
  loss.mat.stand = lapply(loss.mat,log) #lapply(loss.mat,function(j)j/min(colMedians(j)))
  ssnames = c("n = p+1","n = 1.5p", "n = 3p")
  names(loss.mat.stand) = ssnames
  # Convert each matrix to a data frame with row, col, and value
  
  data_frame <- bind_rows(
    lapply(names(loss.mat.stand), function(name) {
      as.data.frame(as.table(loss.mat.stand[[name]])) %>%
        mutate(matrix_name = name)
    })
  )[,-1] 
  colnames(data_frame) = c("Method","Loss","SampleSize")
  data_frame$SampleSize = factor(data_frame$SampleSize,levels = ssnames)
  data_frame$Dimension = paste0("p = ",ps[p.ind])
  
  data_frame = data_frame |> 
    mutate(Method = if_else(Method == "kron.mle","Kron",Method)) |> 
    mutate(Method = if_else(Method == "MR.O","MR.D",Method))
  if (cov.method == "cor9"){
    data_frame = data_frame |> 
      mutate(Method = if_else(Method == "MR.C","MR.I",Method))
  }
  
  grand.df = rbind(grand.df,data_frame)
  
  if (cov.method == "comSym3groups"){
    file.name = paste0("output/SIM_continuous_p",ps[p.ind],"_saveall.Rdata")
    load(file.name)
    
    loss.mat.stand = lapply(loss.mat,log) #lapply(loss.mat,function(j)j/min(colMedians(j)))
    ssnames = c("n = p+1","n = 1.5p", "n = 3p")
    names(loss.mat.stand) = ssnames
    # Convert each matrix to a data frame with row, col, and value
    
    data_frame <- bind_rows(
      lapply(names(loss.mat.stand), function(name) {
        as.data.frame(as.table(loss.mat.stand[[name]])) %>%
          mutate(matrix_name = name)
      })
    )[,-1] 
    colnames(data_frame) = c("Method","Loss","SampleSize")
    data_frame$SampleSize = factor(data_frame$SampleSize,levels = ssnames)
    data_frame$Dimension = paste0("p = ",ps[p.ind])
    
    data_frame = data_frame |> 
      mutate(Method = if_else(Method == "kron.mle","Kron",Method)) |> 
      mutate(Method = if_else(Method == "MR.O","MR.C",Method)) |> 
      filter(Method == "MR.C")
    
    grand.df = rbind(grand.df,data_frame)
  }
  
}


grand.df$Method = factor(grand.df$Method, 
                         levels = c("MR.D","MR.C","MR.I","MLE","Kron","CUSP","SIS")) 

grand.df$Dimension = factor(grand.df$Dimension,levels=paste0("p = ",ps))

if (cov.method == "cor9"){
  pdf(paste0("./plots/Loss_boxplots_",cov.method,plotidentifier,".pdf"), family="Times",
      height=7.6,width=8) 
} else {
  pdf(paste0("./plots/Loss_boxplots_",cov.method,plotidentifier,".pdf"), family="Times",
      height=8,width=8) 
}
ggplot(grand.df) +
  geom_boxplot(aes(y=Loss,color=Method)) +
  facet_grid(Dimension~SampleSize, scales = "free_y") +
  theme_minimal(base_size = 14) +
  theme(
    panel.grid.major = element_blank(),  # Remove major gridlines
    panel.grid.minor = element_blank(),  # Remove minor gridlines
    panel.background = element_blank(), # Ensure white background
    axis.line = element_line(color = "black"), # Add axis lines
    axis.text.x = element_blank(),
    strip.text = element_text(face = "bold"), # Bold facet labels
    panel.spacing = unit(.05, "lines"),
    panel.border = element_rect(color = "black", fill = NA, size = 1), 
    strip.background = element_rect(color = "black", size = 1),
    legend.position = legendposit,
    legend.title = element_blank(),
    text = element_text(size = 18, family = "Times")) +
  ylab("Log Loss") + 
  ggtitle(plottitle) +
  scale_color_manual(values = c("MR.D" = "red","MR.C" = "pink",
                                "MLE" = "blue","CUSP" = "darkgreen",
                                "SIS" = "purple","MR.I" = "orange",
                                "Kron" = "black"))
dev.off()

