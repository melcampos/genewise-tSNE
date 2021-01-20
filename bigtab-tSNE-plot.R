library(RColorBrewer)
library(ggplot2)

# Load the main file
file <- read.table("tsne1-2.txt", sep="\t", header=T)

# t-SNE plot with colours for each chromossomes and inversion regions. 

cols <- c("#e41a1c","#377eb8","#4daf4a","#984ea3","#ff7f00","gold2","#f781bf")
cols_t1<-cols[as.factor(file$chro1)]

plot(file$tsne1, file$tsne2,main="", col=alpha(cols_t1,0.5), cex=0.5, pch=16,xlab="",ylab="",axes=F)
