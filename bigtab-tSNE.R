library(RColorBrewer)
library(tsne)
library(Rtsne)

# take the commandline arguments
args <- commandArgs(trailingOnly = TRUE)

# die if not enough args
if (length(args) < 2) {
  stop("not enough command line arguments")
}

tab_directory <- args[1]
output_directory <- args[2]

tab_tsne <- read.table("pcs1-50.txt", header=T,row.names=1)

tab_tsne2 <- as.matrix(tab_tsne)

tsne <- Rtsne(tab_tsne2,theta=0.0,check_duplicates=F,pca=FALSE,verbose=TRUE, perplexity=500, max_iter=5000)

tab3 <- tsne$Y[,1:2]
write.table(tab3, paste(output_directory,"tsne1-2_p500.txt", sep="/"))

pdf(paste(output_directory, "tsne_plot_p500.pdf",sep="/"), width=9, height=6)
plot(main="tsne",x=tsne$Y[,1], y=tsne$Y[,2], pch=19, cex=0.5, col="black")

print("t-SNE is complete - perplexity 500")
