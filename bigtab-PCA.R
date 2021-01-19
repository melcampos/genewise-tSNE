library("bigmemory.sri")
library("bigmemory")
library("bigalgebra")
library("foreach")
library("DBI")
library("biglm")
library("biganalytics")
library("NCmisc")
library("reader")
library("bigpca")
library(ggplot2)

# take the commandline arguments
args <- commandArgs(trailingOnly = TRUE)

# die if not enough args
if (length(args) < 2) {
  stop("not enough command line arguments")
}

tab_directory <- args[1]
output_directory <- args[2]


gene_vectors <- read.table("huge-table.txt", sep=",")
head(gene_vectors)
dim(gene_vectors)

gene_vectors <- as.big.matrix(gene_vectors)

#gene_vectors <- big.t(gene_vectors)
pc <- big.PCA(gene_vectors, SVD = F)

print("big.PCA finished")

print("Saving it...")

# output 1:
pc.percent <- (pc$Evalues^2 / sum(pc$Evalues^2))*100

print("Saving pc plot")

#pop_code <- read.table("pop.txt", header=T, sep="\t")

# output 2:
pdf(paste(output_directory, "big_pca-big.pdf", sep="/"), width=9, height=6)
plot(pc$PCs[,1], pc$PCs[,2], pch=19, cex=0.5, col="black", main="", 
	xlab= format(pc.percent[1], digits=3), ylab= format(pc.percent[2], digits=3))

# output 3:
tab <- pc$PCs[,1:4]
write.table(tab, paste(output_directory,"pcs_1-4.txt", sep="/"))

# output 3:
tab <- pc$PCs[,1:50]
write.table(tab, paste(output_directory,"pcs_1-50.txt", sep="/"))
