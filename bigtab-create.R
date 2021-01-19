library(ggplot2)
library(gdsfmt)
library(SNPRelate)
library(bigmemory.sri)
library(bigmemory)
library(bigalgebra)
library(foreach)
library(DBI)
library(biglm)
library(biganalytics)
library(NCmisc)
library(reader)
library(bigpca)


# take the commandline arguments
args <- commandArgs(trailingOnly = TRUE)

# die if not enough args
if (length(args) < 2) {
  stop("not enough command line arguments")
}

tab_directory <- args[1]
output_directory <- args[2]

file_pattern <- paste(tab_directory, "*.dist.txt", sep="/")
files <- Sys.glob(file_pattern)
nfiles <- length(files)
gene_vectors <- big.matrix((765*765-765)/2, nfiles, init=0)
print("It is going well")
head(gene_vectors[1:5,1:8])
for (i in 1:nfiles) {
  file <- files[i];
  print(paste("reading ", file));
  dist_mat = as.matrix(read.table(file));
  diag(dist_mat) = NA;
  dist_mat[lower.tri(dist_mat)] = NA;
  gene_vectors[,i] <- dist_mat[which(!is.na(dist_mat))];
  }
print("Almost there")

write.big.matrix(gene_vectors, paste(output_directory,"huge-table.txt", sep="/"))                      
head(gene_vectors[1:5,1:8])
