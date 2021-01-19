# PCA for each gene
## Load the programs
library('SNPRelate')
library('gdsfmt')
library('ade4')
library('adegenet')
library('ggplot2')

# take the commandline arguments
args <- commandArgs(trailingOnly = TRUE)

# die if not enough args
if (length(args) < 2) {
  stop("not enough command line arguments")
}

# first argument = INPUT FILENAME (full path from current directory)
vcf.fn <- args[1]


# second argument = OUTPUT FILE PREFIX (e.g. output/foo/AGAP001234)
# all output filenames will start with this prefix
output.prefix <- args[2]

# Loading the VCF file
# create a temporary filename for the gds intermediate file
# this uses the Linux process ID of the running script itself so that
# we can run it in parallel without each process fighting over the same file
gds.fn <- paste("temp", Sys.getpid(), "gds", sep=".")

# Reformating to GDS : two method options "biallelic-only" or "copy.num.of.ref"
snpgdsVCF2GDS(vcf.fn, gds.fn, method="biallelic.only")
snpgdsSummary(gds.fn)

# Data Analysis
# load gds file into memory
genofile <- snpgdsOpen(gds.fn)
print("read genofile")


# clean up temporary file
file.remove(gds.fn)

# Converting the dataset
genotype_snp <- read.gdsn(index.gdsn(genofile, "genotype"))
write.table(genotype_snp, paste(output.prefix, "output", "txt", sep="."), sep="\t")
print("read genotype")

# Principal Component Analysis
pca <- prcomp(genotype_snp)
tab <- pca$x[,1:2]
metadata <- cbind.data.frame(tab, pop_code)
print("pca done")

# Write txt file with variance proportion (%)
pc.percent <- (pca$sdev^2 / sum(pca$sdev^2))*100
write.table(pc.percent, paste(output.prefix, "pc-percent", "txt", sep="."), sep="\t")

# Make a PC1vsPC2 plot
pdf(paste(output.prefix, "ggplot", "pdf", sep="."))
ggplot(data=metadata, aes(x=PC1, y=PC2, color=metadata$pop_code)) + geom_point(size=2) + 
  labs(x=format(pc.percent[1], digits=2), y=format(pc.percent[2], digits=2)) + 
  labs(colour = "Population")

# Write txt file for PCs %
tab2 <- pca$x[,1:4]
write.table(tab2, paste(output.prefix, "pcs1-4", "txt", sep="."), sep="\t")
