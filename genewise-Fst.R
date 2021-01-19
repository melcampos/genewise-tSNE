# Genewise FST for each gene
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

# Population Information (Metadata)
pop_code <- scan("pop_MS.txt", what=character())

print("read pop file")

# Converting the dataset
genotype_snp <- read.gdsn(index.gdsn(genofile, "genotype"))

# Fst calculation
sample.id <- read.gdsn(index.gdsn(genofile, "sample.id"))
flag <- pop_code %in% c("M","S")
samp.sel <- sample.id[flag]
pop.sel <- pop_code[flag]
fst <- snpgdsFst(genofile, sample.id=samp.sel, population=as.factor(pop.sel),
          method="W&C84", autosome.only=FALSE)
          
write.table(fst, paste(output.prefix, "fst_BFS", "txt", sep="."), sep="\t")
print("fst done")
