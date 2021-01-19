#take the commandline arguments
args <- commandArgs(trailingOnly = TRUE)

# die if not enough args
if (length(args) < 2) {
  stop("not enough command line arguments")
}

# first argument = INPUT FILENAME (full path from current directory)
data <- args[1]

# second argument = OUTPUT FILE PREFIX (e.g. output/foo/AGAP001234)
# all output filenames will start with this prefix
output.prefix <- args[2]
print("Entendi")

library("MASS")
library("lsa")

# Distance
print("Distance Matrix...")
mydata <- read.csv(data, sep="\t", row.names=1, header=TRUE, stringsAsFactors = F)
print("Did")

mydata <- as.matrix(mydata)
print("Reading")

d <- dist(mydata, method="manhattan")
print("Writing the table")

mat <- as.matrix(d)
write.table(mat, paste(output.prefix, "dist", "txt", sep="."), sep="\t")
