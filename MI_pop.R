library(ape)

# take the commandline arguments
args <- commandArgs(trailingOnly = TRUE)

# die if not enough args
if (length(args) < 2) {
  stop("not enough command line arguments")
}

genewise_tsne_dist_file <- args[1] # "/home/genome_dist"
pop_code_file <- args[2] # "/home/pop_pop.txt"
pop <- args[3] #"AOM"

pop <- as.character(pop)

#take in genewise tsne distance matrix
genewise_tsne_dist <- read.csv(genewise_tsne_dist_file, sep="\t", stringsAsFactors = F,header=TRUE)
genewise_tsne_weight <- as.matrix(1/(genewise_tsne_dist))
genewise_tsne_weight[genewise_tsne_weight=="Inf"] <- 1
diag(genewise_tsne_weight) <- 0

#take in vector of individuals population code
pop_code <- scan(pop_code_file, what = "character")
pop_code_key <- unique(pop_code)

pop_code_vec <- (pop_code==pop)*1

MI <- unlist(as.vector(Moran.I(pop_code_vec,genewise_tsne_weight)))
cat(MI,"\n")
