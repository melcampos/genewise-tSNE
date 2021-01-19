# genewise-tSNE

Unravelling population structure heterogeneity within the genome of the malaria vector Anopheles gambiae

## Variant calling files preparation
- `gene-extractor.pl` : Extract individual VCF files for each gene or window;
- `filter_SNPs.sh` : Select exonic SNPs;

## Gene by gene PCA plotting, distance matrix and other metrics
- `metrics.sh` & `metrics2.sh` : Check VCF metrics such as: missing data, SNP frequency, number of SNPs, nucleotide diversity, etc.
- `genewise-Fst.R` : Calculate Fst between M and S forms for each gene;
- `genewise-PCA.R` : Perform a principal component analysis for each gene and convert to VCF to GDS;
- `genewise-dist.R` : Perform distance matrix for each gene;

## Linearize matrices and big table format
- `bigtab-create.R` : Create the big table with all linearized values for each gene;
- `bigtab-PCA.R` : PCA of big table created;
- `bigtab-tSNE.R`: tSNE of PCs. 

