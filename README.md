# genewise-tSNE

Unravelling population structure heterogeneity within the genome of the malaria vector Anopheles gambiae

## Variant calling files preparation
- `gene-extractor.pl` : Extract individual VCF files for each gene or window;
- `filter_SNPs.sh` : Select exonic SNPs;

## Gene by gene PCA plotting, distance matrix and other metrics
- `genewise-Fst.R` : Calculate Fst between M and S forms for each gene;
- `hfst_se_by_group.csv` (results) : Standard error estimates of `hfst_by_group.csv`
- `pi_by_group.csv` (results) : Sequence diversity (&pi;) within each sample group

## Linearize matrices and big table format
- `filter_SNPs.sh` : select exonic SNPs
- `hfst_se_by_group.csv` (results) : Standard error estimates of `hfst_by_group.csv`
- `pi_by_group.csv` (results) : Sequence diversity (&pi;) within each sample group

