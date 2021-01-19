Bash script for vcftools and plotting or nucleotide diversity

vcf_directory=$1        #/home/bigdata/chro_3R/exon_only
file_of_gene_ids=$2     #/home/bigdata/exon_only/chro_3R/filtered_SNP_gene_ids.txt
output_directory=$3     #/home/bigdata/chro_3R/exon_only/missingdata_SNP_genewise

gene_ids=$(cat $file_of_gene_ids)

mkdir -p $output_directory
>$output_directory/$"mean_site-pi.txt"

for gene_id in $gene_ids
do
        input_file=$vcf_directory/$gene_id.vcf
        output_file_prefix=$output_directory/$gene_id

  #Gene_wise nucleotide diversity 
        vcftools --vcf $input_file --keep vcftool-pop-lin_lin.pop --site-pi --out $output_file_prefix
        site_file=$output_file_prefix$".sites.pi"

  #Extract all SNP frequencies
        cat $site_file | awk -F '\t' '{sum += $3} END {print sum /NR }' | sed -e 's/^/'$gene_id'\t/'>> $output_directory/$"mean_site-pi.txt"


done
