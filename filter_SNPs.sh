vcf_directory=$1    # e.g./home/bigdata/chro_3R/original_vcf
file_of_gene_ids=$2 # e.g./home/bigdata/chro_3R/gene_id.txt
output_directory=$3 # e.g./home/bigdata/chro_3R/no_intron

gene_ids=$(cat $file_of_gene_ids)
mkdir -p $output_directory

for gene_id in $gene_ids
do
	vcf_file=$vcf_directory/$gene_id$".vcf"
	
	#For combined filtering steps
  cat $vcf_file | grep -v "intron_variant" | grep -v "intragenic_variant" | grep -v "intergenic_region" | grep -v "3_prime_UTR_variant" | grep -v "5_prime_UTR_variant"  > $output_directory/$gene_id$".vcf"

	#For individual filtering steps
	#cat $vcf_file | grep -v "intron_variant" > $output_directory/$gene_id$".vcf"
  #cat $vcf_file | grep -v "intergenic_region" > $output_directory/$gene_id$".vcf"
  #cat $vcf_file | grep -v "3_prime_UTR_variant" > $output_directory/$gene_id$".vcf"
  #cat $vcf_file | grep -v "intragenic_variant" > $output_directory/$gene_id$".vcf"
  #cat $vcf_file | grep -v "5_prime_UTR_variant" > $output_directory/$gene_id$".vcf"
  
  done
