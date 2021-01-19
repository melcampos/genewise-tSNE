#Bash script for vcftools and plotting or missing data

vcf_directory=$1        #/home/bigdata/chro_3R/exon_only
file_of_gene_ids=$2     #/home/bigdata/exon_only/chro_3R/filtered_SNP_gene_ids.txt
output_directory=$3     #/home/bigdata/chro_3R/exon_only/missingdata_SNP_genewise

gene_ids=$(cat $file_of_gene_ids)

mkdir -p $output_directory
echo GENE_ID$'\t'CHR$'\t'POS$'\t'N_DATA$'\t'N_GENOTYPE_FILTERED$'\t'N_MISS$'\t'F_MISS >$output_directory/$"genewise_maxfreqSNP_missingdata.tsv"
>$output_directory/$"genewise_indfreq_missingdata.tsv"
>$output_directory/$"genewise_indnum_missingdata.tsv"
>$output_directory/$"genewise_SNPnumbers.tsv"
>$output_directory/$"SNPfreq_missingdata.txt"

for gene_id in $gene_ids
do
	input_file=$vcf_directory/$gene_id.vcf
	output_file_prefix=$output_directory/$gene_id

  #Gene_wise missing data SNP frequencies
	vcftools --vcf $input_file --missing-site --out $output_file_prefix
	lmiss_file=$output_file_prefix$".lmiss"
  
  #Extract max SNP frequency
  echo "$(tail -n +2  $lmiss_file)" | sort -nrk6,6 -nrk2,2 - | head -1 | sed -e 's/^/'$gene_id'\t/' >> $output_directory/$"genewise_maxfreqSNP_missingdata.tsv"

  #Extract all SNP frequencies
  cat $lmiss_file | awk '{if(NR>1)print}' | awk '{print $6}' | sed -e 's/^/'$gene_id'\t/'>> $output_directory/$"SNPfreq_missingdata.txt"



  #Gene_wise missing data Individual frequencies
  vcftools --vcf $input_file --missing-indv --out $output_file_prefix
  imiss_file=$output_file_prefix$".imiss"
  
  #Extract frequency data for a single gene and input into a file for all genes 
  awk '{print $5}' $imiss_file | sed "1 s/.*/$gene_id/" - | paste $output_directory/$"genewise_indfreq_missingdata.tsv" - > $output_directory/$"temp.imiss"
  mv $output_directory/$"temp.imiss" $output_directory/$"genewise_indfreq_missingdata.tsv"
  
  #Extract number of missing entries for a single gene and input into a file for all genes 
  awk '{print $4}' $imiss_file | sed "1 s/.*/$gene_id/" - | paste $output_directory/$"genewise_indnum_missingdata.tsv" - > $output_directory/$"temp.imiss"
  mv $output_directory/$"temp.imiss" $output_directory/$"genewise_indnum_missingdata.tsv"
  
  #Extract genewise SNP numbers
  sed -n 2p $imiss_file | awk '{print $2}' - | sed -e 's/^/'$gene_id'\t/' >> $output_directory/$"genewise_SNPnumbers.tsv"


done
