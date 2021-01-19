txt_directory=$1
file_of_gene_ids=$2
output_directory=$3

# run the 'cat' command and put all the output into a new variable, called gene_ids
gene_ids=$(cat $file_of_gene_ids)

# make the output directory (won't complain if it's already there)
mkdir -p $output_directory

# this loops through the 'words' in $gene_ids variable 
for gene_id in $gene_ids
do

 txt_file=$txt_directory/$gene_id.output.txt

 # this is where you run the R script
 # wc -l $txt_file > $output_directory/$gene_id.vcf_length.txt
 module add R/3.2.1
 Rscript genewise-dist.R  $txt_file $output_directory/$gene_id
done

