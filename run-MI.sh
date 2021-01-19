dist_directory=$1        #"/home/genome_dist"
file_of_gene_ids=$2     # "/home/genome_dist/gene_ids.txt"
output_directory=$3     #"/home/MoransI/output"
pop_code_file=$4        

# run the 'cat' command and put all the output into a new variable, called gene_ids
gene_ids=$(cat $file_of_gene_ids)

# make the output directory (won't complain if it's already there)
mkdir -p $output_directory

declare -a arr=("AOM" "BFM" "BFS" "CMS" "GAS" "GNS" "GWA" "KES" "UGS")
# this loops through the 'words' in $gene_ids variable

for pop in "${arr[@]}"
do
 output_file=$output_directory$"/MI_"$pop$".tsv"
 echo started $output_file ...
 # next line empties the output file
 >$output_file
 for gene_id in $gene_ids
  do
  txt_file=$dist_directory/${gene_id}.dist.txt
  gunzip $txt_file.gz

  
  # this is where you run the R script
  Rscript MI_pop.R $txt_file $pop_code_file $pop >> $output_file
  gzip $txt_file
  done
done

