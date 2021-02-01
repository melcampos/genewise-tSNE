#!/bin/env perl
#              -*-  mode: Cperl -*-
#
# RUN THIS AFTER. /distance-stats-merge.pl
#
# prep: make distance files with only short-ish distances in them (<5)
#
# for file in $(cd data ; ls -1 tsne??.txt) ; do perl -nale 'print if ($F[3] < 5)' data/$file > filtered-data/$file ; done
#
# usage: cat gene-stats.txt | parallel --jobs 10 --pipe --N 100 ./add-local-distance-stats.pl --input - filtered-data/tsne??.txt > gene-stats-local.txt
#
# will add new columns - the number of different genes seen in the 30 t-SNEs in the closest 5, 20, 50 neighbours to the gene
#
# nun = "number of unique neighbours"
#
# no error checking for bad files etc!
#
# it's going to take a while..., though it could be parallelised a bit
# it reads in a gene list, and then reads through all the files for each gene,
# pulling out distances for that gene, sorting them, doing the overlap analysis
#
# The final output may need column headings edited by hand
#

use strict;
use warnings;
use PDL;
use PDL::Stats::Basic;
use Getopt::Long;

$| = 1;
my $input;

GetOptions("input=s"=>\$input);

my @cutoffs = (5, 20, 50);
# output headings
# mean_variance	median_variance	nun5	nun20	nun50

my @files = @ARGV;

my %input_lines; # gene => line without newline

open(INPUT, $input) || die;
#my $header = <INPUT>;
#chomp($header);
while (<INPUT>) {
  chomp;
  next if (/mean/); # skip header line
  my ($gene_id) = split "\t", $_;
  $input_lines{$gene_id} = $_;
}
close(INPUT);


# NO HEADERS - running it in parallel
# new output header
# print join("\t", $header, map { "nun$_" } @cutoffs)."\n";

# note: the gene_id includes the double quotes at all times
foreach my $gene_id (sort keys %input_lines) {

  my %closest_genes; # 20 => gene_id => 1
  # where 20 means the closest 20 genes

  foreach my $file (@files) {
    my %distances;  # other_gene_id => distance
    open(FILE, $file) || die;
    my $ignore_header = <FILE>;
    while (<FILE>) {
      chomp;
      my ($ignore_index, $gene1, $gene2, $distance) = split "\t", $_;
      if ($gene1 eq $gene_id) {
	$distances{$gene2} = $distance;
      } elsif ($gene2 eq $gene_id) {
	$distances{$gene1} = $distance;
      }
    }
    my @sorted_genes = sort { $distances{$a} <=> $distances{$b} } keys %distances;

    foreach my $n (@cutoffs) {
      for (my $i=0; $i<$n; $i++) {
	if (defined $sorted_genes[$i]) {
	  $closest_genes{$n}{$sorted_genes[$i]} = 1;
	} else {
	  warn "not enough distances for $gene_id at cutoff $n\n";
	  last;
	}
      }
    }
    close(FILE);
  }

  print $input_lines{$gene_id};
  foreach my $n (@cutoffs) {
    print "\t".scalar(keys %{$closest_genes{$n}});
  }
  print "\n";
}
