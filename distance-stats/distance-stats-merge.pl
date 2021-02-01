#!/bin/env perl

#
# usage: ./distance-stats-merge.pl --gene-wise-output gene-stats.txt data/tsne??.txt > distance-stats.txt
#
# no error checking for bad files etc!
#

use strict;
use warnings;
use PDL;
use PDL::Stats::Basic;
use Getopt::Long;

my $gene_output;

GetOptions("gene-wise-output=s"=>\$gene_output);

my @files = @ARGV;

my %fh; #file handles key is filename

grep { open $fh{$_}, "<$_" } @files;

grep { my $header = readline($fh{$_}) } @files;

print join("\t", '"Gene1"', '"Gene2"', '"dist_mean"', '"dist_variance"', '"dist_n"')."\n";

my $gene_fh;
open ($gene_fh, ">$gene_output") if ($gene_output);

my %gene_variances; # gene_id => [ variances ]

my $row = 1;
while (1) {
  my %distances; #  gene_id_1 => gene_id_2 => [ distances ]

  grep { my $line = readline($fh{$_});
	 if ($line) {
	   chomp($line);
	   my @F = split "\t", $line;
	   push @{$distances{$F[1]}{$F[2]}}, $F[3];
	 }
       } @files;

  last unless (keys %distances);

  # check that there is only one key at each level of distances
  # because gene1 and gene2 should be the same for the line read from each file

  if (keys %distances == 1) {
    my ($gene1) = keys %distances;
    if (keys %{$distances{$gene1}} == 1) {
      my ($gene2) = keys %{$distances{$gene1}};
      my $pdl = pdl(@{$distances{$gene1}{$gene2}});
      my $variance = $pdl->var;
      printf "\"%d\"\t%s\t%s\t%f\t%f\t%d\n",
	$row++, $gene1, $gene2, $pdl->average, $variance, $pdl->nelem;

      if ($gene_output) {
	push @{$gene_variances{$gene1}}, $variance;
	push @{$gene_variances{$gene2}}, $variance;
      }
    }
  }
}

if ($gene_output) {
  print $gene_fh join("\t", "mean_variance", "median_variance")."\n";
  foreach my $gene (sort keys %gene_variances) {
    my $pdl = pdl(@{$gene_variances{$gene}});
    print $gene_fh join("\t", $gene, $pdl->average, $pdl->median)."\n";
  }
  close($gene_fh);
}
