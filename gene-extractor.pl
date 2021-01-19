#!/usr/bin/env perl
use v5.14.00;
use strict;
use warnings;
use Carp;
#use autodie qw(:all);
#use IPC::System::Simple qw(system);
use Readonly;
use Getopt::Long qw(:config no_ignore_case);
use Log::Log4perl qw( :easy ); 
Log::Log4perl->easy_init($WARN); 
my $logger = get_logger(); 
use Try::Tiny;

#use Bio::EnsEMBL::Registry;

###############################################################################
# MAIN
# Get command line args
my %opt = %{ opt_check() };

my $window_size = $opt{window};
my $vcf = $opt{vcf};
my $output_dir = $opt{output_dir};
my $gtf_file = $opt{gtf_file};

my $contigs = get_contigs($vcf);

# Create windows
if ($window_size) {
  $window_size *= 1000;
  my @windows = create_windows($contigs, $window_size);

  # Create as many files
  create_vcf_subfiles($vcf, \@windows, $output_dir);
}
elsif ($gtf_file) {
  my @genes_coords = get_genes_coords_from_gtf($gtf_file, $contigs);

  say "Genes: " . scalar(@genes_coords);

  # Create as many files
  create_vcf_subfiles($vcf, \@genes_coords, $output_dir);
}

sub get_genes_coords_from_gtf {
  my ($gtf, $contigs) = @_;
  
  open my $gtfh, "<", $gtf;
  my @genes;
  while (my $line = readline $gtfh) {
    next if $line =~ /#/;
  
    my ($chrom, $annot, $biotype, $start, $end, $a, $b, $c, $names) = split /\t/, $line;
    
    if ($biotype eq 'gene') {
      next unless $contigs->{$chrom};

      my $gene_id = '';
      if ($names =~ /gene_id "(.+?)"/) {
        $gene_id = $1;
      }
      my $window = "$chrom:$start-$end";
      push @genes, [$window, $gene_id];
    }
  }

  return @genes;
}

###############################################################################
sub get_contigs {
  my ($vcf) = @_;
  
  my $lines = `bcftools index -s $vcf`;
  
  my %contigs;

  for my $line (split /\n/, $lines) {
    my ($chrom, $size, $nsnps) = split "\t", $line;

    $contigs{$chrom} = $size;
  }

  return \%contigs;
}

sub create_windows {
  my ($contigs, $window_size) = @_;
  
  my @windows;
  for my $c (keys %$contigs) {
    my $size = $contigs->{$c};
    
    for(my $i = 0; $i < $size; $i += $window_size) {
      my $from = $i;
      my $to = $i + $window_size;
      my $window = "$c:$from-$to";
      push @windows, [$window];
    }
  }

  return @windows;
}

sub  create_vcf_subfiles {
  my ($vcf, $windows, $output_dir) = @_;

  mkdir $output_dir if not -e $output_dir;
  for my $couple (@$windows) {
    my ($window, $gene_id) = @$couple;
    my $vcf_name = '';
    if ($gene_id) {
      $vcf_name = $gene_id;
    } else {
      $vcf_name = $window;
      $vcf_name =~ s/:/_/;
    }
    $vcf_name .= '.vcf';
    my $out_path = "$output_dir/$vcf_name";
    `bcftools view $vcf -r $window > $out_path`;
  }
}

###############################################################################
# Parameters and usage
sub usage {
  my $error = shift;
  my $help = '';
  if ($error) {
    $help = "[ $error ]\n";
  }
  $help .= <<'EOF';
    
    --vcf <path>

    USE WINDOW
    --window <int>    : window size in kb
    --output_dir <int>: sub vcfs output dir

    OR USE GTF
    --gtf <path>
    
    --help            : show this help message
    --verbose         : show detailed progress
    --debug           : show even more information (for debugging purposes)
EOF
  print STDERR "$help\n";
  exit(1);
}

sub opt_check {
  my %opt = ();
  GetOptions(\%opt,
    "vcf:s",
    "window:s",
    "output_dir:s",
    "gtf_file:s",
    "help",
    "verbose",
    "debug",
  ) or usage();

  usage("Need vcf")                if not $opt{vcf};
  usage("Need window or gtf_file")        if not ($opt{window} xor $opt{gtf_file});
  usage("Need output_dir")                if not $opt{output_dir};
  usage()                if $opt{help};
  Log::Log4perl->easy_init($INFO) if $opt{verbose};
  Log::Log4perl->easy_init($DEBUG) if $opt{debug};
  return \%opt;
}
__END__
