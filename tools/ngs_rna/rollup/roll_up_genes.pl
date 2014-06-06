#!/usr/bin/perl
use strict;
use warnings;
#use POSIX;
use Getopt::Std;


$| = 1;

# Monica Britton August 23, 2012
# Revised June 5, 2014 to include multi-line headers in comment format
#
#	roll_up_genes.pl
#	
#   Rolls up counts by transcript ID into counts by gene ID (does not combine across columns)
#
# Input file should have transcipt ID in first column, and gene IDs in subsequent columns (will us the last one)
#
#TCONS_00001181  XLOC_000786     AAEL014214
#TCONS_00001182  XLOC_000787     
#  
# 
#
# The Bioinformatics Core at UC Davis Genome Center
# http://bioinformatics.ucdavis.edu
# Copyright (c) 2012 The Regents of University of California, Davis Campus.
# All rights reserved.


my $usage = "\nusage: $0 -c <counts table> -g gene ID correspondence table\n".
          "\n           This script will roll up transcript counts into gene counts\n".
			"           Input is a table of counts by transcript and a \n".
			"           transcript to gene ID correspondence table (tab delimited, transcript ID in first column)\n".
			"           can have multiple gene ID columns (last is used).\n".
			"           If there is no gene ID, then transcript ID will be used as gene ID\n".
			"           Any line beginning with # is comment, kept in output\n".
			"           Output is a new table of counts by gene ID\n\n".
			"   -c      counts table\n".
			"   -g      gene ID correspondence table\n\n";

			
our($opt_c, $opt_g );

getopts('c:g:') or die $usage;
if ( not -e $opt_c) {print STDERR "counts table not found\n"; die $usage;}
if ( not -e $opt_g) {print STDERR "gene ID correspondence table not found\n"; die $usage;}

my %gene_table;
my %gene_counts;

open (GENE_FILE, "<", $opt_g );											# open vcf file for reading

while (<GENE_FILE>) {
	chomp;
	
	my @ly = split ( /\t/, $_ );

	my $ly_elements = @ly;
	
	my $gene_loc = $ly_elements - 1 ;
	
	$gene_table{$ly[0]}{"GENE_ID"} = $ly[$gene_loc];					# store correspondence between transcript ID and gene ID
	
#	$gene_counts{$ly[$gene_loc]} = 0 ;									# initialize gene count to zero.
	
#	print "@ly\t" . $gene_table{$ly[0]}{"GENE_ID"} . "\t" . $gene_counts{$ly[$gene_loc]} . "\n";

}

close (GENE_FILE);

my $head_count = 0;
my $hy_samples = 1;

open (COUNT_FILE, "<", $opt_c );										# open count file for reading

while (<COUNT_FILE>) {
	chomp;
	
#	print "$_\n";
	
	if ( $_ =~ m/^#/ ) { print "$_\n"; next; }
	
	if ( $head_count < 1 ) { 
		
#		print "$_\n"; 
		
		$head_count = $head_count + 2;
		
		my @hy = split ( /\t/, $_ );
		
		$hy_samples = @hy - 1;
		
#		print "There are $hy_samples samples\n";
		
		}
		
#	else {	
		
		my @ty = split ( /\t/, $_ );
		
		if ( defined $gene_table{$ty[0]}{"GENE_ID"} ) {					# check that the Gene_ID is in $gene_table
			
			my $gene_id = $gene_table{$ty[0]}{"GENE_ID"};
			
			for (my $i=1; $i<=$hy_samples; $i++) {						# iterating through sample columns
			
				if ( defined $gene_counts{$gene_id}{$i}{"SUM"} ) {		# add count to sum
					
					$gene_counts{$gene_id}{$i}{"SUM"} = $gene_counts{$gene_id}{$i}{"SUM"} + $ty[$i];
					
				}
				
				else { $gene_counts{$gene_id}{$i}{"SUM"} =  $ty[$i]; }
			
			}
			
		}
		
		else { print STDERR "$ty[0] not found in Gene Table!"; die; }
		
#	}
	
}

close (COUNT_FILE);


foreach my $key1 (sort keys %gene_counts) {

	print $key1;
	
	for (my $i=1; $i<=$hy_samples; $i++) {					# iterating through samples

		print "\t" . $gene_counts{$key1}{$i}{"SUM"};

	}
	
	print "\n";

}


exit;
