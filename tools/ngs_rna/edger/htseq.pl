#!/usr/bin/perl

use strict;
use warnings;
use Getopt::Std;
use File::Basename;
$| = 1;

# Grab and set all options
my %OPTIONS = (a => 0, i => "gene_id", m => "intersection-nonempty", s => "no", t => "exon");
getopts('a:cg:i:m:o:r:s:t:', \%OPTIONS);

die qq(
Usage:   HTSeq.pl [OPTIONS] Group1=sample1=<SAM/BAM file> [Group1=sample2=<SAM/BAM file> ... Group2=sampleN=<SAM/BAM file> ...]

OPTIONS:	-a	STR	skip all reads with alignment quality lower than the given minimum value (default: $OPTIONS{a})
			-c		reduce the matrix by removing any feature with no counts
			-g	STR	the features file in the GFF/GTF format
			-i	STR	GFF attribute to be used as feature ID (default: $OPTIONS{i})
			-m	STR	mode to handle reads overlapping more than one feature. Possible values for <mode> are union, intersection-strict and intersection-nonempty (default: $OPTIONS{m})
			-o	STR	output file name for expression matrix
			-r	STR	the name of the output report
			-s	STR	whether the data is from a strand-specific assay (default: $OPTIONS{s})
			-t	STR	feature type (3rd column in GFF file) to be used, all features of other type are ignored (default, suitable for RNA-Seq and Ensembl GTF files: $OPTIONS{t})

) if(@ARGV == 0);

my $sam_out;
my @counts;
my @features;
my %report;
my @samplenames;
my $current_group;
my @groups;
my @files;
my $groupcount = 0;
my %grouphash;

foreach my $input (@ARGV) {
	my ($group, $sample, $input) = split "::", $input;
	if(! defined $grouphash{$group}) {
		$groupcount++;
		$grouphash{$group} = "G${groupcount}:$group";
	}
	push @groups, $grouphash{$group};
	push @samplenames, $sample;
	push @files, $input;
}

for(my $index = 0; $index <= $#files; $index++) {
	my $input_file = $files[$index];
	my $sample = $samplenames[$index];
	
	# run htseq
	my @htseq;
	my $COMM;
	my $file_type = `file $input_file`;
	if(grep /text$/, $file_type ) {
		$COMM = "htseq-count -q -r pos -m $OPTIONS{m} -s $OPTIONS{s} -a $OPTIONS{a} -t $OPTIONS{t} -i $OPTIONS{i} $input_file $OPTIONS{g} | grep -v ^__";
		@htseq = `$COMM`;
	} else {
		$COMM = "samtools view $input_file | htseq-count -q -r pos -m $OPTIONS{m} -s $OPTIONS{s} -a $OPTIONS{a} -t $OPTIONS{t} -i $OPTIONS{i} - $OPTIONS{g} | grep -v ^__";
		@htseq = `$COMM`;
	}
	
	my $row = 0;
	$report{$sample} = "Command Used: $COMM\n";

	for(my $row = 0; $row <= $#htseq; $row++) {
		# store the report is an hash
		if(grep /^no_feature|^ambiguous|^too_low_aQual|^not_aligned|^alignment_not_unique/, $htseq[$row]) {
			$report{$sample} .= $htseq[$row];
		} else {
			# store the counts in a matrix
			chomp $htseq[$row];
			my ($feature, $value) = split "\t", $htseq[$row];
			$counts[$row][$index] = $value;
			if($input_file eq $files[0]) {
				push @features, $feature;
			}
		}
	}
}

# print the matrix
open(MATRIX, ">$OPTIONS{o}") || die "Could Not Create Output File $OPTIONS{o}!\n";
print MATRIX "#\t".join("\t", @groups)."\n";
print MATRIX "#Feature\t".join("\t", @samplenames)."\n";
for(my $row = 0; $row <= $#features; $row++) {
	if(defined $OPTIONS{c}) {
		my $rowsum = 0;
		$rowsum += $_ foreach @{ $counts[$row] };
		if(!$rowsum) {
			next;
		}
	}
	print MATRIX "$features[$row]\t".join("\t", @{ $counts[$row] })."\n";
}
close(MATRIX);

# print the report
open(REPORT, ">$OPTIONS{r}") || die "Could Not Create Output File $OPTIONS{r}!\n";
print REPORT "$groups[$_]:$samplenames[$_]\n$report{$samplenames[$_]}\n" foreach (0..$#samplenames);
close(REPORT);



