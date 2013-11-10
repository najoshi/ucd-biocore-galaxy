#!/usr/bin/perl

($bamfile, $bamindex, $outfile1, $outfileid, $newfilepath, $tempfilepath) = @ARGV;

system ("mkdir -p $tempfilepath");
chdir ("$tempfilepath");

system ("ln -s $bamindex file.bam.bai");
system ("ln -s $bamfile file.bam");
$headers = `/share/apps/samtools/samtools view -H file.bam`;

if ($headers eq "") {
	print STDERR "Error: No headers detected in BAM file. Aborting.\n";
	exit(1);
}

foreach $hd (split (/\n/, $headers)) {
	if (($chr)=$hd=~/^\@SQ\tSN:(.+?)\t/) {
		push (@chrs, $chr);
	}
}

system ("mkdir -p $newfilepath");
open ($tempscript, ">$tempfilepath/tempscript.sh");
foreach $chr (@chrs) {
	print $tempscript "/share/apps/samtools/samtools view -b -o $newfilepath/primary_$outfileid"."_$chr"."_visible_bam file.bam $chr &\n";
}

print $tempscript "wait\n";
close ($tempscript);

system ("sh tempscript.sh");
system ("rm -rf $tempfilepath");
