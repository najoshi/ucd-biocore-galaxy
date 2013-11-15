#!/usr/bin/perl

($bamfile, $bamindex, $outfile1, $outfileid, $newfilepath, $tempfilepath) = @ARGV;

system ("mkdir -p $tempfilepath");
chdir ("$tempfilepath");

system ("ln -s $bamindex file.bam.bai");
system ("ln -s $bamfile file.bam");
$headers = `/share/apps/samtools/samtools view -H file.bam`;

if ($headers eq "") {
	print STDERR "Error: No headers detected in BAM file. Aborting.\n";
    system ("rm -rf $tempfilepath");
	exit(1);
}

foreach $hd (split (/\n/, $headers)) {
	if (($chr)=$hd=~/^\@SQ\tSN:(.+?)\t/) {
		push (@chrs, $chr);
	}
}

open ($outf, ">$outfile1");
print $outf "# Log File\n\nHeaders found:\n$headers\nSplitting BAM file into multiple files.\nNote: You may need to refresh history multiple times to get all the output files.\n";
close ($outf);

system ("mkdir -p $newfilepath");
open ($tempscript, ">$tempfilepath/tempscript.sh");
foreach $chr (@chrs) {
	$chrsub = $chr;
	$chrsub =~ s/_/-/g;
	print $tempscript "/share/apps/samtools/samtools view -b -o $newfilepath/primary_$outfileid"."_$chrsub"."_visible_bam file.bam $chr &\n";
}

print $tempscript "wait\n";
close ($tempscript);

system ("sh tempscript.sh");
system ("rm -rf $tempfilepath");
