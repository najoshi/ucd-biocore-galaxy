#!/usr/bin/perl

($sorp, $infile, $outfile) = @ARGV;

$cnt = 1;

open ($infd, "<$infile");
open ($outfd, ">$outfile");
while ($h1=<$infd>) {
	$seq=<$infd>;
	$h2=<$infd>;
	$q=<$infd>;	

	print $outfd "\@$cnt\#0/1\n$seq+\n$q";

	if ($sorp eq "paired") {
		$g1=<$infd>;
		$seq2=<$infd>;
		$g2=<$infd>;
		$q2=<$infd>;

		print $outfd "\@$cnt\#0/2\n$seq2+\n$q2";
	}

	$cnt++;
}
close ($outfd);
close ($infd);
