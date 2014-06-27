#!/usr/bin/perl

open ($infile, "<$ARGV[0]");
open ($outfile, ">$ARGV[1]");
$qual_value = $ARGV[2];

$seq="";
$first=0;
while ($line=<$infile>) {
    chomp $line;

    if ($line =~ /^>(.+?)\s/ || $line =~ /^>(.+?)$/) {

	if ($first == 0) {$first=1;}
	else {
		print $outfile "\@$id\n$seq\n+\n";
		print $outfile $qual_value x length($seq);
		print $outfile "\n";
		$seq = "";
	}

	$id = $1;
    }

    else {
	$seq .= $line;
    }
}

print $outfile "\@$id\n$seq\n+\n";
print $outfile $qual_value x length($seq);
print $outfile "\n";

close ($infile);
close ($outfile);
