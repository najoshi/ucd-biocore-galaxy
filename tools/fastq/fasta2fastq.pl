#!/usr/bin/perl

open ($infile, "<$ARGV[0]");
open ($outfile, ">$ARGV[1]");
$qual_value = $ARGV[2];

while ($header=<$infile>) {
    chomp $header;
    $seq=<$infile>;
    chomp $seq;

    ($id)=$header=~/^>(.+)/;

    print $outfile "\@$id\n$seq\n+\n";
    print $outfile $qual_value x length($seq);
    print $outfile "\n";
}

close ($infile);
close ($outfile);
