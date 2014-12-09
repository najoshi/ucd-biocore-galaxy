#!/usr/bin/perl

$version = $ARGV[0];


foreach $xmlfile (@ARGV[1 .. $#ARGV]) {

open ($xf, "<$xmlfile");
open ($tmpfile, ">/tmp/tmpfile");

while (<$xf>) {
	if ($_ =~ /^<tool/) {
		if ($_ =~ /^<tool.+version=/) {
			$_ =~ s/version=".+?"/version="$version"/;
		}

		else {
			$_ =~ s/>/ version="$version">/;
		}
	}

	print $tmpfile $_;
}
close($xf);
close($tmpfile);


system ("mv /tmp/tmpfile $xmlfile");

}
