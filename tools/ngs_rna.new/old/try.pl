#!/usr/bin/perl

$output = `Rscript /opt/Bio/galaxy-dist/tools/ngs_rna/deseq3.R /opt/Bio/galaxy-dist/database/files/001/dataset_1789.dat unwounded,unwounded,unwounded,regenerated,regenerated,regenerated unwounded,regenerated /opt/Bio/galaxy-dist/database/files/001/dataset_1915.dat /opt/Bio/galaxy-dist/database/files/001/dataset_1916.dat "/opt/Bio/galaxy-dist/database/job_working_directory/1373/dataset_1916_files" "AC05TPACXX.transcripts.no_rRNA.counts" 2>&1`;

my $exitcode = $? >> 8;

print "OUTPUT: $output\nEXIT CODE: $exitcode\n";
