#!/usr/bin/perl

# deseq4_wrapper.pl $counts $column_types $comparison $top_table $diagnostic_html "$diagnostic_html.files_path" "$counts.name"

($counts, $column_types, $comparison, $top_table, $diagnostic_html, $diagnostic_html_files_path, $counts_name) = @ARGV;

# print "$counts, $column_types, $comparison, $top_table, $diagnostic_html, $diagnostic_html_files_path, $counts_name\n";

$output = `Rscript deseq4.R $counts $column_types $comparison $top_table $diagnostic_html "$diagnostic_html_files_path" "$counts_name" 2>&1`;
$exitcode = $? >> 8;

if ($exitcode != 0 && $output =~ /Error: /) {
	die "OUTPUT: $output\nExit code: $exitcode\n";
}

#exit(0);
