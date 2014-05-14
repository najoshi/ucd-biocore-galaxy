#!/usr/bin/python

"""
sam2count_galaxy_edger.py -- Take SAM files and output a table of counts with column
names that are the filenames, and rowname that are the reference
names.
"""

VERSION = 0.90

import sys
import csv
from os import path
try:
    import pysam
except ImportError:
    sys.exit("pysam not installed; please install it\n")

import argparse

def SAM_file_to_counts(filename, sname, ftype="sam", extra=False, use_all_references=True):
    """
    Take SAM filename, and create a hash of mapped and unmapped reads;
    keys are reference sequences, values are the counts of occurences.

    Also, a hash of qualities (either 0 or >0) of mapped reads
    is output, which is handy for diagnostics.
    """
    counts = dict()
    unique = dict()
    nonunique = dict()
    mode = 'r'
    if ftype == "bam":
        mode = 'rb'
    sf = pysam.Samfile(filename, mode)

    if use_all_references:
        # Make dictionary of all entries in header
        try:
            for sn in sf.header['SQ']:
                if extra:
                    unique[sn['SN']] = 0
                    nonunique[sn['SN']] = 0
                counts[sn['SN']] = 0
        except KeyError:
            print "Sample file of sample " + sname + " does not have header."

    for read in sf:
        if not read.is_unmapped:
            id_name = sf.getrname(read.rname) if read.rname != -1 else 0

            if not use_all_references and not counts.get(id_name, False):
                ## Only make keys based on aligning reads, make empty hash
                if extra:
                    unique[id_name] = 0
                    nonunique[id_name] = 0
                ## initiate entry; even if not mapped, record 0 count
                counts[id_name] = counts.get(id_name, 0)
        
            
            counts[id_name] = counts.get(id_name, 0) + 1

            if extra:
                if read.mapq == 0:
                    nonunique[id_name] = nonunique[id_name] + 1
                else:
                    unique[id_name] = unique[id_name] + 1

    if extra:
        return {'counts':counts, 'unique':unique, 'nonunique':nonunique}

    return {'counts':counts}

def collapsed_nested_count_dict(counts_dict, all_ids, order=None):
    """
    Takes a nested dictionary `counts_dict` and `all_ids`, which is
    built with the `table_dict`. All files (first keys) in
    `counts_dict` are made into columns with order specified by
    `order`.

    Output is a dictionary with keys that are the id's (genes or
    transcripts), with values that are ordered counts. A header will
    be created on the first row from the ordered columns (extracted
    from filenames).
    """
    if order is None:
        col_order = counts_dict.keys()
    else:
        col_order = order

    collapsed_dict = dict()
    for i, filename in enumerate(col_order):
        for id_name in all_ids:
            if not collapsed_dict.get(id_name, False):
                collapsed_dict[id_name] = list()

            # get counts and append
            c = counts_dict[filename].get(id_name, 0)
            collapsed_dict[id_name].append(c)
    return {'table':collapsed_dict, 'header':col_order}


def counts_to_file(table_dict, outfilename, delimiter=','):
    """
    A function for its side-effect of writing `table_dict` (which
    contains a table and header), to `outfilename` with the specified
    `delimiter`.
    """
    writer = csv.writer(open(outfilename, 'a'), delimiter=delimiter, lineterminator='\n')
    table = table_dict['table']
    header = table_dict['header']
    
    #header_row = True
    for id_name, fields in table.items():
        #if header_row:
            #row = ['id'] + header
            #writer.writerow(row)
            #header_row = False

        if id_name == 0:
            continue
        row = [id_name]
        row.extend(fields)
        writer.writerow(row)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='parser for sam2counts')
    parser.add_argument("-d", "--delimiter", help="the delimiter (default: tab)", default='\t')
    parser.add_argument("-o", "--out-file", help="output filename (default: counts.txt)", default='counts.txt')
    parser.add_argument("-u", "--extra-output", help="output extra information on non-unique and unique mappers (default: False)")
    parser.add_argument("-r", "--use-all-references", dest="use_all_references",
                      help="Use all the references from the SAM header (default: True)",
                      default=True, action="store_false")
    parser.add_argument("-f", "--extra-out-files", dest="extra_out_files",
                      help="comma-delimited filenames of unique and non-unique output "
                      "(default: unique.txt,nonunique.txt)",
                      default='unique.txt,nonunique.txt')
    parser.add_argument("-v", "--verbose", dest="verbose",
                      help="enable verbose output")
    parser.add_argument("--bam-file", help="bam file", nargs="+", action="append", required=True)
    parser.add_argument("--group", help="group", nargs="+", action="append", required=True)
    parser.add_argument("--treatment", help="treatment", nargs="+", action="append", required=True)
    parser.add_argument("--sample-name", help="sample name", nargs="+", action="append", required=True)
    parser.add_argument("--file-type", help="file type", nargs="+", action="append", required=True, choices=['sam','bam'])
    args = parser.parse_args()

    args.bam_file = [item for sublist in args.bam_file for item in sublist]
    args.group = [item for sublist in args.group for item in sublist]
    args.treatment = [item for sublist in args.treatment for item in sublist]
    args.sample_name = [item for sublist in args.sample_name for item in sublist]
    args.file_type = [item for sublist in args.file_type for item in sublist]
    #print(args.sample_name)

    if (len(args.sample_name) != len(set(args.sample_name))):
        parser.error("Sample names must be unique.")

    if not(len(args.bam_file) == len(args.group) and len(args.group) == len(args.treatment) and len(args.treatment) == len(args.sample_name) and len(args.sample_name) == len(args.file_type)):
        parser.error("Number of total BAM files, groups, treatments, sample names, and types must be the same.")

    file_counts = dict()
    file_unique_counts = dict()
    file_nonunique_counts = dict()
    all_ids = list()

    ## do a pre-run check that all files exist
    for full_filename in args.bam_file:
        if not path.exists(full_filename):
            parser.error("file '%s' does not exist" % full_filename)

    outf = open(args.out_file, "w")
    outf.write("#")
    for (g,t) in zip(args.group,args.treatment):
        outf.write("\t" + g + ":" + t) 
    outf.write("\n#Feature")
    for s in args.sample_name:
        outf.write("\t" + s)
    outf.write("\n")
    outf.close()

    for (full_filename, sn, ft) in zip(args.bam_file, args.sample_name, args.file_type):
        ## read in SAM file, extract counts, and unpack counts
        tmp = SAM_file_to_counts(full_filename, sn, ftype=ft, extra=args.extra_output,
                                 use_all_references=args.use_all_references)

        if args.extra_output:
            counts, unique, nonunique = tmp['counts'], tmp['unique'], tmp['nonunique']
        else:
            counts = tmp['counts']

        ## save counts, and unique/non-unique counts
        file_counts[sn] = counts

        if args.extra_output:
            file_unique_counts[sn] = unique
            file_nonunique_counts[sn] = nonunique

        ## add all ids encountered in this in this file
        all_ids.extend(file_counts[sn].keys())

    ## Uniquify all_ids, and then take the nested file_counts
    ## dictionary, collapse, and write to file.
    all_ids = set(all_ids)
    table_dict = collapsed_nested_count_dict(file_counts, all_ids, order=args.sample_name)
    counts_to_file(table_dict, args.out_file, delimiter=args.delimiter)

    if args.extra_output:
        unique_fn, nonunique_fn = args.extra_out_files.split(',')
        unique_table_dict = collapsed_nested_count_dict(file_unique_counts, all_ids, order=files)
        nonunique_table_dict = collapsed_nested_count_dict(file_nonunique_counts, all_ids, order=files)
        
        counts_to_file(unique_table_dict, unique_fn, delimiter=args.delimiter)
        counts_to_file(nonunique_table_dict, nonunique_fn, delimiter=args.delimiter)
        
