# sam2counts and DESeq in Galaxy

## About

This is a Galaxy package that wraps sam2counts and DESeq for RNA-Seq analysis using a transcriptome reference.  sam2counts takes SAM files that are created from an alignment to a transcriptome and creates counts of aligned reads for each transcript.  DESeq uses the DESeq package from Bioconductor in R and analyzes the count data from sam2counts. DESeq outputs a toptable of transcripts sorted by adjusted p-value and a page of diagnostic plots.

## Requirements

Python 2.6.5
pysam 0.6 (package for Python)
R 2.15
Bioconductor 2.10 (package for R)
DESeq 1.8.3 (package for R)
aroma.light 1.24.0 (package for R)
lattice 0.20-6 (package for R)

## Installation

stderr_wrapper.py and sam2counts_galaxy.py must be in the path or they can remain in the tools directory with the xml files.  deseq.R must be copied to the "tool-data" directory under the main Galaxy install directory.

## Use

sam2counts needs a SAM file (produced by aligning to a transcriptome) with header information as the input.  The count data produced from this SAM file gets fed into DESeq.
