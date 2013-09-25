cargs <- commandArgs()
cargs <- cargs[(which(cargs == "--args")+1):length(cargs)]

filename = cargs[1]
datatype = cargs[2]
reportfile = cargs[3]
temp.files.path = cargs[4]

sink("/dev/null")

dir.create(temp.files.path, recursive=TRUE)

qualtype = if (datatype=="fastqsanger") "sanger" else if (datatype=="fastqillumina") "illumina" else if (datatype=="fastqsolexa") "solexa" else "sanger"

filetype = if (datatype=="fastqsanger") "fastq" else if (datatype=="fastqillumina") "fastq" else if (datatype=="fastqsolexa") "fastq" else datatype

library(qrqc)

setwd(temp.files.path)
seqdata = readSeqFile(filename, type=filetype, quality=qualtype)
report( seqdata, outdir=".", type="html", hash=TRUE, kmers=TRUE )

# move data around so that Galaxy can see it properly
file.rename ( file.path( temp.files.path, "qrqc-report.html"), reportfile)

sink(NULL)
