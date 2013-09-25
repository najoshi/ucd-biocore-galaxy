## Run DESeq (oversimplified) in Galaxy.
## Format: Rscript deseq.R --args tab-delimited-counts.txt comma,delimited,col,names,except,first,gene,id,col comparison,here outputfile
##
## The incoming data must have the first column be the gene names, and
## the rest raw counts. The next argument is the list of groups
## (i.e. liver,liver,kidney,kidney), and the final argument is the
## comparison, i.e. kidney,liver. This produces a top-table called
## top-table.txt ordered by p-value.

cargs <- commandArgs()
cargs <- cargs[(which(cargs == "--args")+1):length(cargs)]

countstable <- cargs[1]
print(countstable)
conds <- unlist(strsplit(cargs[2], ','))
print(conds)
comparison <- unlist(strsplit(cargs[3], ','))
print(comparison)
outputfile <- cargs[4]

library(DESeq)

d <- read.table(countstable, sep="\t", header=TRUE)
rownames(d) <- d[, 1]
d <- d[, -1]
cds <- newCountDataSet(d, conds)
cds <- estimateSizeFactors(cds)
cds <- estimateVarianceFunctions(cds)
res <- nbinomTest(cds, comparison[1], comparison[2])

write.table(res[order(res$padj), ], file=outputfile, quote=FALSE, row.names=FALSE, sep="\t")

