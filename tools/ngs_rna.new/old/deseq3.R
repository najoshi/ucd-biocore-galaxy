## Run DESeq (oversimplified) in Galaxy.
## Format: Rscript deseq3.R tab-delimited-counts.txt comma,delimited,col,names,except,first,gene,id,col comparison,here outputfile
##
## The incoming data must have the first column be the gene names, and
## the rest raw counts. The next argument is the list of groups
## (i.e. liver,liver,kidney,kidney), and the final argument is the
## comparison, i.e. kidney,liver. This produces a top-table called
## top-table.txt ordered by p-value.

cargs <- commandArgs()
cargs <- cargs[(which(cargs == "--args")+1):length(cargs)]

countstable <- cargs[1]
conds <- unlist(strsplit(cargs[2], ','))
comparison <- unlist(strsplit(cargs[3], ','))
outputfile <- cargs[4]
diag.html = cargs[5]
temp.files.path = cargs[6]
counts.name = cargs[7]

# the comparison must only have two values and the conds must
# be a vector from those values, at least one of each.
if (length(unique(conds)) != 2) {
	#write("Error: You can only have two condition types.", stderr())
	stop("You can only have two columns types: ", cargs[2])
}

if (length(comparison) != 2) {
	#write("Comparison types must be a tuple.", stderr())
	stop("Comparison type must be a tuple: ", cargs[3])
}

if (!identical(sort(comparison), sort(unique(conds)))) {
	#write("Column types must use the two names from Comparison type, and vice versa.  Must have at least one of each in the Column types.", stderr())
	stop("Column types must use the two names from Comparison type, and vice versa.  Must have at least one of each in the Column types.\nColumn types: ", cargs[2], "\n", "Comparison type: ", cargs[3])
}


sink("/dev/null")
dir.create(temp.files.path, recursive=TRUE)
library(DESeq)

d <- read.table(countstable, sep="\t", header=TRUE, row.names=1)
cds <- newCountDataSet(d, conds)
cds <- estimateSizeFactors(cds)
cds <- estimateDispersions( cds )

plotDispEsts <- function( cds )
{
   plot(
       rowMeans( counts( cds, normalized=TRUE ) ),
           fitInfo(cds)$perGeneDispEsts,
           pch = '.', log="xy" )
           xg <- 10^seq( -.5, 5, length.out=300 )
           lines( xg, fitInfo(cds)$dispFun( xg ), col="red" )
        }

temp.disp.est.plot = file.path( temp.files.path, "DispersionEstimatePlot.png" )
png( temp.disp.est.plot, width=500, height=500 )
plotDispEsts( cds )
dev.off()

res <- nbinomTest(cds, comparison[1], comparison[2])

write.table(res[order(res$padj), ], file=outputfile, quote=FALSE, row.names=TRUE, sep="\t")


plotDE <- function( res )
   plot(
      res$baseMean,
      res$log2FoldChange,
      log="x", pch=20, cex=.3,
      col = ifelse( res$padj < .1, "red", "black" ) )

temp.de.plot = file.path( temp.files.path, "DiffExpPlot.png")
png( temp.de.plot, width=500, height=500 )
plotDE( res )
dev.off()


temp.pval.plot = file.path( temp.files.path, "PvalHist.png")
png( temp.pval.plot, width=500, height=500 )
hist(res$pval, breaks=100, col="skyblue", border="slateblue", main="")
dev.off()

cdsBlind <- estimateDispersions( cds, method="blind" )
vsd <- getVarianceStabilizedData( cdsBlind )

temp.heatmap.plot = file.path( temp.files.path, "heatmap.png" )
png( temp.heatmap.plot, width=500, height=500 )
select <- order(res$pval)[1:100]
colors <- colorRampPalette(c("white","darkblue"))(100)
heatmap( vsd[select,], col = colors, scale = "none", Rowv=NULL, main="")
dev.off()


mod_lfc <- (rowMeans( vsd[, conditions(cds)==comparison[2], drop=FALSE] ) - rowMeans( vsd[, conditions(cds)==comparison[1], drop=FALSE] ))
lfc <- res$log2FoldChange
finite <- is.finite(lfc)
table(as.character(lfc[!finite]), useNA="always")
largeNumber <- 10
lfc <- ifelse(finite, lfc, sign(lfc) * largeNumber)

temp.modlfc.plot = file.path( temp.files.path, "modlfc.png" )
png( temp.modlfc.plot, width=500, height=500 )
plot( lfc, mod_lfc, pch=20, cex=.3, col = ifelse( finite, "#80808040", "red" ) )
abline( a=0, b=1, col="#40404040" )
dev.off()

temp.sampclust.plot = file.path( temp.files.path, "sampclust.png" )
png( temp.sampclust.plot, width=500, height=500 )
dists <- dist( t( vsd ) )
heatmap( as.matrix( dists ),
    symm=TRUE, scale="none", margins=c(10,10),
        col = colorRampPalette(c("darkblue","white"))(100),
        labRow = paste( pData(cdsBlind)$condition, pData(cdsBlind)$type ) )
dev.off()


library(aroma.light)
library(lattice)

mdsPlot <- function(x, conds=NULL, cex=1, ...) {
  d <- dist(x)

  mds.fit <- cmdscale(d, eig=TRUE, k=2)

  mds.d <- data.frame(x1=mds.fit$points[, 1],
                      x2=mds.fit$points[, 2],
                      labels=rownames(x))
  if (!is.null(conds))
    mds.d$treatment <- as.factor(conds)

  if (!is.null(conds)) {
    p <- xyplot(x2 ~ x1, group=treatment, data=mds.d, panel=function(x, y, ..., groups, subscripts) {
      panel.text(x, y, mds.d$labels[subscripts], cex=cex, col=trellis.par.get()$superpose.line$col[groups])
    }, ...)
  } else {
    p <- xyplot(x2 ~ x1, data=mds.d, panel=function(x, y, ..., groups, subscripts) {
      panel.text(x, y, mds.d$labels[subscripts], cex=cex)
    }, ...)
  }
  return(p)
}

dn = normalizeQuantileRank(as.matrix(d))
p = mdsPlot(t(log10(dn+1)), conds=conds, xlab="dimension 1", ylab="dimension 2", main="")

temp.mds.plot = file.path( temp.files.path, "mds.png" )
png( temp.mds.plot, width=500, height=500 )
print(p)
dev.off()


file.conn = file(diag.html, open="w")
writeLines( c("<html><body bgcolor='lightgray'>"), file.conn)
writeLines( c("<h2><u>Diagnostic Plots for ", counts.name, "</u></h2>"), file.conn)
writeLines( c("<h3>Dispersion Estimates</h3>"), file.conn)
writeLines( c("<img src='DispersionEstimatePlot.png'><br/><br/>"), file.conn)
writeLines( c("<h3>Differential Expression - Base Mean vs. Log2 Fold Change</h3>"), file.conn)
writeLines( c("<img src='DiffExpPlot.png'><br/><br/>"), file.conn)
writeLines( c("<h3>P-value histogram</h3>"), file.conn)
writeLines( c("<img src='PvalHist.png'><br/><br/>"), file.conn)
writeLines( c("<h3>Top 100 Genes/Transcripts by P-value</h3>"), file.conn)
writeLines( c("<img src='heatmap.png'><br/><br/>"), file.conn)
writeLines( c("<h3>Moderated LFC vs. LFC</h3>"), file.conn)
writeLines( c("<img src='modlfc.png'><br/><br/>"), file.conn)
writeLines( c("<h3>VST Sample Clustering</h3>"), file.conn)
writeLines( c("<img src='sampclust.png'><br/><br/>"), file.conn)
writeLines( c("<h3>MDS Plot</h3>"), file.conn)
writeLines( c("<img src='mds.png'><br/><br/>"), file.conn)
writeLines( c("</body></html>"), file.conn)
close(file.conn)

sink(NULL)
