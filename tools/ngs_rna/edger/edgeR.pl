#!/bin/perl

#EdgeR.pl Version 0.0.3
#Contributors: Monica Britton, Blythe Durbin-Johnson, Joseph Fass, Nikhil Joshi, Alex Mawla

use strict;
use warnings;
use Getopt::Std;
use File::Basename;
use File::Path qw(make_path remove_tree);

$| = 1;

my %OPTIONS = (a => "glm",  d => "tag", f => "BH", r => 5, u => "movingave");

getopts('a:d:e:f:h:lmn:o:r:tu:', \%OPTIONS);

 
die qq(
Usage:   edgeR.pl [OPTIONS] factor::factor1::levels [factor::factor2::levels ...] cp::cont_pred1::values [cp::cont_pred2::values ...] cnt::contrast1 [cnt::contrast2] matrix

OPTIONS:	-a	STR	Type Of Analysis [glm, pw, limma] (default: $OPTIONS{a})
			-d	STR	The dispersion estimate to use for GLM analysis [tag] (default: $OPTIONS{d})
			-e	STR	Path to place additional output files
			-f	STR	False discovery rate adjustment method [BH] (default: $OPTIONS{f})
			-h	STR	Name of html file for additional files
			-l		Output the normalised digital gene expression matrix in log2 format (only applicable when using limma and -n is also specified)
			-m		Perform all pairwise comparisons
			-n	STR	File name to output the normalised digital gene expression matrix (only applicable when usinf glm or limma model)
			-o	STR	File name to output csv file with results
			-r	INT	Common Dispersion Rowsum Filter, ony applicable when 1 factor analysis selected (default: $OPTIONS{r})
			-t		Estimate Tagwise Disp when performing 1 factor analysis
			-u	STR	Method for allowing the prior distribution for the dispersion to be abundance- dependent ["movingave"] (default: $OPTIONS{u})
 
) if(!@ARGV);

my $matrix = pop @ARGV;

make_path($OPTIONS{e});
open(Rcmd,">$OPTIONS{e}/r_script.R") or die "Cannot open $OPTIONS{e}/r_script.R\n\n";
print Rcmd "
zz <- file(\"$OPTIONS{e}/r_script.err\", open=\"wt\")
sink(zz)
sink(zz, type=\"message\")

library(edgeR)
library(limma)

toc <- read.table(\"$matrix\", sep=\"\\t\", comment=\"\", as.is=T)
groups <- sapply(toc[1, -1], strsplit, \":\")
for(i in 1:length(groups)) { g <- make.names(groups[[i]][2]); names(groups)[i] <- g; groups[[i]] <- groups[[i]][-2] }
colnames(toc) <- make.names(toc[2,])
toc[,1] <- gsub(\",\", \".\", toc[,1])
tagnames <- toc[-(1:2), 1]
rownames(toc) <- toc[,1]
toc <- toc[-(1:2), -1]
for(i in colnames(toc)) toc[, i] <- as.numeric(toc[,i])
norm_factors <- calcNormFactors(as.matrix(toc))

pw_tests <- list()
uniq_groups <- unique(names(groups))
for(i in 1:(length(uniq_groups)-1)) for(j in (i+1):length(uniq_groups)) pw_tests[[length(pw_tests)+1]] <- c(uniq_groups[i], uniq_groups[j])
DGE <- DGEList(toc, lib.size=norm_factors*colSums(toc), group=names(groups))
pdf(\"$OPTIONS{e}/MA_plots_normalisation.pdf\", width=14)
for(i in 1:length(pw_tests)) {
	j <- c(which(names(groups) == pw_tests[[i]][1])[1], which(names(groups) == pw_tests[[i]][2])[1])
	par(mfrow = c(1, 2))
	maPlot(toc[, j[1]], toc[, j[2]], normalize = TRUE, pch = 19, cex = 0.2, ylim = c(-10, 10), main=paste(\"MA Plot\", colnames(toc)[j[1]], \"vs\", colnames(toc)[j[2]]))
	grid(col = \"blue\")
	abline(h = log2(norm_factors[j[2]]), col = \"red\", lwd = 4)
	maPlot(DGE\$counts[, j[1]]/DGE\$samples\$lib.size[j[1]], DGE\$counts[, j[2]]/DGE\$samples\$lib.size[j[2]], normalize = FALSE, pch = 19, cex = 0.2, ylim = c(-8, 8), main=paste(\"MA Plot\", colnames(toc)[j[1]], \"vs\", colnames(toc)[j[2]], \"Normalised\"))
	grid(col = \"blue\")
}
dev.off()
pdf(file=\"$OPTIONS{e}/MDSplot.pdf\")
plotMDS(DGE, main=\"MDS Plot\", col=as.numeric(factor(names(groups)))+1, xlim=c(-3,3))
dev.off()
tested <- list()
";

my $all_cont;
my @add_cont;
my @fact;
my @fact_names;
my @cp;
my @cp_names;
if(@ARGV) {
	foreach my $input (@ARGV) {
		my @tmp = split "::", $input;
		if($tmp[0] eq "factor") {
			$tmp[1] =~ s/[ \?\(\)\[\]\/\\=+<>:;\"\',\*\^\|\&-]/./g;
			push @fact_names, $tmp[1];
			$tmp[2] =~ s/:/\", \"/g;
			$tmp[2] = "\"".$tmp[2]."\"";
			push @fact, $tmp[2];
		} elsif($tmp[0] eq "cp") {
			$tmp[1] =~ s/[ \?\(\)\[\]\/\\=+<>:;\"\',\*\^\|\&-]/./g;
			push @cp_names, $tmp[1];
			$tmp[2] =~ s/:/, /g;
			push @cp, $tmp[2];
		} elsif($tmp[0] eq "cnt") {
			push @add_cont, $tmp[1];
		} else {
			die("Unknown Input: $input\n");
		}
	}
}

if($OPTIONS{a} eq "pw") {
	print Rcmd "
disp <- estimateCommonDisp(DGE, rowsum.filter=$OPTIONS{r})
";
	if(defined $OPTIONS{t}) {
		print Rcmd "
disp <- estimateTrendedDisp (disp)
disp <- estimateTagwiseDisp(disp, trend=\"$OPTIONS{u}\")
pdf(file=\"$OPTIONS{e}/Tagwise_Dispersion_vs_Abundance.pdf\")
plotBCV(disp, cex=0.4)
abline(h=disp\$common.dispersion, col=\"firebrick\", lwd=3)
dev.off()
";
	}
	print Rcmd "
for(i in 1:length(pw_tests)) {
	tested[[i]] <- exactTest(disp, pair=pw_tests[[i]])
	names(tested)[i] <- paste(pw_tests[[i]][2], \"-\", pw_tests[[i]][1], sep=\"\")
}
pdf(file=\"$OPTIONS{e}/Smear_Plots.pdf\")
for(i in 1:length(pw_tests)) {
	dt <- decideTestsDGE(tested[[i]], p.value=0.05, adjust.method=\"$OPTIONS{f}\")
	if(sum(dt) > 0) {
		de_tags <- rownames(disp)[which(dt != 0)]
		ttl <- \"Diff. Exp. Genes With adj. Pvalue < 0.05\"
	} else {
		de_tags <- rownames(topTags(tested[[i]], n=100)\$table)
		ttl <- \"Top 100 tags\"
	}

	if(length(dt) < 5000) {
		pointcex = 0.5
	} else {
		pointcex = 0.2
	}
	plotSmear(disp, pair=pw_tests[[i]], de.tags = de_tags, main = paste(\"Smear Plot\", names(tested)[i]), cex=0.5)
	abline(h = c(-1, 1), col = \"blue\")
	legend(\"topright\", c(\"2 Fold Change\", ttl) , lty=c(1, NA), pch=c(NA, 19), pt.cex=0.5, col=c(\"blue\", \"red\"), bty=\"n\")
}
dev.off()
"; 
} 
elsif($OPTIONS{a} eq "glm") {
	for(my $fct = 0; $fct <= $#fact_names; $fct++) {
		print Rcmd "
			$fact_names[$fct] <- c($fact[$fct])
		";
	}
	for(my $fct = 0; $fct <= $#cp_names; $fct++) {
		print Rcmd "
			$cp_names[$fct] <- c($cp[$fct])
		";
	}
	my $all_fact = "";
	if(@fact_names) {
		foreach (@fact_names) {
			$all_fact .= " + factor($_)";
		}
    	}
	my $all_cp = "";
	if(@cp_names) {
		$all_cp = " + ".join(" + ", @cp_names);
	}
	print Rcmd "
		group_fact <- factor(names(groups))
		design <- model.matrix(~ -1 + group_fact${all_fact}${all_cp})
		colnames(design) <- sub(\"group_fact\", \"\", colnames(design))
	";
	foreach my $fct (@fact_names) {
		print Rcmd "
			colnames(design) <- make.names(sub(\"factor.$fct.\", \"\", colnames(design)))
		";
	}
	if($OPTIONS{d} eq "tag") {
		print Rcmd "
			disp <- estimateGLMCommonDisp(DGE, design)
			disp <- estimateGLMTrendedDisp(disp, design)
			disp <- estimateGLMTagwiseDisp(disp, design)
			fit <- glmFit(disp, design)
			pdf(file=\"$OPTIONS{e}/Tagwise_Dispersion_vs_Abundance.pdf\")
			plotBCV(disp, cex=0.4)
			dev.off()
		";
	}
	if(@add_cont) {
		$all_cont = "\"".join("\", \"", @add_cont)."\"";
		print Rcmd "
			cont <- c(${all_cont})
			for(i in uniq_groups)  cont <- gsub(paste(groups[[i]], \"([^0-9])\", sep=\"\"), paste(i, \"\\\\1\", sep=\"\"), cont)
			for(i in uniq_groups)  cont <- gsub(paste(groups[[i]], \"\$\", sep=\"\"), i, cont)
";
	} else {
		print Rcmd "
cont <- NULL
";
	}
	if(defined $OPTIONS{m}) {
		print Rcmd "
for(i in 1:length(pw_tests)) cont <- c(cont, paste(pw_tests[[i]][2], \"-\", pw_tests[[i]][1], sep=\"\"))
";
	}
	if(!defined $OPTIONS{m} && !@add_cont){
		die("No Contrasts have been specified, you must at least either select multiple pairwise comparisons or specify a custom contrast\n");
	}
	print Rcmd "
fit <- glmFit(disp, design)
cont <- makeContrasts(contrasts=cont, levels=design)
for(i in colnames(cont)) tested[[i]] <- glmLRT(fit, contrast=cont[,i])
pdf(file=\"$OPTIONS{e}/Smear_Plots.pdf\")
for(i in colnames(cont)) {
        dt <- decideTestsDGE(tested[[i]], p.value=0.05, adjust.method=\"$OPTIONS{f}\")
        if(sum(dt) > 0) {
                de_tags <- rownames(disp)[which(dt != 0)]
                ttl <- \"Diff. Exp. Genes With adj. Pvalue < 0.05\"
        } else {
                de_tags <- rownames(topTags(tested[[i]], n=100)\$table)
                ttl <- \"Top 100 tags\"
        }

        if(length(dt) < 5000) {
                pointcex = 0.5
        } else {
                pointcex = 0.2
        }
        plotSmear(disp, de.tags = de_tags, main = paste(\"Smear Plot\", i), cex=pointcex)
        abline(h = c(-1, 1), col = \"blue\")
        legend(\"topright\", c(\"2 Fold Change\", ttl) , lty=c(1, NA), pch=c(NA, 19), pt.cex=0.5, col=c(\"blue\", \"red\"), bty=\"n\")
}
dev.off()

	";
	if(defined $OPTIONS{n}) {
		print Rcmd "
			tab <- data.frame(ID=rownames(fit\$fitted.values), fit\$fitted.values, stringsAsFactors=F)
			write.table(tab, \"$OPTIONS{n}\", quote=F, sep=\"\\t\", row.names=F)
		";
	}
} elsif($OPTIONS{a} eq "limma") {
	for(my $fct = 0; $fct <= $#fact_names; $fct++) {
		print Rcmd "
$fact_names[$fct] <- c($fact[$fct])
";
	}
	for(my $fct = 0; $fct <= $#cp_names; $fct++) {
		print Rcmd "
$cp_names[$fct] <- c($cp[$fct])
";
	}
	my $all_fact = "";
	if(@fact_names) {
		foreach (@fact_names) {
			$all_fact .= " + factor($_)";
		}
	}
	my $all_cp = "";
	if(@cp_names) {
		$all_cp = " + ".join(" + ", @cp_names);
	}
	print Rcmd "

group_fact <- factor(names(groups))
design <- model.matrix(~ -1 + group_fact${all_fact}${all_cp})
colnames(design) <- sub(\"group_fact\", \"\", colnames(design))
";
	foreach my $fct (@fact_names) {
		print Rcmd "
colnames(design) <- make.names(sub(\"factor.$fct.\", \"\", colnames(design)))
";
	}
	print Rcmd "
isexpr <- rowSums(cpm(toc)>1) >= 1
toc <- toc[isexpr, ]
pdf(file=\"$OPTIONS{e}/LIMMA_voom.pdf\")
y <- voom(toc, design, plot=TRUE, lib.size=colSums(toc)*norm_factors)
dev.off()

pdf(file=\"$OPTIONS{e}/LIMMA_MDS_plot.pdf\")
plotMDS(y, labels=colnames(toc), col=as.numeric(factor(names(groups)))+1, gene.selection=\"common\")
dev.off()
fit <- lmFit(y, design)
";
	if(defined $OPTIONS{n}) {
		if(defined $OPTIONS{l}) {
			print Rcmd "
tab <- data.frame(ID=rownames(y\$E), y\$E, stringsAsFactors=F)
";
		} else {
			print Rcmd "
tab <- data.frame(ID=rownames(y\$E), 2^y\$E, stringsAsFactors=F)
";
		}
		print Rcmd "
write.table(tab, \"$OPTIONS{n}\", quote=F, sep=\"\\t\", row.names=F)
";
	}
	if(@add_cont) {
		$all_cont = "\"".join("\", \"", @add_cont)."\"";
		print Rcmd "
cont <- c(${all_cont})
for(i in uniq_groups)  cont <- gsub(paste(groups[[i]], \"([^0-9])\", sep=\"\"), paste(i, \"\\\\1\", sep=\"\"), cont)
for(i in uniq_groups)  cont <- gsub(paste(groups[[i]], \"\$\", sep=\"\"), i, cont)
";
	} else {
		print Rcmd "
cont <- NULL
";
	}
	if(defined $OPTIONS{m}) {
		print Rcmd "
for(i in 1:length(pw_tests)) cont <- c(cont, paste(pw_tests[[i]][2], \"-\", pw_tests[[i]][1], sep=\"\"))
";
	}
	if(!defined $OPTIONS{m} && !@add_cont){
		die("No Contrasts have been specified, you must at least either select multiple pairwise comparisons or specify a custom contrast\n");
	}
	print Rcmd "
cont <- makeContrasts(contrasts=cont, levels=design)
fit2 <- contrasts.fit(fit, cont)
fit2 <- eBayes(fit2)
";
} else {
	die("Anaysis type $OPTIONS{a} not found\n");

}
if($OPTIONS{a} ne "limma") {
	print Rcmd "
options(digits = 6)
tab <- NULL
for(i in names(tested)) {
	tab_tmp <- topTags(tested[[i]], n=Inf, adjust.method=\"$OPTIONS{f}\")[[1]]
	colnames(tab_tmp) <- paste(i, colnames(tab_tmp), sep=\":\")
	tab_tmp <- tab_tmp[tagnames,]
	if(is.null(tab)) {
		tab <- tab_tmp
	} else tab <- cbind(tab, tab_tmp)
}
tab <- cbind(Feature=rownames(tab), tab)
";
} else {
	print Rcmd "
tab <- NULL
options(digits = 6)
for(i in colnames(fit2)) {
	tab_tmp <- topTable(fit2, coef=i, n=Inf, sort.by=\"none\", adjust.method=\"$OPTIONS{f}\")
	colnames(tab_tmp)[-1] <- paste(i, colnames(tab_tmp)[-1], sep=\":\")
	if(is.null(tab)) {
		tab <- tab_tmp
	} else tab <- cbind(tab, tab_tmp[,-1])
}
";
}
print Rcmd "
write.table(tab, \"$OPTIONS{o}\", quote=F, sep=\"\\t\", row.names=F)
sink(type=\"message\")
sink()
";
close(Rcmd);
system("R --no-restore --no-save --no-readline < $OPTIONS{e}/r_script.R > $OPTIONS{e}/r_script.out");

open(HTML, ">$OPTIONS{h}");
print HTML "<html><head><title>EdgeR: Empirical analysis of digital gene expression data</title></head><body><h3>EdgeR Additional Files:</h3><p><ul>\n";
print HTML "<li><a href=MA_plots_normalisation.pdf>MA_plots_normalisation.pdf</a></li>\n";
print HTML "<li><a href=MDSplot.pdf>MDSplot.pdf</a></li>\n";
if($OPTIONS{a} eq "pw") {
	if(defined $OPTIONS{t}) {
		print HTML "<li><a href=Tagwise_Dispersion_vs_Abundance.pdf>Tagwise_Dispersion_vs_Abundance.pdf</a></li>\n";
	}
	print HTML "<li><a href=Smear_Plots.pdf>Smear_Plots.pdf</a></li>\n";
} elsif($OPTIONS{a} eq "glm" && $OPTIONS{d} eq "tag") {
	print HTML "<li><a href=Tagwise_Dispersion_vs_Abundance.pdf>Tagwise_Dispersion_vs_Abundance.pdf</a></li>\n";
	print HTML "<li><a href=Smear_Plots.pdf>Smear_Plots.pdf</a></li>\n";
} elsif($OPTIONS{a} eq "limma") {
	print HTML "<li><a href=LIMMA_MDS_plot.pdf>LIMMA_MDS_plot.pdf</a></li>\n";
	print HTML "<li><a href=LIMMA_voom.pdf>LIMMA_voom.pdf</a></li>\n";
}
print HTML "<li><a href=r_script.R>r_script.R</a></li>\n";
print HTML "<li><a href=r_script.out>r_script.out</a></li>\n";
print HTML "<li><a href=r_script.err>r_script.err</a></li>\n";
print HTML "</ul></p>\n";
close(HTML);

