library(optparse)
library(DNAcopy)

option_list <- list(
		make_option("--sample-id", type="character", help="Sample ID"),
		make_option("--lrr", type="character", help="file containing LRR values"),
		make_option("--baf", type="character", help="file containing BAF values"),
		make_option("--plot-output-file", type="character", help="name of output file for chromosome plots (PDF format)"),
		make_option("--loh-segments-output-file", type="character", help="name of output file for LOH segments")
)
opt <- parse_args(OptionParser(option_list=option_list))

#opt <- data.frame('sample-id' = "11221UpT", 'lrr' = "/mnt/projects/zohre/results/cna/11221Up.genome-coverage.lrr.bed", 'baf' = "/mnt/projects/zohre/results/cna/11221UpT.genome-coverage.baf.bed.gz", 'plot-output-file' = "/mnt/projects/zohre/results/cna/allpatients.snp-profile.pdf", 'loh-segments-output-file' = "/mnt/projects/zohre/results/cna/allpatients.loh-segments.tsv", stringsAsFactors=F, check.names=F)   

# read coverage data
lrr <- read.delim(opt$'lrr', header=F, colClasses=c("factor", "integer", "integer", "character", "numeric"))
colnames(lrr) <- c("chr", "start", "end", "name", "ratio")
baf <- read.delim(gzfile(opt$'baf'), header=F, colClasses=c("factor", "integer", "integer", "character", "numeric"))
colnames(baf) <- c("chr", "start", "end", "name", "baf")

# segment coverage
cna.cov <-CNA(genomdat = lrr[,"ratio"], chrom = lrr[,"chr"], maploc = lrr[,"start"], data.type = 'logratio')
cna.smooth <- smooth.CNA(cna.cov)
seg.cov <- segment(cna.smooth, alpha = 0.001, verbose=0, min.width=2)

# segment mirrored AF
baf$maf <- abs(baf$baf - 0.5) + 0.5
cna.maf <-CNA(genomdat = baf$maf, chrom = baf$chr, maploc = baf$start, data.type = 'logratio')
cna.smooth <- smooth.CNA(cna.maf)
seg.maf <- segment(cna.smooth, alpha = 0.001, verbose=0, min.width=2)

# write loh segments
#if (!is.null(opt$'loh-segments-output-file')) {
#	lohs <- seg.maf$output
#	names(lohs)[1] <- "sample"
#	lohs$sample <- opt$'sample-id'
#	write.table(lohs, file=opt$'loh-segments-output-file', row.names=F, col.names=T, sep="\t", quote=F)
#}

# plot function coverage
plot.coverage <- function(data, segs) {
	plot(data[,"start"], data[,"ratio"], ylim=c(-1.5,1.5), xaxt='n', yaxt='n', type="n")
	rect(par("usr")[1], par("usr")[3], par("usr")[2], par("usr")[4], col=rgb(0.99,0.99,0.99))
	points(data[,"start"], data[,"ratio"], ylim=c(-1.5,1.5), cex=0.3, col="darkblue", lwd=0.5)
	axis(2, at=c(-1, 0, 1))
	abline(0, 0)
	for (i in which(abs(segs$seg.mean) >= 0.2)) {
		lines(c(segs$loc.start[i], segs$loc.end[i]),c(segs$seg.mean[i],segs$seg.mean[i]),type="l", col="orange", lty=1, lwd=2, yaxt='n')
	}
	abline(h=c(-1, log2(2/3), log2(3/2), 1), lty=2)
}

# plot function maf
plot.maf <- function(data, segs, xlim, title=title) {
	if (nrow(data) > 0) {
		plot(data[,"start"], data[,"baf"], ylim=c(0,1), xlim=xlim, cex=0.3, xaxt='n', yaxt='n', main=title, type='n')
		rect(par("usr")[1], par("usr")[3], par("usr")[2], par("usr")[4], col=rgb(0.99,0.99,0.99))
		points(data[,"start"], data[,"baf"], ylim=c(0,1), cex=0.3, col="darkblue", lwd=0.5)
	}
	else {
		plot(0, 0, ylim=c(0,1), cex=0.3, xaxt='n', yaxt='n', main=title)
		rect(par("usr")[1], par("usr")[3], par("usr")[2], par("usr")[4], col=rgb(0.97,0.97,0.97))
	}
	for(i in which(segs$seg.mean >= 0.6 & segs$seg.mean <= 0.9)) {
		lines(c(segs$loc.start[i], segs$loc.end[i]),c(segs$seg.mean[i],segs$seg.mean[i]),type="l", col="orange", lty=1, lwd=2)
		lines(c(segs$loc.start[i], segs$loc.end[i]),c(1-segs$seg.mean[i],1-segs$seg.mean[i]),type="l", col="orange", lty=1, lwd=2)
	}
	axis(2, at=c(0, 0.25, 1/3, 0.5, 2/3, 0.75, 1), labels=c("0", "", "", "0.5", "", "", "1"))
	abline(h=0.5)
	abline(h=c(1/3, 0.25, 2/3, 0.75), lty=2)
}

pdf(opt$'plot-output-file', width=15, height=10)
layout(matrix(c(1,3,5,7,9,2,4,6,8,10,11,13,15,17,19,12,14,16,18,20,21,23,25,27,29,22,24,26,28,30,31,33,35,37,39,32,34,36,38,40,41,43,45,47,49,42,44,46,48,50),10,5,byrow=T));
par(oma=c(0, 0, 2, 0))
for (c in c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "X", "Y")) {
	par(mar=c(0,2,1,0.5))
	plot.maf(baf[baf$chr==c,], seg.maf$output[seg.maf$output$chrom==c,], xlim=c(0,max(lrr$start[lrr$chr==c])), title=c)
	par(mar=c(0.2,2,0.2,0.5))
	plot.coverage(lrr[lrr$chr==c,], seg.cov$output[seg.cov$output$chrom==c,])
}
mtext(opt$'sample-id', outer=TRUE, cex=1.5)
dev.off()
