options(warn=1)
library(DNAcopy)
library(optparse)

option_list <- list(
    make_option("--patient", type="character", help="patient ID"),
		make_option("--tumor", type="character", help="tumor coverage data file"),
		make_option("--normal", type="character", help="remission coverage data file"),
		make_option("--gccontent", type="character", help="GC content per bin (optional)"),
		make_option("--output", type="character", help="output file name")
)
opt <- parse_args(OptionParser(option_list=option_list))

if (is.na(opt$patient)) stop("patient ID not specified")
if (is.na(opt$tumor)) stop("tumor coverage file not specified")
if (is.na(opt$normal)) stop("normal coverage file not specified")
if (is.na(opt$output)) stop("output file not specified")

# for test purposes
#opt <- data.frame(patient="11291", tumor="/mnt/projects/zohre/results/cna/11291T.genome-coverage.tsv", normal="/mnt/projects/zohre/results/cna/11291K.genome-coverage.tsv", gccontent="/mnt/projects/zohre/results/cna/gcPerBin.tsv", output="/mnt/projects/zohre/results/cna/11291.genome-coverage.pdf", stringsAsFactors=F)
normal.chrs <- c("2")

t <- read.delim(opt$tumor, header=F, colClasses=c("factor", "integer", "numeric"))
names(t) <- c("chr", "bin", "tumor")
n <- read.delim(opt$normal, header=F, colClasses=c("factor", "integer", "numeric"))
names(n) <- c("chr", "bin", "normal")
m <- merge(n, t, by=c("chr", "bin"))

# remove bins with very low read counts in tumor and normal
#m <- m[m$tumor > 10000 | m$normal > 10000,]

m$ratio <- log2((m[,"tumor"]+0.1 / sum(m[,"tumor"])) / (m[,"normal"]+0.1 / sum(m[,"normal"])))
m$ratio <- m$ratio - ifelse(!is.null(normal.chrs), mean(m$ratio[which(m$chr %in% normal.chrs & is.finite(m$ratio))]), 0)
m$ratio[m$ratio > 1.5] <- 1.5
m$ratio[m$ratio < -1.5] <- -1.5

if (!is.null(opt$gccontent)) {
  gc <- read.delim(opt$gccontent, header=F, colClasses=c("factor", "integer", "numeric"))
  names(gc) <- c("chr", "bin", "gc")
  gc$gc.sq <- gc$gc^2
  gc$gc.cub <- gc$gc^3
  m <- merge(m, gc, by=c("chr", "bin"))
  fit.ratio <- lm(m$ratio ~ gc + gc.sq + gc.cub, data=m) ; summary(fit.ratio)
  m$ratio <- m$ratio - fit.ratio$coeff[2] * (m$gc-mean(m$gc))
  m$ratio <- m$ratio - fit.ratio$coeff[3] * (m$gc.sq-mean(m$gc.sq))
  m$ratio <- m$ratio - fit.ratio$coeff[4] * (m$gc.cub-mean(m$gc.cub))
  m$ratio[m$ratio > 1.5] <- 1.5
  m$ratio[m$ratio < -1.5] <- -1.5
}
m <- m[order(m$chr, m$bin),]

set.seed(25)
CNA.object <-CNA(genomdat = m[,"ratio"], chrom = m[,"chr"], maploc = m[,"bin"], data.type = 'logratio')
CNA.smoothed <- smooth.CNA(CNA.object)
segs <- segment(CNA.smoothed, alpha = 0.01, verbose=0, min.width=2, undo.splits="sdundo", undo.SD=3)
#segs <- segment(CNA.smoothed, alpha = 0.01, verbose=0, min.width=2)

smooth <- function(x,n=5){filter(x,rep(1/n,n), sides=2)}

pdf(opt$output, width=30, paper="A4r")
par(mfrow=c(5,5), mar=c(0.5,2,1.5,1), oma=c(0, 0, 2, 0))
for (c in c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "X", "Y"))
{
	plot(m[m$chr==c,"bin"], m[m$chr==c,"ratio"], ylim=c(-1.5,1.5), main=c, cex=0.3, xaxt='n')
	abline(0, 0)
	for(i in which(segs$output$chrom==c)) {
		lines(c(segs$output$loc.start[i], segs$output$loc.end[i]),c(segs$output$seg.mean[i],segs$output$seg.mean[i]),type="l", col="orange", lty=1, lwd=3)
	}
	gc.smoothed <- smooth((m[m$chr==c,"gc"]-mean(m[,"gc"]))*10, 30)
	lines(m[m$chr==c,"bin"], gc.smoothed, lwd=0.5, col="blue")
}
mtext(opt$patient, outer=TRUE, cex=1.5)
dev.off()
