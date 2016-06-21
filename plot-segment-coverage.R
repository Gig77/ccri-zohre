options(warn=1)
library(DNAcopy)
library(optparse)

option_list <- list(
		make_option("--patient", type="character", help="patient ID"),
		make_option("--tumor", type="character", help="tumor coverage data file"),
		make_option("--normal", type="character", help="remission coverage data file"),
		make_option("--output", type="character", help="output file name")
)
opt <- parse_args(OptionParser(option_list=option_list))

if (is.na(opt$patient)) stop("patient ID not specified")
if (is.na(opt$tumor)) stop("tumor coverage file not specified")
if (is.na(opt$normal)) stop("normal coverage file not specified")
if (is.na(opt$output)) stop("output file not specified")

# for test purposes
#opt <- data.frame(patient="1020540_dia", tumor="../reseq/cnv/segmented_coverage/1020540_Diagnosis.segmented-coverage.tsv", normal="../reseq/cnv/segmented_coverage/1020540_Remission.segmented-coverage.tsv", stringsAsFactors=F)
normal.chrs <- c("2")

t <- read.delim(opt$tumor, header=F, colClasses=c("factor", "integer", "numeric"))
n <- read.delim(opt$normal, header=F, colClasses=c("factor", "integer", "numeric"))
m <- merge(n, t, by=c("V1", "V2"))
m$ratio <- log((m[,"V3.y"] / sum(m[,"V3.y"])) / (m[,"V3.x"] / sum(m[,"V3.x"])), 2)
m$ratio <- m$ratio - ifelse(!is.null(normal.chrs), mean(m$ratio[which(m$V1 %in% normal.chrs & is.finite(m$ratio))]), 0)
m$ratio[m$ratio > 1.5] <- 1.5
m$ratio[m$ratio < -1.5] <- -1.5

set.seed(25)
CNA.object <-CNA(genomdat = m[,"ratio"], chrom = m[,"V1"], maploc = m[,"V2"], data.type = 'logratio')
CNA.smoothed <- smooth.CNA(CNA.object)
segs <- segment(CNA.smoothed, alpha = 0.01, verbose=0, min.width=2, undo.splits="sdundo", undo.SD=3)
#segs <- segment(CNA.smoothed, alpha = 0.01, verbose=0, min.width=2)

pdf(opt$output, width=30, paper="A4r")
par(mfrow=c(5,5), mar=c(0.5,2,1.5,1), oma=c(0, 0, 2, 0))
for (c in c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "X", "Y"))
{
	plot(m[m$V1==c,"V2"], m[m$V1==c,"ratio"], ylim=c(-1.5,1.5), main=c, cex=0.3, xaxt='n')
	abline(0, 0)
	for(i in which(segs$output$chrom==c)) {
		lines(c(segs$output$loc.start[i], segs$output$loc.end[i]),c(segs$output$seg.mean[i],segs$output$seg.mean[i]),type="l", col="orange", lty=1, lwd=3)
	}
}
mtext(opt$patient, outer=TRUE, cex=1.5)
dev.off()
