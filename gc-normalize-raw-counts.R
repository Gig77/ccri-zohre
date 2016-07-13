counts <- read.delim("/mnt/projects/zohre/results/cna/11232K.genome-coverage.250000.tsv", header=F, colClasses=c("factor", "integer", "numeric"))
names(counts) <- c("chr", "bin", "count")

gc <- read.delim("/mnt/projects/zohre/results/cna/gcPerBin.250000.tsv", header=F, colClasses=c("factor", "integer", "numeric"))
names(gc) <- c("chr", "bin", "gc")

counts <- merge(counts, gc, by=c("chr", "bin"))

# remove freak bins and sort
counts <- counts[counts$chr != "Y" & counts$count > 2^16 & counts$gc > 0.35,]
counts <- counts[order(counts$chr, counts$bin),]

counts$gc.sq <- counts$gc^2
counts$gc.cub <- counts$gc^3
counts$gc.forth <- counts$gc^4
counts$gc.fifth <- counts$gc^5
counts$log2 <- log2(counts$count+1)

fit <- lm(log2 ~ gc + gc.sq + gc.cub, data=counts) ; summary(fit)
#fit <- lm(log2 ~ gc + gc.sq, data=counts) ; summary(fit)
#fit <- lm(log2 ~ gc, data=counts) ; summary(fit)

# check fit visually
counts.sorted <- counts[order(counts$gc),]
plot(log2~gc, data=counts.sorted, cex=0.2)
lines(counts.sorted$gc, predict(fit, counts.sorted), col="orange", lwd=5)

counts$log2.corr <- counts$log2
counts$log2.corr <- counts$log2.corr - fit$coeff["gc"] * (counts$gc-mean(counts$gc))
counts$log2.corr <- counts$log2.corr - fit$coeff["gc.sq"] * (counts$gc.sq-mean(counts$gc.sq))
counts$log2.corr <- counts$log2.corr - fit$coeff["gc.cub"] * (counts$gc.cub-mean(counts$gc.cub))

chrom <- "5"
counts.chr <- counts[counts$chr == chrom,]

smooth <- function(x,n=5){filter(x,rep(1/n,n), sides=2)}

counts.chr$log2.smooth <- smooth(counts.chr$log2, 20)
counts.chr$log2.corr.smooth <- smooth(counts.chr$log2.corr, 20)
counts.chr$gc.smooth <- smooth(counts.chr$gc, 20)
counts.chr <- counts.chr[complete.cases(counts.chr),]

plot(1:nrow(counts.chr), counts.chr$log2.smooth-mean(counts.chr$log2.smooth), type="l", lwd=2, ylab="log2(count)", xlab=paste0("bin chr", chrom), col="red", xlim=c(20,nrow(counts.chr)-20))
lines(1:nrow(counts.chr), counts.chr$log2.corr.smooth-mean(counts.chr$log2.corr.smooth), lwd=2, col="blue")
lines(1:nrow(counts.chr), (counts.chr$gc.smooth-mean(counts.chr$gc.smooth))*10, lwd=1, col="black", lty=3)
legend("bottomright", c("counts", "gc"), fill=c("red", "blue"))
