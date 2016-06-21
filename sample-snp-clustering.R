library(VariantAnnotation)

vcf <- readVcf("/mnt/projects/zohre/results/snps/allsamples.noMissing.vcf.gz", genome = "mm10")
gt <- geno(vcf)$GT ; dim(gt)
pl <- geno(vcf)$PL ; dim(pl)

# remove indels, multi-allelic sites, low quality sites, haplotype chromosomes, and non-informative sites (i.e. GT identical in all samples)
keep <- nchar(ref(vcf)) == 1 & nchar(unstrsplit(CharacterList(alt(vcf)))) == 1 
keep <- keep & qual(vcf) >= 900
keep <- keep & as.logical(seqnames(rowRanges(vcf)) %in% c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "X", "Y"))
keep <- keep & apply(gt, 1, function(x) length(unique(x))) > 1
gt <- gt[keep,] ; dim(gt)
pl <- pl[keep,] ; dim(pl)

# keep only sites with exquisite genotype quality in each sample
pl.min <- sapply(pl, FUN=min) ; attr(pl.min, "dim") <- attr(pl, "dim") ; attr(pl.min, "dimnames") <- attr(pl, "dimnames")
pl.med <- sapply(pl, FUN=median) ; attr(pl.med, "dim") <- attr(pl, "dim") ; attr(pl.med, "dimnames") <- attr(pl, "dimnames")
keep <- apply(pl.min, 1, max) == 0 & apply(pl.med, 1, min) > 70
gt <- gt[keep,] ; dim(gt)
pl <- pl[keep,] ; dim(pl)

# convert to numeric matrix
gt <- gsub("0/0", 0, gt)
gt <- gsub("0/1", 1, gt)
gt <- gsub("1/0", 1, gt)
gt <- gsub("1/1", 2, gt)
mode(gt) <- "numeric"

pdf("/mnt/projects/zohre/results/sample-snp-clustering.allchr.pdf")
pca <- prcomp(t(gt), center = TRUE, scale. = TRUE)
varianceExplained <- pca$sdev^2 / sum(pca$sdev^2) 
xspan <- max(pca$x[,1])-min(pca$x[,1])
plot(pca$x[,1], pca$x[,2], 
     cex=0.6, 
     pch=19,
     col=as.factor(gsub("(\\d+).*", "\\1", colnames(gt))),
     xlab=sprintf("PC1 (%.0f%%)", varianceExplained[1]*100), 
     ylab=sprintf("PC2 (%.0f%%)", varianceExplained[2]*100), 
     xlim=c(min(pca$x[,1])-xspan/10, max(pca$x[,1])+xspan/10), 
     main=sprintf("%d SNPs", nrow(gt)))
text(pca$x[,1], pca$x[,2], colnames(gt), cex=0.5)
dev.off()

