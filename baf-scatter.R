library(RColorBrewer)

names.kid <- c("11221K",   "11682K", "11232K", "11291K", "11689K", "11746K",   "11746K",   "9193K",   "9193K")
names.tum <- c("11221UpT", "11682T", "11232T", "11291T", "11689T", "11746_1T", "11746_2T", "9193_1T", "9193_2T")

pdf("/mnt/projects/zohre/results/baf-scatter.pdf", height=10, width=10)
for (i in 1:length(names.kid)) {
  kid <- read.delim(gzfile(paste0("/mnt/projects/zohre/results/cna/", names.kid[i], ".all.baf.bed.gz")))
  colnames(kid) <- c("chr", "start", "end", "name", "baf.kid")
  tum <- read.delim(gzfile(paste0("/mnt/projects/zohre/results/cna/", names.tum[i], ".all.baf.bed.gz")))
  colnames(tum) <- c("chr", "start", "end", "name", "baf.tum")
  
  m <- merge(kid, tum, by=c("chr", "start", "end", "name"))
  m <- m[(m$baf.kid >= 0.05 | m$baf.tum >= 0.05) & (m$baf.kid <= 0.95 | m$baf.tum <= 0.95),]
  
  smoothScatter(
    m[,c("baf.kid", "baf.tum")], 
    nrpoints=0, 
    colramp=function(n) colorRampPalette(c("blue", "red", "yellow", "white"))(n),
    main=paste0(names.kid[i], " vs. ", names.tum[i], "\n", "(n=", nrow(m), ")"),
    xlab=paste("BAF", names.kid[i]),
    ylab=paste("BAF", names.tum[i]),
    xlim=c(0, 1),
    ylim=c(0, 1)
  )
  abline(h=0.5, col="black", lty=2)
  abline(v=0.5, col="black", lty=2)
}
dev.off()
