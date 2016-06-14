# source("https://bioconductor.org/biocLite.R") ; biocLite("rtracklayer")
library(rtracklayer)
library(reshape2)
library(tidyr)
library(GenomicRanges)

# mm9 -> mm10
t <- read.delim("/data_synology/max/zohre/agilent/S0276129_Covered.bed", stringsAsFactors = F, check.names = F, header = F)
names(t) <- c("chr", "start", "end", "ann")
t$chr[t$chr=="MT"] <- "M"
t$chr <- paste0("chr", t$chr)
t <- makeGRangesFromDataFrame(t, keep.extra.columns = T)

# liftover hg18 to hg19
chain <- import.chain("/data_synology/max/zohre/agilent/mm9ToMm10.over.chain")
t.mm10 <- liftOver(t, chain)@unlistData
t.mm10 <- as.data.frame(t.mm10)
t.mm10$seqnames <- gsub("^chr", "", t.mm10$seqnames)
t.mm10$seqnames[t.mm10$seqnames=="M"] <- "MT"

write.table(t.mm10[,c("seqnames", "start", "end", "ann")], "/data_synology/max/zohre/agilent/S0276129_Covered.mm10.bed", col.names = F, row.names = F, quote = F, sep="\t")
