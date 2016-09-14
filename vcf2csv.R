#source("https://bioconductor.org/biocLite.R")
#biocLite("VariantAnnotation")

library(VariantAnnotation)

vcf <- readVcf("/home/anduril/zohre/results/snpeff/test.varscan.dbsnp.snpeff.vcf", genome = "mm10")
samples(header(vcf))

geno(vcf)$FREQ

ann2df <- function(snpEffAnn) {
  df <- data.frame(do.call("rbind", strsplit(snpEffAnn, '\\|')))
  names(df) <- c("Allele", "Annotation", "Annotation_Impact", "Gene_Name", "Gene_ID", "Feature_Type", "Feature_ID",
                 "Transcript_BioType", "Rank", "HGVS.c", "HGVS.p", "cDNA.pos_cDNA.length", "CDS.pos_CDS.length",
                 "AA.pos_AA.length", "Distance", "Info")
  df
}

annotations <- lapply(info(vcf)$ANN, FUN=ann2df)

