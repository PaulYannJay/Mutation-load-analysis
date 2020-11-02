library(SNPRelate, quietly=T)

options(scipen=999)
args = commandArgs(trailingOnly=TRUE)
file=args[1]
chr=args[2]
deb=args[3]
fin=args[4]
out=args[5]
snpgdsVCF2GDS(file, paste(file, ".gds", sep=""),  method="biallelic.only",verbose=F)
genofile <- snpgdsOpen(paste(file, ".gds", sep=""), allow.duplicate=TRUE)
a=snpgdsSummary(genofile)$sample.id
ccm_pca<-snpgdsPCA(genofile, autosome.only=FALSE,  sample.id=a, verbose=F)

table=data.frame(ccm_pca$eigenvect[,1], ccm_pca$eigenvect[,2], ccm_pca$eigenvect[,3],ccm_pca$eigenvect[,4],ccm_pca$eigenvect[,5],ccm_pca$eigenvect[,6],ccm_pca$sample.id, chr, deb, fin)
write.table(table, out, quote=F, row.name=F, append=TRUE, col.name=F)
