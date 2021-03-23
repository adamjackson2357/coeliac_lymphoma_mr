library(snpStats)
setwd("plink_files")

# Extract SNP names
bim=read.table("ukb_imp_merged.bim", stringsAsFactors=FALSE)
mysnps=bim[!duplicated(bim[,2]),2]

# Rename participants
fam=read.table("ukb_imp.fam")
fam[,1]=fam[,2]=1:nrow(fam)
write.table(fam, "ukb_imp_merged.fam", quote=FALSE, row.names=FALSE, col.names=FALSE)

# Read plink data
mydata=read.plink("ukb_imp_merged", select.snps=mysnps, na.strings="-9")

# Extract genotype data
genodata=mydata$genotypes
genodata=data.frame(genodata)

# Rename participants
fam=read.table("ukb_imp.fam")
rownames(genodata)=fam[,1]

# Save genetic data
saveRDS(genodata, "../genetic_data_extracted.rds")
