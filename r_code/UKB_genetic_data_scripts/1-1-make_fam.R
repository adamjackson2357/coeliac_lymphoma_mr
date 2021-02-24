library(data.table)

mydata=readRDS("../Data/Covid_dataset_severity.rds")

case=rep(0,nrow(mydata)) # 0 for excluded
case[which(mydata$covid_infected==1)]=1 # 1 for controls
case[which(mydata$severe_covid==1)]=2 # 2 for cases
names(case)=rownames(mydata)

fam=data.frame(fread("../Data/ukb19266_cal_chr21_v2_s488264.fam"))
rownames(fam)=fam[,1]
ids=intersect(rownames(fam), names(case))
fam[,6]=0
fam[ids,6]=case[ids]
print(table(fam[,6]))

write.table(fam, "../Data/ukb_geno_m1.fam", row.names=FALSE, col.names=FALSE, quote=FALSE)
