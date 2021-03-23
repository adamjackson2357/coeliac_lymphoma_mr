library(data.table)

# Loading the data
covars=readRDS("../Data/Covid_dataset_severity.rds")

# Preparing dummy variables
covars$eid=as.numeric(rownames(covars))
covars$Batch=as.factor(as.character(covars$Batch))
conf_list=c("eid", colnames(covars)[1:which(colnames(covars)=="PC10")])
conf_list=conf_list[!grepl("airpollution", conf_list)]
print(conf_list)
mymodel=model.matrix(as.formula(paste0("~", paste(conf_list, collapse="+"))), data=covars)
print(head(mymodel))
mymodel=mymodel[,-1]
colnames(mymodel)[1]="IID"
mymodel=cbind(FID=mymodel[,1], mymodel)
colnames(mymodel)=gsub(" ", "", colnames(mymodel))

# Saving confounders file
write.table(mymodel, "../Data/confounders.txt", row.names=FALSE, col.names=TRUE, quote=FALSE, sep=" ")




