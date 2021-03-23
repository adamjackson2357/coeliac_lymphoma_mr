# Clear variables and set the path
dev.off()
rm(list=ls())
path=dirname(rstudioapi::getActiveDocumentContext()$path)
setwd(path)

library(data.table)
library(openxlsx)
library(tidyverse)

### Creating the dataset

# read in the data
mydata=data.frame(fread("../data/ukb26390.csv", select=1))

# sample it for testing
mydata <- mydata %>%
  slice(1:100)

rownames(mydata)=mydata[,1]
diseases=c("cancer", "cardiovascular", "hypertension", "diabetes", "respiratory")
names(diseases)=c("Cancer", "Cardiovascular", "Hypertension", "Diabetes", "Respiratory")
x=matrix(0,nrow=nrow(mydata),ncol=length(diseases))
colnames(x)=names(diseases)
mydata=cbind(mydata,x)


### Disease definition

# Cancer (all except benign)

icd10_cancer=c(c(paste0("C0", 0:9), paste0("C", 10:14)), # lip, oral cavity, larynx
               paste0("C", 15:26), # digestive
               paste0("C", 30:39), # respiratory
               paste0("C", 40:41), # bone, cartilage
               paste0("C", 43:44), # skin
               paste0("C", 45:49), # mesothelial
               paste0("C", 50), # breast
               paste0("C", 51:58), # female genital
               paste0("C", 60:63), # male genital
               paste0("C", 64:68), # urinary
               paste0("C", 69:72), # nervous
               paste0("C", 73:75), # endocrin
               paste0("C", 76:80), # other sites
               paste0("C", 81:96), # lymphoid
               paste0("C", 97), # multiple
               paste0("D0", 0:9), # in situ
               paste0("D", 37:48)) # uncertain behaviour

icd9_cancer=c(paste0(140:149), # lip, oral cavity, larynx
              paste0(150:159), # digestive
              paste0(160:165), # respiratory
              paste0(170:176), # bone, cartilage
              paste0(179:189), # genital
              paste0(190:199), # unspecified
              paste0(200:209), # lymphoid
              paste0(230:234), # in situ
              paste0(235:238), # uncertain 
              paste0(239)) # unspecified


# Cardiovascular (all circulatory except hypertension)

icd10_cardiovascular=c(paste0("I0",0:2), # acute rheumatic fever
                       paste0("I0",5:9), # chronic rheumatic
                       paste0("I",20:25), # ischemic heart
                       paste0("I",26:28), # pulmonary
                       paste0("I",30:52), # other heart disease
                       paste0("I",60:69), # cerebrovascular
                       paste0("I",70:76), # arteries -- excluding I77.6 (in autoimmune)
                       paste0("I77", 0:5),
                       paste0("I77", 7:9),
                       paste0("I",78:79), 
                       paste0("I",81:82), # vein thrombosis
                       paste0("I",95:99)) # other circulatory

icd9_cardiovascular=c(paste0(390:392), # acute rheumatic fever
                      paste0(393:398), # chronic rheumatic
                      paste0(410:414), # ischemic heart
                      paste0(415:417), # pulmonary
                      paste0(420:429), # other heart disease
                      paste0(430:438), # cerebrovascular
                      paste0(440:445), # arteries -- excluding 446 and 4476 (in autoimmune)
                      paste0(447, 0:5), 
                      paste0(447, 7:9),
                      paste0(448:449), 
                      "453") #other vein thrombosis

# Hypertension

icd10_hypertension=c(paste0("I", 10:15))

icd9_hypertension=paste0(401:405)

# Diabetes
icd10_diabetes=c("E11", "E12", "E13", "E14")

# icd9_diabetes="250"
icd9_diabetes=c(paste0(paste0("250", 0:9), 0),paste0(paste0("250", 0:9), 2))

# Respiratory (all)
icd10_respiratory=c(paste0("J", 30:39), # other upper
                    paste0("J", 40:47), # chronic lower
                    paste0("J", 60:70), # lung 
                    paste0("J", 80:84), # other interstitium
                    paste0("J", 85:86), # lower
                    paste0("J", 90:94), # other pleura
                    paste0("J", 95:99)) # other respiratory

icd9_respiratory=c(paste0(470:478), # other upper
                   paste0(480:488), # influenza and pneumonia
                   paste0(490:496), # chronic
                   paste0(500:508), # lung
                   paste0(510:519)) # other

### HES: Hospital Episode Statistics

hes=data.frame(fread("../data/hesin_diag.txt"))

# reduce data down to the sampled eids
hes <- hes %>%
  semi_join(mydata, by="eid")

for (d in 1:length(diseases)){
  print(names(diseases)[d])
  disease=diseases[d]
  icd10_list=eval(parse(text=paste0("icd10_",disease)))
  icd9_list=eval(parse(text=paste0("icd9_",disease)))
  myeids=NULL
  
  # ICD10 codes
  pb=txtProgressBar(style=3)
  for (k in 1:length(icd10_list)){
    setTxtProgressBar(pb, k/length(icd10_list))
    tmp=as.character(hes$eid[grepl(paste0("^", icd10_list[k]), hes$diag_icd10)])
    myeids=c(myeids, tmp)
  }
  myeids=unique(myeids)
  cat("\n")
  print(length(myeids))
  
  # ICD9 codes
  pb=txtProgressBar(style=3)
  for (k in 1:length(icd9_list)){
    setTxtProgressBar(pb, k/length(icd9_list) )
    tmp=as.character(hes$eid[grepl(paste0("^", icd9_list[k]), hes$diag_icd9)])
    myeids=c(myeids, tmp)
  }
  myeids=unique(myeids)
  cat("\n")
  print(length(myeids))
  
  print(table(mydata[myeids,names(diseases)[d]]))
  mydata[myeids,names(diseases)[d]]=1
  print(table(mydata[,names(diseases)[d]]))
  cat("\n")
}

print(apply(mydata[,names(diseases)],2,sum))
saveRDS(mydata, "../data/disease_outcomes.rds")
