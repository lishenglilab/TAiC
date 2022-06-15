library(tidyverse)
library(survival)
library(survminer)
args=commandArgs(T)

# list expression files
files=list.files(path=args[2],pattern="*.txt")
files=grep("transcript",files,value = T)
i=grep(args[1],files,value = T)

# import expression data
exp=data.table::fread(paste0(args[2],"/",i),header =T,data.table=F,check.names=F)
rownames(exp)=exp[,1]
exp=exp[,-1]
colnames(exp)=gsub("....$","",colnames(exp))

t=str_extract(i,"(?<=mx_).*(?=_filtered)")
print(paste0("Begin:",t))

# remove duplicated columns
if(sum(table(colnames(exp)) > 1)>0){
  exp=exp[,!duplicated(colnames(exp))]
  print(paste0(t," duplicated"))
}

# import clincial data
pan_cli=data.table::fread("TCGA-CDR-SupplementalTableS1.csv",header =T,data.table=F,check.names=F)

# extract clinical information
cli= pan_cli %>% 
  filter(type==t)
cli=cli[,c("bcr_patient_barcode","ajcc_pathologic_tumor_stage","OS","OS.time")]
colnames(cli)[1]="Id"
colnames(cli)[2]="stage"
colnames(cli)[3]="fustat"
colnames(cli)[4]="futime"

cli=subset(cli,futime>=30)
cli$futime=cli$futime/365
if(sum(is.na(cli$fustat))>0|sum(is.na(cli$futime))>0){
  cli=cli[!is.na(cli$fustat) & !is.na(cli$futime),]
  print(paste0(t," cli NA"))
}

# extract barcodes with both expression and clinical data
copatient=intersect(colnames(exp),cli$Id) 
exp=exp[,copatient]
rownames(cli)=cli$Id
cli=cli[copatient,]

print(paste0("clinical: ",nrow(cli)))

# Cox regression survival analysis
mySurv=with(cli,Surv(futime, fustat))
cox_result_continuous <-apply(exp, 1 , function(gene){
  cox=coxph(mySurv ~ gene)
  coxSummary = summary(cox)

  tmp <- cbind(HR=coxSummary$conf.int[,"exp(coef)"],
               HR.95L=coxSummary$conf.int[,"lower .95"],
               HR.95H=coxSummary$conf.int[,"upper .95"],
               pvalue=coxSummary$coefficients[,"Pr(>|z|)"])

  return(tmp[1,])

})
cox_result_continuous=t(cox_result_continuous)

cox_result_classified <-apply(exp, 1 , function(gene){
  group=ifelse(gene>median(gene),'high','low') 
  survival_dat <- data.frame(group= group,
                             stringsAsFactors = F)
  cox=coxph(mySurv ~ group, data =  survival_dat)
  coxSummary = summary(cox)

  tmp <-cbind(HR=1/coxSummary$conf.int[,"exp(coef)"],
              HR.95L=1/coxSummary$conf.int[,"upper .95"],
              HR.95H=1/coxSummary$conf.int[,"lower .95"],
              pvalue=coxSummary$coefficients[,"Pr(>|z|)"])
  return(tmp[1,])

})
cox_result_classified=t(cox_result_classified)

# kaplan-meier survival analysis
log_rank_p <-apply(exp, 1 , function(gene){
  group=ifelse(gene>median(gene),'high','low') 
  survival_dat <- data.frame(group= group,
                             stringsAsFactors = F)
  data.survdiff=survdiff(mySurv~group,data=survival_dat)
  p.val = 1 - pchisq(data.survdiff$chisq, length(data.survdiff$n) - 1)
  return(p.val)

})
log_rank_result=as.data.frame(cbind(cox_result_classified[,1:3],pvalue=log_rank_p))

colnames(cox_result_continuous)=c("HR_Cox","HR.95L_Cox","HR.95H_Cox","pvalue_Cox")
colnames(log_rank_result)=c("HR_KM","HR.95L_KM","HR.95H_KM","pvalue_KM")
survival_result=cbind(cox_result_continuous,log_rank_result)
write.table(survival_result,file=paste0("./TCGA_survival/",t,"_survival.txt"),sep="\t")
