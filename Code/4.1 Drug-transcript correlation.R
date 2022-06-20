## 

# R version 4.0.2 (2020-06-22)
# tidyverse version=1.3.1

library(tidyverse)
args=commandArgs(T)

## import data
file=paste0("/home/Wei.Hu/Task_RBP_transcript/transcript_split/transcript",args[1],".txt")
exp=read.table(file,header =T,sep="\t",check.names=F)
colnames(exp)=gsub("_.*","",colnames(exp))

auc=read.table("/home/Wei.Hu/Task_RBP_transcript/CTRP_AUC_matrix.txt",sep="\t",header =T,check.names=F)

# correlation analysis
for (i in rownames(auc)){
  print(i)
  auc.i=auc[i,]
  auc.i=auc.i[,!is.na(auc.i)]
  a=intersect(colnames(exp),colnames(auc.i))
  exp1=exp[,a]
  auc.i=auc.i[,a]
  data.i=exp1[apply(exp1,1, function(x) sum(x>=0.1)/ncol(exp1) >= 0.2),]
  rownames <- rownames(data.i)
  dataframes=lapply(rownames,function(x){
    dd=cor.test(as.numeric(data.i[x,]), as.numeric(auc.i), method ="spearman")
    data.frame(drug=i, Gene_2=x, cor=dd$estimate, p.value=dd$p.value)
  })
  genecor <- do.call(rbind, dataframes)
  
  write.table(genecor, paste0("/home/Wei.Hu/Task_RBP_transcript/CTRP-transcript/Cor_transcript_",i,"_",str_extract(file, "\\d+"),".txt"), sep = "\t", quote = F, row.names = F)
}
