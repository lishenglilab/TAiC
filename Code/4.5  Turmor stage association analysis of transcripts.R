## Perform the association analysis of transcripts with the tumor stages in the TCGA cohorts

# R version 4.0.2 (2020-06-22)
# tidyverse version=1.3.1

library(tidyverse)

# list expression files
files=list.files(path="/extraspace/projects/CCLE.RBP.transcript/results/TCGA/expression_mx/filter",pattern="*.txt")
files=grep("transcript",files,value = T)

for(i in files) {
  # import expression data
  exp=data.table::fread(paste0("/extraspace/projects/CCLE.RBP.transcript/results/TCGA/expression_mx/filter/",i),header =T,data.table=F,check.names=F)
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
  
  # extract barcodes with both expression and clinical data
  copatient=intersect(colnames(exp),cli$Id) 
  exp=exp[,copatient]
  rownames(cli)=cli$Id
  cli=cli[copatient,]
  
  data=cbind(cli,t(exp))
  
  index=data$stage %in% c("[Discrepancy]","[Not Applicable]","[Not Available]","[Unknown]","I/II NOS","IS","Stage 0","Stage X")
  data=data[!index,]
  data$stage=gsub("A$|B$|C$",'',data$stage)
  
  # Kruskal Wallis test
  kruskal_test <-apply(data[,5:ncol(data)], 2 , function(gene){
    test=kruskal.test(gene~stage,data)
    return(test$p.value)
  })
  
  kruskal_test=as.data.frame(kruskal_test)
  colnames(kruskal_test)="kruskal_test_pvalue"
  kruskal_test$id=rownames(kruskal_test)
  
  data$stage[data$stage=="Stage I"]=1
  data$stage[data$stage=="Stage II"]=2
  data$stage[data$stage=="Stage III"]=3
  data$stage[data$stage=="Stage IV"]=4
  
  genecor_function <- function(data1,data2){
    y=as.numeric(data1)
    rownames <- rownames(data2)
    dataframes=lapply(rownames,function(x){
      dd=cor.test(as.numeric(data2[x,]), y, method ="spearman")
      data.frame(id=x, cor=dd$estimate, p.value=dd$p.value)
    })
    return(do.call(rbind,dataframes))
  }
  
  # correlation analysis
  genecor <- genecor_function (data1 = data$stage, data2=t(data[,5:ncol(data)]))
  genecor$FDR <- p.adjust(genecor$p.value, method = "fdr")  
  genecor=left_join(genecor,kruskal_test,by="id")
  
  write.table(genecor,file=paste0("./TCGA_stage/",t,"_stage.txt"),sep="\t",row.names = F)
}
