## Perform differential expression analysis of transcripts in 17 cancer types of the TCGA cohorts

# R version 4.0.2 (2020-06-22)
# tidyverse version=1.3.1

library(tidyverse)

#list normal files
files=list.files(path="/home/weihu/pancancer_transcript_project/Data/TCGA/TCGA_exp_normal_all",pattern="*.txt")
files=grep("transcript",files,value = T)

#list tumor files
files_cancer=list.files(path="/mnt/extraspace/projects/CCLE.RBP.transcript/results/TCGA/expression_mx/filter",pattern="*.txt")
files_cancer=grep("transcript",files_cancer,value = T)

for(i in files){
  
  cancer=toupper(str_extract(i,"(?<=mx_).*(?=.txt)"))
  exp_nor=data.table::fread(paste0("/home/weihu/pancancer_transcript_project/Data/TCGA/TCGA_exp_normal_all/",i),
                            header =T,data.table=F,check.names=F)
  
  files_cancer_i <- grep(cancer, files_cancer,value = T)
  exp_tumor=data.table::fread(paste0("/mnt/extraspace/projects/CCLE.RBP.transcript/results/TCGA/expression_mx/filter/",files_cancer_i),
                              header =T,data.table=F,check.names=F)
  
  exp_tumor=left_join(exp_tumor,exp_nor,by="Transcript")
  rownames(exp_tumor)=exp_tumor[,1]
  exp_tumor=exp_tumor[,-1]
  
  group_list=ifelse(as.numeric(substr(colnames(exp_tumor),14,15)) < 10,"Tumor","Normal")
  normal=exp_tumor[,group_list=="Normal"]
  colnames(normal)=gsub("....$","",colnames(normal))
  tumor=exp_tumor[,group_list=="Tumor"]
  colnames(tumor)=gsub("....$","",colnames(tumor))
  
  copatient=intersect(colnames(normal),colnames(tumor))
  normal=normal[,copatient]
  tumor=tumor[,copatient]
  
  round2 = function(x, n) {
    posneg = sign(x)
    z = abs(x)*10^n
    z = z + 0.5
    z = trunc(z)
    z = z/10^n
    z*posneg
  }
  
  # remove low abundance transcripts
  normal_filtered=normal[rowSums(normal >0.1) >= round2(ncol(normal)*0.25,0),]
  tumorl_filtered=tumor[rowSums(tumor >0.1) >= round2(ncol(tumor)*0.25,0),]
  
  cotranscript=union(rownames(normal_filtered),rownames(tumorl_filtered))
  
  normal=normal[cotranscript,]
  tumor=tumor[cotranscript,]
  
  #set TPM < 0.1 as 0
  normal[normal < 0.1] <- 0
  tumor[tumor < 0.1] <- 0
  
  print(paste0(cancer,": ",ncol(normal)," ",ncol(tumor)))
  
  if(identical(rownames(normal),rownames(tumor))){
    
    deg=data.frame()
    for (t in rownames(normal)){
      group1=normal[t,]
      group2=tumor[t,]
      logFC= log2(mean(as.numeric(group2))/mean(as.numeric(group1)))
      p <- t.test(log2(as.numeric(group1)+1),log2(as.numeric(group2)+1),
                  paired =T # paired t test
      )$p.value
      result=data.frame(id=t,normal=sum(group1)/ncol(group1),tumor=sum(group2)/ncol(group2),log2FoldChange=logFC,p.value=p)
      deg=rbind(deg,result)
    }
    deg$FDR <- p.adjust(deg$p.value, method = "fdr")
    
    logFC_t=log2(1.5)
    P.Value_adj = 0.05
    k1 = (deg$FDR < P.Value_adj)&(deg$log2FoldChange < -logFC_t)
    k2 = (deg$FDR < P.Value_adj)&(deg$log2FoldChange > logFC_t)
    Change = ifelse(k1,"Down",ifelse(k2,"Up","Stable"))
    deg$Change <- Change
    
    write.table(deg,file=paste0("/home/weihu/pancancer_transcript_project/A2_Novel_transcript/Data/TCGA_DEG/",cancer,"_DEG.txt"),sep="\t",row.names = F,quote = F)
  }
  
}
