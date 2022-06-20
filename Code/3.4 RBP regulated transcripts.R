## Get RBP regulated transcript from RBP KD-RNA-seq data

# R version 4.0.2 (2020-06-22)
# tidyverse version=1.3.1
# DESeq2 version=1.3.1

## DESeq2 analysis of controls and corresponding knockdown replicates
# use HepG2 as exsample
library(tidyverse)
library(DESeq2)

# import data
exp=data.table::fread("/home/weihu/RBP_transcript_drug_project/4_RBP_regulation/data/count_tpm_matrix/shRNAdeq_count_hepG2.txt",header =T,data.table=F,check.names=F)
rownames(exp)=exp$transcript_id
exp=exp[,-1]

for(i in unique(gsub("_.*","",colnames(exp)[-c(ncol(exp)-1,ncol(exp))]))){
  # project name
  project=i
  
  groupA="HepG2_control"
  groupB=i
  size=2 # sample size
  
  data_count=exp %>% 
    dplyr::select(starts_with(groupA),starts_with(paste(i,"_rep",sep="")))
  group_list=c(rep(groupA,size),rep(groupB,size))
  condition <- factor(group_list,levels = (c(groupA,groupB)))
  
  colData <- data.frame(row.names = colnames(data_count), 
                        condition=condition)
    
  dds <- DESeqDataSetFromMatrix(countData = data_count, colData = colData, design = ~ condition)
  dds <- DESeq(dds)
  res <- results(dds, contrast = c("condition",rev(levels(condition))))
  summary(res)
  resordered <- res[order(res$padj),]
  deg <- as.data.frame(resordered)
  
  logFC_t=log2(1.5)
  P.Value_adj = 0.05
  k1 = (deg$padj < P.Value_adj)&(deg$log2FoldChange < -logFC_t)
  k2 = (deg$padj < P.Value_adj)&(deg$log2FoldChange > logFC_t)
  Change = ifelse(k1,"Down",ifelse(k2,"Up","Stable"))
  deg$Change <- Change
  deg$transcript_id=rownames(deg)
  
  deg=deg %>% 
    dplyr::select(transcript_id,everything()) %>% 
    arrange(desc(abs(log2FoldChange)))
  write.table(deg,file=paste0("/home/weihu/RBP_transcript_drug_project/4_RBP_regulation/data/shRNAseq_DEseq2_HepG2/HepG2-",project,"_deg.txt"),sep="\t",row.names = F)
}

## select significant differentially expressed transcripts with RBP binding sites
# use HepG2 as exsample
files=list.files(path="/home/weihu/RBP_transcript_drug_project/4_RBP_regulation/data/shRNAseq_DEseq2_HepG2",pattern="*.txt")
RBP=gsub("HepG2-","",files)
RBP=gsub("_deg.txt","",RBP)

for(i in RBP){
  print(i)
  deg=data.table::fread(paste0("/home/weihu/RBP_transcript_drug_project/4_RBP_regulation/data/shRNAseq_DEseq2_HepG2/","HepG2-",i,"_deg.txt"),header =T,data.table=F,check.names=F)
  deg=deg[,-c(2, 4, 5)]
  deg=deg[deg$Change %in% c("Up","Down"),]
  
  eclip=data.table::fread(paste0("/home/weihu/RBP_transcript_drug_project/4_RBP_regulation/data/HepG2_intersect_bed/",i,"_exon_eCLIP_intersect"),header =F,data.table=F,check.names=F)
  
  deg_eclip=deg[deg$transcript_id %in% gsub("_.*", "", eclip$V4),]
  deg_eclip$RBP=i
  deg_eclip$cell_type="HepG2"
  
  write.table(deg_eclip,file=paste0("/home/weihu/RBP_transcript_drug_project/4_RBP_regulation/data/DEseq2_eCLIP_HepG2/HepG2_deg_eclip_",i,".txt"),sep="\t",row.names=F,quote = F)
}
