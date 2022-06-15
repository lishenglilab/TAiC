library(tidyverse)
library(rtracklayer)

## generate bed file of exons of transcripts
# import data
AnnoData <- import("/home/weihu/RBP_transcript_drug_project/data/gtf/CCLE_1017_merged_filtered_0.1.gtf")
table(AnnoData$type) # check type

# extract exon information
index <- which(AnnoData$type == "exon") # extract indexes of exons
ccle_exon <- data.frame(Type = AnnoData$type[index],
                        Gene_Id = AnnoData$gene_id[index],
                        Gene_Name = AnnoData$gene_name[index],
                        Transcript_Id = AnnoData$transcript_id[index],
                        AnnoData@seqnames[index],
                        AnnoData@ranges[index],
                        AnnoData@strand[index])
colnames(ccle_exon)[5:9] <- c("Chr","Start","End","Width","Strand")

# generate bed file
ccle_exon_bed <- data.frame(ccle_exon$Chr, ccle_exon$Start, ccle_exon$End, paste0(ccle_exon$Transcript_Id, "_", 
                                                                                  ccle_exon$Chr, ":",
                                                                                  ccle_exon$Start, "|",
                                                                                  ccle_exon$End, ":",
                                                                                  ccle_exon$Strand),
                            0, ccle_exon$Strand)


# export
setwd("/home/weihu/RBP_transcript_drug_project/data/bed")
write.table(ccle_exon_bed, "CCLE_exon.bed", sep = "\t", row.names = F, col.names = F, quote = F)

## DESeq2 analysis of controls and corresponding knockdown replicates
# use HepG2 as exsample
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
