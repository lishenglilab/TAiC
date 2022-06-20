## retrieve the genomic positions and gene types of transcripts

# R version 4.0.2 (2020-06-22)
# tidyverse version=1.3.1
# rtracklayer version=1.5.0

library(tidyverse)
library(rtracklayer)

# import CCLE GTF file
merged = import("./data/gtf/CCLE_1017_merged_filtered.gtf")
index = which(merged$type == "transcript") # extract transcript index
# extract transcript information
merged.info=data.frame(Gene_Id = merged$gene_id[index],
                      Transcript_Id = merged$transcript_id[index],
                      Gene_Name=merged$gene_name[index],
                      Ref_Gene_Id=merged$ref_gene_id[index])

# import gencode.v35 GTF file
gencode.v35 = import("/home/public/reference/gtf/human/gencode.v35.annotation.gtf")
index = which(gencode.v35$type == "transcript") # extract transcript index
# extract transcript information
tmp=data.frame(Transcript_Id = gencode.v35$transcript_id[index],
               Gene_Type=gencode.v35$gene_type[index])

merged.info=left_join(merged.info,tmp,by=c("Transcript_Id"))

# set gene types
match=data.frame(Gene_Type = c("IG_C_gene","IG_D_gene","IG_J_pseudogene","IG_V_gene","lncRNA","misc_RNA","Mt_tRNA","processed_pseudogene",
                               "rRNA","scaRNA","snoRNA","sRNA","TR_C_gene","TR_J_gene","TR_V_gene","transcribed_processed_pseudogene",
                               "transcribed_unprocessed_pseudogene","Unknown","vault_RNA","IG_C_pseudogene","IG_J_gene","IG_pseudogene",
                               "IG_V_pseudogene","miRNA","Mt_rRNA","polymorphic_pseudogene","protein_coding","ribozyme","rRNA_pseudogene",
                               "scRNA","snRNA","TEC","TR_D_gene","TR_J_pseudogene","transcribed_unitary_pseudogene","translated_processed_pseudogene",
                               "unitary_pseudogene","unprocessed_pseudogene","pseudogene","TR_V_pseudogene","translated_unprocessed_pseudogene"
),
Classified_type=c("protein_coding","protein_coding","pseudogene","protein_coding","lncRNA","ncRNA","ncRNA","pseudogene",
                  "ncRNA","ncRNA","ncRNA","ncRNA","protein_coding","protein_coding","protein_coding","pseudogene",
                  "pseudogene","unknown","ncRNA","pseudogene","protein_coding","pseudogene",
                  "pseudogene","ncRNA","ncRNA","pseudogene","protein_coding","ncRNA","pseudogene",
                  "ncRNA","ncRNA","other","protein_coding","pseudogene","pseudogene","pseudogene",
                  "pseudogene","pseudogene","pseudogene","pseudogene","pseudogene"
))

merged.info=left_join(merged.info,match,by="Gene_Type")

# classify new genes
merged.info$Classified_type <- factor(merged.info$Classified_type,levels = c("protein_coding",
                                                                             "lncRNA",
                                                                             "ncRNA",
                                                                             "pseudogene",
                                                                             "other"))
merged.info <- merged.info[order(merged.info$Classified_type),]
merged.info$Classified_type <- as.character(merged.info$Classified_type)

merged_gene_type <- merged.info %>% 
  group_by(Gene_Id) %>% 
  summarise(Gene_Origin = .data$Classified_type[1])

merged.info <- left_join(merged.info, merged_gene_type, by = "Gene_Id")

merged.info_ENSG <- merged.info[grepl("ENSG", merged.info$Gene_Id),]
merged.info_ENST <- merged.info[grepl("ENST", merged.info$Transcript_Id),]
merged.info$Gene_Origin[is.na(merged.info$Gene_Origin)] = "intergenic"

merged.info$Transcript_Origin <- merged.info$Classified_type
merged.info$Transcript_Origin[is.na(merged.info$Transcript_Origin)] <- merged.info$Gene_Origin[is.na(merged.info$Transcript_Origin)]

colnames(merged.info)[5:6] <- c("Gene_Subype", "Gene_Type")

write.table(merged.info,"transcript_gene_origin.txt",sep="\t",row.names = F,quote=F)
