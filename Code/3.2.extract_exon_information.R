# R version 4.0.2 (2020-06-22)
# tidyverse version=1.3.1
# rtracklayer version=1.5.0
library(tidyverse)
library(rtracklayer)

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
