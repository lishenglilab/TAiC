# R version 4.0.2 (2020-06-22)
# tidyverse version=1.3.1
# rtracklayer version=1.5.0

library(tidyverse)
library(rtracklayer)

# import data
ccle <- import("/home/weihu/RBP_transcript_drug_project/data/gtf/CCLE_1017_merged_filtered_0.1.gtf")
genecode_v35 <- import("/home/public/reference/gtf/human/gencode.v35.annotation.gtf")

# extract novel and noncoding transcripts
noncoding_id <- unique(genecode_v35$transcript_id[genecode_v35$transcript_type != "protein_coding"])
novel_id <- unique(grep("MSTRG", ccle$transcript_id, value = T))
ccle_novel_noncoding <- ccle[ccle$transcript_id %in% c(noncoding_id, novel_id)]

# export
setwd("/home/weihu/RBP_transcript_drug_project/data/gtf")
export(ccle_novel_noncoding, "CCLE_1017_merged_novel_noncoding.gtf", format = "gtf")
