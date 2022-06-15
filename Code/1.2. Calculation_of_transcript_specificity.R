# R version 4.0.2 (2020-06-22)
# tidyverse version=1.3.1
## Calculate the expression specificity of transcripts across different cancer cell lineages or cancer types
library(tidyverse)
args=commandArgs(T)

# Import data
exp <- data.table::fread(paste0("transcript_", args[1], ".txt"), header = T, data.table = F, check.names = F)
meta <- data.table::fread("info.ccle-site-type_check_02.txt", header = T, data.table = F, check.names = F)

# Remove primary sites with less than 5 cancer cell lines
summary <- data.frame(table(meta$Site_Primary))
abandon <- as.character(summary$Var1[summary$Freq < 5])
exp <- exp[, !colnames(exp) %in% meta$CCLE_ID[meta$Site_Primary %in% abandon]]
meta <- meta[!meta$Site_Primary %in% abandon, ]

# Calculate average expression of each primary site
expl <- exp %>% 
  pivot_longer(cols = 2:ncol(exp),
               names_to = "Cell_Line",
               values_to = "TPM")
expl <- left_join(expl, meta[, c(1:2)], by = c("Cell_Line" = "CCLE_ID"))
expl_mean <- expl %>% 
  group_by(Transcript_Id, Site_Primary) %>% 
  summarise(Mean = mean(TPM))
  
# Retain transcripts with at least 0.1 average TPM in one primary site
expl_meanw <- expl_mean %>% 
  pivot_wider(names_from = Site_Primary, values_from = Mean )
expl_meanw <- expl_meanw[apply(expl_meanw[, 2:ncol(expl_meanw)], 1, function(x) sum(x >= 0.1) >= 1 ),]

# Calculate specificity score and expression ratio of each transcript
expl_mean <- expl_meanw %>%
  pivot_longer(cols = 2:ncol(expl_meanw),
               names_to = "Site_Primary",
               values_to = "Mean")

N <- length(unique(expl_mean$Site_Primary))
rt <- data.frame()

for(t in unique(expl_mean$Transcript_Id)){
  expl_mean_t <- expl_mean %>% 
    filter(Transcript_Id == t)
  sum_t <- sum(expl_mean_t$Mean)
  
  Pig_t <- expl_mean_t %>% 
    group_by(Site_Primary) %>% 
    summarise(Pig = Mean/sum_t) %>% 
    arrange(desc(Pig))
  
  ratio <- Pig_t$Pig[1]/Pig_t$Pig[2]
  
  Sg <- expl_mean_t %>% 
    group_by(Site_Primary) %>% 
    summarise(int=(Mean/sum_t)*log2(Mean/sum_t))
  Sg$int[Sg$int == "NaN"]=0
  Sgt <- log2(N) + sum(Sg$int)
  Sgt <- data.frame(Transcript_Id = t, Specificity_Score = Sgt, Ratio = ratio, Primary_site = Pig_t$Site_Primary[1])
  rt <- rbind(rt, Sgt)
}

write.table(rt, file = paste0("CCLE_transcript_primary_site_specificity_score_ratio_", args[1], ".txt"), sep = "\t", row.names = F, quote = F)
