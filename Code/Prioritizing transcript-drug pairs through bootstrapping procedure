library(tidyverse)
library(caret)
library(glmnet)
args=commandArgs(T)

# import expression data
exp <- data.table::fread("CCLE_transcript_0.1.txt", header = T, data.table = F, check.names = F)
rownames(exp) <- exp[, 1]
exp <- exp[, -1]
colnames(exp) <- gsub("_.*", "", colnames(exp))

# import drug sensitivity data
auc <- read.table(paste0("CTRP_AUC_",args[1],".txt"), header = T, sep = "\t", check.names = F, row.names = 1)

# set identical id of cell lines
com_id <- intersect(colnames(exp),colnames(auc))
exp1 <- exp[, com_id]
auc1 <- auc[, com_id]

# build result matrix
feature_score <- as.data.frame(matrix(0,ncol=1,nrow=(nrow(exp1)+1)))
rownames(feature_score) <- c("Intercept",rownames(exp1))

# Associating transcript expression and drug response by Elastic Net Regression
for (i in rownames(auc1)){
  print(i)
  
  # import spearman correlation results
  genecor <- data.table::fread(paste0("Cor_transcript-",i,".txt"), header =T , data.table = F, check.names = F)
  
  # remove transcript-drug pairs with correlation index lower than 0.2
  genecor_cor <- subset(genecor, abs(cor) > 0.2 & FDR < 0.05)
  rm(genecor)
  
  if(nrow(genecor_cor) >= 2){
    
    # process data for drug i
    index <- rownames(exp1) %in% genecor_cor$Gene_2
    expi <- exp1[index,]
    data <- as.data.frame(rbind(expi,auc1))
    tdata <- as.data.frame(t(data)) 
    X <- tdata[, rownames(expi)]
    Y <- as.numeric(tdata[, i]) 
    rt <- as.data.frame(cbind(X, Y))
    rt <- subset(rt, Y != "NA")
    
    set.seed(123)
    seeds <- vector(mode = "list", length = 51)
    for(a in 1:50) seeds[[a]] <- sample.int(1000, 25)
    seeds[[51]] <- sample.int(1000, 1)
    
    # Model building: Elastic Net Regression 
    control <- trainControl(method = "repeatedcv", 
                            number = 10, 
                            repeats = 5,
                            search = "random", 
                            verboseIter = TRUE,
                            seeds = seeds) 
    
    # Training ELastic Net Regression model 
    elastic_model <- train( Y ~ .,
                            data = rt, 
                            method = "glmnet", 
                            preProcess = c("center", "scale"), 
                            tuneLength = 25,
                            trControl = control) 
    
    # Extract best tuning parameter
    alpha <- as.numeric(elastic_model$bestTune)[1]
    lambda <- as.numeric(elastic_model$bestTune)[2]
    best_parameter <- data.frame(i,alpha,lambda)
    
    # Bootstrap coef for Elastic Net
    bootstrap_function  <- function(model_data, ndraws) {
      coeff_mtx <- matrix(0, 
                          nrow = ndraws, 
                          ncol = ncol(model_data))
      
      for (i in 1:ndraws) {
        # Bootstrap the data by sampling
        bootstrap_ids <- sample(seq(nrow(model_data)),
                                  nrow(model_data),
                                  replace = TRUE)
        d <- model_data[bootstrap_ids,]
        
        # Estimate the model and save the coefficients
        fit <- glmnet(scale(d[,1:ncol(d)-1]), as.numeric(d[,ncol(d)]), alpha = best_parameter$alpha, lambda = best_parameter$lambda)
        coeff_mtx[i,]  <- as.matrix(coef(fit))
      }
      colnames(coeff_mtx) <- c("Intercept", colnames(rt[,1:ncol(rt)-1]))
      return(coeff_mtx)
    }
    
    # Run bootstrap function
    set.seed(123)
    bootstrap_coef <- bootstrap_function(rt, 1000)
    
    # Calculate prediction score
    score <- apply(bootstrap_coef, 2, function(x) {
      if (sum(x>0) > sum(x<0)){
        sum(x>0)/1000 - sum(x<0)/1000
      } else { 
        sum(x<0)/1000 - sum(x>0)/1000}
    })
    
    # generate results
    feature <- data.frame(score)
    colnames(feature) <- i
    row <- nrow(feature_score) - nrow(feature)
    m <- as.data.frame(matrix(0, ncol = 1, nrow = row))
    colnames(m) <- i
    index <- rownames(exp1) %in% rownames(feature)
    rownames(m) <- rownames(exp1)[!index]
    n <- rbind(feature, m)
    n <- n[match(rownames(feature_score), rownames(n)), ]
    feature_score <- cbind(feature_score,n)
    colnames(feature_score)[ncol(feature_score)] <- i
    
  }else
  {
    print(paste0(i," failed"))
  }
}

feature_score <- feature_score[,-1]
write.table(feature_score, file = paste0("CTRP_transcript_score_",args[1],".txt"), sep = "\t", quote = F)
