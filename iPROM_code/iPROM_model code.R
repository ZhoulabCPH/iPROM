###  iPROM ###
library(xgboost)
library(Matrix)
library(verification)
library(pROC)
library(ggplot2)
library(ggpubr)
library(modEvA)
library(caret)
train <- read.csv("train_data.csv",stringsAsFactors = F,
                  row.names = 1,check.names = F)
test <- read.csv("test_data.csv",stringsAsFactors = F,
                 row.names = 1,check.names = F)
train <- train[,2:52]
test <- test[,2:26]
train1 <- apply(train,2,scale)
rownames(train1) <- rownames(train)
test1 <- apply(test,2,scale)
rownames(test1) <- rownames(test)
train <- train1
test <- test1
imm_label <- read.csv("Immune_consensusCluster_label.csv",
                      stringsAsFactors = F,
                      row.names = 1,check.names = F)
imm_label$mod_lable <- imm_label$imm_cluster
imm_label$mod_lable[which(imm_label$imm_cluster == "imm_hot")] <- 1
imm_label$mod_lable[which(imm_label$imm_cluster == "imm_cold")] <- 0
train_label <- imm_label[colnames(train),]
test_label <- imm_label[colnames(test),]
train <- as.matrix(t(train[c("P31995",
                             "P08637","P07858"),]))
test <- as.matrix(t(test[c("P31995",
                           "P08637","P07858"),]))
traindata1 <- data.matrix(train)  ## 
traindata2 <- Matrix(traindata1,sparse=T)  ## 
traindata3 <- factor(train_label$mod_lable,levels = c(0,1))   ### 
traindata4 <- list(data=traindata2,label=train_label$mod_lable)  ### candidate training data

dtrain <- xgb.DMatrix(data = traindata4$data, label = traindata4$label)
set.seed(3)
mxgb4m <- xgboost(data = dtrain,   
                  objective='binary:logistic',
                  nround= 300, 
                  nfold = 5,  
                  max_depth=2,  
                  subsample = 0.9,
                  colsample_bytree = 0.9,
                  eta=0.3,
                  eval_metric = "rmse")
## training all
train_pred <- predict(mxgb4m,train,type = "response")
train_pred_result <- data.frame(pred = train_pred,
                                name = rownames(train),
                                Cohort = rep("Train",length(train_pred)))
a <- roc(train_label$mod_lable,train_pred)
train.roc <- roc(train_label$mod_lable,train_pred, plot=TRUE, 
                 print.thres=TRUE, ci=TRUE,
                 print.auc=TRUE,legacy.axes = TRUE,col = "#1687A7")
train.roc <- ci.sp(train.roc, sensitivities=seq(0, 1, .01), boot.n=100)  ## 
plot(train.roc, type="shape", col="#8A9EB5",ci=TRUE)  ## 
AUC(obs=train_label$mod_lable,pred=train_pred,curve ="PR", simplif=TRUE, main = "PR curve")
actual <- as.numeric(train_label$mod_lable)
actual[actual == 0] <- 2
train_pred[train_pred > 0.5] <- 1
train_pred[train_pred <= 0.5] <- 2
confusionMatrix(factor(train_pred), factor(actual),
                mode = "everything",positive="1")
## testing cohort
test_pred <- predict(mxgb4m,test,type = "response")
test_pred_result <- data.frame(pred = test_pred,
                               name = rownames(test),
                               Cohort = rep("Test",length(test_pred)))
b <- roc(test_label$mod_lable,test_pred)
test.roc <- roc(test_label$mod_lable,test_pred, plot=TRUE, 
                print.thres=TRUE, ci=TRUE,
                print.auc=TRUE,legacy.axes = TRUE,col = "#1687A7")
test.roc <- ci.sp(test.roc, sensitivities=seq(0, 1, .01), boot.n=100)  ## 
plot(test.roc, type="shape", col="#8A9EB5",ci=TRUE)  ## 
AUC(obs=test_label$mod_lable,pred=test_pred,curve ="PR", 
    simplif=TRUE, main = "PR curve")
actual <- as.numeric(test_label$mod_lable)
actual[actual == 0] <- 2
test_pred[test_pred > 0.5] <- 1
test_pred[test_pred <= 0.5] <- 2
confusionMatrix(factor(test_pred), factor(actual),
                mode = "everything",
                positive="1")



###  George cohort ###
library(xgboost)
library(Matrix)
library(verification)
library(pROC)
library(ggplot2)
library(ggpubr)
library(modEvA)
library(caret)
train <- read.csv("train_data.csv",stringsAsFactors = F,
                  row.names = 1,check.names = F)
train <- train[,2:52]
train1 <- apply(train,2,scale)
rownames(train1) <- rownames(train)
train <- train1
imm_label <- read.csv("Immune_consensusCluster_label.csv",
                      stringsAsFactors = F,
                      row.names = 1,check.names = F)
imm_label$mod_lable <- imm_label$imm_cluster
imm_label$mod_lable[which(imm_label$imm_cluster == "imm_hot")] <- 1
imm_label$mod_lable[which(imm_label$imm_cluster == "imm_cold")] <- 0
train_label <- imm_label[colnames(train),]
train <- as.matrix(t(train[c("P31995",
                             "P08637","P07858"),]))
traindata1 <- data.matrix(train)  ## 
traindata2 <- Matrix(traindata1,sparse=T)  ## 
traindata3 <- factor(train_label$mod_lable,levels = c(0,1))   ### 
traindata4 <- list(data=traindata2,label=train_label$mod_lable)  ### candidate training data
dtrain <- xgb.DMatrix(data = traindata4$data, label = traindata4$label)
set.seed(3)
mxgb4m <- xgboost(data = dtrain,   
                  objective='binary:logistic',
                  nround= 300, 
                  nfold = 5,  
                  max_depth=2,  
                  subsample = 0.9,
                  colsample_bytree = 0.9,
                  eta=0.3,
                  eval_metric = "rmse")
## RNA-seq
data_lcnec <- read.csv("LCNEC_count_RSME.csv",
                       stringsAsFactors = F,row.names = 1,
                       header = T,check.names = F)
data_lcnec <- log2(data_lcnec+1)
data_sclc <- read.table("81tumor_18816gene_data_fpkm.txt",
                        stringsAsFactors = F,row.names = 1,
                        header = T,check.names = F)
data_sclc <- log2(data_sclc+1)
label <- read.csv("RNAseq_consensusCluster_label.csv",row.names = 1,
                  stringsAsFactors = F,check.names = F)
label$lab <- 1:dim(label)[1]
label$lab[which(label$label == "C1")] <- 1
label$lab[which(label$label == "C2")] <- 0
label_lcnec <- label[which(label$type == "LCNEC"),]
label_sclc <- label[which(label$type == "SCLC"),]
data_lcnec <- data_lcnec[c("CTSB","FCGR2C","FCGR3A"),]
linshi <- apply(data_lcnec,1,scale)
rownames(linshi) <- colnames(data_lcnec)
data_lcnec1 <- linshi
colnames(data_lcnec1) <- c("P31995","P08637","P07858")
lcnec_pred <- predict(mxgb4m,data_lcnec1,type = "response")
lcnec_pred_result <- data.frame(pred = lcnec_pred,
                                name = rownames(data_lcnec1),
                                Cohort = rep("LCNEC",length(lcnec_pred)))
label_lcnec <- label_lcnec[rownames(data_lcnec1),]
data_sclc <- data_sclc[c("CTSB","FCGR2C","FCGR3A"),]
linshi <- apply(data_sclc,1,scale)
rownames(linshi) <- colnames(data_sclc)
data_sclc1 <- linshi
colnames(data_sclc1) <- c("P31995","P08637","P07858")
sclc_pred <- predict(mxgb4m,data_sclc1,type = "response")
sclc_pred_result <- data.frame(pred = sclc_pred,
                               name = rownames(data_sclc1),
                               Cohort = rep("SCLC",length(sclc_pred)))
label_sclc <- label_sclc[rownames(data_sclc1),]
###  Lu-NECs
lcnec_pred <- predict(mxgb4m,data_lcnec1,type = "response")
sclc_pred <- predict(mxgb4m,data_sclc1,type = "response")
all_lab <- c(label_lcnec$lab,label_sclc$lab)
all_pred <- c(lcnec_pred,sclc_pred)
all_ctsb <- c(data_lcnec1[,1],data_sclc1[,1])

all.roc <- roc(all_lab,all_pred, plot=TRUE, 
               print.thres=TRUE, ci=TRUE,
               print.auc=TRUE,legacy.axes = TRUE,col = "#1687A7")
all.roc <- ci.sp(all.roc, sensitivities=seq(0, 1, .01), boot.n=100)  ## 
plot(all.roc, type="shape", col="#8A9EB5",ci=TRUE)  ## 
AUC(obs=all_lab,pred=all_pred,curve ="PR", simplif=TRUE, main = "PR curve")
actual <- as.numeric(all_lab)
actual[actual == 0] <- 2
all_pred[all_pred > 0.5] <- 1
all_pred[all_pred <= 0.5] <- 2
confusionMatrix(factor(all_pred), factor(actual),
                mode = "everything",positive="1")
data_mat <- data.frame(score = c(lcnec_pred,sclc_pred),
                       lab = all_lab)
red <- "#D94E48";
blue <- "#A6A6A6";
white <- rgb(255,255,255,maxColorValue = 255)
data_mat <- data_mat[rev(order(data_mat$score)),]
pheatmap(data_mat,fontsize=6,cellwidth = 6,gaps_row = 55,
         color  = colorRampPalette(c(blue,white,red))(100),
         legend_breaks=seq(0,1,by=0.1),
         clustering_method = "ward.D2",
         border_color = "grey60",
         cluster_cols = F, cluster_rows = F,
         show_rownames = T, show_colnames = T)
rownames(data_mat) <- c(rownames(data_lcnec1),rownames(data_sclc1))
#write.csv(data_mat,"George_iPROM_result.csv",quote = F)



# Roper cohort
library(xgboost)
library(Matrix)
library(verification)
library(pROC)
library(ggplot2)
library(ggpubr)
library(modEvA)
library(caret)
train <- read.csv("train_data.csv",stringsAsFactors = F,
                  row.names = 1,check.names = F)
train <- train[,2:52]
train1 <- apply(train,2,scale)
rownames(train1) <- rownames(train)
train <- train1
imm_label <- read.csv("Immune_consensusCluster_label.csv",
                      stringsAsFactors = F,
                      row.names = 1,check.names = F)
imm_label$mod_lable <- imm_label$imm_cluster
imm_label$mod_lable[which(imm_label$imm_cluster == "imm_hot")] <- 1
imm_label$mod_lable[which(imm_label$imm_cluster == "imm_cold")] <- 0
train_label <- imm_label[colnames(train),]
train <- as.matrix(t(train[c("P31995","P08637","P07858"),]))
traindata1 <- data.matrix(train)  ## 
traindata2 <- Matrix(traindata1,sparse=T)  ## 
traindata3 <- factor(train_label$mod_lable,levels = c(0,1))   ### 
traindata4 <- list(data=traindata2,label=train_label$mod_lable)  ### candidate training data
dtrain <- xgb.DMatrix(data = traindata4$data, label = traindata4$label)
set.seed(3)
mxgb4m <- xgboost(data = dtrain,   
                  objective='binary:logistic',
                  nround= 300, nfold = 5,  
                  max_depth=2,  subsample = 0.9,
                  colsample_bytree = 0.9,eta=0.3,
                  eval_metric = "rmse")
setwd("D:\\北京多种肺癌免疫蛋白亚型\\公共数据\\免疫治疗_SCLC")
roper <- read.csv("data_rpkm.csv",stringsAsFactors = F,
                  row.names = 1,check.names = F)
roper_t <- apply(roper[c("FCGR2C","FCGR3A","CTSB"),], 1, scale)
roper_t <- t(roper_t)
colnames(roper_t) <- colnames(roper)
# P31995 : FCGR2C
# P08637 : FCGR3A
# P07858 : CTSB
roper <- t(roper_t)
colnames(roper) <- c("P31995","P08637","P07858")
roper <- as.matrix(roper)
tumor_pred <- predict(mxgb4m,roper,type = "response")
tumor_pred <- data.frame(name = rownames(roper),
                         iPROM_score = tumor_pred)
tumor_pred <- cbind(tumor_pred,roper)
#write.csv(tumor_pred,"Roper_iPROM_result.csv",quote = F)



## IHC cohort
library(xgboost)
library(Matrix)
library(verification)
library(pROC)
library(ggplot2)
library(ggpubr)
library(modEvA)
library(caret)
train <- read.csv("train_data.csv",stringsAsFactors = F,
                  row.names = 1,check.names = F)
train <- train[,2:52]
train1 <- apply(train,2,scale)
rownames(train1) <- rownames(train)
train <- train1
imm_label <- read.csv("Immune_consensusCluster_label.csv",
                      stringsAsFactors = F,
                      row.names = 1,check.names = F)
imm_label$mod_lable <- imm_label$imm_cluster
imm_label$mod_lable[which(imm_label$imm_cluster == "imm_hot")] <- 1
imm_label$mod_lable[which(imm_label$imm_cluster == "imm_cold")] <- 0
train_label <- imm_label[colnames(train),]
train <- as.matrix(t(train[c("P31995","P08637","P07858"),]))
traindata1 <- data.matrix(train)  ## 
traindata2 <- Matrix(traindata1,sparse=T)  ## 
traindata3 <- factor(train_label$mod_lable,levels = c(0,1))   ### 
traindata4 <- list(data=traindata2,label=train_label$mod_lable)  ### candidate training data
dtrain <- xgb.DMatrix(data = traindata4$data, label = traindata4$label)
set.seed(3)
mxgb4m <- xgboost(data = dtrain,   
                  objective='binary:logistic',
                  nround= 300, nfold = 5,  
                  max_depth=2,  subsample = 0.9,
                  colsample_bytree = 0.9,eta=0.3,
                  eval_metric = "rmse")
ihc <- read.csv("IHC_iPROM_result.csv",stringsAsFactors = F,
                row.names = 1,check.names = F)
ihc1 <- apply(ihc[,1:3], 2, scale)
rownames(ihc1) <- rownames(ihc)
# P31995 : FCGR2C
# P08637 : FCGR3A
# P07858 : CTSB
ihc_tumor <- ihc1
colnames(ihc_tumor) <- c("P31995","P08637","P07858")
ihc_tumor <- as.matrix(ihc_tumor)
iPROM_score <- predict(mxgb4m,ihc_tumor,type = "response")

pre_result <- data.frame(iPROM_score = iPROM_score)
rownames(pre_result) <- rownames(ihc)
pre_result <- cbind(ihc,pre_result)
#write.csv(pre_result,"IHC_iPROM_result.csv",quote = F)






































