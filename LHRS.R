#-------------------------------------------------------------------------------
# IHS calculation

library(limma)
library(dplyr)
library(tibble)
library(tidyverse)
library(tidyr)
library(nlme)
library(ConsensusClusterPlus)
library(reshape2)
library(survival)
library(future.apply)
library(ggplot2)
library(ggpubr)
library(forcats)
library(ggrepel)
library(timeROC)
library(randomForestSRC)
library(stringr)
library(doParallel)
library(foreach)


raw_data <- read.table("logTPM_RJ_v5.txt", sep="\t",
                       header=T, check.names=F, row.names=1)
raw_data <- raw_data[, -which(str_detect(colnames(raw_data), "N"))]

# extract PCGs
load("gtf_v22.Rdata")
gene <- intersect(rownames(raw_data), gene_info$SYMBOL)
raw_data <- raw_data[gene,]

# exclude low expressed genes
raw_data <- 2^raw_data-1
raw_data <- raw_data[(rowSums(raw_data>1) >= 0.3*ncol(raw_data)),] %>% as.matrix()
raw_data <- log2(raw_data+1)

sampleN <- 14


# Intra-tumor variability score (ITVS)

data <- raw_data
sampleid <- data.frame(colnames(data))
names(sampleid)[1] <- "id"
sampleid <- separate(data = sampleid, col = id,
                     into = c("score1", "score2"), sep = "_")[,1]
data <- data.frame(id=sampleid, scale(t(data)))

var_result=data.frame()
for(i in colnames(data[,2:ncol(data)])){
  
  a <- data[i]
  a <- unlist(a)
  iris.lme <- lme(a~1, random=~1|id, data=data)
  var <- data.frame(as.character(VarCorr(iris.lme)))[,1] %>% as.numeric()
  
  var_sum <- data.frame(id=i,
                        between_group=var[1],
                        within_group=var[2],
                        std1=var[3],
                        std2=var[4],
                        all=var[1]+var[2],
                        WT=var[2]/(var[1]+var[2]))
  
  var_result <- rbind(var_result, var_sum)
  
}

var_result$id <- gsub('\\.', '-', var_result$id)


# Clustering concordance score (CCS)

data <- scale(t(raw_data))

integ_result <- data.frame()
for(gene in colnames(data)){
  data_cluster <- t(data[,gene])
  
  PGOR <- data.frame()
  for (i in 2:sampleN) {
    
    hc <- hclust(dist(t(data_cluster)), method = "ward.D2")
    cluster <- cutree(hc, k=i) %>% data.frame()
    names(cluster)[1]<- "cluster"
    
    cluster$type <- rownames(cluster)
    cluster <- separate(cluster,type,into = c("type","c"),sep = "_")
    cluster <- cluster[,-3]
    
    function1 <- function(a){
      ifelse(length(unique(a))==1,1,0)
    }
    
    aggdata<-aggregate(cluster[,1],by=list(cluster$type),FUN=function1)
    names(aggdata)[1:2] <- c("id","consistency")
    result_PGOR <- length(which(aggdata$consistency == 1))/length(aggdata$id) %>%
      data.frame()
    PGOR <- rbind(PGOR,result_PGOR)
  }
  
  names(PGOR)[1] <- "proportion"
  row <- data.frame(c(1))
  names(row)[1] <- "proportion"
  PGOR <- rbind(row, PGOR)
  
  x <- as.numeric(row.names(PGOR))
  y <- PGOR$proportion
  
  f <- approxfun(x,y)
  integral <- integrate(f, x[1], x[length(x)])$value
  integral <- data.frame(integral)
  integral$id <- gene
  integ_result <- rbind(integ_result, integral)
}

integ_result <- data.frame(id=integ_result$id,
                           raw_integral=integ_result$integral,
                           norm_integral=1-(integ_result$integral/(sampleN-1)))

all_result <- merge(var_result, integ_result, by="id")
all_result <- data.frame(id=all_result$id, WT=all_result$WT,
                         norm_integral=all_result$norm_integral,
                         ITH_score=sqrt((all_result$WT)*(all_result$norm_integral)))

write.table(all_result, "ITH_score_RJ_new.txt",
            sep="\t", quote=F, row.names=F, col.names=T)


#-------------------------------------------------------------------------------
# filtering step 1

HCC_cohorts <- readRDS("HCC_cohorts_survival.rds")
ITH_score <- read.table("ITH_score_RJ_new.txt", sep="\t", header=T,
                        check.names=F, quote="", na.strings = "")


allcox <- lapply(names(HCC_cohorts), function(i){
  
  cox_input <- HCC_cohorts[[i]]
  cox_input <- cox_input[,-1]
  
  coxf<-function(x){
    mycox <- coxph(Surv(time,status) ~ cox_input[,x], data = cox_input)
  }
  
  plan(multisession, workers = 12)
  
  cox_data <- future_lapply(colnames(cox_input[,3:ncol(cox_input)]), function(a){
    mycox <- coxf(a)
    coxResult <- summary(mycox)
    coxR <- cbind(id = a,HR = coxResult$coefficients[,"exp(coef)"],
                  P = coxResult$coefficients[,"Pr(>|z|)"], 
                  lower_95CI = coxResult$conf.int[,3],
                  upper_95CI = coxResult$conf.int[,4])
  })
  
  plan(sequential)
  
  cox_data <- do.call(rbind, cox_data)
  cox_data <- na.omit(cox_data)
  head(cox_data)
  
  cox_data <- data.frame(cox_data)
  cox_data <- data.frame(cohort=i, id=cox_data$id,
                         pvalue=as.numeric(cox_data$P), HR=cox_data$HR, 
                         FDR=p.adjust(as.numeric(cox_data$P), method="fdr"))
  cox_data
})


allcox <- do.call(rbind, allcox)
allcox$HR <- as.numeric(allcox$HR)


allgene <- ITH_score$id[ITH_score$ITH_score < 0.25]
allcox <- allcox[allcox$id %in% allgene,]

allHR <- data.frame(cohort=allcox$cohort, id=allcox$id, HR=allcox$HR)
allHR <- acast(allHR, id ~ cohort, value.var = "HR", fun.aggregate=mean)
allHR <- allHR[,c("TCGA","CHCC","GSE14520")]

allHR[allHR<1] <- -1
allHR[allHR>1] <- 1
allHR[is.na(allHR)] <- 0
allHR <- data.frame(allHR)
allHR$sum <- rowSums(allHR)

LHR <- rownames(allHR)[allHR$sum==-3]
HHR <- rownames(allHR)[allHR$sum==3]

subcox <- allcox[allcox$pvalue < 0.05,]
subcox_L <- subcox[subcox$id %in% LHR,]
siggene_L <- c(subcox_L$id[subcox_L$cohort=="TCGA"],
               subcox_L$id[subcox_L$cohort=="CHCC"],
               subcox_L$id[subcox_L$cohort=="GSE14520"])
number_L <- data.frame(table(siggene_L))
number_L <- number_L[number_L$Freq>1,]

subcox_H <- subcox[subcox$id %in% HHR,]
siggene_H <- c(subcox_H$id[subcox_H$cohort=="TCGA"],
               subcox_H$id[subcox_H$cohort=="CHCC"],
               subcox_H$id[subcox_H$cohort=="GSE14520"])
number_H <- data.frame(table(siggene_H))
number_H <- number_H[number_H$Freq>1,]

number_L$type <- -1; names(number_L)[1] <- "id"
number_H$type <- 1; names(number_H)[1] <- "id"
number <- rbind(number_L, number_H)

write.table(number, "gene_filter1.txt", sep="\t",
            quote=F, row.names=F, col.names=T)


#-------------------------------------------------------------------------------
# filtering step 2

filter1 <- read.table("gene_filter1.txt", sep="\t", header=T)
training <- c("TCGA", "CHCC", "GSE14520")


allcox <- lapply(training, function(i){
  
  all_cox_input <- HCC_cohorts[[i]]
  common <- intersect(colnames(all_cox_input), filter1$id)
  all_cox_input <- all_cox_input[,c("id", "status", "time", common)]
  
  plan(multisession, workers = 15)
  
  random_cox <- future_lapply(1:1000, function(j){
    random_id <- sample(all_cox_input$id, size=0.7*nrow(all_cox_input))
    cox_input <- all_cox_input[all_cox_input$id %in% random_id,]
    cox_input <- cox_input[,-1]
    
    coxf<-function(x){
      mycox <- coxph(Surv(time,status) ~ cox_input[,x], data = cox_input)
    }
    
    cox_data <- lapply(colnames(cox_input[,3:ncol(cox_input)]), function(a){
      mycox <- coxf(a)
      coxResult <- summary(mycox)
      coxR <- cbind(id = a,HR = coxResult$coefficients[,"exp(coef)"],
                    P = coxResult$coefficients[,"Pr(>|z|)"], 
                    lower_95CI = coxResult$conf.int[,3],
                    upper_95CI = coxResult$conf.int[,4])
    })
    
    cox_data <- do.call(rbind, cox_data)
    cox_data <- na.omit(cox_data)
    head(cox_data)
    
    cox_data <- data.frame(cox_data)
    cox_data <- data.frame(random=j, id=cox_data$id,
                           pvalue=as.numeric(cox_data$P))
    cox_data <- cox_data[cox_data$pvalue<0.05,]
    cox_data
  })
  
  plan(sequential)
  
  random_cox <- do.call(rbind, random_cox)
  frequency <- data.frame(table(random_cox$id))
  frequency$cohort <- i
  frequency
})

allcox <- do.call(rbind, allcox)
result <- acast(allcox, Var1 ~ cohort,
                value.var = "Freq", fun.aggregate=mean)
result[is.na(result)] <- 0

result <- data.frame(result)
result <- data.frame(id=rownames(result),
                     result, sum=rowSums(result))

result <- result[order(result$sum, decreasing=T),]
result$id <- fct_inorder(result$id)

filter2 <- result[result$sum>2000,]
filter2 <- merge(filter2, filter1, by="id")
filter2 <- filter2[,c(1,5,7)]
write.table(filter2, "gene_filter2.txt", sep="\t",
            quote=F, row.names=F, col.names=T)


#-------------------------------------------------------------------------------
# MD-based RSF 

filter2 <- read.table("gene_filter2.txt", sep="\t", header=T)

#--------------------
train1 <- HCC_cohorts[["TCGA"]]
train2 <- HCC_cohorts[["CHCC"]]
train3 <- HCC_cohorts[["GSE14520"]]

common <- intersect(colnames(train1), filter2$id)
common <- intersect(colnames(train2), common)
common <- intersect(colnames(train3), common)

train1 <- train1[,c("status", "time", common)]
train2 <- train2[,c("status", "time", common)]
train3 <- train3[,c("status", "time", common)]

#--------------------
data <- train1

#--------------------
IHS1 <- read.table("ITH_score_RJ_new.txt",sep="\t",header=T)
ITH1 <- read.table("logTPM_RJ_v5.txt", sep="\t", header=T, row.names=1)

#-------------------------------------------------------------------------------
MultiNames <- colnames(data[,3:ncol(data)])
mtry <- sqrt(length(MultiNames))
res.rsf <- rfsrc(Surv(time, status) ~ ., data, proximity=T,
                 tree.err=T, mtry=mtry, forest=T, ntree=500,
                 splitrule="logrank", importance=TRUE)

cl<- makeCluster(detectCores()-1)      
registerDoParallel(cl)
getDoParWorkers()

all_result <- foreach(j=1:1000, 
                      .combine="rbind",
                      .packages = c("timeROC",
                                    "randomForestSRC",
                                    "Hmisc",
                                    "survival",
                                    "survminer",
                                    "tidyverse",
                                    "stringr"),
                      .verbose = T) %dopar% {
                        
                        vars <- var.select(object=res.rsf,
                                           cause=1, method="vh", 
                                           conservative=c("high"), 
                                           ntree=500, K=5, mtry=mtry,
                                           splitrule="logrank", nsplit=3, 
                                           refit=T, fast=T, nrep=10, nstep=5,
                                           prefit=list(action=T, ntree=500, nsplit=3),
                                           verbose=TRUE)
                        
                        signature <- vars$topvars
                        file <- paste("signature/A", j, ".txt", sep="")
                        sig_result <- data.frame(signature)
                        write.table(sig_result, file, sep="\t",
                                    quote=F, row.names=F, col.names=T)
                        
                        LHR <- filter2$id[filter2$type==-1]
                        Lgene <- intersect(LHR, signature)
                        
                        HHR <- filter2$id[filter2$type==1]
                        Hgene <- intersect(HHR, signature)
                        
                        
                        if (length(Lgene)!=0 & length(Hgene)!=0) {
                          
                          #------------------------------
                          # train1
                          train1_score <- data.frame(train1[,c(1:2)],
                                                     score=apply(data.frame(train1[,Hgene]),1,mean)-
                                                       apply(data.frame(train1[,Lgene]),1,mean))
                          ROC1 <- timeROC(T=train1_score$time,
                                          delta=train1_score$status,
                                          marker=train1_score$score,
                                          cause=1, weighting="marginal", 
                                          times=c(12,24,36,48), iid=TRUE)
                          AUC_train1 <- mean(ROC1[["AUC"]])
                          
                          # train2
                          train2_score <- data.frame(train2[,c(1:2)],
                                                     score=apply(data.frame(train2[,Hgene]),1,mean)-
                                                       apply(data.frame(train2[,Lgene]),1,mean))
                          ROC2 <- timeROC(T=train2_score$time,
                                          delta=train2_score$status,
                                          marker=train2_score$score,
                                          cause=1, weighting="marginal", 
                                          times=c(12,24,36,48), iid=TRUE)
                          AUC_train2 <- mean(ROC2[["AUC"]])
                          
                          # train3
                          train3_score <- data.frame(train3[,c(1:2)],
                                                     score=apply(data.frame(train3[,Hgene]),1,mean)-
                                                       apply(data.frame(train3[,Lgene]),1,mean))
                          ROC3 <- timeROC(T=train3_score$time,
                                          delta=train3_score$status,
                                          marker=train3_score$score,
                                          cause=1, weighting="marginal", 
                                          times=c(12,24,36,48), iid=TRUE)
                          AUC_train3 <- mean(ROC3[["AUC"]])
                          
                          #------------------------------
                          result <- data.frame(id=paste("A", j, sep=""),
                                               
                                               AUC_train1=AUC_train1,
                                               AUC_train2=AUC_train2,
                                               AUC_train3=AUC_train3,

                                               all_AUC_train_mean=mean(c(AUC_train1, AUC_train2, AUC_train3)),
                                               all_AUC_train_sd=-log2(sd(c(AUC_train1, AUC_train2, AUC_train3))))
                          
                          result
                        }
                      }

stopCluster(cl)

write.table(all_result, "all_sig_result.txt",
            sep="\t", quote=F, row.names=F, col.names=T)




