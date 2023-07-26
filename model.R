# 06/25/2023
# For publication
# 

setwd("C:/Users/sunid/Box/DuttaLab/deepLearning/")
library(tfdatasets)
library(caret)
library(PepTools) 
library(purrr)
library(readxl)
# library(dplyr)
library(gt)
library(gapminder)

# Read in table of high confidence miRNAs

mir_lvl <- read.table("C:/Users/sunid/Downloads/high_conf_mature_miRs.txt", header = F)


# Name precursor miRNA column

mir_lvl$precursor <- gsub("-5p","", mir_lvl$V1)
mir_lvl$precursor <- gsub("-3p","", mir_lvl$precursor)

# further editing miRNA name to stay consistent across dataframes

mir_lvl$miRNA <- toupper(gsub("-","",mir_lvl$V1))
mir_lvl$miRNA <- toupper(gsub("3P","_3P",mir_lvl$miRNA))
mir_lvl$miRNA <- toupper(gsub("5P","_5P",mir_lvl$miRNA))
mir_lvl$miRNA <- (gsub("HSA","",mir_lvl$miRNA))
head(mir_lvl)

# Read in effect size information from ENCODE shRNA sequencing

all_mirEdited <- readRDS("100821all_mirEdited.rds")

# Change column header names

names(all_mirEdited) <- gsub("Var2", "miRNA", names(all_mirEdited))
names(all_mirEdited) <- gsub("Var1", "query", names(all_mirEdited))

# Merge miRNA names dataframe and ES dataframe so that only high confidence miRNAs are considered 

all_mirEdited <- all_mirEdited[all_mirEdited$miRNA %in% mir_lvl$miRNA,]

# Read in list of gene-miRNA binders

mir_gt <- read.table("C:/Users/sunid/Downloads/matureOnly-final-ENCODE-miRs-HepG2.txt", header = F)

# Edit column names and miRNA naming format

colnames(mir_gt) <- c("Entrez","Gene")
mir_gt$miRNA <- gsub("hsa-","",mir_gt$Gene)
mir_gt$miRNA <- toupper(gsub("-","",mir_gt$miRNA))
mir_gt$miRNA <- toupper(gsub("3P","_3P",mir_gt$miRNA))
mir_gt$miRNA <- toupper(gsub("5P","_5P",mir_gt$miRNA))

# Renaming dataframe

all_mir_PepSeq <- all_mirEdited

# Combining ES and binding dataframes  

all_mir_PepSeq$effector <- paste(all_mir_PepSeq$query, all_mir_PepSeq$miRNA) %in% paste(mir_gt$Entrez, mir_gt$miRNA)

head(all_mir_PepSeq)
dim((all_mir_PepSeq))

all_mir_PepSeq$effector <- gsub("FALSE","notBinder",all_mir_PepSeq$effector)
all_mir_PepSeq$effector <- gsub("TRUE","Binder",all_mir_PepSeq$effector)
head(all_mir_PepSeq)
dim(all_mir_PepSeq)

# Defining binders and effectors

all_mir_PepSeq <- all_mir_PepSeq %>%
  mutate(effector = case_when(effector == "Binder" & padj < 0.05 ~ 'BinderEffector',
                              effector == "Binder" & padj > 0.05 ~ 'BinderNotEffector',
                              effector == "notBinder" & padj < 0.05 ~ 'notBinderEffector',
                              effector == "notBinder" & padj > 0.05 ~ 'notBinderNotEffector',
                              TRUE ~ effector))

# Reading in precursor miRNA sequences & editing miRNA naming format

miRseqs <- read.csv("edited-nkim_precursor_miRNAseqs.csv")[,c(1,2)]
colnames(miRseqs) <- c("miRNA","seq")
head(miRseqs)

miRseqs$miRNA <- gsub("-5p","_5P",miRseqs$miRNA)
miRseqs$miRNA <- gsub("-3p","_3P",miRseqs$miRNA)
miRseqs$miRNA <- gsub("hsa-","",miRseqs$miRNA)
miRseqs$miRNA <- toupper(gsub("-","",miRseqs$miRNA))
miRseqs$miRNA <- toupper(miRseqs$miRNA)
head(miRseqs)

all_mir_PepSeq$miRNAmat <- all_mir_PepSeq$miRNA
all_mir_PepSeq$miRNA <- gsub("_5P","",all_mir_PepSeq$miRNA)
all_mir_PepSeq$miRNA <- gsub("_3P","",all_mir_PepSeq$miRNA)

# Merging binder-effector dataframe and pre-miRNA sequences

all_mir_PepSeq <- unique(merge(all_mir_PepSeq, miRseqs, by = "miRNA"))

# Splitting sequences so that each nucleotide is in its own column

ES_mirSeqSplit <- strsplit((all_mir_PepSeq$seq), "")
t2 <- data.frame(lapply((ES_mirSeqSplit), type.convert), stringsAsFactors=FALSE)

# Numerical encoding for miRNA sequences

ES_mir1hot <- map_df(t2, ~ gsub("A", "1", .x))
ES_mir1hot <- map_df(ES_mir1hot, ~ gsub("T", "2", .x))
ES_mir1hot <- map_df(ES_mir1hot, ~ gsub("C", "3", .x))
ES_mir1hot <- map_df(ES_mir1hot, ~ gsub("G", "4", .x))

ES_mir1hot <- data.matrix(t(ES_mir1hot))

ES_mir1hot <- data.frame(matrix(as.numeric(ES_mir1hot),    # Convert to numeric matrix
                                ncol = ncol(ES_mir1hot)),data.frame(all_mir_PepSeq$query))[,-126]
head(ES_mir1hot)


all_mir_PepSeq <- unique(cbind(all_mir_PepSeq,ES_mir1hot))

head(all_mir_PepSeq)



# Obtaining the amino acid sequences for each gene

library(biomaRt)
mart <- useEnsembl("ensembl", dataset="hsapiens_gene_ensembl")
seq <-  getSequence(id = (as.character(all_mir_PepSeq$query)),
                    type = "hgnc_symbol",
                    seqType = "peptide",
                    mart = mart)
seqNoNA <- subset(seq, peptide != "Sequence unavailable")
library(data.table)
seqNoNA <- seqNoNA[!seqNoNA$peptide %like% "\\*",]
seqNoNA <- seqNoNA[!seqNoNA$peptide %like% "X",]

library(Peptides)
ai <- aIndex(seqNoNA$peptide)
bm <- boman(seqNoNA$peptide)
char <- charge(seqNoNA$peptide)
hydrophoIn <- hydrophobicity(seqNoNA$peptide)
instab <- instaIndex(seqNoNA$peptide)
molwgt <- mw(seqNoNA$peptide)
isoPt <- pI(seqNoNA$peptide)

library(plyr)
seqNoNA$crucianiProp <- ldply(seqNoNA$peptide,function(x){
  matrix(unlist(crucianiProperties(x)), ncol = 3, nrow = 1)
})


seqNoNA$Kidera <- ldply(seqNoNA$peptide,function(x){
  matrix(unlist(kideraFactors(x)), ncol = 10, nrow = 1)
})


seqNoNA$FP <- ldply(seqNoNA$peptide,function(x){
  matrix(unlist(protFP(x)), ncol = 8, nrow = 1)
})
seqNoNA <- seqNoNA %>% dplyr::select(-peptide) 
head(seqNoNA)

seqNoNA <- data.frame(seqNoNA,ai,bm,char,hydrophoIn,instab,molwgt,isoPt)
DF <- as.data.table(seqNoNA)
head(DF)
sapply(DF, class)

pc_prot <- DF %>%
  group_by(hgnc_symbol) %>%
  summarize_all(mean)
pc_prot



# Identifying potentially important domains/Go functions for each protein

ensembl <- useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl")

# reg of translation
transl <- getBM(attributes = c('hgnc_symbol'), 
                filters = 'go', 
                values = 'GO:0006417', 
                mart = ensembl)

all_mir_PepSeq$Translationregulation <- with(all_mir_PepSeq, ifelse(all_mir_PepSeq$query %in% transl$hgnc_symbol, 1, 0))


# response to virus GO term
virus <- getBM(attributes = c('hgnc_symbol'), 
               filters = 'go', 
               values = 'GO:0009615', 
               mart = ensembl)

all_mir_PepSeq$ViralRNAregulation <- with(all_mir_PepSeq, ifelse(all_mir_PepSeq$query %in% virus$hgnc_symbol, 1, 0))


# mirna metabolic process GO term
mirmetproc <- getBM(attributes = c('hgnc_symbol'), 
                    filters = 'go', 
                    values = 'GO:0010586', 
                    mart = ensembl)

all_mir_PepSeq$microRNAprocessing <- with(all_mir_PepSeq, ifelse(all_mir_PepSeq$query %in% mirmetproc$hgnc_symbol, 1, 0))


# dsRNA binding GO term
dsRNA <- getBM(attributes = c('hgnc_symbol'), 
               filters = 'go', 
               values = 'GO:0003725', 
               mart = ensembl)

all_mir_PepSeq$dRBM <- with(all_mir_PepSeq, ifelse(all_mir_PepSeq$query %in% dsRNA$hgnc_symbol, 1, 0))


# Nuclease activity GO term
nuclease <- getBM(attributes = c('hgnc_symbol'), 
                  filters = 'go', 
                  values = 'GO:0004518', 
                  mart = ensembl)

all_mir_PepSeq$Nuclease <- with(all_mir_PepSeq, ifelse(all_mir_PepSeq$query %in% nuclease$hgnc_symbol, 1, 0))


# RNA helicase activity GO term
rnahelicase <- getBM(attributes = c('hgnc_symbol'), 
                     filters = 'go', 
                     values = 'GO:0003724', 
                     mart = ensembl)


all_mir_PepSeq$Helicase <- with(all_mir_PepSeq, ifelse(all_mir_PepSeq$query %in% rnahelicase$hgnc_symbol, 1, 0))

# RNA destablization GO term
rnadest <- getBM(attributes = c('hgnc_symbol'), 
                 filters = 'go', 
                 values = 'GO:0006402', 
                 mart = ensembl)


all_mir_PepSeq$RNAstabilitydecay <- with(all_mir_PepSeq, ifelse(all_mir_PepSeq$query %in% rnadest$hgnc_symbol, 1, 0))


# Identifying proteins that interact with proteins involved in miRNA metabolic processes and adjusting binder-effector class
library(dplyr)
library(readr)
df <- list.files(path="C:/Users/sunid/Box/DuttaLab/deepLearning/", full.names = TRUE,pattern="*.tab3.txt") %>%
  lapply(read.table) %>%
  bind_rows

all_mir_PepSeq$miRProtInt <- with(all_mir_PepSeq, ifelse(all_mir_PepSeq$query %in% df$V1, 1, 0))

all_mir_PepSeq <- all_mir_PepSeq %>%
  mutate(effector = case_when(miRProtInt == 1 & padj < 0.05 & effector == "BinderEffector" ~ 'BinderEffector',
                              miRProtInt == 1 & padj < 0.05 & effector == "notBinderEffector"  ~ 'BinderEffector',
                              miRProtInt == 1 & padj > 0.05 ~ 'BinderNotEffector',
                              microRNAprocessing == 1 & padj < 0.05 ~ 'BinderEffector',
                              microRNAprocessing == 1 & padj > 0.05 ~ 'BindernotEffector',
                              TRUE ~ effector))


# all_mirEditedTreiber <- readRDS("C:/Users/sunid/Box/DuttaLab/deepLearning/Treiber101021all_mirEdited.rds")


names(pc_prot) <- gsub("hgnc_symbol","query", names(pc_prot))


srnaEff2 <- merge(all_mir_PepSeq,pc_prot, "query")
mirs <- unique(srnaEff2$miRNAmat)
# saveRDS(mirs,"top_mirs.rds")

names(srnaEff2)


names(srnaEff2) <- gsub("srnaEff2.effector","effector", names(srnaEff2))
head(srnaEff2)

srnaEff2 <- data.frame(srnaEff2 %>% dplyr::select(-c(query,miRNA,miRNAmat,seq)))
head(srnaEff2)
dim(srnaEff2)

# Define the training, testing, and validation sets (e.g. 80% of the data for training)

set.seed(10)
trainIndex <- createDataPartition(srnaEff2[,1], p = .8, list = F,
                                  times = 1)

train_data <- srnaEff2[rbind(trainIndex,498),]

valIndex <- createDataPartition(srnaEff2[-trainIndex,][,1], p = .5, list = F,
                                times = 1)


val_data <- srnaEff2[-trainIndex,][valIndex,]
test_data <- srnaEff2[-trainIndex,][-valIndex,]

# Define the F1 metric

f1 <- function(data, lev = NULL, model = NULL) {
  f1_val <- MLmetrics::F1_Score(y_pred = data$pred,
                                y_true = data$obs,
                                positive = lev[1])
  c(F1 = f1_val)
} 

###### final ML algorithm, ranger
set.seed(123)
nnetfit = caret::train(y = as.factor(train_data$effector), x = train_data %>% dplyr::select(-c(effector)),
                       method = "ranger",
                       preProc =  c('center','scale'),
                       trControl = trainControl(method = "cv", verboseIter = TRUE, returnData = FALSE,classProbs=TRUE,savePredictions = TRUE, sampling = "up", summaryFunction=f1))
nnetfit

# saveRDS(nnetfit,"multiclass_nnetfit_mirprtints10042022.rds")

nnetfit <- readRDS("multiclass_nnetfit_mirprtints10042022.rds")

#### test_data 
ps <- data.frame(predict(nnetfit, test_data %>% dplyr::select(-c(effector)), type = "prob"))


head(ps)
pred <- ps
pred$effector <- as.character('0')
pred$effector <- (colnames(ps)[apply(ps,1,which.max)])
summary(as.factor(pred$effector))
caret::confusionMatrix(as.factor(pred$effector), as.factor(test_data$effector), mode = "everything")

##########
pred3 <- pred
pred3$effector <- gsub('\\<BinderEffector\\>',1,pred3$effector)
pred3$effector <- gsub('\\<BinderNotEffector\\>',2,pred3$effector)
pred3$effector <- gsub('\\<notBinderEffector\\>',3,pred3$effector)
pred3$effector <- gsub('\\<notBinderNotEffector\\>',4,pred3$effector)

test_data3 <- test_data
test_data3$effector <- gsub('\\<BinderEffector\\>',1,test_data3$effector)
test_data3$effector <- gsub('\\<BinderNotEffector\\>',2,test_data3$effector)
test_data3$effector <- gsub('\\<notBinderEffector\\>',3,test_data3$effector)
test_data3$effector <- gsub('\\<notBinderNotEffector\\>',4,test_data3$effector)

library(pROC)
roc.multi <- multiclass.roc(as.numeric(pred3$effector), as.numeric(test_data3$effector))

roc.multi

rs <- roc.multi[['rocs']]
par(pty="s")
plot.roc(rs[[1]])
sapply(2:length(rs),function(i) lines.roc(rs[[i]],col=i))



# Validation set and confusion matrix
ps <- data.frame(predict(nnetfit, val_data %>% dplyr::select(-c(effector)), type = "prob"))


head(ps)
pred <- ps
pred$effector <- as.character('0')
pred$effector <- (colnames(ps)[apply(ps,1,which.max)])
summary(as.factor(pred$effector))
caret::confusionMatrix(as.factor(pred$effector), as.factor(val_data$effector), mode = "everything")

# ROC curve

pred3 <- pred
pred3$effector <- gsub('\\<BinderEffector\\>',1,pred3$effector)
pred3$effector <- gsub('\\<BinderNotEffector\\>',2,pred3$effector)
pred3$effector <- gsub('\\<notBinderEffector\\>',3,pred3$effector)
pred3$effector <- gsub('\\<notBinderNotEffector\\>',4,pred3$effector)

val_data3 <- val_data
val_data3$effector <- gsub('\\<BinderEffector\\>',1,val_data3$effector)
val_data3$effector <- gsub('\\<BinderNotEffector\\>',2,val_data3$effector)
val_data3$effector <- gsub('\\<notBinderEffector\\>',3,val_data3$effector)
val_data3$effector <- gsub('\\<notBinderNotEffector\\>',4,val_data3$effector)

library(pROC)
roc.multi <- multiclass.roc(as.numeric(pred3$effector), as.numeric(val_data3$effector))

roc.multi

rs <- roc.multi[['rocs']]
par(pty="s")
plot.roc(rs[[1]])
sapply(2:length(rs),function(i) lines.roc(rs[[i]],col=i))