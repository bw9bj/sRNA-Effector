# args[1]==gene of interest
# args[2]==log2 FC file (experimental/control)

suppressPackageStartupMessages({
library(tidyverse)
library(caret)
library(PepTools)
library(Peptides)
library(fgsea)
library(optparse)
})

sewd <- getwd()
args <- commandArgs(trailingOnly = TRUE)

option_list = list(
		     make_option(c("-g", "--gene"), type="character", default=NULL, help="Gene of interest", metavar="character"),
		     make_option(c("-f", "--file"), type="character", default=NULL, help="Log2FC input file name", metavar="character"),
		     make_option(c("-m", "--miRNA"), type="character", default="0", help="MiRNA of interest. Example format: miR-21-5p", metavar="character")
		       ) 
 
opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

# read in sRNA-Effector model
fit3 <- readRDS("multiclass_nnetfit_mirprtints10042022.rds")


# mir_gt2 <- read.table(paste0("C:/Users/sunid/Downloads/",Var1[[1]]), sep = "\t", header = T)

# gene of interest
gene <- opt$gene
file <- opt$file
cat(paste("\n\n**********Identifying miRNAs bound and/or affected by",gene,"**********"))
cat(paste("\n\nReading in file",file))

# read in file with the format gene symbol (column 1), log2FC (experimental/control) (column 2)
mir_gt <- read.table(paste0(sewd,'/',file), header = T)

colnames(mir_gt) <- c("Gene.symbol","logFC")

mir_gmt <- gmtPathways(paste0(sewd,'/',"c3.mir.mirdb.v7.4.symbols.gmt"))

if(nchar(opt$miRNA) > 1){
        miRinput <- opt$miRNA
cat(paste("\n\n Generating CDF plot for", miRinput))
        miRinput <- toupper(gsub("-","",miRinput))
        miRinput <- toupper(gsub("3P","_3P",miRinput))
        miRinput <- toupper(gsub("5P","_5P",miRinput))
        miRinput <- (gsub("HSA","",miRinput))

        mtar <- data.frame(mir_gmt[miRinput])
        colnames(mtar) <- "mRNA"
        ntar <- mir_gt[!mir_gt$Gene.symbol%in%mtar$mRNA,]$logFC
        tar <- mir_gt[mir_gt$Gene.symbol%in%mtar$mRNA,]$logFC
        cent_ntar <- ntar-median(ntar)
        cent_tar <- tar-median(ntar)
               if(length(cent_tar >= 10)){
                        wt <- wilcox.test(cent_ntar,cent_tar)
                        pval <- as.numeric(wt$p.value)
                        Dval <- as.numeric(ks.test(cent_ntar,cent_tar, exact = T)$statistic)
                        ES <- median(cent_tar) - median(cent_ntar)
                        ESmax <- (quantile(cent_tar,.9) - quantile(cent_ntar,.9))[[1]]
                        ESmin <- (quantile(cent_tar,.1) - quantile(cent_ntar,.1))[[1]]
                        ESmean <- mean(cent_tar) - mean(cent_ntar)
                        SDntar <- sd(cent_ntar)
                        SDtar <- sd(cent_tar)
                        pdf(paste0(miRinput,"_targets.pdf"))
                        plot(ecdf(cent_ntar),verticals = TRUE,do.p = FALSE, col=c("black"),ylab="Fraction of genes",xlab="Log fold change (Experimental/Control)",lwd=2,xlim=c(-1.5,1.5), main = paste(opt$miRNA,"Targets"))
                        text(-.25,0.9,labels = paste0("p-value = ",formatC(pval, format = "e", digits = 2)))
                        text(-.25,0.85,labels = paste0("ES = ",round(ES,3)))
                        text(-.25,0.8,labels = paste0("Targets = ",length(tar)),col = "red")
                        plot(ecdf(cent_tar),verticals = TRUE,do.p = FALSE,col=c("red"),add=T,lwd=2)
                        invisible(dev.off())
               }
}


cat("\n\nCalculating effect sizes...")
all_mirsAgo2 <- plyr::ldply(names(mir_gmt), function(x){
			      set.seed(123)	    
			      mtar <- data.frame(mir_gmt[x])
			      colnames(mtar) <- "mRNA"

			      ntar <- mir_gt[!mir_gt$Gene.symbol%in%mtar$mRNA,]$logFC
			      tar <- mir_gt[mir_gt$Gene.symbol%in%mtar$mRNA,]$logFC
			      cent_ntar <- ntar-median(ntar)
			      cent_tar <- tar-median(ntar)
			      

			      if(length(cent_tar >= 10)){
				      wt <- wilcox.test(cent_ntar,cent_tar)
				      pval <- as.numeric(wt$p.value)

				      Dval <- as.numeric(ks.test(cent_ntar,cent_tar, exact = T)$statistic)

				      ES <- median(cent_tar) - median(cent_ntar)

				      ESmax <- (quantile(cent_tar,.9) - quantile(cent_ntar,.9))[[1]]

				      ESmin <- (quantile(cent_tar,.1) - quantile(cent_ntar,.1))[[1]]

				      ESmean <- mean(cent_tar) - mean(cent_ntar)

				      SDntar <- sd(cent_ntar)

				      SDtar <- sd(cent_tar)

				      Vars <- data.frame(gene,x)

				      df <- cbind(pval,Dval,ES,ESmax,ESmin,ESmean,SDntar,SDtar)

				      data.frame(Vars,df)


}}, .progress = "text")

padj <- p.adjust(all_mirsAgo2$pval, method = "fdr")
all_mirsAgo2$padj <- padj

colnames(all_mirsAgo2) <- gsub("x","miRNA", colnames(all_mirsAgo2))


cat(paste("\n\nCurating features for",gene,"\n\n"))

mir_lvl <- read.table(paste0(sewd,"/","high_conf_mature_miRs.txt"), header = F)

high_mir_lvl <- mir_lvl

high_mir_lvl$precursor <- gsub("-5p","", high_mir_lvl$V1)
high_mir_lvl$precursor <- gsub("-3p","", high_mir_lvl$precursor)

high_mir_lvl$miRNA <- toupper(gsub("-","",high_mir_lvl$V1))
high_mir_lvl$miRNA <- toupper(gsub("3P","_3P",high_mir_lvl$miRNA))
high_mir_lvl$miRNA <- toupper(gsub("5P","_5P",high_mir_lvl$miRNA))
high_mir_lvl$miRNA <- (gsub("HSA","",high_mir_lvl$miRNA))

all_mirsAgo2 <- all_mirsAgo2[all_mirsAgo2$miRNA %in% high_mir_lvl$miRNA,]

miRseqs <- read.csv("edited-nkim_precursor_miRNAseqs.csv")[,c(1,2)]
colnames(miRseqs) <- c("miRNA","seq")

miRseqs$miRNA <- gsub("-5p","_5P",miRseqs$miRNA)
miRseqs$miRNA <- gsub("-3p","_3P",miRseqs$miRNA)
miRseqs$miRNA <- gsub("hsa-","",miRseqs$miRNA)
miRseqs$miRNA <- sub("(_[^_]+).*", "", miRseqs$miRNA)
miRseqs$miRNA <- toupper(gsub("-","",miRseqs$miRNA))
miRseqs$miRNA <- toupper(miRseqs$miRNA)

all_mirsAgo2$miRNAmat <- all_mirsAgo2$miRNA
all_mirsAgo2$miRNA <- gsub("_5P","",all_mirsAgo2$miRNA)
all_mirsAgo2$miRNA <- gsub("_3P","",all_mirsAgo2$miRNA)
all_mirsAgo2$miRNA <- gsub("LET_7","LET7",all_mirsAgo2$miRNA)


all_mirsAgo2 <- unique(merge(all_mirsAgo2, miRseqs, by = "miRNA"))


suppressWarnings({
ES_mirSeqSplit <- strsplit((all_mirsAgo2$seq), "")
t2 <- data.frame(lapply((ES_mirSeqSplit), type.convert), stringsAsFactors=FALSE)

ES_mir1hot <- map_df(t2, ~ gsub("A", "1", .x))
ES_mir1hot <- map_df(ES_mir1hot, ~ gsub("T", "2", .x))
ES_mir1hot <- map_df(ES_mir1hot, ~ gsub("C", "3", .x))
ES_mir1hot <- map_df(ES_mir1hot, ~ gsub("G", "4", .x))

ES_mir1hot <- data.matrix(t(ES_mir1hot))

ES_mir1hot <- data.frame(matrix(as.numeric(ES_mir1hot),    # Convert to numeric matrix
				                                ncol = ncol(ES_mir1hot)),data.frame(all_mirsAgo2$miRNA))[,-126]

})
all_mirsAgo2 <- unique(cbind(all_mirsAgo2,ES_mir1hot))

# prot seqs
suppressPackageStartupMessages(library(biomaRt))
ensembl <- useEnsembl(biomart = "ensembl", 
		                            dataset = "hsapiens_gene_ensembl")

seq <-  getSequence(id = (as.character(all_mirsAgo2$gene)),
		                        type = "hgnc_symbol",
					                    seqType = "peptide",
					                    mart = ensembl)
seqNoNA <- subset(seq, peptide != "Sequence unavailable")
suppressPackageStartupMessages(library(data.table))
seqNoNA$peptide <- gsub("\\*","",seqNoNA$peptide)
seqNoNA <- seqNoNA[!seqNoNA$peptide %like% "X",]

suppressPackageStartupMessages(library(Peptides))
ai <- mean(aIndex(seqNoNA$peptide))
bm <- mean(boman(seqNoNA$peptide))
char <- mean(charge(seqNoNA$peptide))
hydrophoIn <- mean(hydrophobicity(seqNoNA$peptide))
instab <- mean(instaIndex(seqNoNA$peptide))
#mem <- membpos(PEPseqs$peptide)$MembPos
molwgt <- mean(mw(seqNoNA$peptide))
isoPt <- mean(pI(seqNoNA$peptide))

seqNoNA$crucianiProp <- plyr::ldply(seqNoNA$peptide,function(x){
				        matrix(unlist(crucianiProperties(x)), ncol = 3, nrow = 1)
							    })
			
seqNoNA$Kidera <- plyr::ldply(seqNoNA$peptide,function(x){
				  matrix(unlist(kideraFactors(x)), ncol = 10, nrow = 1)
							    })

seqNoNA$FP <- plyr::ldply(seqNoNA$peptide,function(x){
			      matrix(unlist(protFP(x)), ncol = 8, nrow = 1)
							    })
seqNoNA <- seqNoNA %>% dplyr::select(-peptide) 

seqNoNA <- data.frame(seqNoNA,ai,bm,char,hydrophoIn,instab,molwgt,isoPt)

DF <- as.data.table(seqNoNA)

pc_prot <- DF %>%
	  group_by(hgnc_symbol) %>%
	    summarize_all(mean)

    # pulling out prots with important domains/Go functions


    # reg of translation
    transl <- getBM(attributes = c('hgnc_symbol'), 
		                    filters = 'go', 
				                    values = 'GO:0006417', 
				                    mart = ensembl)

    all_mirsAgo2$Translationregulation <- with(all_mirsAgo2, ifelse(all_mirsAgo2$gene %in% transl$hgnc_symbol, 1, 0))


    # response to virus
    virus <- getBM(attributes = c('hgnc_symbol'), 
		                  filters = 'go', 
				                 values = 'GO:0009615', 
				                 mart = ensembl)

    all_mirsAgo2$ViralRNAregulation <- with(all_mirsAgo2, ifelse(all_mirsAgo2$gene %in% virus$hgnc_symbol, 1, 0))


    # mirna metabolic process
    mirmetproc <- getBM(attributes = c('hgnc_symbol'), 
			                    filters = 'go', 
					                        values = 'GO:0010586', 
					                        mart = ensembl)

    all_mirsAgo2$microRNAprocessing <- with(all_mirsAgo2, ifelse(all_mirsAgo2$gene %in% mirmetproc$hgnc_symbol, 1, 0))


    # dsRNA binding
    dsRNA <- getBM(attributes = c('hgnc_symbol'), 
		                  filters = 'go', 
				                 values = 'GO:0003725', 
				                 mart = ensembl)

    all_mirsAgo2$dRBM <- with(all_mirsAgo2, ifelse(all_mirsAgo2$gene %in% dsRNA$hgnc_symbol, 1, 0))


    # nuclease activity
    nuclease <- getBM(attributes = c('hgnc_symbol'), 
		                        filters = 'go', 
					                  values = 'GO:0004518', 
					                  mart = ensembl)

    all_mirsAgo2$Nuclease <- with(all_mirsAgo2, ifelse(all_mirsAgo2$gene %in% nuclease$hgnc_symbol, 1, 0))


    # rna helicase activity
    rnahelicase <- getBM(attributes = c('hgnc_symbol'), 
			                      filters = 'go', 
					                           values = 'GO:0003724', 
					                           mart = ensembl)


    all_mirsAgo2$Helicase <- with(all_mirsAgo2, ifelse(all_mirsAgo2$gene %in% rnahelicase$hgnc_symbol, 1, 0))

    # rna catabolic processes
    rnadest <- getBM(attributes = c('hgnc_symbol'), 
		                      filters = 'go', 
				                       values = 'GO:0006402', 
				                       mart = ensembl)


    all_mirsAgo2$RNAstabilitydecay <- with(all_mirsAgo2, ifelse(all_mirsAgo2$gene %in% rnadest$hgnc_symbol, 1, 0))

    # does the prot of interest interact with proteins involved in mirna metabolic processes?
    suppressPackageStartupMessages({
	    library(dplyr)
            library(readr)
    })




    df <- read.table(paste0(sewd,"/","miRmetProts.txt") )

    all_mirsAgo2$miRProtInt <- with(all_mirsAgo2, ifelse(all_mirsAgo2$gene %in% df$V1, 1, 0))

colnames(pc_prot) <- gsub("hgnc_symbol","gene",colnames(pc_prot))

srnaEff2 <- merge(all_mirsAgo2,pc_prot, "gene")

srnaEff2 <- na.omit(srnaEff2)

srnaEff3 <- na.omit(data.frame(srnaEff2 %>% dplyr::select(-c(gene,miRNA,miRNAmat,seq))))
names(srnaEff3) <- gsub("ESmamiRNA","ESmax", names(srnaEff3))


ps <- data.frame(predict(fit3, srnaEff3, type = "prob"))

#thresh_prob <- 0.5
pred <- ps
pred$effector <- as.character('0')

pred$effector <- colnames(ps)[apply(ps,1,which.max)]

summary(as.factor(pred$effector))

colnames(srnaEff2) <- gsub("gene","Gene",colnames(srnaEff2))

# visualization
preds2 <- unique(data.frame(pred$effector,na.omit(srnaEff2 %>% dplyr::select(c(Gene,miRNAmat,ES)))))


ord_ESrank <- preds2[order(preds2$ES),]

write.table(ord_ESrank,"predictions.txt",quote = F, row.names = F)

tab_ord_ESrank <- data.frame(table(ord_ESrank$pred.effector),gene)
colnames(tab_ord_ESrank) <- c("prediction","Freq","GOI")

pdf(paste0(gene,"_stackedbarplot.pdf"))
print(ggplot(tab_ord_ESrank, aes(fill=prediction, y=Freq, x=GOI)) + 
	  geom_bar(position="fill", stat="identity") + xlab("") + ylab("Proportion") )
invisible(dev.off())

ord_ESrankBE <- (subset(data.frame(ord_ESrank), pred.effector == "BinderEffector"))


if(dim(ord_ESrankBE)[1] > 0){
#ord_ESrankBE <- ord_ESrankBE[order(ord_ESrankBE$Freq),]

ord_ESrankBE <- tibble::rowid_to_column(ord_ESrankBE, "ID")


library(ggrepel)

# Identify the top 10 hits
Top_Hits = rbind(head(ord_ESrankBE,5), tail(ord_ESrankBE,5) )
Top_Hits
# Add column label, containing the gene name for the top hits or nothing for all others
ord_ESrankBE$label = if_else(ord_ESrankBE$miRNAmat %in% Top_Hits$miRNAmat,  
			                                  ord_ESrankBE$miRNAmat, NA_character_)

pdf("BErankplot.pdf")
print( ggplot(ord_ESrankBE,aes(ID,ES, label = miRNAmat)) + geom_point() + xlab("rank") + geom_label_repel(data=ord_ESrankBE, label = ord_ESrankBE$label, na.rm=TRUE) + theme_classic() + geom_hline(yintercept = 0, linetype = "dotted") + ylab("Effect Size") + ggtitle("Predicted Binders & Effectors") + theme(text = element_text(size = 20)) )
invisible(dev.off())
}

cat("\n\nDone!\n\n")
