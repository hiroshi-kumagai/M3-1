# libraries
library(tidyverse)
library(tximport)
library(limma)
library(DESeq2)
library(org.Mm.eg.db)
library(org.Hs.eg.db)

setwd("/Volumes/HDD14TB/RNAseq_expression/Dataset/Sarcopenia_Chinese_new")

# directory paths
dir <- "/Volumes/HDD14TB/RNAseq_expression/Dataset/Sarcopenia_Chinese_new"
datadir <- paste(dir, "data", sep = "/")
resultsdir <- paste(dir, "results", sep = "/")

# loading in RSEM data
file <- paste(datadir, "sampleinfo_MCSNP_muscle.csv", sep = "/")
sampleinfo <- read.csv(file, header = TRUE)
sampleinfo$MC_SNP <- factor(sampleinfo$MC_SNP, 
                               levels = c("WT", "K14Q"))

# load in files for tximport
files <- sampleinfo$path
names(files) <- sampleinfo$run
head(files)
all(file.exists(files))

txi <- tximport::tximport(files, type = "rsem", txIn = FALSE, txOut = FALSE)
saveRDS(txi, "txi.MOTScSNP_China.new.rds") # for later use

file <- paste(datadir, "txi.MOTScSNP_China.new.rds", sep = "/")
txi <- readRDS(file)
txi$length[txi$length == 0] <- 1
txi$counts[1:5, 1:5]
txi$length[1:5, 1:5]

# DESeq2
all(sampleinfo$run == colnames(txi$counts))
rownames(sampleinfo) <- colnames(txi$counts)
dds <- DESeq2::DESeqDataSetFromTximport(txi, sampleinfo, ~ MC_SNP)
dds <- DESeq2::DESeq(dds)
dds <- BiocGenerics::estimateSizeFactors(dds)

# get GENENAME and SYMBOL for dds results
res_annotation <- clusterProfiler::bitr(rownames(dds),
                                        fromType = "ENSEMBL",
                                        toType = c("GENENAME", "SYMBOL"),
                                        OrgDb = org.Hs.eg.db,
                                        drop = FALSE)

# group and annotate the results
res_WT_K14Q <- data.frame(DESeq2::results(dds,contrast=c("MC_SNP", "WT", "K14Q")))
write.csv(res_WT_K14Q,"res_WT_K14Q.csv")

# create a list of DESseq2 results
deseq_results <- res_WT_K14Q

# merge with gene annotation data
deseq_results$ENSEMBL <- rownames(deseq_results)
deseq_results <- merge(deseq_results, res_annotation, by = "ENSEMBL")
write.csv(deseq_results,"res_WT_K14Q_genename.csv")

