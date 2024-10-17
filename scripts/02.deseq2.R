####summarizeOverlaps
setwd('/Volumes/HDD14TB/RNAseq_expression/Dataset/biopsy/')
sampleTable<-read.csv("sampleinfo.csv")

filenames <- file.path(getwd(), paste0(sampleTable$Run, ".Aligned.sortedByCoord.out.bam"))
file.exists(filenames)

library("BiocManager")
library("S4Vectors")
library("IRanges")
library("GenomeInfoDb")
library("GenomicRanges")
library("BiocGenerics")
library("BiocParallel")
library("XVector")
library("Biostrings")
library("Rsamtools")
bamfiles<-BamFileList(filenames, yieldSize = 2000000)
seqinfo(bamfiles[1])

library("Biobase")
library("AnnotationDbi")
library("GenomicFeatures")
gtffile <- file.path(getwd(),"Homo_sapiens.GRCh38.100.gtf")
txdb <- makeTxDbFromGFF(gtffile, format = "gtf", circ_seqs = character())
txdb

ebg <- exonsBy(txdb, by="gene")
ebg

library("matrixStats")
library("DelayedArray")
library("SummarizedExperiment")
library("GenomicAlignments")
library("BiocParallel")

se <- summarizeOverlaps(features = ebg,
                        reads = bamfiles,
                        mode = "Union",
                        singleEnd = FALSE,
                        ignore.strand = FALSE,
                        fragments = FALSE,
                        inter.feature = FALSE,
                        preprocess.reads=invertStrand)

saveRDS(se, "biopsy.rds")
se <- readRDS('biopsy.rds')



####DESEQ2
ibrary("DESeq2")
library("dplyr")
colData(se) <- DataFrame(sampleTable)
colData(se)
se$m1382 <- factor(se$m1382)
se$m1382 <- relevel(se$m1382, "c")
se$sex <- factor(se$sex) ## add batch as a factor
se$age <- factor(se$age) ## add age as a factor
dds <- DESeqDataSet(se, design = ~ m1382 + sex + age)
dds <- DESeq(dds)
res <- results(dds)

write.csv(res, "Result_biopsy.csv")
