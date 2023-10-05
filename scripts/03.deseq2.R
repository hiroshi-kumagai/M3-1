#### DESEQ2 ####


# libraries
library(tidyverse)
library(tximport)
library(limma)
library(DESeq2)
library(org.Mm.eg.db)


# directories
dir <- "/Users/Emma/GitHub/MOTSc_cast_immobilization"
datadir <- paste(dir, "data", sep = "/")
resultsdir <- paste(dir, "results", sep = "/")

# files
file <- paste(datadir, "sampleinfo.csv", sep = "/")
sampleinfo <- read.csv(file, header = TRUE)
sampleinfo$condition <- factor(sampleinfo$condition, 
                               levels = c("Non_cast", "Cast_water", "Cast_MOTSc"))

file <- paste(datadir, "txi.motsc.rsem.tpm.3gp.ncw.rds", sep = "/")
txi <- readRDS(file)
txi$length[txi$length == 0] <- 1
txi$counts[1:5, 1:5]

# deseq2
rownames(sampleinfo) <- colnames(txi$counts)
dds <- DESeq2::DESeqDataSetFromTximport(txi, sampleinfo, ~ condition)
dds <- DESeq2::DESeq(dds)
dds <- BiocGenerics::estimateSizeFactors(dds)


# get GENENAME and SYMBOL for dds results
res_annotation <- clusterProfiler::bitr(rownames(dds),
                                        fromType = "ENSEMBL",
                                        toType = c("GENENAME", "SYMBOL"),
                                        OrgDb = org.Mm.eg.db,
                                        drop = FALSE)

res_CastWater_NonCast <- data.frame(
    DESeq2::results(dds,contrast=c("condition", "Cast_water", "Non_cast"))
    )
res_CastMOTSc_CastWater <- data.frame(
    DESeq2::results(dds, contrast=c("condition", "Cast_MOTSc", "Cast_water"))
    )
deseq_results <- list(res_CastWater_NonCast = res_CastWater_NonCast, 
                      res_CastMOTSc_CastWater = res_CastMOTSc_CastWater)
lapply(deseq_results, head, 3)



deseq_results <- lapply(deseq_results, function(dt){
                        dt$ENSEMBL <- rownames(dt)
                        dt <- merge(dt, res_annotation, by = "ENSEMBL")
                        return(dt)
                        })
lapply(deseq_results, nrow)



res_CastMOTSc_CastWater_adjp005 <- subset(deseq_results[[2]], padj < 0.05)
nrow(res_CastMOTSc_CastWater_adjp005)


#| label: volcanoplot
#| warning: false
#| message: false
#| fig-cap: "Figure 4 (B): Volcano plot that compares the gene expression levels between non cast immobilization and cast immobilization control groups."
#| cap-location: margin

# using res_CastWater_NonCast
results <- deseq_results[[1]]
signum <- results %>% 
    dplyr::filter(padj < 0.05 & abs(log2FoldChange) > 0.5) %>% 
    nrow()
name <- "Volcano plot"
subtitle <- "Non cast vs Cast water"
    
EnhancedVolcano::EnhancedVolcano(results,
                                 lab = "",
                                 x = 'log2FoldChange',
                                 y = 'padj',
                                 ylim = c(0, 80),
                                 xlab = bquote(~Log[2]~ 'fold change'),
                                 ylab = bquote(~-Log[10]~adjusted~italic(P)),
                                 pCutoff = 0.05,
                                 FCcutoff = 0.5,
                                 pointSize = 2.0, 
                                 labSize = 3.0,
                                 col = c('black', 'red', 'orange', 'blue'),
                                 legendPosition = '',
                                 legendLabels = c('NS','Log2 FC','adj.p.val','adj.p.val & Log2 FC')) +
ggplot2::labs(title = name, titleLabSize = 10, subtitle = subtitle, caption = "") +
ggplot2::coord_cartesian(xlim=c(-10, 25)) +
ggplot2::scale_x_continuous(breaks=seq(-10, 25, 5)) +
ggplot2::theme(axis.line = element_line(colour = "black"),
               panel.grid.minor = element_blank(),
               panel.border = element_blank(),
               panel.background = element_blank()
               ) 



target_samples <- dds@colData@listData$condition %in% c("Cast_water", "Cast_MOTSc")
target_name <- names(deseq_results[2])
target_genes <- res_CastMOTSc_CastWater_adjp005$ENSEMBL
dds_sub <- dds[target_genes,target_samples]
print(target_name)
print(dds_sub)



vsd <- DESeq2::varianceStabilizingTransformation(dds_sub)
mat <- SummarizedExperiment::assay(vsd)

# matrix
print(mat[1:3,1:5])

# z-scores
mat <- mat - rowMeans(mat)
print(mat[1:3,1:5])


#| warning: false
#| message: false
#| label: heatmap
#| fig-cap: "Figure 4 (C): Heat map of the differentially expressed genes (FDR < 0.05) between cast immobilization control and cast immobilization with MOTS-c treatment groups."
#| cap-location: margin

annotation.df <- data.frame(condition = factor(colData(vsd)[,c("condition")],
                                               levels = c("Cast_water", "Cast_MOTSc")))
rownames(annotation.df) <- colnames(mat)
mycolors <- list(condition = c("red", "blue"))
names(mycolors$condition) <- levels(annotation.df$condition)
breaksList = seq(-1, 1, by = 0.1)

pheatmap::pheatmap(mat, cluster_cols = FALSE, 
                   annotation_col = annotation.df, 
                   annotation_names_col = FALSE, 
                   annotation_colors = mycolors,
                   show_rownames = FALSE, 
                   show_colnames = FALSE, 
                   treeheight_row = 0, 
                   cellwidth = 20,
                   breaks = breaksList,legend_breaks = c(-1, 0, 1),
                   legend_labels = c("-1", "0", "1"),
                   color = grDevices::colorRampPalette(rev(RColorBrewer::brewer.pal(n = 7, name = "RdYlBu")))(length(breaksList)))



all(res_CastMOTSc_CastWater_adjp005$ENSEMBL == rownames(mat))



# background genes
bkgd.genes.entrez <- clusterProfiler::bitr(rownames(dds), 
                                           fromType = "ENSEMBL",
                                           toType = "ENTREZID",
                                           OrgDb = org.Mm.eg.db)
nrow(bkgd.genes.entrez)
head(bkgd.genes.entrez)



# significant genes (adjusted p-value < 0.05)
sig.genes.entrez <- clusterProfiler::bitr(res_CastMOTSc_CastWater_adjp005$ENSEMBL,
                                          fromType = "ENSEMBL",
                                          toType = "ENTREZID",
                                          OrgDb = org.Mm.eg.db)
nrow(sig.genes.entrez)
head(sig.genes.entrez)



# preparing rWikiPathways terms
wp.mm.gmt <- rWikiPathways::downloadPathwayArchive(organism = "Mus musculus", 
                                                   format = "gmt",
                                                   destpath = datadir)
wp2gene <- rWikiPathways::readPathwayGMT(wp.mm.gmt)
wpid2gene <- wp2gene %>% dplyr::select(wpid,gene)
head(wpid2gene) #TERM2GENE
wpid2name <- wp2gene %>% dplyr::select(wpid,name)
head(wpid2name) #TERM2NAME



# get the pathways for significant genes
ewp.sig.wiki005 <- clusterProfiler::enricher(sig.genes.entrez[[2]],
                                             universe = bkgd.genes.entrez[[2]],
                                             # pAdjustMethod = "none",
                                             pvalueCutoff = 0.05,
                                             TERM2GENE = wpid2gene,
                                             TERM2NAME = wpid2name)

ewp.sig.wiki005 <- DOSE::setReadable(ewp.sig.wiki005, 
                                     org.Mm.eg.db, 
                                     keyType = "ENTREZID")
nrow(ewp.sig.wiki005)

# show the pathways
ewp.sig.wiki005$Description


#| label: dotplot
#| fig-height: 2.5
#| fig-cap: "Figure 4 (E): Wikipathway significant terms and related genes that were enriched by MOTS-c treatment-induced gene expression in the skeletal muscle."
#| cap-location: margin
enrichplot::dotplot(ewp.sig.wiki005, 
                    showCategory = 2, 
                    font.size = 12) +
    xlim(0, 0.2) +
    theme(text = element_text(size = 12))


#| warning: false
#| message: false

target <- res_CastMOTSc_CastWater_adjp005
OE_foldchanges <- target$log2FoldChange
names(OE_foldchanges) <- target$SYMBOL
head(OE_foldchanges)


#| warning: false
#| message: false
#| label: cnetplot
#| fig-width: 10
#| fig-height: 8
#| fig-cap: "Figure 4 (F): Wikipathway significant terms and related genes that were enriched by MOTS-c treatment-induced gene expression in the skeletal muscle."
#| cap-location: margin

par(cex = 1.5)
enrichplot::cnetplot(ewp.sig.wiki005, 
                     categorySize="pvalue", 
                     showCategory = 2, 
                     foldChange=OE_foldchanges, 
                     vertex.label.font=12,
                     circular = TRUE, colorEdge = TRUE,
                     cex_label_gene = 1.5,
                     cex_label_category = 1.5,
                     cex_category = 3,
                     cex_gene = 1.5,
                     color_category='blue') +
    ggplot2::scale_color_gradient(name = "fold change", 
                                  low='lightpink', 
                                  high='red') +
    ggplot2::theme(legend.title=element_text(size = 18),
                   legend.text=element_text(size = 18))



sessionInfo()

