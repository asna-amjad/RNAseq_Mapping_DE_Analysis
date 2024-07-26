##### Make Graphs in R #####
############################

INPATH <- "/Users/Asna/Desktop/DE_analysis/"

sampleFiles <- c("E2_R1.counts", "E2_R2.counts", "Untr_R1.counts", "Untr_R2.counts")
sampleCondition <- c("treated", "treated", "untreated", "untreated")
sampleFiles
E2_R1 <- read.table("E2_R1.counts", header = F)
m_rows <- nrow(E2_R1)
m_rows

conditions <- c("treated", "treated", "untreated", "untreated")
samples <- c("E2_R1", "E2_R2", "Untr_R1", "Untr_R2")
counts = data.frame(matrix(ncol = 4, nrow = 51686))
 for (ii in sampleFiles) {
   print(ii)
   temp <- read.table(paste(INPATH, ii, sep = ""), header = F, stringsAsFactors = F)
   counts[, which(sampleFiles == ii)] <- temp[, 2]
 }
rownames(counts) <- temp[, 1]
colnames(counts) <- sampleFiles
head(counts)
cds <- DGEList(counts, group = sampleCondition)
cds
names(cds)
cds <- cds[rowSums(1e+06 * cds$counts/expandAsMatrix(cds$samples$lib.size, dim(cds)) > 1) >= 3, ]

## Normalize

cds <- calcNormFactors(cds)
cds$samples
cds <- estimateCommonDisp(cds)
cds <- estimateTagwiseDisp(cds)
names(cds)

## generate tables

de.tgw <- exactTest(cds, dispersion = cds$tagwise.dispersion,
                    pair = c("untreated", "treated"))

resultsTbl.tgw <- topTags(de.tgw, n = nrow(de.tgw$table), adjust.method = "BH")$table
de.results.tgw <- resultsTbl.tgw[resultsTbl.tgw$FDR <= 0.05, ]
dim(de.results.tgw)

## =====================
## Part 2.1
## =====================

library(edgeR)
counts <- as.matrix(read.csv("E2_gene_matrix.csv", row.names = "gene_id"))
head(counts)

## Calculation of FPKM requires gene lengths from geneCounts files

## Import the files into RStudio

E2R1tab <- read.table("E2_rep1_geneCounts.tab", header=F, stringsAsFactor=T, skip=1)
E2R2tab <- read.table("E2_rep2_geneCounts.tab", header=F, stringsAsFactor=T, skip=1)
UntrR1tab <- read.table("Untr_rep1_geneCounts.tab", header=F, stringsAsFactor=T, skip=1)
UntrR2tab <- read.table("Untr_rep2_geneCounts.tab", header=F, stringsAsFactor=T, skip=1)
 
E2r1tab <- data.frame(gene=E2R1tab[,2], FPKM=E2R1tab[,8])
E2r2tab <- data.frame(gene=E2R2tab[,2], FPKM=E2R2tab[,8])
Untrr1tab <- data.frame(gene=UntrR1tab[,2], FPKM=UntrR1tab[,8])
Untrr2tab <- data.frame(gene=UntrR2tab[,2], FPKM=UntrR2tab[,8])
head(E2r1tab)

E2merge <- merge(E2r1tab, E2r2tab, by="gene", suffixes = c("_E2r1", "_E2r2"))
head(E2merge)
E2uniq <- E2merge[!duplicated(E2merge$gene),]
head(E2uniq)

## ======================
## Part 1.3 
## ======================

qplot(log2(E2uniq[, 2]), log2(Untruniq[, 2]), xlab = "E2", ylab = "Untreated", main = "Untreated and Treated Correlation", pch=20, cex=0.2)
E2UNTR_lm = lm(log2(E2uniq[, 2]) ~ log2(Untruniq[, 2]))
abline(E2UNTR_lm[1], E2UNTR_lm[2], lwd = 2, col = "black")
data = data.frame(E2uniq[, 2], Untruniq[, 2])

Untrmerge <- merge(Untrr1tab, Untrr2tab, by="gene", suffixes = c("_Untrr1", "_Untrr2"))
Untruniq <- Untrmerge[!duplicated(Untrmerge$gene),]
head(Untruniq)

cds <- DGEList(counts, group = conditions)
names(cds)
head(cds$counts)

cds <- cds[rowSums(1e+06 * cds$counts/expandAsMatrix(cds$samples$lib.size, dim(cds)) > 5) >= 4, ]
dim(cds)
cds <- calcNormFactors(cds)
cds$samples
cds$samples$lib.size * cds$samples$norm.factors
cds <- estimateCommonDisp(cds)
names(cds)
cds$common.dispersion
cds <- estimateTagwiseDisp(cds, prior.df = 10)
head(cds)

## ======================
## Part 2.2
## ======================
plotMDS(cds, method = "bcv", col = as.numeric(cds$samples$group))
legend("top", as.character(unique(cds$samples$group)),
                col = 1:3, pch = 20, inset = 0.01, bty = "n")

## ======================
## Part 2.3
## ======================
de.cmn <- exactTest(cds, dispersion = cds$common.dispersion, pair = c("untreated", "treated"))
de.tgw <- exactTest(cds, dispersion = cds$tagwise.dispersion, pair = c("untreated", "treated"))
de.poi <- exactTest(cds, dispersion = 1e-06, pair = c("untreated", "treated"))
names(de.tgw)

resultsTbl.cmn <- topTags(de.cmn, n = nrow(de.tgw$table))$table
resultsTbl.tgw <- topTags(de.tgw, n = nrow(de.tgw$table))$table
resultsTbl.poi <- topTags(de.poi, n = nrow(de.poi$table))$table
head(resultsTbl.tgw)

de.genes.tbl.cmn <- resultsTbl.cmn[resultsTbl.cmn$FDR <= 0.05,]
de.genes.tbl.tgw <- resultsTbl.tgw[resultsTbl.tgw$FDR <= 0.05,]
de.genes.tbl.poi <- resultsTbl.poi[resultsTbl.poi$FDR <= 0.05,]
dim(de.genes.tbl.cmn)

dim(de.genes.tbl.tgw)
dim(de.genes.tbl.poi)

## ========================
## Part 2.3
## MA Plot of changed genes
## ========================

plot(resultsTbl.cmn$logCPM, resultsTbl.cmn$logFC, pch = 20, col = "black",
           cex = 0.2, main = "Common dispersion", xlab = "log(total CPM)",
           ylab = "log(Fold Change)")
points(de.genes.tbl.cmn$logCPM, de.genes.tbl.cmn$logFC, pch = 20,
           col = "violetred3", cex = 0.2)

plot(resultsTbl.tgw$logCPM, resultsTbl.tgw$logFC, pch = 20, col = "black",
           cex = 0.2, main = "Tagwise dispersion", xlab = "log(total CPM)",
           ylab = "log(Fold Change)")
points(de.genes.tbl.tgw$logCPM, de.genes.tbl.tgw$logFC, pch = 20,
           col = "violetred3", cex = 0.2)

plot(resultsTbl.poi$logCPM, resultsTbl.poi$logFC, pch = 20, col = "black",
           cex = 0.2, main = "Poisson model", xlab = "log(mean CPM)",
           ylab = "log(Fold Change)")
points(de.genes.tbl.poi$logCPM, de.genes.tbl.poi$logFC, pch = 20,
           col = "violetred3", cex = 0.2)

## ==========================================
## Part 2.4
## Make a text file for all genes from table
## ==========================================

write.table(resultsTbl.tgw, file = "all_genes_tgw.txt", col.names = T, row.names = T)

resultsTbl.tgw <- read.table("./all_genes_tgw.txt")
E2_up <- de.genes.tbl.tgw[de.genes.tbl.tgw$logFC > 0, ]
dim(E2_up)

E2_down <- de.genes.tbl.tgw[de.genes.tbl.tgw$logFC < 0, ]
dim(E2_down)
head(E2_up)

## Remove transcript suffix
## ===========================

Genes.all <- rownames(resultsTbl.tgw)
E2_up_names <- rownames(E2_up)
E2_down_names <- rownames(E2_down)
#Remove everything before the “|” and will only keep the NM#
Genes.all_split <- sapply(Genes.all, function(x) unlist(strsplit(x, split = "[|]"))[2])
E2_up_split <- sapply(E2_up_names, function(x) unlist(strsplit(x, split = "[|]"))[2])
E2_down_split <- sapply(E2_down_names, function(x) unlist(strsplit(x, split = "[|]"))[2])
head(E2_up_split)

up.vector = as.integer(Genes.all_split %in% E2_up_split)
down.vector = as.integer(Genes.all_split %in% E2_down_split)
head(up.vector)

names(up.vector) <- Genes.all_split
names(down.vector) <- Genes.all_split

head(up.vector)

## Need to convert REFSEQ to ENTREZID
## ====================================

library(AnnotationDbi)
library(org.Hs.eg.db)
Genes.all_splitc <- as.character(Genes.all_split)
E2_up_splitc <- as.character(E2_up_split)
E2_down_splitc <- as.character(E2_down_split)
all_entrezID <- unique(mapIds(x = org.Hs.eg.db, column = "ENTREZID", keytype = "REFSEQ", multiVals = "first", keys = Genes.all_splitc))
up_entrezID <- unique(mapIds(x = org.Hs.eg.db, column = "ENTREZID", keytype = "REFSEQ", multiVals = "first", keys = E2_up_splitc))
down_entrezID <- unique(mapIds(x = org.Hs.eg.db, column = "ENTREZID", keytype = "REFSEQ", multiVals = "first", keys = E2_down_splitc))
head(down_entrezID)

head(E2up_splitc)
E2up_df <- as.data.frame(E2_up_splitc)
E2down_df <- as.data.frame(E2_down_splitc)
all_df <- as.data.frame(Genes.all_splitc)
head(E2up_df)

write.table(E2up_df, file="E2_upreg_genes.txt", col.names=F, row.names=F, quote=F)
write.table(E2up_df, file="E2_downreg_genes.txt", col.names=F, row.names=F, quote=F)
write.table(all_df, file="all_genes.txt", col.names=F, row.names=F, quote=F)

## =========================================
## Part 2.5
## Perform GO analysis using ClusterProfiler
## =========================================
library(clusterProfiler)

## BP
up_go_enrich <- enrichGO(gene = up_entrezID, universe = all_entrezID,
                         OrgDb = "org.Hs.eg.db", keyType = "ENTREZID", readable = T,
                         ont = "BP", pAdjustMethod = "fdr", qvalueCutoff = 0.05)
down_go_enrich <- enrichGO(gene = down_entrezID, universe = all_entrezID,
                           OrgDb = "org.Hs.eg.db", keyType = "ENTREZID", readable = T,
                           ont = "BP", pAdjustMethod = "fdr", qvalueCutoff = 0.05)

barplot(up_go_enrich, showCategory = 20, title = "BP GO Up enrichment",
        font.size = 7.5, label_format = 20)

barplot(down_go_enrich, showCategory = 20, title = "BP GO down enrichment",
        font.size = 7.5, label_format = 20)

## CC
up_go_enrich <- enrichGO(gene = up_entrezID, universe = all_entrezID,
                         OrgDb = "org.Hs.eg.db", keyType = "ENTREZID", readable = T,
                         ont = "CC", pAdjustMethod = "fdr", qvalueCutoff = 0.05)
down_go_enrich <- enrichGO(gene = down_entrezID, universe = all_entrezID,
                           OrgDb = "org.Hs.eg.db", keyType = "ENTREZID", readable = T,
                           ont = "CC", pAdjustMethod = "fdr", qvalueCutoff = 0.05)

barplot(up_go_enrich, showCategory = 20, title = "CC GO Up enrichment",
        font.size = 7.5, label_format = 20)

barplot(down_go_enrich, showCategory = 20, title = "CC GO down enrichment",
        font.size = 7.5, label_format = 20)

## MF
up_go_enrich <- enrichGO(gene = up_entrezID, universe = all_entrezID,
                         OrgDb = "org.Hs.eg.db", keyType = "ENTREZID", readable = T,
                         ont = "MF", pAdjustMethod = "fdr", qvalueCutoff = 0.05)
down_go_enrich <- enrichGO(gene = down_entrezID, universe = all_entrezID,
                           OrgDb = "org.Hs.eg.db", keyType = "ENTREZID", readable = T,
                           ont = "MF", pAdjustMethod = "fdr", qvalueCutoff = 0.05)


barplot(up_go_enrich, showCategory = 20, title = "MF GO Up enrichment",
        font.size = 7.5, label_format = 20)

barplot(down_go_enrich, showCategory = 20, title = "MF GO down enrichment",
        font.size = 7.5, label_format = 20)

## =======================
## Part 3.1
## =======================

sampleFiles <- c("downgenesclosest.txt", "upgenesclosest.txt","nonregulatedgenesclosest.txt")
sampleCondition <- c("down", "up", "non")
downreg <- read.table("downgenesclosest.txt", header = T)
upreg <- read.table("upgenesclosest.txt", header = T)
nonreg <- read.table("nonregulatedgenesclosest.txt", header = T)

## Find the distances to the peaks

upgenes_cdf <- ecdf(upreg[,12])
downgenes_cdf <- ecdf(downreg[,12])
nongenes_cdf <- ecdf(nonreg[,12])

## Create CDF Plot 

plot(upgenes_cdf, col='steelblue', xlim= c(0, 60000),
     main = "CDF of E2_ER distance to TSSs", ylab='fraction', xlab="Distance (bp)")
abline(v=median(upgenes_cdf[,12]), lty=2)
abline(h=.5, lty=2)
lines(downgenes_cdf,col="violetred",lty=2)
lines(nonRgenes_cdf,col="green",lty=2)
legend("topleft", c("Up Reg", "Down Reg", "Non Reg"), inset=0.03, lty=1, col=c("steelblue", "violetred", "green"), cex=1)

## Make boxplots for regulated genes

boxplot(upreg[,12], outline=F, col='steelblue', lwd=2,
        xlab='E2 peaks', ylab='distance (bp)',
        main = "Distance of peaks from UpReg gene TSSs")
boxplot(downreg[,12], outline=F, col='violetred', lwd=2,
        xlab='E2 peaks', ylab='distance (bp)',
        main = "Distance of peaks from DownReg gene TSSs")
boxplot(nonR[,12], outline=F, col='green', lwd=2,
        xlab='E2 peaks', ylab='distance (bp)',
        main = "Distance of peaks from NonReg gene TSSs")

## =============================
