## RNA-seq analysis

# avoid creating Rplots files
if(!interactive()) pdf(NULL)

# load R packages
library("DESeq2")
library("calibrate")
library("genefilter")
library("gplots")
library("EnhancedVolcano")
library("RColorBrewer")

# import featureCounts data
countdata <- read.table("counts_DENV_ZIKV.txt", header=TRUE, row.names=1)

# remove the first five columns (chr, start, end, strand, length)
countdata <- countdata[, 6:ncol(countdata)]

# remove prefix and suffix from filename
colnames(countdata) <- gsub("X.home.laisep.rnaseq_arbovirus.ArbovirusFiocruzBA.83677594.ANALYSIS.ALIGN.", "", colnames(countdata))
colnames(countdata) <- gsub("_Aligned.out.bam", "", colnames(countdata))

# convert to matrix
countdata <- as.matrix(countdata)
head(countdata)

# assign condition
# first "(rep(x)" contains the experiment [exp]
# second "(rep(x)" are controls [ctrl])
(condition <- factor(c(rep("DENV", 8), rep("ZIKV", 12))))

# create a coldata frame
(coldata <- data.frame(row.names=colnames(countdata), condition))

# instantiate the DESeqDataSet including counts and metadata
dds <- DESeqDataSetFromMatrix(countData=countdata, colData=coldata, design=~condition)

# normalize counts
dds <- estimateSizeFactors(dds)
counts(dds, normalized=TRUE)

# idx <- rowSums(counts(dds,normalized=TRUE)>=5)>=N, which "N" is the number of samples in the smallest group [exp or ctrl]
idx <- rowSums(counts(dds, normalized=TRUE) >= 5) >=8
dds <- dds[idx,]

# run the DESeq pipeline
dds <- DESeq(dds)

# plot dispersions
pdf("qc-dispersions_DENV_ZIKV.pdf", 50, 50, pointsize=100)
plotDispEsts(dds, main="dispersion plot")
dev.off()

# transform countdata to log2 scale to decrease the differences between samples with small counts
rld <- rlogTransformation(dds)

# sample distance heatmap
pdf("qc-distance-heatmap_DENV_ZIKV.pdf", w=50, h=50, pointsize=100)
distsRL <- as.matrix(dist(t(assay(rld))))
hmcol <- colorRampPalette(brewer.pal(11,"RdYlGn"))(100)
rownames(distsRL) <- colnames(distsRL)
heatmap.2(distsRL,trace="none",col=rev(hmcol),margin=c(7,7),dendrogram="both",main="sample distance matrix")
dev.off()

# dot plot of the effect of the transformation to rlog
par(mfrow = c( 1, 2))
dds <- estimateSizeFactors(dds)
pdf("qc-rlog_DENV_ZIKV.pdf")
par(mfrow=c(1,2))
plot(log2(1+counts(dds)[,1:2]),col=rgb(0,0,0,.2),pch=16,cex=1.0,main="log2")
plot(assay(rld)[,1:2],col=rgb(0,0,0,.2),pch=16,cex=1.0,main="rlog")
dev.off()

# principal components analysis (PCA)
(mycols <- brewer.pal(8, "Dark2")[1:length(unique(condition))])
rld_pca <- function (rld, intgroup = "condition", ntop = 500, colors=NULL, legendpos="bottomright", main="PCA", textcx=1, ...) {
  require(genefilter)
  require(calibrate)
  require(RColorBrewer)
  rv = rowVars(assay(rld))
  select = order(rv, decreasing=TRUE)[seq_len(min(ntop, length(rv)))]
  pca = prcomp(t(assay(rld)[select, ]))
  fac = factor(apply(as.data.frame(colData(rld)[, intgroup, drop=FALSE]), 1, paste, collapse = " : "))
  if (is.null(colors)) {
    if (nlevels(fac) >= 3) {
      colors = brewer.pal(nlevels(fac), "Paired")}
    else {colors = c("black", "red")}}
  pc1var <- round(summary(pca)$importance[2,1]*100, digits=1)
  pc2var <- round(summary(pca)$importance[2,2]*100, digits=1)
  pc1lab <- paste0("PC1 (",as.character(pc1var),"%)")
  pc2lab <- paste0("PC1 (",as.character(pc2var),"%)")
  plot(PC2~PC1, data=as.data.frame(pca$x), bg=colors[fac], pch=21, xlab=pc1lab, ylab=pc2lab, main=main, ...)
  with(as.data.frame(pca$x), textxy(PC1, PC2, labs=rownames(as.data.frame(pca$x)), cex=textcx))
  legend(legendpos, legend=levels(fac), col=colors, pch=20)}
pdf("qc-pca_DENV_ZIKV.pdf", 50, 50, pointsize=70)
rld_pca(rld, colors=mycols, intgroup="condition", xlim=c(-35, 45))
dev.off()

# decrease the fold change noise with shrinkage function 
resShrink <- lfcShrink(dds,contrast=c("condition","DENV","ZIKV"),type="ashr")

# MA plot of the effect of the shrinkage correction
res <- results(dds)
mcols(res,use.names=TRUE)
pdf("qc-shrinkage-correction_DENV_ZIKV.pdf", 50, 50, pointsize=80)
par(mfrow=c(1,2))
plotMA(res,ylim=c(-7,7),main="unshurunken log2 fold change")
plotMA(resShrink,ylim=c(-7,7),main="shurunken log2 fold change")
dev.off()

# get differential expression results
table(resShrink$padj<0.05)

# order by adjusted p-value
resShrink <- resShrink[order(resShrink$padj), ]

# merge with normalized count data
resShrinkdata <- merge(as.data.frame(resShrink),as.data.frame(counts(dds,normalized=TRUE)),by="row.names",sort=FALSE)
names(resShrinkdata)[1] <- "GeneId"
head(resShrinkdata)

pdf("diffexprShrinkage-volcanoplot_DENV_ZIKV.pdf", 10, 10, pointsize=20)
EnhancedVolcano(resShrink, lab=rownames(resShrink), x='log2FoldChange', y='padj', pCutoff=0.05, FCcutoff=2, pointSize = 5.0, labSize = 5.0)
dev.off()

# write results
write.csv(resShrinkdata, file="diffexprShrinkage-results_DENV_ZIKV.csv")