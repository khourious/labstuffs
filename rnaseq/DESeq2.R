# RNA-seq analysis

# install and load R packages
#BiocManager::install(version = "3.13")
#BiocManager::install("calibrate")
#BiocManager::install("genefilter")
#BiocManager::install("gplots")
#BiocManager::install("DESeq2")
#install.packages("RColorBrewer")
library("calibrate")
library("genefilter")
library("gplots")
library("DESeq2")
library("RColorBrewer")

# clear the workspace
rm(list = ls())

# import featureCounts data
countdata <- read.table("counts.txt", header=TRUE, row.names=1)

# remove the first five columns (chr, start, end, strand, length)
countdata <- countdata[ ,6:ncol(countdata)]

# remove prefix and suffix from filename
colnames(countdata) <- gsub("X.ANALYSIS.ALIGN.", "", colnames(countdata))
colnames(countdata) <- gsub("_Aligned.out.bam", "", colnames(countdata))

# convert to matrix
countdata <- as.matrix(countdata)
head(countdata)

# assign condition
# "(rep("GROUP1", 1) - experimental group 1 and number of replicates
# "(rep("GROUP1", 1) - experimental group 2 and number of replicates
(condition <- factor(c(rep("GROUP1", 1), rep("GROUP2", 1))))

# create a coldata frame and instantiate the DESeqDataSe
(coldata <- data.frame(row.names=colnames(countdata), condition))
dds <- DESeqDataSetFromMatrix(countData=countdata, colData=coldata, design=~condition)
dds <- estimateSizeFactors(dds)
counts(dds, normalized=TRUE)
idx <- rowSums(counts(dds, normalized=TRUE) >= 5) >=3
dds <- dds[idx,]

# run the DESeq pipeline
dds <- DESeq(dds)

# plot dispersions
pdf("qc-dispersions.pdf", 50, 50, pointsize=100)
plotDispEsts(dds, main="Dispersion plot")
dev.off()

# regularized log transformation for clustering/heatmaps, etc
rld <- rlogTransformation(dds)
head(assay(rld))
hist(assay(rld))
par( mfrow = c( 1, 2 ) )
dds <- estimateSizeFactors(dds)
plot( log2( 1 + counts(dds, normalized=TRUE)[ , 1:2] ),
      col=rgb(0,0,0,.2), pch=16, cex=0.3 )
plot( assay(rld)[ , 1:2],
      col=rgb(0,0,0,.2), pch=16, cex=0.3 )

# colors for plots below
library(RColorBrewer)
(mycols <- brewer.pal(8, "Dark2")[1:length(unique(condition))])

# sample distance heatmap
sampleDists <- as.matrix(dist(t(assay(rld))))
library(gplots)
pdf("qc-heatmap-samples.pdf", w=50, h=50, pointsize=100)
heatmap.2(as.matrix(sampleDists), key=F, trace="none",
          col=colorpanel(100, "black", "white"),
          ColSideColors=mycols[condition], RowSideColors=mycols[condition],
          margin=c(10, 10), main="Sample Distance Matrix")
dev.off()

# principal components analysis
rld_pca <- function (rld, intgroup = "condition", ntop = 500, colors=NULL, legendpos="topleft", main="PCA Biplot", textcx=1, ...) {
  require(genefilter)
  require(calibrate)
  require(RColorBrewer)
  rv = rowVars(assay(rld))
  select = order(rv, decreasing=TRUE)[seq_len(min(ntop, length(rv)))]
  pca = prcomp(t(assay(rld)[select, ]))
  fac = factor(apply(as.data.frame(colData(rld)[, intgroup, drop=FALSE]), 1, paste, collapse = " : "))
  if (is.null(colors)) {
    if (nlevels(fac) >= 3) {
      colors = brewer.pal(nlevels(fac), "Paired")
    }   else {
      colors = c("black", "red")
    }
  }
  pc1var <- round(summary(pca)$importance[2,1]*100, digits=1)
  pc2var <- round(summary(pca)$importance[2,2]*100, digits=1)
  pc1lab <- paste0("PC1 (",as.character(pc1var),"%)")
  pc2lab <- paste0("PC1 (",as.character(pc2var),"%)")
  plot(PC2~PC1, data=as.data.frame(pca$x), bg=colors[fac], pch=21, xlab=pc1lab, ylab=pc2lab, main=main, ...)
  with(as.data.frame(pca$x), textxy(PC1, PC2, labs=rownames(as.data.frame(pca$x)), cex=textcx))
  legend(legendpos, legend=levels(fac), col=colors, pch=20)}
pdf("qc-pca.pdf", 50, 50, pointsize=100)
rld_pca(rld, colors=mycols, intgroup="condition", xlim=c(-35, 45))
dev.off()

# get differential expression results
res <- results(dds)
table(res$padj<0.05)

# order by adjusted p-value
res <- res[order(res$padj), ]

# merge with normalized count data
resdata <- merge(as.data.frame(res), as.data.frame(counts(dds, normalized=TRUE)), by="row.names", sort=FALSE)
names(resdata)[1] <- "Gene"
head(resdata)

# write results
write.csv(resdata, file="diffexpr-results.csv")

# examine plot of p-values
hist(res$pvalue, breaks=20, col="grey")

# examine independent filtering
attr(res, "filterThreshold")
plot(attr(res,"filterNumRej"), type="b", xlab="quantiles of baseMean", ylab="number of rejections")

# MA plot
maplot <- function (res, thresh=0.05, labelsig=TRUE, textcx=1, ...) {
  with(res, plot(baseMean, log2FoldChange, pch=20, cex=.5, log="x", ...))
  with(subset(res, padj<thresh), points(baseMean, log2FoldChange, col="red", pch=20, cex=1.5))
  if (labelsig) {
    require(calibrate)
    with(subset(res, padj<thresh), textxy(baseMean, log2FoldChange, labs=Gene, cex=textcx, col=2))
  }
}
pdf("diffexpr-maplot.pdf", 10, 10, pointsize=15)
maplot(resdata, main="MA Plot")
dev.off()

# volcano plot with "significant" genes labeled
volcanoplot <- function (resdata, lfcthresh=2, sigthresh=0.05, main="Volcano Plot", legendpos="topright", labelsig=TRUE, textcx=1, ...) {
  with(res, plot(log2FoldChange, -log10(pvalue), pch=20, main=main, ...))
  with(subset(res, padj<sigthresh ), points(log2FoldChange, -log10(pvalue), pch=20, col="red", ...))
  with(subset(res, abs(log2FoldChange)>lfcthresh), points(log2FoldChange, -log10(pvalue), pch=20, col="orange", ...))
  with(subset(res, padj<sigthresh & abs(log2FoldChange)>lfcthresh), points(log2FoldChange, -log10(pvalue), pch=20, col="green", ...))
  if (labelsig) {
    require(calibrate)
    with(subset(res, padj<sigthresh & abs(log2FoldChange)>lfcthresh), textxy(log2FoldChange, -log10(pvalue), labs=NULL, cex=textcx, ...))
  }
  legend(legendpos, xjust=1, yjust=1, legend=c(paste("FDR<",sigthresh,sep=""), paste("|LogFC|>",lfcthresh,sep=""), "both"), pch=20, col=c("red","orange","green"))
}
pdf("diffexpr-volcanoplot.pdf", 10, 10, pointsize=15)
volcanoplot(resdata, lfcthresh=1, sigthresh=0.05, textcx=.8, xlim=c(-2.3, 2))
dev.off()
