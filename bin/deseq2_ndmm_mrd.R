
# Differential expression analysis based on: 
# https://www.bioconductor.org/packages/devel/workflows/vignettes/rnaseqGene/inst/doc/rnaseqGene.html#summarizedexperiment

library("tximeta")
library("DESeq2")
library("dplyr")
library("ggplot2")
library("AnnotationDbi")
library("genefilter")
library("pheatmap")
library("ggrepel")
library("biomaRt")
library("DEGreport")


# sample_metadata <- read.table("data/raw/samples.csv", sep=',', header=TRUE, stringsAsFactors = FALSE)
# sample_metadata$files <- paste0("data/processed/all_samples","/", sample_metadata$names, "/", sample_metadata$names, "_salmon/quant.sf")


# Parse arguments from a command line
args = commandArgs(trailingOnly=TRUE)


# Read metadata table
sample_metadata <- read.table(args[1], sep=',', header=TRUE, stringsAsFactors = FALSE)



# Ad $files columns with path to the quant.sf files (result from salmon)
sample_metadata$files <- paste0(args[3],"/processed/", sample_metadata$names, "/salmon/quant.sf")
head(sample_metadata)

# Select only fresh (no cryo preserved) sample 
filtered_samples <- subset(sample_metadata, specimen != "RRMM")




suppressPackageStartupMessages(library(tximeta))
makeLinkedTxome(indexDir="/scratch/project/open-27-18/salmon_gencode/salmon_index",
                source="GENCODE",
                organism="Homo sapiens",
                release="13",
                genome="GRCh38",
                fasta="/scratch/project/open-27-18/gencode/gencode.v36.transcripts.fa",
                gtf="/scratch/project/open-27-18/gencode/gencode.v36.annotation.gtf",
                write=TRUE)



# Parse salmon result with tximeta
se <- tximeta(filtered_samples)

#se <- tximport(sample_metadata$files, type = "salmon", tx2gene = tx2gene)

# summarize transcript counts to genes
gse <- summarizeToGene(se)
row.names(gse) = substr(rownames(gse), 1, 15)


# Construction of DESeqDataSet object
dds <- DESeqDataSet(gse, design = ~ patient + specimen)
dds$specimen <- relevel(dds$specimen, ref = "NDMM")



# Remove rows of the DESeqDataSet that have no counts, or only a single count across all samples
keep <- rowSums(counts(dds)) > 1
dds <- dds[keep,]
head(dds)


save(dds, file = paste0(args[3], "/results/", args[2], "/dds.RData"))


# EXPLORATORY ANALYSIS AND VISUALIZATION
########################################
# counts normalization
rld <- rlog(dds, blind = FALSE)

# principal component analysis

pcaData <- plotPCA(rld, intgroup = c( "specimen", "patient"), returnData = TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))

pdf(paste0(args[3], "/results/", args[2], "/pca_rld.pdf"))
ggplot(pcaData, aes(x = PC1, y = PC2, shape = specimen, color = patient)) +
  geom_point(size =3) +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  coord_fixed() +
  ggtitle("PCA with rld data")
dev.off()


# Differential expression analysis
##################################


 # differential expression pipeline
dds <- DESeq(dds, betaPrior=FALSE)



pdf(paste0(args[3], "/results/", args[2], "/dispersion.pdf"))
plotDispEsts(dds)
dev.off()

pdf(paste0(args[3], "/results/", args[2], "/cooks.pdf"))
par(mar=c(8,5,2,2))
boxplot(log10(assays(dds)[["cooks"]]), range=0, las=2)
dev.off()


res <- results(dds, name="specimen_MRD_vs_NDMM", alpha=0.1)
# or to shrink log fold changes association with condition:
# res <- lfcShrink(dds, coef="specimen_NDMM_vs_MRD", type="apeglm")
# res

# out_cooks = assays(dds)[["cooks"]]
# write.csv(out_cooks, file=paste0(args[3], "/results/", args[2], "/cooks.csv"))

# filter our hits with a baseMean lower than 10
res = subset(res, res$baseMean > 10)


mart <- useMart("ENSEMBL_MART_ENSEMBL")
mart <- useDataset("hsapiens_gene_ensembl", mart)

annotLookup <- getBM(
  mart=mart,
  attributes=c("ensembl_gene_id", "gene_biotype", "external_gene_name", "description"),
  filter="ensembl_gene_id",
  values=rownames(res),
  uniqueRows=TRUE)


#############################
############################
# genes with no annotation
problematic_ids = setdiff(rownames(res), annotLookup$ensembl_gene_id)

fileConn <- file(paste0(args[3], "/results/", args[2], "/problematic_ids.out"))
writeLines(problematic_ids, fileConn)    
close(fileConn)

res <- subset(res, !(rownames(res) %in% problematic_ids))
ens = rownames(res)

annotLookup <- data.frame(
  ens[match(annotLookup$ensembl_gene_id, ens)],
  annotLookup)


colnames(annotLookup) <- c( 
      "original_id", 
      c("ensembl_gene_id", "gene_biotype", "external_gene_name", "description")) 

res$symbol <- annotLookup$external_gene_name
res$biotype <- annotLookup$gene_biotype
res$description <- annotLookup$description



# order according to padj
resOrdered <- res[order(res$padj),]
resOrderedDF <- as.data.frame(resOrdered)
write.csv(resOrderedDF, file = paste0(args[3], "/results/", args[2], "/deseq_result_paired.csv"))

counts_ = counts(dds, normalized=T)
write.csv(counts_, file=paste0(args[3], "/results/", args[2], "/norm_counts.csv"))




# plot heatmap with all significant genes
pdf(paste0(args[3], "/results/", args[2], "/significant_genes_heatmap.pdf"),  height=14)
resSig <- subset(res, padj < 0.1)
resSig <- subset(resSig, biotype != "IG_V_gene")
resSig <- subset(resSig, biotype != "IG_C_gene")
# resSig <- subset(resSig, biotype == "protein_coding")
resSig <- subset(resSig, abs(log2FoldChange) > 1)
select <- rownames(resSig[ order(resSig$padj), ], )
mat = assay(rld)[select,]
mat = t(scale(t(mat)))
mat = head(mat, 50)
write.csv(mat, file=paste0(args[3], "/results/", args[2], "/heatmap_matrix.csv"))
anno <- as.data.frame(colData(rld)["specimen"])
pheatmap(mat, cluster_rows=TRUE, show_rownames=TRUE, 
         cluster_cols=TRUE, annotation_col=anno, show_colnames=TRUE)
dev.off()

#heatmap with geneids
pdf(paste0(args[3], "/results/", args[2], "/significant_genes_heatmap_geneid.pdf"),  height=14)
resSig <- subset(res, padj < 0.1)
resSig <- subset(resSig, biotype != "IG_V_gene")
resSig <- subset(resSig, biotype != "IG_C_gene")
resSig <- subset(resSig, abs(log2FoldChange) > 1)
select <- rownames(resSig[ order(resSig$padj), ], )
reorder <- resSig[order(resSig$padj), ]
mat = assay(rld)[select,]
row.names(mat) <- reorder$symbol
row.names(mat) <- ifelse(row.names(mat)=="", rownames(reorder), row.names(mat))
mat = t(scale(t(mat)))
mat = head(mat, 50)
anno <- as.data.frame(colData(rld)["specimen"])
pheatmap(mat, cluster_rows=TRUE, show_rownames=TRUE, 
         cluster_cols=TRUE, annotation_col=anno, show_colnames=TRUE)
dev.off()



pdf(paste0(args[3], "/results/", args[2], "/ma_plot.pdf"))
plotMA(res)
dev.off()

summary(res)