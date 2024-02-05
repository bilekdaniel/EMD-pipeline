# USAGE:
# OVERVIEW: Rscript "script location" "arg1" "arg2" "arg3" "arg4"
# arg1 = path to csv file with samples.
# arg2 = variable for name of subfolder inside results folder.
# arg3 = path to txt file containing genes for exclusion. In this case exclusion of IG genes.
# arg4 = OPTIONAL; path to txt file containing genes to be excluded from the heat map. For it to work you need to uncomment appropriate code snippet at line ~200
# EXAMPLE: Rscript scripts/deseq2_contrast.R data/raw/samples_monocytes.csv test scripts/util/IG_GENES_filtered.txt 

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



# Parse arguments from a command line
args = commandArgs(trailingOnly=TRUE)

#Prepare directories
directory <- getwd()

if (!dir.exists(paste0(directory, "/results/", args[2], "/graphs"))) {
  dir.create(paste0(directory, "/results/", args[2], "/graphs"), recursive = TRUE)
  cat("Folders created","/results/", args[2], "/graphs")
} else {
  cat("Folders already exists:", "/results/", args[2], "/graphs")
}
print(paste0(directory, "/results/", args[2], "/graphs"))

# Read metadata table
sample_metadata <- read.table(args[1], sep=',', header=TRUE, stringsAsFactors = FALSE)


# Ad $files columns with path to the quant.sf files (result from salmon)
sample_metadata$files <- paste0(directory,"/processed/", sample_metadata$names, "/salmon/quant.sf")
head(sample_metadata)

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
se <- tximeta(sample_metadata)

# summarize transcript counts to genes
gse <- summarizeToGene(se)
row.names(gse) = substr(rownames(gse), 1, 15)

# Construction of DESeqDataSet object
dds <- DESeqDataSet(gse, design = ~ patient + specimen)
dds$specimen <- relevel(dds$specimen, ref = "NDMM")

# Removes IG genes
file_path <- args[3]
lines <- readLines(file_path)
cleaned_lines <- gsub('"', '', lines)
gene_vector <- unlist(strsplit(cleaned_lines, ", "))
gene_vector <- as.character(gene_vector)
exclude <- which(rownames(dds) %in% gene_vector)
dds <- dds[-exclude, ]

# Remove rows of the DESeqDataSet that have no counts, or only a single count across all samples
keep <- rowSums(counts(dds)) > 1
dds <- dds[keep,]
head(dds)
save(dds, file = paste0(directory, "/results/", args[2], "/dds.RData"))


# # EXPLORATORY ANALYSIS AND VISUALIZATION
# ########################################
# # counts normalization
 rld <- rlog(dds, blind = FALSE)
 save(rld, file = paste0(directory, "/results/", args[2], "/rld.RData"))
# # principal component analysis

pcaData <- plotPCA(rld, intgroup = c( "specimen", "patient"), returnData = TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))

pdf(paste0(directory, "/results/", args[2], "/graphs/pca_rld.pdf"))
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

pdf(paste0(directory, "/results/", args[2], "/graphs/dispersion.pdf"))
plotDispEsts(dds)
dev.off()

pdf(paste0(directory, "/results/", args[2], "/graphs/cooks.pdf"))
par(mar=c(8,5,2,2))
boxplot(log10(assays(dds)[["cooks"]]), range=0, las=2)
dev.off()

counts_ = counts(dds, normalized=T)
write.csv(counts_, file=paste0(directory, "/results/", args[2], "/graphs/norm_counts.csv"))

# Choose your contrast. In case of first example it takes change in MRD expression in comparison to NDMM
res <- results(dds, contrast=c("specimen","MRD","NDMM"))
folder_name <- "NDMM_VS_MRD"
#res <- results(dds, contrast=c("specimen","RRMM","NDMM"))
#  folder_name <- "NDMM_VS_RRMM"
# res <- results(dds, contrast=c("specimen","RRMM","MRD"))
#    folder_name <- "MRD_VS_RRMM"

#check that specific result sub-folder exists:
if (!dir.exists(paste0(directory, "/results/", args[2], "/graphs/", folder_name))) {
  dir.create(paste0(directory, "/results/", args[2], "/graphs/", folder_name), recursive = TRUE)
  cat("Folder created", folder_name, "\n")
} else {
  cat("Folder already exists:", folder_name, "\n")
}

res = subset(res, res$baseMean > 10)

mart <- useMart("ENSEMBL_MART_ENSEMBL")
mart <- useDataset("hsapiens_gene_ensembl", mart)

annotLookup <- getBM(
  mart=mart,
  attributes=c("ensembl_gene_id", "gene_biotype", "external_gene_name", "description"),
  filter="ensembl_gene_id",
  values=rownames(res),
  uniqueRows=TRUE)

############################
# genes with no annotation
problematic_ids = setdiff(rownames(res), annotLookup$ensembl_gene_id)

fileConn <- file(paste0(directory, "/results/", args[2], "/graphs/", folder_name, "/problematic_ids.out"))
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
write.csv(resOrderedDF, file = paste0(directory, "/results/", args[2], "/graphs/", folder_name, "/deseq_result_paired.csv"))


# plot heatmap with all significant genes
pdf(paste0(directory, "/results/", args[2], "/graphs/", folder_name, "/significant_genes_heatmap.pdf"))
resSig <- subset(res, padj < 0.1)
resSig <- subset(resSig, biotype != "IG_V_gene")
resSig <- subset(resSig, biotype != "IG_C_gene")
# resSig <- subset(resSig, biotype == "protein_coding")
resSig <- subset(resSig, abs(log2FoldChange) > 1)
select <- rownames(resSig[ order(resSig$padj), ], )
mat = assay(rld)[select,]
mat = t(scale(t(mat)))
mat = head(mat, 50)
write.csv(mat, file=paste0(directory, "/results/", args[2], "/heatmap_matrix.csv"))
anno <- as.data.frame(colData(rld)["specimen"])
pheatmap(mat, cluster_rows=TRUE, show_rownames=TRUE, 
         cluster_cols=TRUE, annotation_col=anno, show_colnames=TRUE)
dev.off()



#file_path <- args[4]
#lines <- readLines(file_path)
#cleaned_lines <- gsub('"', '', lines)
#neutro_genes <- unlist(strsplit(cleaned_lines, ", "))
#neutro_genes <- as.character(neutro_genes)

#heatmap with geneids
pdf(paste0(directory, "/results/", args[2], "/graphs/", folder_name, "/significant_genes_heatmap_geneid.pdf"))
resSig <- subset(res, padj < 0.05)
resSig <- subset(resSig, biotype != "IG_V_gene")
resSig <- subset(resSig, biotype != "IG_C_gene")
resSig <- subset(resSig, biotype == "protein_coding")
resSig <- subset(resSig, abs(log2FoldChange) > 1)
#resSig <- resSig[!(resSig$symbol %in% neutro_genes), ]
select <- rownames(resSig[ order(resSig$padj), ], )
reorder <- resSig[order(resSig$padj), ]
mat = assay(rld)[select,]
row.names(mat) <- reorder$symbol
row.names(mat) <- ifelse(row.names(mat)=="", rownames(reorder), row.names(mat))
mat = t(scale(t(mat)))
# mat = head(mat, 50)
anno <- as.data.frame(colData(rld)["specimen"])
pheatmap(mat, cluster_rows=TRUE, show_rownames=TRUE, cluster_cols=FALSE, annotation_col=anno, show_colnames=TRUE, fontsize = 8, cellwidth = 8, treeheight_row = 250)
dev.off()



pdf(paste0(directory, "/results/", args[2], "/graphs/", folder_name, "/ma_plot.pdf"))
plotMA(res)
dev.off()

summary(res)

