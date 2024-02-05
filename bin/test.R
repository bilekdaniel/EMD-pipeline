
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


# Read metadata table
sample_metadata <- read.table(args[1], sep=',', header=TRUE, stringsAsFactors = FALSE)

# Ad $files columns with path to the quant.sf files (result from salmon)
sample_metadata$files <- paste0(args[3],"/processed/", sample_metadata$names, "/salmon/quant.sf")

# Select only fresh (no cryo preserved) sample 
filtered_samples <- subset(sample_metadata, specimen == "NDMM")

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
dds <- DESeqDataSet(gse, design = ~ 1)

# Remove rows of the DESeqDataSet that have no counts, or only a single count across all samples
keep <- rowSums(counts(dds)) > 1
dds <- dds[keep,]

save(dds, file = paste0(args[3], "results/", args[2], "/dds.RData"))

#counts_ = counts(dds, normalized=TRUE)
#write.csv(counts_, file=paste0(args[3], "/results/", args[2], "/norm_counts.csv"))



