# Statistical analysis of differential transcript usage following Salmon quantification
# https://bioconductor.riken.jp/packages/3.10/workflows/vignettes/rnaseqDTU/inst/doc/rnaseqDTU.html


#README
#The names of the samples shouldn't contain a "-" symbol. Although I implemented
#a fix, I can't yet guarantee it works 100% of a time. The script need a little
#tinkering each time it is used on new samples -> row 131 change the name of the
#condition.

library("tximport")
library("DRIMSeq")
library("GenomicFeatures")
library("stageR")
library("ggplot2")

#    ${dir} arg1
#    ${csv} arg2
#    ${gtf} arg3
#    ${run_name} arg4 
#    ${alpha} arg5
#    ${mge} arg6
#    ${mfe} arg7
#    ${mfp} arg8



args = commandArgs(trailingOnly=TRUE)

suppressWarnings(dir.create(paste(args[1], "results", args[4], "DRIMSeq/plots", sep="/")))

# read sample information from csv file
samps <- read.table(args[2], sep=",", header=TRUE)


#samps$files <- paste0(args[1], "processed/", samps$names, "/salmon/quant.sf")

# convert the "specimen" column to a factor
#samps$condition <- factor(samps$specimen)


# create file paths using file.path
files <- paste0(args[1], "processed/", samps$names, "/salmon/quant.sf")


# set file names as sample ids
names(samps)[names(samps) == "names"] <- "sample_id"
names(samps)[names(samps) == "specimen"] <- "condition"
samps$condition <- factor(samps$condition)
names(files) <- samps$sample_id

files

#sanitization
#samps$sample_id <- gsub("-", "", samps$sample_id)
#samps$sample_id <- gsub("_", "", samps$sample_id)
head(samps)

# import transcript quantification data with tximport
txi <- tximport(files, type="salmon", txOut=TRUE,
    countsFromAbundance="scaledTPM")

# extract transcript counts
cts <- txi$counts
cts <- cts[rowSums(cts) > 0,]

# create transcript to gene mapping
txdb <- makeTxDbFromGFF(args[3], format="auto")
txdf <- select(txdb, keys(txdb, "GENEID"), "TXNAME", "GENEID")  
# count the number of transcripts per gene
tab <- table(txdf$GENEID)
txdf$ntx <- tab[match(txdf$GENEID, names(tab))]

cts[1:3,1:3]
head(txdf)

# check if all quantified transcripts are in the database
all(rownames(cts) %in% txdf$TXNAME)

# match quantified transcripts to gene ids
txdf <- txdf[match(rownames(cts),txdf$TXNAME),]

# check if matching is successful
txdf <- txdf[match(rownames(cts),txdf$TXNAME),]
all(rownames(cts) == txdf$TXNAME) 


print("samples_ID")
samps$sample_id


# create a data frame with gene ids, transcript ids and counts
counts <- data.frame(gene_id=txdf$GENEID,
    feature_id=txdf$TXNAME, cts, check.names = TRUE)

colnames(counts)
samps$sample_id <- gsub("-", ".", samps$sample_id)
print("colnames")
colnames(counts)
# create a DTU object from the counts and sample information
d <- dmDSdata(counts=counts, samples=samps)
d
counts(d[1,])[,1:4]
#n.small <- min(aggregate(sample_id ~ condition, samps,
#    function(x) length(unique(x)))$sample_id)
# filter out low-expression features based on minimum number of samples and expression levels
n <- nrow(samps)
n.small <- min(table(samps$condition))
d <- dmFilter(d,
                min_samps_feature_expr=n.small,
                min_feature_expr=args[7],
                min_samps_feature_prop=n.small,
                min_feature_prop=args[8],
                min_samps_gene_expr=n,
                min_gene_expr=args[6])

# bin gene ids based on the number of transcripts
table(table(counts(d)$gene_id))

# create the design matrix for the full model
design_full <- model.matrix(~condition, data=DRIMSeq::samples(d))
colnames(design_full)
#for testing purposes for less computations
d <- d[1:500,]

# Calculate the precision for d using the full design
set.seed(1)
system.time({
  d <- dmPrecision(d, design=design_full)
  d <- dmFit(d, design=design_full)
  d <- dmTest(d, coef = "conditionNDMM")
})


# Extract the results of the DRIMSeq analysis into res
res <- DRIMSeq::results(d)
head(res)
# Extract the results of the DRIMSeq analysis at the transcript level into res.txp
res.txp <- DRIMSeq::results(d, level="feature")
head(res.txp)
# Define the no.na function to replace NAs with 1
no.na <- function(x) ifelse(is.na(x), 1, x)

# Replace NAs in res$pvalue with 1
res$pvalue <- no.na(res$pvalue)

# Replace NAs in res.txp$pvalue with 1
res.txp$pvalue <- no.na(res.txp$pvalue)
print("1")
# Write res as a table to the file specified by xargs$output_TSV1 with "drimseq_genes_DTU.tsv" as suffix
save_location <- paste0(args[1], "results", args[4], "/DRIMSeq/genes_DTU.tsv")
print(save_location)

write.table(res,
    file = paste0(args[1], "results", args[4], "/DRIMSeq/genes_DTU.tsv"),
    row.names=FALSE, na="",col.names=TRUE, sep="\t")
print("2")
# Write res.txp as a table to the file specified by xargs$output_TSV2 with "drimseq_transcripts_proportions.tsv" as suffix
write.table(res.txp,
    file = paste0(args[1], "results", args[4], "/DRIMSeq/transcripts_proportions.tsv"),
    row.names=FALSE, na="",col.names=TRUE, sep="\t")
print("3")
#save(d, file = paste(args[1], "results", args[4], "/DRIMSeq/d.RData"))
print("4")
# Find the indices of the significant results with adjusted p-value less than xargs$alpha
signif_idx <- which(res$adj_pvalue < args[5])

head(rownames(cts))
print("txdfTXNAME")
head(txdf$TXNAME)

#idx <- which(res$adj_pvalue < 0.05)[1]
#res[idx,]
#pdf(paste0(args[1],"results/DRIMSeq/plots/test.pdf"))
#    plotProportions(d, res$gene_id[idx], "condition")
#dev.off()
#Plot the proportion for each significant result and save the plot as pdf
for (i in 1:length(signif_idx)){
    idx = signif_idx[i]
    pdf(paste0(args[1], "results", args[4], "/DRIMSeq/plots/",paste(res$gene_id[idx],"pdf",sep=".")))
    print(plotProportions(d, res$gene_id[which(res$adj_pvalue < args[5])[i]], group_variable="condition", plot_type = "boxplot1"))
    dev.off()
}

# Prepare the data for StageR analysis
pScreen <- res$pvalue

# Define a function to shorten the gene/transcript names
strp <- function(x) substr(x,1,15)

# Shorten the gene names for pScreen
names(pScreen) <- strp(res$gene_id)

# Convert res.txp$pvalue into a matrix pConfirmation
pConfirmation <- matrix(res.txp$pvalue, ncol=1)

# Set row names of "pConfirmation" to cleaned "feature_id" from "res.txp"
rownames(pConfirmation) <- strp(res.txp$feature_id)

# Create data frame "tx2gene" with only "feature_id" and "gene_id" from "res.txp"
tx2gene <- res.txp[,c("feature_id", "gene_id")]

# Clean "feature_id" and "gene_id" in "tx2gene"
for (i in 1:2) tx2gene[,i] <- strp(tx2gene[,i])

# Create "stageRObj" using "stageRTx" on "pScreen", "pConfirmation", and "tx2gene"
# Adjust "stageRObj" with "stageWiseAdjustment" and "dtu" method, using "stageR_alpha" from "xargs"
# Get adjusted p-values from "stageRObj" and store in "drim.padj", suppressing warnings
stageRObj <- stageRTx(pScreen=pScreen, pConfirmation=pConfirmation,
                         pScreenAdjusted=FALSE, tx2gene=tx2gene)
stageRObj <- stageWiseAdjustment(stageRObj, method="dtu", alpha=args[5])
suppressWarnings({
  drim.padj <- getAdjustedPValues(stageRObj, order=FALSE,
                                      onlySignificantGenes=FALSE)
})

# Write "drim.padj" to "StageR_DRIMSeq.tsv"
write.table(drim.padj,
    file = paste0(args[1], "results", args[4], "/DRIMSeq/StageR.tsv"),
    row.names=FALSE, na="",col.names=TRUE, sep="\t")