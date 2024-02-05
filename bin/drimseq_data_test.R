args = commandArgs(trailingOnly=TRUE)


#parser <- ArgumentParser(description = "")
#parser$add_argument("--baseDir", type= "character")
#parser$add_argument("--input1", type= "character")
#parser$add_argument("--out_path", type= "character")
#parser$add_argument("--input2", type= "character")
#parser$add_argument("--alpha", type= "double")
#parser$add_argument("--mge", type= "double")
#parser$add_argument("--mfe", type= "double")
#parser$add_argument("--mfp", type= "double")

#xargs <- parser$parse_args()


suppressWarnings(dir.create(paste(args[1], "results", args[3], "DRIMSeq/plots", sep="/")))

# read sample information from csv file
samps <- read.table(args[2], sep=",", header=TRUE)

samps$files <- paste0(args[1], "processed/", samps$names, "/salmon/quant.sf")

# convert the "specimen" column to a factor
samps$specimen <- factor(samps$specimen)

#sanitization
samps$names <- gsub("-", "", samps$names)
samps$names <- gsub("_", "", samps$names)

head(samps)