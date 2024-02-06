
maf_files <- list.files(path = "/scratch/project/open-27-18/emd_trans/processed/",
                        pattern = "\\.maf$",
                        full.names = TRUE,
                        recursive = TRUE)
                        
mymaf = maftools::merge_mafs(mafs = maf_files, verbose = TRUE)
maftools::write.mafSummary(maf = mymaf, basename = "mymaf")