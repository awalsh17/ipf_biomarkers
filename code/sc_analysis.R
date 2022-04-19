# Need to run the analysis using Haberman scRNAseq data

library(Seurat) # see renv.lock for version information

# data is in the /data directory on my local machine.
# data was downloaded from GEO (GSE135893)

metadata <- readr::read_csv("data/GSE135893_IPF_metadata.csv")
colnames(metadata)[1] <- "barcode"
annotated <- readRDS("data/GSE135893_ILD_annotated_fullsize.rds")

# Run differential expression for all cell clusters IPF v control

# Make dot plots for genes of interest (all samples)

