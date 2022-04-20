# Need to run the analysis using Haberman scRNAseq data

library(Seurat) # see renv.lock for version information
library(MAST)

# data is in the /data directory on my local machine.
# data was downloaded from GEO (GSE135893)

metadata <- readr::read_csv("data/GSE135893_IPF_metadata.csv.gz")
colnames(metadata)[1] <- "barcode"
annotated <- readRDS("data/GSE135893_ILD_annotated_fullsize.rds")

levels(annotated) # what are listed

# Run differential expression for all cell clusters IPF v control
# head(FindMarkers(pbmc, ident.1 = "CD14+ Mono", ident.2 = "FCGR3A+ Mono", test.use = "MAST"))

# immune.combined$celltype.stim <- paste(Idents(immune.combined), immune.combined$stim, sep = "_")
# immune.combined$celltype <- Idents(immune.combined)
# Idents(immune.combined) <- "celltype.stim"
# b.interferon.response <- FindMarkers(immune.combined, ident.1 = "B_STIM", ident.2 = "B_CTRL", verbose = FALSE)

# Combine all diff. expr. into a big table. Write out!


# Make dot plots for genes of interest (all samples)

