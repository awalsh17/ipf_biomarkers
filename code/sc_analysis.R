# Need to run the analysis using Haberman scRNAseq data

library(Seurat) # see renv.lock for version information
library(MAST)
library(dplyr)

# data is in the /data directory on my local machine.
# data was downloaded from GEO (GSE135893)

metadata <- readr::read_csv("data/GSE135893_IPF_metadata.csv.gz")
colnames(metadata)[1] <- "barcode"
# this file is enormous
annotated <- readRDS("data/GSE135893_ILD_annotated_fullsize.rds")

levels(annotated) # what cell types - there are 31
unique(annotated$Diagnosis) # "IPF", "Control", "NSIP", "cHP", "Unclassifiable ILD", "sacroidosis"

annotated$celltype <- Idents(annotated)
# Run differential expression for all cell clusters IPF v control
Idents(annotated) <- "Diagnosis"

# # subset to IPF and Control - not enough memory? crash
# annotated <- subset(annotated, idents =  c("IPF","Control"))

# create new group
annotated$cell.disease <- paste(Idents(annotated), annotated$celltype, sep = "_")
Idents(annotated) <- "cell.disease"

# Test in all cell types
all_comparisons <- 
  tidyr::crossing("disease" = unique(annotated$Diagnosis),
                  "cells" = unique(annotated$celltype)) %>% 
  dplyr::filter(disease %in% c("IPF","Control")) %>% 
  mutate(ident = paste(disease, cells, sep = "_")) %>% 
  tidyr::pivot_wider(names_from = disease, values_from = ident) %>% 
  filter(cells != "HAS1 High Fibroblasts") # not enough

# FindMarkers returns data.frames where genes are the row.names
# should probably move so the gene names are a column
# test <- FindMarkers(annotated, ident.1 = "Control_AT1", ident.2 = "IPF_AT1", test.use = "MAST")

all_diff <- purrr::map2(all_comparisons$Control, all_comparisons$IPF,
                        ~FindMarkers(object = annotated, ident.1 = .x,
                                     ident.2 = .y, test.use = "MAST", verbose = F))

# Combine all diff. expr. into a big table. Write out!
all_diff <- lapply(all_diff, function(x) tibble::rownames_to_column(x, var = "gene"))
names(all_diff) <- all_comparisons$cells
all_diff <- bind_rows(all_diff, .id = "cell")

write.csv(all_diff, "results/results_mast_ipf_v_control.csv", row.names = F)

# Positive avg_logFC - higher in group 1 (Control) than group 2 (IPF)
# Make dot plots or violin plots for genes of interest (all samples)
mygenes <- readLines("gene_symbols.txt")

Idents(annotated) <- "Diagnosis"
plots <- VlnPlot(annotated, features = mygenes,
                 idents = c("Control","IPF"),
                 split.by = "Diagnosis",
                 group.by = "celltype", 
                 pt.size = 0, combine = FALSE)

saveRDS(plots, "results/vlnplots_measuredgenes.Rds")

# intersect results with list of interest

diff_mygenes <- all_diff %>% filter(gene %in% mygenes)

write.csv(diff_mygenes, "results/results_mast_ipf_v_control_selected.csv", row.names = F)

# Can also move files to google drive:
# library(googledrive)

# googledrive::drive_auth()

drive_upload(media = "results/vlnplots_measuredgenes.Rds",
             path = "projects/")
