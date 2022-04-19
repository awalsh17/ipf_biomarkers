# Functions


plot_heatmap <- function(input_df,
                         input_genes) {
  ### For the median data
  mat_data <- input_df %>%
    dplyr::filter(description %in% input_genes) %>%
    # select only median rows
    dplyr::select(adipose_subcutaneous:whole_blood) %>%
    dplyr::select(-cells_cultured_fibroblasts, -cells_ebv_transformed_lymphocytes) %>%
    as.matrix()

  mat_annotation <- input_df %>%
    dplyr::filter(description %in% input_genes) %>%
    # select only annotation rows
    dplyr::select(ensembl_gene_id, description)

  column_groups <- c("adipose",
                 "adipose",
                 "adrenal_gland",
                 "artery",
                 "artery",
                 "artery",
                 "bladder",
                 "brain",
                 "brain",
                 "brain",
                 "brain",
                 "brain",
                 "brain",
                 "brain",
                 "brain",
                 "brain",
                 "brain",
                 "brain",
                 "brain",
                 "brain",
                 "breast",
                 # "cells_cultured_fibroblasts",
                 # "cells_ebv_transformed_lymphocytes",
                 "cervix",
                 "cervix",
                 "colon",
                 "colon",
                 "esophagus",
                 "esophagus",
                 "esophagus",
                 "fallopian_tube",
                 "heart",
                 "heart",
                 "kidney",
                 "kidney",
                 "liver",
                 "lung",
                 "minor_salivary_gland",
                 "muscle_skeletal",
                 "nerve_tibial",
                 "ovary",
                 "pancreas",
                 "pituitary",
                 "prostate",
                 "skin",
                 "skin",
                 "small_intestine_terminal_ileum",
                 "spleen",
                 "stomach",
                 "testis",
                 "thyroid",
                 "uterus",
                 "vagina",
                 "whole_blood")

  mat_col_annotation <- HeatmapAnnotation(tissue = anno_simple(column_groups,
                                         col = c("adipose" = "gray",
                                                 "adipose" = "gray",
                                                 "adrenal_gland" = "gray",
                                                 "artery" = "gray",
                                                 "artery" = "gray",
                                                 "artery" = "gray",
                                                 "bladder" = "gray",
                                                 "brain" = "gray",
                                                 "brain" = "gray",
                                                 "brain" = "gray",
                                                 "brain" = "gray",
                                                 "brain" = "gray",
                                                 "brain" = "gray",
                                                 "brain" = "gray",
                                                 "brain" = "gray",
                                                 "brain" = "gray",
                                                 "brain" = "gray",
                                                 "brain" = "gray",
                                                 "brain" = "gray",
                                                 "brain" = "gray",
                                                 "breast" = "gray",
                                                 # "cells_cultured_fibroblasts" = "gray",
                                                 # "cells_ebv_transformed_lymphocytes" = "gray",
                                                 "cervix" = "gray",
                                                 "cervix" = "gray",
                                                 "colon" = "gray",
                                                 "colon" = "gray",
                                                 "esophagus" = "gray",
                                                 "esophagus" = "gray",
                                                 "esophagus" = "gray",
                                                 "fallopian_tube" = "gray",
                                                 "heart" = "gray",
                                                 "heart" = "gray",
                                                 "kidney" = "gray",
                                                 "kidney" = "gray",
                                                 "liver" = "gray",
                                                 "lung" = "dodgerblue3",
                                                 "minor_salivary_gland" = "gray",
                                                 "muscle_skeletal" = "gray",
                                                 "nerve_tibial" = "gray",
                                                 "ovary" = "gray",
                                                 "pancreas" = "gray",
                                                 "pituitary" = "gray",
                                                 "prostate" = "gray",
                                                 "skin" = "gray",
                                                 "skin" = "gray",
                                                 "small_intestine_terminal_ileum" = "gray",
                                                 "spleen" = "gray",
                                                 "stomach" = "gray",
                                                 "testis" = "gray",
                                                 "thyroid" = "gray",
                                                 "uterus" = "gray",
                                                 "vagina" = "gray",
                                                 "whole_blood" = "gray"
                                         )))
  # create groups based on tissue to separate

  col_split <- c(rep("other", 34),
                 "lung",
                 rep("other", 17)
                 )


  rownames(mat_data) <- mat_annotation$description
  map_colors <- c("white","dodgerblue3")

  lgd <- list(title = "median log2(TPM+1)",
               title_position = "lefttop-rot",
               legend_height = unit(4, "cm"))

  Heatmap(log2(mat_data+1),
          col = map_colors,
          column_split = col_split,
          top_annotation = mat_col_annotation,
          clustering_method_rows = "ward.D2",
          clustering_method_columns = "ward.D2",
          show_column_names = FALSE,
          show_row_names = TRUE,
          heatmap_legend_param = lgd)
}

#' Return omnipathdb results
get_omnipath_receptors <- function(gene_list) {
  read.delim(
    paste0("https://omnipathdb.org/interactions/?genesymbols=1&partners=",
           paste(gene_list, collapse = ","),
           "&datasets=ligrecextra&fields=sources,references"
    )
  ) %>%
    mutate(n_sources = stringr::str_count(string = sources, pattern = ";") + 1)
}

#' Make dotplot for GTEx
plot_gtex_dots <- function(input_data, input_genes) {
  plot_data <- input_data %>%
    dplyr::filter(description %in% input_genes) %>%
    # select only median rows
    dplyr::select(gene = description, adipose_subcutaneous:whole_blood) %>%
    dplyr::select(-cells_cultured_fibroblasts,
                  -cells_ebv_transformed_lymphocytes) %>%
    # pivot to long
    tidyr::pivot_longer(cols = -gene)

  myplot <- plot_data %>%
    group_by(gene) %>%
    mutate(top_tissue = ifelse(value == max(value), name, NA)) %>%
    ungroup() %>%
    mutate(Tissue = ifelse(name == "lung","Lung","Other"),
           size = ifelse(name == "lung", 2, 1)) %>%
    ggplot(aes(x = gene, y = log2(value + 1),
               color = Tissue, size = size, shape = Tissue)) +
    geom_jitter(width = 0.2, height = 0) +
    scale_color_manual(values = c("dodgerblue3","gray")) +
    scale_shape_manual(values = c("Lung" = 15, "Other" = 1)) +
    scale_size(range = c(1,2), guide = "none") +
    # geom_text(aes(label = top_tissue), size = 2, color = "black") +
    geom_hline(yintercept = log2(2), color = "seagreen",
               size = 0.5, linetype = "dashed") +
    labs(x = NULL, y = expression(log[2](TPM + 1))) +
    theme_minimal(base_family = "Avenir") +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5),
          legend.position = "top",
          panel.grid.major.x = element_blank() )
  return(myplot)
}

#' Make a plot of all significant diff expr
plot_sc_diff <- function(input_data) {
  ## Using the Haberman differential expression
  input_data %>%
    mutate(cell_group = case_when(
      cluster %in% c("cDCs","Macrophages","Mast Cells","T Cells",
                     "NK Cells", "Proliferating Macrophages") ~ "immune",
      cluster %in% c("Fibroblasts","Myofibroblasts",
                     "Smooth Muscle Cells") ~ "fibroblasts",
      cluster %in% c("Lymphatic Endothelial Cells",
                     "Endothelial Cells") ~ "endothelial",
      TRUE ~ "epithelial"
    )) %>%
    ggplot(aes(x = gene, y = avg_diff)) +
    # geom_point(aes(color = cell_group, size = log)) +
    geom_jitter(aes(color = cell_group, size = log),
                height = 0, width = 0.1, shape = 1) +
    # geom_text(aes(label = cluster), size = 2, angle = 45) +
    scale_size(range = c(1,3)) +
    coord_flip() +
    labs(x = NULL,
         y = "Difference IPF vs control",
         size = bquote("-"~log[10]*"(adj. p)"),
         color = "Cell type"
    ) +
    theme_minimal(base_family = "Avenir") +
    theme(legend.position = "right",
          panel.grid.major.y = element_blank() )

}
