library(targets)
tar_option_set(packages = c(
  "dplyr",
  "ggplot2",
  "ComplexHeatmap"
))
source("code/functions.R")
elevated_genes <- c("SFTPD",
                    "SFTPA1",
                    "MMP7",
                    "CCL22",
                    "TIMP1",
                    "CXCL13",
                    "CXCL10",
                    "CXCL12",
                    "CHI3L1")

list(
  # Load our saved results
  tar_target(haberman_elevated_file,
             "../../sent_from_bms/diff_expr_table_ligands.csv",
             format = "file"),
  tar_target(haberman_receptor_file,
             "../../sent_from_bms/diff_expr_table.csv",
             format = "file"),
  tar_target(haberman_diff_elevated,
             read.csv(haberman_elevated_file)),
  tar_target(haberman_diff_receptors,
             read.csv(haberman_receptor_file)),

  # Load GTEx data

  tar_target(gtex_median_file,
           "https://storage.googleapis.com/gtex_analysis_v8/rna_seq_data/GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_median_tpm.gct.gz"),

  tar_target(gtex_medians, readr::read_delim(gtex_median_file,
                                           delim = "\t",
                                           skip = 2,
                                           show_col_types = FALSE) %>%
             janitor::clean_names() %>%
             # get rid of these PAR_Y all zero rows
             dplyr::filter(!grepl("PAR_Y", name)) %>%
             tidyr::separate(name, into=c("ensembl_gene_id","ensembl_suffix"),
                             remove=FALSE)),

  # Load gene names
  tar_target(gene_file,
             "gene_symbols.txt",
             format = "file"),
  tar_target(genes_to_plot, readLines(gene_file)),

  # Make heatmap
  tar_target(heatmap,
             plot_heatmap(input_df = gtex_medians,
                          input_genes = genes_to_plot)),
  tar_target(heatmap_up,
             plot_heatmap(input_df = gtex_medians,
                          input_genes = elevated_genes)),

  # Make GTEx jitter plot
  tar_target(gtex_plot,
             plot_gtex_dots(input_data = gtex_medians,
                            input_genes = elevated_genes)),

  # Query omnipath -
  tar_target(receptors,
             get_omnipath_receptors(elevated_genes) %>%
               filter(n_sources > 1)),

  # Make dot plots from Haberman diff expr. with elevated and receptors
  tar_target(plot_haberman_elevated,
             plot_sc_diff(input_data = haberman_diff_elevated)),
  tar_target(plot_haberman_receptors,
             plot_sc_diff(input_data = haberman_diff_receptors))
)
