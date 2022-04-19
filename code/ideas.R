
library(dplyr)
library(ggplot2)

elevated_genes <- c("SFTPD",
                    "SFTPA1",
                    "MMP7",
                    "CCL22",
                    "TIMP1",
                    "CXCL13",
                    "CXCL10",
                    "CXCL12",
                    "CHI3L1")

## Using the Haberman differential expression
diff_expr_table_ligands %>%
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
