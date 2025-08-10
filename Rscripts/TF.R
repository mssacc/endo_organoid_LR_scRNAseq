##################
###DoRothEA TF Script

##Load required libraries
library(dorothea)
library(dplyr)
library(ggplot2)
library(igraph)
library(ggraph)
library(tidygraph)
library(tibble)
library(pheatmap)
library(tidyr)
library(RColorBrewer)
library(patchwork)

data(dorothea_hs, package = "dorothea")
net <- dorothea_hs                   
## Filter rows where confidence is A or B
net = dorothea_hs %>%
  filter(confidence %in% c("A","B","C"))

net <- net[,c('tf','target','mor','confidence')]
colnames(net) <- c('source','target','mor','confidence')


#Load in DE analysis
# Load each CSV and add a cell_type column
pre_unciliated_DE <- read.csv("/data/gpfs/projects/punim1901/flames_v2/seurat_workspace/DEGs/DE_analysis_Pre_Unciliated.csv") %>%
  mutate(cell_type = "Pre-Unciliated")
unciliated_DE <- read.csv("/data/gpfs/projects/punim1901/flames_v2/seurat_workspace/DEGs/DE_analysis_Unciliated.csv") %>%
  mutate(cell_type = "Unciliated")
ciliated_DE <- read.csv("/data/gpfs/projects/punim1901/flames_v2/seurat_workspace/DEGs/DE_analysis_Ciliated.csv") %>%
  mutate(cell_type = "Ciliated")
secretory_DE <- read.csv("/data/gpfs/projects/punim1901/flames_v2/seurat_workspace/DEGs/DE_analysis_Secretory.csv") %>%
  mutate(cell_type = "Secretory")
pre_ciliated_DE <- read.csv("/data/gpfs/projects/punim1901/flames_v2/seurat_workspace/DEGs/DE_analysis_Pre_Ciliated.csv") %>%
  mutate(cell_type = "Pre-Ciliated")
proliferative_DE <- read.csv("/data/gpfs/projects/punim1901/flames_v2/seurat_workspace/DEGs/DE_analysis_Proliferative.csv") %>%
  mutate(cell_type = "Proliferative")

#Combine all into one data frame
combined_DE <- bind_rows(pre_unciliated_DE, unciliated_DE, ciliated_DE,
                         secretory_DE, pre_ciliated_DE, proliferative_DE)


###Print out DEGs that are also TFs
  results <- combined_DE
  results$sig <- ifelse(abs(results$log2FoldChange) >=0.585 & results$padj<0.05,TRUE,FALSE)
  
  #Pull out DEGs that are TFs
  DEGTFs <- results %>%
    filter(gene %in% net$source & sig)
  
  #Pull out DEGs that are targets
  DEG_targets <- results %>%
    filter(gene %in% net$target & sig)



###Create combined dataframe of common TF and TARGET DEGs per cell subtype
results_2 <- combined_DE
results_2 <- results_2 %>%
  mutate(sig = ifelse(abs(log2FoldChange) >= 0.585 & padj < 0.05, TRUE, FALSE))

#Get TFs that are DE and in the DoRothEA regulon
DE_TFs <- results_2 %>%
  filter(gene %in% net$source, sig) %>%
  select(gene, cell_type, log2FoldChange, padj)

# Join regulon with DE TF info
regulon_with_DE_TFs <- net %>%
  inner_join(DE_TFs, by = c("source" = "gene")) %>%
  rename(TF = source, TF_log2FC = log2FoldChange, TF_padj = padj, TF_cell_type = cell_type)

# Now join to find if the target gene is also a DEG in the same cell type
TF_target_DEG_pairs <- regulon_with_DE_TFs %>%
      inner_join(results_2 %>% filter(sig) %>%
      select(gene, log2FoldChange, padj, cell_type),
             by = c("target" = "gene", "TF_cell_type" = "cell_type")) %>%
      rename(target_log2FC = log2FoldChange,
             target_padj = padj,
             cell_type = TF_cell_type)


###Generate some plots
# Count DE TFs per cell type
tf_counts <- DEGTFs %>%
  group_by(cell_type) %>%
  summarise(n_TFs = n_distinct(gene))

ggplot(tf_counts, aes(x = reorder(cell_type, -n_TFs), y = n_TFs, fill = cell_type)) +
  geom_bar(stat = "identity", show.legend = FALSE) +
  theme_minimal(base_size = 14) +
  labs(title = "Number of Differentially Expressed TFs per Cell Type",
       x = "Cell Type", y = "Number of DE TFs") +
  coord_flip()


###Heatmap of TFs
#Create a wide matrix of log2FC per TF per cell type
tf_matrix <- DEGTFs %>%
  select(gene, cell_type, log2FoldChange) %>%
  pivot_wider(names_from = cell_type, values_from = log2FoldChange, values_fill = 0) %>%
  column_to_rownames("gene")

desired_order <- c("Pre-Unciliated", "Unciliated", "Ciliated", "Secretory", "Pre-Ciliated", "Proliferative")
tf_matrix <- tf_matrix[, intersect(desired_order, colnames(tf_matrix))]

#Fine-grained continuous color breaks for smooth gradient
breaks_continuous <- seq(-2, 2, by = 0.1)
#Tick marks every 1 unit on the legend
legend_ticks <- seq(-2, 2, by = 1)
#Smooth red-yellow-blue palette (similar to pheatmap default)
my_colors <- colorRampPalette(rev(brewer.pal(n = 11, name = "RdYlBu")))(length(breaks_continuous) - 1)

pheatmap(as.matrix(tf_matrix),
         cluster_rows = TRUE, cluster_cols = FALSE,
         main = "Log2FC of Differentially Expressed TFs",
         fontsize_row = 10, fontsize_col = 11, angle_col = 0,
         treeheight_row = 0, treeheight_col = 0,
         breaks = breaks_continuous, legend_breaks = legend_ticks,
         color = my_colors)

####All log2FC graphed
#List of TFs to plot
target_TFs <- c("POU2F2", "EGR1", "POU5F1", "BCL11A", "MEIS2", "FOXM1", "SOX9", "TFAP2A")

# Create matrix for selected TFs across all cell types
tf_matrix <- combined_DE %>%
  filter(gene %in% target_TFs) %>%
  select(gene, cell_type, log2FoldChange) %>%
  pivot_wider(names_from = cell_type, values_from = log2FoldChange, values_fill = 0) %>%
  column_to_rownames("gene")

desired_order <- c("Pre-Unciliated", "Unciliated", "Ciliated", "Secretory", "Pre-Ciliated", "Proliferative")
tf_matrix <- tf_matrix[, intersect(desired_order, colnames(tf_matrix))]

pheatmap(as.matrix(tf_matrix),
         cluster_rows = TRUE, cluster_cols = FALSE,
         main = "Log2FC of Differentially Expressed TFs",
         fontsize_row = 10, fontsize_col = 11, angle_col = 0,
         treeheight_row = 0, treeheight_col = 0,
         breaks = breaks_continuous, legend_breaks = legend_ticks,
         color = my_colors)


###Network graph
# Function to generate ggraph plot for a given cell type
plot_tf_network <- function(cell_type_name) {
  network_df <- TF_target_DEG_pairs %>%
    filter(cell_type == cell_type_name) %>%
    select(TF, target, mor)
  if (nrow(network_df) == 0) return(NULL)  # Skip empty plots
  colnames(network_df) <- c("from", "to", "mor")
  graph <- tbl_graph(edges = network_df, directed = TRUE)
  graph <- graph %>%
    mutate(node_type = ifelse(name %in% network_df$from, "TF", "Target"))
  ggraph(graph, layout = "fr") +
    geom_edge_link(aes(color = factor(mor)),
                   arrow = arrow(length = unit(3, 'mm')),
                   end_cap = circle(3, 'mm'),
                   width = 1.2) +
    geom_node_point(aes(fill = node_type), size = 5, shape = 21, color = "black") +
    geom_node_text(aes(label = name), repel = TRUE, size = 3.5) +
    scale_edge_color_manual(values = c("-1" = "#F8766D", "1" = "#849AFF", "0" = "gray"),
                            labels = c("-1" = "Repression", "1" = "Activation", "0" = "Neutral"),
                            name = "Regulation") +
    scale_fill_manual(values = c("TF" = "#0CB702", "Target" = "#FF61CC"),
                      name = "Node Type") +
    theme_graph() +
    theme(plot.title = element_text(hjust = 0.5)) +
    labs(title = paste("TF-target network in", cell_type_name, "cells"))}

p1 <- plot_tf_network("Pre-Unciliated")
p2 <- plot_tf_network("Unciliated")
p3 <- plot_tf_network("Ciliated")
p4 <- plot_tf_network("Secretory")

# Combine and unify legends
(p1 + p2 + p3 + p4) +  plot_layout(guides = "collect") &
  theme(legend.position = "bottom",
        plot.title = element_text(size = 11, hjust = 0.5))



####################
#Shows expression patterns
# Determine global log2FC limits
global_log2fc_vals <- TF_target_DEG_pairs %>%
  select(TF_log2FC, target_log2FC) %>%
  unlist() %>%
  as.numeric()

log2fc_limit <- max(abs(global_log2fc_vals), na.rm = TRUE)  # symmetric around 0
log2fc_breaks <- seq(-log2fc_limit, log2fc_limit, by = 1)

# Updated network plot function with fixed scale
plot_tf_network <- function(cell_type_name) {
  network_df <- TF_target_DEG_pairs %>%
    filter(cell_type == cell_type_name) %>%
    select(TF, target, mor, TF_log2FC, target_log2FC)
  
  if (nrow(network_df) == 0) return(NULL)
  tf_nodes <- network_df %>%
    select(name = TF, log2FC = TF_log2FC) %>%
    mutate(node_type = "TF")
  target_nodes <- network_df %>%
    select(name = target, log2FC = target_log2FC) %>%
    mutate(node_type = "Target")
  node_attrs <- bind_rows(tf_nodes, target_nodes) %>%
    group_by(name) %>%
    slice(1) %>%
    ungroup()
  edges <- network_df %>%
    mutate(mor = factor(mor, levels = c("-1", "1"))) %>%
    select(from = TF, to = target, mor)
  graph <- tbl_graph(nodes = node_attrs, edges = edges, directed = TRUE)
  
  ggraph(graph, layout = "fr") +
    geom_edge_link(aes(color = factor(mor)),
                   arrow = arrow(length = unit(3, 'mm')),
                   end_cap = circle(3, 'mm'), width = 1.2) +
    geom_node_point(aes(fill = log2FC), shape = 21, size = 6, color = "black") +
    geom_node_text(aes(label = name), repel = TRUE, size = 3.5) +
    scale_fill_gradientn(colours = my_colors, limits = c(-2.5, 2.5), breaks = legend_ticks,
                         labels = scales::number_format(accuracy = 1), name = "log2FC") +
    scale_edge_color_manual(values = c("-1" = "#F8766D", "1" = "palegreen3"),
      labels = c("-1" = "Repression", "1" = "Activation"),
      name = "Regulation",
      drop = FALSE) +
    theme_graph() +
    labs(title = paste("TF-target network in", cell_type_name, "cells")) +
    theme(plot.title = element_text(hjust = 0.5))}

# Generate and combine plots with shared legend
p1 <- plot_tf_network("Pre-Unciliated")
p2 <- plot_tf_network("Unciliated")
p3 <- plot_tf_network("Ciliated")
p4 <- plot_tf_network("Secretory")

(p1 + p2 + p3 + p4) + 
  plot_layout(guides = "collect") &
  theme(legend.position = "bottom",
        plot.title = element_text(size = 11, hjust = 0.5))

