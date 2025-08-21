###DESeq2 DEG by Cell Subtype with BEC

#Load required libraries
library(Seurat)
library(DESeq2)
library(tidyverse)
library(ggplot2)
library(limma)
library(EnhancedVolcano)
library(patchwork)
library(dplyr)
library(pheatmap)

#Load Seurat object
all_samples_integrated <- readRDS("/data/gpfs/projects/punim1901/flames_v2/seurat_workspace/all_samples_integrated.rds")

#Define your cell subtypes
cell_types <- c("Pre-Unciliated", "Unciliated", "Ciliated", "Secretory", "Pre-Ciliated", "Proliferative")

#Function to apply limma batch correction
apply_limma_batch_correction <- function(counts_matrix, batch_info) {
  log_counts <- log2(counts_matrix + 1)
  design <- model.matrix(~fertility, data = batch_info)
  corrected_counts <- removeBatchEffect(log_counts,
                                        batch = batch_info$batch,
                                        design = design)
  corrected_counts <- 2^corrected_counts - 1
  corrected_counts[corrected_counts < 0] <- 0
  return(corrected_counts)
}

#Loop through each cell subtype
for (cell in cell_types) {
  message(paste0("\nRunning analysis for: ", cell))
  
  #Subset Seurat object
  subset_obj <- subset(all_samples_integrated, subset = cell_type == cell)
  
  #Aggregate pseudobulk counts
  pseudobulk_counts <- AggregateExpression(subset_obj,
                                           group.by = c("orig.ident", "batch", "fertility"),
                                           assays = "RNA", slot = "counts")$RNA
  
  #Create metadata for the pseudobulk groups
  group_metadata <- data.frame(group_id = colnames(pseudobulk_counts)) %>%
    separate(group_id, into = c("orig.ident", "batch", "fertility"), sep = "_") %>%
    left_join(subset_obj@meta.data %>%
                select(orig.ident, batch, fertility) %>%
                distinct(),
              by = "orig.ident") %>%
    mutate(group_id = colnames(pseudobulk_counts)) %>%
    rename(batch = batch.x, fertility = fertility.x) %>%
    select(-batch.y, -fertility.y) %>%
    drop_na(batch, fertility)
  
  #Apply batch correction
  counts_corrected <- apply_limma_batch_correction(pseudobulk_counts, group_metadata)
  
  #Create DESeq2 dataset
  dds <- DESeqDataSetFromMatrix(round(counts_corrected),
                                colData = group_metadata,
                                design = ~ fertility)
  
  #Filter lowly expressed genes
  keep <- rowSums(counts(dds)) >= 10
  dds <- dds[keep, ]
  
  #Run DESeq2
  dds <- DESeq(dds)
  
  #Get DESeq2 results for Infertile vs Fertile
  res <- results(dds, contrast = c("fertility", "Infertile", "Fertile"))
  
  #Save full results (all genes) in environment
  assign(paste0("res_", gsub("-", "_", cell)), res, envir = .GlobalEnv)
  
  #Add gene column
  res <- res %>%
    as.data.frame() %>%
    rownames_to_column("gene")
  
  #Export res to CSV
  write.csv(res,
            paste0("/data/gpfs/projects/punim1901/flames_v2/seurat_workspace/DEGs/DE_analysis_", gsub("-", "_", cell), ".csv"),
            row.names = FALSE)
  
  #Extract significant DEGs
  sig_res <- res %>%
    filter(padj < 0.05, abs(log2FoldChange) > 0.585)
  
  #Save significant DEGs in environment
  assign(paste0("sig_res_", gsub("-", "_", cell)), sig_res, envir = .GlobalEnv)
  
  #Export sig_res to CSV
  write.csv(sig_res,
            paste0("/data/gpfs/projects/punim1901/flames_v2/seurat_workspace/DEGs/DEGs_", gsub("-", "_", cell), ".csv"),
            row.names = FALSE)
  
  #Create and display volcano plot for all genes
  plot <- EnhancedVolcano(res, lab = rownames(res),
                          x = 'log2FoldChange', y = 'padj',
                          pCutoff = 0.05, FCcutoff = 0.585,
                          pointSize = 3,
                          titleLabSize = 20,
                          title = paste0("DEGs - ", cell, " (Infertile vs Fertile)"),
                          subtitle = NULL,
                          legendPosition = 'bottom',
                          caption = NULL)
  print(plot)
}

#After loop: Combine volcano plots
volcano_plots <- list()

#Generate plots from stored res_* objects
for (cell in cell_types) {
  res_var <- paste0("res_", gsub("-", "_", cell))
  
  if (exists(res_var)) {
    res_data <- get(res_var)
    
    plot <- EnhancedVolcano(res_data, lab = rownames(res_data),
                            x = 'log2FoldChange', y = 'padj',
                            pCutoff = 0.05, FCcutoff = 0.585,
                            pointSize = 2, axisLabSize = 10, labSize = 4,
                            title = paste0(cell, " DEGs - Infertile vs Fertile"),
                            titleLabSize = 14, subtitle = NULL,
                            legendPosition = 'none', caption = NULL)
    
    volcano_plots[[cell]] <- plot
  }
}

#Combine and display all volcano plots with one shared legend
wrap_plots(volcano_plots, ncol = 3) +
  plot_annotation(theme = theme(plot.title = element_text(size = 18, hjust = 0.5))) +
  plot_layout(guides = "collect")  # Collect and display one combined legend


###Combine DEGs into 1 dataframe
#Merge and add a column to denote the source
merged_DEGs <- bind_rows(sig_res_Pre_Unciliated      %>% mutate(source = "Pre-Unciliated"),
                         sig_res_Unciliated          %>% mutate(source = "Unciliated"),
                         sig_res_Ciliated            %>% mutate(source = "Ciliated"),
                         sig_res_Secretory           %>% mutate(source = "Secretory"),
                         sig_res_Pre_Ciliated        %>% mutate(source = "Pre-Ciliated"),
                         sig_res_Proliferative       %>% mutate(source = "Proliferative"))

#Total number of uniques genes across all cell subtypes
unique_genes_count <- merged_DEGs %>%
  summarise(unique_genes = n_distinct(gene))
print(unique_genes_count)

write.csv(merged_DEGs, "/data/gpfs/projects/punim1901/flames_v2/seurat_workspace/DEGs/merged_DEGs.csv")


###Volcano plots per cell subtype
pu <- EnhancedVolcano(res_Pre_Unciliated, lab = rownames(res_Pre_Unciliated),
                      x = 'log2FoldChange', y = 'padj',
                      pCutoff = 0.05, FCcutoff = 0.585,
                      pointSize = 2, labSize = 4, axisLabSize = 14,  
                      title = "Pre-Unciliated DEGs - Infertile vs Fertile", titleLabSize = 18, subtitle = NULL,
                      col = c("#717171","#717171","#717171","#F8766D"),
                      legendPosition = 'none', caption = NULL)
u <- EnhancedVolcano(res_Unciliated, lab = rownames(res_Unciliated),
                      x = 'log2FoldChange', y = 'padj',
                      pCutoff = 0.05, FCcutoff = 0.585,
                      pointSize = 2, labSize = 4, axisLabSize = 14,  
                      title = "Unciliated DEGs - Infertile vs Fertile", titleLabSize = 18, subtitle = NULL,
                      col = c("#717171","#717171","#717171","#ABA300"),
                      legendPosition = 'none', caption = NULL)
c <- EnhancedVolcano(res_Ciliated, lab = rownames(res_Ciliated),
                      x = 'log2FoldChange', y = 'padj',
                      pCutoff = 0.05, FCcutoff = 0.585,
                      pointSize = 2, labSize = 4, axisLabSize = 14,  
                      title = "Ciliated DEGs - Infertile vs Fertile", titleLabSize = 18, subtitle = NULL,
                      col = c("#717171","#717171","#717171","#0CB702"),
                      legendPosition = 'none', caption = NULL)
s <- EnhancedVolcano(res_Secretory, lab = rownames(res_Secretory),
                      x = 'log2FoldChange', y = 'padj',
                      pCutoff = 0.05, FCcutoff = 0.585,
                      pointSize = 2, labSize = 4, axisLabSize = 14,  
                      title = "Secretory DEGs - Infertile vs Fertile", titleLabSize = 18, subtitle = NULL,
                      col = c("#717171","#717171","#717171","#00BFC4"),
                      legendPosition = 'none', caption = NULL)
pc <- EnhancedVolcano(res_Pre_Ciliated, lab = rownames(res_Pre_Ciliated),
                      x = 'log2FoldChange', y = 'padj',
                      pCutoff = 0.05, FCcutoff = 0.585,
                      pointSize = 2, labSize = 4, axisLabSize = 14, 
                      title = "Pre-Ciliated DEGs - Infertile vs Fertile", titleLabSize = 18, subtitle = NULL,
                      col = c("#717171","#717171","#717171","#849AFF"),
                      legendPosition = 'none', caption = NULL)
p <- EnhancedVolcano(res_Proliferative, lab = rownames(res_Proliferative),
                      x = 'log2FoldChange', y = 'padj',
                      pCutoff = 0.05, FCcutoff = 0.585,
                      pointSize = 2, labSize = 4, axisLabSize = 14,  
                      title = "Proliferative DEGs - Infertile vs Fertile", titleLabSize = 18, subtitle = NULL,
                      col = c("#717171","#717171","#717171","#FF61CC"),
                      legendPosition = 'none', caption = NULL)  

(wrap_plots(pu, u, c, s, pc, p, ncol = 3)) +
  plot_annotation(theme = theme(plot.title = element_text(size = 18, hjust = 0.5)))
