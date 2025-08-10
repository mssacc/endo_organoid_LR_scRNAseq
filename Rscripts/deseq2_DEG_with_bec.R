#Load required libraries
library(Seurat)
library(DESeq2)
library(tidyverse)
library(ggplot2)
library(limma)
library(EnhancedVolcano)

#Load RDS object
all_samples_integrated <- readRDS("/data/gpfs/projects/punim1901/flames_v2/seurat_workspace/all_samples_integrated.rds")


#Step 1 - Aggregate pseudobulk counts per group (ensure these columns exist in metadata)
pseudobulk_counts <- AggregateExpression(all_samples_integrated,
                                         group.by = c("orig.ident", "batch", "fertility", "cell_type"),
                                         assays = "RNA", slot = "counts")$RNA


#Step 2 - Create metadata for each pseudobulk group
group_metadata <- data.frame(
  group_id = colnames(pseudobulk_counts)) %>%
  separate(group_id, into = c("orig.ident", "batch", "fertility", "cell_type"), sep = "_") %>%
  left_join(all_samples_integrated@meta.data %>%
      dplyr::select(orig.ident, batch, fertility) %>%
      distinct(),
      by = "orig.ident")

#Fix group_metadata
group_metadata$group_id <- colnames(pseudobulk_counts)
group_metadata <- group_metadata %>%
  rename(batch = batch.x, fertility = fertility.x) %>%
  dplyr::select(-batch.y, -fertility.y) %>%
  drop_na(batch, fertility)


#Step 3 - apply limma batch effect correction
apply_limma_batch_correction <- function(counts_matrix, batch_info) {
  # Convert counts to log2 scale
  log_counts <- log2(counts_matrix + 1)
  # Create design matrix
  design <- model.matrix(~fertility, data = batch_info)
  # Apply limma batch correction
  corrected_counts <- removeBatchEffect(log_counts,
                                        batch = batch_info$batch,
                                        design = design)
  # Convert back to counts scale
  corrected_counts <- 2^corrected_counts - 1
  corrected_counts[corrected_counts < 0] <- 0

  return(corrected_counts)
}
counts_corrected <- apply_limma_batch_correction(pseudobulk_counts, group_metadata)


#Step 4 - Run DESeq2 with fertility as condition and batch as covariate
dds <- DESeqDataSetFromMatrix(round(counts_corrected),
                              colData = group_metadata,
                              design = ~ fertility)
#Filter data
#Remove genes with < 10 reads
keep <- rowSums(counts(dds)) >=10
dds <- dds[keep,]

#Run DESeq2
dds <- DESeq(dds)

#Step 5 - Extract DEGs between Infertile and Fertile
res <- results(dds, contrast = c("fertility", "Infertile", "Fertile"))

#Make gene name
res <- res %>%
  as.data.frame() %>%
  rownames_to_column("gene")

#Save DE analysis with all genes - not just DEGs
write.csv(res, "/data/gpfs/projects/punim1901/flames_v2/seurat_workspace/DEGs/bulk_DE_analysis.csv", row.names = FALSE)


#Step 6 - Filter and plot significant DEGs
sig_res <- res %>%
  filter(padj < 0.05, abs(log2FoldChange) > 0.585)

###Volcano plot
#Define colours for volcano plot
keyvals <- ifelse(
  res$log2FoldChange < -0.585 & res$padj < 0.05, '#F8766D',  
  ifelse(res$log2FoldChange > 0.585 & res$padj < 0.05, 'palegreen3',  
         '#717171'))  
keyvals[is.na(keyvals)] <- '#717171'
keyvals[is.na(keyvals)] <- '#717171'
names(keyvals)[keyvals == 'palegreen3'] <- 'Upregulated'
names(keyvals)[keyvals == '#F8766D'] <- 'Downregulated'
names(keyvals)[keyvals == '#717171'] <- 'Not Significant'
#Generate plot
EnhancedVolcano(res, lab = rownames(res),
                x = 'log2FoldChange', y = 'padj',
                pCutoff = 0.05, FCcutoff = 0.585,
                pointSize = 2, labSize = 4, axisLabSize = 14,  
                title = "Bulk DEGs - Infertile vs Fertile", titleLabSize = 25,
                subtitle = NULL,
                colCustom = keyvals,
                legendPosition = 'right')

#Save DEGs
write.csv(sig_res, "/data/gpfs/projects/punim1901/flames_v2/seurat_workspace/DEGs/bulk_DEGs.csv", row.names = FALSE)


#Step 7 - Run plot
#VST for visualization
vsd <- vst(dds, blind = TRUE)
# Extract the transformed expression matrix
vsd_matrix <- assay(vsd)
# Transpose for PCA (samples as rows, genes as columns -> must flip for prcomp)
pca <- prcomp(t(vsd_matrix))
# Run PCA on corrected VST
pca_df <- as.data.frame(pca$x)
pca_df$group_id <- rownames(pca_df)
#Add metadata
pca_df <- pca_df %>%
  left_join(group_metadata, by = "group_id")


#Step 8 - Plot PCA
ggplot(pca_df, aes(x = PC1, y = PC2, color = fertility, shape = cell_type)) +
  geom_point(size = 3, alpha = 0.8) +
  theme_minimal(base_size = 14) +
  labs(title = "Pseudobulk PCA (fertility)",
       x = paste0("PC1 (", round(summary(pca)$importance[2, 1] * 100, 1), "%)"),
       y = paste0("PC2 (", round(summary(pca)$importance[2, 2] * 100, 1), "%)"))

ggplot(pca_df, aes(x = PC1, y = PC2, color = batch, shape = cell_type)) +
  geom_point(size = 3, alpha = 0.8) +
  theme_minimal(base_size = 14) +
  labs(title = "Pseudobulk PCA (batch)",
       x = paste0("PC1 (", round(summary(pca)$importance[2, 1] * 100, 1), "%)"),
       y = paste0("PC2 (", round(summary(pca)$importance[2, 2] * 100, 1), "%)"))

