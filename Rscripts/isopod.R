### Isopod DTU analysis - loop through cell types

#Load required libraries
library(isopod)
library(Seurat)
library(dplyr)
library(stringr)
library(Matrix)

#Load integrated Seurat object
all_samples_integrated <- readRDS("/data/gpfs/projects/punim1901/flames_v2/seurat_workspace/all_samples_integrated.rds")

#Set isoform assay
DefaultAssay(all_samples_integrated) <- "iso"

#Extract isoform count matrix
iso_counts <- GetAssayData(all_samples_integrated, layer = "counts")

#Convert to data.frame and keep transcript-gene IDs
counts_df <- as.data.frame(Matrix::as.matrix(iso_counts))
counts_df$transcript_gene_id <- rownames(counts_df)

#Split transcript-gene into transcript_id and gene_id
counts_df <- counts_df %>%
  mutate(transcript_id = sub("-[^-]+$", "", transcript_gene_id),
         gene_id = sub("^.*?-", "", transcript_gene_id)) %>%
  select(transcript_id, gene_id, everything(), -transcript_gene_id)

#Reorder columns
counts_df <- counts_df %>% relocate(transcript_id, gene_id)

#Combine fertility and cell_type metadata
all_samples_integrated$fertile_subtype <- paste(all_samples_integrated$fertility, 
                                                all_samples_integrated$cell_type, sep = "_")

#Define cell types to loop through
cell_types_to_test <- c("Pre-Unciliated", "Unciliated", "Ciliated", "Secretory", "Proliferative", "Pre-Ciliated")

#Loop through each cell type
for (cell_type in cell_types_to_test) {
  message("Running DTU analysis for cell type: ", cell_type)
  
  #Define group labels
  group1 <- paste0("Infertile_", cell_type)
  group2 <- paste0("Fertile_", cell_type)
  
  #Subset cells belonging to this cell type
  cells_to_use <- rownames(all_samples_integrated@meta.data)[
    all_samples_integrated$cell_type == cell_type]
  
  #Subset expression matrix to relevant cells
  expr_mat <- counts_df[, c("transcript_id", "gene_id", cells_to_use)]
  expr_cols <- colnames(expr_mat)[!(colnames(expr_mat) %in% c("transcript_id", "gene_id"))]

  ###Filtering
    #Step 1 – Keep isoforms with expression >= 20 in at least one sample
      isoform_pass <- apply(expr_mat[, expr_cols], 1, function(x) any(x >= 20))
      expr_mat <- expr_mat[isoform_pass, ]
  
    #Step 2 – Remove genes with only one isoform
      iso_per_gene <- table(expr_mat$gene_id)
      valid_genes <- names(iso_per_gene[iso_per_gene > 1])
      expr_mat <- expr_mat %>% filter(gene_id %in% valid_genes)
  
    #Step 3 – Remove genes with total expression < 30 across all cells
      gene_expr <- expr_mat %>%
        group_by(gene_id) %>%
        summarise(gene_expr = sum(across(all_of(expr_cols)))) %>%
        filter(gene_expr >= 30)
      expr_mat <- expr_mat %>% filter(gene_id %in% gene_expr$gene_id)


  #Create cell_groups data frame
  cell_groups <- all_samples_integrated@meta.data[cells_to_use, , drop = FALSE] %>%
    mutate(cell_id = rownames(.), cell_group = fertile_subtype) %>%
    select(cell_id, cell_group)

  #Run permutation test
  permutation_results <- get_permutation_pvals(expr_mat, cell_groups, 
                                               transcript_id_colname = 'transcript_id', 
                                               gene_id_colname = 'gene_id',
                                               cell_labels_colname = 'cell_group', 
                                               analysis_group_1 = group1,
                                               analysis_group_2 = group2,
                                               run_on_all_groups = FALSE, cores = 4,
                                               cutoff = 0.05,
                                               do_gene_level_comparisons = TRUE,
                                               report_adjusted_pvalues = TRUE)
  #Save results
  save_path <- paste0("/data/gpfs/projects/punim1901/flames_v2/seurat_workspace/isopod/filtered_results_new/permutation_results_", cell_type, ".rds")
  saveRDS(permutation_results, file = save_path)
}
