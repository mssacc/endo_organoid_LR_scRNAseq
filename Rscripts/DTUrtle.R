###==========================================================================================
###DTUrtle Rscript — Loop through 6 cell subtypes
###==========================================================================================

#Load required libraries
library(DTUrtle)
library(Seurat)
library(stringr)
library(Matrix.utils)
library(GenomicRanges)
library(pheatmap)
library(dplyr)

#Using 4 cores
biocpar <- BiocParallel::MulticoreParam(4)

#Import GTF Annotation to get transcript to gene mapping
tx2gene <- import_gtf(gtf_file = "/data/gpfs/projects/punim1901/flames_v2/seurat_workspace/flames_output_copies/isoform_annotated.gtf")
tx2gene$transcript_id <- sub("\\.\\d+$", "", tx2gene$transcript_id)

#Read gene_symbol mapping from resource table
resource_table <- read.csv("/data/gpfs/projects/punim1901/flames_v2/naming_reference.csv", header = TRUE)

#Move transcript and gene identifier columns to front
tx2gene <- move_columns_to_front(df = tx2gene, columns = c("transcript_id", "gene_id"))

#Add gene_symbol to matrix
tx2gene <- merge(tx2gene, resource_table[, c("transcript_id", "gene_id", "gene_symbol")],
                 by = c("transcript_id", "gene_id"), all.x = TRUE)

#Load integrated Seurat object
all_samples_integrated <- readRDS("/data/gpfs/projects/punim1901/flames_v2/seurat_workspace/all_samples_integrated.rds")

#Define cell subtypes
cell_subtypes <- c("Pre-Unciliated", "Unciliated", "Ciliated",  "Secretory", "Pre-Ciliated", "Proliferative")

#Loop through each cell subtype
for (ctype in cell_subtypes) {
  
  message("===== Processing cell type: ", ctype, " =====")
  
  #Subset object
  obj <- subset(all_samples_integrated, subset = cell_type == ctype)
  
  #Extract isoform counts and clean rownames
  iso_mat <- GetAssayData(obj, assay = "iso", layer = "counts")
  split_names <- str_split_fixed(rownames(iso_mat), "-", n = 2)
  transcript_id <- sub("\\.\\d+$", "", split_names[, 1])
  gene_id <- sub("\\.\\d+$", "", split_names[, 2])
  rownames(iso_mat) <- transcript_id
  obj[["iso"]] <- CreateAssayObject(counts = iso_mat)
  
  #Prepare counts and sample info
  Idents(obj) <- "fertility"
  obj$group <- paste0(obj$orig.ident, "_", obj$fertility)
  
  cts <- as.matrix(iso_mat)
  cts_summed <- aggregate.Matrix(t(cts), groupings = obj$group, fun = "sum")
  cts_summed <- t(cts_summed)
  
  pd <- data.frame(id = colnames(cts_summed),
                   group = gsub(".*_", "", colnames(cts_summed)),
                   stringsAsFactors = FALSE)
  
  #Run DTUrtle
  dturtle <- run_drimseq(counts = cts_summed, tx2gene = tx2gene, pd = pd,
                         id_col = "id", cond_col = "group", filtering_strategy = "bulk",
                         BPPARAM = biocpar)
  
  #Posthoc filtering & StageR
  dturtle <- posthoc_and_stager(dturtle = dturtle, ofdr = 0.05)
  dturtle <- create_dtu_table(dturtle = dturtle)
  
  sig_tx_qvals <- dturtle$FDR_table %>%
    filter(txID %in% dturtle$sig_tx) %>%
    select(txID, geneID, transcript, gene)
  
  #Create output directory
  outdir <- file.path("/data/gpfs/projects/punim1901/flames_v2/seurat_workspace/DTUrtle", gsub(" ", "_", tolower(ctype)))
  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
  
  #Save outputs
  saveRDS(dturtle, file = file.path(outdir, "dturtle_object.rds"))
  saveRDS(cts_summed, file = file.path(outdir, "cts_summed.rds"))
  saveRDS(pd, file = file.path(outdir, "phenotype_data.rds"))
  saveRDS(tx2gene, file = file.path(outdir, "tx2gene_mapping.rds"))
  write.csv(dturtle$dtu_table, file.path(outdir, "dtu_table.csv"), row.names = FALSE)
  write.csv(sig_tx_qvals, file.path(outdir, "sig_transcript_qvals.csv"), row.names = FALSE)
  
  #UMAP coordinate extraction & plotting
  umap_df <- as.data.frame(Embeddings(obj, reduction = "umap"))
  umap_df$group <- obj$group
  reduction_df <- aggregate(. ~ group, data = umap_df, FUN = mean)
  rownames(reduction_df) <- reduction_df$group
  reduction_df <- reduction_df[, c("umap_1", "umap_2")]
  
  plot4 <- plot_dimensional_reduction(dturtle = dturtle,
                                      reduction_df = reduction_df,
                                      savepath = outdir,
                                      add_to_table = "dimensional_reduction",
                                      BPPARAM = biocpar)
                                    
  saveRDS(reduction_df, file = file.path(outdir, "umap_group_means.rds"))
  
  message("===== Finished: ", ctype, " =====\n")}

message("All cell types processed successfully.")
