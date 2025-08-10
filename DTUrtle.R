###==========================================================================================
###DTUrtle Rscript
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

#import gtf Annotation to get transcript to gene mapping
tx2gene <- import_gtf(gtf_file = "/data/gpfs/projects/punim1901/flames_v2/seurat_workspace/flames_output_copies/isoform_annotated.gtf")
  #Remove version number from transcript_id
  tx2gene$transcript_id <- sub("\\.\\d+$", "", tx2gene$transcript_id)

#Read gene_symbol mapping from resource table
resource_table <- read.csv("/data/gpfs/projects/punim1901/flames_v2/naming_reference.csv", header = TRUE)

#Move transcript and gene identifier columns to front
tx2gene <- move_columns_to_front(df = tx2gene, 
                                 columns = c("transcript_id", "gene_id"))
#Add gene_symbol to matrix
tx2gene <- merge(tx2gene, resource_table[, c("transcript_id", "gene_id", "gene_symbol")],
                 by = c("transcript_id", "gene_id"),
                 all.x = TRUE)


#Load integrated Seurat object
all_samples_integrated <- readRDS("/data/gpfs/projects/punim1901/flames_v2/seurat_workspace/all_samples_integrated.rds")
all_samples_integrated <- subset(all_samples_integrated, subset = cell_type == "Proliferative")


#Extract isoform counts from Seurat object and clean rownames
iso_mat <- GetAssayData(all_samples_integrated, assay = "iso", layer = "counts")

#Split transcript_id-gene_id (only on first hyphen to preserve gene names with hyphens)
split_names <- str_split_fixed(rownames(iso_mat), "-", n = 2)
transcript_id <- sub("\\.\\d+$", "", split_names[, 1])
gene_id <- sub("\\.\\d+$", "", split_names[, 2])

#Update rownames in iso matrix and reassign to Seurat object
rownames(iso_mat) <- transcript_id
all_samples_integrated[["iso"]] <- CreateAssayObject(counts = iso_mat)

#Create transcript-to-gene map from Seurat counts (optional: for filtering or checking)
iso_gene_map <- data.frame(transcript_id = transcript_id,
                           gene_id = gene_id,
                           row.names = transcript_id,
                           stringsAsFactors = FALSE)

#Prepare counts and sample sheet for DTUrtle
cts <- as.matrix(iso_mat)


######Figure out which way to do this...maybe subset object by cell subtype and run this analysis separately for each cell???
#Prepare pseudobulk sample info (group by sample + fertility)
Idents(all_samples_integrated) <- "fertility"  # or whatever metadata you want
all_samples_integrated$group <- paste0(all_samples_integrated$orig.ident, "_", all_samples_integrated$fertility)

# Pseudobulk: sum counts across cells by group
cts_summed <- aggregate.Matrix(t(cts), groupings = all_samples_integrated$group, fun = "sum")
cts_summed <- t(cts_summed)

#Create phenotype data
pd <- data.frame(id = colnames(cts_summed),
                 group = gsub(".*_", "", colnames(cts_summed)),  # extract fertility from "sample_fertility"
                 stringsAsFactors = FALSE)

#Run DRIMSeq DTU analysis
dturtle <- run_drimseq(counts = cts_summed, tx2gene = tx2gene, pd = pd,
                       id_col = "id", cond_col = "group", filtering_strategy = "bulk",
                       BPPARAM = biocpar)

#Posthoc filtering and stageR
dturtle <- posthoc_and_stager(dturtle = dturtle, ofdr = 0.05)


#highly flexible function to create a results data frame
dturtle <- create_dtu_table(dturtle = dturtle)

## View results data frame
View(dturtle$dtu_table)

#Extract q-values for significant transcripts
sig_tx_qvals <- dturtle$FDR_table %>%
  filter(txID %in% dturtle$sig_tx) %>%
  select(txID, geneID, transcript, gene)
View(sig_tx_qvals)



# Save DTUrtle object (can be reloaded later for plotting without rerunning everything)
saveRDS(dturtle, file = "/data/gpfs/projects/punim1901/flames_v2/seurat_workspace/DTUrtle/proliferative/dturtle_object.rds")  # === Save object ===

# Save intermediate objects (optional, for reproducibility/debugging)
saveRDS(cts_summed, file = "/data/gpfs/projects/punim1901/flames_v2/seurat_workspace/DTUrtle/proliferative/cts_summed.rds")
saveRDS(pd, file = "/data/gpfs/projects/punim1901/flames_v2/seurat_workspace/DTUrtle/proliferative/phenotype_data.rds")
saveRDS(tx2gene, file = "/data/gpfs/projects/punim1901/flames_v2/seurat_workspace/DTUrtle/proliferative/tx2gene_mapping.rds")

# Create DTU table
dturtle <- create_dtu_table(dturtle = dturtle)

# Save the DTU table as CSV
write.csv(dturtle$dtu_table,
          "/data/gpfs/projects/punim1901/flames_v2/seurat_workspace/DTUrtle/proliferative/dtu_table.csv",
          row.names = FALSE)  # === Save results table ===

# Save sig tx q-values
write.csv(sig_tx_qvals,
          "/data/gpfs/projects/punim1901/flames_v2/seurat_workspace/DTUrtle/proliferative/sig_transcript_qvals.csv",
          row.names = FALSE)  # === Save filtered q-values ===


#Extract UMAP coords for Pre-Unciliated cells
umap_df <- as.data.frame(Embeddings(all_samples_integrated, reduction = "umap"))

#Add group info (used in dturtle pseudobulk step)
umap_df$group <- all_samples_integrated$group  # e.g., E170_Fertile

#Average UMAP coordinates per group (sample)
reduction_df <- aggregate(. ~ group, data = umap_df, FUN = mean)

#Rename columns to match DTUrtle expectations
rownames(reduction_df) <- reduction_df$group
reduction_df <- reduction_df[, c("umap_1", "umap_2")]  # keep only 2 columns

#Now this should work
plot4 <- plot_dimensional_reduction(
  dturtle = dturtle,
  reduction_df = reduction_df,
  savepath = "/data/gpfs/projects/punim1901/flames_v2/seurat_workspace/DTUrtle/proliferative", 
  add_to_table = "dimensional_reduction",
  BPPARAM = biocpar)

#Save reduction dataframe
saveRDS(reduction_df, file = "/data/gpfs/projects/punim1901/flames_v2/seurat_workspace/DTUrtle/proliferative/umap_group_means.rds")











# 
# 
# 
# 
# #change to results folder
# setwd("/data/gpfs/projects/punim1901/flames_v2/seurat_workspace/DTUrtle/proliferative")    
# 
# #create plots, save them to disk and link them in the `dtu_table`.
# plot1 <- plot_proportion_barplot(dturtle = dturtle, 
#                                  savepath = "/data/gpfs/projects/punim1901/flames_v2/seurat_workspace/DTUrtle/proliferative", 
#                                  add_to_table = "barplot",
#                                  BPPARAM = biocpar)
# 
# plot2 <- plot_proportion_pheatmap(dturtle = dturtle, 
#                                     savepath = "/data/gpfs/projects/punim1901/flames_v2/seurat_workspace/DTUrtle/proliferative", 
#                                     include_expression = TRUE,
#                                     add_to_table = "pheatmap",
#                                     BPPARAM = biocpar)


###Again can't do this function because of no exon info in gtf
        #Now create the GRanges object
        # gr <- GRanges(
        #   seqnames = tx2gene$seqnames,
        #   ranges = IRanges(start = as.integer(tx2gene$start), end = as.integer(tx2gene$end)),
        #   strand = tx2gene$strand,
        #   transcript_id = tx2gene$transcript_id,
        #   gene_id = tx2gene$gene_id,
        #   gene_name = tx2gene$gene_symbol,
        #   source = tx2gene$source,
        #   type = tx2gene$type,
        #   score = tx2gene$score,
        #   phase = tx2gene$phase)
        # 
        # mcols(gr) <- mcols(gr)[, c("transcript_id", "gene_id", "gene_name", "type")]
        # 
        # plot3 <- plot_transcripts_view(dturtle = dturtle, 
        #                                  gtf = gr, 
        #                                  genome = 'hg38', 
        #                                  one_to_one = TRUE,
        #                                  savepath = "/data/gpfs/projects/punim1901/flames_v2/seurat_workspace/DTUrtle/proliferative", 
        #                                  add_to_table = "transcript_view",
        #                                  BPPARAM = biocpar)



