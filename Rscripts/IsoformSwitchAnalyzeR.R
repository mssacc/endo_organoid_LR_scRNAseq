# =========================================================================
# Differential Transcript Usage (DTU) Analysis using IsoformSwitchAnalyzeR
# =========================================================================

#Load libraries
library(pfamAnalyzeR)
library(IsoformSwitchAnalyzeR)
library(rtracklayer)
library(Seurat)
library(tidyr)

#Load integrated Seurat object
all_samples_integrated <- readRDS("/data/gpfs/projects/punim1901/flames_v2/seurat_workspace/all_samples_integrated.rds")

#Define cell subtypes of interest
cell_subtypes <- c("Pre-Unciliated", "Unciliated", "Ciliated", "Secretory", "Pre-Ciliated", "Proliferative")

#Load GTF and gene name reference
tx2gene <- import("/data/gpfs/projects/punim1901/flames_v2/seurat_workspace/flames_output_copies/isoform_annotated.gtf")
resource_table <- read.csv("/data/gpfs/projects/punim1901/flames_v2/naming_reference.csv", header = TRUE)
#Output directory
output_dir <- "/data/gpfs/projects/punim1901/flames_v2/seurat_workspace/IsoformSwitchAnalyzeR"


# =========================================================================
#Loop through each cell subtype
for (subtype in cell_subtypes) {
  cat("\n\nProcessing cell subtype:", subtype, "\n")
  
  #Subset Seurat object
  subset_obj <- subset(all_samples_integrated, subset = cell_type == subtype)
  cat("Cells in subset:", ncol(subset_obj), "\n")
  if (ncol(subset_obj) < 25) {
    cat("Skipping", subtype, "- too few cells.\n")
    next}
  
  #Get pseudobulk counts for isoform assay
  pseudo.seurat.isoforms <- AggregateExpression(subset_obj, assays = "iso", return.seurat = FALSE,
                                                group.by = c("orig.ident", "fertility"))
  pseudo.seurat.isoforms.df <- as.data.frame(pseudo.seurat.isoforms)
  
  #Create sample metadata
  sample_data <- as.data.frame(colnames(pseudo.seurat.isoforms.df))
  colnames(sample_data) <- "colnames"
  split_sample_data <- strsplit(sample_data$colnames, "_")
  samps <- data.frame(sampleID = sample_data$colnames,
                      condition = as.factor(sapply(split_sample_data, `[`, 2)))
  
  #Clean transcript IDs in rownames
  rownames_split <- strsplit(rownames(pseudo.seurat.isoforms.df), "-")
  rownames(pseudo.seurat.isoforms.df) <- sapply(rownames_split, `[`, 1)
  
  #Match counts with GTF
  subset_txt2gene <- tx2gene[tx2gene$transcript_id %in% rownames(pseudo.seurat.isoforms.df), ]
  subset_cts <- pseudo.seurat.isoforms.df[rownames(pseudo.seurat.isoforms.df) %in% tx2gene$transcript_id, ]
  subset_txt2gene$isoform_id <- subset_txt2gene$transcript_id
  subset_txt2gene$gene_id <- subset_txt2gene$gene_id
  
  #Create switch list
  aSwitchList <- importRdata(isoformCountMatrix   = subset_cts,
                             isoformRepExpression = subset_cts,
                             designMatrix         = samps,
                             isoformExonAnnoation = "/data/gpfs/projects/punim1901/flames_v2/seurat_workspace/flames_output_copies/isoform_annotated.gtf",
                             isoformNtFasta       = "/data/gpfs/projects/punim1901/flames_v2/seurat_workspace/flames_output_copies/transcript_assembly.fa",
                             showProgress         = TRUE,
                             detectUnwantedEffects = FALSE) ###Added to get pre-ciliated to run (not for the rest)
                             

  #Analyze ORFs
  aSwitchList_ORF <- analyzeORF(switchAnalyzeRlist = aSwitchList, orfMethod = "longest")
  
  #Add gene names to isoformFeatures
  isoform_features <- aSwitchList_ORF[["isoformFeatures"]]
  isoform_features[["gene_name"]] <- resource_table$gene_symbol[match(isoform_features[["gene_id"]], resource_table$gene_id)]
  aSwitchList_ORF[["isoformFeatures"]] <- isoform_features
  
  #Filter low expression
  SwitchListFiltered <- preFilter(switchAnalyzeRlist = aSwitchList_ORF,
                                  geneExpressionCutoff = 30,
                                  isoformExpressionCutoff = 20,
                                  removeSingleIsoformGenes = TRUE)
  
  #Differential isoform usage test
  SwitchListAnalyzed <- isoformSwitchTestDEXSeq(switchAnalyzeRlist = SwitchListFiltered,
                                                reduceToSwitchingGenes = FALSE,
                                                reduceFurtherToGenesWithConsequencePotential = FALSE,
                                                alpha = 0.05,
                                                dIFcutoff = 0.25,
                                                onlySigIsoforms = FALSE)
  
  #Save results as CSV
  iso_feat <- as.data.frame(SwitchListAnalyzed$isoformFeatures)
  csv_file <- file.path(output_dir, paste0("IsoformSwitchResults_", gsub(" ", "_", subtype), ".csv"))
  write.csv(iso_feat, file = csv_file, row.names = FALSE)
  cat("Results saved to:", csv_file, "\n")
  
  #Save object to environment
  analyzed_name <- paste0("SwitchListAnalyzed_", gsub(" ", "_", subtype))
  assign(analyzed_name, SwitchListAnalyzed)
  
  #Save as RDS
  rds_file <- file.path(output_dir, paste0(analyzed_name, ".rds"))
  saveRDS(SwitchListAnalyzed, file = rds_file)
  cat("RDS file saved to:", rds_file, "\n")}

