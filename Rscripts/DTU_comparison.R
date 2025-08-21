###DTU comparison - Rscript

#Load required libraries
library(dplyr)
library(tidyr)
library(stringr)
library(pheatmap)
library(tibble)
library(readxl)

#Helper function to remove version numbers from ENST IDs
strip_version <- function(id) sub("\\.\\d+$", "", id)

#Step 1: Filter DTU transcripts found in ≥2 of 3 methods (ISA, DTUrtle, isopod)
#Load DETs
pre_unciliated_df <- read.csv("/data/gpfs/projects/punim1901/flames_v2/seurat_workspace/DETs/DE_analysis_Pre_Unciliated.csv")
unciliated_df     <- read.csv("/data/gpfs/projects/punim1901/flames_v2/seurat_workspace/DETs/DE_analysis_Unciliated.csv")
ciliated_df       <- read.csv("/data/gpfs/projects/punim1901/flames_v2/seurat_workspace/DETs/DE_analysis_Ciliated.csv")
secretory_df      <- read.csv("/data/gpfs/projects/punim1901/flames_v2/seurat_workspace/DETs/DE_analysis_Secretory.csv")
pre_ciliated_df   <- read.csv("/data/gpfs/projects/punim1901/flames_v2/seurat_workspace/DETs/DE_analysis_Pre_Ciliated.csv")
proliferative_df  <- read.csv("/data/gpfs/projects/punim1901/flames_v2/seurat_workspace/DETs/DE_analysis_Proliferative.csv")

#Merge DETs
DE_results <- bind_rows(pre_unciliated_df %>% mutate(source = "Pre-Unciliated"),
                        unciliated_df     %>% mutate(source = "Unciliated"),
                        ciliated_df       %>% mutate(source = "Ciliated"),
                        secretory_df      %>% mutate(source = "Secretory"),
                        pre_ciliated_df   %>% mutate(source = "Pre-Ciliated"),
                        proliferative_df  %>% mutate(source = "Proliferative"))

#Split into transcript_id and gene_id, remove version number
DE_results_iso <- DE_results %>%
  separate(transcript, into = c("transcript_id", "gene_id"), sep = "-", extra = "merge") %>%
  mutate(transcript_id = strip_version(transcript_id))

#Load DTU data
file_path <- "/data/gpfs/projects/punim1901/SC_paper/DTU_comparison_new.xlsx"

ISA <- read_excel(file_path, sheet = "ISA") %>%
  dplyr::select(cell_type, transcript_id, gene_name) %>%
  distinct() %>%
  mutate(dataset = "ISA",
         transcript_id = strip_version(transcript_id))

DTUrtle <- read_excel(file_path, sheet = "DTUrtle") %>%
  dplyr::select(cell_type, transcript_id, gene_name) %>%
  distinct() %>%
  mutate(dataset = "DTUrtle",
         transcript_id = strip_version(transcript_id))

isopod <- read_excel(file_path, sheet = "IsoPod") %>%
  dplyr::select(cell_type, transcript_id, gene_name) %>%
  distinct() %>%
  mutate(dataset = "isopod",
         transcript_id = strip_version(transcript_id))

#Create join keys
ISA <- ISA %>% mutate(join_key = paste(cell_type, transcript_id, sep = "_"))
DTUrtle <- DTUrtle %>% mutate(join_key = paste(cell_type, transcript_id, sep = "_"))
isopod <- isopod %>% mutate(join_key = paste(cell_type, transcript_id, sep = "_"))
DE_results_iso <- DE_results_iso %>%
  mutate(join_key = paste(source, transcript_id, sep = "_"))

#Mark DTUs that are also DEIs
ISA <- ISA %>%
  mutate(DET = if_else(join_key %in% DE_results_iso$join_key, "Yes", "No")) %>%
  dplyr::select(-join_key)
DTUrtle <- DTUrtle %>%
  mutate(DET = if_else(join_key %in% DE_results_iso$join_key, "Yes", "No")) %>%
  dplyr::select(-join_key)
isopod <- isopod %>%
  mutate(DET = if_else(join_key %in% DE_results_iso$join_key, "Yes", "No")) %>%
  dplyr::select(-join_key)

#Combine DTU results
dtu_combined <- bind_rows(ISA, DTUrtle, isopod) %>%
  filter(DET == "Yes") %>%
  dplyr::select(cell_type, transcript_id, dataset)

#Count number of datasets per (cell_type, transcript_id) pair
dtu_freq <- dtu_combined %>%
  distinct(cell_type, transcript_id, dataset) %>%
  group_by(cell_type, transcript_id) %>%
  summarise(n_methods = n(), .groups = "drop")

#Keep DTUs found in at least 2 datasets
DTU_DET <- dtu_freq %>%
  filter(n_methods >= 2)

#Step 2: Filter DE_results_iso and join transcript names
#Join DTU_DET to DE_results_iso using both transcript_id and source/cell_type
filtered_DE <- DE_results_iso %>%
  inner_join(DTU_DET, by = c("transcript_id", "source" = "cell_type"))

#Join to get full transcript names
filtered_DE <- filtered_DE %>%
  left_join(DE_results %>%
      mutate(transcript_id_match = strip_version(str_extract(transcript, "^[^-]+"))) %>%
      dplyr::select(transcript_id_match, source, transcript),
    by = c("transcript_id" = "transcript_id_match", "source"))


#Step 3: Prepare heatmap matrix
logFC_mat_df <- filtered_DE %>%
  group_by(transcript, source) %>%
  summarise(log2FoldChange = mean(log2FoldChange, na.rm = TRUE), .groups = "drop")

logFC_wide <- logFC_mat_df %>%
  pivot_wider(names_from = source, values_from = log2FoldChange, values_fill = 0) %>%
  column_to_rownames(var = "transcript")

#Reorder columns
desired_order <- c("Pre-Unciliated", "Unciliated", "Ciliated", "Secretory", "Pre-Ciliated", "Proliferative")
logFC_wide <- logFC_wide[, intersect(desired_order, colnames(logFC_wide))]

#Order rows by average log2FC
annotation_row <- data.frame(avg_log2FC = rowMeans(logFC_wide, na.rm = TRUE))
rownames(annotation_row) <- rownames(logFC_wide)
logFC_wide <- logFC_wide[order(-annotation_row$avg_log2FC), ]

#Remove one iso that keeps coming up even though it's only in 1 DTU method
logFC_wide <- logFC_wide[!rownames(logFC_wide) %in% "ENST00000567998.5-SULT1A1", ]


#Step 4: Plot heatmap
my_breaks <- seq(-8, 8, length.out = 100)
legend_ticks <- seq(-8, 8, by = 2)

pheatmap(logFC_wide,
         show_rownames = TRUE,
         cluster_rows = FALSE, cluster_cols = FALSE,
         treeheight_col = 0, treeheight_row = 0,
         angle_col = 0,
         breaks = my_breaks, legend_breaks = legend_ticks,
         main = "Isoform Expression of DTUs Found in ≥2 Methods")
