#Generating reference files for SC analysis

#Load required libraries
library(rtracklayer)
library(gridExtra)
library(data.table)
library(BiocParallel)
library(celda)
library(stringr)
library(cowplot)
library(grid)
library(patchwork)
library(tidyverse)

#Set working directory
setwd("/data/gpfs/projects/punim1901/flames_v2/flames_output/")


#Import the first GTF file (transcripts GTF)
gtf1 <- import("/data/gpfs/projects/punim1901/flames_v2/flames_output/isoform_annotated.gtf")
gtf1_df <- as.data.frame(gtf1)
  
#Select relevant columns from the first GTF
selected_columns1 <- gtf1_df[, c("transcript_id", "gene_id")]
unique_selected_cols <- unique(selected_columns1)
  
#Import the second GTF file (reference GTF with gene symbols)
gtf2 <- read.csv("/home/mssacc/genomes/gene_info.csv")
gtf2_df <- as.data.frame(gtf2)

#Select relevant columns from the second GTF
selected_columns2 <- gtf2_df[, c("gene_symbol", "gene_id")]
unique_gene_symbol <- unique(selected_columns2)
  
#Merge the two data frames on 'gene_id'
combined_data <- merge(unique_selected_cols, 
                         unique_gene_symbol, 
                         by = "gene_id", 
                         all.x = TRUE)
  
#If 'gene_symbol' is missing, replace it with 'gene_id'
combined_data$gene_symbol <- ifelse(is.na(combined_data$gene_symbol), 
                                      combined_data$gene_id, 
                                      combined_data$gene_symbol)
  
#Select relevant columns
final_combined_data <- combined_data[, c("transcript_id", "gene_id", "gene_symbol")]
  
#Write to a CSV file
write.csv(final_combined_data, file = file.path("/data/gpfs/projects/punim1901/flames_v2/naming_reference.csv"), row.names = FALSE)
  
