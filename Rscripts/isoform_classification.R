#Isoform classification

#Load required libraries
library(Seurat)
library(dplyr)
library(tidyr)
library(ggplot2)
library(biomaRt)
library(ggrepel)
library(scales)

#Import datasets
all_samples_integrated <- readRDS("/data/gpfs/projects/punim1901/flames_v2/seurat_workspace/all_samples_integrated.rds")
SQANTI <- read.csv("/data/gpfs/projects/punim1901/flames_v2/flames_output/remove_unknownstrand_classification.txt", sep ='\t')

###Analyze isoforms by all samples
#Aggregate the expression data by cell type 
counts <- AggregateExpression(all_samples_integrated, 
                              assays = "iso", 
                              return.seurat = FALSE,
                              group.by="group")
as.data.frame(counts) -> df
row.names(df) -> df$gene
  
#Split transcript ids into gene and transcript id
pseudobulk_data <- df %>% separate(gene, into = c("transcript_id", "gene_id"), sep = "-",  extra = "merge") 
            
#We can use the pseudobulk_data counts calculated previously 
#Select cell type columns and gather into long format
long_data <- pseudobulk_data %>%
    pivot_longer(cols = starts_with("iso"),
                 names_to = "all_cells",
                 values_to = "expression")
            
#Filter for non-zero expression to identify genes and isoforms expressed in each cell type
filtered_data <- long_data %>%
  filter(expression > 0) %>%
    dplyr::select(all_cells, transcript_id, gene_id) %>%
    distinct()
            
#Calculate unique gene and isoform counts for each cell type
cell_type_summary <- filtered_data %>%
    group_by(all_cells) %>%
    summarise(num_genes = n_distinct(gene_id),
              num_isoforms = n_distinct(transcript_id)) %>%
    pivot_longer(cols = c(num_genes, num_isoforms), 
                 names_to = "Category", 
                 values_to = "Count")
                        
#Plot the data
p1 <- ggplot(cell_type_summary, aes(x = all_cells, y = Count, fill = Category)) +
      geom_bar(stat = "identity", position = "dodge", color = "black", width = 0.8) +
      theme_minimal() +
      labs(title = "Number of genes and isoforms", y = "Count") +
      theme(plot.title = element_text(size = 20, face = "bold", hjust = 0.5),
            axis.title.x = element_text(size = 14),
            axis.title.y = element_text(size = 14),
            axis.text.x = element_text(size = 12, angle = 45, hjust = 1),
            axis.text.y = element_text(size = 12),
            legend.title = element_blank())
p1            

#Plot the structural category per cell type
merged_data <- merge(pseudobulk_data, SQANTI, 
                     by.x = "transcript_id", 
                     by.y = "isoform", 
                     all.x = TRUE)

#Pivot pseudobulk data to long format for cell types
long_data <- merged_data %>%
              pivot_longer(cols = starts_with("iso"),
                           names_to = "all_cells",
                           values_to = "expression")

#Generate our new df that we can use for plotting attributes defined by cell type
filtered_data <- long_data %>%
   filter(expression > 0) %>%
  distinct()

#Calculate unique counts of genes and isoforms per structural category and cell type
category_summary <- filtered_data %>%
    group_by(all_cells, structural_category, subcategory, coding) %>%
    summarise(num_genes = n_distinct(gene_id),
              num_isoforms = n_distinct(transcript_id),
              .groups = "drop")
            
#Plot number of isoforms per structural category
p2 <- ggplot(category_summary, aes(x = all_cells, y = num_isoforms, fill = structural_category)) +
        geom_bar(stat = "identity", position = "stack", width = 0.8) +
        theme_minimal() +
        labs(title = "Isoforms by structural category",
             y = "Number of Isoforms", fill = "Structural Category") +
        theme(plot.title = element_text(size = 20, face = "bold", hjust = 0.5),
              axis.title.x = element_text(size = 14),
              axis.title.y = element_text(size = 14),
              axis.text.x = element_text(size = 12, angle = 45, hjust = 1),
              axis.text.y = element_text(size = 12),
              legend.title = element_text(size = 14),
              legend.text = element_text(size = 12))

#Plot isoforms by structural subcategory            
p3 <- ggplot(category_summary, aes(x = all_cells, y = num_isoforms, fill = subcategory)) +
       geom_bar(stat = "identity", position = "stack", width = 0.8) +
       theme_minimal() +
       labs(title = "Isoforms by structural subcategory",
            y = "Number of Isoforms", fill = "subcategory") +
      theme(plot.title = element_text(size = 20, face = "bold", hjust = 0.5),
            axis.title.x = element_text(size = 14),
            axis.title.y = element_text(size = 14),
            axis.text.x = element_text(size = 12, angle = 45, hjust = 1),
            axis.text.y = element_text(size = 12),
            legend.title = element_text(size = 14),
            legend.text = element_text(size = 12))

#Plot coding vs non-coding isoforms  
p4 <- ggplot(category_summary, aes(x = all_cells, y = num_isoforms, fill = coding)) +
        geom_bar(stat = "identity", position = "fill", width = 0.8) +
        theme_minimal() +
        labs(title = "Proportion of coding isoforms",
             y = "Proportion of Isoforms", fill = "Structural Category") +
        scale_y_continuous(labels = scales::percent) +  # This will show y-axis as percentages
        theme(plot.title = element_text(size = 20, face = "bold", hjust = 0.5),
              axis.title.x = element_text(size = 14),
              axis.title.y = element_text(size = 14),
              axis.text.x = element_text(size = 12, angle = 45, hjust = 1),
              axis.text.y = element_text(size = 12),
              legend.title = element_text(size = 14),
              legend.text = element_text(size = 12))
            
cowplot::plot_grid(p1, p2, p3, p4, ncol = 2)


#Remove version numbers from transcript IDs
filtered_data$transcript_id <- sub("\\..*$", "", filtered_data$transcript_id)

#Retrieve gene biotype from Ensembl
mart <- useEnsembl(biomart = "ENSEMBL_MART_ENSEMBL",
                   dataset = "hsapiens_gene_ensembl",
                   host = "https://jan2020.archive.ensembl.org")
biotype_info_transcript <- getBM(attributes = c("ensembl_transcript_id","transcript_biotype","transcript_length"),
                                 filters = 'ensembl_transcript_id',
                                 values = filtered_data$transcript_id, mart = mart)

#Merge BioMart data into filtered expression data
merged_biomart_data <- merge(filtered_data, biotype_info_transcript,
                             by.x = "transcript_id",
                             by.y = "ensembl_transcript_id",
                             all.x = TRUE)

#Reassign transcript_biotype to "BambuTx" if transcript_id starts with "BambuTx"
merged_biomart_data$transcript_id <- trimws(as.character(merged_biomart_data$transcript_id))
merged_biomart_data <- merged_biomart_data %>%
  mutate(transcript_biotype = ifelse(grepl("^BambuTx", transcript_id), "BambuTx", transcript_biotype))

#Prepare the summary data
category_summary <- merged_biomart_data %>%
  group_by(transcript_biotype) %>%
  summarise(num_isoforms = n(), .groups = "drop") %>%
  mutate(proportion = num_isoforms / sum(num_isoforms)) %>%
  ungroup()

#Group transcript_biotypes with proportion < 0.5% into "Other"
category_summary <- category_summary %>%
  mutate(transcript_biotype = ifelse(proportion < 0.005, "Other", transcript_biotype)) %>%
  group_by(transcript_biotype) %>%
  summarise(num_isoforms = sum(num_isoforms), .groups = "drop") %>%
  mutate(proportion = num_isoforms / sum(num_isoforms)) %>%
  ungroup()

#Ensure consistent factor order
category_summary <- category_summary %>%
  mutate(transcript_biotype = factor(transcript_biotype, levels = c(
    "protein_coding", "lncRNA", "retained_intron", "processed_transcript",
    "nonsense_mediated_decay", "BambuTx", "processed_pseudogene", "TEC", "Other")))

#Custom labels
biotype_labels <- c(protein_coding = "Protein Coding",
                    lncRNA = "lncRNA",
                    retained_intron = "Retained Intron",
                    processed_transcript = "Processed Transcript",
                    nonsense_mediated_decay = "NMD",
                    BambuTx = "BambuTx",
                    processed_pseudogene = "Processed Pseudogene",
                    TEC = "TEC",
                    Other = "Other")

#Create pie chart of transcript subtypes
ggplot(category_summary, aes(x = "", y = proportion, fill = transcript_biotype)) +
  geom_bar(stat = "identity", width = 1, colour = "black", size = 0.2) +  
  coord_polar(theta = "y", start = pi / 6) +
  scale_fill_discrete(labels = biotype_labels) +
  labs(title = "BioMart Transcript Biotypes", fill = "Transcript Biotypes") +
  theme_void() +
  theme(plot.title = element_text(size = 20, face = "bold", hjust = 0.5, margin = margin(b = 10)),
        plot.margin = margin(10, 100, 10, 10),
        legend.title = element_text(size = 13, face = "bold"),
        legend.text = element_text(size = 12),
        strip.text = element_text(size = 14, face = "bold", margin = margin(b = 0)),
        panel.grid = element_blank())

#Count the number of unique BambuTx isoforms
num_bambu_isoforms <- merged_biomart_data %>%
  filter(transcript_biotype == "BambuTx") %>%
  distinct(transcript_id) %>%
  nrow()

cat("Number of unique BambuTx isoforms:", num_bambu_isoforms, "\n")
