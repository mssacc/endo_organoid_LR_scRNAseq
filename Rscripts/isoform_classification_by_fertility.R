###Isoform classification divided by fertility

#Load required libraries
library(Seurat)
library(dplyr)
library(tidyr)
library(biomaRt)
library(ggplot2)
library(ggrepel)
library(scales)

#Import datasets
all_samples_integrated <- readRDS("/data/gpfs/projects/punim1901/flames_v2/seurat_workspace/all_samples_integrated.rds")
SQANTI <- read.csv("/data/gpfs/projects/punim1901/flames_v2/flames_output/remove_unknownstrand_classification.txt", sep ='\t')

# ======================================================================================================================
#Aggregate the expression data by cell type 
counts <- AggregateExpression(all_samples_integrated, 
                              assays = "iso", 
                              return.seurat = FALSE,
                              group.by="fertility")
as.data.frame(counts) -> df
row.names(df) -> df$gene

#Split transcript ids into gene and transcript id
pseudobulk_data <- df %>% separate(gene, into = c("transcript_id", "gene_id"), sep = "-",  extra = "merge") 

#Use the pseudobulk_data counts calculated previously 
#Select cell type columns and gather into long format
long_data <- pseudobulk_data %>%
  pivot_longer(cols = starts_with("iso."),
               names_to = "fertility",
               values_to = "expression")

#Filter for non-zero expression to identify genes and isoforms expressed in each cell type
filtered_data <- long_data %>%
                 filter(expression > 0) %>%
                 dplyr::select(fertility, transcript_id, gene_id) %>%
                 distinct()

#Calculate unique gene and isoform counts for each cell type
cell_type_summary <- filtered_data %>%
  group_by(fertility) %>%
  summarise(num_genes = n_distinct(gene_id),
            num_isoforms = n_distinct(transcript_id)) %>%
  pivot_longer(cols = c(num_genes, num_isoforms), 
               names_to = "Category", 
               values_to = "Count")

# ======================================================================================================================
#Plot gene and isoform count
# Plot gene and isoform count with counts on top
p1 <- ggplot(cell_type_summary, aes(x = fertility, y = Count, fill = Category)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.7), width = 0.6) +
  #geom_text(aes(label = Count), position = position_dodge(width = 0.8), vjust = -0.5, size = 4) +  #Adds number counts over bars
  scale_x_discrete(labels = c("iso.Fertile" = "Fertile", "iso.Infertile" = "Infertile")) +
  scale_fill_manual(values = c("num_genes" = "#849AFF", "num_isoforms" = "#CBD2FF"),
                    labels = c("num_genes" = "Genes", "num_isoforms" = "Isoforms")) +
  theme_minimal() +
  labs(title = "Total Number of Genes\n and Isoforms by Fertility",
       y = "Count") +
  theme(plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
        axis.title.x = element_text(size = 0),
        axis.title.y = element_text(size = 14),
        axis.text.x = element_text(size = 14, angle = 0, hjust = 0.5, face = "bold"),
        axis.text.y = element_text(size = 12),
        legend.title = element_blank(),
        legend.text = element_text(size = 12))
p1


# ======================================================================================================================
#Plot the structural category across fertility
merged_data <- merge(pseudobulk_data, SQANTI, 
                     by.x = "transcript_id", 
                     by.y = "isoform", 
                     all.x = TRUE)

#Pivot pseudobulk data to long format
long_data <- merged_data %>%
  pivot_longer(cols = starts_with("iso."),
               names_to = "fertility",
               values_to = "expression")

#Generate our new df that we can use for plotting attributes defined by fertility
filtered_data <- long_data %>%
                 filter(expression > 0) %>%
                 distinct()

#Calculate unique counts of genes and isoforms per structural category and fertility
category_summary <- filtered_data %>%
  group_by(fertility, structural_category, subcategory, coding) %>%
  summarise(num_genes = n_distinct(gene_id),
            num_isoforms = n_distinct(transcript_id),
            .groups = "drop")

#Plot number of isoforms per structural category for fertility
p2 <- ggplot(category_summary, aes(x = fertility, y = num_isoforms, fill = structural_category)) +
  geom_bar(stat = "identity", position = "stack", width = 0.8) +
  theme_minimal() +
  scale_x_discrete(labels = c("iso.Fertile" = "Fertile", "iso.Infertile" = "Infertile")) +
  labs(title = "Isoforms by Structural Category \n by Fertility",
       x = "Fertility", y = "Number of Isoforms", fill = "Structural Category") +
  theme(plot.title = element_text(size = 20, face = "bold", hjust = 0.5),
        axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 14),
        axis.text.x = element_text(size = 14, angle = 0, hjust = 0.5),
        axis.text.y = element_text(size = 14),
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 12))
p2

# ======================================================================================================================
#Plot 3 - Isoform structural subcategories
desired_order <- c("reference_match",
                   "mono-exon", 
                   "multi-exon", 
                   "at_least_one_novel_splicesite",
                   "combination_of_known_junctions",
                   "combination_of_known_splicesites", 
                   "intron_retention",
                   "NA")
category_summary$subcategory <- factor(category_summary$subcategory, 
                                       levels = desired_order)
levels(category_summary$subcategory) <- c("Reference Match","Mono-Exon","Multi-Exon",
                                          "At Least 1 Novel Splicesite","Combination of Known Junctions",
                                          "Combination of Known Splicesites","Intron Retention","Other")
levels(category_summary$subcategory)

#Plot
p3 <- ggplot(category_summary, aes(x = fertility, y = num_isoforms, fill = subcategory)) +
  geom_bar(stat = "identity", position = "stack", width = 0.8) +
  scale_x_discrete(labels = c("iso.Fertile" = "Fertile", "iso.Infertile" = "Infertile")) +
  theme_minimal() +
  labs(title = "Isoform Structural\n Subcategory by Fertility",
       y = "Number of Isoforms", fill = "Structural Subcategory") +
  theme(plot.title = element_text(size = 20, face = "bold", hjust = 0.5),
        axis.title.x = element_text(size = 0),
        axis.title.y = element_text(size = 14),
        axis.text.x = element_text(size = 14, angle = 0, hjust = 0.5, face = "bold"),
        axis.text.y = element_text(size = 14),
        legend.title = element_text(size = 14, face = "bold"),
        legend.text = element_text(size = 12))
p3

# ======================================================================================================================
#Plot coding vs non-coding isoforms
p4 <- ggplot(category_summary %>% filter(!is.na(coding)),  # Remove NA values
  aes(x = fertility, y = num_isoforms, fill = coding)) +
  geom_bar(stat = "identity", position = "fill", width = 0.8) +
  scale_fill_manual(values = c("coding" = "#00BFC4", "non_coding" ="#A1E6E8"),
                    labels = c("coding" = "Coding", "non_coding" = "Non-Coding")) +
  theme_minimal() +
  labs(title = "Proportion of Coding \n Isoforms Split by Fertility",
       y = "Percentage of Isoforms", fill = "Coding Status") +
  scale_y_continuous(labels = scales::percent) +  # Show y-axis as percentages
  scale_x_discrete(labels = c("iso.Fertile" = "Fertile", "iso.Infertile" = "Infertile")) +  # Rename x-axis labels
  theme(plot.title = element_text(size = 20, face = "bold", hjust = 0.5),
        axis.title.x = element_text(size = 0),
        axis.title.y = element_text(size = 14),
        axis.text.x = element_text(size = 14, angle = 0, hjust = 0.5, face = "bold"),
        axis.text.y = element_text(size = 14),
        legend.title = element_text(size = 14, face = "bold"),
        legend.text = element_text(size = 12))
p4


cowplot::plot_grid(p1, p2, p3, p4, ncol = 2)


# ======================================================================================================================
###Biomart analysis
#Remove version number form isoform
filtered_data$transcript_id <- sub("\\..*$", "", filtered_data$transcript_id)

#Retrieve gene biotype from Ensembl
mart <- useEnsembl(biomart = "ENSEMBL_MART_ENSEMBL",
                   dataset = "hsapiens_gene_ensembl",
                   host = "https://jan2020.archive.ensembl.org") 
biotype_info_transcript <- getBM(attributes = c('ensembl_transcript_id', 'transcript_biotype', 'transcript_length'),
                                 filters = 'ensembl_transcript_id',
                                 values = filtered_data$transcript_id, mart = mart)

#Merge exclusive isoforms and biomart data together 
merged_biomart_data <- merge(filtered_data, biotype_info_transcript, 
                             by.x = "transcript_id", by.y = "ensembl_transcript_id", 
                             all.x = TRUE)

#Reassign transcript_biotype to "BambuTx" if transcript_id starts with "BambuTx"
merged_biomart_data$transcript_id <- trimws(as.character(merged_biomart_data$transcript_id))
merged_biomart_data <- merged_biomart_data %>%
  mutate(transcript_biotype = ifelse(grepl("^BambuTx", transcript_id), "BambuTx", transcript_biotype))

#Prepare the summary data
category_summary_2 <- merged_biomart_data %>%
  group_by(fertility, transcript_biotype) %>%
  summarise(num_isoforms = n(), .groups = "drop") %>%
  group_by(fertility) %>%
  mutate(proportion = num_isoforms / sum(num_isoforms)) %>%
  ungroup()

#Group transcript_biotypes with proportion < 0.5% into "Other"
category_summary_2 <- category_summary_2 %>%
  mutate(transcript_biotype = ifelse(proportion < 0.005, "Other", transcript_biotype)) %>%
  group_by(fertility, transcript_biotype) %>%
  summarise(num_isoforms = sum(num_isoforms), .groups = "drop") %>%
  group_by(fertility) %>%
  mutate(proportion = num_isoforms / sum(num_isoforms)) %>%
  ungroup()

#Ensure consistent factor order
category_summary_2 <- category_summary_2 %>%
  mutate(fertility = factor(fertility, levels = unique(fertility)),
         transcript_biotype = factor(transcript_biotype, levels = c(
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

#Pie chart of structural categories
ggplot(category_summary_2, aes(x = "", y = proportion, fill = transcript_biotype)) +
  geom_bar(stat = "identity", width = 1, colour="black", size=0.2) +  
  coord_polar(theta = "y") +  # Convert to pie chart
  facet_wrap(~fertility, labeller = as_labeller(c("iso.Fertile" = "Fertile", "iso.Infertile" = "Infertile"))) +  
  theme_minimal() +
  scale_fill_discrete(labels = biotype_labels) +
  labs(title = "BioMart Structural Categories of Fertility Unique Isoforms",
       fill = "BioMart Structural Category") +
  theme(plot.title = element_text(size = 25, face = "bold", hjust = 0.5),
        axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        legend.title = element_text(size = 20),
        legend.text = element_text(size = 16),
        strip.text = element_text(size = 18, face = "bold"),
        panel.grid = element_blank())



# ======================================================================================================================
###Novel genes
#Filter novel genes
novel_genes <- pseudobulk_data %>%
  filter(grepl("^BambuGene\\d+$", gene_id))

#Determine detection status in each condition
novel_status <- novel_genes %>%
  mutate(detected_fertile = iso.Fertile > 0,
         detected_infertile = iso.Infertile > 0) %>%
  dplyr::select(gene_id, detected_fertile, detected_infertile) %>%
  as_tibble()

#Classify genes by detection status
gene_classification <- novel_status %>%
  mutate(status = case_when(
      detected_fertile & detected_infertile ~ "Common",
      detected_fertile & !detected_infertile ~ "Unique to Fertile",
      !detected_fertile & detected_infertile ~ "Unique to Infertile",
      TRUE ~ NA_character_)) %>%
  filter(!is.na(status))

#Reshape to long format for plotting: one row per (gene, fertility condition)
plot_data <- gene_classification %>%
  pivot_longer(cols = c(detected_fertile, detected_infertile),
               names_to = "fertility",
               values_to = "detected") %>%
  filter(detected) %>%
  mutate(fertility = ifelse(fertility == "detected_fertile", "Fertile", "Infertile")) %>%
  mutate(category = case_when(
      status == "Common" ~ "Common",
      status == "Unique to Fertile" & fertility == "Fertile" ~ "Unique to Fertile",
      status == "Unique to Infertile" & fertility == "Infertile" ~ "Unique to Infertile",
      TRUE ~ NA_character_)) %>%
  filter(!is.na(category))

#Count genes per fertility and category
counts <- plot_data %>%
  group_by(fertility, category) %>%
  summarise(n = n_distinct(gene_id), .groups = "drop") %>%
  mutate(category = factor(category, levels = c("Common", "Unique to Fertile", "Unique to Infertile")))

#Calculate percentages within each fertility group
counts <- counts %>%
  group_by(fertility) %>%
  mutate(percent = n / sum(n) * 100,
    ypos = cumsum(n) - 0.5 * n) %>%  # for label positioning in the middle of the bar segment
  ungroup()

#Reverse factor levels to put unique on top in the stacked bar
counts$category <- factor(counts$category, levels = rev(levels(counts$category)))

#Plot stacked bar with percentages
ggplot(counts, aes(x = fertility, y = n, fill = category)) +
  geom_bar(stat = "identity") +
  geom_text(aes(y = ypos, label = paste0(round(percent, 1), "%")),
            color = "black", size = 4) +
  scale_fill_manual(values = c(
    "Unique to Fertile" = "palegreen3",
    "Common" = "lightblue",
    "Unique to Infertile" = "#F8766D")) +
  labs(title = "Novel Genes", x = "Fertility Condition", y = "Number of Novel Genes", fill = "Category") +
  theme_minimal(base_size = 14)


# ======================================================================================================================
###Number of isoforms by fertility
#Aggregate isoform expression by cell_type and fertility
counts <- AggregateExpression(all_samples_integrated, 
                              assays = "iso", 
                              return.seurat = FALSE,
                              group.by = c("cell_type", "fertility"))

#Convert to dataframe
df <- as.data.frame(counts)
df$gene <- rownames(df)

#Split transcript-gene info
pseudobulk_data <- df %>%
  separate(gene, into = c("transcript_id", "gene_id"), sep = "-", extra = "merge")

#Pivot longer to get fertility group info
long_df <- pseudobulk_data %>%
  pivot_longer(cols = -c(transcript_id, gene_id),
               names_to = "celltype_fertility",
               values_to = "expression") %>%
  #Extract fertility info from the column name
  separate(celltype_fertility, into = c("cell_type", "fertility"), sep = "_")

#Filter for expressed isoforms (optional threshold, e.g. expression > 0)
long_df_filtered <- long_df %>%
  filter(expression > 0)

#Count isoforms per gene per fertility group, and assign category
isoform_count_per_gene <- long_df_filtered %>%
  dplyr::group_by(gene_id, fertility) %>%
  dplyr::summarise(n_isoforms = n_distinct(transcript_id), .groups = "drop") %>%
  dplyr::mutate(
    isoform_category = dplyr::case_when(
      n_isoforms == 1 ~ "1",
      n_isoforms >= 2 & n_isoforms <= 3 ~ "2-3",
      n_isoforms >= 4 & n_isoforms <= 5 ~ "4-5",
      n_isoforms >= 6 & n_isoforms <= 9 ~ "6-9",
      n_isoforms >= 10 ~ "10+")) %>%
  dplyr::ungroup() %>%
  as_tibble()

#Calculate percentage of genes per bin per fertility group
isoform_count_summary <- isoform_count_per_gene %>%
  dplyr::count(fertility, isoform_category) %>%
  dplyr::group_by(fertility) %>%
  dplyr::mutate(percent = (n / sum(n)) * 100) %>%
  dplyr::ungroup()

#Ensure consistent order for x-axis
isoform_count_summary$isoform_category <- factor(
  isoform_count_summary$isoform_category,
  levels = c("1", "2-3", "4-5", "6-9", "10+"))

ggplot(isoform_count_summary, aes(x = isoform_category, y = percent, fill = fertility)) +
  geom_bar(stat = "identity", 
           position = position_dodge(width = 0.8), 
           color = "black", width = 0.7) +  
  labs(title = "Number of Isoforms per\n Gene by Fertility", 
       x = "Isoforms per Gene", y = "Gene Percentage", fill = "Fertility") +
  scale_fill_manual(values = c("Fertile" = "palegreen3", "Infertile" = "#F8766D")) +
  theme_minimal() +
  theme(plot.title = element_text(size = 16, hjust = 0.5, face = "bold"),
        axis.text.y = element_text(size = 12),
        axis.text.x = element_text(size = 12),
        axis.title.y = element_text(size = 13),
        axis.title.x = element_text(size = 13),
        legend.title = element_text(size = 13),
        legend.text = element_text(size = 12),
        legend.position = "bottom")
