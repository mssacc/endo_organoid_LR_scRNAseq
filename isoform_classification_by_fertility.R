###Isoform classification by fertility
# ======================================================================================================================

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
       x = "Fertility",
       y = "Number of Isoforms",
       fill = "Structural Category") +
  theme(plot.title = element_text(size = 20, face = "bold", hjust = 0.5),
        axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 14),
        axis.text.x = element_text(size = 14, angle = 0, hjust = 0.5),
        axis.text.y = element_text(size = 14),
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 12))
p2

###Pie chart of biotypes
#Calculate the total number of isoforms across all categories for each fertility level
category_summary_percent <- category_summary %>%
  group_by(fertility, structural_category) %>%
  summarise(num_isoforms = sum(num_isoforms)) %>%
  ungroup() %>%
  group_by(fertility) %>%
  mutate(percentage = num_isoforms / sum(num_isoforms)) %>%
  ungroup()

#Define the labeller function for fertility levels
fertility_labeller <- c("iso.Fertile" = "Fertile", "iso.Infertile" = "Infertile")

desired_order <- c("full-splice_match",
                   "incomplete-splice_match", 
                   "novel_in_catalog", 
                   "novel_not_in_catalog",
                   "antisense",
                   "fusion", 
                   "genic", 
                   "genic_intron", 
                   "NA")
category_summary_percent$structural_category <- factor(category_summary_percent$structural_category, 
                                                       levels = desired_order)
levels(category_summary_percent$structural_category) <- c("FSM", "ISM", "NIC", "NNC", "Antisense",
                                                  "Fusion","Genic", "Genic Intron", "Other")
levels(category_summary_percent$structural_category)

# Create pie chart with percentage labels for each fertility status
p2_pie <- ggplot(category_summary_percent, aes(x = "", y = percentage, fill = structural_category)) +
          geom_bar(stat = "identity", width = 1, colour="black", linewidth=0.1) +  
          coord_polar(theta = "y", start = pi / 6) +
          facet_wrap(~ fertility, labeller = labeller(fertility = fertility_labeller)) +
          labs(title = "Isoforms by Structural Category Across Fertility",
               fill = "Structural Category") +
          theme_void() +
          theme(plot.title = element_text(size = 25, face = "bold", hjust = 0.5, margin = margin(b = 20)),
                strip.text = element_text(size = 18),
                legend.title = element_text(size = 20),
                legend.text = element_text(size = 16))
p2_pie

p2_bar <- ggplot(category_summary_percent, aes(x = fertility, y = percentage, fill = structural_category)) +
  geom_bar(stat = "identity", position = "stack", width = 0.8) +
  theme_minimal() +
  scale_x_discrete(labels = c("iso.Fertile" = "Fertile", "iso.Infertile" = "Infertile")) +
  scale_y_continuous(labels = scales::percent) +  # Show y-axis as percentages
  labs(title = "Isoform Structural\n Category by Fertility",
       y = "Percentage of Isoforms",
       fill = "Structural Category") +
  theme(plot.title = element_text(size = 20, face = "bold", hjust = 0.5),
        axis.title.x = element_text(size = 0),
        axis.title.y = element_text(size = 14),
        axis.text.x = element_text(size = 14, angle = 0, hjust = 0.5, face = "bold"),
        axis.text.y = element_text(size = 14),
        legend.title = element_text(size = 14, face = "bold"),
        legend.text = element_text(size = 12))
p2_bar
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
       y = "Number of Isoforms",
       fill = "Structural Subcategory") +
  theme(plot.title = element_text(size = 20, face = "bold", hjust = 0.5),
        axis.title.x = element_text(size = 0),
        axis.title.y = element_text(size = 14),
        axis.text.x = element_text(size = 14, angle = 0, hjust = 0.5, face = "bold"),
        axis.text.y = element_text(size = 14),
        legend.title = element_text(size = 14, face = "bold"),
        legend.text = element_text(size = 12))
p3

# ======================================================================================================================
p4 <- ggplot(category_summary %>% filter(!is.na(coding)),  # Remove NA values
  aes(x = fertility, y = num_isoforms, fill = coding)) +
  geom_bar(stat = "identity", position = "fill", width = 0.8) +
  scale_fill_manual(values = c("coding" = "#00BFC4", "non_coding" ="#A1E6E8"),
                    labels = c("coding" = "Coding", "non_coding" = "Non-Coding")) +
  theme_minimal() +
  labs(title = "Proportion of Coding \n Isoforms Split by Fertility",
       y = "Percentage of Isoforms",
       fill = "Coding Status") +
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

#Plots
cowplot::plot_grid(p1, p2, p3, p4, ncol = 2)

p2_bar + p4 + p3



# ======================================================================================================================
# ======================================================================================================================
# ======================================================================================================================

###Isoforms unique to fertility status
#We define unique as an isoform that has 5+ counts in one cell type and ≤ 1 count in all other cell types
#Filter the data based on the expression cutoff
#Set expression cutoff threshold
expression_cutoff <- 5
max_expression_in_all_other_cells_types = 1
#Filter isoforms based on the new criteria
exclusive_isoforms <- long_data %>%
  group_by(transcript_id) %>%
  filter(sum(expression > expression_cutoff) <= 1,     # Only one cell type where expression > cutoff
    # Ensure all other cell types have 0 expression
    all(expression < max_expression_in_all_other_cells_types | expression > expression_cutoff)) %>%
  ungroup()

#Count the unique isoforms per cell type
#Adjust counting logic
exclusive_isoforms_count <- exclusive_isoforms %>%
  group_by(fertility) %>%
  summarise(unique_isoforms = sum(expression > expression_cutoff)) %>%  # Count only isoforms with valid expression
  ungroup()

#Plot the results
ggplot(exclusive_isoforms_count, aes(x = fertility, y = unique_isoforms, fill = fertility)) +
  geom_bar(stat = "identity", color = "black", width = 0.8) +
  theme_minimal() +
  labs(title = "Number of Isoforms Unique \n to Fertility",
       x = "Fertility", y = "Number of Unique Isoforms") +
  scale_x_discrete(labels = c("iso.Fertile" = "Fertile", "iso.Infertile" = "Infertile")) +  # Rename x-axis labels
  scale_fill_discrete(labels = c("iso.Fertile" = "Fertile", "iso.Infertile" = "Infertile")) +  # Change legend labels
  theme(plot.title = element_text(size = 20, face = "bold", hjust = 0.5),
        axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 14),
        axis.text.x = element_text(size = 12, angle = 0, hjust = 0.5),
        axis.text.y = element_text(size = 12),
        legend.title = element_blank(),
        legend.text = element_text(size = 14))
    
#We could extract these isoforms for further analysis
exclusive_isoforms_uniq <- exclusive_isoforms %>%
  filter(expression >= max_expression_in_all_other_cells_types) %>%
  arrange(desc(expression))
exclusive_isoforms_uniq$transcript_id <- sub("\\..*", "", exclusive_isoforms_uniq$transcript_id)


#Retrieve gene biotype from Ensembl
mart <- useMart(biomart = "ensembl", dataset = "hsapiens_gene_ensembl") 
biotype_info_transcript <- getBM(attributes = c("ensembl_transcript_id", "transcript_is_canonical", 'transcript_biotype', 'transcript_length'),
                                 filters = 'ensembl_transcript_id',
                                 values = exclusive_isoforms_uniq$transcript_id, mart = mart)

#Merge exclusive isoforms and biomart data together 
merged_biomart_data <- merge(exclusive_isoforms_uniq, biotype_info_transcript, 
                             by.x = "transcript_id", by.y = "ensembl_transcript_id", 
                             all.x = TRUE)

#Print out unique isoforms
write.csv(merged_biomart_data, file="/data/gpfs/projects/punim1901/flames_v2/seurat_workspace/isoforms_unique_to_fertility.csv")

#Prepare the summary data
category_summary_2 <- merged_biomart_data %>%
  group_by(fertility, transcript_biotype) %>%
  summarise(num_isoforms = n(), .groups = "drop") %>%
  group_by(fertility) %>%
  mutate(proportion = num_isoforms / sum(num_isoforms)) %>%
  ungroup()

#Identify biotypes present in your data (use to update code below)
unique(merged_biomart_data$transcript_biotype)

#Ensure consistent factor order for cell_type and transcript_biotype
category_summary_2 <- category_summary_2 %>%
  mutate(fertility = factor(fertility, levels = unique(fertility)),
         transcript_biotype = factor(transcript_biotype, levels = c(
            "protein_coding", "lncRNA", "retained_intron", "protein_coding_CDS_not_defined",
            "nonsense_mediated_decay", "processed_pseudogene", "unprocessed_pseudogene",
            "transcribed_unprocessed_pseudogene", "unitary_pseudogene", "misc_RNA", "TEC", "snRNA", "artifact")))

#Create the plot 
ggplot(category_summary_2, aes(x = fertility, y = proportion, fill = transcript_biotype)) +
  geom_bar(stat = "identity", position = "fill", width = 0.8) +
  geom_text(aes(label = scales::percent(proportion, accuracy = 0.1)), 
            position = position_fill(vjust = 0.5), size = 3) +
  theme_minimal() +
  scale_y_continuous(labels = scales::percent) +
  scale_x_discrete(labels = c("iso.Fertile" = "Fertile", "iso.Infertile" = "Infertile")) +  # Rename x-axis labels
  theme(plot.title = element_text(size = 22, face = "bold", hjust = 0.5),
        axis.title.x = element_text(size = 0),
        axis.title.y = element_text(size = 15),
        axis.text.x = element_text(size = 15, angle = 0, hjust = 0.5),
        axis.text.y = element_text(size = 15),
        legend.title = element_text(size = 18),
        legend.text = element_text(size = 15))

#Pie chart of transcript subtypes
ggplot(category_summary_2, aes(x = "", y = proportion, fill = transcript_biotype)) +
  geom_bar(stat = "identity", width = 1, colour="black", size=0.2) +  
  coord_polar(theta = "y") +  # Convert to pie chart
  facet_wrap(~fertility, labeller = as_labeller(c("iso.Fertile" = "Fertile", "iso.Infertile" = "Infertile"))) +  
  theme_minimal() +
  labs(title = "BioMart Structural Categories of Fertility Unique Isoforms",
       fill = "BioMart Structural Category") +
  theme(plot.title = element_text(size = 25, face = "bold", hjust = 0.5),
        axis.title = element_blank(),  # Remove x & y axis labels
        axis.text = element_blank(),   # Remove axis text
        axis.ticks = element_blank(),  # Remove axis ticks
        legend.title = element_text(size = 20),
        legend.text = element_text(size = 16),
        strip.text = element_text(size = 18, face = "bold"),
        panel.grid = element_blank())


#Get info for types of biotypes in all cells
#Remove the ".number" from the 'isoform' column
filtered_data$transcript_id <- sub("\\..*", "", filtered_data$transcript_id)
mart <- useMart(biomart = "ensembl", 
                dataset = "hsapiens_gene_ensembl")

biotype_info_all <- getBM(attributes = c("ensembl_transcript_id", "transcript_is_canonical", 'transcript_biotype', 'transcript_length'),
                          filters = 'ensembl_transcript_id',
                          values = filtered_data$transcript_id, mart = mart)

#Reassign transcript_biotype for BambuTx IDs
merged_data <- merged_data %>%
  mutate(transcript_biotype = ifelse(grepl("^BambuTx", transcript_id), "BambuTx", transcript_biotype))

#Print out unique isoforms
write.csv(merged_data, file="/data/gpfs/projects/punim1901/flames_v2/seurat_workspace/biomart_split_by_fertility.csv")

#Prepare the summary data
category_summary_3 <- merged_data %>%
  group_by(fertility, transcript_biotype) %>%
  summarise(num_isoforms = n(), .groups = "drop") %>%
  group_by(fertility) %>%
  mutate(proportion = num_isoforms / sum(num_isoforms)) %>%
  ungroup()

#Identify biotypes present in your data (use to update code below)
unique(merged_data$transcript_biotype)

#Merge categories less than 0.1 into 1 category
category_summary_3 <- category_summary_3 %>%
  mutate(transcript_biotype = ifelse(proportion < 0.001, "Other", as.character(transcript_biotype))) %>%
  group_by(fertility, transcript_biotype) %>%
  summarise(num_isoforms = sum(num_isoforms), proportion = sum(proportion), .groups = "drop") %>%
  ungroup()
category_summary_3$transcript_biotype <- factor(category_summary_3$transcript_biotype)


#Ensure consistent factor order for cell_type and transcript_biotype
category_summary_3 <- category_summary_3 %>%
  mutate(fertility = factor(fertility, levels = unique(fertility)),
         transcript_biotype = factor(transcript_biotype, levels = c(
            "protein_coding", "lncRNA", "retained_intron", "protein_coding_CDS_not_defined",
            "nonsense_mediated_decay", "processed_pseudogene", "TEC",
            "transcribed_unprocessed_pseudogene", "misc_RNA", "transcribed_processed_pseudogene",
            "unprocessed_pseudogene", "snRNA", "snoRNA", "non_stop_decay", "NA", "Other")))

#Create the plot
#Pie chart of transcript subtypes
ggplot(category_summary_3, aes(x = "", y = proportion, fill = transcript_biotype)) +
  geom_bar(stat = "identity", width = 1, colour="black", size=0.2) +  
  coord_polar(theta = "y") +  # Convert to pie chart
  facet_wrap(~fertility, labeller = as_labeller(c("iso.Fertile" = "Fertile", "iso.Infertile" = "Infertile"))) +  
  theme_minimal() +
  labs(title = "BioMart Structural Categories by Fertility",
        fill = "BioMart Structural Category") +
  theme(plot.title = element_text(size = 25, face = "bold", hjust = 0.5),
        axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        legend.title = element_text(size = 20),
        legend.text = element_text(size = 16),
        strip.text = element_text(size = 18, face = "bold"),
        panel.grid = element_blank())




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
  pivot_longer(
    cols = c(detected_fertile, detected_infertile),
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
  mutate(
    category = factor(category, levels = c("Common", "Unique to Fertile", "Unique to Infertile")))

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
  labs(title = "Novel Genes",
       x = "Fertility Condition", y = "Number of Novel Genes",
       fill = "Category") +
  theme_minimal(base_size = 14)







#lets aggeragte the expresstion data by cell type 
counts <- AggregateExpression(
  all_samples_integrated, 
  assays = "iso", 
  return.seurat = FALSE,
  group.by = "cell_type")

as.data.frame(counts) -> df
row.names(df) -> df$gene

#split transcript ids into gene and transcript id
pseudobulk_data <- df %>% separate(gene, into = c("transcript_id", "gene_id"), sep = "-",  extra = "merge") 
#df$transcript_id <- sub("\\..*", "", df$transcript_id)

# 2. Count the number of isoforms per gene
isoform_count_per_gene <- pseudobulk_data %>%
  group_by(gene_id) %>%
  summarise(n_isoforms = n_distinct(transcript_id))

# 3. count isforms per category 
isoform_count_per_gene <- isoform_count_per_gene %>%
  mutate(isoform_category = case_when(
    n_isoforms == 1 ~ "1",
    n_isoforms >= 2 & n_isoforms <= 3 ~ "2-3",
    n_isoforms >= 4 & n_isoforms <= 5 ~ "4-5",
    n_isoforms >= 6 ~ "6+"))

# 4. Calculate the percentage of genes in each bin
isoform_count_summary <- isoform_count_per_gene %>%
  dplyr::count(isoform_category) %>%
  mutate(percent = (n / sum(n)) * 100)

ggplot(isoform_count_summary, aes(x = isoform_category, y = percent)) +
  geom_bar(stat = "identity", fill = "#849AFF", color = "black") +  
  labs(title = "Number of Isoforms per Gene", 
       x = "Isoforms per Gene", 
       y = "Genes Percentage") +
  theme_minimal() +
  theme(legend.position = "none", 
        plot.title = element_text(hjust = 0.5, face = "bold"),
        axis.text.y = element_text(size = 11),
        axis.text.x = element_text(size = 11),
        axis.title.y = element_text(size = 12),
        axis.title.x = element_text(size = 12))


####Nmber of isoforms by fertility
#Aggregate isoform expression by cell_type and fertility
counts <- AggregateExpression(
  all_samples_integrated, 
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
  pivot_longer(
    cols = -c(transcript_id, gene_id),
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
       x = "Isoforms per Gene", 
       y = "Gene Percentage", 
       fill = "Fertility") +
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
